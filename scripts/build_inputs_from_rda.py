"""
Build Troy-format ERA5 input files for one month directly from NCAR's RDA mirror
of ERA5 (no CDS download needed).

RDA dataset paths used:
  /glade/campaign/collections/rda/data/d633000/e5.oper.an.pl/<YYYYMM>/
    e5.oper.an.pl.128_130_t.ll025sc.<YYYYMMDD>00_<YYYYMMDD>23.nc  (one per day, full 37 lvls)
    e5.oper.an.pl.128_131_u.ll025uv.<YYYYMMDD>00_<YYYYMMDD>23.nc
    e5.oper.an.pl.128_132_v.ll025uv.<YYYYMMDD>00_<YYYYMMDD>23.nc
    e5.oper.an.pl.128_133_q.ll025sc.<YYYYMMDD>00_<YYYYMMDD>23.nc
  /glade/campaign/collections/rda/data/d633000/e5.oper.an.sfc/<YYYYMM>/
    e5.oper.an.sfc.128_134_sp.ll025sc.<YYYYMM>0100_<YYYYMM>{31}23.nc  (one per month)

Method:
  1. RDA stores at native 0.25° (1440 lon × 721 lat × 37 levels).
  2. Sub-sample to 3.0° by taking every 12th point (matches CDS `--grid 3.0 3.0`,
     which is point-sampling, not area-averaging — verified bit-identical
     against current CDS at multiple test points).
  3. Sub-set the 28 levels Troy uses (20-975, drops 1/2/3/5/7/10/15/1000).
  4. Write era_5_m{M}_y{YYYY}_full.nc and era_5_m{M}_y{YYYY}_surface_p.nc
     with the dim order Fortran expects (time, level, latitude, longitude)
     in CDL terms / (longitude, latitude, level, time) in Fortran terms.

Output:
  inputs/<year>/era_5_m{M}_y{YYYY}_full.nc
  inputs/<year>/era_5_m{M}_y{YYYY}_surface_p.nc

Usage:
  python build_inputs_from_rda.py --year 2014 --month 1
  python build_inputs_from_rda.py --year 2014 --month 1 --skip-if-exists
"""
from __future__ import annotations

import argparse
import calendar
import os
import sys
import time as _time

import numpy as np
import netCDF4 as nc

# Force unbuffered output so progress is visible in long-running tasks
sys.stdout.reconfigure(line_buffering=True)
sys.stderr.reconfigure(line_buffering=True)


WORKDIR = os.environ.get("LUCIE_WORKDIR",
    "/glade/derecho/scratch/mdarman/sst_data/era5_data_download/gen_lucie_3d")
RDA_PL = "/glade/campaign/collections/rda/data/d633000/e5.oper.an.pl"      # NCAR RDA mirror
RDA_SFC = "/glade/campaign/collections/rda/data/d633000/e5.oper.an.sfc"    # NCAR RDA mirror

# Variables and their RDA grid suffixes
PL_VARS = [
    ("t", "128_130", "ll025sc"),
    ("u", "128_131", "ll025uv"),
    ("v", "128_132", "ll025uv"),
    ("q", "128_133", "ll025sc"),
]

# Troy's 28 pressure levels (subset of ERA5's 37 native levels). Must be in same
# order as how Troy's archived files store them — checked: ascending by pressure.
TROY_LEVELS = [20, 30, 50, 70, 100, 125, 200, 225, 250, 300, 350, 400, 450, 500, 550,
               600, 650, 700, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975]

# Sub-sampling stride from 0.25° → 3.0°: 12 (since 1440/120 = 12)
SUBSAMPLE_STRIDE = 12


def rda_pl_filename(year: int, month: int, day: int, var_code: str, grid: str) -> str:
    yyyymm = f"{year}{month:02d}"
    yyyymmdd = f"{year}{month:02d}{day:02d}"
    return os.path.join(RDA_PL, yyyymm,
                        f"e5.oper.an.pl.{var_code}.{grid}.{yyyymmdd}00_{yyyymmdd}23.nc")


def rda_sfc_filename(year: int, month: int, var_code: str, grid: str) -> str:
    yyyymm = f"{year}{month:02d}"
    last_day = calendar.monthrange(year, month)[1]
    return os.path.join(RDA_SFC, yyyymm,
                        f"e5.oper.an.sfc.{var_code}.{grid}.{yyyymm}0100_{yyyymm}{last_day:02d}23.nc")


def find_var_name(ds: nc.Dataset, candidates: list[str]) -> str:
    for c in candidates:
        if c in ds.variables:
            return c
    raise KeyError(f"None of {candidates} found in dataset (have {list(ds.variables.keys())})")


def _read_one_file(task: tuple) -> tuple:
    """Worker for multiprocessing.Pool. Reads one daily file and returns the
    sub-sampled (24, n_levels, 61, 120) array as float64.

    task = (day, var_name, file_path, level_idx_list, subsample_stride)
    Returns: (day, var_name, ndarray)
    """
    day, var_name, fp, level_idx_list, stride = task
    level_idx = np.array(level_idx_list, dtype=int)
    with nc.Dataset(fp, "r") as ds:
        vname = find_var_name(ds, [var_name.upper(), var_name])
        # Strided slicing: avoids reading the full 721×1440 grid.
        arr = ds.variables[vname][:, level_idx, ::stride, ::stride]
    return day, var_name, arr.astype(np.float64)


def build_full_nc(year: int, month: int, out_dir: str) -> str:
    """Read T/U/V/Q from RDA daily files, sub-sample, and write a Troy-format
    monthly _full.nc."""
    out_path = os.path.join(out_dir, f"era_5_m{month}_y{year}_full.nc")
    if os.path.exists(out_path):
        os.remove(out_path)

    last_day = calendar.monthrange(year, month)[1]
    n_hours = last_day * 24

    print(f"  Building {out_path} ({last_day} days = {n_hours} hourly steps)")

    # Discover grid metadata from the first file
    first_t = rda_pl_filename(year, month, 1, "128_130_t", "ll025sc")
    if not os.path.exists(first_t):
        raise FileNotFoundError(f"RDA T file missing: {first_t}")

    with nc.Dataset(first_t, "r") as ds:
        rda_lat_full = ds.variables["latitude"][:]   # 721 N→S
        rda_lon_full = ds.variables["longitude"][:]  # 1440 0..359.75
        rda_levels_full = ds.variables["level"][:]   # 37 levels in hPa

    # Identify level indices for the 28 Troy levels
    level_idx = []
    for L in TROY_LEVELS:
        idx = int(np.where(rda_levels_full == L)[0][0])
        level_idx.append(idx)
    level_idx = np.array(level_idx, dtype=int)
    print(f"  Selecting {len(level_idx)} levels: {[int(rda_levels_full[i]) for i in level_idx]}")

    # Lat/lon sub-sampling indices
    lat_idx = np.arange(0, 721, SUBSAMPLE_STRIDE)  # 61 points
    lon_idx = np.arange(0, 1440, SUBSAMPLE_STRIDE)  # 120 points
    sub_lat = rda_lat_full[lat_idx].astype(np.float64)
    sub_lon = rda_lon_full[lon_idx].astype(np.float64)
    sub_lev = rda_levels_full[level_idx].astype(np.float64)
    nlat = len(sub_lat)
    nlon = len(sub_lon)
    nlev = len(sub_lev)
    print(f"  Sub-sampled grid: lon={nlon}, lat={nlat}, level={nlev}")

    # Allocate output arrays
    out_data = {}
    for var_name, _, _ in PL_VARS:
        # Shape (time, level, lat, lon) — same as Python C-order in the output file
        out_data[var_name] = np.empty((n_hours, nlev, nlat, nlon), dtype=np.float64)

    # Read 124 files (4 vars × 31 days) in parallel via multiprocessing.
    # Each file is independent. Decompression is single-threaded per file, so
    # parallelizing across files gives near-linear speedup for the cold-cache
    # case (the bottleneck is /glade/campaign network I/O + HDF5 chunk decompress).
    from multiprocessing import Pool
    n_workers = int(os.environ.get("RDA_BUILD_WORKERS", "4"))
    print(f"  Reading 124 files with {n_workers} worker processes...")

    tasks = []
    for day in range(1, last_day + 1):
        for var_name, var_code_suffix, grid in PL_VARS:
            full_code = f"{var_code_suffix}_{var_name}"
            fp = rda_pl_filename(year, month, day, full_code, grid)
            if not os.path.exists(fp):
                raise FileNotFoundError(f"RDA file missing: {fp}")
            tasks.append((day, var_name, fp, level_idx.tolist(), SUBSAMPLE_STRIDE))

    t_start = _time.monotonic()
    n_done = 0
    n_total = len(tasks)
    with Pool(n_workers) as pool:
        for result in pool.imap_unordered(_read_one_file, tasks, chunksize=1):
            day, var_name, arr = result
            hour_offset = (day - 1) * 24
            out_data[var_name][hour_offset:hour_offset + 24, :, :, :] = arr
            n_done += 1
            if n_done % 8 == 0 or n_done == n_total:
                elapsed = _time.monotonic() - t_start
                rate = n_done / elapsed if elapsed > 0 else 0
                eta = (n_total - n_done) / rate if rate > 0 else 0
                print(f"    {n_done:3d}/{n_total} files done ({elapsed:.0f}s elapsed, ETA {eta:.0f}s)")

    print("  Writing NetCDF output...")

    # Write Fortran-expected format. Python (C order) declares dim tuple as
    # ("time", "level", "latitude", "longitude") so Fortran (column-major)
    # sees (longitude, latitude, level, time) — what regrid_era.f90 expects.
    with nc.Dataset(out_path, "w", format="NETCDF4") as out:
        out.createDimension("longitude", nlon)
        out.createDimension("latitude", nlat)
        out.createDimension("level", nlev)
        out.createDimension("time", n_hours)

        v_lon = out.createVariable("longitude", "f8", ("longitude",))
        v_lat = out.createVariable("latitude", "f8", ("latitude",))
        v_lev = out.createVariable("level", "f8", ("level",))
        v_tim = out.createVariable("time", "f8", ("time",))
        v_lon[:] = sub_lon
        v_lat[:] = sub_lat
        v_lev[:] = sub_lev
        v_tim[:] = np.arange(n_hours, dtype=np.float64)
        v_lon.units = "degrees_east"
        v_lat.units = "degrees_north"
        v_lev.units = "millibars"
        v_tim.units = "hours since 1900-01-01"

        for var_name in ("t", "u", "v", "q"):
            v = out.createVariable(
                var_name, "f8", ("time", "level", "latitude", "longitude"),
                zlib=False,
            )
            v[:] = out_data[var_name]
            # Identity packing attrs (Fortran's read_netcdf_4d_dp_era requires them)
            v.scale_factor = 1.0
            v.add_offset = 0.0
            v.units = {"t": "K", "u": "m s**-1", "v": "m s**-1", "q": "kg kg**-1"}[var_name]

    print(f"  -> {out_path}")
    return out_path


def build_surface_p_nc(year: int, month: int, out_dir: str) -> str:
    out_path = os.path.join(out_dir, f"era_5_m{month}_y{year}_surface_p.nc")
    if os.path.exists(out_path):
        os.remove(out_path)

    fp = rda_sfc_filename(year, month, "128_134_sp", "ll025sc")
    if not os.path.exists(fp):
        raise FileNotFoundError(f"RDA SP file missing: {fp}")

    print(f"  Reading {fp}")
    with nc.Dataset(fp, "r") as ds:
        rda_lat_full = ds.variables["latitude"][:]
        rda_lon_full = ds.variables["longitude"][:]
        sp_full = ds.variables["SP"][:]   # (n_hours, 721, 1440)

    lat_idx = np.arange(0, 721, SUBSAMPLE_STRIDE)
    lon_idx = np.arange(0, 1440, SUBSAMPLE_STRIDE)
    sp_sub = sp_full[:, lat_idx, :][:, :, lon_idx].astype(np.float64)
    sub_lat = rda_lat_full[lat_idx].astype(np.float64)
    sub_lon = rda_lon_full[lon_idx].astype(np.float64)
    n_hours = sp_sub.shape[0]
    nlat = len(sub_lat)
    nlon = len(sub_lon)
    print(f"  Sub-sampled SP: lon={nlon}, lat={nlat}, time={n_hours}")

    with nc.Dataset(out_path, "w", format="NETCDF4") as out:
        out.createDimension("longitude", nlon)
        out.createDimension("latitude", nlat)
        out.createDimension("time", n_hours)
        v_lon = out.createVariable("longitude", "f8", ("longitude",))
        v_lat = out.createVariable("latitude", "f8", ("latitude",))
        v_tim = out.createVariable("time", "f8", ("time",))
        v_lon[:] = sub_lon
        v_lat[:] = sub_lat
        v_tim[:] = np.arange(n_hours, dtype=np.float64)
        v_sp = out.createVariable("sp", "f8", ("time", "latitude", "longitude"))
        v_sp[:] = sp_sub
        v_sp.scale_factor = 1.0
        v_sp.add_offset = 0.0
        v_sp.units = "Pa"

    print(f"  -> {out_path}")
    return out_path


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--year", type=int, required=True)
    p.add_argument("--month", type=int, required=True)
    p.add_argument("--inputs-dir", default=os.path.join(WORKDIR, "inputs"))
    p.add_argument("--skip-if-exists", action="store_true")
    args = p.parse_args()

    out_dir = os.path.join(args.inputs_dir, f"{args.year}")
    os.makedirs(out_dir, exist_ok=True)

    full_path = os.path.join(out_dir, f"era_5_m{args.month}_y{args.year}_full.nc")
    sp_path = os.path.join(out_dir, f"era_5_m{args.month}_y{args.year}_surface_p.nc")

    if args.skip_if_exists and os.path.exists(full_path) and os.path.exists(sp_path):
        print(f"[skip-all] Both files exist:\n  {full_path}\n  {sp_path}")
        return

    print(f"=== Build inputs from RDA for y={args.year} m={args.month} ===")
    build_full_nc(args.year, args.month, out_dir)
    build_surface_p_nc(args.year, args.month, out_dir)
    print("\n=== Done ===")
    print(f"  {full_path}")
    print(f"  {sp_path}")


if __name__ == "__main__":
    main()
