"""
Convert Troy's existing yearly NetCDF (`era_5_y{YYYY}_regridded_mpi.nc` or
`_fixed_var_gcc.nc`) directly into per-hour LUCIE-format h5 files.

For years 1950-2021 we already have Troy's bit-perfect output at:
  /glade/derecho/scratch/mdarman/ERA5_hr_haiwen/ERA5_T30_1950_TO_2021/<year>/era_5_y<year>_regridded_mpi.nc

These ARE the gold standard (verified bit-identical to LUCIE training files
in session 3). So for those years, we don't need to re-run the regridder at
all — we just slice the existing yearly file at each hour and write per-hour
h5 files.

This is MUCH faster than running the regridder from scratch (~1-2 minutes per
month vs ~13 hours from RDA cold-cache).

Usage:
  python extract_h5_from_yearly.py --year 2014
  python extract_h5_from_yearly.py --year 2014 --month 1     # only Jan
  python extract_h5_from_yearly.py --year 2014 --skip-existing
  python extract_h5_from_yearly.py --year-range 1990 2014    # multiple years

Each output:
  h5_out/{YYYY}_{HHHH:04d}.h5 — LUCIE-format file matching the training data
  for the prognostic variables (T_0..7, U_0..7, V_0..7, Q_0..7, logp).
  Plus boundary conditions (sst, tisr, orography, land_sea_mask) when available.
"""
from __future__ import annotations

import argparse
import calendar
import os
import sys
import time as _time
from datetime import datetime, timezone

import numpy as np
import h5py
import netCDF4 as nc

# Force unbuffered output
sys.stdout.reconfigure(line_buffering=True)
sys.stderr.reconfigure(line_buffering=True)

WORKDIR = os.environ.get("LUCIE_WORKDIR",
    "/glade/derecho/scratch/mdarman/sst_data/era5_data_download/gen_lucie_3d")

# Where Troy's yearly files live — change these for your setup.
# See RECIPE.md for how to obtain these files.
TROY_YEARLY_BASE = os.environ.get("LUCIE_TROY_YEARLY",
    "/glade/derecho/scratch/mdarman/ERA5_hr_haiwen/ERA5_T30_1950_TO_2021")

# Boundary-condition sources — change these for your setup
LUCIE_REF_BASE = os.environ.get("LUCIE_REF_BASE",
    "/glade/derecho/scratch/mdarman/ERA5_hr_haiwen/LUCIE_3D/1_step_1hr_h5df_test")
SST_DIR = os.environ.get("LUCIE_SST_DIR",
    "/glade/derecho/scratch/mdarman/sst_data/era5_data_download/sst_2014-2023")
TISR_DIR = os.environ.get("LUCIE_TISR_DIR",
    "/glade/derecho/scratch/mdarman/sst_data/era5_data_download/tisr_2015-2023")

STATIC_REF_H5_CANDIDATES = [
    os.path.join(LUCIE_REF_BASE, "test", "2013_0000.h5"),
    os.path.join(LUCIE_REF_BASE, "val", "2010_0000.h5"),
    os.path.join(LUCIE_REF_BASE, "train", "1990_0000.h5"),
]


def find_yearly_file(year: int) -> str:
    """Find Troy's yearly NC file for the given year."""
    candidates = [
        os.path.join(TROY_YEARLY_BASE, str(year), f"era_5_y{year}_regridded_mpi.nc"),
        os.path.join(TROY_YEARLY_BASE, str(year),
                     f"era_5_y{year}_regridded_mpi_fixed_var_gcc.nc"),
    ]
    for p in candidates:
        if os.path.exists(p):
            return p
    raise FileNotFoundError(
        f"No yearly file for {year}. Looked in:\n  " + "\n  ".join(candidates)
    )


def find_static_ref() -> str:
    for p in STATIC_REF_H5_CANDIDATES:
        if os.path.exists(p):
            return p
    raise FileNotFoundError(
        f"None of the static reference h5 files exist: {STATIC_REF_H5_CANDIDATES}"
    )


def find_lucie_ref_h5(year: int, hour_idx: int) -> str | None:
    fname = f"{year}_{hour_idx:04d}.h5"
    for sub in ("train", "val", "test"):
        p = os.path.join(LUCIE_REF_BASE, sub, fname)
        if os.path.exists(p):
            return p
    return None


def find_sst_h5(year: int, hour_idx: int) -> str | None:
    p = os.path.join(SST_DIR, f"{year}_{hour_idx:04d}.h5")
    return p if os.path.exists(p) else None


def find_tisr_h5(year: int, hour_idx: int) -> str | None:
    p = os.path.join(TISR_DIR, f"{year}_{hour_idx:04d}.h5")
    return p if os.path.exists(p) else None


def load_field_from_h5(path: str, field: str) -> np.ndarray:
    with h5py.File(path, "r") as f:
        return np.array(f[f"input/{field}"], dtype=np.float32)


def discover_bc(year: int, hour_idx: int) -> dict:
    """Auto-discover sst & tisr for one hour."""
    bc = {"sst": None, "tisr": None}
    ref = find_lucie_ref_h5(year, hour_idx)
    if ref is not None:
        try:
            with h5py.File(ref, "r") as f:
                if "sst" in f["input"]:
                    bc["sst"] = np.array(f["input/sst"], dtype=np.float32)
                if "tisr" in f["input"]:
                    bc["tisr"] = np.array(f["input/tisr"], dtype=np.float32)
        except Exception:
            pass
    if bc["sst"] is None:
        p = find_sst_h5(year, hour_idx)
        if p:
            try:
                bc["sst"] = load_field_from_h5(p, "sst")
            except Exception:
                pass
    if bc["tisr"] is None:
        p = find_tisr_h5(year, hour_idx)
        if p:
            try:
                bc["tisr"] = load_field_from_h5(p, "tisr")
            except Exception:
                pass
    return bc


def extract_year(year: int, output_dir: str, month_filter: int | None = None,
                 skip_existing: bool = False) -> tuple[int, int]:
    """Extract hourly h5 files from Troy's yearly NC for the given year.
    Optionally restrict to a single month. Returns (n_written, n_skipped)."""
    yearly_path = find_yearly_file(year)
    print(f"[y={year}] reading {yearly_path}")

    # Determine year length
    is_leap = calendar.isleap(year)
    n_hours_year = 8784 if is_leap else 8760

    # Determine which hours to extract
    if month_filter is not None:
        if not (1 <= month_filter <= 12):
            raise ValueError(f"month_filter must be 1-12, got {month_filter}")
        # Compute hour range for this month
        start_dt = datetime(year, month_filter, 1, 0, tzinfo=timezone.utc)
        if month_filter == 12:
            end_dt = datetime(year + 1, 1, 1, 0, tzinfo=timezone.utc)
        else:
            end_dt = datetime(year, month_filter + 1, 1, 0, tzinfo=timezone.utc)
        h_start = int((start_dt - datetime(year, 1, 1, tzinfo=timezone.utc)).total_seconds() // 3600)
        h_end = int((end_dt - datetime(year, 1, 1, tzinfo=timezone.utc)).total_seconds() // 3600)
        print(f"[y={year}] month={month_filter} → hours [{h_start}, {h_end})")
    else:
        h_start, h_end = 0, n_hours_year
        print(f"[y={year}] full year → hours [0, {n_hours_year})")

    statics = {}
    static_ref = find_static_ref()
    print(f"[y={year}] static fields source: {static_ref}")
    statics["orography"] = load_field_from_h5(static_ref, "orography")
    statics["land_sea_mask"] = load_field_from_h5(static_ref, "land_sea_mask")

    os.makedirs(output_dir, exist_ok=True)

    # Open the yearly NC and read just the slice we need
    print(f"[y={year}] loading regridded variables for hours [{h_start}, {h_end})")
    t0 = _time.monotonic()
    with nc.Dataset(yearly_path, "r") as ds:
        T = ds.variables["Temperature"][h_start:h_end]   # (n, 8, 48, 96)
        U = ds.variables["U-wind"][h_start:h_end]
        V = ds.variables["V-wind"][h_start:h_end]
        Q = ds.variables["Specific_Humidity"][h_start:h_end]
        logp = ds.variables["logp"][h_start:h_end]
    print(f"[y={year}] read {T.shape[0]} timesteps in {_time.monotonic()-t0:.1f}s")

    n_written = 0
    n_skipped = 0
    n_no_sst = 0
    n_no_tisr = 0
    t1 = _time.monotonic()

    for local_idx in range(T.shape[0]):
        hour_idx = h_start + local_idx
        out_filename = f"{year}_{hour_idx:04d}.h5"
        out_path = os.path.join(output_dir, out_filename)

        if skip_existing and os.path.exists(out_path):
            n_skipped += 1
            continue

        bc = discover_bc(year, hour_idx)
        if bc["sst"] is None:
            n_no_sst += 1
        if bc["tisr"] is None:
            n_no_tisr += 1

        with h5py.File(out_path, "w") as f:
            grp = f.create_group("input")
            for k in range(8):
                grp.create_dataset(f"Temperature_{k}",
                                   data=T[local_idx, k].astype(np.float32))
                grp.create_dataset(f"U-wind_{k}",
                                   data=U[local_idx, k].astype(np.float32))
                grp.create_dataset(f"V-wind_{k}",
                                   data=V[local_idx, k].astype(np.float32))
                grp.create_dataset(f"Specific_Humidity_{k}",
                                   data=Q[local_idx, k].astype(np.float32))
            grp.create_dataset("logp", data=logp[local_idx].astype(np.float32))
            grp.create_dataset("orography", data=statics["orography"])
            grp.create_dataset("land_sea_mask", data=statics["land_sea_mask"])
            if bc["sst"] is not None:
                grp.create_dataset("sst", data=bc["sst"])
            if bc["tisr"] is not None:
                grp.create_dataset("tisr", data=bc["tisr"])
            grp.create_dataset("time", data=str(hour_idx).encode())

        n_written += 1
        if n_written % 200 == 0:
            elapsed = _time.monotonic() - t1
            rate = n_written / elapsed if elapsed > 0 else 0
            print(f"[y={year}]   {n_written}/{T.shape[0]} written ({rate:.1f}/s)")

    elapsed = _time.monotonic() - t1
    print(f"[y={year}] done: {n_written} written, {n_skipped} skipped in {elapsed:.1f}s")
    if n_no_sst:
        print(f"[y={year}]   warning: {n_no_sst} files with no SST source")
    if n_no_tisr:
        print(f"[y={year}]   warning: {n_no_tisr} files with no TISR source")
    return n_written, n_skipped


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    grp = p.add_mutually_exclusive_group(required=True)
    grp.add_argument("--year", type=int, help="Single year to process")
    grp.add_argument("--year-range", type=int, nargs=2, metavar=("Y1", "Y2"),
                     help="Inclusive year range")
    p.add_argument("--month", type=int, default=None,
                   help="Optional: restrict to one month (only with --year)")
    p.add_argument("--output-dir", default=os.path.join(WORKDIR, "h5_out"))
    p.add_argument("--skip-existing", action="store_true")
    args = p.parse_args()

    if args.year_range and args.month:
        sys.exit("--month can only be used with --year, not --year-range")

    years = [args.year] if args.year else list(range(args.year_range[0], args.year_range[1] + 1))

    total_written = 0
    total_skipped = 0
    for y in years:
        try:
            w, s = extract_year(y, args.output_dir, args.month, args.skip_existing)
            total_written += w
            total_skipped += s
        except FileNotFoundError as e:
            print(f"[y={y}] ERROR: {e}")
            continue

    print()
    print(f"=== ALL DONE: {total_written} written, {total_skipped} skipped ===")


if __name__ == "__main__":
    main()
