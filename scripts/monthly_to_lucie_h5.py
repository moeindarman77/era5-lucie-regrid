"""
Convert a monthly Fortran-regridded NetCDF (the output of run_month.sh) into
per-hour LUCIE-format HDF5 files matching the training data format.

For each of the 744 (or 696, or 720) timesteps in a monthly file, write:
  h5_out/{YYYY}_{HHHH:04d}.h5

where HHHH is the hour-of-year index (0..8759 or 0..8783 for leap years).

Each h5 file has group `input/` containing:
  Temperature_0..7    (48, 96) float32  ← from Fortran regridder
  U-wind_0..7         (48, 96) float32  ← from Fortran regridder
  V-wind_0..7         (48, 96) float32  ← from Fortran regridder
  Specific_Humidity_0..7 (48, 96) float32 ← from Fortran regridder
  logp                (48, 96) float32  ← from Fortran regridder
  orography           (48, 96) float32  ← from STATIC_REF_H5
  land_sea_mask       (48, 96) float32  ← from STATIC_REF_H5
  sst                 (48, 96) float32  ← auto-discovered (see below)
  tisr                (48, 96) float32  ← auto-discovered (see below)
  tp6hr               (48, 96) float32  ← optional, from --tp-source if given
  time                bytes scalar      ← string of HHHH

Auto-discovery for sst/tisr (in priority order):
  1. LUCIE reference (train/val/test) per-hour h5 file at the same hour-of-year
     — most reliable for years 1981-2014.
  2. /glade/derecho/scratch/mdarman/sst_data/era5_data_download/sst_2014-2023/{YYYY}_{HHHH}.h5
  3. /glade/derecho/scratch/mdarman/sst_data/era5_data_download/tisr_2015-2023/{YYYY}_{HHHH}.h5
  4. Skip variable (warning) if no source available.

NOTE: This script ONLY processes the 5 prognostic variables produced by the
regridder. The full LUCIE training format also has q_con, ocean_ptemp_*, ohc
variables — those are added by separate post-processing pipelines and are NOT
included here. If you need them, copy them from the reference (for years
where reference exists).

Lat orientation: Fortran's regrid_var ends up in S→N order (matches LUCIE
training data convention) due to two flips canceling. NO flip needed here.

Sigma level convention: Fortran writes Sigma_Level dim with index 0 = σ=0.025
(top of atmosphere) and index 7 = σ=0.95 (near surface). LUCIE training data
uses the same convention. So Temperature_0 in our h5 = T at σ=0.025.

Usage:
  python monthly_to_lucie_h5.py --year 2014 --month 1
  python monthly_to_lucie_h5.py --year 2014 --month 1 --output-dir h5_out_test
  python monthly_to_lucie_h5.py --year 2014 --month 1 --skip-existing
"""
from __future__ import annotations

import argparse
import calendar
import os
import sys
from datetime import datetime, timezone

import numpy as np
import h5py
import netCDF4 as nc

WORKDIR = os.environ.get("LUCIE_WORKDIR",
    "/glade/derecho/scratch/mdarman/sst_data/era5_data_download/gen_lucie_3d")

# Sources for boundary-condition fields — change these for your setup.
# See RECIPE.md for how to obtain these files.
LUCIE_REF_BASE = os.environ.get("LUCIE_REF_BASE",
    "/glade/derecho/scratch/mdarman/ERA5_hr_haiwen/LUCIE_3D/1_step_1hr_h5df_test")
SST_DIR = os.environ.get("LUCIE_SST_DIR",
    "/glade/derecho/scratch/mdarman/sst_data/era5_data_download/sst_2014-2023")
TISR_DIR = os.environ.get("LUCIE_TISR_DIR",
    "/glade/derecho/scratch/mdarman/sst_data/era5_data_download/tisr_2015-2023")

# A reference h5 file for static fields (orography, land_sea_mask) — these are constant
# in the LUCIE training data. Pick any reference file that exists.
STATIC_REF_H5_CANDIDATES = [
    os.path.join(LUCIE_REF_BASE, "test", "2013_0000.h5"),
    os.path.join(LUCIE_REF_BASE, "val", "2010_0000.h5"),
    os.path.join(LUCIE_REF_BASE, "train", "1990_0000.h5"),
]


def find_static_ref() -> str:
    for p in STATIC_REF_H5_CANDIDATES:
        if os.path.exists(p):
            return p
    raise FileNotFoundError(
        f"None of the static reference h5 files exist: {STATIC_REF_H5_CANDIDATES}"
    )


def hour_of_year(year: int, month: int, day: int, hour: int) -> int:
    """Convert a UTC datetime to the hour-of-year index used in LUCIE filenames."""
    start = datetime(year, 1, 1, tzinfo=timezone.utc)
    target = datetime(year, month, day, hour, tzinfo=timezone.utc)
    return int((target - start).total_seconds() // 3600)


def month_start_hour(year: int, month: int) -> int:
    return hour_of_year(year, month, 1, 0)


def find_lucie_ref_h5(year: int, hour_idx: int) -> str | None:
    """Find a LUCIE reference h5 file for the given (year, hour) if it exists."""
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


def load_static_fields() -> dict[str, np.ndarray]:
    """Load orography + land_sea_mask from a reference h5 — these never change."""
    ref = find_static_ref()
    print(f"  static fields source: {ref}")
    return {
        "orography": load_field_from_h5(ref, "orography"),
        "land_sea_mask": load_field_from_h5(ref, "land_sea_mask"),
    }


def discover_bc_for_hour(year: int, hour_idx: int) -> dict[str, np.ndarray | None]:
    """Auto-discover sst and tisr for one (year, hour_idx). Returns None for any
    field that couldn't be sourced."""
    bc: dict[str, np.ndarray | None] = {"sst": None, "tisr": None}

    # Priority 1: LUCIE reference h5 (has both sst & tisr)
    ref = find_lucie_ref_h5(year, hour_idx)
    if ref is not None:
        with h5py.File(ref, "r") as f:
            keys = set(f["input"].keys())
        try:
            if "sst" in keys:
                bc["sst"] = load_field_from_h5(ref, "sst")
            if "tisr" in keys:
                bc["tisr"] = load_field_from_h5(ref, "tisr")
        except Exception:
            pass

    # Priority 2: dedicated sst dir
    if bc["sst"] is None:
        p = find_sst_h5(year, hour_idx)
        if p is not None:
            try:
                bc["sst"] = load_field_from_h5(p, "sst")
            except Exception:
                pass

    # Priority 3: dedicated tisr dir
    if bc["tisr"] is None:
        p = find_tisr_h5(year, hour_idx)
        if p is not None:
            try:
                bc["tisr"] = load_field_from_h5(p, "tisr")
            except Exception:
                pass

    return bc


def convert_month(year: int, month: int, output_dir: str,
                  skip_existing: bool = False, max_warnings: int = 5) -> None:
    fort_nc = os.path.join(WORKDIR, "outputs", str(year),
                           f"era_5_m{month}_y{year}_regridded_spectral_mpi.nc")
    if not os.path.exists(fort_nc):
        raise FileNotFoundError(f"Fortran NC not found: {fort_nc}")

    print(f"=== Reading {fort_nc} ===")
    with nc.Dataset(fort_nc, "r") as ds:
        T = ds.variables["Temperature"][:]   # (ntimes, 8, 48, 96)
        U = ds.variables["U-wind"][:]
        V = ds.variables["V-wind"][:]
        Q = ds.variables["Specific_Humidity"][:]
        logp = ds.variables["logp"][:]       # (ntimes, 48, 96)
        ntimes = T.shape[0]
    print(f"  ntimes={ntimes}")

    # Sanity-check ntimes against the expected number for the month
    expected_ntimes = calendar.monthrange(year, month)[1] * 24
    if ntimes != expected_ntimes:
        print(f"  WARNING: expected {expected_ntimes} timesteps, got {ntimes}")

    statics = load_static_fields()
    orography = statics["orography"]
    land_sea_mask = statics["land_sea_mask"]

    base_hour = month_start_hour(year, month)
    print(f"  Month start hour-of-year: {base_hour}")

    os.makedirs(output_dir, exist_ok=True)

    n_written = 0
    n_skipped = 0
    n_no_sst = 0
    n_no_tisr = 0
    warned = 0

    for t_idx in range(ntimes):
        hour_idx = base_hour + t_idx
        out_filename = f"{year}_{hour_idx:04d}.h5"
        out_path = os.path.join(output_dir, out_filename)

        if skip_existing and os.path.exists(out_path):
            n_skipped += 1
            continue

        bc = discover_bc_for_hour(year, hour_idx)
        if bc["sst"] is None:
            n_no_sst += 1
            if warned < max_warnings:
                print(f"  WARNING: no SST source for hour {hour_idx}")
                warned += 1
        if bc["tisr"] is None:
            n_no_tisr += 1

        with h5py.File(out_path, "w") as f:
            grp = f.create_group("input")

            # Prognostic atmospheric variables (5 of them, 8 sigma levels each)
            for k in range(8):
                grp.create_dataset(f"Temperature_{k}",
                                   data=T[t_idx, k].astype(np.float32))
                grp.create_dataset(f"U-wind_{k}",
                                   data=U[t_idx, k].astype(np.float32))
                grp.create_dataset(f"V-wind_{k}",
                                   data=V[t_idx, k].astype(np.float32))
                grp.create_dataset(f"Specific_Humidity_{k}",
                                   data=Q[t_idx, k].astype(np.float32))

            grp.create_dataset("logp", data=logp[t_idx].astype(np.float32))

            # Static fields
            grp.create_dataset("orography", data=orography)
            grp.create_dataset("land_sea_mask", data=land_sea_mask)

            # Boundary conditions (only write if available)
            if bc["sst"] is not None:
                grp.create_dataset("sst", data=bc["sst"])
            if bc["tisr"] is not None:
                grp.create_dataset("tisr", data=bc["tisr"])

            # Time
            grp.create_dataset("time", data=str(hour_idx).encode())

        n_written += 1
        if n_written % 100 == 0:
            print(f"  ...{n_written}/{ntimes} written")

    print()
    print(f"=== Done: {n_written} written, {n_skipped} skipped ===")
    if n_no_sst:
        print(f"  WARNING: {n_no_sst} files missing SST")
    if n_no_tisr:
        print(f"  WARNING: {n_no_tisr} files missing TISR")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--year", type=int, required=True)
    p.add_argument("--month", type=int, required=True)
    p.add_argument("--output-dir", default=os.path.join(WORKDIR, "h5_out"),
                   help="Where to write {YYYY}_{HHHH}.h5 files (default: %(default)s)")
    p.add_argument("--skip-existing", action="store_true",
                   help="Skip writing files that already exist.")
    args = p.parse_args()

    convert_month(args.year, args.month, args.output_dir,
                  skip_existing=args.skip_existing)


if __name__ == "__main__":
    main()
