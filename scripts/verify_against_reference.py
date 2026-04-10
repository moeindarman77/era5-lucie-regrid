"""
Compare a generated LUCIE-format h5 file against a reference h5 file.
Reports per-variable RMSE, MaxAbsErr, ranges, and a PASS/FAIL flag.

Usage:
  python verify_against_reference.py \
      --generated h5_out/2013_0000.h5 \
      --reference /glade/derecho/scratch/mdarman/ERA5_hr_haiwen/LUCIE_3D/1_step_1hr_h5df_test/test/2013_0000.h5
"""
from __future__ import annotations

import argparse
import os

import h5py
import numpy as np


def verify(generated_path: str, reference_path: str) -> bool:
    print(f"\n{'='*78}")
    print(f"VERIFICATION: {os.path.basename(generated_path)} vs {os.path.basename(reference_path)}")
    print(f"  generated: {generated_path}")
    print(f"  reference: {reference_path}")
    print(f"{'='*78}")

    with h5py.File(generated_path, "r") as fg, h5py.File(reference_path, "r") as fr:
        gen_vars = set(fg["input"].keys())
        ref_vars = set(fr["input"].keys())

        common = sorted(gen_vars & ref_vars)
        only_gen = sorted(gen_vars - ref_vars)
        only_ref = sorted(ref_vars - gen_vars)

        if only_gen:
            print(f"\n  Only in generated: {only_gen}")
        if only_ref:
            print(f"  Only in reference (not checked): {len(only_ref)} vars: {only_ref[:10]}{'...' if len(only_ref) > 10 else ''}")

        print(f"\n  {'Variable':<25} {'RMSE':>14} {'MaxAbsErr':>14} {'GenRange':>26} {'RefRange':>26} {'Status'}")
        print(f"  {'-'*114}")

        all_pass = True
        for var in common:
            g = fg["input"][var]
            r = fr["input"][var]

            if g.shape == () or r.shape == ():
                continue
            try:
                gd = np.array(g, dtype=np.float64)
                rd = np.array(r, dtype=np.float64)
            except Exception:
                continue

            if gd.shape != rd.shape:
                print(f"  {var:<25} SHAPE MISMATCH: {gd.shape} vs {rd.shape}")
                all_pass = False
                continue

            diff = gd - rd
            rmse = float(np.sqrt(np.mean(diff ** 2)))
            mxe = float(np.abs(diff).max())
            gr = f"[{gd.min():.4g}, {gd.max():.4g}]"
            rr = f"[{rd.min():.4g}, {rd.max():.4g}]"

            mag = max(abs(rd.max()), abs(rd.min()), 1e-10)
            threshold = max(1e-5, 1e-6 * mag)
            status = "PASS" if rmse < threshold else "FAIL"
            if status == "FAIL":
                all_pass = False
            print(f"  {var:<25} {rmse:>14.6g} {mxe:>14.6g} {gr:>26} {rr:>26} {status}")

        print(f"\n  {'='*78}")
        if all_pass:
            print(f"  ALL VARIABLES PASS")
        else:
            print(f"  SOME VARIABLES DEVIATE FROM REFERENCE — see RMSE column above")
        print(f"  {'='*78}")

    return all_pass


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--generated", required=True)
    p.add_argument("--reference", required=True)
    args = p.parse_args()
    verify(args.generated, args.reference)


if __name__ == "__main__":
    main()
