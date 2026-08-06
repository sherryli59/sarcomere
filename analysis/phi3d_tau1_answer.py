#!/usr/bin/env python3
"""Does a constant, load-INDEPENDENT tau1 recover the s1s2 phi_3d boost?

The old ("s1s2") kinetics stabilized status-1 bonds by lambda*f_load -- a
LOAD-DEPENDENT lifetime the user flagged as unphysical. At lambda=1e-3 that
bought +0.33 in legacy myosin phi_3d over s2only (status-1 unstabilized).

This sweep replaces that with a constant tau1 (base_lifetime), load-independent,
catch stabilization on status 2 only. For each (tau1, lambda) we ask where its
phi_3d lands between the s2only FLOOR (tau1->0) and the s1s2 CEILING (old
load-dependent). If some finite tau1 reaches the ceiling, the boost was just
"bonds live longer" (physical, keepable). If none do, the boost required the
load-dependence -- and s2only is the right model.

phi_3d = legacy myosin-only definition (max over 3D Fourier modes), identical to
compare_phi3d_s1s2_vs_s2only.py, so the s1s2/s2only reference numbers are
directly comparable.

For each completed (>=min_step) run we average phi_3d over its own last `window`
frames (its steady state). Seeds are 97000+rep across all three sweeps, so the
tau1 / s1s2 / s2only cells at a given (lambda, rep) share initial conditions.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import h5py
import numpy as np

from extract_structure_3d import (
    fourier_order_3d_frames,
    print_compute_device_info,
    resolve_compute_device,
)

SCRATCH = Path("/scratch/users/xsli/sarcomere")
TAU1_DIR = SCRATCH / "data_a2000_tau1_sweep"
S1S2_DIR = SCRATCH / "data_a2000_s1s2_lambda"
S2ONLY_DIR = SCRATCH / "data_a2000_status2_lambda"

# tau1 label -> float, matching the sbatch generator
TAU_LABELS = {
    "0": 0.0, "1e-4.5": 10 ** -4.5, "1e-4": 1e-4,
    "1e-3.5": 10 ** -3.5, "1e-3": 1e-3, "1e-2.5": 10 ** -2.5,
    # lambda=1e-2 extension rungs (a2000_tau1_ext_lt1e-2.sbatch)
    "1e-2": 1e-2, "1e-1.5": 10 ** -1.5,
}
LT_LABELS = {"1e-3": 1e-3, "1e-2": 1e-2}


def step_of(h):
    for k in ("/resume/current_step", "/resume_prev/current_step"):
        if k in h:
            return int(np.asarray(h[k]).reshape(-1)[-1])
    return int(np.asarray(h["/state/current_step"])[-1, 0])


def phi_late(path, window, stride, device, min_step, batch=16):
    """Legacy myosin phi_3d over the last `window` frames; None if too young."""
    if not path.exists():
        return None
    with h5py.File(path, "r", locking=False) as h:
        if step_of(h) < min_step or "/myosin/center" not in h:
            return None
        centers = h["/myosin/center"]
        n = centers.shape[0]
        idx = list(range(max(0, n - window), n, stride))
        if not idx:
            return None
        res = []
        for s in range(0, len(idx), batch):
            coords = np.asarray(centers[idx[s:s + batch]], dtype=np.float64)
            res.extend(fourier_order_3d_frames(coords, device))
    return float(np.mean([r["phi_3d_max"] for r in res]))


def mean_late(directory, stem, window, stride, device, min_step):
    """Average over whatever replicas of this stem are complete."""
    vals = []
    for rep in range(2):
        p = directory / stem.format(rep=rep)
        v = phi_late(p, window, stride, device, min_step)
        if v is not None:
            vals.append(v)
    return (float(np.mean(vals)), len(vals)) if vals else (None, 0)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--window", type=int, default=300)
    p.add_argument("--stride", type=int, default=10)
    p.add_argument("--min-step", type=int, default=290000)
    p.add_argument("--device", choices=("auto", "cuda", "cpu", "numpy"),
                   default="auto")
    p.add_argument("--output", type=Path,
                   default=Path("phi3d_tau1_answer.csv"))
    args = p.parse_args()

    device = resolve_compute_device(args.device)
    if device is None:
        print("Using numpy.", flush=True)
    else:
        print_compute_device_info(device)

    # Reference floors/ceilings per lambda (own-late, completed replicas).
    ref = {}
    for lt in LT_LABELS:
        s1, n1 = mean_late(S1S2_DIR, f"A2000_M1000_s1s2_lt{lt}_r{{rep:02d}}_L12x2x2.h5",
                           args.window, args.stride, device, args.min_step)
        s2, n2 = mean_late(S2ONLY_DIR, f"A2000_M1000_s2only_lt{lt}_r{{rep:02d}}_L12x2x2.h5",
                           args.window, args.stride, device, args.min_step)
        ref[lt] = (s1, n1, s2, n2)
        print(f"[ref] lt={lt}: s1s2={s1} (n={n1})  s2only={s2} (n={n2})", flush=True)

    rows = []
    for lt in LT_LABELS:
        s1, _n1, s2, _n2 = ref[lt]
        for tau in TAU_LABELS:
            stem = f"A2000_M1000_s2catch_tau{tau}_lt{lt}_r{{rep:02d}}_L12x2x2.h5"
            phi, n = mean_late(TAU1_DIR, stem, args.window, args.stride,
                               device, args.min_step)
            if phi is None:
                print(f"tau={tau:>7} lt={lt}: no complete replica", flush=True)
                continue
            # where does tau1 land in the [s2only, s1s2] span? 0=floor, 1=ceiling
            frac = None
            if s1 is not None and s2 is not None and abs(s1 - s2) > 1e-9:
                frac = (phi - s2) / (s1 - s2)
            rows.append({
                "lambda": LT_LABELS[lt], "lt": lt, "tau1": TAU_LABELS[tau],
                "tau1_label": tau, "n_rep": n, "phi_tau1": phi,
                "phi_s2only_floor": s2, "phi_s1s2_ceiling": s1,
                "recovered_frac": frac,
            })
            fs = f"{frac:+.2f}" if frac is not None else "  n/a"
            print(f"tau={tau:>7} lt={lt}: phi={phi:.3f} (n={n})  "
                  f"floor={s2:.3f} ceil={s1:.3f}  recovered={fs}", flush=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else
                           ["lambda"])
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {len(rows)} rows to {args.output}")

    print("\n=== phi_3d vs constant tau1  (floor=s2only, ceiling=s1s2) ===")
    for lt in LT_LABELS:
        s1, _a, s2, _b = ref[lt]
        print(f"\nlambda={LT_LABELS[lt]:g}   floor(s2only)={s2:.3f}  "
              f"ceiling(s1s2)={s1:.3f}  gap={s1-s2:+.3f}")
        print(f"  {'tau1':>9} {'phi':>7} {'recovered_frac':>15}")
        for r in sorted([r for r in rows if r["lt"] == lt],
                        key=lambda r: r["tau1"]):
            fr = f"{r['recovered_frac']:+.2f}" if r["recovered_frac"] is not None else "n/a"
            print(f"  {r['tau1']:>9.3g} {r['phi_tau1']:>7.3f} {fr:>15}")


if __name__ == "__main__":
    main()
