#!/usr/bin/env python3
"""Legacy myosin Phi_3D for the s1s2 sweep, compared to s2only at MATCHED time.

Why a matched window is mandatory here
--------------------------------------
The s2only sweep is complete (~300k steps, ~1500 frames). The s1s2 sweep is
still running and much shorter. Comparing "each run's own last 300 frames"
would compare s1s2's transient against s2only's steady state and manufacture a
difference that is purely an age gap -- the same trap that made the raw
status-1 counts look divergent at high lambda.

So for each (lambda, replica) we take the number of frames the S1S2 run has
actually committed, and evaluate BOTH sweeps over the SAME absolute frame
window. Both sweeps use save_every=200 and identical nsteps, so frame index
maps to simulation step identically across them -- matched frames == matched
steps.

Phi_3D is the legacy myosin-only definition (max over 3D Fourier modes), the
same one behind legacy_status2_directional_summary.csv, so numbers here are
directly comparable to that file.

Reported per (lambda, replica):
  phi_matched_s1s2 / phi_matched_s2only  -- like-for-like, the honest comparison
  phi_own_late_s2only                    -- s2only at its own steady state, for context
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
S1S2_DIR = SCRATCH / "data_a2000_s1s2_lambda"
S2ONLY_DIR = SCRATCH / "data_a2000_status2_lambda"

LABELS = ("0", "1e-5.5", "1e-5", "1e-4.5", "1e-4", "1e-3.5", "1e-3",
          "1e-2.5", "1e-2", "1e-1.5", "1e-1", "1e-0.5", "1e0")

FIELDS = ("lt", "lambda", "rep", "n_frames_s1s2", "n_frames_s2only",
          "win_lo", "win_hi", "step_lo", "step_hi", "n_sampled",
          "phi_matched_s1s2", "phi_matched_s2only", "phi_own_late_s2only",
          "nx_s1s2", "nx_s2only")


def lt_to_float(lt):
    lt = str(lt).strip()
    if lt in ("0", "lt0", ""):
        return 0.0
    mant, exp = lt.split("e")
    return float(mant) * 10.0 ** float(exp)


def committed(handle, key="/myosin/center"):
    if "/meta/committed_frames" in handle:
        c = int(np.asarray(handle["/meta/committed_frames"]).reshape(-1)[0])
        if c > 0:
            return min(c, handle[key].shape[0])
    return handle[key].shape[0]


def phi_over(path, lo, hi, stride, device, batch=16):
    """Legacy myosin Phi_3D averaged over frames [lo,hi) with the given stride."""
    idx = list(range(lo, hi, stride))
    if not idx:
        return None, None, 0
    with h5py.File(path, "r", locking=False) as handle:
        if "/myosin/center" not in handle:
            return None, None, 0
        centers = handle["/myosin/center"]
        idx = [i for i in idx if i < centers.shape[0]]
        if not idx:
            return None, None, 0
        results = []
        for s in range(0, len(idx), batch):
            coords = np.asarray(centers[idx[s:s + batch]], dtype=np.float64)
            results.extend(fourier_order_3d_frames(coords, device))
    phi = float(np.mean([r["phi_3d_max"] for r in results]))
    modes = [(r["nx"], r["ny"], r["nz"]) for r in results]
    nx = Counter(modes).most_common(1)[0][0][0]
    return phi, nx, len(results)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--output", type=Path,
                   default=Path("phi3d_s1s2_vs_s2only_matched.csv"))
    p.add_argument("--window", type=int, default=300,
                   help="frames back from the s1s2 frontier to average over")
    p.add_argument("--stride", type=int, default=10)
    p.add_argument("--min-frames", type=int, default=120,
                   help="skip a lambda whose s1s2 run is younger than this")
    p.add_argument("--save-every", type=int, default=200)
    p.add_argument("--device", choices=("auto", "cuda", "cpu", "numpy"),
                   default="auto")
    args = p.parse_args()

    device = resolve_compute_device(args.device)
    if device is None:
        print("Using numpy for legacy Phi_3D.", flush=True)
    else:
        print_compute_device_info(device)

    rows, skipped = [], []
    for lt in LABELS:
        for rep in range(2):
            stem = f"_lt{lt}_r{rep:02d}_L12x2x2.h5"
            f1 = S1S2_DIR / f"A2000_M1000_s1s2{stem}"
            f2 = S2ONLY_DIR / f"A2000_M1000_s2only{stem}"
            if not f1.exists() or not f2.exists():
                skipped.append((lt, rep, "missing file"))
                continue
            try:
                with h5py.File(f1, "r", locking=False) as h:
                    n1 = committed(h)
                with h5py.File(f2, "r", locking=False) as h:
                    n2 = committed(h)
            except Exception as exc:  # noqa: BLE001
                skipped.append((lt, rep, f"{type(exc).__name__}"))
                continue

            if n1 < args.min_frames:
                skipped.append((lt, rep, f"s1s2 only {n1} frames"))
                continue

            # Evaluate both sweeps over the SAME absolute frames, anchored on
            # the s1s2 frontier (the shorter run sets the pace).
            hi = min(n1, n2)
            lo = max(0, hi - args.window)
            pm1, nx1, k1 = phi_over(f1, lo, hi, args.stride, device)
            pm2, nx2, _k2 = phi_over(f2, lo, hi, args.stride, device)
            # s2only at its own late time, for context only.
            lo2 = max(0, n2 - args.window)
            po2, _nx, _k = phi_over(f2, lo2, n2, args.stride, device)
            if pm1 is None or pm2 is None:
                skipped.append((lt, rep, "no frames in window"))
                continue

            rows.append({
                "lt": lt, "lambda": lt_to_float(lt), "rep": f"r{rep:02d}",
                "n_frames_s1s2": n1, "n_frames_s2only": n2,
                "win_lo": lo, "win_hi": hi,
                "step_lo": lo * args.save_every, "step_hi": hi * args.save_every,
                "n_sampled": k1,
                "phi_matched_s1s2": pm1, "phi_matched_s2only": pm2,
                "phi_own_late_s2only": po2,
                "nx_s1s2": nx1, "nx_s2only": nx2,
            })
            print(f"lt={lt:>7} r{rep:02d} frames s1s2={n1:4d}/s2only={n2:4d} "
                  f"win=[{lo},{hi}) steps=[{lo*args.save_every},{hi*args.save_every}) "
                  f"| phi s1s2={pm1:.3f} s2only={pm2:.3f} "
                  f"(s2only own-late={po2:.3f})", flush=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as s:
        w = csv.DictWriter(s, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {len(rows)} rows to {args.output}")

    # ---- aggregate over replicas ----
    by_lt = defaultdict(list)
    for r in rows:
        by_lt[r["lt"]].append(r)
    print("\n=== legacy Phi_3D, MATCHED frame window (s1s2 sets the window) ===")
    print(f"{'lambda':>9} {'n':>2} {'steps':>16} | {'s1s2':>7} {'s2only':>7} "
          f"{'diff':>7} | {'s2only own-late':>15}")
    print("-" * 76)
    for lt in sorted(by_lt, key=lt_to_float):
        rs = by_lt[lt]
        a = float(np.mean([r["phi_matched_s1s2"] for r in rs]))
        b = float(np.mean([r["phi_matched_s2only"] for r in rs]))
        c = float(np.mean([r["phi_own_late_s2only"] for r in rs]))
        lo = min(r["step_lo"] for r in rs)
        hi = max(r["step_hi"] for r in rs)
        print(f"{lt_to_float(lt):>9.3g} {len(rs):>2} {f'{lo}-{hi}':>16} | "
              f"{a:>7.3f} {b:>7.3f} {a-b:>+7.3f} | {c:>15.3f}")

    if skipped:
        print(f"\n{len(skipped)} (lambda, rep) skipped:")
        for lt, rep, why in skipped:
            print(f"  lt={lt} r{rep:02d}: {why}")


if __name__ == "__main__":
    main()
