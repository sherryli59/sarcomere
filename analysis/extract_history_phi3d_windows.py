#!/usr/bin/env python3
"""Extract full-3D order in selected step windows of one history trajectory."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import h5py
import numpy as np

from extract_sarcomeric_order import (
    build_3d_modes,
    fourier_order_3d,
    shared_sarcomeric_order,
)


FIELDS = (
    "filename", "window_start", "window_end", "frame", "step",
    "status1_fraction", "status2_fraction", "phi_3d_myosin",
    "phi_3d_actin_cb2", "shared_phi_3d_myosin_corrected",
    "shared_phi_3d_actin_cb2_corrected", "shared_joint_band_amplitude",
    "shared_nx_3d", "shared_ny_3d", "shared_nz_3d", "q_sarcomere",
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--window", action="append", required=True,
        help="Inclusive simulation-step window START:END.",
    )
    parser.add_argument("--samples", type=int, default=3)
    args = parser.parse_args()

    windows = [tuple(map(int, item.split(":"))) for item in args.window]
    k_matrix, mode_matrix = build_3d_modes(12.0, 2.0, 2.0, 30, 5)
    rows = []
    with h5py.File(args.input, "r", locking=False) as handle:
        n_frames = min(
            handle["/myosin/center"].shape[0],
            handle["/actin/center"].shape[0],
            handle["/actin/cb_status"].shape[0],
            handle["/state/current_step"].shape[0],
        )
        steps = np.asarray(handle["/state/current_step"][:n_frames, 0], dtype=int)
        for start, end in windows:
            candidates = np.flatnonzero((steps >= start) & (steps <= end))
            if not len(candidates):
                continue
            select = candidates[
                np.unique(np.linspace(0, len(candidates) - 1, args.samples, dtype=int))
            ]
            for frame in select:
                myosin = np.asarray(handle["/myosin/center"][frame], dtype=float)
                actin = np.asarray(handle["/actin/center"][frame], dtype=float)
                status = np.asarray(handle["/actin/cb_status"][frame]).reshape(-1)
                actin_cb2 = actin[status == 2]
                phi_m = fourier_order_3d(myosin, k_matrix, mode_matrix)
                phi_z = fourier_order_3d(actin_cb2, k_matrix, mode_matrix)
                shared = shared_sarcomeric_order(
                    myosin, actin_cb2, k_matrix, mode_matrix, 2.0, 2.0, 2, 2
                )
                rows.append({
                    "filename": args.input.name,
                    "window_start": start,
                    "window_end": end,
                    "frame": int(frame),
                    "step": int(steps[frame]),
                    "status1_fraction": float(np.mean(status == 1)),
                    "status2_fraction": float(np.mean(status == 2)),
                    "phi_3d_myosin": phi_m["phi"],
                    "phi_3d_actin_cb2": phi_z["phi"],
                    "shared_phi_3d_myosin_corrected": shared["phi_m_corrected"],
                    "shared_phi_3d_actin_cb2_corrected": shared["phi_z_corrected"],
                    "shared_joint_band_amplitude": shared["joint"],
                    "shared_nx_3d": shared["nx"],
                    "shared_ny_3d": shared["ny"],
                    "shared_nz_3d": shared["nz"],
                    "q_sarcomere": shared["q_sarcomere"],
                })
                print(f"{args.input.name}: step {steps[frame]}", flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} rows to {args.output}", flush=True)


if __name__ == "__main__":
    main()
