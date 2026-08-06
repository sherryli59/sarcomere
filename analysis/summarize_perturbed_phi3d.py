#!/usr/bin/env python3
"""Summarize late-time full-3D order for the catch-bond perturbation runs."""

from __future__ import annotations

import argparse
import csv
import os
import re
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import h5py
import numpy as np

from extract_sarcomeric_order import (
    build_3d_modes,
    fourier_order_3d,
    shared_sarcomeric_order,
)


SCRATCH = Path("/scratch/users/xsli/sarcomere")
CAMPAIGNS = {
    "status2_directional": (SCRATCH / "data_a2000_status2_lambda", 300_000),
    "status2_nondirectional": (
        SCRATCH / "data_a2000_status2_nondir_lambda",
        300_000,
    ),
    # Original (pre-split) kinetics on the current code: lambda*f_load
    # stabilizes status 1 and status 2 alike. The control that separates the
    # kinetics change from the -O0 -> -O2 / refactor speedups.
    "status1_and_status2": (SCRATCH / "data_a2000_s1s2_lambda", 300_000),
    # Constant, load-independent tau1 with catch stabilization on status 2 only.
    "constant_tau1": (SCRATCH / "data_a2000_tau1_sweep", 300_000),
    "lifetime_match": (SCRATCH / "data_status2_lifetime_match", 150_000),
    "constant_bracket_a1500": (SCRATCH / "data_intervention", 300_000),
    "history_current": (SCRATCH / "data_history_status2", 400_000),
}

FIELDS = (
    "campaign", "filename", "frame", "step", "checkpoint_step", "complete",
    "status1_fraction", "status2_fraction", "n_myosin", "n_actin_cb2",
    "phi_3d_myosin", "nx_3d_myosin", "ny_3d_myosin", "nz_3d_myosin",
    "phi_3d_actin_cb2", "nx_3d_actin_cb2", "ny_3d_actin_cb2",
    "nz_3d_actin_cb2", "shared_nx_3d", "shared_ny_3d", "shared_nz_3d",
    "shared_phi_3d_myosin_corrected", "shared_phi_3d_actin_cb2_corrected",
    "shared_joint_band_amplitude", "shared_registration_pi", "q_sarcomere",
)


def checkpoint_step(handle: h5py.File) -> int:
    for key in ("/resume/current_step", "/resume_prev/current_step"):
        if key in handle:
            return int(np.asarray(handle[key]).reshape(-1)[-1])
    return int(np.asarray(handle["/state/current_step"])[-1, 0])


def selected_frames(n_frames: int, samples: int, window: int) -> np.ndarray:
    start = max(0, n_frames - 1 - window)
    return np.unique(np.linspace(start, n_frames - 1, min(samples, n_frames), dtype=int))


def analyze(path: Path, campaign: str, target: int, args, k_matrix, mode_matrix):
    records = []
    with h5py.File(path, "r", locking=False) as handle:
        current = checkpoint_step(handle)
        complete = current >= target
        if not complete and not args.include_incomplete:
            return records
        required = ("/myosin/center", "/actin/center", "/actin/cb_status")
        if not all(key in handle for key in required):
            return records
        n_frames = min(handle[key].shape[0] for key in required)
        frames = selected_frames(n_frames, args.late_samples, args.late_window_frames)
        state_steps = handle.get("/state/current_step")
        for frame in frames:
            myosin = np.asarray(handle["/myosin/center"][frame], dtype=float)
            actin = np.asarray(handle["/actin/center"][frame], dtype=float)
            status = np.asarray(handle["/actin/cb_status"][frame]).reshape(-1)
            actin_cb2 = actin[status == 2]
            phi_m = fourier_order_3d(myosin, k_matrix, mode_matrix)
            phi_z = fourier_order_3d(actin_cb2, k_matrix, mode_matrix)
            shared = shared_sarcomeric_order(
                myosin, actin_cb2, k_matrix, mode_matrix,
                args.ly, args.lz, args.yz_bins, args.yz_bins,
            )
            step = (
                int(np.asarray(state_steps[frame]).reshape(-1)[0])
                if state_steps is not None and frame < state_steps.shape[0]
                else frame * args.save_every
            )
            records.append({
                "campaign": campaign,
                "filename": path.name,
                "frame": int(frame),
                "step": step,
                "checkpoint_step": current,
                "complete": int(complete),
                "status1_fraction": float(np.mean(status == 1)),
                "status2_fraction": float(np.mean(status == 2)),
                "n_myosin": len(myosin),
                "n_actin_cb2": len(actin_cb2),
                "phi_3d_myosin": phi_m["phi"],
                "nx_3d_myosin": phi_m["nx"],
                "ny_3d_myosin": phi_m["ny"],
                "nz_3d_myosin": phi_m["nz"],
                "phi_3d_actin_cb2": phi_z["phi"],
                "nx_3d_actin_cb2": phi_z["nx"],
                "ny_3d_actin_cb2": phi_z["ny"],
                "nz_3d_actin_cb2": phi_z["nz"],
                "shared_nx_3d": shared["nx"],
                "shared_ny_3d": shared["ny"],
                "shared_nz_3d": shared["nz"],
                "shared_phi_3d_myosin_corrected": shared["phi_m_corrected"],
                "shared_phi_3d_actin_cb2_corrected": shared["phi_z_corrected"],
                "shared_joint_band_amplitude": shared["joint"],
                "shared_registration_pi": shared["registration_pi"],
                "q_sarcomere": shared["q_sarcomere"],
            })
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--campaign", action="append", choices=CAMPAIGNS)
    parser.add_argument("--include-incomplete", action="store_true")
    parser.add_argument("--late-samples", type=int, default=3)
    parser.add_argument("--late-window-frames", type=int, default=50)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--lx", type=float, default=12.0)
    parser.add_argument("--ly", type=float, default=2.0)
    parser.add_argument("--lz", type=float, default=2.0)
    parser.add_argument("--yz-bins", type=int, default=2)
    parser.add_argument("--n-max-x", type=int, default=30)
    parser.add_argument("--n-max-perp", type=int, default=5)
    parser.add_argument("--save-every", type=int, default=100)
    args = parser.parse_args()
    if not 0 <= args.shard_index < args.shard_count:
        parser.error("require 0 <= shard-index < shard-count")

    selected = args.campaign or list(CAMPAIGNS)
    files = []
    for campaign in selected:
        directory, target = CAMPAIGNS[campaign]
        files.extend((campaign, target, path) for path in sorted(directory.glob("*.h5")))
    files = files[args.shard_index::args.shard_count]
    k_matrix, mode_matrix = build_3d_modes(
        args.lx, args.ly, args.lz, args.n_max_x, args.n_max_perp
    )
    rows = []
    for index, (campaign, target, path) in enumerate(files, start=1):
        records = analyze(path, campaign, target, args, k_matrix, mode_matrix)
        rows.extend(records)
        print(
            f"[{index}/{len(files)}] {path.name}: {len(records)} late frames",
            flush=True,
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} rows to {args.output}", flush=True)


if __name__ == "__main__":
    main()
