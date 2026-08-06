#!/usr/bin/env python3
"""Apply the legacy myosin Phi_3D definition to completed perturbation runs."""

from __future__ import annotations

import argparse
import csv
import os
from collections import Counter
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import h5py
import numpy as np

from extract_structure_3d import (
    fourier_order_3d_frames,
    print_compute_device_info,
    resolve_compute_device,
)
from summarize_perturbed_phi3d import CAMPAIGNS, checkpoint_step


FIELDS = (
    "campaign", "filename", "checkpoint_step", "n_sampled_frames",
    "Phi_3D", "Theta_deg", "nx", "ny", "nz",
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--campaign", action="append", choices=CAMPAIGNS)
    parser.add_argument("--eval-window", type=int, default=300)
    parser.add_argument("--stride", type=int, default=10)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument(
        "--device", choices=("auto", "cuda", "cpu", "numpy"), default="auto"
    )
    args = parser.parse_args()
    if not 0 <= args.shard_index < args.shard_count:
        parser.error("require 0 <= shard-index < shard-count")

    device = resolve_compute_device(args.device)
    if device is None:
        print("Using numpy for legacy Phi_3D.", flush=True)
    else:
        print_compute_device_info(device)

    rows = []
    selected = args.campaign or [name for name in CAMPAIGNS if name != "history_current"]
    for campaign in selected:
        directory, target = CAMPAIGNS[campaign]
        paths = sorted(directory.glob("*.h5"))[
            args.shard_index::args.shard_count
        ]
        for path in paths:
            with h5py.File(path, "r", locking=False) as handle:
                current = checkpoint_step(handle)
                if current < target or "/myosin/center" not in handle:
                    continue
                centers = handle["/myosin/center"]
                n_frames = centers.shape[0]
                frame_indices = list(
                    range(max(0, n_frames - args.eval_window), n_frames, args.stride)
                )
                results = []
                for start in range(0, len(frame_indices), args.batch_size):
                    indices = frame_indices[start:start + args.batch_size]
                    coords = np.asarray(centers[indices], dtype=np.float64)
                    results.extend(fourier_order_3d_frames(coords, device))
            modes = [(r["nx"], r["ny"], r["nz"]) for r in results]
            nx, ny, nz = Counter(modes).most_common(1)[0][0]
            rows.append({
                "campaign": campaign,
                "filename": path.name,
                "checkpoint_step": current,
                "n_sampled_frames": len(results),
                "Phi_3D": float(np.mean([r["phi_3d_max"] for r in results])),
                "Theta_deg": float(np.mean([r["theta_deg"] for r in results])),
                "nx": nx,
                "ny": ny,
                "nz": nz,
            })
            print(
                f"{campaign}: {path.name}: Phi_3D={rows[-1]['Phi_3D']:.4f} "
                f"over {len(results)} frames",
                flush=True,
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} runs to {args.output}", flush=True)


if __name__ == "__main__":
    main()
