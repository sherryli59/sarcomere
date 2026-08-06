#!/usr/bin/env python3
"""Legacy all-myosin Phi_3D for the A2100 tau--lambda sweep."""

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

from extract_structure_3d import fourier_order_3d_frames, print_compute_device_info, resolve_compute_device


PATTERN = re.compile(
    r"^A2100_M1000_s2tau_tau(?P<tau_label>[^_]+)_lt(?P<lambda_label>[^_]+)"
    r"_r(?P<replica>\d+)_L12x2x2\.h5$"
)
RUN_FIELDS = (
    "tau_label", "tau_value", "lambda_label", "lambda_value", "replica", "filename",
    "checkpoint_step", "n_sampled_frames", "Phi_3D", "Theta_deg", "nx", "ny", "nz",
)
MEAN_FIELDS = (
    "tau_label", "tau_value", "lambda_label", "lambda_value", "n_replicates",
    "Phi_3D_mean", "Phi_3D_std", "Theta_deg_mean",
)


def value(label: str) -> float:
    return float(label)


def checkpoint_step(handle: h5py.File) -> int:
    for key in ("/resume/current_step", "/resume_prev/current_step"):
        if key in handle:
            return int(np.asarray(handle[key]).reshape(-1)[-1])
    return int(np.asarray(handle["/state/current_step"])[-1, 0])


def run_row(path: Path, args: argparse.Namespace, device) -> dict[str, object] | None:
    match = PATTERN.match(path.name)
    if match is None:
        return None
    with h5py.File(path, "r", locking=False) as handle:
        step = checkpoint_step(handle)
        if step < args.min_step or "/myosin/center" not in handle:
            return None
        centers = handle["/myosin/center"]
        indices = list(range(max(0, centers.shape[0] - args.eval_window), centers.shape[0], args.stride))
        results = []
        for start in range(0, len(indices), args.batch_size):
            coords = np.asarray(centers[indices[start:start + args.batch_size]], dtype=np.float64)
            results.extend(fourier_order_3d_frames(coords, device))

    modes = [(item["nx"], item["ny"], item["nz"]) for item in results]
    return {
        "tau_label": match["tau_label"],
        "tau_value": value(match["tau_label"]),
        "lambda_label": match["lambda_label"],
        "lambda_value": value(match["lambda_label"]),
        "replica": int(match["replica"]),
        "filename": path.name,
        "checkpoint_step": step,
        "n_sampled_frames": len(results),
        "Phi_3D": float(np.mean([item["phi_3d_max"] for item in results])),
        "Theta_deg": float(np.mean([item["theta_deg"] for item in results])),
        "nx": Counter(modes).most_common(1)[0][0][0],
        "ny": Counter(modes).most_common(1)[0][0][1],
        "nz": Counter(modes).most_common(1)[0][0][2],
    }


def write_csv(path: Path, fields: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def aggregate(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["tau_label"]), str(row["lambda_label"]))].append(row)
    output = []
    for (tau, lam), group in sorted(grouped.items(), key=lambda item: (value(item[0][0]), value(item[0][1]))):
        phi = np.asarray([float(row["Phi_3D"]) for row in group])
        theta = np.asarray([float(row["Theta_deg"]) for row in group])
        output.append({
            "tau_label": tau,
            "tau_value": value(tau),
            "lambda_label": lam,
            "lambda_value": value(lam),
            "n_replicates": len(group),
            "Phi_3D_mean": float(np.mean(phi)),
            "Phi_3D_std": float(np.std(phi, ddof=1)) if len(phi) > 1 else 0.0,
            "Theta_deg_mean": float(np.mean(theta)),
        })
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("/scratch/users/xsli/sarcomere/data_A2100_tau_lambda_s2only"))
    parser.add_argument("--per-run-output", type=Path, required=True)
    parser.add_argument("--mean-output", type=Path, required=True)
    parser.add_argument("--eval-window", type=int, default=300)
    parser.add_argument("--stride", type=int, default=10)
    parser.add_argument("--min-step", type=int, default=300000)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--device", choices=("auto", "cuda", "cpu", "numpy"), default="auto")
    args = parser.parse_args()

    device = resolve_compute_device(args.device)
    print_compute_device_info(device)
    rows = []
    for path in sorted(args.input_dir.glob("*.h5")):
        row = run_row(path, args, device)
        if row is not None:
            rows.append(row)
            print(f"tau={row['tau_label']} lambda={row['lambda_label']} r{row['replica']:02d}: Phi_3D={row['Phi_3D']:.4f}", flush=True)
    rows.sort(key=lambda row: (float(row["tau_value"]), float(row["lambda_value"]), int(row["replica"])))
    write_csv(args.per_run_output, RUN_FIELDS, rows)
    means = aggregate(rows)
    write_csv(args.mean_output, MEAN_FIELDS, means)
    print(f"wrote {len(rows)} complete trajectories and {len(means)} cells", flush=True)


if __name__ == "__main__":
    main()
