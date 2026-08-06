#!/usr/bin/env python3
"""Check pulse and diastolic statistics across integration timesteps."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


PATTERN = re.compile(r"dt-(?P<dt>[0-9.e+-]+)_rep-(?P<rep>\d+)")


def flattened_sum(handle: h5py.File, name: str) -> np.ndarray:
    values = np.asarray(handle[name], dtype=float)
    return values.reshape(values.shape[0], -1).sum(axis=1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    args = parser.parse_args()
    rows = []
    for path in args.inputs:
        match = PATTERN.search(path.stem)
        if match is None:
            raise ValueError(path)
        dt = float(match.group("dt"))
        with h5py.File(path, "r") as handle:
            time = np.asarray(handle["state/current_step"][:, 0], dtype=float) * dt
            q = flattened_sum(handle, "actin/am_q")
            strong = flattened_sum(handle, "actin/am_strong_count")
        tail = time >= time[-1] - 0.1
        rows.append({
            "path": str(path), "dt": dt, "replicate": int(match.group("rep")),
            "q_peak": float(q.max()), "q_integral": float(np.trapezoid(q, time)),
            "q_diastolic_tail": float(np.mean(q[tail])),
            "strong_diastolic_tail": float(np.mean(strong[tail])),
        })
    rows.sort(key=lambda row: (-row["dt"], row["replicate"]))
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    with args.output_prefix.with_suffix(".csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

    metrics = ("q_peak", "q_integral", "q_diastolic_tail", "strong_diastolic_tail")
    aggregate = []
    for dt in sorted({float(row["dt"]) for row in rows}, reverse=True):
        subset = [row for row in rows if float(row["dt"]) == dt]
        result: dict[str, float] = {"dt": dt}
        for metric in metrics:
            values = np.array([row[metric] for row in subset], dtype=float)
            result[f"{metric}_mean"] = float(values.mean())
            result[f"{metric}_std"] = float(values.std(ddof=1))
        aggregate.append(result)
    aggregate_path = args.output_prefix.with_name(
        args.output_prefix.name + "_aggregate").with_suffix(".csv")
    with aggregate_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(aggregate[0]))
        writer.writeheader(); writer.writerows(aggregate)

    fig, axes = plt.subplots(2, 2, figsize=(7.6, 5.8))
    for axis, metric in zip(axes.flat, metrics):
        x = np.arange(len(aggregate))
        axis.errorbar(x, [row[f"{metric}_mean"] for row in aggregate],
                      yerr=[row[f"{metric}_std"] for row in aggregate], fmt="o-")
        axis.set_xticks(x, [f'{row["dt"]:g}' for row in aggregate])
        axis.set_xlabel("dt")
        axis.set_ylabel(metric.replace("_", " "))
    fig.tight_layout()
    fig.savefig(args.output_prefix.with_suffix(".png"), dpi=200)


if __name__ == "__main__":
    main()
