#!/usr/bin/env python3
"""Aggregate calcium-only and strong-head-feedback activation sweeps."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


PATTERN = re.compile(r"kc-(?P<kc>[0-9.]+)_sb-(?P<sb>[0-9.]+)_rep-(?P<rep>\d+)")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    args = parser.parse_args()
    rows = []
    for path in args.inputs:
        match = PATTERN.fullmatch(path.stem)
        if match is None:
            raise ValueError(path)
        with h5py.File(path, "r") as handle:
            step = np.asarray(handle["state/current_step"][:, 0], dtype=float)
            q = np.asarray(handle["actin/am_q"], dtype=float).reshape(len(step), -1).mean(axis=1)
            strong = np.asarray(handle["actin/am_strong_count"], dtype=float).reshape(len(step), -1).sum(axis=1)
        time = step * 1e-4
        tail = time >= time[-1] - 0.1
        rows.append({
            "path": str(path), "k_calcium": float(match.group("kc")),
            "k_strong": float(match.group("sb")), "replicate": int(match.group("rep")),
            "q_peak": float(q.max()), "q_integral": float(np.trapezoid(q, time)),
            "q_diastolic_tail": float(q[tail].mean()),
            "strong_integral": float(np.trapezoid(strong, time)),
        })
    rows.sort(key=lambda row: (row["k_calcium"], row["k_strong"], row["replicate"]))
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    with args.output_prefix.with_suffix(".csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    aggregate = []
    for kc in sorted({row["k_calcium"] for row in rows}):
        for sb in sorted({row["k_strong"] for row in rows}):
            subset = [row for row in rows if row["k_calcium"] == kc and row["k_strong"] == sb]
            out = {"k_calcium": kc, "k_strong": sb}
            for metric in ("q_peak", "q_integral", "q_diastolic_tail", "strong_integral"):
                values = np.array([row[metric] for row in subset], dtype=float)
                out[f"{metric}_mean"] = float(values.mean())
                out[f"{metric}_std"] = float(values.std(ddof=1))
            aggregate.append(out)
    aggregate_path = args.output_prefix.with_name(args.output_prefix.name + "_aggregate").with_suffix(".csv")
    with aggregate_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(aggregate[0])); writer.writeheader(); writer.writerows(aggregate)
    fig, axes = plt.subplots(1, 3, figsize=(9.5, 3.0))
    for axis, metric in zip(axes, ("q_peak", "q_integral", "q_diastolic_tail")):
        for kc in sorted({row["k_calcium"] for row in aggregate}):
            subset = [row for row in aggregate if row["k_calcium"] == kc]
            axis.errorbar([row["k_strong"] for row in subset],
                          [row[f"{metric}_mean"] for row in subset],
                          yerr=[row[f"{metric}_std"] for row in subset], marker="o", label=f"kC={kc:g}")
        axis.set_xlabel("strong-feedback rate")
        axis.set_ylabel(metric.replace("_", " "))
    axes[0].legend(frameon=False)
    fig.tight_layout(); fig.savefig(args.output_prefix.with_suffix(".png"), dpi=200)


if __name__ == "__main__":
    main()
