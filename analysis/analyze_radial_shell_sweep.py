#!/usr/bin/env python3
"""Aggregate six-actin preferred-gap and radial-stiffness sweeps."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


PATTERN = re.compile(r"gap-(?P<gap>[0-9.]+)_kr-(?P<kr>[0-9.]+)_rep-(?P<rep>\d+)")


def material_vector(handle: h5py.File, prefix: str) -> np.ndarray:
    return np.stack([np.asarray(handle[f"myosin/{prefix}_{axis}"])[:, 0, 0]
                     for axis in "xyz"], axis=1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    args = parser.parse_args()
    rows = []
    for path in args.inputs:
        match = PATTERN.fullmatch(path.stem)
        if match is None: raise ValueError(path)
        gap = float(match.group("gap")); target = 0.025 + 0.015 + gap
        with h5py.File(path, "r") as handle:
            ac = np.asarray(handle["actin/center"], dtype=float)
            mc = np.asarray(handle["myosin/center"], dtype=float)[:, 0, :]
            axis = np.asarray(handle["myosin/direction"], dtype=float)[:, 0, :]
            e1 = material_vector(handle, "material_e1")
            e2 = material_vector(handle, "material_e2")
            bound = np.asarray(handle["state/am_head_bound_actin_current"], dtype=int)[:, 0, :]
        displacement = ac - mc[:, None, :]
        transverse = displacement - np.sum(displacement * axis[:, None, :], axis=2)[:, :, None] * axis[:, None, :]
        y = np.sum(transverse * e1[:, None, :], axis=2)
        z = np.sum(transverse * e2[:, None, :], axis=2)
        radius = np.sqrt(y * y + z * z)
        theta = np.mod(np.arctan2(z, y), 2 * np.pi)
        psi6 = np.abs(np.mean(np.exp(6j * theta), axis=1))
        ordered = np.sort(theta, axis=1)
        gaps = np.diff(np.concatenate([ordered, ordered[:, :1] + 2 * np.pi], axis=1), axis=1)
        gap_cv = gaps.std(axis=1) / gaps.mean(axis=1)
        distinct = np.array([len(set(frame[frame >= 0].tolist())) for frame in bound])
        tail = slice(max(1, int(0.8 * len(psi6))), None)
        rows.append({
            "path": str(path), "surface_gap": gap,
            "radial_stiffness": float(match.group("kr")), "replicate": int(match.group("rep")),
            "target_centerline": target, "psi6_initial": float(psi6[0]),
            "psi6_tail": float(psi6[tail].mean()),
            "gap_cv_tail": float(gap_cv[tail].mean()),
            "radius_tail": float(radius[tail].mean()),
            "radius_error_tail": float(np.abs(radius[tail] - target).mean()),
            "radius_cv_tail": float((radius[tail].std(axis=1) / radius[tail].mean(axis=1)).mean()),
            "distinct_bound_tail": float(distinct[tail].mean()),
        })
    rows.sort(key=lambda row: (row["surface_gap"], row["radial_stiffness"], row["replicate"]))
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    with args.output_prefix.with_suffix(".csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    aggregate = []
    for gap in sorted({row["surface_gap"] for row in rows}):
        for stiffness in sorted({row["radial_stiffness"] for row in rows}):
            subset = [row for row in rows if row["surface_gap"] == gap and row["radial_stiffness"] == stiffness]
            out = {"surface_gap": gap, "radial_stiffness": stiffness}
            for metric in ("psi6_tail", "gap_cv_tail", "radius_tail", "radius_error_tail", "radius_cv_tail", "distinct_bound_tail"):
                values = np.array([row[metric] for row in subset], dtype=float)
                out[f"{metric}_mean"] = float(values.mean()); out[f"{metric}_std"] = float(values.std(ddof=1))
            aggregate.append(out)
    aggregate_path = args.output_prefix.with_name(args.output_prefix.name + "_aggregate").with_suffix(".csv")
    with aggregate_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(aggregate[0])); writer.writeheader(); writer.writerows(aggregate)
    fig, axes = plt.subplots(1, 3, figsize=(9.5, 3.0))
    for axis_plot, metric in zip(axes, ("psi6_tail", "gap_cv_tail", "radius_error_tail")):
        for gap in sorted({row["surface_gap"] for row in aggregate}):
            subset = [row for row in aggregate if row["surface_gap"] == gap]
            axis_plot.errorbar([row["radial_stiffness"] for row in subset],
                               [row[f"{metric}_mean"] for row in subset],
                               yerr=[row[f"{metric}_std"] for row in subset], marker="o", label=f"gap={gap:g}")
        axis_plot.set_xscale("log"); axis_plot.set_xlabel("radial stiffness"); axis_plot.set_ylabel(metric.replace("_", " "))
    axes[0].legend(frameon=False); fig.tight_layout(); fig.savefig(args.output_prefix.with_suffix(".png"), dpi=200)


if __name__ == "__main__":
    main()
