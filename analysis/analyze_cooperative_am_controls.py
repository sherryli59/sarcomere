#!/usr/bin/env python3
"""Aggregate thin/thick feedback factorial and timestep-convergence runs."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


NAME = re.compile(
    r"(?P<load>isometric|lowload)_thin-(?P<thin>on|off)_thick-(?P<thick>on|off)"
    r"(?:_rep-(?P<rep>\d+))?"
)


def series(handle: h5py.File, name: str) -> np.ndarray:
    values = np.asarray(handle[name], dtype=float)
    return values.reshape(values.shape[0], -1).sum(axis=1)


def summarize(path: Path, dt: float) -> dict[str, float | str]:
    match = NAME.search(path.stem)
    if match is None:
        raise ValueError(f"Unrecognized factorial filename: {path}")
    with h5py.File(path, "r") as handle:
        step = np.asarray(handle["state/current_step"][:, 0], dtype=float)
        time = step * dt
        q = series(handle, "actin/am_q")
        strong = series(handle, "actin/am_strong_count")
        active_force = series(handle, "actin/am_motor_force")
        on = series(handle, "myosin/am_on_count_minus") + series(
            handle, "myosin/am_on_count_plus")
        off = series(handle, "myosin/am_off_count_minus") + series(
            handle, "myosin/am_off_count_plus")
        tension = series(handle, "myosin/am_tension_minus") + series(
            handle, "myosin/am_tension_plus")
    total_heads = np.maximum(on + off, 1.0)
    tail = time >= max(0.0, time[-1] - 0.1)
    labels = match.groupdict()
    labels["rep"] = int(labels["rep"] or 0)
    return {
        "path": str(path),
        **labels,
        "q_peak": float(q.max()),
        "q_integral": float(np.trapezoid(q, time)),
        "q_diastolic_tail": float(np.mean(q[tail])),
        "strong_peak": float(strong.max()),
        "strong_integral": float(np.trapezoid(strong, time)),
        "force_peak": float(active_force.max()),
        "force_integral": float(np.trapezoid(active_force, time)),
        "on_fraction_mean": float(np.mean(on / total_heads)),
        "on_fraction_peak": float(np.max(on / total_heads)),
        "tension_peak": float(tension.max()),
    }


def write_rows(path: Path, rows: list[dict[str, float | str]]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--dt", type=float, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    args = parser.parse_args()
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    rows = sorted(
        [summarize(path, args.dt) for path in args.inputs],
        key=lambda row: (row["load"], row["thin"], row["thick"], row["rep"]),
    )
    write_rows(args.output_prefix.with_suffix(".csv"), rows)

    metrics = ["q_integral", "strong_integral", "force_integral", "on_fraction_mean"]
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 6.2))
    labels = ["neither", "thin", "thick", "both"]
    states = [("off", "off"), ("on", "off"), ("off", "on"), ("on", "on")]
    for axis, metric in zip(axes.flat, metrics):
        x = np.arange(4)
        for offset, load in zip((-0.18, 0.18), ("isometric", "lowload")):
            values = []
            errors = []
            for state in states:
                samples = np.array([
                    row[metric] for row in rows
                    if row["load"] == load and
                    (row["thin"], row["thick"]) == state
                ], dtype=float)
                values.append(samples.mean())
                errors.append(samples.std(ddof=1) if len(samples) > 1 else 0.0)
            axis.bar(x + offset, values, yerr=errors, width=0.34, label=load)
        axis.set_xticks(x, labels, rotation=20)
        axis.set_ylabel(metric.replace("_", " "))
    axes[0, 0].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(args.output_prefix.with_suffix(".png"), dpi=200)

    synergy_rows = []
    for load in ("isometric", "lowload"):
        for rep in sorted({int(row["rep"]) for row in rows if row["load"] == load}):
            lookup = {
                (row["thin"], row["thick"]): row for row in rows
                if row["load"] == load and int(row["rep"]) == rep
            }
            for metric in metrics:
                synergy = (float(lookup[("on", "on")][metric])
                           - float(lookup[("on", "off")][metric])
                           - float(lookup[("off", "on")][metric])
                           + float(lookup[("off", "off")][metric]))
                synergy_rows.append({"load": load, "rep": rep, "metric": metric,
                                     "delta_synergy": synergy})
    write_rows(args.output_prefix.with_name(args.output_prefix.name + "_synergy").with_suffix(".csv"), synergy_rows)


if __name__ == "__main__":
    main()
