#!/usr/bin/env python3
"""Audit tension-dependent myosin OFF/SRX <-> detached-ON recruitment."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


def analyze(path: Path, dt: float, tail_fraction: float) -> dict[str, float | str]:
    with h5py.File(path, "r") as handle:
        meta = handle["meta"].attrs
        tension = float(meta["am_thick_force_clamp"])
        recruit0 = float(meta["am_thick_recruit_rate"])
        fold0 = float(meta["am_thick_fold_rate"])
        scale = float(meta["am_thick_tension_scale"])
        steps = np.asarray(handle["state/current_step"][:, 0], dtype=float)
        on_minus = np.asarray(
            handle["myosin/am_on_count_minus"][:], dtype=float).reshape(len(steps), -1)
        on_plus = np.asarray(
            handle["myosin/am_on_count_plus"][:], dtype=float).reshape(len(steps), -1)
        off_minus = np.asarray(
            handle["myosin/am_off_count_minus"][:], dtype=float).reshape(len(steps), -1)
        off_plus = np.asarray(
            handle["myosin/am_off_count_plus"][:], dtype=float).reshape(len(steps), -1)
        events = np.asarray(handle["actin_myo/thick_transitions"], dtype=float)

    on = on_minus.sum(axis=1) + on_plus.sum(axis=1)
    off = off_minus.sum(axis=1) + off_plus.sum(axis=1)
    total = on + off
    start = max(1, int((1.0 - tail_fraction) * len(steps)))
    tail_time_start = steps[start] * dt
    tail_time_end = steps[-1] * dt
    tail_events = events[events[:, 3] * dt >= tail_time_start]
    recruit_events = int(np.count_nonzero(
        (tail_events[:, 4] == 3) & (tail_events[:, 5] == 0)))
    fold_events = int(np.count_nonzero(
        (tail_events[:, 4] == 0) & (tail_events[:, 5] == 3)))

    time = steps * dt
    off_exposure = float(np.trapezoid(off[start:], time[start:]))
    on_exposure = float(np.trapezoid(on[start:], time[start:]))
    observed_recruit = recruit_events / off_exposure if off_exposure > 0 else math.nan
    observed_fold = fold_events / on_exposure if on_exposure > 0 else math.nan
    expected_recruit = recruit0 * math.exp(min(20.0, tension / scale))
    expected_fold = fold0 * math.exp(-min(20.0, tension / scale))
    expected_on = expected_recruit / (expected_recruit + expected_fold)

    return {
        "path": str(path),
        "tension": tension,
        "tail_start_time": tail_time_start,
        "tail_end_time": tail_time_end,
        "observed_on_fraction": float(np.mean(on[start:] / total[start:])),
        "expected_on_fraction": expected_on,
        "observed_recruit_rate": observed_recruit,
        "expected_recruit_rate": expected_recruit,
        "observed_fold_rate": observed_fold,
        "expected_fold_rate": expected_fold,
        "recruit_events": recruit_events,
        "fold_events": fold_events,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--dt", type=float, required=True)
    parser.add_argument("--tail-fraction", type=float, default=0.5)
    parser.add_argument("--output-prefix", type=Path, required=True)
    args = parser.parse_args()

    rows = sorted(
        (analyze(path, args.dt, args.tail_fraction) for path in args.inputs),
        key=lambda row: float(row["tension"]),
    )
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_prefix.with_suffix(".csv")
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    tension = np.array([row["tension"] for row in rows], dtype=float)
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.1))
    comparisons = (
        ("observed_on_fraction", "expected_on_fraction", "ON fraction"),
        ("observed_recruit_rate", "expected_recruit_rate", "OFF→ON rate"),
        ("observed_fold_rate", "expected_fold_rate", "ON→OFF rate"),
    )
    for axis, (observed, expected, ylabel) in zip(axes, comparisons):
        axis.plot(tension, [row[expected] for row in rows], "k--", label="dictated")
        axis.plot(tension, [row[observed] for row in rows], "o-", label="measured")
        axis.set_xlabel("half-filament tension clamp")
        axis.set_ylabel(ylabel)
    axes[0].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(args.output_prefix.with_suffix(".png"), dpi=200)


if __name__ == "__main__":
    main()
