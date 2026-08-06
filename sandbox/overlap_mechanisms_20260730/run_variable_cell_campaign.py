#!/usr/bin/env python3
"""Convergence campaign with a dynamic periodic repeat length."""

from __future__ import annotations

import csv
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from run_campaign import tail_mean, tail_slope
from variable_cell_model import (
    default_variable_state,
    simulate_variable_cell,
    variable_cell_parameters,
)


OUT = Path(__file__).resolve().parent / "results"
MECHANISMS = (
    "baseline_turnover",
    "johann_entropic",
    "titin",
    "entropic_plus_titin",
    "register_plus_titin",
    "entropic_plus_series",
)
INITIALS = {
    "reference": (0.10, 4.30),
    "aa_low": (0.04, 4.30),
    "aa_high": (0.16, 4.30),
    "short_repeat": (0.10, 3.90),
    "long_repeat": (0.10, 4.70),
}
SEEDS = range(3)


def run_one(task):
    name, initial_name, y0, length0, seed = task
    cp = variable_cell_parameters(name)
    result = simulate_variable_cell(
        cp,
        seed,
        default_variable_state(y0, length0),
    )
    y = np.mean(result["overlaps"][:, :2], axis=1)
    x = np.mean(result["overlaps"][:, 2:], axis=1)
    length = result["length"]
    row = {
        "mechanism": name,
        "initial": initial_name,
        "seed": seed,
        "aa_tail_um": tail_mean(y),
        "am_tail_um": tail_mean(x),
        "length_tail_um": tail_mean(length),
        "aa_tail_velocity_um_s": tail_slope(result["time"], y),
        "am_tail_velocity_um_s": tail_slope(result["time"], x),
        "length_tail_velocity_um_s": tail_slope(
            result["time"], length
        ),
    }
    representative = result if initial_name == "reference" and seed == 0 else None
    return row, representative


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    representative = {}
    tasks = [
        (name, initial_name, y0, length0, seed)
        for name in MECHANISMS
        for initial_name, (y0, length0) in INITIALS.items()
        for seed in SEEDS
    ]
    with ProcessPoolExecutor(max_workers=4) as pool:
        for row, result in pool.map(run_one, tasks):
            rows.append(row)
            if result is not None:
                representative[row["mechanism"]] = result
            print(
                row["mechanism"], row["initial"], row["seed"], flush=True
            )

    with (OUT / "variable_cell_runs.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    summary = []
    for name in MECHANISMS:
        subset = [row for row in rows if row["mechanism"] == name]
        group_means = []
        group_velocities = []
        for initial_name in INITIALS:
            group = [row for row in subset if row["initial"] == initial_name]
            group_means.append(
                [
                    np.mean([row["aa_tail_um"] for row in group]),
                    np.mean([row["am_tail_um"] for row in group]),
                    np.mean([row["length_tail_um"] for row in group]),
                ]
            )
            group_velocities.append(
                [
                    np.mean([row["aa_tail_velocity_um_s"] for row in group]),
                    np.mean([row["am_tail_velocity_um_s"] for row in group]),
                    np.mean(
                        [row["length_tail_velocity_um_s"] for row in group]
                    ),
                ]
            )
        group_means = np.asarray(group_means)
        group_velocities = np.asarray(group_velocities)
        row = {
            "mechanism": name,
            "mean_aa_um": np.mean(group_means[:, 0]),
            "mean_am_um": np.mean(group_means[:, 1]),
            "mean_length_um": np.mean(group_means[:, 2]),
            "aa_ic_spread_um": np.ptp(group_means[:, 0]),
            "am_ic_spread_um": np.ptp(group_means[:, 1]),
            "length_ic_spread_um": np.ptp(group_means[:, 2]),
            "max_abs_group_tail_velocity_um_s": np.max(
                np.abs(group_velocities)
            ),
        }
        row["strict_pass"] = (
            row["aa_ic_spread_um"] < 0.025
            and row["am_ic_spread_um"] < 0.025
            and row["length_ic_spread_um"] < 0.05
            and row["max_abs_group_tail_velocity_um_s"] < 0.01
            and 0.06 < row["mean_aa_um"] < 0.14
            and 0.02 < row["mean_am_um"] < 0.655
        )
        summary.append(row)

    with (OUT / "variable_cell_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)

    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    for name, result in representative.items():
        time = result["time"]
        axes[0].plot(
            time, np.mean(result["overlaps"][:, :2], axis=1), label=name
        )
        axes[1].plot(
            time, np.mean(result["overlaps"][:, 2:], axis=1), label=name
        )
        axes[2].plot(time, result["length"], label=name)
    axes[0].axhline(0.10, color="black", ls=":")
    axes[0].set_ylabel("mean A–A overlap (µm)")
    axes[1].set_ylabel("mean A–M overlap (µm)")
    axes[2].set_ylabel("repeat length (µm)")
    axes[2].set_xlabel("time (s)")
    axes[0].legend(ncol=3, fontsize=8)
    fig.suptitle("Dynamic-repeat convergence under the same external load")
    fig.tight_layout()
    fig.savefig(OUT / "variable_cell_timeseries.png", dpi=180)

    # The decisive comparison: different starting repeat lengths under the
    # same boundary condition for the two combined mechanisms.
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    for col, name in enumerate(("entropic_plus_titin", "register_plus_titin")):
        cp = variable_cell_parameters(name)
        for label, (_, length0) in INITIALS.items():
            if label not in ("reference", "short_repeat", "long_repeat"):
                continue
            result = simulate_variable_cell(
                cp, 100, default_variable_state(0.10, length0)
            )
            axes[0, col].plot(
                result["time"], result["length"], label=f"{length0:.1f} µm"
            )
            axes[1, col].plot(
                result["time"],
                np.mean(result["overlaps"][:, 2:], axis=1),
                label=f"{length0:.1f} µm",
            )
        axes[0, col].set_title(name.replace("_", " "))
        axes[0, col].set_ylabel("repeat length (µm)")
        axes[1, col].set_ylabel("mean A–M overlap (µm)")
        axes[1, col].set_xlabel("time (s)")
        axes[0, col].legend()
    fig.tight_layout()
    fig.savefig(OUT / "variable_cell_length_convergence.png", dpi=180)


if __name__ == "__main__":
    main()
