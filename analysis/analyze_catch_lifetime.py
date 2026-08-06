#!/usr/bin/env python3
"""Estimate ACTN2 bond persistence with censoring-aware exposure/break rates."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


NAME = re.compile(
    r"heads(?P<heads>\d+)_memory(?P<memory>\d+(?:p\d+)?)_seed(?P<seed>\d+)_aa_bond_episodes\.csv"
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    parser.add_argument("--duration", type=float, default=0.2)
    args = parser.parse_args()

    runs: list[dict] = []
    for path in sorted(args.directory.glob("*_aa_bond_episodes.csv")):
        match = NAME.fullmatch(path.name)
        if match is None:
            continue
        with path.open() as stream:
            episodes = list(csv.DictReader(stream))
        exposure = sum(float(row["duration_time"]) for row in episodes)
        breaks = sum(int(row["censored"]) == 0 for row in episodes)
        status2_time = sum(float(row["status2_time"]) for row in episodes)
        runs.append({
            "heads": int(match["heads"]),
            "memory": float(match["memory"].replace("p", ".")),
            "seed": int(match["seed"]),
            "bonded_exposure": exposure,
            "breaks": breaks,
            "bond_occupancy": exposure / args.duration,
            "status2_fraction": status2_time / exposure if exposure else np.nan,
        })
    with (args.directory / "catch_lifetime_runs.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(runs[0]))
        writer.writeheader()
        writer.writerows(runs)

    groups: dict[tuple[int, float], list[dict]] = {}
    for row in runs:
        groups.setdefault((row["heads"], row["memory"]), []).append(row)
    summary: list[dict] = []
    for (heads, memory), group in sorted(groups.items()):
        exposure = sum(row["bonded_exposure"] for row in group)
        breaks = sum(row["breaks"] for row in group)
        # For an exponential hazard this exposure/event estimator naturally
        # includes right-censored episodes.
        lifetime = exposure / breaks if breaks else np.inf
        summary.append({
            "heads": heads,
            "memory": memory,
            "n_runs": len(group),
            "breaks": breaks,
            "bonded_exposure": exposure,
            "lifetime_mle": lifetime,
            "bond_occupancy_mean": np.mean([row["bond_occupancy"] for row in group]),
            "status2_fraction_mean": np.mean([row["status2_fraction"] for row in group]),
        })
    with (args.directory / "catch_lifetime_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(summary[0]))
        writer.writeheader()
        writer.writerows(summary)

    figure, axes = plt.subplots(1, 2, figsize=(9.5, 4.0))
    styles = {0.0: ("#3b528b", "no memory"), 0.005: ("#e07a28", "5 ms")}
    for memory in sorted({row["memory"] for row in summary}):
        selected = [row for row in summary if row["memory"] == memory]
        heads = np.asarray([row["heads"] for row in selected])
        color, label = styles.get(memory, (None, f"{memory:g}"))
        axes[0].plot(heads, [row["lifetime_mle"] for row in selected],
                     marker="o", color=color, label=label)
        axes[1].plot(heads, [row["bond_occupancy_mean"] for row in selected],
                     marker="o", color=color, label=label)
    for axis in axes:
        axis.set_xscale("log", base=2)
        axis.set_xticks([1, 8, 35], ["1", "8", "35"])
        axis.set_xlabel("available head slots per AM pair")
        axis.spines[["top", "right"]].set_visible(False)
    axes[0].set_yscale("log")
    axes[0].set_ylabel("bond lifetime estimate (exposure / breaks)")
    axes[1].set_ylabel("ACTN2 bond occupancy")
    axes[1].set_ylim(0.65, 1.02)
    axes[0].legend(frameon=False)
    figure.tight_layout()
    figure.savefig(args.directory / "catch_lifetime.png", dpi=220)
    figure.savefig(args.directory / "catch_lifetime.svg")
    print(f"runs={len(runs)} summary={args.directory / 'catch_lifetime_summary.csv'}")


if __name__ == "__main__":
    main()
