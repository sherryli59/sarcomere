#!/usr/bin/env python3
"""Plot AA occupancy and status-2 duty across the AM status-memory sweep."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


PATTERN = re.compile(r"(?P<geometry>free|fixed)_memory(?P<memory>0|1ms|5ms|20ms)$")
MEMORY_MS = {"0": 0.0, "1ms": 1.0, "5ms": 5.0, "20ms": 20.0}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("summary_csv", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    grouped: dict[str, list[dict[str, float]]] = {"free": [], "fixed": []}
    with args.summary_csv.open(newline="") as stream:
        for row in csv.DictReader(stream):
            match = PATTERN.fullmatch(row["condition"])
            if match is None:
                continue
            grouped[match.group("geometry")].append({
                "memory": MEMORY_MS[match.group("memory")],
                "bond": float(row["mean_aa_bond_fraction"]),
                "bond_sem": float(row["sem_aa_bond_fraction"]),
                "status2": float(row["mean_status2_fraction_when_bonded"]),
                "status2_sem": float(row["sem_status2_fraction_when_bonded"]),
            })

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.5), sharex=True, sharey=True)
    colors = {"free": "#0072B2", "fixed": "#D55E00"}
    for geometry in ("free", "fixed"):
        rows = sorted(grouped[geometry], key=lambda row: row["memory"])
        x = np.asarray([row["memory"] for row in rows])
        for axis, field, sem_field in (
            (axes[0], "bond", "bond_sem"),
            (axes[1], "status2", "status2_sem"),
        ):
            axis.errorbar(
                x, [row[field] for row in rows],
                yerr=[row[sem_field] for row in rows],
                marker="o", linewidth=1.8, capsize=3,
                color=colors[geometry], label=geometry,
            )
    axes[0].set_title("AA bond occupancy")
    axes[1].set_title("Status-2 duty when bonded")
    axes[0].set_ylabel("Fraction of time")
    for axis in axes:
        axis.set_xlabel("Status-2 memory (ms)")
        axis.set_xticks([0, 1, 5, 20])
        axis.set_ylim(-0.02, 1.02)
        axis.spines[["top", "right"]].set_visible(False)
    axes[0].legend(frameon=False)
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=220, bbox_inches="tight")
    print(args.output)


if __name__ == "__main__":
    main()
