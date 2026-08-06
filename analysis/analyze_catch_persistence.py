#!/usr/bin/env python3
"""Aggregate ACTN2 status-2 persistence across explicit-head/memory probes."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


NAME = re.compile(
    r"heads(?P<heads>\d+)_memory(?P<memory>\d+(?:p\d+)?)_seed(?P<seed>\d+)_aa_flicker_summary\.csv"
)


def read_key_value(path: Path) -> dict[str, float]:
    with path.open() as stream:
        return {row["metric"]: float(row["value"]) for row in csv.DictReader(stream)}


def read_status_episodes(path: Path) -> list[dict]:
    with path.open() as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--zero-dir", type=Path, required=True)
    parser.add_argument("--memory-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    run_rows: list[dict] = []
    for directory in (args.zero_dir, args.memory_dir):
        for path in sorted(directory.glob("*_aa_flicker_summary.csv")):
            match = NAME.fullmatch(path.name)
            if match is None:
                continue
            memory = float(match["memory"].replace("p", "."))
            # The pre-fix directory contains placeholder nonzero-memory runs;
            # retain only its valid zero-memory controls.
            if directory == args.zero_dir and memory != 0.0:
                continue
            metrics = read_key_value(path)
            episode_path = Path(str(path).replace(
                "_aa_flicker_summary.csv", "_aa_status_episodes.csv"))
            episodes = read_status_episodes(episode_path)
            status2 = [row for row in episodes if int(row["status"]) == 2]
            observed_status2_dwell = (
                float(np.mean([float(row["duration_time"]) for row in status2]))
                if status2 else np.nan
            )
            run_rows.append({
                "heads": int(match["heads"]),
                "memory": memory,
                "seed": int(match["seed"]),
                "status2_fraction": metrics["status2_time_fraction"],
                "n_internal_status1_gaps": metrics["n_internal_status1_gaps"],
                "mean_internal_status1_gap": metrics["mean_internal_status1_gap"],
                "mean_observed_status2_dwell": observed_status2_dwell,
                "n_status2_episodes": len(status2),
            })

    run_fields = tuple(run_rows[0])
    with (args.output_dir / "catch_persistence_runs.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=run_fields)
        writer.writeheader()
        writer.writerows(run_rows)

    summary_rows: list[dict] = []
    groups: dict[tuple[int, float], list[dict]] = {}
    for row in run_rows:
        groups.setdefault((row["heads"], row["memory"]), []).append(row)
    for (heads, memory), group in sorted(groups.items()):
        fractions = np.asarray([row["status2_fraction"] for row in group])
        gaps = np.asarray([row["n_internal_status1_gaps"] for row in group])
        dwells = np.asarray([row["mean_observed_status2_dwell"] for row in group])
        summary_rows.append({
            "heads": heads,
            "memory": memory,
            "n": len(group),
            "status2_fraction_mean": np.mean(fractions),
            "status2_fraction_sd": np.std(fractions, ddof=1),
            "internal_gaps_mean": np.mean(gaps),
            "mean_observed_status2_dwell": np.nanmean(dwells),
        })
    with (args.output_dir / "catch_persistence_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)

    figure, axes = plt.subplots(1, 2, figsize=(9.8, 4.0))
    colors = {0.0: "#3b528b", 0.001: "#21918c", 0.005: "#e07a28"}
    labels = {0.0: "no memory", 0.001: "1 ms", 0.005: "5 ms"}
    for memory in sorted({row["memory"] for row in summary_rows}):
        selected = [row for row in summary_rows if row["memory"] == memory]
        x = np.asarray([row["heads"] for row in selected], dtype=float)
        y = np.asarray([row["status2_fraction_mean"] for row in selected])
        sd = np.asarray([row["status2_fraction_sd"] for row in selected])
        axes[0].errorbar(x, y, yerr=sd, marker="o", capsize=3,
                         color=colors[memory], label=labels[memory])
        gaps = np.asarray([row["internal_gaps_mean"] for row in selected])
        axes[1].plot(x, gaps, marker="o", color=colors[memory], label=labels[memory])
    for axis in axes:
        axis.set_xscale("log", base=2)
        axis.set_xticks([1, 8, 35], ["1", "8", "35"])
        axis.set_xlabel("available head slots per AM pair")
        axis.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("fraction of bonded time in status 2")
    axes[0].set_ylim(0.35, 1.03)
    axes[0].legend(frameon=False)
    axes[1].set_ylabel("internal status-1 interruptions per 0.2 time")
    figure.tight_layout()
    figure.savefig(args.output_dir / "catch_persistence.png", dpi=220)
    figure.savefig(args.output_dir / "catch_persistence.svg")
    print(f"runs={len(run_rows)} output={args.output_dir}")


if __name__ == "__main__":
    main()
