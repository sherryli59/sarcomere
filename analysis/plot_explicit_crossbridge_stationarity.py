#!/usr/bin/env python3
"""Visualize free, crosslinked, and pinned minimal cross-bridge mechanics."""

from __future__ import annotations

import argparse
import csv
import os
import re
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/sarcomere_mplconfig")
import h5py
import matplotlib.pyplot as plt
import numpy as np


PATTERN = re.compile(r"am_xlink_(?P<condition>.+)_r(?P<rep>\d+)\.h5$")
CONDITIONS = (
    "free_noaa", "free_persistentaa", "fixed_noaa", "fixed_persistentaa",
)
STYLE = {
    "free_noaa": ("#e69f00", "--", "Free, no AA"),
    "free_persistentaa": ("#0072b2", "-", "Free + persistent AA"),
    "fixed_noaa": ("#777777", "--", "Pinned, no AA"),
    "fixed_persistentaa": ("#7e57c2", "-", "Pinned + persistent AA"),
}


def run_metrics(path: Path, dt: float, box_x: float) -> dict[str, np.ndarray]:
    with h5py.File(path, "r", locking=False) as handle:
        committed = int(np.asarray(handle["/meta/committed_frames"]).reshape(-1)[0])
        steps = np.asarray(handle["/state/current_step"][:committed]).reshape(-1)
        centers = np.asarray(handle["/actin/center"][:committed], dtype=float)
        directions = np.asarray(handle["/actin/direction"][:committed], dtype=float)
        active_force = np.asarray(
            handle["/actin/am_motor_force"][:committed], dtype=float
        ).reshape((committed, -1))
        pair_counts = np.asarray(
            handle["/state/aa_pairs_current_count"][:committed]
        ).reshape(-1)
        pair_status = handle["/state/aa_status_current"]
        status2 = np.zeros(committed, dtype=float)
        bonded = pair_counts > 0
        for frame, count in enumerate(pair_counts):
            if count:
                status2[frame] = float(np.any(
                    np.asarray(pair_status[frame, : int(count), 0]) == 2
                ))

    displacement = np.diff(centers, axis=0)
    displacement[..., 0] -= box_x * np.rint(displacement[..., 0] / box_x)
    productive_step = np.sum(displacement * directions[:-1], axis=2)
    productive_path = np.concatenate((
        np.zeros((1, centers.shape[1])),
        np.cumsum(np.abs(productive_step), axis=0),
    ), axis=0).mean(axis=1)

    com_step = np.mean(displacement, axis=1)
    com_path = np.concatenate((
        np.zeros(1), np.cumsum(np.linalg.norm(com_step, axis=1))
    ))

    relative_step = displacement[:, 0] - displacement[:, 1]
    relative_path = np.concatenate((
        np.zeros(1), np.cumsum(np.linalg.norm(relative_step, axis=1))
    ))
    return {
        "time": steps * dt,
        "productive_path": productive_path,
        "com_path": com_path,
        "relative_path": relative_path,
        "active_force": np.mean(active_force, axis=1),
        "bonded": bonded.astype(float),
        "status2": status2,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--dt", type=float, default=1e-5)
    parser.add_argument("--box-x", type=float, default=4.0)
    args = parser.parse_args()

    grouped: defaultdict[str, list[dict[str, np.ndarray]]] = defaultdict(list)
    for path in sorted(args.input_dir.glob("am_xlink_*_r*.h5")):
        match = PATTERN.match(path.name)
        if match and match.group("condition") in CONDITIONS:
            grouped[match.group("condition")].append(
                run_metrics(path, args.dt, args.box_x)
            )
    missing = [condition for condition in CONDITIONS if not grouped[condition]]
    if missing:
        parser.error(f"missing conditions: {', '.join(missing)}")

    fields = ("productive_path", "com_path", "relative_path", "active_force",
              "bonded", "status2")
    summary: dict[str, dict[str, np.ndarray]] = {}
    csv_rows: list[dict] = []
    for condition in CONDITIONS:
        runs = grouped[condition]
        n_frames = min(len(run["time"]) for run in runs)
        result = {"time": runs[0]["time"][:n_frames]}
        for field in fields:
            values = np.stack([run[field][:n_frames] for run in runs])
            result[f"mean_{field}"] = np.mean(values, axis=0)
            result[f"sem_{field}"] = (
                np.std(values, axis=0, ddof=1) / np.sqrt(values.shape[0])
            )
        summary[condition] = result
        for frame in range(n_frames):
            row = {"condition": condition, "frame": frame,
                   "time": result["time"][frame], "n_runs": len(runs)}
            for field in fields:
                row[f"mean_{field}"] = result[f"mean_{field}"][frame]
                row[f"sem_{field}"] = result[f"sem_{field}"][frame]
            csv_rows.append(row)

    prefix = args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    with Path(f"{prefix}_timeseries.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(csv_rows[0]))
        writer.writeheader()
        writer.writerows(csv_rows)

    fig, axes = plt.subplots(2, 3, figsize=(13.8, 7.4), sharex=True)
    panels = (
        (axes[0, 0], "productive_path", "Cumulative axial path / actin"),
        (axes[0, 1], "relative_path", "Cumulative relative path"),
        (axes[0, 2], "com_path", "Cumulative pair COM path"),
        (axes[1, 0], "active_force", "Mean active AM force / actin"),
        (axes[1, 1], "bonded", "Probability AA pair is bonded"),
        (axes[1, 2], "status2", "Probability pair is status 2"),
    )
    for axis, field, ylabel in panels:
        for condition in CONDITIONS:
            color, linestyle, label = STYLE[condition]
            result = summary[condition]
            time = result["time"]
            mean = result[f"mean_{field}"]
            sem = result[f"sem_{field}"]
            axis.plot(time, mean, color=color, ls=linestyle, lw=2, label=label)
            axis.fill_between(time, mean - sem, mean + sem, color=color, alpha=0.14)
        axis.set_ylabel(ylabel)
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(alpha=0.2)
    for axis in axes[1]:
        axis.set_xlabel("Simulation time")
    axes[1, 1].set_ylim(-0.03, 1.03)
    axes[1, 2].set_ylim(-0.03, 1.03)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(Path(f"{prefix}_mechanics.png"), dpi=240, bbox_inches="tight")
    fig.savefig(Path(f"{prefix}_mechanics.svg"), bbox_inches="tight")
    plt.close(fig)

    # A compact endpoint view makes the partial rescue by AA easy to compare.
    endpoint_fields = ("productive_path", "relative_path", "active_force", "status2")
    fig, axes = plt.subplots(1, 4, figsize=(12.0, 3.4))
    x = np.arange(len(CONDITIONS))
    colors = [STYLE[condition][0] for condition in CONDITIONS]
    for axis, field in zip(axes, endpoint_fields):
        means = [summary[c][f"mean_{field}"][-1] for c in CONDITIONS]
        sems = [summary[c][f"sem_{field}"][-1] for c in CONDITIONS]
        axis.bar(x, means, yerr=sems, color=colors, alpha=0.85, capsize=3)
        axis.set_xticks(x, ("Free\nno AA", "Free\n+ AA", "Pinned\nno AA",
                           "Pinned\n+ AA"), rotation=25, ha="right")
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", alpha=0.2)
    axes[0].set_ylabel("Final axial path")
    axes[1].set_ylabel("Final relative path")
    axes[2].set_ylabel("Final active force")
    axes[3].set_ylabel("Final status-2 probability")
    axes[3].set_ylim(0.0, 1.05)
    fig.tight_layout()
    fig.savefig(Path(f"{prefix}_endpoints.png"), dpi=240, bbox_inches="tight")
    fig.savefig(Path(f"{prefix}_endpoints.svg"), bbox_inches="tight")
    plt.close(fig)
    print(f"conditions={len(grouped)} runs={sum(map(len, grouped.values()))} prefix={prefix}")


if __name__ == "__main__":
    main()
