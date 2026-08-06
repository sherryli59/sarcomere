#!/usr/bin/env python3
"""Analyze the tilted one-actin/one-myosin explicit-head mechanics probe."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


NAME = re.compile(r"heads(?P<heads>\d+)_torque(?P<torque>true|false)_seed(?P<seed>\d+)\.h5")


def scalar_series(handle: h5py.File, path: str) -> np.ndarray:
    return np.asarray(handle[path], dtype=float).reshape((handle[path].shape[0], -1))[:, 0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    parser.add_argument("--dt", type=float, default=1e-5)
    parser.add_argument("--save-every", type=int, default=200)
    args = parser.parse_args()

    traces: list[dict] = []
    rows: list[dict] = []
    for path in sorted(args.directory.glob("*.h5")):
        match = NAME.fullmatch(path.name)
        if match is None:
            continue
        heads = int(match["heads"])
        seed = int(match["seed"])
        torque = match["torque"] == "true"
        with h5py.File(path, "r", locking=False) as handle:
            actin = np.asarray(handle["/actin/direction"], dtype=float)[:, 0, :]
            myosin = np.asarray(handle["/myosin/direction"], dtype=float)[:, 0, :]
            cosine = np.clip(np.abs(np.sum(actin * myosin, axis=1)), 0.0, 1.0)
            angle = np.degrees(np.arccos(cosine))
            strong = scalar_series(handle, "/actin/am_strong_count")
            weak = scalar_series(handle, "/actin/am_weak_count")
            force = scalar_series(handle, "/actin/am_motor_force")
        time = np.arange(angle.size) * args.save_every * args.dt
        late = time >= 0.75 * time[-1]
        traces.append({
            "heads": heads, "seed": seed, "torque": torque,
            "time": time, "angle": angle, "strong": strong,
        })
        rows.append({
            "heads": heads,
            "torque": int(torque),
            "seed": seed,
            "initial_angle_deg": angle[0],
            "final_angle_deg": angle[-1],
            "late_angle_deg": float(np.mean(angle[late])),
            "mean_strong_heads": float(np.mean(strong)),
            "late_strong_heads": float(np.mean(strong[late])),
            "mean_weak_heads": float(np.mean(weak)),
            "mean_motor_force": float(np.mean(force)),
        })

    fields = tuple(rows[0])
    with (args.directory / "multihead_alignment_runs.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    groups: dict[tuple[int, bool], list[dict]] = {}
    for row in rows:
        groups.setdefault((row["heads"], bool(row["torque"])), []).append(row)
    summary: list[dict] = []
    for (heads, torque), group in sorted(groups.items()):
        summary.append({
            "heads": heads,
            "torque": int(torque),
            "n": len(group),
            "final_angle_mean_deg": np.mean([x["final_angle_deg"] for x in group]),
            "final_angle_sd_deg": np.std([x["final_angle_deg"] for x in group], ddof=1) if len(group) > 1 else 0.0,
            "late_angle_mean_deg": np.mean([x["late_angle_deg"] for x in group]),
            "mean_strong_heads": np.mean([x["mean_strong_heads"] for x in group]),
            "late_strong_heads": np.mean([x["late_strong_heads"] for x in group]),
        })
    with (args.directory / "multihead_alignment_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(summary[0]))
        writer.writeheader()
        writer.writerows(summary)

    figure, axes = plt.subplots(1, 2, figsize=(10.2, 4.1))
    palette = plt.cm.viridis(np.linspace(0.08, 0.92, 5))
    torque_on_heads = sorted({x["heads"] for x in traces if x["torque"]})
    for color, heads in zip(palette, torque_on_heads):
        selected = [x for x in traces if x["torque"] and x["heads"] == heads]
        stack = np.stack([x["angle"] for x in selected])
        time = selected[0]["time"]
        mean = np.mean(stack, axis=0)
        sd = np.std(stack, axis=0, ddof=1)
        axes[0].plot(time, mean, color=color, label=f"{heads} slots")
        axes[0].fill_between(time, mean - sd, mean + sd, color=color, alpha=0.16, linewidth=0)
    control = [x for x in traces if not x["torque"]]
    if control:
        stack = np.stack([x["angle"] for x in control])
        time = control[0]["time"]
        axes[0].plot(time, np.mean(stack, axis=0), color="0.25", linestyle="--", label="35, torque off")
    axes[0].set(xlabel="time", ylabel="actin–myosin angle (deg)", ylim=(0, 25))
    axes[0].legend(frameon=False, fontsize=8)

    on_summary = [x for x in summary if x["torque"]]
    x = np.asarray([x["heads"] for x in on_summary], dtype=float)
    angle = np.asarray([x["final_angle_mean_deg"] for x in on_summary], dtype=float)
    angle_sd = np.asarray([x["final_angle_sd_deg"] for x in on_summary], dtype=float)
    strong = np.asarray([x["mean_strong_heads"] for x in on_summary], dtype=float)
    axes[1].errorbar(x, angle, yerr=angle_sd, color="#31688e", marker="o", capsize=3)
    axes[1].set_xscale("log", base=2)
    axes[1].set_xticks(x, [str(int(value)) for value in x])
    axes[1].set(xlabel="available head slots", ylabel="final angle (deg)", ylim=(0, 25))
    twin = axes[1].twinx()
    twin.plot(x, strong, color="#b43c39", marker="s", linestyle="--")
    twin.set_ylabel("mean strong heads", color="#b43c39")
    twin.tick_params(axis="y", colors="#b43c39")
    figure.tight_layout()
    figure.savefig(args.directory / "multihead_alignment.png", dpi=220)
    figure.savefig(args.directory / "multihead_alignment.svg")
    print(f"analyzed {len(rows)} runs; wrote {args.directory / 'multihead_alignment_summary.csv'}")


if __name__ == "__main__":
    main()
