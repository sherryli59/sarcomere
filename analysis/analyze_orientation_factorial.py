#!/usr/bin/env python3
"""Analyze global/local orientational response in the AM torque factorial."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


DT = 1.0e-6
SAVE_EVERY = 200
BOX = np.array([6.0, np.sqrt(2.0), np.sqrt(2.0)])


def q_metrics(u: np.ndarray) -> tuple[float, np.ndarray, float]:
    q = 1.5 * np.einsum("ni,nj->ij", u, u) / len(u) - 0.5 * np.eye(3)
    values, vectors = np.linalg.eigh(q)
    director = vectors[:, -1]
    angle_x = np.degrees(np.arccos(np.clip(abs(director[0]), 0.0, 1.0)))
    return float(values[-1]), director, float(angle_x)


def local_pair_order(center: np.ndarray, u: np.ndarray, cutoff: float = 1.0) -> float:
    delta = center[:, None, :] - center[None, :, :]
    delta -= BOX * np.rint(delta / BOX)
    distance2 = np.einsum("ijk,ijk->ij", delta, delta)
    keep = np.triu((distance2 < cutoff * cutoff) & (distance2 > 0.0), 1)
    i, j = np.nonzero(keep)
    if len(i) == 0:
        return float("nan")
    dot = np.einsum("ij,ij->i", u[i], u[j])
    return float(np.mean(0.5 * (3.0 * dot * dot - 1.0)))


def cap_fraction(dataset: h5py.Dataset, n_particles: int, n_steps: int) -> float:
    count = float(np.sum(np.asarray(dataset)[-1]))
    return count / (n_particles * n_steps)


def analyze(path: Path) -> tuple[list[dict], dict]:
    condition, seed_text = path.stem.rsplit("_seed", 1)
    seed = int(seed_text)
    with h5py.File(path, "r", locking=False) as handle:
        ac = np.asarray(handle["/actin/center"])
        au = np.asarray(handle["/actin/direction"])
        mc = np.asarray(handle["/myosin/center"])
        mu = np.asarray(handle["/myosin/direction"])
        mt = np.asarray(handle["/myosin/torque"])
        strong = np.asarray(handle["/actin/am_strong_count"])
        weak = np.asarray(handle["/actin/am_weak_count"])
        n_steps = (len(au) - 1) * SAVE_EVERY
        cap = {
            "actin_force_cap_fraction": cap_fraction(
                handle["/actin/cumulative_force_cap_events"], len(au[0]), n_steps),
            "actin_torque_cap_fraction": cap_fraction(
                handle["/actin/cumulative_torque_cap_events"], len(au[0]), n_steps),
            "myosin_force_cap_fraction": cap_fraction(
                handle["/myosin/cumulative_force_cap_events"], len(mu[0]), n_steps),
            "myosin_torque_cap_fraction": cap_fraction(
                handle["/myosin/cumulative_torque_cap_events"], len(mu[0]), n_steps),
        }

    rows = []
    au0, mu0 = au[0], mu[0]
    for frame in range(len(au)):
        sa, da, aa = q_metrics(au[frame])
        sm, dm, am = q_metrics(mu[frame])
        adisp = np.degrees(np.arccos(np.clip(np.abs(np.einsum("ij,ij->i", au0, au[frame])), 0, 1)))
        mdisp = np.degrees(np.arccos(np.clip(np.abs(np.einsum("ij,ij->i", mu0, mu[frame])), 0, 1)))
        rows.append({
            "condition": condition, "seed": seed, "frame": frame,
            "time_s": frame * DT * SAVE_EVERY,
            "actin_nematic": sa, "myosin_nematic": sm,
            "actin_director_angle_x_deg": aa,
            "myosin_director_angle_x_deg": am,
            "am_director_misalignment_deg": float(np.degrees(
                np.arccos(np.clip(abs(np.dot(da, dm)), 0.0, 1.0)))),
            "actin_local_order": local_pair_order(ac[frame], au[frame]),
            "myosin_local_order": local_pair_order(mc[frame], mu[frame]),
            "actin_angular_displacement_deg": float(np.mean(adisp)),
            "myosin_angular_displacement_deg": float(np.mean(mdisp)),
            "myosin_torque_rms": float(np.sqrt(np.mean(np.sum(mt[frame] ** 2, axis=1)))),
            "weak_heads": int(np.sum(weak[frame])),
            "strong_heads": int(np.sum(strong[frame])),
        })

    tail = rows[max(0, int(0.8 * len(rows))):]
    summary = {"condition": condition, "seed": seed}
    for key in (
        "actin_nematic", "myosin_nematic", "actin_local_order", "myosin_local_order",
        "actin_director_angle_x_deg", "myosin_director_angle_x_deg",
        "am_director_misalignment_deg",
        "actin_angular_displacement_deg", "myosin_angular_displacement_deg",
        "myosin_torque_rms", "weak_heads", "strong_heads",
    ):
        summary[f"initial_{key}"] = rows[0][key]
        summary[f"final_{key}"] = rows[-1][key]
        summary[f"tail_{key}"] = float(np.mean([row[key] for row in tail]))
    summary.update(cap)
    return rows, summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    all_rows, summaries = [], []
    for path in sorted(args.input.glob("*.h5")):
        with h5py.File(path, "r", locking=False) as handle:
            if handle["/actin/direction"].shape[0] != 101:
                print(f"Skipping incomplete trajectory: {path}")
                continue
        rows, summary = analyze(path)
        all_rows.extend(rows)
        summaries.append(summary)
    if not summaries:
        raise SystemExit(f"No HDF5 trajectories found in {args.input}")

    for name, records in (("timeseries.csv", all_rows), ("summary.csv", summaries)):
        with (args.output / name).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)

    conditions = sorted({row["condition"] for row in all_rows})
    colors = plt.cm.tab10(np.linspace(0, 1, len(conditions)))
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.2), sharex=True)
    metrics = (
        ("actin_nematic", "Global actin nematic"),
        ("myosin_nematic", "Global myosin nematic"),
        ("actin_local_order", "Local actin pair order"),
        ("myosin_local_order", "Local myosin pair order"),
    )
    for condition, color in zip(conditions, colors):
        group = [row for row in all_rows if row["condition"] == condition]
        seeds = sorted({row["seed"] for row in group})
        for axis, (metric, label) in zip(axes.flat, metrics):
            curves = []
            for seed in seeds:
                run = [row for row in group if row["seed"] == seed]
                curves.append([row[metric] for row in run])
            values = np.asarray(curves)
            time = np.asarray([row["time_s"] for row in run])
            axis.plot(time, np.mean(values, axis=0), color=color, lw=1.8, label=condition)
            axis.set_ylabel(label)
            axis.grid(alpha=0.18)
            axis.spines[["top", "right"]].set_visible(False)
    axes[1, 0].set_xlabel("Production time (s)")
    axes[1, 1].set_xlabel("Production time (s)")
    axes[0, 1].legend(frameon=False, fontsize=7, bbox_to_anchor=(1.03, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(args.output / "nematic_factorial.png", dpi=220)
    fig.savefig(args.output / "nematic_factorial.svg")

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.2), sharex=True)
    metrics = (
        ("am_director_misalignment_deg", "AM director misalignment (deg)"),
        ("myosin_angular_displacement_deg", "Myosin angular displacement (deg)"),
        ("myosin_torque_rms", "Myosin torque RMS"),
        ("strong_heads", "Strong AM heads"),
    )
    for condition, color in zip(conditions, colors):
        group = [row for row in all_rows if row["condition"] == condition]
        seeds = sorted({row["seed"] for row in group})
        for axis, (metric, label) in zip(axes.flat, metrics):
            curves = []
            for seed in seeds:
                run = [row for row in group if row["seed"] == seed]
                curves.append([row[metric] for row in run])
            values = np.asarray(curves)
            time = np.asarray([row["time_s"] for row in run])
            axis.plot(time, np.mean(values, axis=0), color=color, lw=1.8,
                      label=condition)
            axis.set_ylabel(label)
            axis.grid(alpha=0.18)
            axis.spines[["top", "right"]].set_visible(False)
    axes[1, 0].set_xlabel("Production time (s)")
    axes[1, 1].set_xlabel("Production time (s)")
    axes[0, 1].legend(frameon=False, fontsize=7, bbox_to_anchor=(1.03, 1),
                      loc="upper left")
    fig.tight_layout()
    fig.savefig(args.output / "orientation_mechanics.png", dpi=220)
    fig.savefig(args.output / "orientation_mechanics.svg")

    tilted = [condition for condition in conditions
              if condition.startswith("alignedA_tiltedM_")]
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2))
    final_changes = []
    for condition, color in zip(tilted, plt.cm.tab10(np.linspace(0, 1, len(tilted)))):
        group = [row for row in all_rows if row["condition"] == condition]
        seeds = sorted({row["seed"] for row in group})
        curves = []
        for seed in seeds:
            run = [row for row in group if row["seed"] == seed]
            curve = np.asarray([row["am_director_misalignment_deg"] for row in run])
            curves.append(curve - curve[0])
        values = np.asarray(curves)
        time = np.asarray([row["time_s"] for row in run])
        mean = np.mean(values, axis=0)
        sd = np.std(values, axis=0, ddof=1)
        axes[0].plot(time, mean, color=color, lw=2, label=condition)
        axes[0].fill_between(time, mean - sd, mean + sd, color=color,
                             alpha=0.12, linewidth=0)
        final_changes.append((condition, mean[-1], sd[-1], color))
    axes[0].axhline(0, color="0.45", lw=1, ls="--")
    axes[0].set_xlabel("Production time (s)")
    axes[0].set_ylabel("Change in AM misalignment (deg)")
    axes[0].grid(alpha=0.18)
    labels = [item[0].removeprefix("alignedA_tiltedM_") for item in final_changes]
    x = np.arange(len(labels))
    axes[1].bar(x, [item[1] for item in final_changes],
                yerr=[item[2] for item in final_changes],
                color=[item[3] for item in final_changes], capsize=3)
    axes[1].axhline(0, color="0.45", lw=1, ls="--")
    axes[1].set_xticks(x, labels, rotation=35, ha="right")
    axes[1].set_ylabel("Final change in AM misalignment (deg)")
    axes[1].grid(axis="y", alpha=0.18)
    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
    axes[0].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(args.output / "tilted_alignment_comparison.png", dpi=220)
    fig.savefig(args.output / "tilted_alignment_comparison.svg")


if __name__ == "__main__":
    main()
