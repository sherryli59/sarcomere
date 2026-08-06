#!/usr/bin/env python3
"""Audit signed per-head AM alignment drive and half-filament cancellation."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


DT = 1e-6
SAVE_EVERY = 200
BETA = 233.644
MYOSIN_TUMBLE_DIFFUSION = 0.0589838
N_MYOSINS = 250


def director(u: np.ndarray) -> np.ndarray:
    return np.linalg.eigh(u.T @ u / len(u))[1][:, -1]


def values(handle: h5py.File, name: str) -> np.ndarray:
    return np.asarray(handle[f"/myosin/{name}"])[:, :, 0]


def analyze(path: Path) -> tuple[list[dict], dict]:
    condition, seed_text = path.stem.rsplit("_seed", 1)
    seed = int(seed_text)
    with h5py.File(path, "r", locking=False) as handle:
        au = np.asarray(handle["/actin/direction"])
        mu = np.asarray(handle["/myosin/direction"])
        data = {}
        for metric in (
            "drive", "abs_drive", "active_drive", "tether_drive",
            "strong_drive", "head_count", "positive_count", "negative_count",
        ):
            for half in ("minus", "plus"):
                data[f"{metric}_{half}"] = values(
                    handle, f"am_align_{metric}_{half}")

    initial_misalignment = None
    rows = []
    for frame in range(len(au)):
        da, dm = director(au[frame]), director(mu[frame])
        misalignment = np.degrees(np.arccos(np.clip(abs(da @ dm), 0, 1)))
        if initial_misalignment is None:
            initial_misalignment = misalignment
        summed = {name: float(np.sum(array[frame])) for name, array in data.items()}
        total_drive = summed["drive_minus"] + summed["drive_plus"]
        total_abs = summed["abs_drive_minus"] + summed["abs_drive_plus"]
        half_abs = abs(summed["drive_minus"]) + abs(summed["drive_plus"])
        total_heads = summed["head_count_minus"] + summed["head_count_plus"]
        minus_is_dominant = summed["head_count_minus"] >= summed["head_count_plus"]
        dominant_drive = summed["drive_minus"] if minus_is_dominant else summed["drive_plus"]
        nondominant_drive = summed["drive_plus"] if minus_is_dominant else summed["drive_minus"]
        dominant_heads = summed["head_count_minus"] if minus_is_dominant else summed["head_count_plus"]
        nondominant_heads = summed["head_count_plus"] if minus_is_dominant else summed["head_count_minus"]
        positive = summed["positive_count_minus"] + summed["positive_count_plus"]
        negative = summed["negative_count_minus"] + summed["negative_count_plus"]
        rows.append({
            "condition": condition, "seed": seed, "frame": frame,
            "time_s": frame * DT * SAVE_EVERY,
            "am_misalignment_deg": misalignment,
            "delta_am_misalignment_deg": misalignment - initial_misalignment,
            "drive_minus": summed["drive_minus"],
            "drive_plus": summed["drive_plus"],
            "drive_total": total_drive,
            "dominant_half_drive": dominant_drive,
            "nondominant_half_drive": nondominant_drive,
            "dominant_half_heads": dominant_heads,
            "nondominant_half_heads": nondominant_heads,
            "drive_abs_total": total_abs,
            "active_drive_total": summed["active_drive_minus"] + summed["active_drive_plus"],
            "tether_drive_total": summed["tether_drive_minus"] + summed["tether_drive_plus"],
            "strong_drive_total": summed["strong_drive_minus"] + summed["strong_drive_plus"],
            "head_count": total_heads,
            "aligning_head_fraction": positive / total_heads if total_heads else np.nan,
            "opposing_head_fraction": negative / total_heads if total_heads else np.nan,
            "global_cancellation": 1.0 - abs(total_drive) / total_abs if total_abs else np.nan,
            "between_half_cancellation": 1.0 - abs(total_drive) / half_abs if half_abs else np.nan,
            "predicted_alignment_rate_deg_s": (
                BETA * MYOSIN_TUMBLE_DIFFUSION * total_drive / N_MYOSINS
                * 180.0 / np.pi),
        })

    predicted_delta = 0.0
    rows[0]["predicted_delta_am_misalignment_deg"] = 0.0
    frame_dt = DT * SAVE_EVERY
    for previous, current in zip(rows[:-1], rows[1:]):
        predicted_delta -= 0.5 * frame_dt * (
            previous["predicted_alignment_rate_deg_s"] +
            current["predicted_alignment_rate_deg_s"])
        current["predicted_delta_am_misalignment_deg"] = predicted_delta

    tail = rows[int(0.8 * len(rows)):]
    summary = {"condition": condition, "seed": seed}
    for metric in (
        "delta_am_misalignment_deg", "drive_minus", "drive_plus", "drive_total",
        "dominant_half_drive", "nondominant_half_drive",
        "dominant_half_heads", "nondominant_half_heads",
        "drive_abs_total", "active_drive_total", "tether_drive_total",
        "strong_drive_total", "head_count", "aligning_head_fraction",
        "opposing_head_fraction", "global_cancellation", "between_half_cancellation",
        "predicted_alignment_rate_deg_s",
    ):
        x = np.asarray([row[metric] for row in tail], dtype=float)
        summary[f"tail_{metric}"] = float(np.nanmean(x))
    summary["final_delta_am_misalignment_deg"] = rows[-1]["delta_am_misalignment_deg"]
    summary["final_predicted_delta_am_misalignment_deg"] = rows[-1][
        "predicted_delta_am_misalignment_deg"]
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
        raise SystemExit("No complete trajectories")

    for filename, records in (("timeseries.csv", all_rows), ("summary.csv", summaries)):
        with (args.output / filename).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)

    conditions = sorted({row["condition"] for row in all_rows})
    colors = dict(zip(conditions, plt.cm.tab10(np.linspace(0, 1, len(conditions)))))
    metrics = (
        ("delta_am_misalignment_deg", "Change in AM misalignment (deg)"),
        ("predicted_delta_am_misalignment_deg", "Drive-predicted AM change (deg)"),
        ("drive_total", "Signed alignment drive"),
        ("active_drive_total", "Active-motor alignment drive"),
        ("tether_drive_total", "Passive-tether alignment drive"),
        ("drive_abs_total", "Absolute per-head alignment drive"),
        ("global_cancellation", "Global head cancellation"),
        ("between_half_cancellation", "Between-half cancellation"),
        ("aligning_head_fraction", "Fraction of heads aligning"),
        ("head_count", "Applied bound heads"),
    )
    fig, axes = plt.subplots(5, 2, figsize=(11, 14), sharex=True)
    for condition in conditions:
        group = [row for row in all_rows if row["condition"] == condition]
        seeds = sorted({row["seed"] for row in group})
        for axis, (metric, label) in zip(axes.flat, metrics):
            curves = []
            for seed in seeds:
                run = [row for row in group if row["seed"] == seed]
                curves.append([row[metric] for row in run])
            x = np.asarray(curves, dtype=float)
            time = np.asarray([row["time_s"] for row in run])
            axis.plot(time, np.nanmean(x, axis=0), lw=2,
                      color=colors[condition], label=condition)
            axis.set_ylabel(label)
            axis.grid(alpha=0.18)
            axis.spines[["top", "right"]].set_visible(False)
    axes[-1, 0].set_xlabel("Production time (s)")
    axes[-1, 1].set_xlabel("Production time (s)")
    axes[0, 1].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(args.output / "alignment_projection_audit.png", dpi=220)
    fig.savefig(args.output / "alignment_projection_audit.svg")

    torque_conditions = [c for c in conditions if c.startswith("head_torque")]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for condition in torque_conditions:
        group = [row for row in all_rows if row["condition"] == condition]
        seeds = sorted({row["seed"] for row in group})
        time = np.asarray([row["time_s"] for row in group if row["seed"] == seeds[0]])
        for metric, ls, label in (("dominant_half_drive", "-", "occupied/dominant half"),
                                  ("nondominant_half_drive", "--", "nondominant half")):
            curves = [[row[metric] for row in group if row["seed"] == seed]
                      for seed in seeds]
            axes[0].plot(time, np.mean(curves, axis=0), ls=ls, lw=2,
                         color=colors[condition], label=f"{condition}: {label}")
        plus = np.concatenate([[row["nondominant_half_drive"] for row in group if row["seed"] == seed]
                               for seed in seeds])
        minus = np.concatenate([[row["dominant_half_drive"] for row in group if row["seed"] == seed]
                                for seed in seeds])
        axes[1].scatter(minus, plus, s=8, alpha=0.35,
                        color=colors[condition], label=condition)
    axes[0].axhline(0, color="0.5", ls=":")
    axes[0].set_xlabel("Production time (s)")
    axes[0].set_ylabel("Signed half-filament alignment drive")
    axes[0].legend(frameon=False, fontsize=7)
    axes[1].axhline(0, color="0.5", ls=":")
    axes[1].axvline(0, color="0.5", ls=":")
    axes[1].set_xlabel("Dominant-half drive")
    axes[1].set_ylabel("Nondominant-half drive")
    axes[1].legend(frameon=False, fontsize=8)
    for axis in axes:
        axis.grid(alpha=0.18)
        axis.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(args.output / "half_filament_cancellation.png", dpi=220)
    fig.savefig(args.output / "half_filament_cancellation.svg")


if __name__ == "__main__":
    main()
