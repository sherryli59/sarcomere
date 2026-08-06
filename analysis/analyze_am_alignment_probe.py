#!/usr/bin/env python3
"""Analyze the tilted single-crossbridge off-center-torque probe."""

import csv
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent / "am_alignment_probe_20260720"
DT = 1.0e-5
SAVE_EVERY = 20
PROBE_SEEDS = {9101, 9102, 9103}


def load_run(path: Path):
    with h5py.File(path, "r") as handle:
        actin = np.asarray(handle["/actin/direction"])[:, 0, :]
        myosin = np.asarray(handle["/myosin/direction"])[:, 0, :]
        strong = np.asarray(handle["/actin/am_strong_count"])[:, 0, 0]
        weak = np.asarray(handle["/actin/am_weak_count"])[:, 0, 0]
        motor_force = np.asarray(handle["/actin/am_motor_force"])[:, 0, 0]

    dot = np.einsum("ij,ij->i", actin, myosin)
    dot /= np.linalg.norm(actin, axis=1) * np.linalg.norm(myosin, axis=1)
    return {
        "condition": path.stem.split("_seed", 1)[0],
        "seed": int(path.stem.rsplit("seed", 1)[1]),
        "time": np.arange(len(dot)) * DT * SAVE_EVERY,
        "nematic_angle_deg": np.degrees(
            np.arccos(np.clip(np.abs(dot), 0.0, 1.0))
        ),
        "abs_axis_dot": np.abs(dot),
        "weak": weak,
        "strong": strong,
        "motor_force": motor_force,
    }


runs = [
    load_run(path)
    for path in sorted(ROOT.glob("*.h5"))
    if int(path.stem.rsplit("seed", 1)[1]) in PROBE_SEEDS
]

timeseries_fields = [
    "condition", "seed", "frame", "time", "nematic_angle_deg",
    "abs_axis_dot", "weak", "strong", "motor_force",
]
with (ROOT / "am_alignment_probe_timeseries.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=timeseries_fields)
    writer.writeheader()
    for run in runs:
        for frame in range(len(run["time"])):
            writer.writerow({
                "condition": run["condition"],
                "seed": run["seed"],
                "frame": frame,
                "time": run["time"][frame],
                "nematic_angle_deg": run["nematic_angle_deg"][frame],
                "abs_axis_dot": run["abs_axis_dot"][frame],
                "weak": run["weak"][frame],
                "strong": run["strong"][frame],
                "motor_force": run["motor_force"][frame],
            })

summary = []
for run in runs:
    tail_start = int(0.8 * len(run["time"]))
    summary.append({
        "condition": run["condition"],
        "seed": run["seed"],
        "initial_angle_deg": run["nematic_angle_deg"][0],
        "final_angle_deg": run["nematic_angle_deg"][-1],
        "late_mean_angle_deg": np.mean(run["nematic_angle_deg"][tail_start:]),
        "angle_change_deg": (
            run["nematic_angle_deg"][-1] - run["nematic_angle_deg"][0]
        ),
        "strong_fraction": np.mean(run["strong"]),
        "mean_motor_force": np.mean(run["motor_force"]),
    })

summary_fields = list(summary[0])
with (ROOT / "am_alignment_probe_summary.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=summary_fields)
    writer.writeheader()
    writer.writerows(summary)

colors = {"torque": "#2a6fbb", "center": "#d06b36"}
labels = {
    "torque": "Off-center force (torque on)",
    "center": "Same force at center (torque off)",
}

fig, axes = plt.subplots(2, 1, figsize=(7.2, 6.4), sharex=True)
for condition in ["torque", "center"]:
    group = [run for run in runs if run["condition"] == condition]
    time = group[0]["time"]
    angles = np.stack([run["nematic_angle_deg"] for run in group], axis=1)
    forces = np.stack([run["motor_force"] for run in group], axis=1)
    color = colors[condition]

    for column in range(angles.shape[1]):
        axes[0].plot(time, angles[:, column], color=color, alpha=0.22, lw=0.8)
    mean_angle = np.mean(angles, axis=1)
    sd_angle = np.std(angles, axis=1, ddof=1)
    axes[0].plot(time, mean_angle, color=color, lw=2.2, label=labels[condition])
    axes[0].fill_between(time, mean_angle - sd_angle, mean_angle + sd_angle,
                         color=color, alpha=0.14, linewidth=0)

    mean_force = np.mean(forces, axis=1)
    sd_force = np.std(forces, axis=1, ddof=1)
    axes[1].plot(time, mean_force, color=color, lw=2.0)
    axes[1].fill_between(time, mean_force - sd_force, mean_force + sd_force,
                         color=color, alpha=0.14, linewidth=0)

axes[0].axhline(20.0, color="0.55", ls="--", lw=1.0, label="Initial angle")
axes[0].set_ylabel("Nematic AM angle (degrees)")
axes[0].set_ylim(bottom=0.0)
axes[0].legend(frameon=False, fontsize=9)
axes[0].set_title("Tilted single-crossbridge mechanics probe (3 replicates)")
axes[1].set_ylabel("Active motor force")
axes[1].set_xlabel("Simulation time")
axes[1].set_ylim(bottom=0.0)
for axis in axes:
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(alpha=0.18, linewidth=0.6)
fig.tight_layout()
fig.savefig(ROOT / "am_alignment_probe.png", dpi=220)
fig.savefig(ROOT / "am_alignment_probe.svg")

for row in summary:
    print(" ".join(f"{key}={value}" for key, value in row.items()))
print("\nCondition means:")
for condition in ["torque", "center"]:
    group = [row for row in summary if row["condition"] == condition]
    print(
        condition,
        "initial=", np.mean([row["initial_angle_deg"] for row in group]),
        "final=", np.mean([row["final_angle_deg"] for row in group]),
        "final_sd=", np.std([row["final_angle_deg"] for row in group], ddof=1),
        "late_mean=", np.mean([row["late_mean_angle_deg"] for row in group]),
        "change=", np.mean([row["angle_change_deg"] for row in group]),
        "strong_fraction=", np.mean([row["strong_fraction"] for row in group]),
        "motor_force=", np.mean([row["mean_motor_force"] for row in group]),
    )
