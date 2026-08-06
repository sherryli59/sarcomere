#!/usr/bin/env python3
"""Analyze matched local-ACTN2 motor-load and bundle-translation controls."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


PATTERN = re.compile(
    r"(?P<condition>common_linked|common_noaa|opposed_active|opposed_catch|opposed_free|opposed_nostrong)"
    r"_r(?P<rep>\d+)\.h5$"
)


def local_state_exposure(
    handle: h5py.File, dt: float
) -> tuple[float, float, int, int, int, int]:
    """Return exact status exposure, transitions and detachments."""
    if "/catch_bond/local_link_frame_offsets" not in handle:
        return 0.0, 0.0, 0, 0, 0, 0
    offsets = np.asarray(handle["/catch_bond/local_link_frame_offsets"], dtype=int)
    identity = np.asarray(handle["/catch_bond/local_link_identity"], dtype=int)
    active: dict[int, tuple[int, int]] = {}
    if len(offsets):
        _, offset, count = offsets[0]
        for row in identity[offset : offset + count]:
            link_id = (int(np.uint32(row[0])) << 32) | int(np.uint32(row[1]))
            active[link_id] = (int(row[4]), int(offsets[0, 0]))
    events = (
        np.asarray(handle["/catch_bond/local_link_events"], dtype=float)
        if "/catch_bond/local_link_events" in handle
        else np.empty((0, 11))
    )
    exposure = {1: 0.0, 2: 0.0}
    promotions = 0
    demotions = 0
    detachments = {1: 0, 2: 0}
    for row in events[np.argsort(events[:, 3], kind="stable")] if len(events) else []:
        link_id = int(row[0])
        step = int(row[3])
        code = int(row[4])
        old_status = int(row[5])
        new_status = int(row[6])
        if link_id in active:
            status, start = active[link_id]
            exposure[status] += max(0, step - start) * dt
        if code == 1:
            active[link_id] = (new_status, step)
        elif code == 2:
            if old_status in detachments:
                detachments[old_status] += 1
            active.pop(link_id, None)
        elif code == 3:
            if old_status == 1 and new_status == 2:
                promotions += 1
            if old_status == 2 and new_status == 1:
                demotions += 1
            active[link_id] = (new_status, step)
    final_step = int(np.asarray(handle["/resume/current_step"]).reshape(-1)[0])
    for status, start in active.values():
        exposure[status] += max(0, final_step - start) * dt
    return (
        exposure[1], exposure[2], promotions, demotions,
        detachments[1], detachments[2],
    )


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("--dt", type=float, default=1e-5)
    parser.add_argument("--lx", type=float, default=4.0)
    parser.add_argument("--output-prefix", type=Path)
    args = parser.parse_args()
    prefix = args.output_prefix or args.input_dir / "local_aa_coupling"

    per_run: list[dict] = []
    for path in sorted(args.input_dir.glob("*_r*.h5")):
        match = PATTERN.fullmatch(path.name)
        if match is None:
            continue
        with h5py.File(path, "r", locking=False) as handle:
            centers = np.asarray(handle["/actin/center"], dtype=float)
            dx = np.diff(centers[:, :, 0], axis=0)
            dx -= args.lx * np.rint(dx / args.lx)
            displacement = np.vstack(
                [np.zeros((1, centers.shape[1])), np.cumsum(dx, axis=0)]
            )
            com = displacement.mean(axis=1)
            relative = displacement[:, 1] - displacement[:, 0]
            motor_force = np.asarray(
                handle["/actin/am_motor_force"], dtype=float
            ).reshape((len(centers), -1))
            strong = np.asarray(
                handle["/actin/am_strong_count"], dtype=float
            ).reshape((len(centers), -1))
            tension = np.asarray(
                handle["/actin/aa_local_tension_max"], dtype=float
            ).reshape((len(centers), -1))
            force_caps = np.asarray(
                handle["/actin/cumulative_force_cap_events"], dtype=float
            ).reshape((len(centers), -1))
            torque_caps = np.asarray(
                handle["/actin/cumulative_torque_cap_events"], dtype=float
            ).reshape((len(centers), -1))
            final_step = int(
                np.asarray(handle["/resume/current_step"]).reshape(-1)[0]
            )
            (
                status1_time, status2_time, promotions, demotions,
                status1_detachments, status2_detachments,
            ) = local_state_exposure(handle, args.dt)
        total_link_time = status1_time + status2_time
        updates = max(1, final_step * centers.shape[1])
        per_run.append(
            {
                "condition": match.group("condition"),
                "replicate": int(match.group("rep")),
                "bundle_com_displacement": float(com[-1]),
                "bundle_com_path": float(np.sum(np.abs(np.diff(com)))),
                "relative_displacement": float(relative[-1]),
                "relative_path": float(np.sum(np.abs(np.diff(relative)))),
                "status1_exposure": status1_time,
                "status2_exposure": status2_time,
                "status2_fraction": (
                    status2_time / total_link_time if total_link_time else np.nan
                ),
                "status12_promotions": promotions,
                "status21_demotions": demotions,
                "status1_detachments": status1_detachments,
                "status2_detachments": status2_detachments,
                "status1_off_rate": (
                    status1_detachments / status1_time if status1_time else np.nan
                ),
                "status2_off_rate": (
                    status2_detachments / status2_time if status2_time else np.nan
                ),
                "mean_motor_force": float(np.mean(motor_force)),
                "mean_strong_heads": float(np.mean(strong)),
                "mean_link_tension": float(np.mean(np.max(tension, axis=1))),
                "actin_force_cap_fraction": float(np.sum(force_caps[-1]) / updates),
                "actin_torque_cap_fraction": float(np.sum(torque_caps[-1]) / updates),
            }
        )

    numeric = [key for key in per_run[0] if key not in {"condition", "replicate"}]
    grouped: defaultdict[str, list[dict]] = defaultdict(list)
    for row in per_run:
        grouped[row["condition"]].append(row)
    summary: list[dict] = []
    for condition, rows in sorted(grouped.items()):
        result: dict[str, float | int | str] = {
            "condition": condition,
            "n_runs": len(rows),
        }
        for field in numeric:
            values = np.asarray([row[field] for row in rows], dtype=float)
            values = values[np.isfinite(values)]
            result[f"mean_{field}"] = float(np.mean(values)) if len(values) else np.nan
            result[f"sem_{field}"] = (
                float(np.std(values, ddof=1) / np.sqrt(len(values)))
                if len(values) > 1
                else np.nan
            )
        total_status1_time = sum(float(row["status1_exposure"]) for row in rows)
        total_status2_time = sum(float(row["status2_exposure"]) for row in rows)
        total_status1_detach = sum(int(row["status1_detachments"]) for row in rows)
        total_status2_detach = sum(int(row["status2_detachments"]) for row in rows)
        result["pooled_status1_off_rate"] = (
            total_status1_detach / total_status1_time
            if total_status1_time else np.nan
        )
        result["pooled_status2_off_rate"] = (
            total_status2_detach / total_status2_time
            if total_status2_time else np.nan
        )
        summary.append(result)

    write_csv(Path(f"{prefix}_per_run.csv"), per_run)
    write_csv(Path(f"{prefix}_summary.csv"), summary)

    lookup = {row["condition"]: row for row in summary}
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.4))
    panels = [
        ("mean_bundle_com_displacement", "Bundle COM displacement (µm)"),
        ("mean_relative_path", "Internal relative path (µm)"),
        ("mean_status2_fraction", "Status-2 link-time fraction"),
    ]
    order = [
        "common_linked", "common_noaa", "opposed_active", "opposed_catch",
        "opposed_free", "opposed_nostrong",
    ]
    labels = [
        "Common +AA", "Common no AA", "Opp. λ=0.1",
        "Opp. λ=1", "Opp. free M", "Opp. no strong",
    ]
    for axis, (field, ylabel) in zip(axes, panels):
        values = [lookup[name][field] for name in order]
        sems = [lookup[name][field.replace("mean_", "sem_")] for name in order]
        axis.bar(
            range(6), values, yerr=sems,
            color=["#277da1", "#90be6d", "#f94144", "#d00000", "#f9844a", "#f9c74f"],
            capsize=3,
        )
        axis.set_xticks(range(6), labels, rotation=28, ha="right")
        axis.set_ylabel(ylabel)
        axis.axhline(0.0, color="black", linewidth=0.7)
    fig.tight_layout()
    fig.savefig(f"{prefix}.png", dpi=220)
    plt.close(fig)
    print(f"analyzed {len(per_run)} runs -> {prefix}_summary.csv")


if __name__ == "__main__":
    main()
