#!/usr/bin/env python3
"""Aggregate the four-arm minimal crosslink diagnostic."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import h5py

from summarize_aa_flicker import read_events as read_aa_events, reconstruct
from summarize_am_transitions import build_dwells, read_events as read_am_events, read_last_step


PATTERN = re.compile(
    r"am_xlink_(?P<condition>[A-Za-z0-9_.-]+)_r(?P<rep>\d+)\.h5$"
)


def write_csv(path: Path, fields: tuple[str, ...], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("--dt", type=float, default=1e-5)
    parser.add_argument("--actin-length", type=float, default=1.0)
    parser.add_argument("--antiparallel-overlap", type=float, default=0.1)
    parser.add_argument("--lx", type=float, default=4.0)
    parser.add_argument("--output-prefix", type=Path)
    args = parser.parse_args()
    prefix = args.output_prefix or args.input_dir / "minimal_crosslink"

    per_run: list[dict] = []
    for path in sorted(args.input_dir.glob("am_xlink_*_r*.h5")):
        match = PATTERN.search(path.name)
        if match is None:
            continue
        final_step = read_last_step(path)
        am_dwells = build_dwells(read_am_events(path), args.dt, final_step)
        strong_exposure = sum(
            float(row["dwell_time"]) for row in am_dwells if row["state"] == "strong"
        )
        strong_detachments = sum(
            row["state"] == "strong" and row["exit_reason"] == "strong_detachment"
            for row in am_dwells
        )
        aa_events, aa_last_step = read_aa_events(path)
        status_rows, bond_rows = reconstruct(aa_events, aa_last_step, args.dt)
        status_time = sum(float(row["duration_time"]) for row in status_rows)
        status2_time = sum(
            float(row["duration_time"]) for row in status_rows if row["status"] == 2
        )
        bond_time = sum(float(row["duration_time"]) for row in bond_rows)
        observation_time = max(final_step * args.dt, args.dt)
        gaps = [row for row in status_rows if row["internal_status1_gap"] == 1]
        with h5py.File(path, "r", locking=False) as handle:
            committed = int(np.asarray(handle["/meta/committed_frames"]).reshape(-1)[0])
            centers = np.asarray(handle["/actin/center"][:committed], dtype=float)
            directions = np.asarray(handle["/actin/direction"][:committed], dtype=float)
            motor_force = np.asarray(
                handle["/actin/am_motor_force"][:committed], dtype=float
            ).reshape((committed, -1))
            actin_force_caps = np.asarray(
                handle["/actin/cumulative_force_cap_events"][:committed], dtype=float
            ).reshape((committed, -1))
            actin_torque_caps = np.asarray(
                handle["/actin/cumulative_torque_cap_events"][:committed], dtype=float
            ).reshape((committed, -1))
            myosin_force_caps = np.asarray(
                handle["/myosin/cumulative_force_cap_events"][:committed], dtype=float
            ).reshape((committed, -1))
            myosin_torque_caps = np.asarray(
                handle["/myosin/cumulative_torque_cap_events"][:committed], dtype=float
            ).reshape((committed, -1))
        delta = np.diff(centers, axis=0)
        # The minimal boxes are periodic only along x. Undo only unambiguous
        # one-frame crossings; the diagnostic displacement is then projected
        # onto the barbed-to-pointed actin direction.
        delta[..., 0] -= args.lx * np.rint(delta[..., 0] / args.lx)
        axial_steps = np.sum(delta * directions[:-1], axis=2)
        net_axial = np.sum(axial_steps, axis=0)
        path_axial = np.sum(np.abs(axial_steps), axis=0)

        # AA registry is defined by the two barbed ends, as in aa_energy().
        # Track the pair-internal error separately from collective translation.
        barbed = centers - 0.5 * args.actin_length * directions
        barbed_delta = barbed[:, 1] - barbed[:, 0]
        barbed_delta[:, 0] -= args.lx * np.rint(barbed_delta[:, 0] / args.lx)
        ds0 = np.sum(barbed_delta * directions[:, 0], axis=1)
        ds1 = np.sum(barbed_delta * directions[:, 1], axis=1)
        polarity = np.sum(directions[:, 0] * directions[:, 1], axis=1)
        target = args.antiparallel_overlap * 0.5 * (1.0 - polarity)
        registry_error = np.sqrt(
            0.5 * ((ds0 - target) ** 2 + (ds1 + target) ** 2))
        n_steps = max(final_step, 1)
        n_actin_updates = n_steps * centers.shape[1]
        n_myosin_updates = n_steps * myosin_force_caps.shape[1]
        per_run.append({
            "condition": match.group("condition"),
            "replicate": int(match.group("rep")),
            "strong_detachments": strong_detachments,
            "strong_exposure": strong_exposure,
            "strong_detachment_rate": strong_detachments / strong_exposure if strong_exposure else np.nan,
            "aa_bond_fraction": min(1.0, bond_time / observation_time),
            "status2_fraction_when_bonded": status2_time / status_time if status_time else np.nan,
            "n_internal_status1_gaps": len(gaps),
            "mean_internal_status1_gap": float(np.mean([row["duration_time"] for row in gaps])) if gaps else np.nan,
            "mean_abs_net_axial_displacement": float(np.mean(np.abs(net_axial))),
            "mean_axial_path_length": float(np.mean(path_axial)),
            "initial_barbed_registry_error": float(registry_error[0]),
            "final_barbed_registry_error": float(registry_error[-1]),
            "barbed_registry_error_path": float(np.sum(np.abs(np.diff(registry_error)))),
            "actin_force_cap_fraction": float(np.sum(actin_force_caps[-1]) / n_actin_updates),
            "actin_torque_cap_fraction": float(np.sum(actin_torque_caps[-1]) / n_actin_updates),
            "myosin_force_cap_fraction": float(np.sum(myosin_force_caps[-1]) / n_myosin_updates),
            "myosin_torque_cap_fraction": float(np.sum(myosin_torque_caps[-1]) / n_myosin_updates),
            "mean_saved_active_motor_force": float(np.mean(motor_force)),
            "final_saved_active_motor_force": float(np.mean(motor_force[-1])),
        })

    numeric_fields = (
        "strong_detachment_rate", "aa_bond_fraction",
        "status2_fraction_when_bonded", "n_internal_status1_gaps",
        "mean_internal_status1_gap",
        "mean_abs_net_axial_displacement", "mean_axial_path_length",
        "initial_barbed_registry_error", "final_barbed_registry_error",
        "barbed_registry_error_path", "actin_force_cap_fraction",
        "actin_torque_cap_fraction", "myosin_force_cap_fraction",
        "myosin_torque_cap_fraction",
        "mean_saved_active_motor_force", "final_saved_active_motor_force",
    )
    grouped: defaultdict[str, list[dict]] = defaultdict(list)
    for row in per_run:
        grouped[row["condition"]].append(row)
    summary: list[dict] = []
    for condition, rows in sorted(grouped.items()):
        result = {"condition": condition, "n_runs": len(rows)}
        for field in numeric_fields:
            values = np.asarray([row[field] for row in rows], dtype=float)
            finite = values[np.isfinite(values)]
            result[f"mean_{field}"] = float(np.mean(finite)) if finite.size else np.nan
            result[f"sem_{field}"] = (
                float(np.std(finite, ddof=1) / np.sqrt(finite.size))
                if finite.size > 1 else np.nan
            )
        total_detach = sum(row["strong_detachments"] for row in rows)
        total_exposure = sum(row["strong_exposure"] for row in rows)
        result["pooled_strong_detachment_rate"] = (
            total_detach / total_exposure if total_exposure else np.nan
        )
        summary.append(result)

    per_run_fields = (
        "condition", "replicate", "strong_detachments", "strong_exposure",
        "strong_detachment_rate", "aa_bond_fraction", "status2_fraction_when_bonded",
        "n_internal_status1_gaps", "mean_internal_status1_gap",
        "mean_abs_net_axial_displacement", "mean_axial_path_length",
        "initial_barbed_registry_error", "final_barbed_registry_error",
        "barbed_registry_error_path", "actin_force_cap_fraction",
        "actin_torque_cap_fraction", "myosin_force_cap_fraction",
        "myosin_torque_cap_fraction",
        "mean_saved_active_motor_force", "final_saved_active_motor_force",
    )
    summary_fields = tuple(summary[0].keys()) if summary else ("condition", "n_runs")
    write_csv(Path(f"{prefix}_per_run.csv"), per_run_fields, per_run)
    write_csv(Path(f"{prefix}_summary.csv"), summary_fields, summary)
    print(f"runs={len(per_run)} conditions={len(summary)}")


if __name__ == "__main__":
    main()
