#!/usr/bin/env python3
"""Reconstruct AA bond/status episodes and quantify status-2 flicker."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np


COLUMNS = (
    "actin_i", "actin_j", "step", "old_bond", "new_bond",
    "old_status", "new_status", "motor_force_i", "motor_force_j",
    "strong_count_i", "strong_count_j", "load_i", "load_j",
)


def read_events(path: Path, include_pending: bool = True) -> tuple[np.ndarray, int]:
    chunks: list[np.ndarray] = []
    last_step = 0
    with h5py.File(path, "r") as handle:
        if "/catch_bond/state_transitions" in handle:
            saved = np.asarray(handle["/catch_bond/state_transitions"], dtype=float)
            if saved.size:
                chunks.append(saved.reshape((-1, len(COLUMNS))))
        resume_group = None
        if include_pending:
            resume_group = "/resume" if "/resume" in handle else (
                "/resume_prev" if "/resume_prev" in handle else None
            )
        if resume_group is not None:
            pending_path = f"{resume_group}/aa_state_events_pending"
            if pending_path in handle:
                pending = np.asarray(handle[pending_path], dtype=float).reshape(-1)
                complete = pending.size - pending.size % len(COLUMNS)
                if complete:
                    chunks.append(pending[:complete].reshape((-1, len(COLUMNS))))
            step_path = f"{resume_group}/current_step"
            if step_path in handle:
                last_step = int(np.asarray(handle[step_path]).reshape(-1)[0])
        if last_step == 0 and "/state/current_step" in handle:
            values = np.asarray(handle["/state/current_step"]).reshape(-1)
            if values.size:
                last_step = int(values[-1])
    if not chunks:
        return np.empty((0, len(COLUMNS)), dtype=float), last_step
    events = np.concatenate(chunks, axis=0)
    events = events[np.lexsort((events[:, 1], events[:, 0], events[:, 2]))]
    return events, max(last_step, int(np.max(events[:, 2])))


def write_rows(path: Path, fields: tuple[str, ...], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def reconstruct(events: np.ndarray, last_step: int, dt: float) -> tuple[list[dict], list[dict]]:
    active: dict[tuple[int, int], dict] = {}
    episode_counter: defaultdict[tuple[int, int], int] = defaultdict(int)
    status_rows: list[dict] = []
    bond_rows: list[dict] = []

    def close_status(record: dict, end_step: int, context: dict, censored: int) -> None:
        duration = max(0, end_step - record["status_start"])
        row = {
            "actin_i": record["actin_i"], "actin_j": record["actin_j"],
            "bond_episode": record["episode"], "status": record["status"],
            "start_step": record["status_start"], "end_step": end_step,
            "duration_steps": duration, "duration_time": duration * dt,
            "entry_motor_force_i": record["entry_motor_force_i"],
            "entry_motor_force_j": record["entry_motor_force_j"],
            "exit_motor_force_i": context["motor_force_i"],
            "exit_motor_force_j": context["motor_force_j"],
            "entry_strong_count_i": record["entry_strong_count_i"],
            "entry_strong_count_j": record["entry_strong_count_j"],
            "exit_strong_count_i": context["strong_count_i"],
            "exit_strong_count_j": context["strong_count_j"],
            "censored": censored,
            "internal_status1_gap": 0,
        }
        status_rows.append(row)
        record["status_time"][record["status"]] += duration

    def close_bond(record: dict, end_step: int, censored: int) -> None:
        duration = max(0, end_step - record["bond_start"])
        bond_rows.append({
            "actin_i": record["actin_i"], "actin_j": record["actin_j"],
            "bond_episode": record["episode"], "start_step": record["bond_start"],
            "end_step": end_step, "duration_steps": duration,
            "duration_time": duration * dt, "initial_status": record["initial_status"],
            "final_status": record["status"], "n_status_changes": record["n_status_changes"],
            "status1_time": record["status_time"][1] * dt,
            "status2_time": record["status_time"][2] * dt,
            "censored": censored, "left_censored": record["left_censored"],
        })

    for event in events:
        context = dict(zip(COLUMNS, event))
        for name in COLUMNS[:7]:
            context[name] = int(round(context[name]))
        key = tuple(sorted((context["actin_i"], context["actin_j"])))
        old_bond, new_bond = context["old_bond"], context["new_bond"]
        old_status, new_status = context["old_status"], context["new_status"]
        step = context["step"]
        record = active.get(key)

        if old_bond == 0 and new_bond == 1:
            episode_counter[key] += 1
            active[key] = {
                "actin_i": key[0], "actin_j": key[1], "episode": episode_counter[key],
                "bond_start": step, "status_start": step, "status": new_status,
                "initial_status": new_status, "n_status_changes": 0,
                "status_time": defaultdict(int), "left_censored": 0,
                "entry_motor_force_i": context["motor_force_i"],
                "entry_motor_force_j": context["motor_force_j"],
                "entry_strong_count_i": context["strong_count_i"],
                "entry_strong_count_j": context["strong_count_j"],
            }
            continue

        if record is None and old_bond == 1:
            episode_counter[key] += 1
            record = {
                "actin_i": key[0], "actin_j": key[1], "episode": episode_counter[key],
                "bond_start": step, "status_start": step, "status": old_status,
                "initial_status": old_status, "n_status_changes": 0,
                "status_time": defaultdict(int), "left_censored": 1,
                "entry_motor_force_i": context["motor_force_i"],
                "entry_motor_force_j": context["motor_force_j"],
                "entry_strong_count_i": context["strong_count_i"],
                "entry_strong_count_j": context["strong_count_j"],
            }
            active[key] = record

        if record is None:
            continue
        if old_status != new_status or old_bond != new_bond:
            close_status(record, step, context, 0)
        if old_bond == 1 and new_bond == 0:
            close_bond(record, step, 0)
            active.pop(key, None)
        elif old_bond == 1 and new_bond == 1 and old_status != new_status:
            record["n_status_changes"] += 1
            record["status"] = new_status
            record["status_start"] = step
            record["entry_motor_force_i"] = context["motor_force_i"]
            record["entry_motor_force_j"] = context["motor_force_j"]
            record["entry_strong_count_i"] = context["strong_count_i"]
            record["entry_strong_count_j"] = context["strong_count_j"]

    terminal_context = {
        "motor_force_i": np.nan, "motor_force_j": np.nan,
        "strong_count_i": np.nan, "strong_count_j": np.nan,
    }
    for record in active.values():
        close_status(record, last_step, terminal_context, 1)
        close_bond(record, last_step, 1)

    grouped: defaultdict[tuple[int, int, int], list[int]] = defaultdict(list)
    for index, row in enumerate(status_rows):
        grouped[(row["actin_i"], row["actin_j"], row["bond_episode"])].append(index)
    for indices in grouped.values():
        ordered = sorted(indices, key=lambda idx: status_rows[idx]["start_step"])
        for position in range(1, len(ordered) - 1):
            previous = status_rows[ordered[position - 1]]
            current = status_rows[ordered[position]]
            following = status_rows[ordered[position + 1]]
            if previous["status"] == 2 and current["status"] == 1 and following["status"] == 2:
                current["internal_status1_gap"] = 1
    return status_rows, bond_rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("trajectory", type=Path)
    parser.add_argument("--dt", type=float, default=1e-5)
    parser.add_argument("--tau1", type=float, default=1e-3)
    parser.add_argument("--output-prefix", type=Path)
    parser.add_argument(
        "--saved-only", action="store_true",
        help="Ignore mutable pending events in /resume when reading an active writer",
    )
    args = parser.parse_args()
    if not args.dt > 0 or args.tau1 < 0:
        parser.error("--dt must be positive and --tau1 non-negative")
    prefix = args.output_prefix or args.trajectory.with_suffix("")

    events, last_step = read_events(
        args.trajectory, include_pending=not args.saved_only
    )
    status_rows, bond_rows = reconstruct(events, last_step, args.dt)
    transition_rows = []
    for event in events:
        row = dict(zip(COLUMNS, event))
        for name in COLUMNS[:7]:
            row[name] = int(round(row[name]))
        transition_rows.append(row)

    transition_fields = COLUMNS
    status_fields = (
        "actin_i", "actin_j", "bond_episode", "status", "start_step", "end_step",
        "duration_steps", "duration_time", "entry_motor_force_i", "entry_motor_force_j",
        "exit_motor_force_i", "exit_motor_force_j", "entry_strong_count_i",
        "entry_strong_count_j", "exit_strong_count_i", "exit_strong_count_j",
        "censored", "internal_status1_gap",
    )
    bond_fields = (
        "actin_i", "actin_j", "bond_episode", "start_step", "end_step",
        "duration_steps", "duration_time", "initial_status", "final_status",
        "n_status_changes", "status1_time", "status2_time", "censored", "left_censored",
    )
    gaps = [row for row in status_rows if row["internal_status1_gap"] == 1]
    gap_times = np.asarray([row["duration_time"] for row in gaps], dtype=float)
    total_status_time = sum(row["duration_time"] for row in status_rows)
    status2_time = sum(row["duration_time"] for row in status_rows if row["status"] == 2)
    metrics = [
        {"metric": "n_state_transitions", "value": len(transition_rows)},
        {"metric": "n_bond_episodes", "value": len(bond_rows)},
        {"metric": "n_status_episodes", "value": len(status_rows)},
        {"metric": "n_internal_status1_gaps", "value": len(gaps)},
        {"metric": "mean_internal_status1_gap", "value": float(np.mean(gap_times)) if gap_times.size else np.nan},
        {"metric": "median_internal_status1_gap", "value": float(np.median(gap_times)) if gap_times.size else np.nan},
        {"metric": "p_internal_gap_gt_tau1", "value": float(np.mean(gap_times > args.tau1)) if gap_times.size else np.nan},
        {"metric": "status2_time_fraction", "value": status2_time / total_status_time if total_status_time else np.nan},
        {"metric": "tau1", "value": args.tau1},
    ]

    write_rows(Path(f"{prefix}_aa_state_transitions.csv"), transition_fields, transition_rows)
    write_rows(Path(f"{prefix}_aa_status_episodes.csv"), status_fields, status_rows)
    write_rows(Path(f"{prefix}_aa_bond_episodes.csv"), bond_fields, bond_rows)
    write_rows(Path(f"{prefix}_aa_status1_gaps.csv"), status_fields, gaps)
    write_rows(Path(f"{prefix}_aa_flicker_summary.csv"), ("metric", "value"), metrics)
    print(
        f"events={len(transition_rows)} bonds={len(bond_rows)} "
        f"status_episodes={len(status_rows)} internal_gaps={len(gaps)}"
    )


if __name__ == "__main__":
    main()
