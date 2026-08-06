#!/usr/bin/env python3
"""Summarize dynamic actin-myosin transition events from a trajectory."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np


COLUMNS = (
    "actin_id",
    "myosin_id",
    "step",
    "old_state",
    "new_state",
    "transition_reason",
    "activation",
    "axial_strain",
    "motor_force",
    "preceding_cb_status",
    "head_id",
)
LEGACY_WIDTH = len(COLUMNS) - 1
STATE_NAMES = {0: "detached", 1: "weak", 2: "strong"}
REASON_NAMES = {
    0: "attachment",
    1: "weak_detachment",
    2: "powerstroke",
    3: "strong_detachment",
    4: "geometric_rupture",
    5: "neighbor_loss",
    6: "capacity_limit",
}
RUPTURE_COLUMNS = (
    "actin_id", "myosin_id", "step", "head_id", "old_state",
    "transition_reason", "retention_failure_mask", "segment_distance",
    "rupture_distance", "actin_myosin_dot", "expected_half", "head_half",
    "productive_axial_overlap", "binding_ratio", "partial_ratio",
)
RETENTION_FAILURE_BITS = {
    1: "not_in_neighbor_list",
    2: "beyond_rupture_distance",
    4: "wrong_productive_half",
    8: "no_productive_axial_overlap",
    16: "association_gate_closed",
    32: "unexpected_candidate_loss",
}


def decode_retention_failure(mask: int) -> str:
    if mask == 0:
        return "none"
    return "|".join(
        name for bit, name in RETENTION_FAILURE_BITS.items() if mask & bit
    ) or f"unknown_mask_{mask}"


def read_forced_ruptures(path: Path, include_pending: bool = True) -> np.ndarray:
    chunks: list[np.ndarray] = []
    with h5py.File(path, "r") as handle:
        if "/actin_myo/forced_ruptures" in handle:
            saved = np.asarray(handle["/actin_myo/forced_ruptures"], dtype=float)
            if saved.size:
                chunks.append(saved.reshape((-1, len(RUPTURE_COLUMNS))))
        if include_pending:
            for group in ("/resume", "/resume_prev"):
                pending_path = f"{group}/am_forced_rupture_events_pending"
                if pending_path not in handle:
                    continue
                pending = np.asarray(handle[pending_path], dtype=float).reshape(-1)
                complete = pending.size - pending.size % len(RUPTURE_COLUMNS)
                if complete:
                    chunks.append(
                        pending[:complete].reshape((-1, len(RUPTURE_COLUMNS)))
                    )
                break
    if not chunks:
        return np.empty((0, len(RUPTURE_COLUMNS)), dtype=float)
    rows = np.concatenate(chunks, axis=0)
    return rows[np.lexsort((rows[:, 3], rows[:, 1], rows[:, 0], rows[:, 2]))]


def read_events(path: Path, include_pending: bool = True) -> np.ndarray:
    chunks: list[np.ndarray] = []
    with h5py.File(path, "r") as handle:
        if "/actin_myo/head_transitions" in handle:
            saved = np.asarray(handle["/actin_myo/head_transitions"], dtype=float)
            if saved.size:
                chunks.append(saved.reshape((-1, len(COLUMNS))))
        elif "/actin_myo/transitions" in handle:
            saved = np.asarray(handle["/actin_myo/transitions"], dtype=float)
            if saved.size:
                legacy = saved.reshape((-1, LEGACY_WIDTH))
                chunks.append(np.column_stack((legacy, np.zeros(len(legacy)))))
        # Events after the last trajectory-frame commit live in the resume
        # snapshot and have not yet been appended to /actin_myo/transitions.
        if include_pending:
            resume_path = None
            if "/resume/am_transition_events_pending" in handle:
                resume_path = "/resume/am_transition_events_pending"
            elif "/resume_prev/am_transition_events_pending" in handle:
                resume_path = "/resume_prev/am_transition_events_pending"
            if resume_path is not None:
                pending = np.asarray(handle[resume_path], dtype=float).reshape(-1)
                model_version = int(handle["/meta"].attrs.get("am_model_version", 2))
                width = len(COLUMNS) if model_version >= 3 else LEGACY_WIDTH
                complete = pending.size - pending.size % width
                if complete:
                    values = pending[:complete].reshape((-1, width))
                    if width == LEGACY_WIDTH:
                        values = np.column_stack((values, np.zeros(len(values))))
                    chunks.append(values)
    if not chunks:
        return np.empty((0, len(COLUMNS)), dtype=float)
    events = np.concatenate(chunks, axis=0)
    # The simulator merges thread-local buffers, so explicitly restore temporal
    # order while retaining deterministic pair ordering within a step.
    return events[np.lexsort((events[:, 10], events[:, 1], events[:, 0], events[:, 2]))]


def read_last_step(path: Path, include_resume: bool = True) -> int:
    with h5py.File(path, "r") as handle:
        if include_resume:
            for group in ("/resume", "/resume_prev"):
                step_path = f"{group}/current_step"
                if step_path in handle:
                    return int(np.asarray(handle[step_path]).reshape(-1)[-1])
        if "/state/current_step" in handle:
            values = np.asarray(handle["/state/current_step"]).reshape(-1)
            if values.size:
                return int(values[-1])
    return 0


def write_rows(path: Path, fieldnames: tuple[str, ...], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def build_dwells(events: np.ndarray, dt: float, final_step: int | None = None) -> list[dict]:
    active: dict[tuple[int, int], dict] = {}
    rows: list[dict] = []
    for event in events:
        values = dict(zip(COLUMNS, event))
        actin_id = int(round(values["actin_id"]))
        myosin_id = int(round(values["myosin_id"]))
        head_id = int(round(values["head_id"]))
        step = int(round(values["step"]))
        old_state = int(round(values["old_state"]))
        new_state = int(round(values["new_state"]))
        reason = int(round(values["transition_reason"]))
        cb_status = int(round(values["preceding_cb_status"]))
        key = (actin_id, myosin_id, head_id)

        entry = active.get(key)
        if entry is not None and entry["state"] == old_state and old_state != 0:
            dwell_steps = max(0, step - entry["start_step"])
            rows.append(
                {
                    "actin_id": actin_id,
                    "myosin_id": myosin_id,
                    "head_id": head_id,
                    "state": STATE_NAMES.get(old_state, str(old_state)),
                    "start_step": entry["start_step"],
                    "end_step": step,
                    "dwell_steps": dwell_steps,
                    "dwell_time": dwell_steps * dt,
                    "exit_state": STATE_NAMES.get(new_state, str(new_state)),
                    "exit_reason": REASON_NAMES.get(reason, str(reason)),
                    "entry_cb_status": entry["cb_status"],
                    "exit_cb_status": cb_status,
                    "entry_force": entry["force"],
                    "exit_force": values["motor_force"],
                    "entry_activation": entry["activation"],
                    "exit_activation": values["activation"],
                }
            )

        if new_state == 0:
            active.pop(key, None)
        else:
            active[key] = {
                "state": new_state,
                "start_step": step,
                "cb_status": cb_status,
                "force": values["motor_force"],
                "activation": values["activation"],
            }
    if final_step is not None:
        for (actin_id, myosin_id, head_id), entry in active.items():
            dwell_steps = max(0, final_step - entry["start_step"])
            rows.append({
                "actin_id": actin_id,
                "myosin_id": myosin_id,
                "head_id": head_id,
                "state": STATE_NAMES.get(entry["state"], str(entry["state"])),
                "start_step": entry["start_step"],
                "end_step": final_step,
                "dwell_steps": dwell_steps,
                "dwell_time": dwell_steps * dt,
                "exit_state": "censored",
                "exit_reason": "right_censored",
                "entry_cb_status": entry["cb_status"],
                "exit_cb_status": entry["cb_status"],
                "entry_force": entry["force"],
                "exit_force": np.nan,
                "entry_activation": entry["activation"],
                "exit_activation": np.nan,
            })
    return rows


def summarize_dwells(dwells: list[dict]) -> list[dict]:
    groups: dict[tuple[str, int, str], list[float]] = defaultdict(list)
    for row in dwells:
        groups[(row["state"], int(row["entry_cb_status"]), row["exit_reason"])].append(
            float(row["dwell_time"])
        )
    rows: list[dict] = []
    for (state, cb_status, exit_reason), samples in sorted(groups.items()):
        values = np.asarray(samples, dtype=float)
        positive = values[values > 0.0]
        mean = float(np.mean(values)) if values.size else np.nan
        rows.append(
            {
                "state": state,
                "entry_cb_status": cb_status,
                "exit_reason": exit_reason,
                "n_completed": int(values.size),
                "mean_dwell": mean,
                "median_dwell": float(np.median(values)),
                "std_dwell": float(np.std(values, ddof=1)) if values.size > 1 else 0.0,
                "p10_dwell": float(np.quantile(values, 0.1)),
                "p90_dwell": float(np.quantile(values, 0.9)),
                "event_rate_estimate": 1.0 / float(np.mean(positive)) if positive.size else np.nan,
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("trajectory", type=Path)
    parser.add_argument("--dt", type=float, default=1e-5)
    parser.add_argument("--output-prefix", type=Path)
    parser.add_argument(
        "--saved-only", action="store_true",
        help="Ignore mutable pending events in /resume when reading an active writer",
    )
    args = parser.parse_args()
    if not args.dt > 0.0:
        parser.error("--dt must be positive")

    prefix = args.output_prefix or args.trajectory.with_suffix("")
    events = read_events(args.trajectory, include_pending=not args.saved_only)
    transition_rows = []
    for event in events:
        row = dict(zip(COLUMNS, event))
        for key in (
            "actin_id", "myosin_id", "step", "old_state", "new_state",
            "transition_reason", "preceding_cb_status",
            "head_id",
        ):
            row[key] = int(round(row[key]))
        row["old_state_name"] = STATE_NAMES.get(row["old_state"], str(row["old_state"]))
        row["new_state_name"] = STATE_NAMES.get(row["new_state"], str(row["new_state"]))
        row["transition_reason_name"] = REASON_NAMES.get(
            row["transition_reason"], str(row["transition_reason"])
        )
        transition_rows.append(row)

    rupture_rows = []
    for event in read_forced_ruptures(
        args.trajectory, include_pending=not args.saved_only
    ):
        row = dict(zip(RUPTURE_COLUMNS, event))
        for key in (
            "actin_id", "myosin_id", "step", "head_id", "old_state",
            "transition_reason", "retention_failure_mask", "expected_half",
            "head_half",
        ):
            row[key] = int(round(row[key]))
        row["old_state_name"] = STATE_NAMES.get(
            row["old_state"], str(row["old_state"])
        )
        row["transition_reason_name"] = REASON_NAMES.get(
            row["transition_reason"], str(row["transition_reason"])
        )
        row["retention_failure_name"] = decode_retention_failure(
            row["retention_failure_mask"]
        )
        rupture_rows.append(row)

    transition_fields = COLUMNS + (
        "old_state_name", "new_state_name", "transition_reason_name",
    )
    rupture_fields = RUPTURE_COLUMNS + (
        "old_state_name", "transition_reason_name", "retention_failure_name",
    )
    dwell_fields = (
        "actin_id", "myosin_id", "head_id", "state", "start_step", "end_step",
        "dwell_steps", "dwell_time", "exit_state", "exit_reason", "entry_cb_status",
        "exit_cb_status", "entry_force", "exit_force", "entry_activation",
        "exit_activation",
    )
    summary_fields = (
        "state", "entry_cb_status", "exit_reason", "n_completed", "mean_dwell",
        "median_dwell", "std_dwell", "p10_dwell", "p90_dwell",
        "event_rate_estimate",
    )
    dwells = build_dwells(
        events, args.dt,
        read_last_step(args.trajectory, include_resume=not args.saved_only),
    )
    summaries = summarize_dwells(dwells)
    write_rows(Path(f"{prefix}_am_transitions.csv"), transition_fields, transition_rows)
    write_rows(
        Path(f"{prefix}_am_forced_ruptures.csv"), rupture_fields, rupture_rows
    )
    write_rows(Path(f"{prefix}_am_dwells.csv"), dwell_fields, dwells)
    write_rows(Path(f"{prefix}_am_dwell_summary.csv"), summary_fields, summaries)
    print(
        f"events={len(transition_rows)} forced_ruptures={len(rupture_rows)} "
        f"completed_dwells={len(dwells)} "
        f"groups={len(summaries)}"
    )


if __name__ == "__main__":
    main()
