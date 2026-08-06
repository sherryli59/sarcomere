#!/usr/bin/env python3
"""Summarize order/occupancy responses in the lambda history-switch experiment.

Input is the per-frame CSV produced by ``extract_sarcomeric_order.py`` for the
history-switch HDF5 directory.  Each protocol is split into its kinetic phases;
the tail of every phase is averaged to avoid scoring the instantaneous switch
transient as a plateau.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

import numpy as np


HISTORY_RE = re.compile(
    r"_history_(medium_control|high_control|medium_high|high_medium_high)_(r\d+)_L"
)


def read_records(path: Path) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    with path.open() as stream:
        for row in csv.DictReader(stream):
            match = HISTORY_RE.search(row["filename"])
            if not match:
                continue
            records.append(
                {
                    "filename": row["filename"],
                    "protocol": match.group(1),
                    "replica": match.group(2),
                    "step": int(float(row["step"])),
                    "q_sarcomere": float(row["q_sarcomere"]),
                    "joint_band_amplitude": float(row["shared_joint_band_amplitude"]),
                    "n_actin_cb2": float(row["n_actin_cb2"]),
                }
            )
    if not records:
        raise ValueError(f"no history-switch rows found in {path}")
    return records


def phase_ranges(
    protocol: str,
    total_steps: int,
    medium_high_switch: int,
    pulse_start: int,
    pulse_end: int,
) -> list[tuple[str, int, int, str]]:
    if protocol == "medium_control":
        return [("medium_control", 0, total_steps, "medium")]
    if protocol == "high_control":
        return [("high_control", 0, total_steps, "high")]
    if protocol == "medium_high":
        return [
            ("medium_assembly", 0, medium_high_switch, "medium"),
            ("high_hold", medium_high_switch, total_steps, "high"),
        ]
    if protocol == "high_medium_high":
        return [
            ("high_initial", 0, pulse_start, "high"),
            ("medium_rescue", pulse_start, pulse_end, "medium"),
            ("high_rehold", pulse_end, total_steps, "high"),
        ]
    raise ValueError(f"unknown protocol {protocol}")


def phase_tail(rows: list[dict[str, object]], start: int, end: int, tail_fraction: float):
    phase = [row for row in rows if start <= int(row["step"]) <= end]
    if not phase:
        return []
    cutoff = start + (1.0 - tail_fraction) * (end - start)
    tail = [row for row in phase if int(row["step"]) >= cutoff]
    return tail or [phase[-1]]


def summarize(
    records: list[dict[str, object]],
    total_steps: int,
    medium_high_switch: int,
    pulse_start: int,
    pulse_end: int,
    tail_fraction: float,
) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in records:
        grouped[(str(row["protocol"]), str(row["replica"]))].append(row)

    summaries: list[dict[str, object]] = []
    for (protocol, replica), rows in sorted(grouped.items()):
        rows.sort(key=lambda row: int(row["step"]))
        for phase, start, end, lambda_state in phase_ranges(
            protocol, total_steps, medium_high_switch, pulse_start, pulse_end
        ):
            tail = phase_tail(rows, start, end, tail_fraction)
            if not tail:
                continue
            summaries.append(
                {
                    "protocol": protocol,
                    "replica": replica,
                    "phase": phase,
                    "lambda_state": lambda_state,
                    "phase_start_step": start,
                    "phase_end_step": end,
                    "n_tail_samples": len(tail),
                    "q_sarcomere": float(np.mean([row["q_sarcomere"] for row in tail])),
                    "joint_band_amplitude": float(
                        np.mean([row["joint_band_amplitude"] for row in tail])
                    ),
                    "n_actin_cb2": float(np.mean([row["n_actin_cb2"] for row in tail])),
                }
            )
    return summaries


def phase_means(rows: list[dict[str, object]]) -> dict[str, float]:
    grouped: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        grouped[str(row["phase"])].append(float(row["q_sarcomere"]))
    return {phase: float(np.mean(values)) for phase, values in grouped.items()}


def contrasts(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    q = phase_means(rows)
    definitions = (
        (
            "high_accessibility_gap",
            "medium_control",
            "high_control",
            "Positive means high lambda cannot assemble as well from disorder.",
        ),
        (
            "high_preservation_gap",
            "high_hold",
            "high_control",
            "Positive means high lambda preserves medium-assembled order better than it creates it.",
        ),
        (
            "medium_rescue_gain",
            "medium_rescue",
            "high_initial",
            "Positive means the medium-lambda pulse rescues a high-lambda arrested state.",
        ),
        (
            "post_rescue_retention",
            "high_rehold",
            "medium_rescue",
            "Near zero means restored high lambda retains the rescued order.",
        ),
    )
    output: list[dict[str, object]] = []
    for name, positive, negative, interpretation in definitions:
        if positive not in q or negative not in q:
            continue
        output.append(
            {
                "contrast": name,
                "positive_phase": positive,
                "negative_phase": negative,
                "q_difference": q[positive] - q[negative],
                "interpretation": interpretation,
            }
        )
    return output


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"no rows to write to {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("order_csv", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--total-steps", type=int, default=400000)
    parser.add_argument("--medium-high-switch", type=int, default=200000)
    parser.add_argument("--pulse-start", type=int, default=100000)
    parser.add_argument("--pulse-end", type=int, default=300000)
    parser.add_argument("--phase-tail-fraction", type=float, default=0.25)
    args = parser.parse_args()

    if not 0.0 < args.phase_tail_fraction <= 1.0:
        parser.error("--phase-tail-fraction must lie in (0, 1]")
    records = read_records(args.order_csv)
    summaries = summarize(
        records,
        args.total_steps,
        args.medium_high_switch,
        args.pulse_start,
        args.pulse_end,
        args.phase_tail_fraction,
    )
    response = contrasts(summaries)
    write_csv(args.output_dir / "history_phase_summary.csv", summaries)
    write_csv(args.output_dir / "history_contrasts.csv", response)
    print(f"wrote history summaries to {args.output_dir}")
    for row in response:
        print(f"  {row['contrast']}: {float(row['q_difference']):+.4f}")


if __name__ == "__main__":
    main()
