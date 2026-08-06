#!/usr/bin/env python3
"""Aggregate fixed single-motor force-clamp trajectories."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

import numpy as np

from summarize_am_transitions import build_dwells, read_events, read_last_step


PATTERN = re.compile(r"am_single_F(?P<force>0|0p5|1|2)_r(?P<rep>\d+)\.h5$")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("--dt", type=float, default=1e-5)
    parser.add_argument("--k-strong-off", type=float, default=80.0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = args.output or args.input_dir / "single_motor_force_clamp_summary.csv"

    pooled: defaultdict[float, list[float]] = defaultdict(list)
    exposure: defaultdict[float, float] = defaultdict(float)
    n_files: defaultdict[float, int] = defaultdict(int)
    for path in sorted(args.input_dir.glob("am_single_F*_r*.h5")):
        match = PATTERN.search(path.name)
        if match is None:
            continue
        force_ratio = float(match.group("force").replace("p", "."))
        dwells = build_dwells(read_events(path), args.dt, read_last_step(path))
        strong_kinetic = [
            float(row["dwell_time"]) for row in dwells
            if row["state"] == "strong" and row["exit_reason"] == "strong_detachment"
            and float(row["dwell_time"]) > 0.0
        ]
        strong_exposure = sum(
            float(row["dwell_time"]) for row in dwells if row["state"] == "strong"
        )
        pooled[force_ratio].extend(strong_kinetic)
        exposure[force_ratio] += strong_exposure
        n_files[force_ratio] += 1

    fields = (
        "force_ratio", "n_files", "n_strong_detachments", "mean_strong_dwell",
        "sem_strong_dwell", "measured_detachment_rate", "expected_detachment_rate",
        "rate_ratio_measured_expected",
    )
    rows = []
    for force_ratio in sorted(pooled):
        values = np.asarray(pooled[force_ratio], dtype=float)
        total_exposure = exposure[force_ratio]
        measured_rate = values.size / total_exposure if total_exposure > 0.0 else np.nan
        mean = 1.0 / measured_rate if measured_rate > 0.0 else np.nan
        expected_rate = args.k_strong_off * np.exp(-force_ratio)
        rows.append({
            "force_ratio": force_ratio,
            "n_files": n_files[force_ratio],
            "n_strong_detachments": int(values.size),
            "mean_strong_dwell": mean,
            "sem_strong_dwell": mean / np.sqrt(values.size) if values.size else np.nan,
            "measured_detachment_rate": measured_rate,
            "expected_detachment_rate": expected_rate,
            "rate_ratio_measured_expected": measured_rate / expected_rate if expected_rate > 0 else np.nan,
        })
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"conditions={len(rows)} output={output}")


if __name__ == "__main__":
    main()
