#!/usr/bin/env python3
"""Audit fixed-population actin renewal events and realized Poisson rates."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import h5py
import numpy as np


def analyze(path: Path) -> dict[str, object]:
    with h5py.File(path, "r", locking=False) as handle:
        meta = handle["meta"].attrs
        rate = float(meta.get("actin_turnover_rate", 0.0))
        start = float(meta.get("actin_turnover_start_time", 0.0))
        dt = float(meta["dt"])
        n_actins = int(handle["actin/center"].shape[1])
        final_step = int(handle["state/current_step"][-1, 0])
        active_time = max(0.0, final_step * dt - start)
        events = (
            np.asarray(handle["actin/turnover_events"], dtype=float)
            if "actin/turnover_events" in handle
            else np.empty((0, 20), dtype=float)
        )
        success = events[:, 3] > 0.5 if len(events) else np.zeros(0, bool)
        expected = n_actins * rate * active_time
        observed = int(np.count_nonzero(success))
        successful_actins = (
            np.unique(events[success, 0].astype(int)) if observed else np.empty(0, int)
        )
        variance = max(expected, 1.0)
        z_score = (observed - expected) / math.sqrt(variance)
        centers_finite = bool(np.isfinite(handle["actin/center"][:]).all())
        directions = np.asarray(handle["actin/direction"][:], dtype=float)
        direction_norm_error = float(
            np.max(np.abs(np.linalg.norm(directions, axis=2) - 1.0))
        )
        cumulative = np.asarray(
            handle["actin/cumulative_turnover_events"][-1], dtype=float
        ).reshape(-1)
    return {
        "path": str(path),
        "n_actins": n_actins,
        "rate_per_s": rate,
        "active_time_s": active_time,
        "expected_events": expected,
        "successful_events": observed,
        "failed_placements": int(len(events) - observed),
        "poisson_z_score": z_score,
        "realized_events_per_actin": observed / max(1, n_actins),
        "realized_unique_fraction": len(successful_actins) / max(1, n_actins),
        "expected_unique_fraction": 1.0 - math.exp(-rate * active_time),
        "mean_placement_attempts": float(np.mean(events[:, 4])) if len(events) else 0.0,
        "max_placement_attempts": int(np.max(events[:, 4])) if len(events) else 0,
        "released_am_heads": int(np.sum(events[:, 5])) if len(events) else 0,
        "released_actn2": int(np.sum(events[:, 6])) if len(events) else 0,
        "released_titin": int(np.sum(events[:, 7])) if len(events) else 0,
        "cumulative_feature_events": int(np.rint(cumulative.sum())),
        "finite_centers": centers_finite,
        "max_direction_norm_error": direction_norm_error,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = [analyze(path) for path in args.inputs]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        if row["successful_events"] != row["cumulative_feature_events"]:
            raise SystemExit(f"turnover count mismatch: {row['path']}")
        if not row["finite_centers"] or row["max_direction_norm_error"] > 1e-10:
            raise SystemExit(f"invalid renewed geometry: {row['path']}")
    print(args.output)


if __name__ == "__main__":
    main()

