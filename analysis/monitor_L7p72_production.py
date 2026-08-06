#!/usr/bin/env python3
"""Read-only health/timing snapshot for the A1280/M640 production run."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from datetime import datetime, timezone

import h5py
import numpy as np


def scalar_sum(handle: h5py.File, path: str, frame: int) -> float:
    if path not in handle:
        return math.nan
    return float(np.asarray(handle[path][frame], dtype=float).sum())


def nematic(directions: np.ndarray) -> float:
    u = np.asarray(directions, dtype=float).reshape(-1, 3)
    norms = np.linalg.norm(u, axis=1)
    u = u[norms > 0] / norms[norms > 0, None]
    if len(u) == 0:
        return math.nan
    q = 1.5 * np.einsum("ni,nj->ij", u, u) / len(u) - 0.5 * np.eye(3)
    return float(np.linalg.eigvalsh(q)[-1])


def append_csv(path: str, row: dict[str, object]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    exists = os.path.exists(path)
    with open(path, "a", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(row))
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("trajectory")
    parser.add_argument("--csv", required=True)
    parser.add_argument("--job-id", default="")
    parser.add_argument("--elapsed-seconds", type=float, default=math.nan)
    parser.add_argument("--target-steps", type=int, default=1_000_000)
    parser.add_argument("--ramp-steps", type=int, default=20_000)
    args = parser.parse_args()

    now = datetime.now(timezone.utc).isoformat()
    try:
        with h5py.File(args.trajectory, "r", locking=False) as h:
            nframes = int(h["/actin/center"].shape[0])
            if nframes == 0:
                raise RuntimeError("trajectory contains no committed frames")
            last = nframes - 1
            steps = np.asarray(h["/state/current_step"][:, 0], dtype=np.int64)
            step = int(steps[last])
            centers_a = np.asarray(h["/actin/center"][last], dtype=float)
            centers_m = np.asarray(h["/myosin/center"][last], dtype=float)
            dirs_a = np.asarray(h["/actin/direction"][last], dtype=float)
            dirs_m = np.asarray(h["/myosin/direction"][last], dtype=float)
            finite = bool(
                np.isfinite(centers_a).all()
                and np.isfinite(centers_m).all()
                and np.isfinite(dirs_a).all()
                and np.isfinite(dirs_m).all()
            )

            force_paths = (
                "/actin/cumulative_force_cap_events",
                "/myosin/cumulative_force_cap_events",
            )
            torque_paths = (
                "/actin/cumulative_torque_cap_events",
                "/myosin/cumulative_torque_cap_events",
            )
            force_caps = sum(scalar_sum(h, path, last) for path in force_paths)
            torque_caps = sum(scalar_sum(h, path, last) for path in torque_paths)
            weak = scalar_sum(h, "/actin/am_weak_count", last)
            strong = scalar_sum(h, "/actin/am_strong_count", last)
            status1 = 0.5 * scalar_sum(h, "/actin/aa_local_status1_count", last)
            status2 = 0.5 * scalar_sum(h, "/actin/aa_local_status2_count", last)

            interval_force_fractions: list[float] = []
            interval_torque_fractions: list[float] = []
            first = max(1, nframes - 4)
            for frame in range(first, nframes):
                delta_step = int(steps[frame] - steps[frame - 1])
                if delta_step <= 0:
                    continue
                denom = 1920.0 * delta_step
                f0 = sum(scalar_sum(h, path, frame - 1) for path in force_paths)
                f1 = sum(scalar_sum(h, path, frame) for path in force_paths)
                t0 = sum(scalar_sum(h, path, frame - 1) for path in torque_paths)
                t1 = sum(scalar_sum(h, path, frame) for path in torque_paths)
                interval_force_fractions.append((f1 - f0) / denom)
                interval_torque_fractions.append((t1 - t0) / denom)

            latest_force_fraction = (
                interval_force_fractions[-1] if interval_force_fractions else math.nan
            )
            latest_torque_fraction = (
                interval_torque_fractions[-1] if interval_torque_fractions else math.nan
            )
            accelerating = (
                len(interval_force_fractions) >= 3
                and interval_force_fractions[-3]
                < interval_force_fractions[-2]
                < interval_force_fractions[-1]
            )
            runaway_caps = bool(
                step >= args.ramp_steps + 3_000
                and accelerating
                and latest_force_fraction > 1.0e-3
            )
            rate = (
                step / args.elapsed_seconds
                if step > 0
                and math.isfinite(args.elapsed_seconds)
                and args.elapsed_seconds > 0
                else math.nan
            )
            eta_hours = (
                max(0, args.target_steps - step) / rate / 3600
                if math.isfinite(rate) and rate > 0
                else math.nan
            )
            row = {
                "observed_utc": now,
                "job_id": args.job_id,
                "elapsed_seconds": args.elapsed_seconds,
                "frames": nframes,
                "step": step,
                "sim_time_s": step * 1.0e-6,
                "steps_per_wall_s": rate,
                "eta_hours": eta_hours,
                "force_caps": force_caps,
                "torque_caps": torque_caps,
                "interval_force_cap_fraction": latest_force_fraction,
                "interval_torque_cap_fraction": latest_torque_fraction,
                "force_caps_accelerating": int(accelerating),
                "weak_heads": weak,
                "strong_heads": strong,
                "status1_links": status1,
                "status2_links": status2,
                "actin_nematic": nematic(dirs_a),
                "myosin_nematic": nematic(dirs_m),
                "finite": int(finite),
                "runaway_caps": int(runaway_caps),
            }
    except (OSError, KeyError, RuntimeError, ValueError) as error:
        print(json.dumps({"observed_utc": now, "read_error": str(error)}))
        return 0

    append_csv(args.csv, row)
    print(json.dumps(row, allow_nan=True))
    if not finite:
        return 42
    if runaway_caps:
        return 43
    return 0


if __name__ == "__main__":
    sys.exit(main())
