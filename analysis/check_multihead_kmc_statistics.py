#!/usr/bin/env python3
"""Validate independent explicit-head KMC hazards in fixed AM geometry."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    parser.add_argument("--dt", type=float, default=1e-5)
    parser.add_argument("--k-on", type=float, default=100.0)
    parser.add_argument("--k-weak-off", type=float, default=200.0)
    parser.add_argument("--k-powerstroke", type=float, default=200.0)
    parser.add_argument("--k-strong-off", type=float, default=80.0)
    parser.add_argument("--force-scale", type=float, default=1.0)
    args = parser.parse_args()

    exposure = np.zeros(3, dtype=float)
    transitions: defaultdict[tuple[int, int, int], int] = defaultdict(int)
    positive_gates: list[float] = []
    strong_exit_forces: list[float] = []
    files = sorted(args.directory.glob("kmc_seed*.h5"))
    if not files:
        parser.error("no kmc_seed*.h5 files")

    for path in files:
        with h5py.File(path, "r", locking=False) as handle:
            events = np.asarray(handle["/actin_myo/head_transitions"], dtype=float)
            nheads = int(handle["/meta"].attrs["am_heads_per_pair"])
            final_step = int(np.asarray(handle["/resume/current_step"]).reshape(-1)[0])
            gate = np.asarray(handle["/actin/partial_binding_ratio"], dtype=float).reshape(-1)
            positive_gates.extend(gate[gate > 0].tolist())

        by_head: defaultdict[tuple[int, int, int], list[np.ndarray]] = defaultdict(list)
        for event in events:
            key = (int(round(event[0])), int(round(event[1])), int(round(event[10])))
            by_head[key].append(event)
        # This diagnostic has exactly one geometrically allowed AM pair.
        for head in range(nheads):
            state = 0
            last_step = 0
            for event in sorted(by_head[(0, 0, head)], key=lambda row: row[2]):
                step = int(round(event[2]))
                old_state = int(round(event[3]))
                new_state = int(round(event[4]))
                reason = int(round(event[5]))
                if old_state != state:
                    raise RuntimeError(
                        f"state discontinuity in {path.name}, head {head}, step {step}: "
                        f"tracked {state}, event old {old_state}"
                    )
                exposure[state] += max(0, step - last_step) * args.dt
                transitions[(old_state, new_state, reason)] += 1
                if old_state == 2 and new_state == 0 and reason == 3:
                    strong_exit_forces.append(float(event[8]))
                state = new_state
                last_step = step
            exposure[state] += max(0, final_step - last_step) * args.dt

    attach = transitions[(0, 1, 0)]
    weak_off = transitions[(1, 0, 1)]
    powerstroke = transitions[(1, 2, 2)]
    strong_off = transitions[(2, 0, 3)]
    gate = float(np.median(positive_gates))
    mean_force = float(np.mean(strong_exit_forces))
    rows = []

    def hazard_row(name: str, count: int, state: int, expected: float) -> None:
        measured = count / exposure[state]
        rows.append({
            "test": name,
            "events": count,
            "exposure": exposure[state],
            "measured": measured,
            "expected": expected,
            "ratio": measured / expected,
            "poisson_relative_se": 1.0 / math.sqrt(count),
        })

    hazard_row("attachment", attach, 0, args.k_on * gate)
    hazard_row("weak_total_exit", weak_off + powerstroke, 1,
               args.k_weak_off + args.k_powerstroke)
    measured_branch = powerstroke / (weak_off + powerstroke)
    expected_branch = args.k_powerstroke / (args.k_weak_off + args.k_powerstroke)
    rows.append({
        "test": "powerstroke_branch",
        "events": weak_off + powerstroke,
        "exposure": np.nan,
        "measured": measured_branch,
        "expected": expected_branch,
        "ratio": measured_branch / expected_branch,
        "poisson_relative_se": np.nan,
    })
    hazard_row("strong_detachment", strong_off, 2,
               args.k_strong_off * math.exp(-mean_force / args.force_scale))

    output = args.directory / "multihead_kmc_statistics.csv"
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"files={len(files)} heads/file={nheads} gate={gate:.8g} mean_strong_force={mean_force:.8g}")
    for row in rows:
        print(f"{row['test']:22s} measured={row['measured']:.5g} "
              f"expected={row['expected']:.5g} ratio={row['ratio']:.4f} n={row['events']}")


if __name__ == "__main__":
    main()
