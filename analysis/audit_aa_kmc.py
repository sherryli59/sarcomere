#!/usr/bin/env python3
"""Audit ACTN2 state exposure and all-cause exits from transition events."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import h5py
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("trajectory", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--dt", type=float, default=1e-5)
    parser.add_argument("--base-lifetime", type=float, default=1e-3)
    parser.add_argument("--lifetime-coeff", type=float, default=1e-2)
    args = parser.parse_args()

    with h5py.File(args.trajectory, "r", locking=False) as handle:
        events = np.asarray(handle["/catch_bond/state_transitions"])
        final_step = int(np.asarray(handle["/resume/current_step"]).reshape(-1)[0])
    events = events[np.argsort(events[:, 2], kind="stable")]

    active: dict[tuple[int, int], tuple[int, int]] = {}
    exposure = {1: 0.0, 2: 0.0}
    exits: dict[tuple[int, int], int] = {}
    endpoint_loads: dict[tuple[int, int], list[float]] = {}
    mismatches = 0
    for event in events:
        i, j, step, _, _, old_status, new_status = map(int, np.rint(event[:7]))
        key = (min(i, j), max(i, j))
        if old_status in (1, 2):
            entry = active.get(key)
            if entry is None or entry[0] != old_status:
                mismatches += 1
            else:
                exposure[old_status] += (step - entry[1]) * args.dt
        if old_status in (1, 2) and new_status != old_status:
            transition = (old_status, new_status)
            exits[transition] = exits.get(transition, 0) + 1
            endpoint_loads.setdefault(transition, []).append(
                max(0.0, min(float(event[11]), float(event[12])))
            )
        if new_status in (1, 2):
            active[key] = (new_status, step)
        else:
            active.pop(key, None)
    for status, start_step in active.values():
        exposure[status] += (final_step - start_step) * args.dt

    rows = []
    for status in (1, 2):
        n_detach = exits.get((status, 0), 0)
        n_switch = sum(n for (old, new), n in exits.items() if old == status and new != 0)
        all_cause_rate = n_detach / exposure[status]
        if status == 1:
            prescribed = 1.0 / args.base_lifetime
            note = "all-cause detach includes geometric loss; prescribed KMC rate is a lower bound"
            endpoint_load = np.nan
        else:
            loads = np.asarray(endpoint_loads.get((2, 0), []), dtype=float)
            prescribed = float(np.mean(1.0 / (args.base_lifetime + args.lifetime_coeff * loads)))
            endpoint_load = float(np.mean(loads))
            note = "prescribed estimate uses exit-time load; geometric loss is included in measured exits"
        rows.append({
            "status": status, "exposure": exposure[status], "n_detach": n_detach,
            "n_status_switch": n_switch, "all_cause_detach_rate": all_cause_rate,
            "prescribed_kmc_rate_or_endpoint_estimate": prescribed,
            "rate_ratio": all_cause_rate / prescribed,
            "mean_load_at_detach": endpoint_load, "note": note,
        })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"events={len(events)} mismatches={mismatches} active_final={len(active)} output={args.output}")
    for row in rows:
        print(
            f"status={row['status']} exposure={row['exposure']:.6g} "
            f"detach={row['n_detach']} all_cause_rate={row['all_cause_detach_rate']:.6g} "
            f"prescribed={row['prescribed_kmc_rate_or_endpoint_estimate']:.6g}"
        )


if __name__ == "__main__":
    main()
