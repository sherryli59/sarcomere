#!/usr/bin/env python3
"""Compare historical and current Cell A trajectories for exact compatibility."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np


CORE_PATHS = (
    "state/current_step",
    "actin/center",
    "actin/direction",
    "actin/velocity",
    "actin/force",
    "actin/torque",
    "actin/cb_status",
    "actin/aa_local_status1_count",
    "actin/aa_local_status2_count",
    "actin/am_weak_count",
    "actin/am_strong_count",
    "actin/am_q",
    "myosin/center",
    "myosin/direction",
    "myosin/velocity",
    "myosin/force",
    "myosin/torque",
    "state/am_head_state_current",
    "state/am_head_bound_actin_current",
    "state/am_head_target_axial_current",
    "state/am_head_binding_azimuth_current",
    "state/am_head_stroke_progress_current",
    "state/titin_bound_actin_current",
    "state/titin_force_current",
)

# These datasets exist only when the sparse finite-valency myomesin treatment
# is active.  If either trajectory contains one, require the other trajectory
# to contain and exactly match it.  This lets the same harness prove both the
# historical k_mm=0 path and preservation of the active treatment path.
OPTIONAL_TREATMENT_PATHS = (
    "myomesin/bridge_state",
    "myomesin/events",
    "myomesin/frame_offsets",
    "myomesin/next_id_by_frame",
)


def compare(old_path: Path, new_path: Path) -> list[dict]:
    results: list[dict] = []
    with (
        h5py.File(old_path, "r", locking=False) as old,
        h5py.File(new_path, "r", locking=False) as new,
    ):
        paths = list(CORE_PATHS)
        paths.extend(
            path for path in OPTIONAL_TREATMENT_PATHS
            if path in old or path in new
        )
        for path in paths:
            if path not in old or path not in new:
                results.append({
                    "path": path,
                    "equal": False,
                    "reason": f"missing old={path in old} new={path in new}",
                })
                continue
            first = old[path][:]
            second = new[path][:]
            if first.shape != second.shape:
                results.append({
                    "path": path,
                    "equal": False,
                    "reason": f"shape {first.shape} != {second.shape}",
                })
                continue
            unequal = first != second
            count = int(np.count_nonzero(unequal))
            record = {"path": path, "equal": count == 0, "unequal": count}
            if count:
                frame_axis = np.flatnonzero(
                    np.any(unequal.reshape((unequal.shape[0], -1)), axis=1)
                )
                record["first_frame"] = int(frame_axis[0])
                if np.issubdtype(first.dtype, np.number):
                    delta = np.abs(first.astype(float) - second.astype(float))
                    finite = np.isfinite(delta)
                    record["max_abs_difference"] = (
                        float(np.max(delta[finite])) if np.any(finite) else None
                    )
            results.append(record)
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("old", type=Path)
    parser.add_argument("new", type=Path)
    args = parser.parse_args()

    results = compare(args.old, args.new)
    failed = [record for record in results if not record["equal"]]
    for record in results:
        if record["equal"]:
            print(f"PASS {record['path']}")
        else:
            details = ", ".join(
                f"{key}={value}"
                for key, value in record.items()
                if key not in {"path", "equal"}
            )
            print(f"FAIL {record['path']}: {details}")
    print(
        f"summary: {len(results) - len(failed)}/{len(results)} exact; "
        f"{len(failed)} mismatched"
    )
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
