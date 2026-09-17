#!/usr/bin/env python3
"""Periodically append health records while a simulation process is alive."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import threading

from check_live_production_hdf5 import inspect


def process_exists(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def record(
    trajectory: Path,
    output: Path,
    actin_force_cap_fraction_max: float,
    myosin_force_cap_fraction_max: float,
) -> None:
    try:
        result = inspect(
            trajectory,
            actin_force_cap_fraction_max,
            myosin_force_cap_fraction_max,
        )
    except (OSError, RuntimeError, KeyError, ValueError) as error:
        result = {
            "trajectory": str(trajectory),
            "healthy": False,
            "transient_read_error": str(error),
        }
    result["checked_utc"] = datetime.now(timezone.utc).isoformat()
    with output.open("a") as stream:
        stream.write(json.dumps(result, separators=(",", ":")) + "\n")
        stream.flush()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trajectory", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--pid", type=int, required=True)
    parser.add_argument("--initial-seconds", type=float, default=60.0)
    parser.add_argument("--interval-seconds", type=float, default=300.0)
    parser.add_argument(
        "--actin-force-cap-fraction-max", type=float, default=0.05
    )
    parser.add_argument(
        "--myosin-force-cap-fraction-max", type=float, default=0.05
    )
    args = parser.parse_args()

    waiter = threading.Event()
    delay = max(0.0, args.initial_seconds)
    while process_exists(args.pid):
        waiter.wait(delay)
        if not process_exists(args.pid):
            break
        if args.trajectory.is_file() and args.trajectory.stat().st_size > 0:
            record(
                args.trajectory,
                args.output,
                args.actin_force_cap_fraction_max,
                args.myosin_force_cap_fraction_max,
            )
        delay = max(1.0, args.interval_seconds)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
