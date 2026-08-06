#!/usr/bin/env python3
"""Return success only when the reduced sparse-ACTN2 pilot is healthy."""

import argparse
import csv
import math
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", type=Path)
    args = parser.parse_args()
    with args.csv.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    failures = []
    if len(rows) != 15:
        failures.append(f"expected 15 runs, found {len(rows)}")
    if any(int(float(row["finite"])) != 1 for row in rows):
        failures.append("non-finite trajectory")
    cap_fields = (
        "actin_force_cap_fraction", "actin_torque_cap_fraction",
        "myosin_force_cap_fraction", "myosin_torque_cap_fraction",
    )
    for field in cap_fields:
        worst = max((float(row[field]) for row in rows), default=math.inf)
        if worst >= 1e-4:
            failures.append(f"{field}={worst:.3g} >= 1e-4")
    enabled = [row for row in rows if row["condition"] != "aa_off"]
    if not enabled or any(int(float(row["associations"])) == 0 for row in enabled):
        failures.append("AA-enabled run without association")
    if not enabled or max(float(row["mean_strong_heads"]) for row in enabled) <= 0:
        failures.append("no strong AM occupancy")
    status2 = [row for row in rows if row["condition"] == "two_status2"]
    if sum(int(float(row["promotions"])) for row in status2) == 0:
        failures.append("no motor-load status-2 promotion")
    rates = [float(row["status1_off_rate"]) for row in enabled]
    rates = [value for value in rates if math.isfinite(value)]
    if not rates or not 50 <= sum(rates) / len(rates) <= 200:
        failures.append("status-1 off-rate outside broad 100/s acceptance band")
    if failures:
        print("REDUCED PILOT FAILED PROMOTION GATE")
        for failure in failures:
            print(f"- {failure}")
        raise SystemExit(1)
    print("REDUCED PILOT PASSED PROMOTION GATE")


if __name__ == "__main__":
    main()
