#!/usr/bin/env python3
"""Run the short, non-production full-simulator overlap diagnostic matrix."""

from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import os
from pathlib import Path
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parent
BINARY = ROOT / "build" / "axial_overlap_diagnostic"
OUT = ROOT / "results" / "campaign"
MECHANISMS = (
    "turnover",
    "titin",
    "register_titin",
    "opposing_restore_titin",
    "opposing_assist_titin",
    "bidirectional_titin",
)
INITIAL_CONDITIONS = (
    (4.3, 0.04),
    (4.3, 0.10),
    (4.3, 0.16),
    (3.9, 0.10),
    (4.7, 0.10),
)
SEEDS = (1001, 1002, 1003, 1004)
NSTEPS = 300_000
SAVE_EVERY = 1_500


def stem(mechanism, length, overlap, seed):
    return (
        f"{mechanism}__L{length:.1f}__AA{overlap:.2f}__s{seed}"
        .replace(".", "p")
    )


def run_one(task):
    mechanism, length, overlap, seed = task
    prefix = OUT / stem(mechanism, length, overlap, seed)
    command = [
        str(BINARY), mechanism, str(length), str(overlap), str(seed),
        str(NSTEPS), str(SAVE_EVERY), str(prefix),
    ]
    environment = dict(os.environ)
    environment.update(
        OMP_NUM_THREADS="1",
        OMP_DYNAMIC="FALSE",
        HDF5_USE_FILE_LOCKING="FALSE",
    )
    started = time.time()
    completed = subprocess.run(
        command, env=environment, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    )
    return {
        "task": task,
        "returncode": completed.returncode,
        "elapsed": time.time() - started,
        "output": completed.stdout[-2000:],
    }


def linear_slope(xs, ys):
    n = len(xs)
    if n < 2:
        return 0.0
    xm = sum(xs) / n
    ym = sum(ys) / n
    denominator = sum((x - xm) ** 2 for x in xs)
    if denominator == 0.0:
        return 0.0
    return sum((x - xm) * (y - ym) for x, y in zip(xs, ys)) / denominator


def summarize_one(task):
    mechanism, length, overlap, seed = task
    path = OUT / f"{stem(mechanism, length, overlap, seed)}.csv"
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    tail = rows[max(1, int(0.8 * len(rows))):]
    times = [float(row["time"]) for row in tail]

    def values(name):
        return [float(row[name]) for row in tail]

    result = {
        "mechanism": mechanism,
        "initial_Lx": length,
        "initial_AA": overlap,
        "seed": seed,
    }
    for name in ("Lx", "aa_overlap", "am_overlap"):
        data = values(name)
        result[f"tail_{name}"] = sum(data) / len(data)
        result[f"slope_{name}"] = linear_slope(times, data)
    for name in (
        "weak_heads", "strong_heads", "aa_status1", "aa_status2",
        "actin_tilt", "myosin_tilt", "q_total",
    ):
        data = values(name)
        result[f"tail_{name}"] = sum(data) / len(data)
    return result


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def aggregate(rows):
    summary = []
    for mechanism in MECHANISMS:
        selected = [row for row in rows if row["mechanism"] == mechanism]
        record = {"mechanism": mechanism}
        for name in ("Lx", "aa_overlap", "am_overlap"):
            values = [row[f"tail_{name}"] for row in selected]
            slopes = [abs(row[f"slope_{name}"]) for row in selected]
            record[f"mean_{name}"] = sum(values) / len(values)
            record[f"spread_{name}"] = max(values) - min(values)
            record[f"max_abs_slope_{name}"] = max(slopes)
        for name in (
            "weak_heads", "strong_heads", "aa_status1", "aa_status2",
            "actin_tilt", "myosin_tilt",
        ):
            values = [row[f"tail_{name}"] for row in selected]
            record[f"mean_{name}"] = sum(values) / len(values)
        record["target_AA_error"] = abs(record["mean_aa_overlap"] - 0.10)
        record["finite_AM"] = (
            0.35 <= record["mean_am_overlap"] <= 0.70
        )
        record["pass"] = (
            record["target_AA_error"] <= 0.025
            and record["spread_aa_overlap"] <= 0.05
            and record["spread_am_overlap"] <= 0.05
            and record["spread_Lx"] <= 0.10
            and record["max_abs_slope_aa_overlap"] <= 0.05
            and record["max_abs_slope_am_overlap"] <= 0.05
            and record["max_abs_slope_Lx"] <= 0.10
            and record["finite_AM"]
            and record["mean_strong_heads"] >= 1.0
            and record["mean_aa_status1"] + record["mean_aa_status2"] >= 1.0
        )
        summary.append(record)
    return summary


def main():
    if not BINARY.exists():
        raise SystemExit(f"missing diagnostic binary: {BINARY}")
    OUT.mkdir(parents=True, exist_ok=True)
    tasks = [
        (mechanism, length, overlap, seed)
        for mechanism in MECHANISMS
        for length, overlap in INITIAL_CONDITIONS
        for seed in SEEDS
    ]
    workers = min(4, os.cpu_count() or 1)
    failures = []
    started = time.time()
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(run_one, task): task for task in tasks}
        completed_count = 0
        for future in as_completed(futures):
            result = future.result()
            completed_count += 1
            if result["returncode"] != 0:
                failures.append(result)
            if completed_count % 10 == 0 or completed_count == len(tasks):
                print(
                    f"completed {completed_count}/{len(tasks)} "
                    f"in {time.time() - started:.1f}s",
                    flush=True,
                )
    if failures:
        for failure in failures:
            print(f"FAILED {failure['task']}\n{failure['output']}", file=sys.stderr)
        return 1
    rows = [summarize_one(task) for task in tasks]
    write_csv(OUT / "run_summary.csv", rows)
    write_csv(OUT / "mechanism_summary.csv", aggregate(rows))
    print(f"campaign complete in {time.time() - started:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
