#!/usr/bin/env python3
"""Timing across a campaign's runs: wall clock, section breakdown, ACTN2 phases.

Reads the last profile line of each run log, so it works on a campaign that is
still in flight -- the numbers are cumulative to whatever step the run reached.
"""
from __future__ import annotations
import argparse, re
from pathlib import Path

def last_line(path: Path, prefix: str) -> str | None:
    found = None
    with path.open(errors="replace") as fh:
        for line in fh:
            if line.startswith(prefix):
                found = line
    return found

def kv(line: str, suffix: str) -> dict[str, float]:
    return {m.group(1): float(m.group(2))
            for m in re.finditer(rf"(\w+?){suffix}=([0-9.eE+-]+)", line)}

ap = argparse.ArgumentParser()
ap.add_argument("--root", type=Path, required=True)
ap.add_argument("--top", type=int, default=8)
ap.add_argument("--arms", default="",
                help="comma-separated arm names in array order, e.g. "
                     "'fullon,intermediate,hcm'; the log name carries only the "
                     "array index, so this is how a log maps to an arm")
ap.add_argument("--seeds-per-arm", type=int, default=3)
args = ap.parse_args()

logs = sorted((args.root / "logs").glob("*.[lo][ou][gt]"))
runs = []
for log in logs:
    rp = last_line(log, "runtime_profile")
    if not rp:
        continue
    prof = kv(rp, "_wall_s")
    steps_m = re.search(r"steps=(\d+)", rp)
    steps = int(steps_m.group(1)) if steps_m else 0
    sect_line = last_line(log, "parallel_section_profile:") or ""
    threads_m = re.search(r"threads=(\d+)", sect_line)
    sect = kv(sect_line, "_mean_s")
    phase = kv(last_line(log, "aa_phase_profile:") or "", "_mean_s")
    idx_m = re.search(r"_(\d+)$", log.stem)
    idx = int(idx_m.group(1)) if idx_m else -1
    arm_names = [a for a in args.arms.split(",") if a]
    arm = (arm_names[idx // args.seeds_per_arm]
           if arm_names and 0 <= idx // args.seeds_per_arm < len(arm_names)
           else log.stem)
    runs.append({"tag": log.stem, "arm": arm, "steps": steps, "prof": prof,
                 "sect": sect, "phase": phase,
                 "threads": int(threads_m.group(1)) if threads_m else 0})

if not runs:
    raise SystemExit(f"no runtime_profile lines under {args.root}/logs")

threads = {r["threads"] for r in runs if r["threads"]}
print(f"threads: {sorted(threads)}\n")
print(f"{'run':44s}{'steps':>8}{'elapsed_s':>11}{'ms/step':>9}"
      f"{'update':>9}{'parallel':>10}{'neighbor':>10}{'sample':>9}")
for r in sorted(runs, key=lambda r: r["tag"]):
    p, n = r["prof"], max(r["steps"], 1)
    print(f"{(r['arm']+' '+r['tag'].split('_')[-1])[:44]:44s}{r['steps']:>8}{p.get('elapsed',0):>11.1f}"
          f"{1000*p.get('elapsed',0)/n:>9.2f}"
          f"{p.get('update_system',0):>9.1f}{p.get('parallel_update',0):>10.1f}"
          f"{p.get('neighbor',0):>10.1f}{p.get('sample',0):>9.1f}")

def aggregate(key: str, label: str) -> None:
    keys: dict[str, float] = {}
    for r in runs:
        for k, v in r[key].items():
            keys[k] = keys.get(k, 0.0) + v
    if not keys:
        return
    total = sum(keys.values())
    print(f"\n{label} (summed over {len(runs)} runs, mean over threads)")
    for k, v in sorted(keys.items(), key=lambda kv: -kv[1])[:args.top]:
        print(f"  {k:26s}{v:10.1f} s{100*v/total:>7.1f}%")

aggregate("sect", "parallel_update sections")
aggregate("phase", "actn2_chemistry phases")

# per-arm normalised cost, so arms can be compared at equal step count
print("\nper-arm ms/step (elapsed)")
by_arm: dict[str, list[float]] = {}
for r in runs:
    arm = r["arm"]
    if r["steps"]:
        by_arm.setdefault(arm, []).append(1000 * r["prof"].get("elapsed", 0) / r["steps"])
for arm, vals in sorted(by_arm.items()):
    mean = sum(vals) / len(vals)
    print(f"  {arm:20s}{mean:8.2f}  (n={len(vals)}, "
          f"{', '.join(f'{v:.2f}' for v in vals)})")
