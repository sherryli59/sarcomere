#!/usr/bin/env python3
"""Check AM KMC event statistics against the configured continuous-time rates."""

from __future__ import annotations

import argparse
import csv
import glob
import math
import re
import sys
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from summarize_am_transitions import build_dwells, read_events, read_last_step


FORCE_RE = re.compile(r"am_single_F(0|0p5|1|2)_r\d+\.h5$")


def force_ratio(path: str) -> float:
    match = FORCE_RE.search(path)
    if match is None:
        raise ValueError(path)
    return float(match.group(1).replace("p", "."))


def detached_exposure(events: np.ndarray, final_step: int, dt: float) -> tuple[int, float]:
    start = 0
    exposure = 0.0
    attachments = 0
    detached = True
    for event in events:
        step, old_state, new_state = map(int, np.rint(event[[2, 3, 4]]))
        if old_state == 0 and new_state == 1:
            if detached:
                exposure += max(0, step - start) * dt
            attachments += 1
            detached = False
        elif new_state == 0:
            start = step
            detached = True
    if detached:
        exposure += max(0, final_step - start) * dt
    return attachments, exposure


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--dt", type=float, default=1e-5)
    parser.add_argument("--k-on", type=float, default=1000.0)
    parser.add_argument("--k-weak-off", type=float, default=200.0)
    parser.add_argument("--k-powerstroke", type=float, default=200.0)
    parser.add_argument("--k-strong-off", type=float, default=80.0)
    args = parser.parse_args()

    files = sorted(glob.glob(str(args.input_dir / "am_single_F*_r*.h5")))
    if not files:
        parser.error("no force-clamp trajectories found")

    events_by_force: dict[float, list[np.ndarray]] = {}
    dwells_by_force: dict[float, list[dict]] = {}
    attach_count = 0
    attach_exposure = 0.0
    positive_gates: list[float] = []
    for filename in files:
        ratio = force_ratio(filename)
        events = read_events(Path(filename))
        final_step = read_last_step(Path(filename))
        events_by_force.setdefault(ratio, []).append(events)
        dwells_by_force.setdefault(ratio, []).extend(
            build_dwells(events, args.dt, final_step)
        )
        count, exposure = detached_exposure(events, final_step, args.dt)
        attach_count += count
        attach_exposure += exposure
        with h5py.File(filename, "r", locking=False) as handle:
            gate = np.asarray(handle["/actin/partial_binding_ratio"][:, 0, 0])
            positive_gates.extend(gate[gate > 0.0].tolist())

    rows: list[dict] = []
    weak_events = np.concatenate(
        [events for chunks in events_by_force.values() for events in chunks], axis=0
    )
    weak_events = weak_events[
        (weak_events[:, 3] == 1) & np.isin(weak_events[:, 5], [1, 2])
    ]
    all_dwells = [row for rows in dwells_by_force.values() for row in rows]
    weak_exposure = sum(row["dwell_time"] for row in all_dwells if row["state"] == "weak")
    weak_measured = len(weak_events) / weak_exposure
    weak_expected = args.k_weak_off + args.k_powerstroke
    rows.append({"test": "weak_total_exit", "condition": "all", "n_events": len(weak_events),
                 "exposure": weak_exposure, "measured": weak_measured,
                 "expected": weak_expected, "ratio": weak_measured / weak_expected})

    n_powerstroke = int(np.count_nonzero(weak_events[:, 5] == 2))
    branch_measured = n_powerstroke / len(weak_events)
    branch_expected = args.k_powerstroke / weak_expected
    rows.append({"test": "weak_powerstroke_branch", "condition": "all",
                 "n_events": len(weak_events), "exposure": math.nan,
                 "measured": branch_measured, "expected": branch_expected,
                 "ratio": branch_measured / branch_expected})

    gate = float(np.median(positive_gates))
    attach_measured = attach_count / attach_exposure
    attach_expected = args.k_on * gate
    rows.append({"test": "attachment", "condition": f"geometry_gate={gate:.8g}",
                 "n_events": attach_count, "exposure": attach_exposure,
                 "measured": attach_measured, "expected": attach_expected,
                 "ratio": attach_measured / attach_expected})

    strong_plot = []
    for ratio in sorted(events_by_force):
        events = np.concatenate(events_by_force[ratio], axis=0)
        dwells = dwells_by_force[ratio]
        n_detach = int(np.count_nonzero((events[:, 3] == 2) & (events[:, 5] == 3)))
        exposure = sum(row["dwell_time"] for row in dwells if row["state"] == "strong")
        measured = n_detach / exposure
        expected = args.k_strong_off * math.exp(-ratio)
        rows.append({"test": "strong_detachment", "condition": f"F/F0={ratio:g}",
                     "n_events": n_detach, "exposure": exposure, "measured": measured,
                     "expected": expected, "ratio": measured / expected})
        strong_plot.append((ratio, n_detach, measured, expected))

    prefix = args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    csv_path = Path(f"{prefix}.csv")
    fields = ("test", "condition", "n_events", "exposure", "measured", "expected", "ratio")
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    x = np.asarray([row[0] for row in strong_plot])
    count = np.asarray([row[1] for row in strong_plot])
    measured = np.asarray([row[2] for row in strong_plot])
    expected = np.asarray([row[3] for row in strong_plot])
    fig, axis = plt.subplots(figsize=(5.5, 4.0))
    axis.errorbar(x, measured, yerr=measured / np.sqrt(count), fmt="o", capsize=3,
                  color="#3b528b", label="measured (Poisson SE)")
    axis.plot(x, expected, "--", color="#e76f51", label=r"$80e^{-F/F_0}$")
    axis.set(xlabel=r"Resisting load $F/F_0$", ylabel="Strong detachment rate")
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(frameon=False)
    fig.tight_layout()
    plot_path = Path(f"{prefix}_strong_detachment.png")
    fig.savefig(plot_path, dpi=220, bbox_inches="tight")
    plt.close(fig)

    print(f"files={len(files)} csv={csv_path} plot={plot_path}")
    for row in rows:
        print(f"{row['test']:25s} {row['condition']:22s} ratio={row['ratio']:.4f} n={row['n_events']}")


if __name__ == "__main__":
    main()
