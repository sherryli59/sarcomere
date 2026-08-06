#!/usr/bin/env python3
"""Plot the three minimal motor-validation scenarios from diagnostic CSVs."""

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_csv(path: Path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    return {
        key: np.asarray([float(row[key]) for row in rows], dtype=float)
        for key in rows[0]
    }


def only(directory: Path, pattern: str) -> Path:
    matches = sorted(directory.glob(pattern))
    if len(matches) != 1:
        raise RuntimeError(f"expected one {pattern} in {directory}, found {len(matches)}")
    return matches[0]


def block_fv(data, width_s=0.002):
    t = data["time_s"]
    aa = data["aa_overlap_um"]
    load = data["actn2_load_pn_per_actin"]
    status2 = data["status2_links"]
    # Diagnostic CSV cadence is configurable. Ensure each F--V block spans at
    # least two sampled intervals instead of silently returning an empty plot
    # when a long run is saved more coarsely than the historical 2-ms blocks.
    if len(t) > 1:
        width_s = max(width_s, 2.1 * float(np.median(np.diff(t))))
    edges = np.arange(t[0], t[-1] + width_s, width_s)
    force, velocity, state = [], [], []
    for left, right in zip(edges[:-1], edges[1:]):
        mask = (t >= left) & (t <= right)
        indices = np.flatnonzero(mask)
        if len(indices) < 2:
            continue
        first, last = indices[0], indices[-1]
        force.append(np.mean(load[indices]))
        velocity.append(-(aa[last] - aa[first]) / (t[last] - t[first]))
        state.append(np.mean(status2[indices]))
    return np.asarray(force), np.asarray(velocity), np.asarray(state)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("suite_dir", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    free = read_csv(only(args.suite_dir, "free_*.csv"))
    linked = read_csv(only(args.suite_dir, "linked_*.csv"))
    balanced = read_csv(only(args.suite_dir, "balanced_*.csv"))
    bf, bv, bs = block_fv(balanced)

    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)

    ax = axes[0, 0]
    t_ms = 1e3 * balanced["time_s"]
    ax.plot(t_ms, balanced["aa_overlap_um"], color="black", lw=2, label="A–A")
    for index, color in enumerate(("#0072B2", "#56B4E9", "#D55E00", "#E69F00")):
        ax.plot(t_ms, balanced[f"am{index}_um"], color=color, lw=1.2,
                label=f"A{index}–owner M")
    ax.axhline(0.10, color="black", ls="--", lw=0.8)
    ax.set(xlabel="time (ms)", ylabel="axial overlap (µm)",
           title="Balanced cell: A–A register and A–M overlaps remain bounded")
    ax.legend(frameon=False, fontsize=8, ncol=2)

    ax = axes[0, 1]
    points = ax.scatter(bf, bv, c=bs, cmap="viridis", s=34,
                        edgecolor="none")
    ax.axhline(0.0, color="0.4", lw=0.8)
    ax.set(xlabel="measured ACTN2 load (pN/actin)",
           ylabel="A–A extension velocity (µm/s)",
           title="Endogenous force–velocity samples (2 ms windows)")
    fig.colorbar(points, ax=ax, label="mean status-2 links")

    ax = axes[1, 0]
    ax.plot(1e3 * free["time_s"], free["q0_um"], lw=2,
            color="#0072B2", label="free actin")
    ax.plot(1e3 * linked["time_s"], linked["q0_um"], lw=1.5,
            color="#D55E00", label="linked, motor-bound actin")
    ax.plot(1e3 * linked["time_s"], -linked["q1_um"], lw=1.5,
            color="#E69F00", label="linked, unbound follower")
    ax.set(xlabel="time (ms)", ylabel="directed travel (µm)",
           title="Free transport and one-sided follower transport")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 1]
    ax.plot(t_ms, balanced["actn2_load_pn_per_actin"], color="#CC79A7",
            lw=1.5, label="ACTN2 load")
    ax.set(xlabel="time (ms)", ylabel="ACTN2 load (pN/actin)",
           title="Load-bearing catch-state formation")
    state_axis = ax.twinx()
    state_axis.step(t_ms, balanced["status2_links"], where="mid",
                    color="#009E73", lw=1.0, label="status 2")
    state_axis.set_ylabel("status-2 molecules")
    lines = ax.lines + state_axis.lines
    ax.legend(lines, [line.get_label() for line in lines], frameon=False,
              fontsize=8, loc="upper right")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)


if __name__ == "__main__":
    main()
