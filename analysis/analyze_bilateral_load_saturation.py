#!/usr/bin/env python3
"""Audit whether ACTN2 bilateral arm load approaches motor force balance."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/sarcomere_matplotlib")

import h5py
import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trajectory", type=Path)
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--plot", type=Path, required=True)
    args = parser.parse_args()

    with h5py.File(args.trajectory, "r", locking=False) as handle:
        dt = float(handle["/meta"].attrs["dt"])
        offsets = np.asarray(
            handle["/catch_bond/local_link_frame_offsets"], dtype=int)
        identity = np.asarray(
            handle["/catch_bond/local_link_identity"], dtype=int)
        mechanics = np.asarray(
            handle["/catch_bond/local_link_mechanics"], dtype=float)
        calcium = np.asarray(handle["/actin/am_calcium"], dtype=float)
        calcium = calcium.reshape(len(offsets), -1).mean(axis=1)
        state = np.asarray(handle["/state/am_head_state_current"], dtype=int)
        bound = np.asarray(
            handle["/state/am_head_bound_actin_current"], dtype=int)
        strong = np.count_nonzero((bound >= 0) & (state == 2), axis=(1, 2))

    time = offsets[:, 0] * dt
    bilateral = np.zeros(len(time))
    pair_bilateral = np.zeros(len(time))
    support = np.zeros(len(time))
    max_extension = np.zeros(len(time))
    raw_bilateral = np.zeros(len(time))
    allocated_bilateral = np.zeros(len(time))
    aa2 = np.zeros(len(time), dtype=int)
    for frame, (_, start, count) in enumerate(offsets):
        rows = mechanics[start:start + count]
        statuses = identity[start:start + count, 4]
        if count:
            bilateral[frame] = np.max(rows[:, 7])
            support[frame] = np.max(rows[:, 8])
            max_extension[frame] = np.max(rows[:, 4])
            aa2[frame] = np.count_nonzero(statuses == 2)
            if rows.shape[1] >= 16:
                raw_bilateral[frame] = np.max(
                    np.maximum(0.0, np.minimum(rows[:, 12], rows[:, 13])))
                allocated_bilateral[frame] = np.max(
                    np.maximum(0.0, np.minimum(rows[:, 14], rows[:, 15])))
            else:
                raw_bilateral[frame] = bilateral[frame]
                allocated_bilateral[frame] = bilateral[frame]
            pair_totals: dict[tuple[int, int], float] = {}
            for ident, row in zip(identity[start:start + count], rows):
                pair = tuple(sorted((int(ident[2]), int(ident[3]))))
                pair_totals[pair] = pair_totals.get(pair, 0.0) + row[7]
            pair_bilateral[frame] = max(pair_totals.values(), default=0.0)

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("time_s", "calcium", "strong_heads", "status2_links",
                         "max_single_link_bilateral_arm_load_pN",
                         "max_pair_summed_bilateral_arm_load_pN",
                         "max_productive_support_pN",
                         "max_link_extension_um", "max_raw_spring_bilateral_pN",
                         "max_allocated_motor_bilateral_pN"))
        writer.writerows(zip(time, calcium, strong, aa2, bilateral,
                             pair_bilateral, support, max_extension,
                             raw_bilateral, allocated_bilateral))

    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True,
                             constrained_layout=True)
    axes[0].plot(1e3 * time, calcium, color="#777777", label="calcium")
    head_axis = axes[0].twinx()
    head_axis.plot(1e3 * time, strong, color="#5e3c99", label="strong AM heads")
    axes[0].set_ylabel("calcium")
    head_axis.set_ylabel("strong heads")
    axes[1].plot(1e3 * time, bilateral, color="#18a558",
                 label="maximum effective single-link catch load")
    axes[1].plot(1e3 * time, pair_bilateral, color="#007f5f", alpha=.8,
                 label="maximum pair-summed arm load")
    axes[1].plot(1e3 * time, support, color="#5e3c99", alpha=.85,
                 label="max productive AM support")
    if np.any(np.abs(raw_bilateral - bilateral) > 1e-12):
        axes[1].plot(1e3 * time, raw_bilateral, color="#ef5a29", alpha=.55,
                     label="raw spring bilateral load (not used by KMC)")
    axes[1].set_ylabel("load (pN)")
    axes[1].legend(loc="upper left")
    axes[2].plot(1e3 * time, 1e3 * max_extension, color="#ef5a29",
                 label="maximum ACTN2 extension")
    status_axis = axes[2].twinx()
    status_axis.step(1e3 * time, aa2, where="post", color="#18a558",
                     alpha=.65, label="status-2 links")
    axes[2].set_ylabel("extension (nm)")
    status_axis.set_ylabel("status-2 links")
    axes[2].set_xlabel("time (ms)")
    for axis in axes:
        axis.grid(alpha=.2)
    args.plot.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.plot, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
