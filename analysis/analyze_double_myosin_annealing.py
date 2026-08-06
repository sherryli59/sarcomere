#!/usr/bin/env python3
"""Measure axial register annealing in the four-actin/four-myosin rig."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


INTENDED_PAIRS = ((0, 0), (0, 2), (1, 1), (1, 3),
                  (2, 0), (2, 2), (3, 1), (3, 3))


def minimum_image(displacement, box, periodic):
    result = np.asarray(displacement, dtype=float).copy()
    result[..., periodic] -= (
        box[periodic] * np.rint(result[..., periodic] / box[periodic])
    )
    return result


def axial_register(first_center, first_direction, second_center,
                   second_direction, box, periodic):
    displacement = minimum_image(second_center - first_center, box, periodic)
    mean_axis = first_direction + second_direction
    mean_axis /= np.linalg.norm(mean_axis, axis=1)[:, None]
    return np.sum(displacement * mean_axis, axis=1)


def projected_overlap(actin_center, actin_direction, myosin_center,
                      myosin_direction, actin_length, myosin_length,
                      box, periodic):
    displacement = minimum_image(actin_center - myosin_center, box, periodic)
    actin_midpoint = np.sum(displacement * myosin_direction, axis=1)
    actin_half_extent = 0.5 * actin_length * np.abs(
        np.sum(actin_direction * myosin_direction, axis=1)
    )
    lower = np.maximum(-0.5 * myosin_length,
                       actin_midpoint - actin_half_extent)
    upper = np.minimum(0.5 * myosin_length,
                       actin_midpoint + actin_half_extent)
    return np.maximum(0.0, upper - lower)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("trajectory", type=Path)
    parser.add_argument("output_prefix", type=Path)
    args = parser.parse_args()

    with h5py.File(args.trajectory, "r", locking=False) as handle:
        steps = handle["state/current_step"][:, 0]
        dt = float(handle["meta"].attrs["dt"])
        actin_center = handle["actin/center"][:]
        actin_direction = handle["actin/direction"][:]
        myosin_center = handle["myosin/center"][:]
        myosin_direction = handle["myosin/direction"][:]
        head_actin = handle["state/am_head_bound_actin_current"][:]
        head_state = handle["state/am_head_state_current"][:]
        cb_status = handle["actin/cb_status"][:, :, 0]
        box = np.array([
            float(handle["meta"].attrs["box_x"]),
            float(handle["meta"].attrs["box_y"]),
            float(handle["meta"].attrs["box_z"]),
        ])
        periodic = np.array([
            bool(handle["meta"].attrs["periodic_x"]),
            bool(handle["meta"].attrs["periodic_y"]),
            bool(handle["meta"].attrs["periodic_z"]),
        ])
        actin_length = float(handle["meta"].attrs["actin_length"])
        myosin_length = float(handle["meta"].attrs["myosin_length"])

    time_ms = steps * dt * 1e3
    right_register = axial_register(
        myosin_center[:, 0], myosin_direction[:, 0],
        myosin_center[:, 2], myosin_direction[:, 2], box, periodic)
    left_register = axial_register(
        myosin_center[:, 1], myosin_direction[:, 1],
        myosin_center[:, 3], myosin_direction[:, 3], box, periodic)
    mean_error = 0.5 * (np.abs(right_register) + np.abs(left_register))

    overlaps = np.column_stack([
        projected_overlap(
            actin_center[:, actin], actin_direction[:, actin],
            myosin_center[:, myosin], myosin_direction[:, myosin],
            actin_length, myosin_length, box, periodic)
        for actin, myosin in INTENDED_PAIRS
    ])

    support = np.zeros((len(steps), 4), dtype=int)
    for frame in range(len(steps)):
        for actin in range(4):
            support[frame, actin] = sum(
                np.any((head_actin[frame, myosin] == actin) &
                       np.isin(head_state[frame, myosin], (1, 2)))
                for myosin in range(4)
            )

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    table = np.column_stack((
        time_ms, 1e3 * right_register, 1e3 * left_register,
        1e3 * mean_error, overlaps,
        support, np.count_nonzero(cb_status == 2, axis=1),
    ))
    header = ["time_ms", "right_register_nm", "left_register_nm",
              "mean_abs_register_nm"]
    header += [f"overlap_A{a}_M{m}_um" for a, m in INTENDED_PAIRS]
    header += [f"attached_myosins_A{a}" for a in range(4)]
    header += ["status2_actins"]
    np.savetxt(args.output_prefix.with_suffix(".csv"), table,
               delimiter=",", header=",".join(header), comments="")

    fig, axes = plt.subplots(3, 1, figsize=(8.2, 8.4), sharex=True)
    axes[0].plot(time_ms, 1e3 * np.abs(right_register), lw=2,
                 label="right myosin pair")
    axes[0].plot(time_ms, 1e3 * np.abs(left_register), lw=2,
                 label="left myosin pair")
    axes[0].plot(time_ms, 1e3 * mean_error, color="black", lw=2.6,
                 label="pair mean")
    axes[0].set_ylabel("myosin register error (nm)")
    axes[0].legend(frameon=False, ncol=3, fontsize=9)

    for column, (actin, myosin) in enumerate(INTENDED_PAIRS):
        axes[1].plot(time_ms, overlaps[:, column], lw=1.2, alpha=0.7,
                     label=f"A{actin}–M{myosin}")
    axes[1].plot(time_ms, overlaps.mean(axis=1), color="black", lw=2.5,
                 label="mean")
    axes[1].set_ylabel("axial A–M overlap (µm)")
    axes[1].legend(frameon=False, ncol=5, fontsize=7)

    for actin in range(4):
        axes[2].step(time_ms, support[:, actin], where="post", lw=1.5,
                     label=f"A{actin}")
    axes[2].set_ylim(-0.1, 2.25)
    axes[2].set_yticks((0, 1, 2))
    axes[2].set_ylabel("occupied myosin neighbors")
    axes[2].set_xlabel("time (ms)")
    axes[2].legend(frameon=False, ncol=4)
    fig.suptitle("Axial annealing of two myosins sharing each actin")
    fig.tight_layout()
    fig.savefig(args.output_prefix.with_suffix(".png"), dpi=200)
    plt.close(fig)

    final = mean_error[-1] * 1e3
    initial = mean_error[0] * 1e3
    print(f"mean register error: {initial:.1f} -> {final:.1f} nm")
    print(f"reduction: {(1.0 - final / initial) * 100:.1f}%")
    print(f"two-myosin occupancy: {np.mean(support >= 2) * 100:.1f}% of actin-frames")
    print(f"final status-2 actins: {np.count_nonzero(cb_status[-1] == 2)}/4")


if __name__ == "__main__":
    main()
