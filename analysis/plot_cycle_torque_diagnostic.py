#!/usr/bin/env python3
"""Resolve apparent heartbeat-boundary torque events in a minimal trajectory."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/sarcomere_matplotlib")

import h5py
import matplotlib.pyplot as plt
import numpy as np


def director_angle(reference: np.ndarray, values: np.ndarray) -> np.ndarray:
    dot = np.abs(np.sum(values * reference, axis=-1))
    return np.degrees(np.arccos(np.clip(dot, 0.0, 1.0)))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trajectory", type=Path)
    parser.add_argument("--start-ms", type=float, default=270.0)
    parser.add_argument("--stop-ms", type=float, default=310.0)
    parser.add_argument("--pulse-ms", type=float, default=300.0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with h5py.File(args.trajectory, "r", locking=False) as handle:
        dt = float(handle["meta"].attrs["dt"])
        offsets = np.asarray(handle["catch_bond/local_link_frame_offsets"], dtype=int)
        time_ms = 1e3 * offsets[:, 0] * dt
        actin_direction = np.asarray(handle["actin/direction"], dtype=float)
        myosin_direction = np.asarray(handle["myosin/direction"], dtype=float)
        actin_torque = np.asarray(handle["actin/torque"], dtype=float)
        myosin_torque = np.asarray(handle["myosin/torque"], dtype=float)
        actin_tangent = actin_torque - actin_direction * np.sum(
            actin_torque * actin_direction, axis=-1)[..., None]
        myosin_tangent = myosin_torque - myosin_direction * np.sum(
            myosin_torque * myosin_direction, axis=-1)[..., None]
        state = np.asarray(handle["state/am_head_state_current"], dtype=int)
        bound = np.asarray(handle["state/am_head_bound_actin_current"], dtype=int)
        occupied = bound >= 0
        strong = np.count_nonzero(occupied & (state == 2), axis=(1, 2))
        weak = np.count_nonzero(occupied & (state == 1), axis=(1, 2))
        identity = np.asarray(handle["catch_bond/local_link_identity"], dtype=int)
        aa1 = np.zeros(len(time_ms), dtype=int)
        aa2 = np.zeros(len(time_ms), dtype=int)
        for frame, (_, start, count) in enumerate(offsets):
            status = identity[start:start + count, 4]
            aa1[frame] = np.count_nonzero(status == 1)
            aa2[frame] = np.count_nonzero(status == 2)
        ruptures = np.asarray(
            handle.get("actin_myo/forced_ruptures_v2", np.empty((0, 19))),
            dtype=float,
        )

    selected = (time_ms >= args.start_ms) & (time_ms <= args.stop_ms)
    indices = np.flatnonzero(selected)
    first = indices[0]
    t = time_ms[selected]
    fig, axes = plt.subplots(4, 1, figsize=(10, 9), sharex=True)
    axes[0].plot(t, strong[selected], color="#54278f", label="strong AM")
    axes[0].plot(t, weak[selected], "--", color="#8073ac", label="weak AM")
    axes[0].set_ylabel("heads")
    twin = axes[0].twinx()
    twin.plot(t, aa2[selected], color="#009e73", label="status 2")
    twin.plot(t, aa1[selected], "--", color="#e69f00", label="status 1")
    twin.set_ylabel("ACTN2 links")
    axes[0].legend(frameon=False, loc="upper left", ncol=2)
    twin.legend(frameon=False, loc="upper right", ncol=2)

    for i, color in ((2, "#0072b2"), (3, "#d55e00")):
        axes[1].plot(t, np.linalg.norm(actin_torque[selected, i], axis=1),
                     color=color, alpha=0.25, label=f"actin {i} stored norm")
        axes[1].plot(t, np.linalg.norm(actin_tangent[selected, i], axis=1),
                     color=color, lw=1.8, label=f"actin {i} tangent")
    axes[1].set_ylabel("generalized torque\n(pN µm)")
    axes[1].legend(frameon=False, ncol=2, fontsize=8)

    for i, color in ((0, "#cc79a7"), (1, "#56b4e9")):
        axes[2].plot(t, np.linalg.norm(myosin_tangent[selected, i], axis=1),
                     color=color, label=f"myosin {i} tangent")
    axes[2].set_ylabel("myosin tangent\ntorque (pN µm)")
    axes[2].legend(frameon=False, ncol=2)

    for i, color in ((2, "#0072b2"), (3, "#d55e00")):
        axes[3].plot(t, director_angle(
            actin_direction[first, i], actin_direction[selected, i]),
            color=color, label=f"actin {i}")
    axes[3].plot(t, director_angle(
        myosin_direction[first, 0], myosin_direction[selected, 0]),
        color="#cc79a7", label="myosin 0")
    axes[3].set(xlabel="time (ms)", ylabel=f"angle from {time_ms[first]:g} ms (deg)")
    axes[3].legend(frameon=False, ncol=3)

    for axis in axes:
        axis.axvline(args.pulse_ms, color="#222222", ls=":", lw=1.2)
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", color="#eeeeee", lw=0.7)
    for row in ruptures:
        event_ms = 1e3 * row[2] * dt
        if args.start_ms <= event_ms <= args.stop_ms:
            for axis in axes:
                axis.axvline(event_ms, color="#e41a1c", ls="--", lw=0.9)
    axes[0].set_title(
        "End of cycle 3: dotted = next calcium pulse; dashed red = forced head rupture")
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=220, bbox_inches="tight")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
