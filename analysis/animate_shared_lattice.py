#!/usr/bin/env python3
"""Render the corrected shared 4-thick/8-thin lattice with live KMC counts."""
import argparse
import os
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/sarcomere_matplotlib")

import h5py
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib.patches import Circle, Rectangle
import numpy as np


def pbc_delta(delta, box):
    result = np.asarray(delta, dtype=float).copy()
    result -= box * np.round(result / box)
    return result


def rod(axis, center, direction, length, color, width):
    ends = np.vstack((
        center - 0.5 * length * direction,
        center + 0.5 * length * direction))
    axis.plot(*ends.T, color=color, lw=width, solid_capstyle="round")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("trajectory", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--stride", type=int, default=4)
    parser.add_argument("--fps", type=int, default=18)
    args = parser.parse_args()
    with h5py.File(args.trajectory, "r", locking=False) as h:
        ac = np.asarray(h["actin/center"])
        au = np.asarray(h["actin/direction"])
        mc = np.asarray(h["myosin/center"])
        mu = np.asarray(h["myosin/direction"])
        steps = np.asarray(h["state/current_step"])[:, 0]
        states = np.asarray(h["state/am_head_state_current"])
        bound = np.asarray(h["state/am_head_bound_actin_current"])
        pair_count = np.asarray(h["state/aa_pairs_current_count"])[:, 0]
        pairs = np.asarray(h["state/aa_pairs_current"])
        pair_status = np.asarray(h["state/aa_status_current"])[:, :, 0]
        meta = h["meta"].attrs
        dt = float(meta["dt"])
        box = np.array([meta["box_x"], meta["box_y"], meta["box_z"]],
                       dtype=float)
        al = float(meta["actin_length"])
        ml = float(meta["myosin_length"])
        ar = float(meta["actin_radius"])
        mr = float(meta.get(
            "myosin_radius",
            meta["am_preferred_centerline"] -
            meta["actin_radius"] - meta["am_surface_gap"]))
        if "actin/aa_local_status1_count" in h:
            s1 = 0.5 * np.asarray(
                h["actin/aa_local_status1_count"]).reshape(len(steps), 16).sum(1)
            s2 = 0.5 * np.asarray(
                h["actin/aa_local_status2_count"]).reshape(len(steps), 16).sum(1)
        else:
            s1 = np.array([
                np.count_nonzero(pair_status[f, :int(pair_count[f])] == 1)
                for f in range(len(steps))])
            s2 = np.array([
                np.count_nonzero(pair_status[f, :int(pair_count[f])] == 2)
                for f in range(len(steps))])
    time = steps * dt
    weak = np.count_nonzero((states == 1) & (bound >= 0), axis=(1, 2))
    strong = np.count_nonzero((states == 2) & (bound >= 0), axis=(1, 2))

    colors = ("#0072B2", "#D55E00")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(14, 9), facecolor="white")
    writer = FFMpegWriter(
        fps=args.fps, codec="libx264",
        extra_args=("-pix_fmt", "yuv420p", "-crf", "18", "-threads", "1"))
    frames = range(0, len(time), max(1, args.stride))
    with writer.saving(fig, str(args.output), dpi=120):
        for f in frames:
            fig.clear()
            grid = fig.add_gridspec(
                2, 3, height_ratios=(1.0, 0.46),
                width_ratios=(1.35, 1.0, 1.0), hspace=0.3, wspace=0.28)
            ax3 = fig.add_subplot(grid[0, 0], projection="3d")
            ax3.set_proj_type("ortho")
            for m in range(4):
                rod(ax3, mc[f, m], mu[f, m], ml, "#263238", 7)
            for a in range(16):
                rod(ax3, ac[f, a], au[f, a], al,
                    colors[0 if a < 8 else 1], 3)
            ax3.set(
                xlim=(-1.1, 1.1),
                ylim=(-0.55 * box[1], 0.55 * box[1]),
                zlim=(-0.55 * box[2], 0.55 * box[2]),
                xlabel="x (µm)", ylabel="y (µm)", zlabel="z (µm)")
            ax3.set_box_aspect((3.0, 1.0, 1.0))
            ax3.view_init(elev=18, azim=-66)
            ax3.set_title("force-balanced sarcomere cell")

            for half in range(2):
                axis = fig.add_subplot(grid[0, 1 + half])
                axis.add_patch(Rectangle(
                    (-0.5e3 * box[1], -0.5e3 * box[2]),
                    1e3 * box[1], 1e3 * box[2],
                    fc="none", ec="#9E9E9E", lw=1.2))
                for m in range(4):
                    y = 1e3 * ((mc[f, m, 1] + 0.5 * box[1]) %
                              box[1] - 0.5 * box[1])
                    z = 1e3 * ((mc[f, m, 2] + 0.5 * box[2]) %
                              box[2] - 0.5 * box[2])
                    axis.add_patch(Circle(
                        (y, z), 1e3 * mr, fc="#37474F",
                        ec="white", lw=0.7, zorder=3))
                offset = 8 * half
                wrapped = {}
                for a in range(offset, offset + 8):
                    y = ((ac[f, a, 1] + 0.5 * box[1]) %
                         box[1] - 0.5 * box[1])
                    z = ((ac[f, a, 2] + 0.5 * box[2]) %
                         box[2] - 0.5 * box[2])
                    wrapped[a] = np.array([y, z])
                for slot in range(int(pair_count[f])):
                    i, j = map(int, pairs[f, slot])
                    if i not in wrapped or j not in wrapped:
                        continue
                    start = wrapped[i]
                    delta = wrapped[j] - start
                    delta -= box[1:] * np.round(delta / box[1:])
                    status_color = (
                        "#009E73" if pair_status[f, slot] == 2 else "#E69F00")
                    axis.plot(
                        1e3 * np.array([start[0], start[0] + delta[0]]),
                        1e3 * np.array([start[1], start[1] + delta[1]]),
                        color=status_color, lw=1.2, alpha=0.75, zorder=1)
                for a, point in wrapped.items():
                    axis.add_patch(Circle(
                        tuple(1e3 * point), 1e3 * ar,
                        fc=colors[half], ec="white", lw=0.7, zorder=4))
                axis.set(
                    xlim=(-0.55e3 * box[1], 0.55e3 * box[1]),
                    ylim=(-0.55e3 * box[2], 0.55e3 * box[2]),
                    aspect="equal", xlabel="y (nm)", ylabel="z (nm)",
                    title=("-x half: + polarity" if half == 0
                           else "+x half: − polarity"))
                axis.grid(alpha=0.15)

            counts = fig.add_subplot(grid[1, :])
            counts.plot(1e3 * time, weak, color="#56B4E9", label="weak AM")
            counts.plot(1e3 * time, strong, color="#D55E00", label="strong AM")
            counts.plot(1e3 * time, s1, color="#E69F00",
                        label="status-1 ACTN2")
            counts.plot(1e3 * time, s2, color="#009E73",
                        label="status-2 ACTN2")
            counts.axvline(1e3 * time[f], color="black", lw=1)
            counts.set(
                xlim=(0, 1e3 * time[-1]), xlabel="time (ms)",
                ylabel="molecular count")
            counts.legend(ncol=4, loc="upper right")
            counts.grid(alpha=0.2)
            fig.suptitle(
                "Shared 4-thick/8-thin lattice "
                f"(each thin→3 thick, each thick→6 thin) — "
                f"t={1e3*time[f]:.2f} ms",
                fontsize=14, weight="bold")
            writer.grab_frame(facecolor="white")
    plt.close(fig)
    print(args.output)


if __name__ == "__main__":
    main()
