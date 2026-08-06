#!/usr/bin/env python3
"""Render only current cb2 actins and myosins carrying heads on them."""

from __future__ import annotations

import argparse
import os
import subprocess
import tempfile
from itertools import product
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/sarcomere-mplconfig")

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


WEAK, STRONG = 1, 2


def minimum_image(delta: np.ndarray, box: np.ndarray) -> np.ndarray:
    return delta - box * np.rint(delta / box)


def clipped_periodic_segments(
    starts: np.ndarray, ends: np.ndarray, box: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    lower, upper = -0.5 * box, 0.5 * box
    clipped_starts, clipped_ends = [], []
    for start, end in zip(np.asarray(starts), np.asarray(ends)):
        choices = []
        for axis in range(3):
            offsets = [0.0]
            if min(start[axis], end[axis]) < lower[axis]:
                offsets.append(box[axis])
            if max(start[axis], end[axis]) > upper[axis]:
                offsets.append(-box[axis])
            choices.append(offsets)
        for offset in product(*choices):
            first = start + np.asarray(offset)
            delta = end - start
            enter, leave = 0.0, 1.0
            valid = True
            for axis in range(3):
                if abs(delta[axis]) < 1e-14:
                    if first[axis] < lower[axis] or first[axis] > upper[axis]:
                        valid = False
                    continue
                t0 = (lower[axis] - first[axis]) / delta[axis]
                t1 = (upper[axis] - first[axis]) / delta[axis]
                enter = max(enter, min(t0, t1))
                leave = min(leave, max(t0, t1))
            if valid and leave > enter + 1e-12:
                clipped_starts.append(first + enter * delta)
                clipped_ends.append(first + leave * delta)
    return (np.asarray(clipped_starts).reshape((-1, 3)),
            np.asarray(clipped_ends).reshape((-1, 3)))


def draw_segments(ax, starts, ends, color, width, alpha=1.0, zorder=2):
    starts, ends = np.asarray(starts), np.asarray(ends)
    for start, end in zip(starts, ends):
        ax.plot(
            [start[0], end[0]], [start[1], end[1]], [start[2], end[2]],
            color=color, linewidth=width, alpha=alpha, solid_capstyle="round",
            zorder=zorder,
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("trajectory", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--fps", type=int, default=8)
    parser.add_argument("--start-frame", type=int)
    parser.add_argument("--stop-frame", type=int)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(args.trajectory, "r", locking=False) as handle:
        meta = handle["meta"].attrs
        actin_center = handle["actin/center"][:]
        actin_direction = handle["actin/direction"][:]
        myosin_center = handle["myosin/center"][:]
        myosin_direction = handle["myosin/direction"][:]
        cb_status = handle["actin/cb_status"][:, :, 0]
        steps = handle["state/current_step"][:, 0]
        head_state = handle["state/am_head_state_current"][:]
        head_actin = handle["state/am_head_bound_actin_current"][:]
        head_target = handle["state/am_head_target_axial_current"][:]
        head_azimuth = handle["state/am_head_binding_azimuth_current"][:]
        material_e1 = handle["state/myosin_material_e1_current"][:]
        material_e2 = handle["state/myosin_material_e2_current"][:]
        link_offsets = handle["catch_bond/local_link_frame_offsets"][:]
        link_identity = handle["catch_bond/local_link_identity"][:]
        link_coords = handle["catch_bond/local_link_coordinates"][:]
        dt = float(meta["dt"])
        actin_length = float(meta["actin_length"])
        myosin_length = float(meta["myosin_length"])
        myosin_radius = float(meta["myosin_radius"])
        bare_fraction = float(meta["am_bare_zone_fraction"])
        box = np.array([meta["box_x"], meta["box_y"], meta["box_z"]], dtype=float)

    cb2_count = np.sum(cb_status == 2, axis=1)
    first_cb2 = int(np.flatnonzero(cb2_count > 0)[0]) if np.any(cb2_count > 0) else 0
    first = args.start_frame if args.start_frame is not None else first_cb2
    last = args.stop_frame if args.stop_frame is not None else len(steps)
    frames = list(range(max(0, first), min(last, len(steps))))
    if not frames:
        raise SystemExit("no selected frames")

    colors = {
        "actin_plus": "#79b9ee",
        "actin_minus": "#164f86",
        "myosin": "#ef7f35",
        "aa2": "#d43f8d",
        "weak": "#e3b124",
        "strong": "#d84b35",
        "barbed": "#2a9d68",
    }

    with tempfile.TemporaryDirectory(prefix="cb2_video_", dir=args.output.parent) as temporary:
        directory = Path(temporary)
        for output_index, frame in enumerate(frames):
            centers_a = actin_center[frame]
            directions_a = actin_direction[frame]
            centers_m = myosin_center[frame]
            directions_m = myosin_direction[frame]
            visible_actin = cb_status[frame] == 2

            states = head_state[frame]
            bound = head_actin[frame]
            attached_rows = np.argwhere(
                np.isin(states, (WEAK, STRONG))
                & (bound >= 0)
                & visible_actin[np.maximum(bound, 0)]
            )
            visible_myosin_ids = np.unique(attached_rows[:, 0]) if len(attached_rows) else np.empty(0, int)
            visible_myosin = np.zeros(len(centers_m), dtype=bool)
            visible_myosin[visible_myosin_ids] = True

            fig = plt.figure(figsize=(14, 7.6), facecolor="white")
            ax = fig.add_subplot(111, projection="3d")
            ax.set_facecolor("white")

            # Filament rods.
            a_half = 0.5 * actin_length * directions_a
            a_starts, a_ends = centers_a - a_half, centers_a + a_half
            for polarity, color in ((True, colors["actin_plus"]), (False, colors["actin_minus"])):
                mask = visible_actin & ((directions_a[:, 0] >= 0) == polarity)
                starts, ends = clipped_periodic_segments(a_starts[mask], a_ends[mask], box)
                draw_segments(ax, starts, ends, color, 4.4, 0.98, 4)
            m_half = 0.5 * myosin_length * directions_m
            m_starts, m_ends = clipped_periodic_segments(
                (centers_m - m_half)[visible_myosin],
                (centers_m + m_half)[visible_myosin], box,
            )
            draw_segments(ax, m_starts, m_ends, colors["myosin"], 6.0, 0.88, 2)

            # Barbed ends are s=-L/2 in the simulator.
            barbed = centers_a[visible_actin] - 0.5 * actin_length * directions_a[visible_actin]
            barbed = (barbed + 0.5 * box) % box - 0.5 * box
            if len(barbed):
                ax.scatter(barbed[:, 0], barbed[:, 1], barbed[:, 2], s=24,
                           color=colors["barbed"], edgecolors="white", linewidths=.5,
                           depthshade=False, zorder=6)

            # Mature ACTN2 molecules connecting two currently visible cb2 actins.
            _, offset, count = map(int, link_offsets[frame])
            identity = link_identity[offset : offset + count]
            coords = link_coords[offset : offset + count]
            if len(identity):
                keep = (
                    (identity[:, 4] == 2)
                    & visible_actin[identity[:, 2]]
                    & visible_actin[identity[:, 3]]
                )
                identity, coords = identity[keep], coords[keep]
                if len(identity):
                    ai, aj = identity[:, 2], identity[:, 3]
                    aa_start = centers_a[ai] + coords[:, :1] * directions_a[ai]
                    raw_end = centers_a[aj] + coords[:, 1:] * directions_a[aj]
                    aa_end = aa_start + minimum_image(raw_end - aa_start, box)
                    aa_start, aa_end = clipped_periodic_segments(aa_start, aa_end, box)
                    draw_segments(ax, aa_start, aa_end, colors["aa2"], 2.2, 0.92, 5)

            # Material AM heads on the selected cb2 actins.
            weak_segments, strong_segments = [], []
            heads_per_half = states.shape[1] // 2
            bare_half = 0.5 * bare_fraction * myosin_length
            head_region = 0.5 * myosin_length - bare_half
            for myosin_id, head_id in attached_rows:
                actin_id = int(bound[myosin_id, head_id])
                half_index, local_head = divmod(int(head_id), heads_per_half)
                half_sign = -1 if half_index == 0 else 1
                fraction = (local_head + 0.5) / heads_per_half
                anchor_axial = half_sign * (bare_half + fraction * head_region)
                angle = head_azimuth[frame, myosin_id, head_id]
                anchor = (
                    centers_m[myosin_id]
                    + anchor_axial * directions_m[myosin_id]
                    + myosin_radius * (
                        np.cos(angle) * material_e1[frame, myosin_id]
                        + np.sin(angle) * material_e2[frame, myosin_id]
                    )
                )
                actin_center_image = anchor + minimum_image(centers_a[actin_id] - anchor, box)
                target = float(head_target[frame, myosin_id, head_id])
                site = actin_center_image + target * directions_a[actin_id]
                segment = (anchor, site)
                (weak_segments if states[myosin_id, head_id] == WEAK else strong_segments).append(segment)
            for segments, color, width in (
                (weak_segments, colors["weak"], 0.9),
                (strong_segments, colors["strong"], 1.5),
            ):
                if segments:
                    starts, ends = zip(*segments)
                    starts, ends = clipped_periodic_segments(np.asarray(starts), np.asarray(ends), box)
                    draw_segments(ax, starts, ends, color, width, 0.82, 5)

            ax.set_xlim(-0.5 * box[0], 0.5 * box[0])
            ax.set_ylim(-0.5 * box[1], 0.5 * box[1])
            ax.set_zlim(-0.5 * box[2], 0.5 * box[2])
            ax.set_box_aspect(box)
            ax.view_init(elev=20, azim=-67)
            ax.set_xlabel("x (µm)", labelpad=8)
            ax.set_ylabel("y (µm)", labelpad=8)
            ax.set_zlabel("z (µm)", labelpad=8)
            ax.grid(True, alpha=.18)
            ax.set_title(
                f"t = {steps[frame] * dt * 1e3:.0f} ms   |   "
                f"cb2 actins {np.count_nonzero(visible_actin)}   |   "
                f"attached myosins {np.count_nonzero(visible_myosin)}   |   "
                f"weak/strong heads {len(weak_segments)}/{len(strong_segments)}",
                pad=14,
            )
            legend = [
                Line2D([0], [0], color=colors["actin_plus"], lw=4, label="actin polarity +x"),
                Line2D([0], [0], color=colors["actin_minus"], lw=4, label="actin polarity −x"),
                Line2D([0], [0], color=colors["myosin"], lw=5, label="attached myosin"),
                Line2D([0], [0], color=colors["aa2"], lw=2, label="status-2 ACTN2"),
                Line2D([0], [0], color=colors["weak"], lw=1.5, label="weak head"),
                Line2D([0], [0], color=colors["strong"], lw=2, label="strong head"),
                Line2D([0], [0], marker="o", color="none", markerfacecolor=colors["barbed"],
                       markeredgecolor="white", markersize=7, label="barbed end"),
            ]
            ax.legend(handles=legend, loc="upper center", bbox_to_anchor=(0.5, -0.02),
                      ncol=4, frameon=False)
            fig.subplots_adjust(left=.02, right=.98, top=.91, bottom=.13)
            destination = directory / f"frame_{output_index:05d}.png"
            fig.savefig(destination, dpi=135, facecolor="white")
            if output_index == len(frames) - 1:
                fig.savefig(args.output.with_suffix(".png"), dpi=180, facecolor="white")
            plt.close(fig)
            print(f"rendered {output_index + 1}/{len(frames)}", flush=True)

        subprocess.run([
            "ffmpeg", "-y", "-framerate", str(args.fps),
            "-i", str(directory / "frame_%05d.png"),
            "-c:v", "libx264", "-preset", "slow", "-crf", "18",
            "-pix_fmt", "yuv420p", "-movflags", "+faststart", str(args.output),
        ], check=True)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
