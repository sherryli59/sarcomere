#!/usr/bin/env python3
"""Focused replay of A-M bare-zone penetration after bilateral titin perturbation."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/sarcomere_matplotlib")

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter


MYOSIN = "#2878b5"
BARE = "#d9e2ec"
ACTIN = ("#ef6c31", "#d83a8b", "#ef6c31", "#d83a8b")
TITIN = "#0aa6ad"
PENETRATION = "#d62728"
OVERLAP = "#5aa0d6"


def minimum_image(value: float, box: float) -> float:
    return value - box * np.rint(value / box)


def load(path: Path):
    with h5py.File(path, "r", locking=False) as h:
        data = {
            "ac": h["actin/center"][:],
            "ad": h["actin/direction"][:],
            "mc": h["myosin/center"][:],
            "md": h["myosin/direction"][:],
            "weak": h["actin/am_weak_count"][:, :, 0],
            "strong": h["actin/am_strong_count"][:, :, 0],
            "bound": h["state/titin_bound_actin_current"][:],
            "titin_force": h["state/titin_force_current"][:],
            "material_e1": h["state/myosin_material_e1_current"][:],
            "material_e2": h["state/myosin_material_e2_current"][:],
            "offsets": h["catch_bond/local_link_mechanics_v2_frame_offsets"][:],
            "dt": float(h["meta"].attrs["dt"]),
            "box": float(h["meta"].attrs["box_x"]),
            "actin_length": float(h["meta"].attrs["actin_length"]),
            "myosin_length": float(h["meta"].attrs["myosin_length"]),
            "bare_fraction": float(h["meta"].attrs["am_bare_zone_fraction"]),
            "links_per_half": int(h["meta"].attrs["titin_links_per_half"]),
            "myosin_radius": float(h["meta"].attrs["myosin_radius"]),
            "titin_k": float(h["meta"].attrs["titin_k"]),
            "titin_rest": float(h["meta"].attrs["titin_rest_length"]),
        }
    if data["bound"].ndim == 3:
        data["bound"] = data["bound"][:, :, 0]
    if data["titin_force"].ndim == 3:
        data["titin_force"] = data["titin_force"][:, :, 0]
    data["time"] = data["offsets"][:, 0] * data["dt"] * 1e3
    data["time"] -= data["time"][1]
    data["time"][0] = -0.5
    return data


def pair_geometry(data, frame: int, actin: int, myosin: int, half: int):
    axis = data["md"][frame, myosin]
    axis = axis / np.linalg.norm(axis)
    displacement = data["ac"][frame, actin] - data["mc"][frame, myosin]
    displacement[0] = minimum_image(displacement[0], data["box"])
    center = float(np.dot(displacement, axis))
    direction = float(np.dot(data["ad"][frame, actin], axis))
    barbed = center - 0.5 * data["actin_length"] * direction
    pointed = center + 0.5 * data["actin_length"] * direction
    lo, hi = sorted((barbed, pointed))
    myosin_half = 0.5 * data["myosin_length"]
    overlap = max(0.0, min(myosin_half, hi) - max(-myosin_half, lo))
    bare_half = 0.5 * data["bare_fraction"] * data["myosin_length"]
    penetration = (
        max(0.0, pointed + bare_half)
        if half == 0 else max(0.0, bare_half - pointed)
    )
    force = 0.0
    for titin_index, candidate in enumerate(data["bound"][frame]):
        if int(candidate) != actin:
            continue
        per_half = data["links_per_half"]
        anchor_myosin, local = divmod(titin_index, 2 * per_half)
        anchor_half, track = divmod(local, per_half)
        anchor_axis = data["md"][frame, anchor_myosin]
        anchor_axis = anchor_axis / np.linalg.norm(anchor_axis)
        phase = np.pi / per_half if per_half > 1 and anchor_half == 0 else 0.0
        azimuth = phase + 2.0 * np.pi * track / per_half
        radial = (
            np.cos(azimuth) * data["material_e1"][frame, anchor_myosin]
            + np.sin(azimuth) * data["material_e2"][frame, anchor_myosin]
        )
        half_sign = -1.0 if anchor_half == 0 else 1.0
        thick_anchor = (
            data["mc"][frame, anchor_myosin]
            + half_sign * 0.5 * data["myosin_length"] * anchor_axis
            + data["myosin_radius"] * radial
        )
        actin_barbed = (
            data["ac"][frame, actin]
            - 0.5 * data["actin_length"] * data["ad"][frame, actin]
        )
        connector = thick_anchor - actin_barbed
        connector[0] = minimum_image(connector[0], data["box"])
        distance = np.linalg.norm(connector)
        force = data["titin_k"] * max(0.0, distance - data["titin_rest"])
        break
    return {
        "center": center,
        "direction": direction,
        "barbed": barbed,
        "pointed": pointed,
        "lo": lo,
        "hi": hi,
        "overlap": overlap,
        "penetration": penetration,
        "force": force,
    }


def all_geometry(data, frame: int):
    owners = (0, 1, 0, 1)
    halves = (0, 1, 1, 0)
    return [
        pair_geometry(data, frame, actin, owners[actin], halves[actin])
        for actin in range(4)
    ]


def draw_myosin_panel(ax, data, frame: int, myosin: int, actins):
    myosin_half = 0.5 * data["myosin_length"]
    bare_half = 0.5 * data["bare_fraction"] * data["myosin_length"]
    ax.set_xlim(-2.08, 2.08)
    ax.set_ylim(-0.42, 0.42)
    ax.axvspan(-bare_half, bare_half, color=BARE, alpha=0.75, zorder=0)
    ax.plot([-myosin_half, myosin_half], [0, 0], color=MYOSIN, lw=13,
            solid_capstyle="round", zorder=1)
    ax.plot([-bare_half, bare_half], [0, 0], color="#f0f3f6", lw=13,
            solid_capstyle="butt", zorder=2)
    ax.text(0, 0.06, f"M{myosin} bare zone", ha="center", va="bottom",
            fontsize=9, color="#4a5865")
    geometry = all_geometry(data, frame)
    halves = (0, 1, 1, 0)
    for lane, actin in zip((0.23, -0.23), actins):
        g = geometry[actin]
        color = ACTIN[actin]
        ax.plot([g["lo"], g["hi"]], [lane, lane], color=color, lw=9,
                solid_capstyle="round", zorder=4)
        ax.scatter([g["barbed"]], [lane], s=72, facecolor="white",
                   edgecolor=color, linewidth=2.2, zorder=7)
        ax.scatter([g["pointed"]], [lane], s=80, marker="|", color=color,
                   linewidth=3, zorder=7)
        direction = np.sign(g["direction"])
        ax.annotate("", xy=(g["center"] + 0.22 * direction, lane),
                    xytext=(g["center"] - 0.05 * direction, lane),
                    arrowprops=dict(arrowstyle="-|>", color="white", lw=1.5),
                    zorder=8)
        overlap_lo = max(-myosin_half, g["lo"])
        overlap_hi = min(myosin_half, g["hi"])
        if overlap_hi > overlap_lo:
            ax.plot([overlap_lo, overlap_hi], [lane - 0.055, lane - 0.055],
                    color=OVERLAP, lw=5, alpha=0.9, solid_capstyle="butt")
        half = halves[actin]
        titin_anchor = -myosin_half if half == 0 else myosin_half
        taut = g["force"] > 1e-9
        ax.plot([titin_anchor, g["barbed"]], [0, lane], color=TITIN,
                lw=3.2 if taut else 2.0, ls="-" if taut else "--", zorder=5)
        if g["penetration"] > 0:
            boundary = -bare_half if half == 0 else bare_half
            ax.plot([boundary, g["pointed"]], [lane + 0.075, lane + 0.075],
                    color=PENETRATION, lw=6, solid_capstyle="butt", zorder=6)
            ax.text(0.5 * (boundary + g["pointed"]), lane + 0.12,
                    f"penetration {1e3*g['penetration']:.0f} nm",
                    ha="center", va="bottom", fontsize=9,
                    color=PENETRATION, weight="bold")
        ax.text(-2.02, lane, f"A{actin}", color=color, fontsize=11,
                weight="bold", va="center")
        ax.text(2.02, lane,
                f"AM {1e3*g['overlap']:.0f} nm  |  titin {g['force']:.1f} pN  |  "
                f"W{int(data['weak'][frame, actin])}/S{int(data['strong'][frame, actin])}",
                ha="right", va="center", fontsize=8.5, color="#3f4850")
    ax.set_title(f"Myosin {myosin}: coordinates relative to the myosin center",
                 fontsize=11, weight="bold")
    ax.set_yticks([])
    ax.set_xlabel("axial coordinate relative to myosin center (µm)")
    ax.grid(axis="x", alpha=0.16)
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)


def render(fig, data, frame: int, series):
    fig.clear()
    gs = fig.add_gridspec(2, 6, height_ratios=(1.15, 1.0), hspace=0.42,
                          wspace=0.70)
    ax0 = fig.add_subplot(gs[0, :3])
    ax1 = fig.add_subplot(gs[0, 3:])
    penetration_ax = fig.add_subplot(gs[1, :3])
    overlap_ax = fig.add_subplot(gs[1, 3:])
    draw_myosin_panel(ax0, data, frame, 0, (0, 2))
    draw_myosin_panel(ax1, data, frame, 1, (3, 1))

    time = data["time"]
    penetration_ax.plot(time, 1e3 * series["max_penetration"],
                        color=PENETRATION, lw=2.7,
                        label="maximum pointed-end penetration")
    penetration_ax.plot(time, 1e3 * series["midpoint_error"],
                        color=TITIN, lw=2.3,
                        label="maximum local midpoint error")
    penetration_ax.axhline(0, color="#555", lw=1)
    penetration_ax.set_ylabel("distance (nm)")
    penetration_ax.set_ylim(-8, 585)
    penetration_ax.legend(loc="upper right", fontsize=8)

    for actin in range(4):
        overlap_ax.plot(time, 1e3 * series["overlap"][:, actin],
                        color=ACTIN[actin], lw=2.1,
                        ls="-" if actin < 2 else "--", label=f"A{actin}")
    overlap_ax.axhline(340, color="#555", ls=":", lw=1.2,
                       label="centered start: 340 nm")
    overlap_ax.set_ylabel("A–M overlap (nm)")
    overlap_ax.set_ylim(-15, 940)
    overlap_ax.legend(loc="upper right", fontsize=8, ncol=3)

    for ax in (penetration_ax, overlap_ax):
        ax.axvline(time[frame], color="#111", lw=1.4)
        ax.set_xlim(time[0], time[-1])
        ax.set_xlabel("time after 550 nm displacement (ms)")
        ax.grid(alpha=0.18)

    geometry = all_geometry(data, frame)
    max_penetration = max(g["penetration"] for g in geometry)
    mean_overlap = np.mean([g["overlap"] for g in geometry])
    if frame == 0:
        phase = "CENTERED — bilateral titin already bound"
    elif frame == 1:
        phase = "DISPLACED — A0/A1 penetrate the bare zone"
    else:
        phase = "RELEASED — far-side titin tension removes bare-zone penetration"
    fig.suptitle(
        f"A–M burial diagnostic — {phase}\n"
        f"t = {time[frame]:.2f} ms  |  max penetration {1e3*max_penetration:.0f} nm"
        f"  |  mean A–M overlap {1e3*mean_overlap:.0f} nm",
        fontsize=14, weight="bold", y=0.98,
    )
    fig.text(
        0.02, 0.015,
        "open circle: barbed end  |  short bar: pointed end  |  cyan: titin "
        "(solid tension / dashed slack)  |  blue underline: A–M overlap  |  red: bare-zone penetration",
        fontsize=9, color="#444",
    )
    fig.subplots_adjust(left=0.055, right=0.975, bottom=0.105, top=0.86)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trajectory", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--fps", type=int, default=15)
    args = parser.parse_args()
    data = load(args.trajectory)
    all_frames = range(len(data["ac"]))
    geometry = [all_geometry(data, frame) for frame in all_frames]
    max_penetration = np.array([
        max(g["penetration"] for g in frame_geometry)
        for frame_geometry in geometry
    ])
    overlap = np.array([
        [g["overlap"] for g in frame_geometry]
        for frame_geometry in geometry
    ])
    midpoint_error = []
    for frame in all_frames:
        errors = []
        for myosin in range(2):
            items = [g for a, g in enumerate(geometry[frame])
                     if (0, 1, 0, 1)[a] == myosin]
            errors.append(abs(0.5 * sum(
                (item["barbed"] for item in items), start=0.0)))
        midpoint_error.append(max(errors))
    series = {
        "max_penetration": max_penetration,
        "overlap": overlap,
        "midpoint_error": np.asarray(midpoint_error),
    }

    # Hold the centered and displaced states long enough to inspect them, then
    # show the otherwise sub-frame 0--2 ms recoil in slow motion.
    movie_frames = [0] * 18 + [1] * 30
    movie_frames += [frame for frame in range(2, 7) for _ in range(6)]
    movie_frames += [frame for frame in range(7, len(data["ac"]))]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(16, 9), facecolor="white")
    writer = FFMpegWriter(
        fps=args.fps, codec="libx264",
        extra_args=("-pix_fmt", "yuv420p", "-crf", "18", "-threads", "1"),
    )
    with writer.saving(fig, str(args.output), dpi=120):
        for number, frame in enumerate(movie_frames, 1):
            render(fig, data, frame, series)
            writer.grab_frame(facecolor="white")
            print(f"rendered {number}/{len(movie_frames)}", flush=True)
    plt.close(fig)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
