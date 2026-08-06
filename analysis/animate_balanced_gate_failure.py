#!/usr/bin/env python3
"""Focused replay of the two A--A pairs in the minimal balanced gate."""

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


PAIR_COLORS = ("#177b57", "#d34a32")
ACTIN_COLORS = ("#ef6c31", "#d83a8b", "#ef6c31", "#d83a8b")
MYOSIN_COLOR = "#2777b8"
WEAK_COLOR = "#f2b84b"
STRONG_COLOR = "#8c2d72"
STATUS1_COLOR = "#e3a008"
STATUS2_COLOR = "#138f62"
TITIN_COLOR = "#13a7ad"


def minimum_image(value: float, box: float) -> float:
    return value - box * np.rint(value / box)


def load(path: Path):
    with h5py.File(path, "r", locking=False) as h:
        d = {
            "ac": h["actin/center"][:],
            "ad": h["actin/direction"][:],
            "mc": h["myosin/center"][:],
            "md": h["myosin/direction"][:],
            "weak": h["actin/am_weak_count"][:, :, 0],
            "strong": h["actin/am_strong_count"][:, :, 0],
            "active_x": h["actin/am_active_force_x"][:, :, 0],
            "productive_overlap": h["actin/am_overlap_productive_weighted"][:, :, 0],
            "titin": h["actin/titin_bound"][:, :, 0],
            "identity": h["catch_bond/local_link_identity"][:],
            "coordinates": h["catch_bond/local_link_coordinates"][:],
            "mechanics": h["catch_bond/local_link_mechanics_v2"][:],
            "offsets": h["catch_bond/local_link_mechanics_v2_frame_offsets"][:],
            "dt": float(h["meta"].attrs["dt"]),
            "box": float(h["meta"].attrs["box_x"]),
            "actin_length": float(h["meta"].attrs["actin_length"]),
            "myosin_length": float(h["meta"].attrs["myosin_length"]),
        }
    d["time"] = d["offsets"][:, 0] * d["dt"] * 1e3
    frames = len(d["time"])
    d["overlap"] = np.zeros((frames, 2))
    d["registry"] = np.full((frames, 2), np.nan)
    d["support"] = np.full((frames, 2), np.nan)
    d["status1"] = np.zeros((frames, 2), dtype=int)
    d["status2"] = np.zeros((frames, 2), dtype=int)
    for frame, (_, identity_start, mechanics_start, count) in enumerate(d["offsets"]):
        identity_start, mechanics_start, count = map(
            int, (identity_start, mechanics_start, count)
        )
        ids = d["identity"][identity_start:identity_start + count]
        mechanics = d["mechanics"][mechanics_start:mechanics_start + count]
        for pair, (a, b) in enumerate(((0, 1), (2, 3))):
            dx = minimum_image(d["ac"][frame, b, 0] - d["ac"][frame, a, 0], d["box"])
            d["overlap"][frame, pair] = max(0.0, d["actin_length"] - abs(dx))
            mask = np.array(
                [{int(row[2]), int(row[3])} == {a, b} for row in ids], dtype=bool
            )
            if not np.any(mask):
                continue
            statuses = ids[mask, 4]
            selected = mechanics[mask]
            d["registry"][frame, pair] = np.median(selected[:, 17])
            d["support"][frame, pair] = np.max(selected[:, 8])
            d["status1"][frame, pair] = np.count_nonzero(statuses == 1)
            d["status2"][frame, pair] = np.count_nonzero(statuses == 2)
    return d


def pair_local_x(d, frame, pair):
    a, b = ((0, 1), (2, 3))[pair]
    length = d["actin_length"]
    barbed_a = d["ac"][frame, a, 0] - 0.5 * length * d["ad"][frame, a, 0]
    barbed_b_raw = d["ac"][frame, b, 0] - 0.5 * length * d["ad"][frame, b, 0]
    barbed_b = barbed_a + minimum_image(barbed_b_raw - barbed_a, d["box"])
    origin = 0.5 * (barbed_a + barbed_b)
    centers = {}
    for actin in (a, b):
        barbed = barbed_a if actin == a else barbed_b
        centers[actin] = barbed + 0.5 * length * d["ad"][frame, actin, 0] - origin
    myosins = {}
    owners = {0: 0, 1: 1, 2: 0, 3: 1}
    for actin in (a, b):
        owner = owners[actin]
        raw = d["mc"][frame, owner, 0]
        myosins[owner] = centers[actin] + minimum_image(
            raw - d["ac"][frame, actin, 0], d["box"]
        )
    return origin, centers, myosins


def active_links(d, frame, pair):
    _, identity_start, _, count = d["offsets"][frame]
    identity_start, count = int(identity_start), int(count)
    a, b = ((0, 1), (2, 3))[pair]
    result = []
    for row in range(identity_start, identity_start + count):
        ident = d["identity"][row]
        if {int(ident[2]), int(ident[3])} != {a, b}:
            continue
        result.append((int(ident[2]), int(ident[3]), int(ident[4]),
                       *d["coordinates"][row]))
    return result


def draw_pair(ax, d, frame, pair):
    a, b = ((0, 1), (2, 3))[pair]
    _, centers, myosins = pair_local_x(d, frame, pair)
    actin_y = {a: 0.16, b: -0.16}
    ax.set_xlim(-1.15, 1.15)
    ax.set_ylim(-0.43, 0.43)
    ax.axvline(0.0, color="#d7dde4", lw=1.2, zorder=0)
    for owner, center in myosins.items():
        ax.plot(
            [center - 0.5 * d["myosin_length"], center + 0.5 * d["myosin_length"]],
            [0.0, 0.0], color=MYOSIN_COLOR, lw=11, alpha=0.24,
            solid_capstyle="round", zorder=1,
        )
        ax.text(center, 0.02, f"M{owner}", color=MYOSIN_COLOR,
                fontsize=8, ha="center", va="bottom", clip_on=True)
    for actin in (a, b):
        center = centers[actin]
        direction = d["ad"][frame, actin, 0]
        lo, hi = sorted((center - 0.5 * d["actin_length"] * direction,
                         center + 0.5 * d["actin_length"] * direction))
        y = actin_y[actin]
        ax.plot([lo, hi], [y, y], color=ACTIN_COLORS[actin], lw=8,
                solid_capstyle="round", zorder=3)
        barbed = center - 0.5 * d["actin_length"] * direction
        ax.scatter([barbed], [y], s=80, marker="o", facecolor="white",
                   edgecolor=ACTIN_COLORS[actin], linewidth=2.2, zorder=6)
        ax.annotate("", xy=(center + 0.40 * direction, y),
                    xytext=(center + 0.18 * direction, y),
                    arrowprops=dict(arrowstyle="-|>", color="white", lw=1.5),
                    zorder=7)
        ax.text(lo if direction > 0 else hi, y + (0.07 if y > 0 else -0.08),
                f"A{actin}", color=ACTIN_COLORS[actin], weight="bold",
                fontsize=10, ha="center")
        weak = int(d["weak"][frame, actin])
        strong = int(d["strong"][frame, actin])
        ax.text(0.98, 0.91 if y > 0 else 0.09,
                f"A{actin}: weak {weak}   strong {strong}",
                transform=ax.transAxes, ha="right", va="center", fontsize=9,
                color=STRONG_COLOR if strong else "#5d6670",
                bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.82))
        force = d["active_x"][frame, actin]
        if abs(force) > 0.05:
            dx = np.sign(force) * min(0.32, 0.06 + 0.035 * abs(force))
            ax.annotate("", xy=(center + dx, y), xytext=(center, y),
                        arrowprops=dict(arrowstyle="->", lw=2.5,
                                        color=STRONG_COLOR), zorder=8)
    for i, j, status, s_i, s_j in active_links(d, frame, pair):
        x_i = centers[i] + s_i * d["ad"][frame, i, 0]
        x_j = centers[j] + s_j * d["ad"][frame, j, 0]
        color = STATUS2_COLOR if status == 2 else STATUS1_COLOR
        ax.plot([x_i, x_j], [actin_y[i], actin_y[j]], color=color,
                lw=3.0, alpha=0.95, zorder=5)
    overlap_nm = 1e3 * d["overlap"][frame, pair]
    registry_nm = 1e3 * d["registry"][frame, pair]
    support = d["support"][frame, pair]
    locked = d["status2"][frame, pair] > 0
    state = "LOCKED" if locked else "status 1 / mobile"
    ax.set_title(
        f"A{a}–A{b}: {state}\n"
        f"overlap {overlap_nm:.0f} nm  |  registry error {registry_nm:.0f} nm  |  support {support:.2f} pN",
        color=PAIR_COLORS[pair], fontsize=11, weight="bold",
    )
    ax.set_yticks([])
    ax.set_xlabel("local axial coordinate (µm)")
    ax.grid(axis="x", alpha=0.16)
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)


def render(fig, d, frame):
    fig.clear()
    gs = fig.add_gridspec(2, 6, height_ratios=(1.15, 1.0), hspace=0.48, wspace=0.82)
    pair0 = fig.add_subplot(gs[0, :3])
    pair1 = fig.add_subplot(gs[0, 3:])
    overlap = fig.add_subplot(gs[1, :2])
    registry = fig.add_subplot(gs[1, 2:4])
    support = fig.add_subplot(gs[1, 4:])
    draw_pair(pair0, d, frame, 0)
    draw_pair(pair1, d, frame, 1)
    time = d["time"]
    for pair, label in enumerate(("A0–A1", "A2–A3")):
        color = PAIR_COLORS[pair]
        overlap.plot(time, 1e3 * d["overlap"][:, pair], color=color, lw=2.4, label=label)
        registry.plot(time, 1e3 * d["registry"][:, pair], color=color, lw=2.4, label=label)
        support.plot(time, d["support"][:, pair], color=color, lw=2.4, label=label)
    overlap.axhline(100, color="#555", ls="--", lw=1.2, label="preferred 100 nm")
    overlap.set_ylabel("A–A overlap (nm)")
    overlap.set_ylim(-5, 175)
    registry.axhspan(0, 40, color="#dff3e8", alpha=0.9, label="locking window ≤40 nm")
    registry.axhline(40, color="#555", ls="--", lw=1.2)
    registry.set_ylabel("registry error (nm)")
    registry.set_ylim(0, 145)
    support.axhline(0.25, color="#555", ls="--", lw=1.2, label="support gate 0.25 pN")
    support.set_ylabel("productive opposed support (pN)")
    support.set_ylim(-0.1, max(5.6, np.nanmax(d["support"]) * 1.08))
    a3 = support.twinx()
    a3.plot(time, d["strong"][:, 3], color=STRONG_COLOR, lw=1.8,
            alpha=0.75, label="A3 strong heads")
    a3.set_ylim(-0.2, max(7, np.max(d["strong"][:, 3]) + 1))
    a3.set_ylabel("A3 strong heads", color=STRONG_COLOR)
    for ax in (overlap, registry, support):
        ax.axvline(time[frame], color="#111", lw=1.4)
        ax.set_xlim(time[0], time[-1])
        ax.set_xlabel("time (ms)")
        ax.grid(alpha=0.18)
        ax.legend(loc="best", fontsize=7.5)
    if time[frame] < 33.5:
        message = "A2 has motor force; A3 still has no force-producing head → no opposed support"
    elif time[frame] < 67:
        message = "A2–A3 has left the 40 nm register window; status-1 links do not restore axial register"
    else:
        message = "A3 finally develops force and opposed support appears—but the pair is already out of register"
    fig.suptitle(f"Balanced-gate failure replay — t = {time[frame]:.0f} ms\n{message}",
                 fontsize=15, weight="bold", y=0.98)
    fig.text(0.015, 0.012,
             "orange connectors: ACTN2 status 1  |  green: status 2  |  purple arrows: active motor force  |  open circles: barbed ends",
             fontsize=9, color="#4a4f55")
    fig.subplots_adjust(left=0.06, right=0.94, bottom=0.10, top=0.86)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trajectory", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--fps", type=int, default=12)
    args = parser.parse_args()
    d = load(args.trajectory)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(16, 9), facecolor="white")
    writer = FFMpegWriter(
        fps=args.fps, codec="libx264",
        extra_args=("-pix_fmt", "yuv420p", "-crf", "18", "-threads", "1"),
    )
    with writer.saving(fig, str(args.output), dpi=120):
        for frame in range(len(d["time"])):
            render(fig, d, frame)
            writer.grab_frame(facecolor="white")
            print(f"rendered {frame + 1}/{len(d['time'])}", flush=True)
    plt.close(fig)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
