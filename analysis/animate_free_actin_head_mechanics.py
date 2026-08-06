#!/usr/bin/env python3
"""Render a one-actin/one-myosin trajectory with explicit head mechanics."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/sarcomere_matplotlib")

import h5py
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
import numpy as np


DETACHED, WEAK, STRONG, OFF = 0, 1, 2, 3
MYOSIN = "#2474b5"
BARE = "#d9e2ec"
ACTIN = "#ef6c32"
WEAK_C = "#f2bd31"
STRONG_C = "#df3b2f"
TETHER = "#00a8b5"
NET = "#27272a"
GHOST = "#a1a1aa"


def minimum_image(value, box, periodic):
    result = np.asarray(value, dtype=float).copy()
    valid = periodic & (box > 0)
    result[valid] -= box[valid] * np.rint(result[valid] / box[valid])
    return result


def load(path: Path):
    with h5py.File(path, "r", locking=False) as h:
        meta = h["/meta"].attrs
        data = {
            "time_ms": h["/state/current_step"][:, 0] * float(meta["dt"]) * 1e3,
            "ac": h["/actin/center"][:],
            "ad": h["/actin/direction"][:],
            "mc": h["/myosin/center"][:],
            "md": h["/myosin/direction"][:],
            "actin_force": h["/actin/force"][:, 0, :],
            "myosin_force": h["/myosin/force"][:, 0, :],
            "state": h["/state/am_head_state_current"][:, 0, :],
            "bound": h["/state/am_head_bound_actin_current"][:, 0, :],
            "target": h["/state/am_head_target_axial_current"][:, 0, :],
            "azimuth": h["/state/am_head_binding_azimuth_current"][:, 0, :],
            "stroke": h["/state/am_head_stroke_progress_current"][:, 0, :],
            "e1": h["/state/myosin_material_e1_current"][:, 0, :],
            "e2": h["/state/myosin_material_e2_current"][:, 0, :],
            "dt": float(meta["dt"]),
            "actin_length": float(meta.get("actin_length", 1.0)),
            "myosin_length": float(meta.get("myosin_length", 1.5)),
            "myosin_radius": float(meta.get("myosin_radius", 0.025)),
            "actin_radius": float(meta.get("actin_radius", 0.015)),
            "surface_gap": float(meta.get("am_surface_gap", 0.015)),
            "bare_fraction": float(meta.get("am_bare_zone_fraction", 0.10)),
            "motor_k": float(meta.get("am_motor_k", 1000.0)),
            "weak_ratio": float(meta.get("am_weak_stiffness_ratio", 0.1)),
            "stall": float(meta.get("am_stall_force", 4.5)),
            "head_k": float(meta.get("am_head_transverse_k", 0.0)),
            "slack": float(meta.get("am_tether_slack", 0.002)),
            "align_k": float(meta.get("am_alignment_k", 0.0)),
            "box": np.array([meta.get("box_x", 0.0), meta.get("box_y", 0.0),
                             meta.get("box_z", 0.0)], dtype=float),
            "periodic": np.array([meta.get("periodic_x", 0),
                                  meta.get("periodic_y", 0),
                                  meta.get("periodic_z", 0)], dtype=bool),
        }
    return data


def anchor_geometry(data, frame):
    states = data["state"][frame]
    count = len(states)
    heads_per_half = count // 2
    crown_count = max(1, (heads_per_half + 2) // 3)
    anchors = np.zeros((count, 3))
    axial_positions = np.zeros(count)
    bare_half = 0.5 * data["bare_fraction"] * data["myosin_length"]
    head_region = 0.5 * (1.0 - data["bare_fraction"]) * data["myosin_length"]
    axis = data["md"][frame, 0]
    for head in range(count):
        half_index, local = divmod(head, heads_per_half)
        sign = -1.0 if half_index == 0 else 1.0
        crown, spoke = divmod(local, 3)
        axial = sign * (bare_half + (crown + 0.5) / crown_count * head_region)
        crown_angle = 2.0 * np.pi * spoke / 3.0 + crown * np.deg2rad(40.0)
        angle = data["azimuth"][frame, head] if states[head] in (WEAK, STRONG) else crown_angle
        radial = (np.cos(angle) * data["e1"][frame] +
                  np.sin(angle) * data["e2"][frame])
        anchors[head] = (data["mc"][frame, 0] + axial * axis +
                         data["myosin_radius"] * radial)
        axial_positions[head] = axial
    return anchors, axial_positions


def bound_mechanics(data, frame):
    anchors, axial_positions = anchor_geometry(data, frame)
    states = data["state"][frame]
    axis_a = data["ad"][frame, 0]
    preferred = data["actin_radius"] + data["surface_gap"]
    rows = []
    for head, state in enumerate(states):
        if state not in (WEAK, STRONG) or data["bound"][frame, head] != 0:
            continue
        raw_site = (data["ac"][frame, 0] +
                    data["target"][frame, head] * axis_a)
        site = anchors[head] + minimum_image(
            raw_site - anchors[head], data["box"], data["periodic"])
        connector = site - anchors[head]
        z = float(np.dot(connector, axis_a))
        transverse = connector - z * axis_a
        rho = float(np.linalg.norm(transverse))
        radial_extension = rho - preferred
        loaded = 0.0
        if radial_extension > data["slack"]:
            loaded = radial_extension - data["slack"]
        elif radial_extension < -data["slack"]:
            loaded = radial_extension + data["slack"]
        f_perp = np.zeros(3)
        if data["head_k"] > 0.0 and rho > 1e-12 and loaded != 0.0:
            f_perp = -data["head_k"] * loaded * transverse / rho
        if state == STRONG:
            strain = data["stroke"][frame, head] - z
            active = np.clip(data["motor_k"] * strain, 0.0, data["stall"])
            f_axial = active * axis_a
        else:
            strain = -z
            f_axial = -(data["weak_ratio"] * data["motor_k"] * z) * axis_a
        transverse_unit = transverse / rho if rho > 1e-12 else np.zeros(3)
        preferred_point = anchors[head] + z * axis_a + preferred * transverse_unit
        rows.append({
            "head": head, "state": int(state), "anchor": anchors[head],
            "site": site, "z": z, "rho": rho, "loaded": loaded,
            "f_perp": f_perp, "f_axial": f_axial,
            "preferred_point": preferred_point, "strain": strain,
            "target": data["target"][frame, head],
        })
    return anchors, axial_positions, rows


def precompute(data):
    frames = len(data["time_ms"])
    strong_force = np.zeros(frames)
    weak_force = np.zeros(frames)
    transverse_force = np.zeros(frames)
    max_tether_extension = np.zeros(frames)
    relative_slide = np.zeros(frames)
    angle = np.zeros(frames)
    initial_axis = data["ad"][0, 0]
    for frame in range(frames):
        _, _, rows = bound_mechanics(data, frame)
        for row in rows:
            axial = float(np.dot(row["f_axial"], data["ad"][frame, 0]))
            if row["state"] == STRONG:
                strong_force[frame] += axial
            else:
                weak_force[frame] += axial
            transverse_force[frame] += np.linalg.norm(row["f_perp"])
            max_tether_extension[frame] = max(
                max_tether_extension[frame], abs(row["loaded"]))
        if frame:
            da = minimum_image(data["ac"][frame, 0] - data["ac"][frame - 1, 0],
                               data["box"], data["periodic"])
            dm = minimum_image(data["mc"][frame, 0] - data["mc"][frame - 1, 0],
                               data["box"], data["periodic"])
            relative_slide[frame] = relative_slide[frame - 1] + np.dot(
                da - dm, data["ad"][frame, 0])
        angle[frame] = np.degrees(np.arccos(np.clip(
            abs(np.dot(data["ad"][frame, 0], data["md"][frame, 0])), 0.0, 1.0)))
    occupied = np.any(np.isin(data["state"], (WEAK, STRONG)), axis=1)
    occupied_frames = np.flatnonzero(occupied)
    last_bound_frame = int(occupied_frames[-1]) if len(occupied_frames) else -1
    return {
        "strong_force": strong_force,
        "weak_force": weak_force,
        "transverse_force": transverse_force,
        "max_tether_extension": 1e3 * max_tether_extension,
        "relative_slide": 1e3 * relative_slide,
        "angle": angle,
        "weak_count": np.sum(data["state"] == WEAK, axis=1),
        "strong_count": np.sum(data["state"] == STRONG, axis=1),
        "net_axial": np.sum(data["actin_force"] * data["ad"][:, 0, :], axis=1),
        "initial_axis": initial_axis,
        "last_bound_frame": last_bound_frame,
    }


def arrow3d(ax, start, force, color, scale, label=None, width=1.6):
    norm = float(np.linalg.norm(force))
    if norm < 1e-8:
        return
    vector = force * scale
    ax.quiver(*start, *vector, color=color, linewidth=width,
              arrow_length_ratio=0.28, normalize=False, label=label)


def draw_cross_section(ax, data, frame, selected):
    ax.set_aspect("equal")
    ax.set_xlim(-100, 100)
    ax.set_ylim(-100, 100)
    ax.set_xticks((-50, 0, 50))
    ax.set_yticks((-50, 0, 50))
    ax.tick_params(labelsize=7)
    ax.set_xlabel("local transverse 1 (nm)", fontsize=7)
    ax.set_ylabel("local transverse 2 (nm)", fontsize=7)
    ax.set_title("head cross-section", fontsize=8)
    ax.add_patch(Circle((0, 0), 1e3 * data["myosin_radius"],
                        facecolor=MYOSIN, edgecolor="none", alpha=0.85))
    axis_m = data["md"][frame, 0]
    delta = minimum_image(data["ac"][frame, 0] - data["mc"][frame, 0],
                          data["box"], data["periodic"])
    c1 = 1e3 * np.dot(delta, data["e1"][frame])
    c2 = 1e3 * np.dot(delta, data["e2"][frame])
    ax.add_patch(Circle((c1, c2), 1e3 * data["actin_radius"],
                        facecolor=ACTIN, edgecolor="none", alpha=0.9))
    if selected is not None:
        def project(point):
            rel = point - data["mc"][frame, 0]
            return 1e3 * np.array((np.dot(rel, data["e1"][frame]),
                                   np.dot(rel, data["e2"][frame])))
        anchor = project(selected["anchor"])
        site = project(selected["site"])
        preferred = project(selected["preferred_point"])
        color = STRONG_C if selected["state"] == STRONG else WEAK_C
        ax.plot(*np.vstack((anchor, site)).T, color=color, lw=2.0)
        ax.plot(*np.vstack((preferred, site)).T, color=TETHER, lw=4.0,
                alpha=0.8, solid_capstyle="round")
        ax.scatter(*anchor, color=color, s=22, zorder=4)
        ax.scatter(*site, facecolor="white", edgecolor=color, s=28, zorder=4)
        if np.linalg.norm(selected["f_perp"]) > 1e-8:
            transverse_force = np.array((
                np.dot(selected["f_perp"], data["e1"][frame]),
                np.dot(selected["f_perp"], data["e2"][frame])))
            ax.arrow(*site, *(8.0 * transverse_force), color=TETHER,
                     width=0.8, head_width=5, length_includes_head=True)
    ax.text(0.02, 0.02, "cyan = loaded part of passive tether",
            transform=ax.transAxes, fontsize=7, color="#444444")


def render(fig, data, traces, frame, end_frame):
    fig.clear()
    outer = fig.add_gridspec(4, 2, width_ratios=(1.6, 1.0),
                             height_ratios=(1.1, 0.9, 0.9, 0.8),
                             left=0.035, right=0.985, bottom=0.07, top=0.92,
                             wspace=0.24, hspace=0.46)
    ax = fig.add_subplot(outer[:, 0], projection="3d")
    walk = fig.add_subplot(outer[0, 1])
    motion = fig.add_subplot(outer[1, 1])
    forces = fig.add_subplot(outer[2, 1])
    counts = fig.add_subplot(outer[3, 1])
    time = data["time_ms"]
    now = time[frame]
    actin_length = data["actin_length"]
    myosin_length = data["myosin_length"]

    visible = slice(0, end_frame + 1)
    all_x = np.concatenate((
        data["ac"][visible, 0, 0] - 0.5 * actin_length,
        data["ac"][visible, 0, 0] + 0.5 * actin_length,
        data["mc"][visible, 0, 0] - 0.5 * myosin_length,
        data["mc"][visible, 0, 0] + 0.5 * myosin_length,
    ))
    ax.set_xlim(np.min(all_x) - 0.08, np.max(all_x) + 0.08)
    yz = 0.20
    ax.set_ylim(-yz, yz)
    ax.set_zlim(-yz, yz)
    ax.set_box_aspect((2.4, 1.0, 1.0))
    ax.set_proj_type("ortho")
    ax.view_init(elev=22, azim=-70)
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("y (µm)")
    ax.set_zlabel("z (µm)")
    ax.grid(False)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_alpha(0.025)

    # Starting positions remain as transparent references so recoil and actin
    # sliding are distinguishable without changing the simulated coordinates.
    for center, direction, length in (
        (data["ac"][0, 0], data["ad"][0, 0], actin_length),
        (data["mc"][0, 0], data["md"][0, 0], myosin_length),
    ):
        ends = np.vstack((center - 0.5 * length * direction,
                          center + 0.5 * length * direction))
        ax.plot(*ends.T, color=GHOST, lw=3.0, alpha=0.18)

    mc = data["mc"][frame, 0]
    md = data["md"][frame, 0]
    outer_half = 0.5 * myosin_length
    bare_half = 0.5 * data["bare_fraction"] * myosin_length
    for lo, hi in ((-outer_half, -bare_half), (bare_half, outer_half)):
        points = np.vstack((mc + lo * md, mc + hi * md))
        ax.plot(*points.T, color=MYOSIN, lw=11, solid_capstyle="round")
    bare = np.vstack((mc - bare_half * md, mc + bare_half * md))
    ax.plot(*bare.T, color=BARE, lw=11, solid_capstyle="butt")
    ax.scatter(*mc, color="#1f3349", s=22, depthshade=False)

    ac = data["ac"][frame, 0]
    ad = data["ad"][frame, 0]
    aends = np.vstack((ac - 0.5 * actin_length * ad,
                       ac + 0.5 * actin_length * ad))
    ax.plot(*aends.T, color=ACTIN, lw=8, solid_capstyle="round")
    ax.quiver(*(ac - 0.10 * ad), *(0.20 * ad), color=ACTIN,
              linewidth=2.2, arrow_length_ratio=0.35)

    anchors, _, rows = bound_mechanics(data, frame)
    states = data["state"][frame]
    detached = np.isin(states, (DETACHED, OFF))
    if np.any(detached):
        ax.scatter(*anchors[detached].T, color="#8d949b", s=8,
                   alpha=0.48, depthshade=False)
    selected = max(rows, key=lambda row: np.linalg.norm(row["f_perp"]),
                   default=None)
    for row in rows:
        color = STRONG_C if row["state"] == STRONG else WEAK_C
        ax.scatter(*row["anchor"], color=color, s=22, depthshade=False)
        ax.plot(*np.vstack((row["anchor"], row["site"])).T,
                color=color, lw=2.7 if row["state"] == STRONG else 1.8,
                ls="-" if row["state"] == STRONG else "--", alpha=0.95)
        ax.scatter(*row["site"], facecolor="white", edgecolor=color,
                   linewidth=1.4, s=35, depthshade=False)
        if abs(row["loaded"]) > 1e-12:
            ax.plot(*np.vstack((row["preferred_point"], row["site"])).T,
                    color=TETHER, lw=4.2, alpha=0.72,
                    solid_capstyle="round")
        arrow3d(ax, row["site"], row["f_axial"], STRONG_C, 0.018, width=1.5)
        arrow3d(ax, row["site"], row["f_perp"], TETHER, 0.018, width=1.5)
    net_force = data["actin_force"][frame]
    arrow3d(ax, ac, net_force, NET, 0.022, width=2.4)

    detail = "no bound heads"
    if selected is not None:
        label = "strong" if selected["state"] == STRONG else "weak"
        detail = (f"selected h{selected['head']}  {label}\n"
                  f"axial strain {1e3 * selected['strain']:+.1f} nm   "
                  f"tether excess {1e3 * selected['loaded']:+.1f} nm\n"
                  f"|F axial| {np.linalg.norm(selected['f_axial']):.2f} pN   "
                  f"|F tether| {np.linalg.norm(selected['f_perp']):.2f} pN")
    if traces["last_bound_frame"] >= 0 and frame > traces["last_bound_frame"]:
        detail = ("WALKED OFF PRODUCTIVE MYOSIN HALF\n"
                  f"no weak/strong heads since "
                  f"{time[traces['last_bound_frame']]:.1f} ms")
    ax.text2D(0.02, 0.02,
              "gray dots detached anchors · yellow dashed weak · red strong\n"
              "white rings binding sites · cyan loaded tether · black net force\n" +
              detail,
              transform=ax.transAxes, fontsize=8.5, color="#303038")
    ax.set_title("Free actin sliding on one free bipolar myosin\n"
                 f"all 70 head anchors · k⊥={data['head_k']:.0f} pN/µm · "
                 f"k_align={data['align_k']:.0f}", fontsize=12)

    cross = ax.inset_axes([0.015, 0.61, 0.28, 0.28])
    draw_cross_section(cross, data, frame, selected)

    # The head-target plot is the walk: each horizontal segment is one head
    # holding a fixed actin material site; gaps are detached intervals and a
    # new vertical level is the same head rebinding farther along the actin.
    for state, color, label in ((WEAK, WEAK_C, "weak"),
                                (STRONG, STRONG_C, "strong")):
        mask = (data["state"][:frame + 1] == state)
        ti, hi = np.nonzero(mask)
        if len(ti):
            walk.scatter(time[ti], 1e3 * data["target"][:frame + 1][ti, hi],
                         s=5, color=color, alpha=0.72, label=label)
    walk.axvline(now, color=NET, lw=1.0, alpha=0.7)
    walk.set_ylabel("bound actin site (nm)")
    walk.set_title("The walk: bind → hold site → detach → rebind", fontsize=9)
    handles, labels = walk.get_legend_handles_labels()
    if handles:
        walk.legend(handles, labels, loc="upper left", ncol=2,
                    fontsize=7, frameon=False)

    motion.plot(time[:frame + 1], traces["relative_slide"][:frame + 1],
                color=ACTIN, lw=2.1, label="relative slide")
    motion.set_ylabel("inward A−M slide (nm)", color=ACTIN)
    motion.tick_params(axis="y", labelcolor=ACTIN)
    angle_ax = motion.twinx()
    angle_ax.plot(time[:frame + 1], traces["angle"][:frame + 1],
                  color=MYOSIN, lw=1.8, label="A–M angle")
    angle_ax.set_ylabel("misalignment (°)", color=MYOSIN)
    angle_ax.tick_params(axis="y", labelcolor=MYOSIN)

    forces.plot(time[:frame + 1], traces["strong_force"][:frame + 1],
                color=STRONG_C, lw=1.8, label="strong axial")
    forces.plot(time[:frame + 1], traces["weak_force"][:frame + 1],
                color=WEAK_C, lw=1.5, label="weak axial")
    forces.plot(time[:frame + 1], traces["transverse_force"][:frame + 1],
                color=TETHER, lw=1.8, label="Σ |passive tether|")
    forces.plot(time[:frame + 1], traces["net_axial"][:frame + 1],
                color=NET, lw=1.3, alpha=0.8, label="net actin axial")
    forces.axhline(0.0, color=GHOST, lw=0.7)
    forces.set_ylabel("force (pN)")
    forces.legend(loc="upper left", ncol=2, fontsize=6.7, frameon=False)

    counts.step(time[:frame + 1], traces["weak_count"][:frame + 1],
                where="post", color=WEAK_C, lw=1.8, label="weak")
    counts.step(time[:frame + 1], traces["strong_count"][:frame + 1],
                where="post", color=STRONG_C, lw=1.8, label="strong")
    ext_ax = counts.twinx()
    ext_ax.plot(time[:frame + 1], traces["max_tether_extension"][:frame + 1],
                color=TETHER, lw=1.5, label="max tether excess")
    counts.set_ylabel("bound heads")
    ext_ax.set_ylabel("max tether excess (nm)", color=TETHER)
    ext_ax.tick_params(axis="y", labelcolor=TETHER)
    lines = counts.lines + ext_ax.lines
    counts.legend(lines, [line.get_label() for line in lines],
                  loc="upper left", ncol=3, fontsize=6.7, frameon=False)

    if 0 <= traces["last_bound_frame"] <= end_frame:
        detach_time = time[traces["last_bound_frame"]]
        for panel in (walk, motion, forces, counts):
            panel.axvline(detach_time, color="#6b7280", lw=1.0, ls=":")
            panel.axvspan(detach_time, time[end_frame], color="#d1d5db",
                         alpha=0.12, lw=0)
        counts.text(detach_time, 0.97, " sustained unbound",
                    transform=counts.get_xaxis_transform(), va="top",
                    fontsize=7, color="#4b5563")

    for panel in (walk, motion, forces, counts):
        panel.set_xlim(time[0], time[end_frame])
        panel.grid(alpha=0.18)
        panel.tick_params(labelsize=7)
    for panel in (walk, motion, forces):
        panel.tick_params(labelbottom=False)
    counts.set_xlabel("simulation time (ms)")
    fig.suptitle(f"t = {now:.1f} ms", fontsize=15, fontweight="bold")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trajectory", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--stride", type=int, default=4)
    parser.add_argument("--fps", type=int, default=25)
    parser.add_argument("--poster", type=Path)
    parser.add_argument("--stop-ms", type=float)
    args = parser.parse_args()
    data = load(args.trajectory)
    traces = precompute(data)
    end_frame = len(data["time_ms"]) - 1
    if args.stop_ms is not None:
        end_frame = min(end_frame, int(np.searchsorted(
            data["time_ms"], args.stop_ms, side="right") - 1))
        end_frame = max(0, end_frame)
    frames = list(range(0, end_frame + 1, max(1, args.stride)))
    if frames[-1] != end_frame:
        frames.append(end_frame)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(15.4, 8.7), facecolor="white")
    writer = FFMpegWriter(fps=args.fps, codec="libx264",
                          extra_args=("-pix_fmt", "yuv420p", "-crf", "19",
                                      "-threads", "1", "-movflags", "+faststart"))
    with writer.saving(fig, str(args.output), dpi=120):
        for index, frame in enumerate(frames, 1):
            render(fig, data, traces, frame, end_frame)
            writer.grab_frame(facecolor="white")
            if index % 25 == 0 or index == len(frames):
                print(f"rendered {index}/{len(frames)}", flush=True)
    if args.poster:
        render(fig, data, traces, end_frame, end_frame)
        fig.savefig(args.poster, dpi=160, facecolor="white")
    plt.close(fig)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
