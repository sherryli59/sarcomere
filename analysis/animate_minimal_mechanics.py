#!/usr/bin/env python3
"""Render sparse AM/ACTN2 mechanics trajectories with diagnostic overlays."""

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


ACTIN_COLORS = ("#ef5a29", "#d7308f")
ACTIN_RING_COLORS = ("#0072B2", "#D55E00", "#56B4E9", "#E69F00")
MYOSIN_COLOR = "#2171b5"
BARE_COLOR = "#d9e2ec"
STATUS1_COLOR = "#f0ad00"
STATUS2_COLOR = "#18a558"
WEAK_COLOR = "#8e6bbf"
STRONG_COLOR = "#5e3c99"
TITIN_COLOR = "#00a6b2"
MYOMESIN_COLOR = "#d81b60"


def endpoints(center: np.ndarray, direction: np.ndarray, length: float):
    return center - 0.5 * length * direction, center + 0.5 * length * direction


def closest_segment_points(p0, p1, q0, q1):
    """Closest points for two finite, nonperiodic line segments."""
    u, v, w = p1 - p0, q1 - q0, p0 - q0
    a, b, c = np.dot(u, u), np.dot(u, v), np.dot(v, v)
    d, e = np.dot(u, w), np.dot(v, w)
    denominator = a * c - b * b
    s = 0.0 if denominator < 1e-14 else np.clip((b * e - c * d) / denominator, 0, 1)
    t = np.clip((b * s + e) / max(c, 1e-14), 0, 1)
    s = np.clip((b * t - d) / max(a, 1e-14), 0, 1)
    return p0 + s * u, q0 + t * v


def load(path: Path, save_every: int):
    with h5py.File(path, "r", locking=False) as h:
        data = {
            "ac": h["actin/center"][:],
            "ad": h["actin/direction"][:],
            "mc": h["myosin/center"][:],
            "md": h["myosin/direction"][:],
            "status": h["actin/cb_status"][:, :, 0],
            "weak": h["actin/am_weak_count"][:, :, 0],
            "strong": h["actin/am_strong_count"][:, :, 0],
            "calcium": h["actin/am_calcium"][:, :, 0],
            "q": h["actin/am_q"][:, :, 0],
            "motor_force": h["actin/am_motor_force"][:, :, 0],
            "am_bonds": h["actin_myo/bonds"][:],
            "dt": float(h["meta"].attrs["dt"]),
            # New trajectories carry the simulator constant explicitly. The
            # fallback is the production value for older diagnostic files.
            "bare_fraction": float(
                h["meta"].attrs.get("am_bare_zone_fraction", 0.10)
            ),
            "actin_length": float(h["meta"].attrs.get("actin_length", 1.0)),
            "myosin_length": float(h["meta"].attrs.get("myosin_length", 1.5)),
            "box": np.array(
                [
                    float(h["meta"].attrs.get("box_x", 0.0)),
                    float(h["meta"].attrs.get("box_y", 0.0)),
                    float(h["meta"].attrs.get("box_z", 0.0)),
                ],
                dtype=float,
            ),
            "periodic": np.array(
                [
                    bool(h["meta"].attrs.get("periodic_x", 0)),
                    bool(h["meta"].attrs.get("periodic_y", 0)),
                    bool(h["meta"].attrs.get("periodic_z", 0)),
                ],
                dtype=bool,
            ),
            "myosin_radius": float(h["meta"].attrs.get("myosin_radius", 0.025)),
        }
        # Conserved head states are the authoritative AM population.  The
        # per-actin feature counters are force-pass diagnostics and are still
        # zero in the initialization frame, even when minimal_catch() has
        # already seeded weak heads.  Reconstructing by bound-actin ownership
        # keeps the first movie frame honest and also catches counter drift.
        if ("state/am_head_state_current" in h and
                "state/am_head_bound_actin_current" in h):
            head_state = h["state/am_head_state_current"][:]
            head_actin = h["state/am_head_bound_actin_current"][:]
            data["head_state"] = head_state
            data["head_actin"] = head_actin
            for key, path in (
                ("head_target", "state/am_head_target_axial_current"),
                ("head_azimuth", "state/am_head_binding_azimuth_current"),
                ("head_stroke", "state/am_head_stroke_progress_current"),
                ("material_e1", "state/myosin_material_e1_current"),
                ("material_e2", "state/myosin_material_e2_current"),
            ):
                if path in h:
                    data[key] = h[path][:]
            weak = np.zeros_like(data["weak"])
            strong = np.zeros_like(data["strong"])
            for actin_index in range(weak.shape[1]):
                owned = head_actin == actin_index
                weak[:, actin_index] = np.count_nonzero(
                    owned & (head_state == 1), axis=(1, 2)
                )
                strong[:, actin_index] = np.count_nonzero(
                    owned & (head_state == 2), axis=(1, 2)
                )
            data["weak"] = weak
            data["strong"] = strong
        if "catch_bond/local_link_frame_offsets" in h:
            data["link_offsets"] = h["catch_bond/local_link_frame_offsets"][:]
            data["link_identity"] = h["catch_bond/local_link_identity"][:]
            data["link_coordinates"] = h["catch_bond/local_link_coordinates"][:]
            if "catch_bond/local_link_mechanics" in h:
                data["link_mechanics"] = h["catch_bond/local_link_mechanics"][:]
        if "state/titin_bound_actin_current" in h:
            data["titin_bound_actin"] = h["state/titin_bound_actin_current"][:]
            data["titin_force"] = h["state/titin_force_current"][:]
            data["titin_links_per_half"] = int(
                h["meta"].attrs.get("titin_links_per_half", 1)
            )
        if "myomesin/frame_offsets" in h:
            data["myomesin_offsets"] = h["myomesin/frame_offsets"][:]
            data["myomesin_state"] = h["myomesin/bridge_state"][:]
            if "myomesin/events" in h:
                data["myomesin_events"] = h["myomesin/events"][:]
    if "link_offsets" in data:
        data["time_ms"] = data["link_offsets"][:, 0] * data["dt"] * 1e3
    else:
        data["time_ms"] = np.arange(len(data["ac"])) * save_every * data["dt"] * 1e3
    data["aa1"] = np.zeros(len(data["ac"]), dtype=int)
    data["aa2"] = np.zeros(len(data["ac"]), dtype=int)
    data["bilateral_load"] = np.zeros(len(data["ac"]), dtype=float)
    data["productive_support"] = np.zeros(len(data["ac"]), dtype=float)
    if "link_offsets" in data:
        for frame, (_, start, count) in enumerate(data["link_offsets"]):
            start, count = int(start), int(count)
            statuses = data["link_identity"][start:start + count, 4]
            data["aa1"][frame] = np.count_nonzero(statuses == 1)
            data["aa2"][frame] = np.count_nonzero(statuses == 2)
            if count and "link_mechanics" in data:
                mechanics = data["link_mechanics"][start:start + count]
                data["bilateral_load"][frame] = np.max(mechanics[:, 7])
                data["productive_support"][frame] = np.max(mechanics[:, 8])
    if "titin_bound_actin" in data:
        data["titin_bound_count"] = np.count_nonzero(
            data["titin_bound_actin"] >= 0, axis=1
        )
        data["titin_total_force"] = np.sum(data["titin_force"], axis=1)
    if "myomesin_offsets" in data:
        data["myomesin_count"] = data["myomesin_offsets"][:, 2].astype(int)
        data["myomesin_max_force"] = np.zeros(len(data["ac"]), dtype=float)
        data["myomesin_max_axial_nm"] = np.zeros(len(data["ac"]), dtype=float)
        pairs = sorted({
            tuple(sorted((int(row[2]), int(row[3]))))
            for row in data["myomesin_state"]
        })
        data["myomesin_pair_axial_nm"] = {
            pair: np.full(len(data["ac"]), np.nan) for pair in pairs
        }
        for frame, (_, start, count) in enumerate(data["myomesin_offsets"]):
            rows = data["myomesin_state"][int(start):int(start + count)]
            if len(rows):
                data["myomesin_max_force"][frame] = np.max(rows[:, 7])
                data["myomesin_max_axial_nm"][frame] = 1e3 * np.max(
                    np.abs(rows[:, 5])
                )
                for row in rows:
                    pair = tuple(sorted((int(row[2]), int(row[3]))))
                    data["myomesin_pair_axial_nm"][pair][frame] = (
                        1e3 * abs(float(row[5]))
                    )
        data["myomesin_formation_axial_nm"] = {}
        for row in data.get("myomesin_events", []):
            if int(row[1]) != 1:
                continue
            pair = tuple(sorted((int(row[3]), int(row[4]))))
            data["myomesin_formation_axial_nm"].setdefault(
                pair, 1e3 * abs(float(row[5]))
            )
    if data["ac"].shape[1] == 4 and data["mc"].shape[1] == 2:
        owners = np.array((0, 1, 0, 1), dtype=int)
        am_dot = np.stack(
            [np.sum(data["ad"][:, a] * data["md"][:, owners[a]], axis=1)
             for a in range(4)],
            axis=1,
        )
        data["max_am_angle"] = np.max(
            np.degrees(np.arccos(np.clip(np.abs(am_dot), -1.0, 1.0))),
            axis=1,
        )
        aa_dot = np.stack(
            [np.sum(data["ad"][:, a] * data["ad"][:, b], axis=1)
             for a, b in ((0, 1), (2, 3))],
            axis=1,
        )
        data["max_aa_angle_error"] = np.max(
            np.degrees(np.arccos(np.clip(-aa_dot, -1.0, 1.0))),
            axis=1,
        )
        mm_dot = np.sum(data["md"][:, 0] * data["md"][:, 1], axis=1)
        data["myosin_angle"] = np.degrees(
            np.arccos(np.clip(np.abs(mm_dot), -1.0, 1.0))
        )
    return data


def plot_rod(ax, center, direction, length, color, width, alpha=1.0):
    left, right = endpoints(center, direction, length)
    ax.plot(*np.vstack((left, right)).T, color=color, lw=width,
            alpha=alpha, solid_capstyle="round")


def minimum_image(displacement, box, periodic):
    """Return the nearest periodic image of a displacement vector."""
    result = np.asarray(displacement, dtype=float).copy()
    valid = periodic & (box > 0)
    result[valid] -= box[valid] * np.rint(result[valid] / box[valid])
    return result


def plot_myosin(ax, center, direction, length, bare_fraction):
    # Read the simulator's head-free fraction from trajectory metadata rather
    # than maintaining an independent visualization constant.
    bare_half = 0.5 * bare_fraction * length
    outer = 0.5 * length
    for lo, hi in ((-outer, -bare_half), (bare_half, outer)):
        points = np.vstack((center + lo * direction, center + hi * direction))
        ax.plot(*points.T, color=MYOSIN_COLOR, lw=12, solid_capstyle="round")
    bare = np.vstack((center - bare_half * direction, center + bare_half * direction))
    ax.plot(*bare.T, color=BARE_COLOR, lw=12, solid_capstyle="butt")
    ax.scatter(*center, color="#16324f", s=18, depthshade=False)


def active_local_links(data, frame):
    if "link_offsets" not in data:
        return []
    _, start, count = data["link_offsets"][frame]
    result = []
    for row in range(int(start), int(start + count)):
        ident = data["link_identity"][row]
        result.append((int(ident[2]), int(ident[3]), int(ident[4]),
                       *data["link_coordinates"][row]))
    return result


def explicit_heads(data, frame):
    required = ("head_state", "head_actin", "head_target", "head_azimuth",
                "material_e1", "material_e2")
    if any(key not in data for key in required):
        return None
    states = data["head_state"][frame]
    bound = data["head_actin"][frame]
    targets = data["head_target"][frame]
    azimuths = data["head_azimuth"][frame]
    e1 = data["material_e1"][frame]
    e2 = data["material_e2"][frame]
    heads_per_myosin = states.shape[1]
    heads_per_half = heads_per_myosin // 2
    bare_half = 0.5 * data["bare_fraction"] * data["myosin_length"]
    head_region = 0.5 * (1.0 - data["bare_fraction"]) * data["myosin_length"]
    crown_count = max(1, (heads_per_half + 2) // 3)
    anchors, sites, visible_states, strains = [], [], [], []
    anchor_states = []
    for myosin_index in range(states.shape[0]):
        for head_index in range(heads_per_myosin):
            half_index, local_head = divmod(head_index, heads_per_half)
            half_sign = -1.0 if half_index == 0 else 1.0
            crown, spoke = divmod(local_head, 3)
            fraction = (crown + 0.5) / crown_count
            axial = half_sign * (bare_half + fraction * head_region)
            state = int(states[myosin_index, head_index])
            anchor_angle = 2.0 * np.pi * spoke / 3.0 + crown * np.deg2rad(40.0)
            angle = (azimuths[myosin_index, head_index]
                     if state in (1, 2) else anchor_angle)
            anchor = (
                data["mc"][frame, myosin_index]
                + axial * data["md"][frame, myosin_index]
                + data["myosin_radius"] * (
                    np.cos(angle) * e1[myosin_index]
                    + np.sin(angle) * e2[myosin_index]
                )
            )
            anchors.append(anchor)
            anchor_states.append(state)
            actin_index = int(bound[myosin_index, head_index])
            if state not in (1, 2) or not 0 <= actin_index < len(data["ac"][frame]):
                continue
            raw_site = (
                data["ac"][frame, actin_index]
                + targets[myosin_index, head_index]
                * data["ad"][frame, actin_index]
            )
            site = anchor + minimum_image(
                raw_site - anchor, data["box"], data["periodic"]
            )
            sites.append(site)
            visible_states.append(state)
            strains.append(np.linalg.norm(site - anchor))
    return (np.asarray(anchors), np.asarray(anchor_states, dtype=int),
            np.asarray(sites).reshape((-1, 3)),
            np.asarray(visible_states, dtype=int), np.asarray(strains))


def explicit_titins(data, frame):
    required = ("titin_bound_actin", "titin_force", "material_e1",
                "material_e2", "titin_links_per_half")
    if any(key not in data for key in required):
        return None
    per_half = data["titin_links_per_half"]
    bound = data["titin_bound_actin"][frame]
    forces = data["titin_force"][frame]
    anchors, sites, visible_forces = [], [], []
    for index, actin_index_value in enumerate(bound):
        actin_index = int(actin_index_value)
        if not 0 <= actin_index < len(data["ac"][frame]):
            continue
        myosin_index, local = divmod(index, 2 * per_half)
        half, track = divmod(local, per_half)
        if myosin_index >= len(data["mc"][frame]):
            continue
        axis = data["md"][frame, myosin_index]
        phase = np.pi / per_half if per_half > 1 and half == 0 else 0.0
        azimuth = phase + 2.0 * np.pi * track / per_half
        radial = (
            np.cos(azimuth) * data["material_e1"][frame, myosin_index]
            + np.sin(azimuth) * data["material_e2"][frame, myosin_index]
        )
        half_sign = -1.0 if half == 0 else 1.0
        anchor = (
            data["mc"][frame, myosin_index]
            + half_sign * 0.5 * data["myosin_length"] * axis
            + data["myosin_radius"] * radial
        )
        raw_site = (
            data["ac"][frame, actin_index]
            - 0.5 * data["actin_length"] * data["ad"][frame, actin_index]
        )
        site = anchor + minimum_image(
            raw_site - anchor, data["box"], data["periodic"]
        )
        anchors.append(anchor)
        sites.append(site)
        visible_forces.append(float(forces[index]))
    return (np.asarray(anchors).reshape((-1, 3)),
            np.asarray(sites).reshape((-1, 3)),
            np.asarray(visible_forces))


def active_myomesin_bridges(data, frame):
    if "myomesin_offsets" not in data:
        return []
    _, start, count = data["myomesin_offsets"][frame]
    result = []
    for row in data["myomesin_state"][int(start):int(start + count)]:
        result.append({
            "i": int(row[2]),
            "j": int(row[3]),
            "axial_nm": 1e3 * abs(float(row[5])),
            "radial_nm": 1e3 * float(row[6]),
            "force": float(row[7]),
        })
    return result


def configure_geometry(ax, data, mode):
    ax.set_proj_type("ortho")
    if mode == "catch":
        balanced = data["ac"].shape[1] > 2
        ax.set_xlim(-2.85 if balanced else -2.15,
                    2.85 if balanced else 2.15)
        dense_lattice = data["ac"].shape[1] >= 24
        transverse_limit = 0.16 if dense_lattice or balanced else 0.42
        ax.set_ylim(-transverse_limit, transverse_limit)
        ax.set_zlim(-transverse_limit, transverse_limit)
        ax.view_init(elev=20, azim=-75)
    else:
        ax.set_xlim(-1.10, 1.25)
        ax.set_ylim(-0.45, 0.45)
        ax.set_zlim(-0.45, 0.45)
        ax.view_init(elev=23, azim=-67)
    # Keep axial distances physical and use only a small, explicitly labelled
    # transverse magnification. The former fixed (2.6, 1, 1) box exaggerated
    # transverse coordinates by roughly threefold in the balanced cell and
    # made a 2-degree AM error look like a large kink.
    xlim, ylim, zlim = ax.get_xlim(), ax.get_ylim(), ax.get_zlim()
    transverse_display_scale = 1.75 if mode == "catch" and balanced else 1.0
    ax.set_box_aspect((xlim[1] - xlim[0],
                       transverse_display_scale * (ylim[1] - ylim[0]),
                       transverse_display_scale * (zlim[1] - zlim[0])))
    ax.set_xlabel("x (µm)")
    if transverse_display_scale > 1.0:
        ax.set_yticks([])
        ax.set_zticks([])
        ax.set_ylabel("")
        ax.set_zlabel("")
        ax.text2D(0.01, 0.96, "transverse display ×1.75",
                  transform=ax.transAxes, fontsize=8, color="#555555")
    else:
        ax.set_ylabel("y (µm)")
        ax.set_zlabel("z (µm)")
    ax.grid(False)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_alpha(0.03)


def render_frame(fig, data, frame, mode, actin_length, myosin_length):
    fig.clear()
    grid = fig.add_gridspec(3, 2, width_ratios=(1.35, 1.0), hspace=0.42)
    ax = fig.add_subplot(grid[:, 0], projection="3d")
    metrics = fig.add_subplot(grid[0, 1])
    activation = fig.add_subplot(grid[1, 1])
    counts = fig.add_subplot(grid[2, 1])
    configure_geometry(ax, data, mode)
    time_ms = data["time_ms"]

    # Ghosted starting state makes nanometre-scale translation visible without
    # altering the actual rendered coordinates.
    dense_lattice = len(data["ac"][frame]) >= 24
    for i in range(len(data["ac"][0])):
        plot_rod(ax, data["ac"][0, i], data["ad"][0, i], actin_length,
                 "#a8a8a8", 2.5 if dense_lattice else 5, alpha=0.18)
    for m in range(len(data["mc"][frame])):
        plot_myosin(
            ax,
            data["mc"][frame, m],
            data["md"][frame, m],
            myosin_length,
            data["bare_fraction"],
        )
    for i in range(len(data["ac"][frame])):
        color = (ACTIN_RING_COLORS[min(i // 6, 3)]
                 if len(data["ac"][frame]) >= 24
                 else ACTIN_COLORS[i % len(ACTIN_COLORS)])
        plot_rod(ax, data["ac"][frame, i], data["ad"][frame, i], actin_length,
                 color, 4 if dense_lattice else 8)
        c, d = data["ac"][frame, i], data["ad"][frame, i]
        ax.quiver(*(c - 0.10 * d), *(0.20 * d), color=color,
                  arrow_length_ratio=0.35,
                  linewidth=1 if dense_lattice else 2)

    # Exact occupied M-band graph.  Each magenta connector joins the two
    # thick-filament bare-zone centers that consume a finite myomesin site.
    # Minimum-image placement keeps a periodic bridge local on screen.
    myomesin_bridges = active_myomesin_bridges(data, frame)
    bridge_labels = []
    for bridge in myomesin_bridges:
        first = data["mc"][frame, bridge["i"]]
        second = first + minimum_image(
            data["mc"][frame, bridge["j"]] - first,
            data["box"], data["periodic"],
        )
        midpoint = 0.5 * (first + second)
        ax.plot(*np.vstack((first, second)).T, color=MYOMESIN_COLOR,
                lw=6.0, alpha=0.95, solid_capstyle="round", zorder=20)
        ax.scatter(*first, facecolor="white", edgecolor=MYOMESIN_COLOR,
                   linewidth=2.0, s=58, depthshade=False, zorder=21)
        ax.scatter(*second, facecolor="white", edgecolor=MYOMESIN_COLOR,
                   linewidth=2.0, s=58, depthshade=False, zorder=21)
        ax.scatter(*midpoint, marker="D", color=MYOMESIN_COLOR,
                   s=46, depthshade=False, zorder=22)
        bridge_labels.append(
            f"M{bridge['i']}–M{bridge['j']}: Δx {bridge['axial_nm']:.0f} nm, "
            f"r⊥ {bridge['radial_nm']:.0f} nm, F {bridge['force']:.1f} pN"
        )
    if bridge_labels:
        ax.text2D(
            0.01, 0.885,
            "coarse myomesin M–M module: magenta; diamonds = module midpoints\n"
            + "\n".join(bridge_labels),
            transform=ax.transAxes, fontsize=8, color=MYOMESIN_COLOR,
        )

    # Explicit titin molecules connect a material domain at each thick-filament
    # tip to one individual actin barbed endpoint. Dashed links are bound but
    # slack; solid links are carrying tension. This is reconstructed from the
    # checkpointed molecular identity rather than from an owner heuristic.
    titins = explicit_titins(data, frame)
    if titins is not None:
        titin_anchors, titin_sites, titin_forces = titins
        for anchor, site, force in zip(titin_anchors, titin_sites, titin_forces):
            taut = force > 1e-9
            ax.plot(*np.vstack((anchor, site)).T, color=TITIN_COLOR,
                    lw=3.2 if taut else 2.0,
                    ls="-" if taut else "--", alpha=0.95 if taut else 0.65,
                    solid_capstyle="round")
            ax.scatter(*anchor, marker="^", color=TITIN_COLOR,
                       s=34, depthshade=False)
            ax.scatter(*site, marker="s", facecolor="white",
                       edgecolor=TITIN_COLOR, linewidth=1.5,
                       s=38, depthshade=False)

    # Every myosin head anchor is visible. Detached/off anchors are neutral;
    # weak and strong anchors inherit the connector color. Bound connectors
    # terminate at the exact saved actin material coordinate, so their length
    # is the actual geometric strain and the endpoint is the binding site.
    heads = explicit_heads(data, frame)
    if heads is not None:
        anchors, anchor_states, sites, site_states, strains = heads
        head_weak, head_strong = "#ffd166", "#ff5a36"
        anchor_palette = {0: "#9aa0a6", 1: head_weak,
                          2: head_strong, 3: "#d0d4d8"}
        for state, color in anchor_palette.items():
            selected = anchor_states == state
            if np.any(selected):
                ax.scatter(*anchors[selected].T, color=color,
                           s=8 if state in (0, 3) else 18,
                           alpha=0.45 if state in (0, 3) else 0.95,
                           depthshade=False)
        bound_index = 0
        for state in site_states:
            anchor_candidates = np.flatnonzero(
                np.isin(anchor_states, (1, 2)))
            anchor = anchors[anchor_candidates[bound_index]]
            site = sites[bound_index]
            color = head_strong if state == 2 else head_weak
            ax.plot(*np.vstack((anchor, site)).T, color=color,
                    lw=2.8 if state == 2 else 1.8,
                    alpha=0.95, solid_capstyle="round")
            ax.scatter(*site, facecolor="white", edgecolor=color,
                       linewidth=1.4, s=32, depthshade=False)
            bound_index += 1
        max_strain_nm = 1e3 * np.max(strains) if len(strains) else 0.0
        ax.text2D(
            0.01, 0.01,
            "heads: gray detached  |  yellow weak  |  red strong\n"
            "white-ring endpoint = actin binding site; connector length = strain\n"
            f"max bound connector length: {max_strain_nm:.1f} nm\n"
            + (f"titin: cyan, solid taut / dashed slack; "
               f"bound {len(titin_anchors)}/{len(data['titin_bound_actin'][frame])}; "
               f"tension {np.sum(titin_forces):.1f} pN"
               if titins is not None else ""),
            transform=ax.transAxes, fontsize=8, color="#333333",
        )

    # Render exact local ACTN2 endpoints and their current status.
    for i, j, status, si, sj in active_local_links(data, frame):
        pi = data["ac"][frame, i] + si * data["ad"][frame, i]
        pj = data["ac"][frame, j] + sj * data["ad"][frame, j]
        # Cross-boundary links are local in the minimum-image convention.
        pj = pi + minimum_image(
            pj - pi, data["box"], data["periodic"]
        )
        color = STATUS2_COLOR if status == 2 else STATUS1_COLOR
        ax.plot(*np.vstack((pi, pj)).T, color=color,
                lw=2.5 if dense_lattice else 6)
        ax.scatter(*((pi + pj) / 2), color=color,
                   s=12 if dense_lattice else 42, depthshade=False)

    # Older trajectories have no explicit head state. Keep the pair-level
    # schematic only as a backwards-compatible fallback.
    for pair in ([] if heads is not None else
                 data["am_bonds"][frame].reshape(-1, 2)):
        i, m = pair.astype(int)
        if i < 0 or m < 0:
            continue
        a0, a1 = endpoints(data["ac"][frame, i], data["ad"][frame, i], actin_length)
        # Use the nearest myosin image for local AM attachments.
        myosin_center = data["ac"][frame, i] + minimum_image(
            data["mc"][frame, m] - data["ac"][frame, i],
            data["box"],
            data["periodic"],
        )
        m0, m1 = endpoints(
            myosin_center, data["md"][frame, m], myosin_length
        )
        pa, pm = closest_segment_points(a0, a1, m0, m1)
        color = STRONG_COLOR if data["strong"][frame, i] > 0 else WEAK_COLOR
        ax.plot(*np.vstack((pa, pm)).T, color=color,
                lw=0.8 if dense_lattice else 2.5, ls="--",
                alpha=0.35 if dense_lattice else 1.0)

    if mode == "catch":
        metrics.plot(
            time_ms[:frame + 1], data["bilateral_load"][:frame + 1],
            color=STATUS2_COLOR, lw=2.5, label="max bilateral ACTN2 arm load",
        )
        metrics.plot(
            time_ms[:frame + 1], data["productive_support"][:frame + 1],
            color=STRONG_COLOR, lw=2.2, label="productive opposed AM support",
        )
        if "titin_total_force" in data:
            metrics.plot(
                time_ms[:frame + 1], data["titin_total_force"][:frame + 1],
                color=TITIN_COLOR, lw=2.0, label="total titin tension",
            )
        if "myomesin_max_force" in data:
            metrics.plot(
                time_ms[:frame + 1],
                data["myomesin_max_force"][:frame + 1],
                color=MYOMESIN_COLOR, lw=2.0,
                label="max myomesin module force",
            )
        # In the balanced periodic cell, actins 0/1 are the central link under
        # test and 2/3 supply the physical counterload at the image Z-disc.
        central_count = min(2, data["status"].shape[1])
        status = int(np.min(data["status"][frame, :central_count]))
        state = "status 2 — load stabilized" if status == 2 else "status 1 — encounter"
        geometry_name = (
            "Two-myosin / four-half / 24-actin periodic lattice"
            if data["ac"].shape[1] >= 24 else
            ("Shared-actin balanced cell: two myosins per actin"
             if (data["ac"].shape[1] == 4 and
                 data["mc"].shape[1] == 4) else
            ("Force-balanced periodic sarcomere cell"
             if data["ac"].shape[1] > 2 else
             ("One-sided motor + ACTN2-linked follower"
              if data["ac"].shape[1] == 2 else
              "Single actomyosin pair"))))
        ax.set_title(f"{geometry_name}\ncentral ACTN2: {state}", fontsize=13)
        metrics.set_ylabel("load (pN)")
        metrics.set_ylim(
            -1,
            max(5, np.max((data["bilateral_load"],
                           data["productive_support"],
                           data.get("titin_total_force",
                                    np.zeros_like(data["bilateral_load"])),
                           data.get("myomesin_max_force",
                                    np.zeros_like(data["bilateral_load"])))) * 1.12),
        )
    else:
        initial_axis = data["ad"][0, 0]
        displacement = 1e3 * np.sum(
            (data["ac"][:, 0] - data["ac"][0, 0]) * initial_axis, axis=1)
        alignment = np.abs(np.sum(data["ad"][:, 0] * data["md"][:, 0], axis=1))
        angle = np.degrees(np.arccos(np.clip(alignment, -1, 1)))
        metrics.plot(time_ms[:frame + 1], displacement[:frame + 1], color=ACTIN_COLORS[0],
                     lw=2.5, label="M-line-directed actin displacement")
        angle_axis = metrics.twinx()
        angle_axis.plot(time_ms[:frame + 1], angle[:frame + 1], color=MYOSIN_COLOR,
                        lw=2.5, label="AM misalignment angle")
        angle_axis.set_ylabel("misalignment (degrees)", color=MYOSIN_COLOR)
        angle_axis.tick_params(axis="y", labelcolor=MYOSIN_COLOR)
        angle_axis.set_ylim(0, 21)
        weak = data["weak"][frame, 0] > 0
        strong = data["strong"][frame, 0] > 0
        occupied = data["weak"][:, 0] + data["strong"][:, 0] > 0
        occupied_frames = np.flatnonzero(occupied)
        final_attached_frame = int(occupied_frames[-1]) if len(occupied_frames) else -1
        permanently_detached = frame > final_attached_frame
        state = "strong, pulling" if strong else (
            "weakly attached" if weak else (
                "detached — productive overlap exhausted"
                if permanently_detached else "between head attachments"
            )
        )
        ax.set_title(
            f"One myosin + one actin with physical rod diffusion\nAM state: {state}",
            fontsize=13,
        )
        metrics.set_ylabel("M-line-directed actin displacement (nm)", color=ACTIN_COLORS[0])
        metrics.tick_params(axis="y", labelcolor=ACTIN_COLORS[0])
        metrics.set_ylim(-1, max(12, np.max(displacement) * 1.15))
        if final_attached_frame + 1 < len(time_ms):
            detach_time = time_ms[final_attached_frame + 1]
            metrics.axvline(detach_time, color="#555555", ls=":", lw=1.5)
            metrics.axvspan(detach_time, time_ms[-1], color="#777777", alpha=0.08)
            metrics.text(
                detach_time, 0.97, "permanent detachment",
                transform=metrics.get_xaxis_transform(), ha="left", va="top",
                color="#555555", fontsize=9,
            )
        lines = [
            line for line in metrics.lines + angle_axis.lines
            if not line.get_label().startswith("_")
        ]
        metrics.legend(lines, [line.get_label() for line in lines], loc="lower right")

    metrics.set_xlabel("time (ms)")
    metrics.set_xlim(time_ms[0], time_ms[-1])
    metrics.grid(alpha=0.22)
    if mode == "catch":
        metrics.legend(loc="upper left")
    if mode == "catch" and "myomesin_pair_axial_nm" in data:
        pair_colors = (MYOMESIN_COLOR, "#7b1fa2", "#ef6c00", "#00838f")
        upper = 0.0
        for index, (pair, values) in enumerate(
                data["myomesin_pair_axial_nm"].items()):
            color = pair_colors[index % len(pair_colors)]
            activation.plot(
                time_ms[:frame + 1], values[:frame + 1],
                color=color, lw=2.4,
                label=f"M{pair[0]}–M{pair[1]} mismatch",
            )
            finite = values[np.isfinite(values)]
            if len(finite):
                upper = max(upper, float(np.max(finite)))
            formation = data["myomesin_formation_axial_nm"].get(pair)
            if formation is not None:
                activation.axhline(
                    formation, color=color, ls=":", lw=1.2, alpha=0.75,
                    label=f"M{pair[0]}–M{pair[1]} at association",
                )
                upper = max(upper, formation)
        activation.set_ylim(0.0, max(10.0, 1.10 * upper))
        activation.set_ylabel("M–M axial mismatch (nm)")
    elif mode == "catch" and "max_am_angle" in data:
        activation.plot(
            time_ms[:frame + 1], data["max_am_angle"][:frame + 1],
            color=STRONG_COLOR, lw=2.2, label="max A–M angle",
        )
        activation.plot(
            time_ms[:frame + 1], data["max_aa_angle_error"][:frame + 1],
            color=STATUS2_COLOR, lw=2.2, label="max A–A antiparallel error",
        )
        activation.plot(
            time_ms[:frame + 1], data["myosin_angle"][:frame + 1],
            color="#777777", lw=1.8, label="M–M angle",
        )
        angle_max = max(
            5.0,
            float(np.max((data["max_am_angle"],
                          data["max_aa_angle_error"],
                          data["myosin_angle"]))) * 1.12,
        )
        activation.set_ylim(-0.1, angle_max)
        activation.set_ylabel("misalignment (degrees)")
    else:
        activation.plot(
            time_ms[:frame + 1],
            np.mean(data["calcium"], axis=1)[:frame + 1],
            color="#777777", lw=2.2, label="calcium",
        )
        activation.plot(
            time_ms[:frame + 1], np.mean(data["q"], axis=1)[:frame + 1],
            color="#222222", lw=2.2, label="mean thin activation q",
        )
        activation.set_ylim(-0.02, 1.02)
        activation.set_ylabel("activation")
    activation.set_xlim(time_ms[0], time_ms[-1])
    activation.grid(alpha=0.22)
    activation.legend(loc="upper right", fontsize=8)
    total_weak = np.sum(data["weak"], axis=1)
    total_strong = np.sum(data["strong"], axis=1)
    counts.plot(time_ms[:frame + 1], total_weak[:frame + 1],
                color=WEAK_COLOR, lw=2.2, label="weak AM heads")
    counts.plot(time_ms[:frame + 1], total_strong[:frame + 1],
                color=STRONG_COLOR, lw=2.2, label="strong AM heads")
    counts.plot(time_ms[:frame + 1], data["aa1"][:frame + 1],
                color=STATUS1_COLOR, lw=2.2, label="status-1 ACTN2")
    counts.plot(time_ms[:frame + 1], data["aa2"][:frame + 1],
                color=STATUS2_COLOR, lw=2.2, label="status-2 ACTN2")
    if "titin_bound_count" in data:
        counts.plot(time_ms[:frame + 1], data["titin_bound_count"][:frame + 1],
                    color=TITIN_COLOR, lw=2.2, label="bound titin")
    if "myomesin_count" in data:
        counts.plot(time_ms[:frame + 1], data["myomesin_count"][:frame + 1],
                    color=MYOMESIN_COLOR, lw=2.2,
                    label="bound M–M modules")
    counts.set_xlim(time_ms[0], time_ms[-1])
    counts.set_ylim(bottom=-0.1)
    counts.set_xlabel("time (ms)")
    counts.set_ylabel("number attached")
    counts.grid(alpha=0.22)
    counts.legend(loc="upper right", fontsize=8, ncol=2)
    fig.suptitle(f"t = {time_ms[frame]:.2f} ms", fontsize=15, weight="bold")
    fig.subplots_adjust(left=0.03, right=0.97, bottom=0.10, top=0.88, wspace=0.22)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trajectory", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--mode", choices=("catch", "single"), required=True)
    parser.add_argument("--save-every", type=int, required=True)
    parser.add_argument("--actin-length", type=float)
    parser.add_argument("--myosin-length", type=float)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--stop", type=int)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--fps", type=int, default=15)
    args = parser.parse_args()

    data = load(args.trajectory, args.save_every)
    actin_length = args.actin_length if args.actin_length is not None else data["actin_length"]
    myosin_length = args.myosin_length if args.myosin_length is not None else data["myosin_length"]
    stop = min(args.stop if args.stop is not None else len(data["ac"]), len(data["ac"]))
    frames = list(range(args.start, stop, max(1, args.stride)))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    # At 130 dpi this is exactly 1820x832; yuv420p requires even dimensions.
    fig = plt.figure(figsize=(14, 7.6), facecolor="white")
    writer = FFMpegWriter(fps=args.fps, codec="libx264",
                          extra_args=("-pix_fmt", "yuv420p", "-crf", "18",
                                      "-threads", "1"))
    with writer.saving(fig, str(args.output), dpi=130):
        for number, frame in enumerate(frames, 1):
            render_frame(fig, data, frame, args.mode,
                         actin_length, myosin_length)
            writer.grab_frame(facecolor="white")
            print(f"rendered {number}/{len(frames)}", flush=True)
    plt.close(fig)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
