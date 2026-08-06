#!/usr/bin/env python3
"""Streamed 3-D trajectory renderer for mechanical and KMC debugging.

The renderer understands the explicit AM-head and sparse local-ACTN2 schemas.
It intentionally uses batched VTK line geometry instead of one cylinder actor
per filament, which makes full 2000/1000 trajectories practical to inspect.
"""

from __future__ import annotations

import argparse
import ctypes.util
import os
import subprocess
import tempfile
from itertools import product
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/sarcomere_matplotlib")

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv


STATE_OFF, STATE_DETACHED, STATE_WEAK, STATE_STRONG = 3, 0, 1, 2


def line_mesh(starts: np.ndarray, ends: np.ndarray) -> pv.PolyData:
    starts = np.asarray(starts, dtype=float).reshape((-1, 3))
    ends = np.asarray(ends, dtype=float).reshape((-1, 3))
    if len(starts) == 0:
        return pv.PolyData()
    points = np.empty((2 * len(starts), 3), dtype=float)
    points[0::2], points[1::2] = starts, ends
    lines = np.column_stack(
        (np.full(len(starts), 2), 2 * np.arange(len(starts)), 2 * np.arange(len(starts)) + 1)
    ).ravel()
    return pv.PolyData(points, lines=lines)


def minimum_image(delta: np.ndarray, box: np.ndarray) -> np.ndarray:
    return delta - box * np.rint(delta / box)


def rods(centers: np.ndarray, directions: np.ndarray, length: float) -> tuple[np.ndarray, np.ndarray]:
    half = 0.5 * length * directions
    return centers - half, centers + half


def periodic_clip_segments(starts: np.ndarray, ends: np.ndarray,
                           box: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Split periodic segments into images clipped to the central box."""
    lower, upper = -.5 * box, .5 * box
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
    return (np.asarray(clipped_starts, dtype=float).reshape((-1, 3)),
            np.asarray(clipped_ends, dtype=float).reshape((-1, 3)))


def nematic(directions: np.ndarray) -> float:
    q = 1.5 * directions.T @ directions / max(1, len(directions)) - 0.5 * np.eye(3)
    return float(np.linalg.eigvalsh(q)[-1])


class Trajectory:
    def __init__(self, path: Path, box: np.ndarray,
                 actin_length: float | None, myosin_length: float | None):
        self.handle = h5py.File(path, "r", locking=False)
        self.box = box
        self.meta = self.handle["/meta"].attrs if "/meta" in self.handle else {}
        stored_actin = float(self.meta.get("actin_hydrodynamic_length", 0.0))
        stored_myosin = float(self.meta.get("myosin_hydrodynamic_length", 0.0))
        self.actin_length = actin_length if actin_length is not None else (
            stored_actin if stored_actin > 0 else 1.0)
        self.myosin_length = myosin_length if myosin_length is not None else (
            stored_myosin if stored_myosin > 0 else 1.5)
        required = ("/actin/center", "/actin/direction", "/myosin/center", "/myosin/direction")
        missing = [name for name in required if name not in self.handle]
        if missing:
            raise ValueError(f"missing trajectory datasets: {missing}")
        self.nframes = min(self.handle[name].shape[0] for name in required)
        if "/state/current_step" in self.handle:
            self.nframes = min(self.nframes, self.handle["/state/current_step"].shape[0])

    def close(self) -> None:
        self.handle.close()

    def diagnostic_series(self, frames: list[int]) -> dict[str, np.ndarray]:
        """Small global traces used beside every movie frame."""
        h = self.handle
        result = {
            "frame": np.asarray(frames, dtype=int),
            "time": np.zeros(len(frames)),
            "weak": np.zeros(len(frames)),
            "strong": np.zeros(len(frames)),
            "calcium": np.zeros(len(frames)),
            "q": np.zeros(len(frames)),
            "aa1": np.zeros(len(frames)),
            "aa2": np.zeros(len(frames)),
            "actin_nematic": np.zeros(len(frames)),
            "myosin_nematic": np.zeros(len(frames)),
        }
        dt = float(self.meta.get("dt", 1.0))
        offsets = h.get("/catch_bond/local_link_frame_offsets")
        identity = h.get("/catch_bond/local_link_identity")
        for output_index, frame in enumerate(frames):
            step = frame
            if "/state/current_step" in h:
                step = int(np.asarray(h["/state/current_step"][frame]).reshape(-1)[0])
            result["time"][output_index] = step * dt
            for key, path in (("weak", "/actin/am_weak_count"),
                              ("strong", "/actin/am_strong_count")):
                if path in h:
                    result[key][output_index] = np.sum(h[path][frame])
            for key, path in (("calcium", "/actin/am_calcium"),
                              ("q", "/actin/am_q")):
                if path in h:
                    result[key][output_index] = np.mean(h[path][frame])
            if offsets is not None and identity is not None and frame < len(offsets):
                _, start, count = np.asarray(offsets[frame], dtype=int)
                status = np.asarray(identity[start:start + count, 4], dtype=int)
                result["aa1"][output_index] = np.count_nonzero(status == 1)
                result["aa2"][output_index] = np.count_nonzero(status == 2)
            result["actin_nematic"][output_index] = nematic(
                np.asarray(h["/actin/direction"][frame], dtype=float))
            result["myosin_nematic"][output_index] = nematic(
                np.asarray(h["/myosin/direction"][frame], dtype=float))
        return result

    def frame(self, index: int) -> dict:
        h = self.handle
        data = {
            "actin_center": np.asarray(h["/actin/center"][index], dtype=float),
            "actin_direction": np.asarray(h["/actin/direction"][index], dtype=float),
            "myosin_center": np.asarray(h["/myosin/center"][index], dtype=float),
            "myosin_direction": np.asarray(h["/myosin/direction"][index], dtype=float),
        }
        for key, path in (
            ("actin_status", "/actin/cb_status"),
            ("actin_load", "/actin/f_load"),
            ("actin_force", "/actin/force"),
            ("myosin_force", "/myosin/force"),
            ("strong_count", "/actin/am_strong_count"),
            ("weak_count", "/actin/am_weak_count"),
        ):
            if path in h and index < h[path].shape[0]:
                data[key] = np.asarray(h[path][index]).squeeze()
        if "/myosin/bonds" in h and index < h["/myosin/bonds"].shape[0]:
            data["myosin_bonds"] = np.asarray(h["/myosin/bonds"][index], dtype=int)
        state_path = "/state/am_head_state_current"
        if state_path in h and index < h[state_path].shape[0]:
            states = np.asarray(h[state_path][index], dtype=int)
            data["myosin_bound_head_count"] = np.count_nonzero(
                np.isin(states, (STATE_WEAK, STATE_STRONG)), axis=1)
            data["myosin_strong_head_count"] = np.count_nonzero(
                states == STATE_STRONG, axis=1)
        data["step"] = index
        if "/state/current_step" in h:
            data["step"] = int(np.asarray(h["/state/current_step"][index]).reshape(-1)[0])
        data["local_links"] = self.local_links(index, data)
        data["heads"] = self.bound_heads(index, data)
        return data

    def local_links(self, frame: int, data: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        h = self.handle
        names = (
            "/catch_bond/local_link_frame_offsets",
            "/catch_bond/local_link_identity",
            "/catch_bond/local_link_coordinates",
        )
        if any(name not in h for name in names) or frame >= h[names[0]].shape[0]:
            return (np.empty((0, 3)), np.empty((0, 3)),
                    np.empty(0, dtype=int), np.empty((0, 2), dtype=int))
        _, offset, count = np.asarray(h[names[0]][frame], dtype=int)
        identity = np.asarray(h[names[1]][offset : offset + count], dtype=int)
        coords = np.asarray(h[names[2]][offset : offset + count], dtype=float)
        if not len(identity):
            return (np.empty((0, 3)), np.empty((0, 3)),
                    np.empty(0, dtype=int), np.empty((0, 2), dtype=int))
        ai, aj = identity[:, 2], identity[:, 3]
        starts = data["actin_center"][ai] + coords[:, :1] * data["actin_direction"][ai]
        raw_ends = data["actin_center"][aj] + coords[:, 1:] * data["actin_direction"][aj]
        ends = starts + minimum_image(raw_ends - starts, self.box)
        return starts, ends, identity[:, 4], identity[:, 2:4]

    def bound_heads(self, frame: int, data: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        h = self.handle
        state_path = "/state/am_head_state_current"
        bound_path = "/state/am_head_bound_actin_current"
        azimuth_path = "/state/am_head_binding_azimuth_current"
        e1_path = "/state/myosin_material_e1_current"
        e2_path = "/state/myosin_material_e2_current"
        required = (state_path, bound_path, azimuth_path, e1_path, e2_path)
        if any(name not in h for name in required) or frame >= h[state_path].shape[0]:
            return (np.empty((0, 3)), np.empty((0, 3)), np.empty(0, dtype=int),
                    np.empty(0, dtype=int), np.empty(0, dtype=int))
        states = np.asarray(h[state_path][frame], dtype=int)
        bound = np.asarray(h[bound_path][frame], dtype=int)
        azimuth = np.asarray(h[azimuth_path][frame], dtype=float)
        e1 = np.asarray(h[e1_path][frame], dtype=float)
        e2 = np.asarray(h[e2_path][frame], dtype=float)
        heads_per_myosin = states.shape[1]
        heads_per_half = heads_per_myosin // 2
        myosin_radius = float(self.meta.get("myosin_radius", 0.025))
        starts, ends, visible_states, visible_myosins, visible_actins = [], [], [], [], []
        for myosin_index, head_index in np.argwhere((bound >= 0) & np.isin(states, (STATE_WEAK, STATE_STRONG))):
            actin_index = int(bound[myosin_index, head_index])
            if not 0 <= actin_index < len(data["actin_center"]):
                continue
            half_index, local_head = divmod(int(head_index), heads_per_half)
            half_sign = -1 if half_index == 0 else 1
            fraction = (local_head + 0.5) / heads_per_half
            bare_half = 0.5 * float(
                self.meta.get("am_bare_zone_fraction", 0.10)) * self.myosin_length
            head_region = 0.5 * self.myosin_length - bare_half
            axial = half_sign * (bare_half + fraction * head_region)
            axis_m = data["myosin_direction"][myosin_index]
            angle = azimuth[myosin_index, head_index]
            anchor = (
                data["myosin_center"][myosin_index] + axial * axis_m
                + myosin_radius * (np.cos(angle) * e1[myosin_index] + np.sin(angle) * e2[myosin_index])
            )
            axis_a = data["actin_direction"][actin_index]
            center_image = anchor + minimum_image(data["actin_center"][actin_index] - anchor, self.box)
            coordinate = np.clip(np.dot(anchor - center_image, axis_a),
                                 -0.5 * self.actin_length, 0.5 * self.actin_length)
            attachment = center_image + coordinate * axis_a
            starts.append(anchor); ends.append(attachment)
            visible_states.append(states[myosin_index, head_index])
            visible_myosins.append(myosin_index); visible_actins.append(actin_index)
        return (np.asarray(starts).reshape((-1, 3)),
                np.asarray(ends).reshape((-1, 3)),
                np.asarray(visible_states, dtype=int),
                np.asarray(visible_myosins, dtype=int),
                np.asarray(visible_actins, dtype=int))


def add_segments(plotter: pv.Plotter, starts: np.ndarray, ends: np.ndarray,
                 color: str, width: float, opacity: float = 1.0,
                 box: np.ndarray | None = None) -> None:
    if box is not None and len(starts):
        starts, ends = periodic_clip_segments(starts, ends, box)
    if len(starts):
        plotter.add_mesh(line_mesh(starts, ends), color=color, line_width=width,
                         opacity=opacity, render_lines_as_tubes=True,
                         lighting=True)


def add_force_arrows(plotter: pv.Plotter, centers: np.ndarray, forces: np.ndarray,
                     color: str, max_arrows: int) -> None:
    centers = np.asarray(centers, dtype=float).reshape((-1, 3))
    forces = np.asarray(forces, dtype=float).reshape((-1, 3))
    magnitude = np.linalg.norm(forces, axis=1)
    valid = np.flatnonzero(magnitude > 1e-12)
    if len(valid) > max_arrows:
        valid = valid[np.argsort(magnitude[valid])[-max_arrows:]]
    if not len(valid):
        return
    scale = 0.25 / max(np.percentile(magnitude[valid], 95), 1e-12)
    plotter.add_arrows(centers[valid], forces[valid], mag=scale, color=color)


def add_dashboard(destination: Path, series: dict[str, np.ndarray],
                  output_index: int) -> None:
    """Append synchronized KMC and order traces to a rendered 3-D frame."""
    image = plt.imread(destination)
    fig = plt.figure(figsize=(16, 9), facecolor="white")
    grid = fig.add_gridspec(3, 2, width_ratios=(3.2, 1.0), hspace=0.36, wspace=0.16)
    visual = fig.add_subplot(grid[:, 0])
    activation = fig.add_subplot(grid[0, 1])
    counts = fig.add_subplot(grid[1, 1])
    order = fig.add_subplot(grid[2, 1])
    visual.imshow(image)
    visual.axis("off")
    stop = output_index + 1
    time_ms = 1e3 * series["time"]
    activation.plot(time_ms[:stop], series["calcium"][:stop], color="#777777",
                    lw=2, label="calcium")
    activation.plot(time_ms[:stop], series["q"][:stop], color="#222222",
                    lw=2, label="mean q")
    activation.set_ylabel("activation")
    activation.set_ylim(-.02, 1.02)
    activation.legend(fontsize=8, ncol=2)
    activation.grid(alpha=.2)
    counts.plot(time_ms[:stop], series["weak"][:stop], color="#8e6bbf",
                lw=2, label="weak AM")
    counts.plot(time_ms[:stop], series["strong"][:stop], color="#5e3c99",
                lw=2, label="strong AM")
    counts.plot(time_ms[:stop], series["aa1"][:stop], color="#f0ad00",
                lw=2, label="status-1 ACTN2")
    counts.plot(time_ms[:stop], series["aa2"][:stop], color="#18a558",
                lw=2, label="status-2 ACTN2")
    counts.set_ylabel("number attached")
    counts.legend(fontsize=8, ncol=2)
    counts.grid(alpha=.2)
    order.plot(time_ms[:stop], series["actin_nematic"][:stop], color="#4477aa",
               lw=2, label="actin nematic")
    order.plot(time_ms[:stop], series["myosin_nematic"][:stop], color="#d9822b",
               lw=2, label="myosin nematic")
    order.set_xlabel("time (ms)")
    order.set_ylabel("nematic order")
    order.set_ylim(0, 1.02)
    order.legend(fontsize=8)
    order.grid(alpha=.2)
    for axis in (activation, counts, order):
        axis.set_xlim(time_ms[0], time_ms[-1] if len(time_ms) > 1 else time_ms[0] + 1)
        axis.axvline(time_ms[output_index], color="#333333", ls=":", lw=1)
    fig.savefig(destination, dpi=100, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def render_frame(traj: Trajectory, index: int, destination: Path,
                 args: argparse.Namespace, output_index: int | None = None) -> None:
    data = traj.frame(index)
    plotter = pv.Plotter(off_screen=True, window_size=(args.width, args.height))
    clean = args.style == "clean"
    plotter.set_background("white" if clean else "#10151c")
    a0, a1 = rods(data["actin_center"], data["actin_direction"], traj.actin_length)
    m0, m1 = rods(data["myosin_center"], data["myosin_direction"], traj.myosin_length)
    if args.crop is None:
        actin_visible = np.ones(len(a0), dtype=bool)
        myosin_visible = np.ones(len(m0), dtype=bool)
        crop_center = np.zeros(3)
        crop_extent = traj.box
    else:
        bounds = np.asarray(args.crop, dtype=float).reshape((3, 2))
        actin_visible = np.all(
            (data["actin_center"] >= bounds[:, 0]) &
            (data["actin_center"] <= bounds[:, 1]), axis=1)
        myosin_visible = np.all(
            (data["myosin_center"] >= bounds[:, 0]) &
            (data["myosin_center"] <= bounds[:, 1]), axis=1)
        crop_center = bounds.mean(axis=1)
        crop_extent = bounds[:, 1] - bounds[:, 0]
    status = np.asarray(data.get("actin_status", np.zeros(len(a0))), dtype=int).reshape(-1)
    weak_count = np.asarray(data.get("weak_count", np.zeros(len(a0))), dtype=float).reshape(-1)
    strong_count = np.asarray(data.get("strong_count", np.zeros(len(a0))), dtype=float).reshape(-1)
    myosin_bound = np.asarray(
        data.get("myosin_bound_head_count", np.zeros(len(m0))), dtype=int).reshape(-1)
    myosin_strong = np.asarray(
        data.get("myosin_strong_head_count", np.zeros(len(m0))), dtype=int).reshape(-1)
    head_start, head_end, head_state, head_myosin, head_actin = data["heads"]
    if args.selection == "engaged":
        actin_visible &= (status > 0) | (weak_count + strong_count > 0)
        myosin_visible &= myosin_bound > 0
    elif args.selection == "loaded":
        actin_visible &= (status == 2) | (strong_count > 0)
        myosin_visible &= myosin_strong > 0
    elif args.selection == "cb2-attached":
        actin_visible &= status == 2
        attached_myosin = np.zeros(len(m0), dtype=bool)
        if len(head_myosin):
            keep = actin_visible[head_actin]
            attached_myosin[np.unique(head_myosin[keep])] = True
        myosin_visible &= attached_myosin
    elif args.selection == "legacy":
        load = np.asarray(data.get("actin_load", np.zeros(len(a0))), dtype=float).reshape(-1)
        actin_visible &= status * np.maximum(load, 0.0) >= args.load_threshold
        legacy_myosin = np.zeros(len(m0), dtype=bool)
        bonds = np.asarray(data.get("myosin_bonds", np.empty((0, 2))), dtype=int).reshape((-1, 2))
        valid = bonds[bonds[:, 0] >= 0] if len(bonds) else bonds
        if len(valid):
            indices = np.unique(valid)
            legacy_myosin[indices[(indices >= 0) & (indices < len(m0))]] = True
        myosin_visible &= legacy_myosin

    if args.actin_color == "status" and "actin_status" in data:
        palette = (
            ((0, "#777777"), (1, "#4c78a8"), (2, "#8f3f71")) if clean else
            ((0, "#4d8cc9"), (1, "#36c9c6"), (2, "#ff4da6"))
        )
        for value, color in palette:
            mask = (status == value) & actin_visible
            add_segments(plotter, a0[mask], a1[mask], color, args.actin_width, .92, traj.box)
    elif args.actin_color in ("polarity", "polarity-blue"):
        plus = (data["actin_direction"][:, 0] >= 0) & actin_visible
        plus_color = "#78b7eb" if args.actin_color == "polarity-blue" else "#ee6677"
        minus_color = "#174f86" if args.actin_color == "polarity-blue" else "#4477aa"
        add_segments(plotter, a0[plus], a1[plus], plus_color, args.actin_width, .94, traj.box)
        minus = (data["actin_direction"][:, 0] < 0) & actin_visible
        add_segments(plotter, a0[minus], a1[minus], minus_color, args.actin_width, .94, traj.box)
    else:
        add_segments(plotter, a0[actin_visible], a1[actin_visible], "#5ba4e6", args.actin_width, .9, traj.box)
    add_segments(plotter, m0[myosin_visible], m1[myosin_visible],
                 "#d9822b" if clean else "#f2c14e", args.myosin_width, .82, traj.box)

    if args.show_polarity:
        barbed = a1[actin_visible]
        barbed = (barbed + .5 * traj.box) % traj.box - .5 * traj.box
        points = pv.PolyData(barbed)
        plotter.add_mesh(points, color="#f8f8f2", point_size=args.polarity_size,
                         render_points_as_spheres=True)

    link_start, link_end, link_status, link_actins = data["local_links"]
    if args.selection == "cb2-attached" and len(link_start):
        keep = (actin_visible[link_actins[:, 0]] &
                actin_visible[link_actins[:, 1]] & (link_status == 2))
        link_start, link_end, link_status, link_actins = (
            link_start[keep], link_end[keep], link_status[keep], link_actins[keep])
    if args.crop is not None and len(link_start):
        bounds = np.asarray(args.crop, dtype=float).reshape((3, 2))
        midpoint = .5 * (link_start + link_end)
        keep = np.all((midpoint >= bounds[:, 0]) & (midpoint <= bounds[:, 1]), axis=1)
        link_start, link_end, link_status = link_start[keep], link_end[keep], link_status[keep]
    if args.show_links:
        add_segments(plotter, link_start[link_status == 1], link_end[link_status == 1],
                     "#39d0c8", args.link_width, .9, traj.box)
        add_segments(plotter, link_start[link_status == 2], link_end[link_status == 2],
                     "#ff3b9d", 1.4 * args.link_width, 1.0, traj.box)

    if args.selection == "cb2-attached" and len(head_start):
        keep = actin_visible[head_actin] & myosin_visible[head_myosin]
        head_start, head_end, head_state, head_myosin, head_actin = (
            head_start[keep], head_end[keep], head_state[keep],
            head_myosin[keep], head_actin[keep])
    if args.crop is not None and len(head_start):
        bounds = np.asarray(args.crop, dtype=float).reshape((3, 2))
        midpoint = .5 * (head_start + head_end)
        keep = np.all((midpoint >= bounds[:, 0]) & (midpoint <= bounds[:, 1]), axis=1)
        head_start, head_end, head_state = head_start[keep], head_end[keep], head_state[keep]
    if len(head_start) > args.max_heads:
        choice = np.linspace(0, len(head_start) - 1, args.max_heads, dtype=int)
        head_start, head_end, head_state = head_start[choice], head_end[choice], head_state[choice]
    if args.show_heads:
        add_segments(plotter, head_start[head_state == STATE_WEAK], head_end[head_state == STATE_WEAK],
                     "#ffd166", args.head_width, .75, traj.box)
        add_segments(plotter, head_start[head_state == STATE_STRONG], head_end[head_state == STATE_STRONG],
                     "#ff5a36", 1.35 * args.head_width, 1.0, traj.box)

    if args.show_forces in ("actin", "both") and "actin_force" in data:
        add_force_arrows(plotter, data["actin_center"][actin_visible],
                         np.asarray(data["actin_force"]).reshape((-1, 3))[actin_visible],
                         "#70d6ff", args.max_arrows)
    if args.show_forces in ("myosin", "both") and "myosin_force" in data:
        add_force_arrows(plotter, data["myosin_center"][myosin_visible],
                         np.asarray(data["myosin_force"]).reshape((-1, 3))[myosin_visible],
                         "#ff9770", args.max_arrows)

    box = pv.Box(bounds=(-traj.box[0]/2, traj.box[0]/2, -traj.box[1]/2,
                         traj.box[1]/2, -traj.box[2]/2, traj.box[2]/2))
    plotter.add_mesh(box, style="wireframe", color="#777777" if clean else "#a0a7b0",
                     opacity=.22 if clean else .32, line_width=1)
    dt = float(traj.meta.get("dt", 1.0))
    text = (
        f"frame {index}/{traj.nframes-1}   step {data['step']}   t={data['step']*dt:.5g} s\n"
        f"S_actin={nematic(data['actin_direction']):.3f}   S_myosin={nematic(data['myosin_direction']):.3f}   "
        f"AA1={np.count_nonzero(link_status==1)}   AA2={np.count_nonzero(link_status==2)}   "
        f"AMweak={np.count_nonzero(head_state==STATE_WEAK)}   AMstrong={np.count_nonzero(head_state==STATE_STRONG)}"
    )
    plotter.add_text(text, position="upper_left", font_size=args.font_size,
                     color="black" if clean else "white", shadow=not clean)
    if args.camera == "axial":
        plotter.camera_position = [crop_center + (1.2 * max(crop_extent), 0, 0), crop_center, (0, 0, 1)]
    elif args.camera == "transverse":
        plotter.camera_position = [crop_center + (0, 0, 2.2 * max(crop_extent)), crop_center, (0, 1, 0)]
    else:
        plotter.camera_position = [crop_center + (0.8*crop_extent[0], 1.8*crop_extent[1], 1.8*crop_extent[2]),
                                   crop_center, (0, 0, 1)]
    plotter.camera.parallel_projection = args.parallel_projection
    plotter.show_axes()
    plotter.screenshot(str(destination))
    plotter.close()
    if output_index is not None and hasattr(args, "dashboard_series"):
        add_dashboard(destination, args.dashboard_series, output_index)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trajectory", type=Path)
    parser.add_argument("output", type=Path, help=".png for one frame or .mp4 for an animation")
    parser.add_argument("--box", type=float, nargs=3, default=(12, 2, 2), metavar=("LX", "LY", "LZ"))
    parser.add_argument("--actin-length", type=float,
                        help="display length; default reads trajectory metadata")
    parser.add_argument("--myosin-length", type=float,
                        help="display length; default reads trajectory metadata")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--stop", type=int)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--frame", type=int, default=-1, help="frame for PNG output; -1 means last")
    parser.add_argument("--fps", type=int, default=15)
    parser.add_argument("--style", choices=("clean", "debug"), default="clean",
                        help="clean suppresses mechanics clutter; debug exposes all filaments")
    parser.add_argument("--selection", choices=("all", "engaged", "loaded", "cb2-attached", "legacy"),
                        help="filament subset; defaults to engaged for clean and all for debug")
    parser.add_argument("--load-threshold", type=float, default=0.01,
                        help="legacy cb_status*f_load visibility threshold")
    parser.add_argument("--camera", choices=("isometric", "axial", "transverse"), default="isometric")
    parser.add_argument("--actin-color", choices=("status", "polarity", "polarity-blue", "plain"), default="status")
    parser.add_argument("--show-polarity", action="store_true")
    parser.add_argument("--show-links", action="store_true",
                        help="draw individual local ACTN2 connectors")
    parser.add_argument("--show-heads", action="store_true",
                        help="draw individual weak/strong AM head connectors")
    parser.add_argument("--show-forces", choices=("none", "actin", "myosin", "both"), default="none")
    parser.add_argument("--max-arrows", type=int, default=150)
    parser.add_argument("--max-heads", type=int, default=5000)
    parser.add_argument("--actin-width", type=float)
    parser.add_argument("--myosin-width", type=float)
    parser.add_argument("--link-width", type=float)
    parser.add_argument("--head-width", type=float)
    parser.add_argument("--polarity-size", type=float, default=5)
    parser.add_argument("--font-size", type=int, default=10)
    parser.add_argument("--width", type=int, default=1600)
    parser.add_argument("--height", type=int, default=900)
    parser.add_argument("--parallel-projection", action="store_true")
    parser.add_argument("--no-dashboard", action="store_true",
                        help="render only the 3-D scene without global diagnostic traces")
    parser.add_argument("--crop", type=float, nargs=6,
                        metavar=("XMIN", "XMAX", "YMIN", "YMAX", "ZMIN", "ZMAX"),
                        help="render only objects centered in this debugging region")
    args = parser.parse_args()
    if args.selection is None:
        args.selection = "loaded" if args.style == "clean" else "all"
    defaults = {
        "clean": (2.0, 4.0, 2.0, 1.5),
        "debug": (4.0, 8.0, 3.0, 2.0),
    }[args.style]
    for name, value in zip(
            ("actin_width", "myosin_width", "link_width", "head_width"), defaults):
        if getattr(args, name) is None:
            setattr(args, name, value)
    if args.stride < 1:
        parser.error("--stride must be positive")
    if (not os.environ.get("DISPLAY") and
            not Path("/dev/nvidiactl").exists() and
            ctypes.util.find_library("OSMesa") is None):
        parser.error(
            "no DISPLAY, GPU EGL device, or OSMesa backend is available; "
            "run on a GPU node (for example the gpu_testjob allocation)"
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    traj = Trajectory(args.trajectory, np.asarray(args.box, dtype=float),
                      args.actin_length, args.myosin_length)
    try:
        if args.output.suffix.lower() == ".png":
            frame = args.frame if args.frame >= 0 else traj.nframes - 1
            render_frame(traj, min(frame, traj.nframes - 1), args.output, args)
            print(f"wrote {args.output}")
            return
        if args.output.suffix.lower() != ".mp4":
            parser.error("output must end in .png or .mp4")
        stop = min(args.stop if args.stop is not None else traj.nframes, traj.nframes)
        frames = list(range(max(0, args.start), stop, args.stride))
        if not frames:
            parser.error("selected frame range is empty")
        if not args.no_dashboard:
            args.dashboard_series = traj.diagnostic_series(frames)
        with tempfile.TemporaryDirectory(prefix="sarcomere_frames_", dir=args.output.parent) as raw:
            directory = Path(raw)
            for output_index, frame in enumerate(frames):
                render_frame(traj, frame, directory / f"frame_{output_index:06d}.png",
                             args, output_index)
                print(f"rendered {output_index+1}/{len(frames)} (trajectory frame {frame})", flush=True)
            subprocess.run([
                "ffmpeg", "-y", "-framerate", str(args.fps),
                "-i", str(directory / "frame_%06d.png"),
                "-c:v", "libx264", "-preset", "slow", "-crf", "18",
                "-pix_fmt", "yuv420p", "-movflags", "+faststart", str(args.output),
            ], check=True)
        print(f"wrote {args.output}")
    finally:
        traj.close()


if __name__ == "__main__":
    main()
