#!/usr/bin/env python3
"""Animation wrapper faithful to Sherry's clean ``visualize_traj3d`` view.

The visual choices are intentionally conservative: cylinders, the original
orange/light-blue/dark-blue palette, white background, fixed camera, and no
diagnostic overlays by default. Only HDF5 streaming and MP4 assembly are new.
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
import numpy as np
import pyvista as pv


def boundary_offsets(center: np.ndarray, direction: np.ndarray,
                     length: float, box: np.ndarray) -> list[np.ndarray]:
    half_extent = .5 * length * np.abs(direction)
    choices = []
    for coordinate, extent, side in zip(center, half_extent, box):
        offsets = [0.0]
        if coordinate - extent < -.5 * side:
            offsets.append(side)
        if coordinate + extent > .5 * side:
            offsets.append(-side)
        choices.append(offsets)
    return [np.asarray(value) for value in product(*choices)]


def add_filaments(plotter: pv.Plotter, centers: np.ndarray,
                  directions: np.ndarray, length: float, radius: float,
                  box: np.ndarray, color: str, opacity: float,
                  resolution: int) -> None:
    for center, direction in zip(centers, directions):
        for offset in boundary_offsets(center, direction, length, box):
            cylinder = pv.Cylinder(
                center=center + offset, direction=direction, height=length,
                radius=radius, resolution=resolution,
            )
            plotter.add_mesh(cylinder, color=color, opacity=opacity)


def add_box(plotter: pv.Plotter, box: np.ndarray) -> None:
    lx, ly, lz = box
    corners = np.array([
        [-lx/2, -ly/2, -lz/2], [lx/2, -ly/2, -lz/2],
        [lx/2, ly/2, -lz/2], [-lx/2, ly/2, -lz/2],
        [-lx/2, -ly/2, lz/2], [lx/2, -ly/2, lz/2],
        [lx/2, ly/2, lz/2], [-lx/2, ly/2, lz/2],
    ])
    faces = np.hstack([
        [4, 0, 1, 2, 3], [4, 4, 5, 6, 7], [4, 0, 1, 5, 4],
        [4, 2, 3, 7, 6], [4, 1, 2, 6, 5], [4, 3, 0, 4, 7],
    ])
    plotter.add_mesh(pv.PolyData(corners, faces), style="wireframe",
                     color="black", line_width=1, opacity=.4)


class Trajectory:
    def __init__(self, path: Path):
        self.handle = h5py.File(path, "r", locking=False)
        self.meta = self.handle["/meta"].attrs if "/meta" in self.handle else {}
        required = ("/actin/center", "/actin/direction",
                    "/myosin/center", "/myosin/direction")
        missing = [name for name in required if name not in self.handle]
        if missing:
            raise ValueError(f"missing datasets: {missing}")
        self.nframes = min(self.handle[name].shape[0] for name in required)

    def close(self) -> None:
        self.handle.close()

    def frame(self, index: int) -> dict[str, np.ndarray]:
        h = self.handle
        data = {
            "actin_center": np.asarray(h["/actin/center"][index], dtype=float),
            "actin_direction": np.asarray(h["/actin/direction"][index], dtype=float),
            "myosin_center": np.asarray(h["/myosin/center"][index], dtype=float),
            "myosin_direction": np.asarray(h["/myosin/direction"][index], dtype=float),
            "cb_status": np.asarray(h["/actin/cb_status"][index], dtype=int).reshape(-1),
        }
        if "/actin_myo/bonds" in h and index < h["/actin_myo/bonds"].shape[0]:
            data["am_bonds"] = np.asarray(h["/actin_myo/bonds"][index], dtype=int)
        if "/myosin/bonds" in h and index < h["/myosin/bonds"].shape[0]:
            data["myosin_bonds"] = np.asarray(h["/myosin/bonds"][index], dtype=int)
        return data


def selected_myosins(data: dict[str, np.ndarray], mode: str) -> tuple[np.ndarray, np.ndarray]:
    n = len(data["myosin_center"])
    am = data.get("am_bonds", np.empty((0, 2), dtype=int)).reshape((-1, 2))
    am = am[am[:, 0] >= 0] if len(am) else am
    highlighted = np.empty(0, dtype=int)
    if len(am):
        cb_actins = np.flatnonzero(data["cb_status"] == 2)
        highlighted = np.unique(am[np.isin(am[:, 0], cb_actins), 1])
        highlighted = highlighted[(highlighted >= 0) & (highlighted < n)]
    if mode == "all":
        displayed = np.arange(n)
    elif mode == "cb_attached":
        displayed = highlighted
    else:
        bonds = data.get("myosin_bonds", np.empty((0, 2), dtype=int)).reshape((-1, 2))
        bonds = bonds[bonds[:, 0] >= 0] if len(bonds) else bonds
        displayed = np.unique(bonds) if len(bonds) else np.empty(0, dtype=int)
        displayed = displayed[(displayed >= 0) & (displayed < n)]
    return displayed, highlighted


def render(traj: Trajectory, frame: int, output: Path, args: argparse.Namespace) -> None:
    data = traj.frame(frame)
    actin_display_length = args.actin_length * args.display_length_scale
    myosin_display_length = args.myosin_length * args.display_length_scale
    plotter = pv.Plotter(off_screen=True, window_size=(args.width, args.height))
    plotter.set_background("white")
    if args.actin_display == "cb":
        actin_mask = data["cb_status"] == 2
    else:
        actin_mask = np.ones(len(data["actin_center"]), dtype=bool)
    add_filaments(
        plotter, data["actin_center"][actin_mask],
        data["actin_direction"][actin_mask], actin_display_length,
        args.actin_radius, args.box, args.actin_color, 1.0, args.resolution,
    )

    displayed, highlighted = selected_myosins(data, args.myosin_display)
    highlight_mask = np.isin(displayed, highlighted)
    add_filaments(
        plotter, data["myosin_center"][displayed[~highlight_mask]],
        data["myosin_direction"][displayed[~highlight_mask]],
        myosin_display_length, args.myosin_radius, args.box,
        args.nonbonded_myosin_color, args.myosin_opacity, args.resolution,
    )
    add_filaments(
        plotter, data["myosin_center"][displayed[highlight_mask]],
        data["myosin_direction"][displayed[highlight_mask]],
        myosin_display_length, args.myosin_radius, args.box,
        args.bonded_myosin_color, args.myosin_opacity, args.resolution,
    )
    if args.show_box:
        add_box(plotter, args.box)
    plotter.set_focus((0, 0, 0))
    # Sherry's reference view uses (10, 10, 10) for an Lx=12 box. Preserve
    # that framing for reduced boxes instead of leaving them tiny on canvas.
    camera_component = (args.camera_distance if args.camera_distance is not None
                        else 10.0 * args.box[0] / 12.0)
    plotter.camera_position = [
        (camera_component, camera_component, camera_component),
        (0, 0, 0), (0, 0, 1),
    ]
    plotter.screenshot(str(output))
    plotter.close()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trajectory", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--box", type=float, nargs=3, default=(12, 2, 2),
                        metavar=("LX", "LY", "LZ"))
    parser.add_argument("--actin-length", type=float, default=1.0)
    parser.add_argument("--myosin-length", type=float, default=1.5)
    parser.add_argument(
        "--display-length-scale", type=float, default=1.0,
        help=("display-only multiplier for both filament lengths; use 0.25 "
              "for an Lx=3 reduced test to match the proportions of the "
              "original Lx=12 visualization"),
    )
    parser.add_argument("--actin-radius", type=float, default=.01)
    parser.add_argument("--myosin-radius", type=float, default=.025)
    parser.add_argument("--actin-display", choices=("cb", "all"), default="cb")
    parser.add_argument("--myosin-display", choices=("all", "bonded", "cb_attached"), default="all")
    parser.add_argument("--actin-color", default="#F95D1A")
    parser.add_argument("--bonded-myosin-color", default="#0077be")
    parser.add_argument("--nonbonded-myosin-color", default="#89cff0")
    parser.add_argument("--myosin-opacity", type=float, default=1.0)
    parser.add_argument(
        "--camera-distance", type=float,
        help="camera coordinate along x, y and z (default: 10*Lx/12)",
    )
    parser.add_argument("--show-box", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--resolution", type=int, default=24)
    parser.add_argument("--frame", type=int, default=-1)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--stop", type=int)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--fps", type=int, default=15)
    parser.add_argument("--width", type=int, default=1600)
    parser.add_argument("--height", type=int, default=1100)
    args = parser.parse_args()
    args.box = np.asarray(args.box, dtype=float)
    if args.display_length_scale <= 0:
        parser.error("--display-length-scale must be positive")
    if (not os.environ.get("DISPLAY") and not Path("/dev/nvidiactl").exists()
            and ctypes.util.find_library("OSMesa") is None):
        parser.error("render on a GPU/EGL node, X display, or OSMesa node")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    traj = Trajectory(args.trajectory)
    try:
        if args.output.suffix.lower() == ".png":
            frame = args.frame if args.frame >= 0 else traj.nframes - 1
            render(traj, min(frame, traj.nframes - 1), args.output, args)
            print(f"wrote {args.output}")
            return
        if args.output.suffix.lower() != ".mp4":
            parser.error("output must be .png or .mp4")
        stop = min(args.stop if args.stop is not None else traj.nframes, traj.nframes)
        frames = list(range(max(0, args.start), stop, max(1, args.stride)))
        with tempfile.TemporaryDirectory(prefix="visualize_traj3d_", dir=args.output.parent) as raw:
            directory = Path(raw)
            for number, frame in enumerate(frames):
                render(traj, frame, directory / f"frame_{number:06d}.png", args)
                print(f"rendered {number+1}/{len(frames)}", flush=True)
            subprocess.run([
                "ffmpeg", "-y", "-framerate", str(args.fps),
                "-i", str(directory / "frame_%06d.png"), "-c:v", "libx264",
                "-crf", "18", "-pix_fmt", "yuv420p", "-movflags", "+faststart",
                str(args.output),
            ], check=True)
        print(f"wrote {args.output}")
    finally:
        traj.close()


if __name__ == "__main__":
    main()
