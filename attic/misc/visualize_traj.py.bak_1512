import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import math
import joblib
from joblib import Parallel, delayed
import os
import shutil
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import sys
import argparse
from mpl_toolkits.mplot3d import Axes3D  # required for 3D plotting
from itertools import product
import pyvista as pv

# Prefer native off-screen rendering; only attempt Xvfb if it is available.
os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")


def configure_headless_rendering():
    """Configure PyVista for headless rendering without requiring Xvfb."""
    pv.OFF_SCREEN = True

    # Try Xvfb only when installed. If unavailable or fails, continue using
    # native off-screen rendering (EGL/OSMesa) so PNG/HTML export still works.
    if shutil.which("Xvfb") is None:
        print("Xvfb not found; using PyVista native off-screen rendering.")
        return

    try:
        pv.start_xvfb()
        print("Started Xvfb for PyVista off-screen rendering.")
    except Exception as exc:
        print(f"Warning: could not start Xvfb ({exc}); continuing with native off-screen rendering.")


configure_headless_rendering()


def should_use_matplotlib_backend(render_backend="auto"):
    """Decide whether to use matplotlib fallback rendering."""
    if render_backend == "matplotlib":
        return True
    if render_backend == "pyvista":
        return False

    # auto mode: if we're headless and Xvfb is unavailable, prefer matplotlib
    # to avoid VTK/EGL/OSMesa crashes on cluster nodes.
    headless = not bool(os.environ.get("DISPLAY"))
    has_xvfb = shutil.which("Xvfb") is not None
    return headless and (not has_xvfb)


import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from itertools import product


# (Include your existing functions: get_offsets_for_filament, plot_filaments_3d, plot_myosin_velocity, etc.)

def plot_bonded_myosins_for_frame(data, frame, myosin_length, Lx, Ly, Lz):
    bonds = data["/myosin/bonds"][frame]
    valid_pairs = bonds[bonds[:, 0] >= 0].astype(int)
    if valid_pairs.size == 0:
        print(f"No bonded myosins in frame {frame}")
        return

    print(f"\nPlotted Frame {frame}: {len(valid_pairs)} bonded pairs.")

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection='3d')

    cmap = plt.cm.get_cmap('tab20', len(valid_pairs))
    norm = mcolors.Normalize(vmin=0, vmax=len(valid_pairs) - 1)

    for pid, (iA, iB) in enumerate(valid_pairs):
        color = cmap(norm(pid))
        for idx in (iA, iB):
            ctr = data["/myosin/center"][frame][idx:idx+1]
            dirv = data["/myosin/direction"][frame][idx:idx+1]
            plot_filaments_3d(ctr, dirv, myosin_length, Lx, Ly, Lz, ax, color=color)

    ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.set_zlabel("Z")
    ax.set_title(f"Bonded Myosin Pairs (Frame {frame})")
    plt.tight_layout()
    plt.savefig(f"bonded_myosins_frame_{frame}.png", dpi=300)
    plt.close(fig)

def print_and_plot_last_frame(filename, myosin_length, Lx, Ly, Lz):
    with h5py.File(filename, 'r') as traj:
        data = hdf5_to_dict(traj)
    last_frame = data["/actin/center"].shape[0] - 1
    print_bonded_myosin_info(data, start_frame=last_frame, end_frame=last_frame+1)
    #plot_bonded_myosins_for_frame(data, last_frame, myosin_length, Lx, Ly, Lz)



def get_offsets_for_filament(coord, d, L):
    """
    Given a coordinate (e.g. x), half-length d along that axis, and box length L,
    compute the set of translations needed so that the filament endpoints fall within the box.
    The endpoints are defined as: p1 = coord - d and p2 = coord + d.
      - If p1 < -L/2, then add an offset of +L.
      - If p2 > L/2, then add an offset of -L.
    Always include 0.
    """
    offsets = {0}
    p1 = coord - d
    p2 = coord + d
    if p1 < -L/2:
        offsets.add(L)
    if p2 > L/2:
        offsets.add(-L)
    return offsets


def plot_filaments_3d(plotter, center, direction, radius, l, Lx, Ly, Lz, color='k', color_spectrum=None, opacity=1.0):
    """
    Plot actin filaments in 3D as cylinders using PyVista.
    - center: (N,3) array of filament centers.
    - direction: (N,3) array of unit direction vectors.
    - radius: cylinder radius (scalar or array).
    - l: filament length (scalar or array).
    - plotter: a pyvista.Plotter instance.
    """
    if np.isscalar(l):
        l = np.ones(center.shape[0]) * l
    if np.isscalar(radius):
        radius = np.ones(center.shape[0]) * radius

    color_values = None
    norm = None
    if color_spectrum is not None and color_spectrum.size > 0:
        color_values = np.asarray(color_spectrum, dtype=float).copy()
        color_values = np.clip(color_values, 0, None)
        max_val = np.max(color_values)
        if max_val > 0:
            color_values /= max_val
        color_values = np.sqrt(color_values)
        norm = plt.Normalize(0, 1)

    for i in range(center.shape[0]):
        if l[i] < 0.01:
            continue

        this_color = color
        if color_values is not None:
            this_color = plt.cm.Blues(norm(color_values[i]))

        # PyVista expects cylinder defined by center, direction, height, and radius
        cyl = pv.Cylinder(center=center[i],
                          direction=direction[i],
                          height=l[i],
                          radius=radius[i],
                          resolution=24)
        
        # Add main cylinder
        plotter.add_mesh(cyl, color=this_color, opacity=opacity)

        # Periodic boundary duplicates
        dx = 0.5 * l[i] * np.array(direction[i])
        ox_set = get_offsets_for_filament(center[i][0], abs(dx[0]), Lx)
        oy_set = get_offsets_for_filament(center[i][1], abs(dx[1]), Ly)
        oz_set = get_offsets_for_filament(center[i][2], abs(dx[2]), Lz)
        
        for ox, oy, oz in product(ox_set, oy_set, oz_set):
            if ox == oy == oz == 0:
                continue
            duplicate_center = center[i] + np.array([ox, oy, oz])
            dup_cyl = pv.Cylinder(center=duplicate_center,
                                  direction=direction[i],
                                  height=l[i],
                                  radius=radius[i],
                                  resolution=24)
            plotter.add_mesh(dup_cyl, color=this_color,opacity=opacity)


def add_vector_arrows(plotter, centers, vectors, color, scale=1.0,
                      shaft_radius=0.01, tip_radius=0.02,
                      thickness_scale=1.0, min_magnitude=1e-6,
                      min_display_length=1e-6, normalize=False,
                      arrow_opacity=1.0):
    """
    Draw vector arrows originating from filament centers (forces, velocities, etc.).
    Arrow length reflects the vector magnitude; optionally normalized per-frame.
    """
    if centers.size == 0 or vectors is None or vectors.size == 0:
        return

    magnitudes = np.linalg.norm(vectors, axis=1)
    if normalize:
        max_mag = np.max(magnitudes)
        norm_factor = max(max_mag, min_magnitude)
    else:
        norm_factor = None

    for center, raw_vec, magnitude in zip(centers, vectors, magnitudes):
        if magnitude < min_magnitude:
            continue
        direction_unit = raw_vec / magnitude
        if normalize:
            arrow_length = max((magnitude / norm_factor) * scale, min_display_length)
        else:
            arrow_length = max(magnitude * scale, min_display_length)
        arrow = pv.Arrow(
            start=center,
            direction=direction_unit * arrow_length,
            scale=1.0,
            tip_length=0.2,
            tip_radius=tip_radius * thickness_scale,
            shaft_radius=shaft_radius * thickness_scale,
        )
        plotter.add_mesh(
            arrow,
            color=color,
            opacity=arrow_opacity,
            smooth_shading=True,
            style='surface',
            show_edges=False,
            lighting=True
        )


def plot_system(frame, data, myosin_length, actin_length, Lx, Ly, Lz,
                myosin_radius, myosin_display="all", actin_display="cb",
                show_actin_force=False, show_myosin_force=False,
                show_actin_velocity=False, show_myosin_velocity=False,
                actin_force_scale=1.0, myosin_force_scale=1.0,
                actin_velocity_scale=1.0, myosin_velocity_scale=1.0,
                actin_force_available=False, myosin_force_available=False,
                actin_velocity_available=False, myosin_velocity_available=False,
                actin_force_thickness=3.0, myosin_force_thickness=3.0,
                actin_velocity_thickness=1.5, myosin_velocity_thickness=1.5,
                myosin_opacity=1.0, render_backend=None, **kwargs):
    """Render the system for a single frame."""
    plotter = pv.Plotter(off_screen=True)

    # ------------------------------------------------------------------
    # Actin filaments
    # ------------------------------------------------------------------
    actin_center = data["/actin/center"][frame]
    actin_direction = data["/actin/direction"][frame]
    f_load = data["/actin/f_load"][frame].flatten()
    cb_status = data["/actin/cb_status"][frame].flatten()
    actin_force = data["/actin/force"][frame] if actin_force_available else None
    actin_velocity = data["/actin/velocity"][frame] if actin_velocity_available else None

    load_metric = np.maximum(f_load, 0) * (cb_status == 2)
    if actin_display == "cb":
        mask = cb_status == 2
        print(f"Frame {frame}: {np.sum(mask)} actin filaments in catch-bond state")
    else:
        mask = np.ones_like(cb_status, dtype=bool)
    actin_center = actin_center[mask]
    actin_direction = actin_direction[mask]
    actin_color_values = load_metric[mask]

    plot_filaments_3d(
        center=actin_center,
        direction=actin_direction,
        radius=0.01,
        l=actin_length,
        Lx=Lx, Ly=Ly, Lz=Lz,
        plotter=plotter,
        color='blue',
        color_spectrum=actin_color_values
    )
    if show_actin_force and actin_force_available:
        add_vector_arrows(
            plotter,
            centers=actin_center,
            vectors=actin_force[mask],
            color='#1f77b4',
            scale=actin_force_scale,
            shaft_radius=0.01 * myosin_radius,
            tip_radius=0.02 * myosin_radius,
            thickness_scale=actin_force_thickness,
            normalize=True,
            arrow_opacity=0.6,
        )
    if show_actin_velocity and actin_velocity_available:
        add_vector_arrows(
            plotter,
            centers=actin_center,
            vectors=actin_velocity[mask],
            color='#2ca02c',
            scale=actin_velocity_scale,
            shaft_radius=0.015 * myosin_radius,
            tip_radius=0.03 * myosin_radius,
            thickness_scale=actin_velocity_thickness,
            normalize=True,
            arrow_opacity=0.6,
        )

    # ------------------------------------------------------------------
    # Myosin filaments
    # ------------------------------------------------------------------
    myosin_centers_frame = data["/myosin/center"][frame]
    myosin_dirs_frame = data["/myosin/direction"][frame]
    n_myosins = myosin_centers_frame.shape[0]
    actin_myo_bonds_ds = data.get("/actin_myo/bonds")
    actin_myo_bonds_frame = None
    actin_myo_bonds_available = False
    if actin_myo_bonds_ds is None:
        print("Warning: /actin_myo/bonds dataset not found; skipping bond-based highlights.")
    else:
        shape = getattr(actin_myo_bonds_ds, "shape", None)
        total_am_frames = shape[0] if shape and len(shape) > 0 else 0
        if frame < total_am_frames:
            actin_myo_bonds_frame = actin_myo_bonds_ds[frame]
            actin_myo_bonds_available = True
        else:
            print(
                f"Warning: /actin_myo/bonds has {total_am_frames} frame(s), "
                f"but frame {frame} was requested; skipping bond-based highlights."
            )
    highlight_indices = np.array([], dtype=int)
    if actin_myo_bonds_available:
        valid_pairs = actin_myo_bonds_frame[actin_myo_bonds_frame[:, 0] >= 0]
        strong_actins = np.where(cb_status == 2)[0]
        if valid_pairs.size > 0 and strong_actins.size > 0:
            valid_pairs = valid_pairs.astype(int)
            mask_strong = np.isin(valid_pairs[:, 0], strong_actins)
            if np.any(mask_strong):
                highlight_indices = np.unique(valid_pairs[mask_strong, 1])

    if myosin_display == "bonded":
        myo_bonds = data["/myosin/bonds"][frame]
        valid_pairs = myo_bonds[myo_bonds[:, 0] >= 0].astype(int)
        bonded_indices = np.unique(valid_pairs.flatten())
        display_indices = bonded_indices if bonded_indices.size > 0 else np.empty(0, dtype=int)
    elif myosin_display == "cb_attached":
        if not actin_myo_bonds_available:
            print("Warning: /actin_myo/bonds dataset not found; displaying all myosins.")
            display_indices = np.arange(n_myosins)
        else:
            valid_pairs = actin_myo_bonds_frame[actin_myo_bonds_frame[:, 0] >= 0]
            cb_indices = np.where(cb_status > 1)[0]
            if valid_pairs.size == 0 or cb_indices.size == 0:
                display_indices = np.empty(0, dtype=int)
            else:
                valid_pairs = valid_pairs.astype(int)
                mask_cb_pairs = np.isin(valid_pairs[:, 0], cb_indices)
                cb_pairs = valid_pairs[mask_cb_pairs]
                if cb_pairs.size == 0:
                    display_indices = np.empty(0, dtype=int)
                else:
                    display_indices = np.unique(cb_pairs[:, 1])
    else:  # plot all myosins
        display_indices = np.arange(n_myosins)

    if display_indices.size > 0:
        myosin_center = myosin_centers_frame[display_indices]
        myosin_direction = myosin_dirs_frame[display_indices]
        myosin_force = data["/myosin/force"][frame][display_indices] if myosin_force_available else None
        myosin_velocity = data["/myosin/velocity"][frame][display_indices] if myosin_velocity_available else None
        highlight_mask = np.isin(display_indices, highlight_indices)
        base_center = myosin_center[~highlight_mask]
        base_direction = myosin_direction[~highlight_mask]
        highlight_center = myosin_center[highlight_mask]
        highlight_direction = myosin_direction[highlight_mask]

        if base_center.size > 0:
            plot_filaments_3d(
                center=base_center,
                direction=base_direction,
                radius=myosin_radius,
                l=myosin_length,
                Lx=Lx, Ly=Ly, Lz=Lz,
                plotter=plotter,
                color='lemon_chiffon',
                opacity=myosin_opacity
            )
        if highlight_center.size > 0:
            plot_filaments_3d(
                center=highlight_center,
                direction=highlight_direction,
                radius=myosin_radius,
                l=myosin_length,
                Lx=Lx, Ly=Ly, Lz=Lz,
                plotter=plotter,
                color='#f5a45b',  # light orange
                opacity=myosin_opacity
            )
        if show_myosin_force and myosin_force_available:
            add_vector_arrows(
                plotter,
                centers=myosin_center,
                vectors=myosin_force,
                color='#d62728',
                scale=myosin_force_scale,
                shaft_radius=0.015 * myosin_radius,
                tip_radius=0.03 * myosin_radius,
                thickness_scale=myosin_force_thickness,
                normalize=True,
                arrow_opacity=0.6,
            )
        if show_myosin_velocity and myosin_velocity_available:
            add_vector_arrows(
                plotter,
                centers=myosin_center,
                vectors=myosin_velocity,
                color='#ff7f0e',
                scale=myosin_velocity_scale,
                shaft_radius=0.02 * myosin_radius,
                tip_radius=0.04 * myosin_radius,
                thickness_scale=myosin_velocity_thickness,
                normalize=True,
                arrow_opacity=0.6,
            )

    plotter.set_background("white")
    plotter.set_focus((0, 0, 0))
    # Example: zoom closer to the myosin bundle
    plotter.camera_position = [
        (10, 10, 10),    # camera location
        (0, 0, 0),    # focal point
        (0, 0, 1)     # view-up direction
    ]
    #plotter.camera.zoom(1.5)   # zoom in by a factor
    # Draw the simulation box as a wireframe cube for reference.
    box_corners = np.array([
        [-0.5 * Lx, -0.5 * Ly, -0.5 * Lz],
        [ 0.5 * Lx, -0.5 * Ly, -0.5 * Lz],
        [ 0.5 * Lx,  0.5 * Ly, -0.5 * Lz],
        [-0.5 * Lx,  0.5 * Ly, -0.5 * Lz],
        [-0.5 * Lx, -0.5 * Ly,  0.5 * Lz],
        [ 0.5 * Lx, -0.5 * Ly,  0.5 * Lz],
        [ 0.5 * Lx,  0.5 * Ly,  0.5 * Lz],
        [-0.5 * Lx,  0.5 * Ly,  0.5 * Lz],
    ])
    faces = np.hstack([
        [4, 0, 1, 2, 3],
        [4, 4, 5, 6, 7],
        [4, 0, 1, 5, 4],
        [4, 2, 3, 7, 6],
        [4, 1, 2, 6, 5],
        [4, 3, 0, 4, 7],
    ])
    box_mesh = pv.PolyData(box_corners, faces)
    plotter.add_mesh(box_mesh, style='wireframe', color='black', line_width=1.0, opacity=0.4)
    return plotter


def plot_system_matplotlib(frame, data, myosin_length, actin_length, Lx, Ly, Lz,
                           myosin_radius, myosin_display="all", actin_display="cb",
                           myosin_opacity=1.0, **kwargs):
    """Fallback renderer using matplotlib only (headless-safe)."""
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    def draw_filaments(centers, directions, length, color, linewidth=1.0, alpha=1.0):
        if centers.size == 0:
            return
        lengths = np.ones(centers.shape[0]) * length if np.isscalar(length) else np.asarray(length)
        for c, d, l in zip(centers, directions, lengths):
            if l < 0.01:
                continue
            p1 = c - 0.5 * l * d
            p2 = c + 0.5 * l * d
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]],
                    color=color, linewidth=linewidth, alpha=alpha)

    # Actin
    actin_center = data["/actin/center"][frame]
    actin_direction = data["/actin/direction"][frame]
    cb_status = data["/actin/cb_status"][frame].flatten()
    if actin_display == "cb":
        mask = cb_status == 2
        print(f"Frame {frame}: {np.sum(mask)} actin filaments in catch-bond state")
    else:
        mask = np.ones_like(cb_status, dtype=bool)
    actin_center = actin_center[mask]
    actin_direction = actin_direction[mask]
    draw_filaments(actin_center, actin_direction, actin_length, color='tab:blue', linewidth=1.2, alpha=0.9)

    # Myosin selection logic
    myosin_centers_frame = data["/myosin/center"][frame]
    myosin_dirs_frame = data["/myosin/direction"][frame]
    n_myosins = myosin_centers_frame.shape[0]

    actin_myo_bonds_ds = data.get("/actin_myo/bonds")
    actin_myo_bonds_frame = None
    actin_myo_bonds_available = False
    if actin_myo_bonds_ds is not None:
        shape = getattr(actin_myo_bonds_ds, "shape", None)
        total_am_frames = shape[0] if shape and len(shape) > 0 else 0
        if frame < total_am_frames:
            actin_myo_bonds_frame = actin_myo_bonds_ds[frame]
            actin_myo_bonds_available = True

    if myosin_display == "bonded":
        myo_bonds = data["/myosin/bonds"][frame]
        valid_pairs = myo_bonds[myo_bonds[:, 0] >= 0].astype(int)
        bonded_indices = np.unique(valid_pairs.flatten()) if valid_pairs.size > 0 else np.empty(0, dtype=int)
        display_indices = bonded_indices
    elif myosin_display == "cb_attached":
        if not actin_myo_bonds_available:
            display_indices = np.arange(n_myosins)
        else:
            valid_pairs = actin_myo_bonds_frame[actin_myo_bonds_frame[:, 0] >= 0]
            cb_indices = np.where(cb_status > 1)[0]
            if valid_pairs.size == 0 or cb_indices.size == 0:
                display_indices = np.empty(0, dtype=int)
            else:
                valid_pairs = valid_pairs.astype(int)
                mask_cb_pairs = np.isin(valid_pairs[:, 0], cb_indices)
                cb_pairs = valid_pairs[mask_cb_pairs]
                display_indices = np.unique(cb_pairs[:, 1]) if cb_pairs.size > 0 else np.empty(0, dtype=int)
    else:
        display_indices = np.arange(n_myosins)

    if display_indices.size > 0:
        myosin_center = myosin_centers_frame[display_indices]
        myosin_direction = myosin_dirs_frame[display_indices]
        draw_filaments(myosin_center, myosin_direction, myosin_length,
                       color='#f5a45b', linewidth=max(1.0, myosin_radius * 20), alpha=myosin_opacity)

    # Box and axes
    ax.set_xlim(-0.5 * Lx, 0.5 * Lx)
    ax.set_ylim(-0.5 * Ly, 0.5 * Ly)
    ax.set_zlim(-0.5 * Lz, 0.5 * Lz)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'Frame {frame}')
    ax.view_init(elev=25, azim=35)
    plt.tight_layout()
    return fig


def plot(ind, nworkers, frame_indices, **kwargs):
    total = len(frame_indices)
    if total == 0:
        return
    start_idx = int(ind * total / nworkers)
    end_idx = int((ind + 1) * total / nworkers)
    end_idx = min(end_idx, total)
    use_matplotlib = should_use_matplotlib_backend(kwargs.get("render_backend", "auto"))
    for frame in frame_indices[start_idx:end_idx]:
        png_path = file_format.format(frame_dir, frame)
        html_path = html_format.format(frame_dir, frame)

        if use_matplotlib:
            fig = plot_system_matplotlib(frame=frame, **kwargs)
            fig.savefig(png_path, dpi=220)
            plt.close(fig)
            png_name = os.path.basename(png_path)
            with open(html_path, "w", encoding="utf-8") as f:
                f.write(
                    "<!doctype html>\n"
                    "<html><head><meta charset='utf-8'>"
                    f"<title>Frame {frame}</title>"
                    "<style>body{margin:0;background:#fff;display:flex;justify-content:center;}"
                    "img{max-width:100vw;max-height:100vh;object-fit:contain;}</style>"
                    "</head><body>"
                    f"<img src='{png_name}' alt='frame {frame}'>"
                    "</body></html>\n"
                )
        else:
            plotter = plot_system(frame=frame, **kwargs)

            # Always render PNG first.
            plotter.screenshot(
                png_path,
                window_size=(2400, 2000)   # or higher
            )

            # Try interactive HTML export; if trame is unavailable, write a
            # lightweight HTML wrapper that embeds the PNG.
            try:
                # plotter.export_vtksz('test.vtkjs')
                plotter.export_html(html_path)
            except Exception as exc:
                print(f"Warning: interactive HTML export failed for frame {frame}: {exc}")
                print("Writing static HTML wrapper around PNG instead.")
                png_name = os.path.basename(png_path)
                with open(html_path, "w", encoding="utf-8") as f:
                    f.write(
                        "<!doctype html>\n"
                        "<html><head><meta charset='utf-8'>"
                        f"<title>Frame {frame}</title>"
                        "<style>body{margin:0;background:#fff;display:flex;justify-content:center;}"
                        "img{max-width:100vw;max-height:100vh;object-fit:contain;}</style>"
                        "</head><body>"
                        f"<img src='{png_name}' alt='frame {frame}'>"
                        "</body></html>\n"
                    )
            plotter.close()


def hdf5_to_dict(hdf5_file):
    """
    Recursively converts an HDF5 file/group into a dictionary.
    This version catches conversion errors and attempts to force datasets into float64.
    """
    data_dict = {}
    def _hdf5_to_dict_recursive(group, group_name):
        for key in group.keys():
            current_path = f"{group_name}/{key}"
            if isinstance(group[key], h5py.Dataset):
                try:
                    data_dict[current_path] = np.array(group[key])
                except ValueError as e:
                    try:
                        data_dict[current_path] = np.array(group[key], dtype=np.float64)
                    except Exception as e2:
                        print(f"Failed to convert dataset {current_path}: {e2}")
            elif isinstance(group[key], h5py.Group):
                _hdf5_to_dict_recursive(group[key], current_path)
    _hdf5_to_dict_recursive(hdf5_file, "")
    return data_dict

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--filename", type=str, default='data/traj3d.h5',
                        help="HDF5 file with 3D trajectory data")
    parser.add_argument("--frame_dir", type=str, default="frames")
    parser.add_argument("--Lx", type=float, default=10)
    parser.add_argument("--Ly", type=float, default=10)
    parser.add_argument("--Lz", type=float, default=10)
    parser.add_argument("--myosin_radius", type=float, default=0.2)
    parser.add_argument("--actin_length", type=float, default=1)
    parser.add_argument("--myosin_length", type=float, default=1.5)
    parser.add_argument("--myosin_opacity", type=float, default=1.0,
                        help="Opacity for myosin filaments (1.0=opaque, <1.0 → semi-transparent).")
    parser.add_argument("--sample_every", type=int, default=1,
                        help="Only render every Nth frame (default 1 renders all frames).")
    parser.add_argument("--print_frame", type=int, default=10000)
    parser.add_argument("--start_frame", type=int, default=0,
                        help="Start frame (inclusive)")
    parser.add_argument("--end_frame", type=int, default=None,
                        help="End frame (exclusive); default is the last frame")
    parser.add_argument(
        "--myosin_display",
        choices=["all", "bonded", "cb_attached"],
        default="all",
        help="Display all myosins, only myosin–myosin bonds, or those attached to catch-bonded actins.",
    )
    parser.add_argument(
        "--actin_display",
        choices=["all", "cb"],
        default="cb",
        help="Display all actins or only those with catch-bond status == 2.",
    )
    parser.add_argument("--show_actin_force", action="store_true",
                        help="Overlay actin force vectors originating at filament centers.")
    parser.add_argument("--show_myosin_force", action="store_true",
                        help="Overlay myosin force vectors originating at filament centers.")
    parser.add_argument("--show_actin_velocity", action="store_true",
                        help="Overlay actin velocity vectors originating at filament centers.")
    parser.add_argument("--show_myosin_velocity", action="store_true",
                        help="Overlay myosin velocity vectors originating at filament centers.")
    parser.add_argument("--actin_force_scale", type=float, default=1.0,
                        help="Scale factor applied to actin force vectors.")
    parser.add_argument("--myosin_force_scale", type=float, default=1.0,
                        help="Scale factor applied to myosin force vectors.")
    parser.add_argument("--actin_velocity_scale", type=float, default=1.0,
                        help="Scale factor applied to actin velocity vectors.")
    parser.add_argument("--myosin_velocity_scale", type=float, default=1.0,
                        help="Scale factor applied to myosin velocity vectors.")
    parser.add_argument("--actin_force_thickness", type=float, default=3.0,
                        help="Multiplier for actin force arrow radii (higher → thicker).")
    parser.add_argument("--myosin_force_thickness", type=float, default=3.0,
                        help="Multiplier for myosin force arrow radii (higher → thicker).")
    parser.add_argument("--actin_velocity_thickness", type=float, default=1.5,
                        help="Multiplier for actin velocity arrow radii (higher → thicker).")
    parser.add_argument("--myosin_velocity_thickness", type=float, default=1.5,
                        help="Multiplier for myosin velocity arrow radii (higher → thicker).")
    parser.add_argument(
        "--render_backend",
        choices=["auto", "pyvista", "matplotlib"],
        default="auto",
        help="Rendering backend: auto chooses matplotlib on headless nodes without Xvfb.",
    )
    return parser.parse_args()




def analyze_myosin_pairs(data, Lx, Ly, Lz, frame_range=None):
    bonds = data["/myosin/bonds"]
    centers = data["/myosin/center"]
    dirs = data["/myosin/direction"]
    nframes = bonds.shape[0]
    box = np.array([Lx, Ly, Lz])
    if frame_range is None:
        frame_range = range(nframes)
    
    def classify_angle(angle):
        if angle < 30:
            return "parallel"
        elif angle > 150:
            return "anti-parallel"
        else:
            return "intermediate"

    for frame in frame_range:
        valid = bonds[frame]
        valid = valid[valid[:, 0] >= 0].astype(int)
        if valid.size == 0:
            continue
        print(f"\nFrame {frame}: analyzing {len(valid)} bonded pairs")
        
        for iA, iB in valid:
            cA = centers[frame, iA]
            cB = centers[frame, iB]
            dA = dirs[frame, iA]
            dB = dirs[frame, iB]
            
            # Normalize direction vectors
            dA_u = dA / np.linalg.norm(dA)
            dB_u = dB / np.linalg.norm(dB)

            # Angle between myosin directions
            cosang = np.clip(np.dot(dA_u, dB_u), -1.0, 1.0)
            angle_dirs = np.degrees(np.arccos(cosang))
            rel_dirs = classify_angle(angle_dirs)

            # Compute center-to-center vector with PBC
            delta = cB - cA
            delta -= np.round(delta / box) * box
            dist = np.linalg.norm(delta)
            delta_u = delta / (dist + 1e-12)

            # Angle between delta vector and each myosin
            angle_A = np.degrees(np.arccos(np.clip(np.dot(delta_u, dA_u), -1.0, 1.0)))
            angle_B = np.degrees(np.arccos(np.clip(np.dot(delta_u, dB_u), -1.0, 1.0)))
            rel_A = classify_angle(angle_A)
            rel_B = classify_angle(angle_B)

            print(
                f" Pair ({iA},{iB}):\n"
                f"   dir-dir angle = {angle_dirs:.1f}° → {rel_dirs}, "
                f"dist = {dist:.3f}\n"
                f"   center→A-dir = {angle_A:.1f}° → {rel_A}\n"
                f"   center→B-dir = {angle_B:.1f}° → {rel_B}"
            )


def print_bonded_myosin_info(data, start_frame=0, end_frame=None):
    """
    For each frame in the HDF5 data, prints:
      - indices of bonded myosin pairs,
      - their center positions and direction vectors.
    """
    myo_bonds = data["/myosin/bonds"]         # shape (n_frames, n_bonds_per_frame, 2)
    myo_centers = data["/myosin/center"]      # shape (n_frames, n_myosins, 3)
    myo_dirs = data["/myosin/direction"]      # shape (n_frames, n_myosins, 3)

    n_frames = myo_bonds.shape[0]
    if end_frame is None or end_frame > n_frames:
        end_frame = n_frames

    for frame in range(start_frame, end_frame):
        bonds_f = myo_bonds[frame]
        # Filter valid bonds where first index >= 0
        valid = bonds_f[bonds_f[:, 0] >= 0].astype(int)
        if valid.size == 0:
            print(f"Frame {frame}: no bonded myosins.")
            continue

        print(f"\nFrame {frame}: {len(valid)} bonded myosin pairs:")
        for idxA, idxB in valid:
            centerA = myo_centers[frame, idxA]
            dirA = myo_dirs[frame, idxA]
            centerB = myo_centers[frame, idxB]
            dirB = myo_dirs[frame, idxB]
            print(f"  Pair ({idxA}, {idxB}):")
            print(f"    Myosin {idxA} — Center: {centerA}, Direction: {dirA}")
            print(f"    Myosin {idxB} — Center: {centerB}, Direction: {dirB}")

if __name__ == "__main__":
    args = parse_args()
    filename = args.filename
    frame_dir = args.frame_dir
    Lx = args.Lx
    Ly = args.Ly
    Lz = args.Lz
    myosin_radius = args.myosin_radius
    actin_length = args.actin_length
    myosin_length = args.myosin_length
    myosin_opacity = np.clip(args.myosin_opacity, 0.0, 1.0)
    myosin_display = args.myosin_display
    actin_display = args.actin_display
    show_actin_force = args.show_actin_force
    show_myosin_force = args.show_myosin_force
    show_actin_velocity = args.show_actin_velocity
    show_myosin_velocity = args.show_myosin_velocity
    actin_force_scale = args.actin_force_scale
    myosin_force_scale = args.myosin_force_scale
    actin_velocity_scale = args.actin_velocity_scale
    myosin_velocity_scale = args.myosin_velocity_scale
    actin_force_thickness = args.actin_force_thickness
    myosin_force_thickness = args.myosin_force_thickness
    actin_velocity_thickness = args.actin_velocity_thickness
    myosin_velocity_thickness = args.myosin_velocity_thickness
    render_backend = args.render_backend

    # Open the HDF5 file and convert to dictionary.
    traj = h5py.File(filename, 'r')
    data = hdf5_to_dict(traj)
    actin_force_available = "/actin/force" in data
    myosin_force_available = "/myosin/force" in data
    actin_velocity_available = "/actin/velocity" in data
    myosin_velocity_available = "/myosin/velocity" in data
    if show_actin_force and not actin_force_available:
        print("Warning: --show_actin_force requested but /actin/force dataset is missing.")
    if show_myosin_force and not myosin_force_available:
        print("Warning: --show_myosin_force requested but /myosin/force dataset is missing.")
    if show_actin_velocity and not actin_velocity_available:
        print("Warning: --show_actin_velocity requested but /actin/velocity dataset is missing.")
    if show_myosin_velocity and not myosin_velocity_available:
        print("Warning: --show_myosin_velocity requested but /myosin/velocity dataset is missing.")
    if should_use_matplotlib_backend(render_backend):
        print("Using matplotlib fallback renderer (headless-safe).")
    else:
        print("Using PyVista renderer.")
    last_frame = data["/actin/center"].shape[0] - 1
    print_and_plot_last_frame(
        args.filename,
        args.myosin_length,
        args.Lx, args.Ly, args.Lz
    )
    print_bonded_myosin_info(data, start_frame=last_frame, end_frame=last_frame + 1)
    analyze_myosin_pairs(data, Lx, Ly, Lz, frame_range=range(last_frame, last_frame + 1))
    nframes = data["/actin/center"].shape[0]
    nparticles = data["/actin/center"].shape[1]
    print(f"Number of particles: {nparticles}")
    print(f"Number of frames: {nframes}")
    start_frame = max(0, args.start_frame)
    end_frame = args.end_frame if args.end_frame is not None else nframes
    end_frame = min(end_frame, nframes)
    sample_every = max(1, args.sample_every)
    frame_indices = list(range(start_frame, end_frame, sample_every))
    if args.print_frame < nframes:
        actin_center = data["/actin/center"][args.print_frame]
        cb_strength = data["/actin/cb_status"][args.print_frame]
        f_load = data["/actin/f_load"][args.print_frame]
        for i in range(actin_center.shape[0]):
            if cb_strength[i] > 1:
                print(f"Actin filament {i}: catch bond status: {cb_strength[i]}")
                print(f"Actin filament {i}: f load: {f_load[i]}")
                print(f"Actin filament {i}: center: {actin_center[i]}")
        #myosin_center = data["/myosin/center"][args.print_frame]
        # for i in range(myosin_center.shape[0]):
        #     print(f"Myosin filament {i}: center: {myosin_center[i]}")

    cpu_workers = joblib.cpu_count()
    print(f"Using {cpu_workers} CPU workers for parallel processing.")
    n_digits = len(str(nframes))
    file_format = "{}/frame_{:0" + str(n_digits) + "d}.png"
    html_format = "{}/frame_{:0" + str(n_digits) + "d}.html"
    if not os.path.exists(frame_dir):
        os.mkdir(frame_dir)
    Parallel(n_jobs=cpu_workers)(
        delayed(plot)(i, cpu_workers,
                      frame_indices=frame_indices,
                      data=data,
                      myosin_radius=myosin_radius,
                      actin_length=actin_length,
                      myosin_length=myosin_length,
                      Lx=Lx, Ly=Ly, Lz=Lz,
                      myosin_display=myosin_display,
                      actin_display=actin_display,
                      show_actin_force=show_actin_force,
                      show_myosin_force=show_myosin_force,
                      show_actin_velocity=show_actin_velocity,
                      show_myosin_velocity=show_myosin_velocity,
                      actin_force_scale=actin_force_scale,
                      myosin_force_scale=myosin_force_scale,
                      actin_velocity_scale=actin_velocity_scale,
                      myosin_velocity_scale=myosin_velocity_scale,
                      actin_force_available=actin_force_available,
                      myosin_force_available=myosin_force_available,
                      actin_velocity_available=actin_velocity_available,
                      myosin_velocity_available=myosin_velocity_available,
                      actin_force_thickness=actin_force_thickness,
                      myosin_force_thickness=myosin_force_thickness,
                      actin_velocity_thickness=actin_velocity_thickness,
                      myosin_velocity_thickness=myosin_velocity_thickness,
                      myosin_opacity=myosin_opacity,
                      render_backend=render_backend)
        for i in range(cpu_workers)
    )
