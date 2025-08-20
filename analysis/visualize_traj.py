import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import math
import joblib
from joblib import Parallel, delayed
import os
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import sys
import argparse
from mpl_toolkits.mplot3d import Axes3D  # required for 3D plotting
from itertools import product
import pyvista as pv
pv.start_xvfb()  # Start Xvfb for off-screen rendering


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

    if color_spectrum is not None:
        color_spectrum = np.sqrt(color_spectrum)
        colormap = plt.cm.Blues
        custom_colormap = plt.cm.ScalarMappable(cmap=plt.cm.Blues)
        norm = plt.Normalize(0, 1)

    for i in range(center.shape[0]):
        if l[i] < 0.01:
            continue

        this_color = color
        if color_spectrum is not None and color_spectrum[i] >= 0.01:
            this_color = plt.cm.Blues(norm(color_spectrum[i]))
        elif color_spectrum is not None:
            continue

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




def plot_system(frame, data, myosin_length, actin_length, Lx, Ly, Lz, myosin_radius):
    plotter = pv.Plotter(off_screen=True)
    # Actin
    actin_center = data["/actin/center"][frame]
    actin_direction = data["/actin/direction"][frame]
    cb_strength = data["/actin/cb_strength"][frame].flatten()
    f_load = data["/actin/f_load"][frame].flatten()
    plot_filaments_3d(
        center=actin_center,
        direction=actin_direction,
        radius=0.02,
        l=actin_length,
        Lx=Lx, Ly=Ly, Lz=Lz,
        plotter=plotter,
        color='blue',
        color_spectrum=cb_strength * f_load
    )

    # Myosins engaged in catch bonds
    bonds = data["/myosin/bonds"][frame]
    valid_pairs = bonds[bonds[:, 0] >= 0].astype(int)
    bonded_indices = np.unique(valid_pairs.flatten())
    if bonded_indices.size > 0:
        myosin_center = data["/myosin/center"][frame][bonded_indices]
        myosin_direction = data["/myosin/direction"][frame][bonded_indices]

        plot_filaments_3d(
            center=myosin_center,
            direction=myosin_direction,
            radius=myosin_radius,
            l=myosin_length,
            Lx=Lx, Ly=Ly, Lz=Lz,
            plotter=plotter,
            color='lemon_chiffon',
        )

    plotter.set_background("white")
    plotter.set_scale(xscale=Lx, yscale=Ly, zscale=Lz)
    plotter.set_focus((0, 0, 0))
    plotter.camera_position = 'xy'

    return plotter


def plot(ind, nworkers, frame_range, **kwargs):
    frame_start = frame_range.start + int(ind * len(frame_range) / nworkers)
    frame_end = frame_range.start + int((ind + 1) * len(frame_range) / nworkers)
    frame_end = min(frame_end, frame_range.stop)
    for frame in range(frame_start, frame_end):
        plotter = plot_system(frame=frame, **kwargs)
        #plotter.export_vtksz('test.vtkjs')
        plotter.export_html(html_format.format(frame_dir, frame))
        plotter.screenshot(file_format.format(frame_dir, frame))
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
    parser.add_argument("--print_frame", type=int, default=10000)
    parser.add_argument("--start_frame", type=int, default=0,
                        help="Start frame (inclusive)")
    parser.add_argument("--end_frame", type=int, default=None,
                        help="End frame (exclusive); default is the last frame")
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

    # Open the HDF5 file and convert to dictionary.
    traj = h5py.File(filename, 'r')
    data = hdf5_to_dict(traj)
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
    frame_range = range(start_frame, end_frame)
    if args.print_frame < nframes:
        actin_center = data["/actin/center"][args.print_frame]
        cb_strength = data["/actin/cb_strength"][args.print_frame]
        f_load = data["/actin/f_load"][args.print_frame]
        for i in range(actin_center.shape[0]):
            if cb_strength[i] > 0.01:
                print(f"Actin filament {i}: catch bond strength: {cb_strength[i]}")
                print(f"Actin filament {i}: f load: {f_load[i]}")
                print(f"Actin filament {i}: center: {actin_center[i]}")
        myosin_center = data["/myosin/center"][args.print_frame]
        for i in range(myosin_center.shape[0]):
            print(f"Myosin filament {i}: center: {myosin_center[i]}")

    cpu_workers = joblib.cpu_count()
    print(f"Using {cpu_workers} CPU workers for parallel processing.")
    n_digits = len(str(nframes))
    file_format = "{}/frame_{:0" + str(n_digits) + "d}.png"
    html_format = "{}/frame_{:0" + str(n_digits) + "d}.html"
    if not os.path.exists(frame_dir):
        os.mkdir(frame_dir)
    Parallel(n_jobs=cpu_workers)(
        delayed(plot)(i, cpu_workers,
                      frame_range=frame_range,
                      data=data,
                      myosin_radius=myosin_radius,
                      actin_length=actin_length,
                      myosin_length=myosin_length,
                      Lx=Lx, Ly=Ly, Lz=Lz)
        for i in range(cpu_workers)
    )



