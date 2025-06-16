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

def plot_filaments_3d(center, direction, l, Lx, Ly, Lz, ax, color='k', color_spectrum=None, **kwargs):
    """
    Plot actin filaments in 3D using unit direction vectors.
    - center: (N,3) array of filament centers.
    - direction: (N,3) array of unit direction vectors.
    - l: filament length (scalar or array).
    - Lx, Ly, Lz: periodic box dimensions.
    - ax: a matplotlib 3D axis.
    - color_spectrum: optional array for color mapping.
    """
    if np.isscalar(l):
        l = np.ones(center.shape[0]) * l

    if color_spectrum is not None:
        color_spectrum = np.sqrt(color_spectrum)
        colormap = plt.cm.Blues
        custom_colormap = mcolors.LinearSegmentedColormap.from_list(
            "custom_blues", colormap(np.linspace(0.3, 1, 256))
        )
        norm = plt.Normalize(0, 1)
        sm = plt.cm.ScalarMappable(cmap=custom_colormap, norm=norm)
        sm.set_array([])

    for i in range(center.shape[0]):
        if l[i] < 0.01:
            continue

        this_color = color
        if color_spectrum is not None: # and color_spectrum[i] >= 0.01:
            this_color = custom_colormap(norm(color_spectrum[i]))
            kwargs['alpha'] = 0.3 + 0.7 * (color_spectrum[i] > 0)
        # elif color_spectrum is not None:
        #     continue

        x, y, z = center[i]
        ux, uy, uz = direction[i]
        dx = 0.5 * l[i] * ux
        dy = 0.5 * l[i] * uy
        dz = 0.5 * l[i] * uz

        ax.quiver(x - dx, y - dy, z - dz, 2*dx, 2*dy, 2*dz, arrow_length_ratio=0.1,
                  color=this_color, **kwargs)

        # PBC duplicates
        ox_set = get_offsets_for_filament(x, abs(dx), Lx)
        oy_set = get_offsets_for_filament(y, abs(dy), Ly)
        oz_set = get_offsets_for_filament(z, abs(dz), Lz)
        for ox, oy, oz in product(ox_set, oy_set, oz_set):
            if ox == oy == oz == 0:
                continue
            ax.quiver(x - dx + ox, y - dy + oy, z - dz + oz,
                      2*dx, 2*dy, 2*dz, arrow_length_ratio=0.1,
                      color=this_color, **kwargs)

    if color_spectrum is not None:
        cbar = plt.colorbar(sm, ax=ax, orientation='vertical', fraction=0.02, pad=0.04)
        cbar.set_label('Catch bond strength')




def plot_myosin_velocity(center, velocity, Lx, Ly, Lz, ax, scale=1.0, color='hotpink'):
    """
    Plot velocity vectors as arrows originating from myosin centers.
    """
    for i in range(center.shape[0]):
        x, y, z = center[i]
        vx, vy, vz = velocity[i] * scale  # apply optional scale factor

        ax.quiver(x, y, z, vx, vy, vz,
                  color=color, arrow_length_ratio=0.2, linewidth=1.2)

def compute_alternative_torque(actin_center, actin_direction, myosin_center, myosin_direction,
                                actin_length, myosin_length, box, kappa_am, myosin_radius):
    """
    Compute alternative myosin torques based on angular misalignment with nearby actins.
    Returns a (Nm, 3) array of torques.
    """
    Nm = myosin_center.shape[0]
    Na = actin_center.shape[0]
    torque_alt = np.zeros_like(myosin_direction)

    for i in range(Nm):
        ui = myosin_direction[i]

        for j in range(Na):
            uj = actin_direction[j]

            # Segment-segment distance filtering (approximate; could replace with fast check)
            r1 = myosin_center[i]
            r2 = actin_center[j]
            dr = r2 - r1

            # Apply minimum image convention (PBC)
            for d in range(3):
                if dr[d] > 0.5 * box[d]:
                    dr[d] -= box[d]
                elif dr[d] < -0.5 * box[d]:
                    dr[d] += box[d]

            # dist = np.linalg.norm(dr)
            # if dist > myosin_radius:
            #     print(f"Skipping myosin {i} and actin {j}: distance {dist} exceeds radius {myosin_radius}")
            #     continue  # Only consider nearby actins

            dot = np.clip(np.dot(ui, uj), -1.0, 1.0)
            if np.abs(dot) >= 0.999:
                continue  # Already aligned or anti-aligned

            theta = np.arccos(dot)
            prefactor = kappa_am * theta / np.sqrt(1 - dot ** 2)
            grad = prefactor * uj  # ∂E/∂ui
            torque_contrib = np.cross(grad, ui)  # project to tangent space
            print("torque contrib:",torque_contrib)
            torque_alt[i] += torque_contrib

    return torque_alt

def plot_system(frame, data, myosin_length, actin_length, Lx, Ly, Lz, myosin_radius):
    """
    Create a 3D plot for a given frame.
    Assumes:
      - Actin center is stored in "/actin/center" as an (N,3) array.
      - Actin direction is stored in "/actin/direction" as an (N,3) array.
      - Myosin data follows a similar convention.
    """
    fig = plt.figure(dpi=500)
    ax = fig.add_subplot(111, projection='3d')

    # Actin data.
    actin_center = data["/actin/center"][frame]
    actin_direction = data["/actin/direction"][frame]
    cb_strength = data["/actin/cb_strength"][frame].flatten()
    actin_direction = data["/actin/direction"][frame]  # shape (N, 3)
    plot_filaments_3d(actin_center, actin_direction, actin_length,
                  Lx, Ly, Lz, ax, color='k', color_spectrum=cb_strength)
    
    # Myosin data.
    myosin_center = data["/myosin/center"][frame]
    myosin_direction = data["/myosin/direction"][frame]
    plot_filaments_3d(myosin_center, myosin_direction, myosin_length,
                   Lx, Ly, Lz, ax, color='C1')
    
    myosin_velocity = data["/myosin/velocity"][frame]
    # #print any non-zero velocities
    # if np.any(myosin_velocity > 0.0001):
    #     print(f"Frame {frame}: Myosin velocities:")
    #     for i in range(myosin_velocity.shape[0]):
    #         if np.linalg.norm(myosin_velocity[i]) > 0.0001:
    #             print(f"Myosin {i}: velocity: {myosin_velocity[i]}")
    #             print(f"Myosin {i}: center: {myosin_center[i]}")
    plot_myosin_velocity(myosin_center, myosin_velocity, Lx, Ly, Lz, ax, scale=5, color='hotpink')

    myosin_torque = data["/myosin/torque"][frame]
    print(myosin_torque)
    plot_myosin_velocity(myosin_center, myosin_torque, Lx, Ly, Lz, ax, scale=0.5, color='C2')
    box = np.array([Lx, Ly, Lz])
    myosin_torque_alt = compute_alternative_torque(
    actin_center, actin_direction,
    myosin_center, myosin_direction,
    actin_length, myosin_length,
    box, kappa_am=1.0,  # or your simulation value
    myosin_radius=myosin_radius
    )
    #plot_myosin_velocity(myosin_center, myosin_torque_alt, Lx, Ly, Lz, ax, scale=0.5, color='lime')
    ax.view_init(elev=20, azim=30)
    ax.set_xlim(-Lx/2, Lx/2)
    ax.set_ylim(-Ly/2, Ly/2)
    ax.set_zlim(-Lz/2, Lz/2)
    ax.set_box_aspect([Lx, Ly, Lz])
    return fig, ax

def plot(ind, nframes, nworkers, **kwargs):
    frame_start = int(ind * nframes / nworkers)
    frame_end = max(int((ind + 1) * nframes / nworkers), 1)
    for frame in range(frame_start, frame_end):
        fig, ax = plot_system(frame=frame, **kwargs)
        plt.savefig(file_format.format(frame_dir, frame))
        plt.close(fig)

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
    return parser.parse_args()

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

    nframes = data["/actin/center"].shape[0]
    nparticles = data["/actin/center"].shape[1]
    print(f"Number of particles: {nparticles}")
    print(f"Number of frames: {nframes}")
    if args.print_frame < nframes:
        actin_center = data["/actin/center"][args.print_frame]
        cb_strength = data["/actin/cb_strength"][args.print_frame]
        for i in range(actin_center.shape[0]):
            if cb_strength[i] > 0.01:
                print(f"Actin filament {i}: catch bond strength: {cb_strength[i]}")
                print(f"Actin filament {i}: center: {actin_center[i]}")
        myosin_center = data["/myosin/center"][args.print_frame]
        for i in range(myosin_center.shape[0]):
            print(f"Myosin filament {i}: center: {myosin_center[i]}")

    #cpu_workers = joblib.cpu_count()
    cpu_workers = 1
    n_digits = len(str(nframes))
    file_format = "{}/frame_{:0" + str(n_digits) + "d}.png"
    if not os.path.exists(frame_dir):
        os.mkdir(frame_dir)
    Parallel(n_jobs=cpu_workers)(
        delayed(plot)(i, nframes, cpu_workers,
                      data=data,
                      myosin_radius=myosin_radius,
                      actin_length=actin_length,
                      myosin_length=myosin_length,
                      Lx=Lx, Ly=Ly, Lz=Lz)
        for i in range(cpu_workers)
    )