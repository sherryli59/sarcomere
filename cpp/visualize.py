import h5py
import numpy as np
import ovito
from ovito.data import DataCollection, Particles, Bonds
from ovito.io import export_file
import os
import argparse
import joblib
from joblib import Parallel, delayed


def create_ovito_visualization(frame, data, myosin_radius, actin_length, myosin_length, Lx, Ly, Lz, output_dir):
    """
    Generates a 3D OVITO visualization where actin and myosin are properly rendered as 3D **filaments**.
    """
    
    # Load actin and myosin data for the given frame
    actin_center = data["/actin/center"][frame]
    actin_theta = data["/actin/theta"][frame]  # Orientation in xy-plane
    actin_phi = data["/actin/phi"][frame]  # Orientation in z-plane
    cb_strength = data["/actin/cb_strength"][frame].flatten()

    myosin_center = data["/myosin/center"][frame]
    myosin_theta = data["/myosin/theta"][frame]
    myosin_phi = data["/myosin/phi"][frame]

    # Create OVITO data collection
    data_collection = DataCollection()
    particles = Particles()
    positions = []
    bonds = Bonds()

    # Set up color mapping for catch-bond strength
    color_map = np.clip(cb_strength, 0, 1)

    ### Add Actin Filaments ###
    actin_positions = []
    actin_colors = []
    actin_bonds = []
    actin_indices = {}

    for i in range(actin_center.shape[0]):
        x, y, z = actin_center[i]
        theta = actin_theta[i]
        phi = actin_phi[i]
        
        # Compute filament endpoints using theta (xy-plane) and phi (z-plane)
        delta_x = 0.5 * actin_length * np.cos(theta) * np.cos(phi)
        delta_y = 0.5 * actin_length * np.sin(theta) * np.cos(phi)
        delta_z = 0.5 * actin_length * np.sin(phi)

        # Two endpoints per filament
        p1 = [x - delta_x, y - delta_y, z - delta_z]
        p2 = [x + delta_x, y + delta_y, z + delta_z]
        
        actin_positions.extend([p1, p2])
        actin_indices[i] = (len(actin_positions) - 2, len(actin_positions) - 1)

        # Color based on catch bond strength (blue-red gradient)
        actin_colors.extend([[cb_strength[i], 0, 1 - cb_strength[i]], [cb_strength[i], 0, 1 - cb_strength[i]]])

        # Add bond connecting the two endpoints
        actin_bonds.append((actin_indices[i][0], actin_indices[i][1]))

    ### Add Myosin Filaments ###
    myosin_positions = []
    myosin_colors = []
    myosin_bonds = []
    myosin_indices = {}

    for i in range(myosin_center.shape[0]):
        x, y, z = myosin_center[i]
        theta = myosin_theta[i]
        phi = myosin_phi[i]
        
        # Compute filament endpoints using theta (xy-plane) and phi (z-plane)
        delta_x = 0.5 * myosin_length * np.cos(theta) * np.cos(phi)
        delta_y = 0.5 * myosin_length * np.sin(theta) * np.cos(phi)
        delta_z = 0.5 * myosin_length * np.sin(phi)

        # Two endpoints per filament
        p1 = [x - delta_x, y - delta_y, z - delta_z]
        p2 = [x + delta_x, y + delta_y, z + delta_z]

        myosin_positions.extend([p1, p2])
        myosin_indices[i] = (len(myosin_positions) - 2, len(myosin_positions) - 1)

        # Myosin is always orange
        myosin_colors.extend([[1, 0.5, 0], [1, 0.5, 0]])

        # Add bond connecting the two endpoints
        myosin_bonds.append((myosin_indices[i][0], myosin_indices[i][1]))

    # Convert positions into NumPy arrays
    all_positions = np.vstack((actin_positions, myosin_positions))
    all_colors = np.vstack((actin_colors, myosin_colors))
    all_bonds = np.array(actin_bonds + myosin_bonds)

    # Add particles (filament endpoints)
    particles.create_property('Position', data=all_positions)
    particles.create_property('Color', data=all_colors)

    # Add bonds between endpoints to make them **filaments**
    bonds.create_property('Topology', data=all_bonds)

    # Add to data collection
    data_collection.add(particles)
    data_collection.add(bonds)

    # Export OVITO file
    output_file = os.path.join(output_dir, f"frame_{frame:04d}.ovito")
    export_file(data_collection, output_file, "ovito")
    print(f"Saved OVITO frame: {output_file}")


# Argument Parser
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--filename", type=str, default='data/traj.h5')
    parser.add_argument("--output_dir", type=str, default="frames_3D")
    parser.add_argument("--Lx", type=float, default=10)
    parser.add_argument("--Ly", type=float, default=10)
    parser.add_argument("--Lz", type=float, default=10)  # 3D space
    parser.add_argument("--myosin_radius", type=float, default=0.2)
    parser.add_argument("--actin_length", type=float, default=1)
    parser.add_argument("--myosin_length", type=float, default=1.5)
    parser.add_argument("--n_jobs", type=int, default=-1)  # Use all CPUs
    return parser.parse_args()


# Main Execution
if __name__ == "__main__":
    args = parse_args()
    filename = args.filename
    output_dir = args.output_dir

    # Load trajectory file
    traj = h5py.File(filename, 'r')
    data = {key: np.array(traj[key]) for key in traj.keys()}

    nframes = data["/actin/center"].shape[0]
    print(f"Total frames: {nframes}")

    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Parallel processing with Joblib
    Parallel(n_jobs=args.n_jobs)(
        delayed(create_ovito_visualization)(
            frame, data, args.myosin_radius, args.actin_length, args.myosin_length, args.Lx, args.Ly, args.Lz, output_dir
        ) for frame in range(nframes)
    )
