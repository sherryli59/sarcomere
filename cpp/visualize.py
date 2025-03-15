import argparse
import os
import h5py
import numpy as np
from ovito.data import DataCollection, Particles, Bonds
from ovito.io import export_file
from ovito.pipeline import Pipeline, PythonScriptSource


# For constructing our custom pipeline node.
from ovito.pipeline import PipelineNode, Pipeline

# --- Custom dynamic sources that returns the pre-built DataCollection for a given frame ---
class MyDynamicSource(PipelineNode):
    def __init__(self, data_list):
        # Initialize the base class with pipeline=None
        super().__init__(pipeline=None)
        self.data_list = data_list

    def compute(self, frame, **kwargs):
        # Convert frame to integer index.
        frame_index = int(frame)
        if frame_index < 0 or frame_index >= len(self.data_list):
            raise IndexError("Frame index out of range.")
        return self.data_list[frame_index]


# --- Function to build a DataCollection for one frame ---
def build_data_collection(frame, data, myosin_radius, actin_length, myosin_length, Lx, Ly, Lz):
    """
    Build a DataCollection for the given frame.
    """
    # Retrieve data for this frame.
    actin_center  = data["/actin/center"][frame]
    actin_theta   = data["/actin/theta"][frame]  # Orientation in xy-plane.
    actin_phi     = data["/actin/phi"][frame]      # Orientation in z-plane.
    cb_strength   = data["/actin/cb_strength"][frame].flatten()

    myosin_center = data["/myosin/center"][frame]
    myosin_theta  = data["/myosin/theta"][frame]
    myosin_phi    = data["/myosin/phi"][frame]

    # Create a new DataCollection.
    dc = DataCollection()
    particles = Particles()
    bonds = Bonds()

    ### Process Actin Filaments ###
    actin_positions = []
    actin_colors    = []
    actin_bonds     = []
    actin_indices   = {}

    for i in range(actin_center.shape[0]):
        x, y, z = actin_center[i]
        theta = actin_theta[i]
        phi   = actin_phi[i]

        # Compute endpoints.
        delta_x = 0.5 * actin_length * np.cos(theta) * np.cos(phi)
        delta_y = 0.5 * actin_length * np.sin(theta) * np.cos(phi)
        delta_z = 0.5 * actin_length * np.sin(phi)

        p1 = [x - delta_x, y - delta_y, z - delta_z]
        p2 = [x + delta_x, y + delta_y, z + delta_z]
        actin_positions.extend([p1, p2])
        actin_indices[i] = (len(actin_positions) - 2, len(actin_positions) - 1)

        # Use a blue–red gradient based on catch-bond strength.
        color = [cb_strength[i], 0, 1 - cb_strength[i]]
        actin_colors.extend([color, color])
        actin_bonds.append((actin_indices[i][0], actin_indices[i][1]))

    ### Process Myosin Filaments ###
    myosin_positions = []
    myosin_colors    = []
    myosin_bonds     = []
    myosin_indices   = {}

    for i in range(myosin_center.shape[0]):
        x, y, z = myosin_center[i]
        theta = myosin_theta[i]
        phi   = myosin_phi[i]

        # Compute endpoints.
        delta_x = 0.5 * myosin_length * np.cos(theta) * np.cos(phi)
        delta_y = 0.5 * myosin_length * np.sin(theta) * np.cos(phi)
        delta_z = 0.5 * myosin_length * np.sin(phi)

        p1 = [x - delta_x, y - delta_y, z - delta_z]
        p2 = [x + delta_x, y + delta_y, z + delta_z]
        myosin_positions.extend([p1, p2])
        myosin_indices[i] = (len(myosin_positions) - 2, len(myosin_positions) - 1)

        # Myosin is always orange.
        myosin_colors.extend([[1, 0.5, 0], [1, 0.5, 0]])
        myosin_bonds.append((myosin_indices[i][0], myosin_indices[i][1]))

    # Combine actin and myosin data.
    all_positions = np.vstack((actin_positions, myosin_positions))
    all_colors    = np.vstack((actin_colors, myosin_colors))
    all_bonds     = np.array(actin_bonds + myosin_bonds)

    # Add properties.
    particles.create_property('Position', data=all_positions)
    particles.create_property('Color', data=all_colors)
    bonds.create_property('Topology', data=all_bonds)

    dc.objects.append(particles)
    dc.objects.append(bonds)
    return dc


# --- Helper to load all HDF5 data into a dictionary ---
def hdf5_to_dict(hdf5_file):
    data_dict = {}
    def _hdf5_to_dict_recursive(group, group_name):
        for key in group.keys():
            if isinstance(group[key], h5py.Dataset):
                data_dict[f"{group_name}/{key}"] = np.array(group[key].astype(np.float64)).squeeze()
            elif isinstance(group[key], h5py.Group):
                _hdf5_to_dict_recursive(group[key], f"{group_name}/{key}")
    _hdf5_to_dict_recursive(hdf5_file, "")
    return data_dict


# --- Argument parser ---
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--filename", type=str, default='data/traj.h5',
                        help="Path to the HDF5 trajectory file.")
    parser.add_argument("--output_file", type=str, default="output.ovito",
                        help="Name for the output OVITO session file.")
    parser.add_argument("--Lx", type=float, default=10)
    parser.add_argument("--Ly", type=float, default=10)
    parser.add_argument("--Lz", type=float, default=10)
    parser.add_argument("--myosin_radius", type=float, default=0.2)
    parser.add_argument("--actin_length", type=float, default=1)
    parser.add_argument("--myosin_length", type=float, default=1.5)
    return parser.parse_args()


# --- Main execution ---
if __name__ == "__main__":
    args = parse_args()

    # Load HDF5 data.
    with h5py.File(args.filename, "r") as f:
        data = hdf5_to_dict(f)
    nframes = data["/actin/center"].shape[0]
    print(f"Total frames: {nframes}")
    #pring actin center
    last_frame = data["/myosin/phi"][-1]
    print(last_frame)
    # Build a list of DataCollections – one for each frame.
    data_list = []
    for frame in range(nframes):
        dc = build_data_collection(frame, data, args.myosin_radius,
                                   args.actin_length, args.myosin_length,
                                   args.Lx, args.Ly, args.Lz)
        data_list.append(dc)

    # Define the create function for the PythonScriptSource
    def create(frame, data_collection):
        dc = build_data_collection(frame, data, args.myosin_radius, args.actin_length, args.myosin_length,args.Lx, args.Ly, args.Lz)
        data_collection.objects.extend(dc.objects)

    # Create a pipeline with the PythonScriptSource
    pipeline = Pipeline(source=PythonScriptSource(function=create))

    # Create a scene and add the pipeline.
    scene = ovito.scene.Scene()
    scene.add_pipeline(pipeline)

    # Export the entire scene as an OVITO session file.
    # The output format identifier is "ovito/session".
    export_file(scene, args.output_file, "ovito/session", multiple_frames=True)
    print(f"Saved OVITO session to '{args.output_file}'")
