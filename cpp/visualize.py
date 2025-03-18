import argparse
import os
import h5py
import numpy as np
from ovito.data import DataCollection, Particles, Bonds, SimulationCell
from ovito.pipeline import Pipeline, PythonScriptSource
import ovito
import joblib
from joblib import Parallel, delayed
from ovito.vis import Viewport, OSPRayRenderer,TachyonRenderer, CoordinateTripodOverlay
from ovito.qt_compat import QtCore


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
    actin_theta   = data["/actin/theta"][frame]
    actin_phi     = data["/actin/phi"][frame]
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
    actin_colors = []
    actin_bonds = []
    actin_indices = {}

    for i in range(actin_center.shape[0]):
        x, y, z = actin_center[i]
        theta = actin_theta[i]
        phi = actin_phi[i]

        delta_x = 0.5 * actin_length * np.cos(theta) * np.sin(phi)
        delta_y = 0.5 * actin_length * np.sin(theta) * np.sin(phi)
        delta_z = 0.5 * actin_length * np.cos(phi)


        p1 = [x - delta_x, y - delta_y, z - delta_z]
        p2 = [x + delta_x, y + delta_y, z + delta_z]
        actin_positions.extend([p1, p2])
        actin_indices[i] = (len(actin_positions) - 2, len(actin_positions) - 1)

        color = [cb_strength[i], 0, 1 - cb_strength[i]]
        actin_colors.extend([color, color])
        actin_bonds.append((actin_indices[i][0], actin_indices[i][1]))

    ### Process Myosin Filaments ###
    myosin_positions = []
    myosin_colors = []
    myosin_bonds = []
    myosin_indices = {}

    for i in range(myosin_center.shape[0]):
        x, y, z = myosin_center[i]
        theta = myosin_theta[i]
        phi = myosin_phi[i]

        delta_x = 0.5 * myosin_length * np.cos(theta) * np.sin(phi)
        delta_y = 0.5 * myosin_length * np.sin(theta) * np.sin(phi)
        delta_z = 0.5 * myosin_length * np.cos(phi)

        p1 = [x - delta_x, y - delta_y, z - delta_z]
        p2 = [x + delta_x, y + delta_y, z + delta_z]
        myosin_positions.extend([p1, p2])
        myosin_indices[i] = (len(myosin_positions) - 2, len(myosin_positions) - 1)

        myosin_colors.extend([[1, 0.5, 0], [1, 0.5, 0]])
        myosin_bonds.append((myosin_indices[i][0], myosin_indices[i][1]))

        # Calculate offset for myosin bonds
    n_actin_particles = len(actin_positions)
    
    # Adjust myosin bond indices for global particle list
    myosin_bonds_adjusted = [
        (idx1 + n_actin_particles, idx2 + n_actin_particles)
        for (idx1, idx2) in myosin_bonds
    ]
    
    # # Create final bond list
    # all_bonds = np.array(actin_bonds + myosin_bonds_adjusted)

    # # Set radii explicitly
    # actin_bond_radii = np.full(len(actin_bonds), 0.01) 
    # myosin_bond_radii = np.full(len(myosin_bonds), myosin_radius)  
    # all_bond_radii = np.concatenate([actin_bond_radii, myosin_bond_radii])


    all_bonds = np.array(myosin_bonds_adjusted)

    # Set radii explicitly
    actin_bond_radii = np.full(len(actin_bonds), 0.01) 
    myosin_bond_radii = np.full(len(myosin_bonds), myosin_radius)  
    all_bond_radii = myosin_bond_radii  

    # Combine data
    all_positions = np.vstack((actin_positions, myosin_positions))
    all_colors = np.vstack((actin_colors, myosin_colors))
    # Create radii arrays
    actin_radii = np.full((len(actin_positions),), 0.01)
    myosin_radii = np.full((len(myosin_positions),), myosin_radius)
    all_radii = np.concatenate((actin_radii, myosin_radii))

    # Add properties to particles and bonds
    particles.create_property('Position', data=all_positions)
    particles.create_property('Color', data=all_colors)
    particles.create_property('Radius', data=all_radii)
    bonds.create_property('Topology', data=all_bonds)
    bonds.create_property('Radius', data=all_bond_radii)  
    bonds.vis.enabled = True
    # bonds.vis.shading = BondsVis.Shading.Flat
    # bonds.vis.radius_mapping.enabled = True  # <-- Critical for using Radius property
    # bonds.vis.radius_mapping.radius = 1.0    # <-- Scaling factor for your radii
    # Assign bonds to particles
    particles.bonds = bonds

    # Create simulation cell
    cell = SimulationCell()
    cell.matrix = [
        [Lx, 0, 0, -Lx/2],  # Center box at origin
        [0, Ly, 0, -Ly/2],
        [0, 0, Lz, -Lz/2]
    ]
    cell.pbc = (True, True, True)  # Set periodic boundary conditions as needed

    dc.objects.append(cell)
    dc.objects.append(particles)

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

def render_image(data, data_path, ind, n_workers):
    nframes = data["/actin/center"].shape[0]
    frame_start = int(ind*(nframes-1)/n_workers)
    frame_end = int((ind+1)*(nframes-1)/n_workers)
    if ind == n_workers-1:
        frame_end = nframes
    # Build a list of DataCollections – one for each frame.
    data_list = []
    for frame in range(frame_start, frame_end):
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
    pipeline.compute()
    pipeline.add_to_scene()
    # cell_vis = pipeline.source.data.cell.vis
    # cell_vis.enabled = False
    vp = Viewport(type = Viewport.Type.Ortho)
    vp.camera_pos = (-10, -15, 15)
    vp.camera_dir = (2, 3, -3)
    vp.fov = 3.5
    vp.zoom_all()

    # Create the overlay.
    tripod = CoordinateTripodOverlay()
    tripod.size = 0.07
    tripod.alignment = QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignBottom
    vp.overlays.append(tripod)
    for frame in range(frame_start, frame_end):
        #make directory if it doesn't exist
        vp.render_image(filename=data_path+f"frame_{frame:04d}.png", 
                        size=(480,240), background=(1,1,1), frame=frame,
                        renderer=TachyonRenderer())

# --- Main execution ---
if __name__ == "__main__":
    args = parse_args()

    # Load HDF5 data.
    with h5py.File(args.filename, "r") as f:
        data = hdf5_to_dict(f)
    nframes = data["/actin/center"].shape[0]
    print(f"Total frames: {nframes}")
    #pring actin center
    myosin_centers = data["/myosin/center"][5:]
    # compute center-center distance by doing the unsqueeze operation
    center_center_distance = np.linalg.norm(np.expand_dims(myosin_centers, axis=-2) - np.expand_dims(myosin_centers, axis=-3), axis=-1)
    center_center_distance = center_center_distance[center_center_distance<0.4]
    center_center_distance = center_center_distance[center_center_distance>0]
    print(center_center_distance)
    cb_strength   = data["/actin/cb_strength"]
    print(np.average(cb_strength, axis=-1))
    exit()
    resolution = (1920, 1080)
    output_dir = "frames/"
    os.makedirs(output_dir, exist_ok=True)
    cpu_workers = joblib.cpu_count()
    Parallel(n_jobs=cpu_workers)(delayed(render_image)(data, output_dir,ind, cpu_workers) for ind in range(cpu_workers))

