import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import joblib
from joblib import Parallel, delayed
import os
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import argparse
from itertools import product

# Attempt to import OVITO
try:
    from ovito.data import DataCollection, Particles, Bonds, SimulationCell
    from ovito.vis import Viewport, TachyonRenderer, OpenGLRenderer
    from ovito.pipeline import Pipeline, StaticSource
except ImportError:
    print("Warning: ovito is not installed in this environment.")
    print("To run this script, install ovito (e.g., via conda-forge: micromamba install -c conda-forge ovito)")


def hdf5_to_dict(hdf5_file):
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

def get_offsets_for_filament(coord, d, L):
    offsets = {0}
    p1 = coord - d
    p2 = coord + d
    if p1 < -L/2:
        offsets.add(L)
    if p2 > L/2:
        offsets.add(-L)
    return offsets

def hex_to_rgb(hex_color):
    hex_color = hex_color.lstrip('#')
    return [int(hex_color[i:i+2], 16) / 255.0 for i in (0, 2, 4)]

def build_scene(frame, data, myosin_length, actin_length, Lx, Ly, Lz,
                myosin_radius, myosin_display="all", actin_display="cb"):
    
    positions = []
    topology = []
    colors = []
    widths = []

    pos_idx = 0

    def add_filament(center, direction, length, radius, color):
        nonlocal pos_idx
        
        dx = 0.5 * length * np.array(direction)
        ox_set = get_offsets_for_filament(center[0], abs(dx[0]), Lx)
        oy_set = get_offsets_for_filament(center[1], abs(dx[1]), Ly)
        oz_set = get_offsets_for_filament(center[2], abs(dx[2]), Lz)
        
        for ox, oy, oz in product(ox_set, oy_set, oz_set):
            duplicate_center = center + np.array([ox, oy, oz])
            p1 = duplicate_center - dx
            p2 = duplicate_center + dx
            
            positions.append(p1)
            positions.append(p2)
            topology.append([pos_idx, pos_idx + 1])
            colors.append(color)
            widths.append(radius)
            pos_idx += 2

    # ------------------ Actin ------------------
    actin_center = data["/actin/center"][frame]
    actin_direction = data["/actin/direction"][frame]
    f_load = data["/actin/f_load"][frame].flatten()
    cb_status = data["/actin/cb_status"][frame].flatten()

    load_metric = np.maximum(f_load, 0) * (cb_status == 2)
    if actin_display == "cb":
        mask = cb_status == 2
        print(f"Frame {frame}: {np.sum(mask)} actin filaments in catch-bond state")
    else:
        mask = np.ones_like(cb_status, dtype=bool)

    actin_center = actin_center[mask]
    actin_direction = actin_direction[mask]
    actin_color_values = load_metric[mask]

    # Precompute norm for colors
    max_val = np.max(actin_color_values) if actin_color_values.size > 0 else 0
    if max_val > 0:
        norm_color_vals = np.sqrt(actin_color_values / max_val)
    else:
        norm_color_vals = np.zeros_like(actin_color_values)

    cmap = plt.cm.Blues

    for i in range(actin_center.shape[0]):
        c_val = norm_color_vals[i]
        color = cmap(c_val)[:3]  # RGB
        add_filament(actin_center[i], actin_direction[i], actin_length, 0.01, color)

    # ------------------ Myosin ------------------
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
        display_indices = np.unique(valid_pairs.flatten()) if valid_pairs.size > 0 else np.empty(0, dtype=int)
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

    lemon_chiffon = hex_to_rgb('#FFFACD')
    light_orange = hex_to_rgb('#F5A45B')

    if display_indices.size > 0:
        myosin_center = myosin_centers_frame[display_indices]
        myosin_direction = myosin_dirs_frame[display_indices]
        
        for i, m_idx in enumerate(display_indices):
            color = light_orange if m_idx in highlight_indices else lemon_chiffon
            add_filament(myosin_center[i], myosin_direction[i], myosin_length, myosin_radius, color)

    # Convert to OVITO arrays
    positions = np.array(positions, dtype=np.float32) if len(positions) > 0 else np.empty((0, 3), dtype=np.float32)
    topology = np.array(topology, dtype=int) if len(topology) > 0 else np.empty((0, 2), dtype=int)
    colors = np.array(colors, dtype=np.float32) if len(colors) > 0 else np.empty((0, 3), dtype=np.float32)
    widths = np.array(widths, dtype=np.float32) if len(widths) > 0 else np.empty((0,), dtype=np.float32)
    
    return positions, topology, colors, widths


def plot_frame_ovito(frame, frame_dir, filename_format, data, myosin_length, actin_length, Lx, Ly, Lz, myosin_radius, myosin_display, actin_display):
    try:
        from ovito.data import DataCollection, Particles, Bonds, SimulationCell
        from ovito.vis import Viewport, TachyonRenderer, OpenGLRenderer
        from ovito.pipeline import Pipeline, StaticSource
        from ovito.modifiers import ConstructSurfaceModifier
    except ImportError:
        return
        
    positions, topology, colors, widths = build_scene(
        frame=frame, data=data, myosin_length=myosin_length,
        actin_length=actin_length, Lx=Lx, Ly=Ly, Lz=Lz,
        myosin_radius=myosin_radius, myosin_display=myosin_display,
        actin_display=actin_display
    )

    dc = DataCollection()
    
    # Setup cell
    cell = SimulationCell(pbc=(True, True, True))
    # Box matrix in OVITO: column-major, [x, y, z, origin]
    cell_matrix = np.array([
        [Lx, 0, 0, -Lx/2],
        [0, Ly, 0, -Ly/2],
        [0, 0, Lz, -Lz/2]
    ])
    cell.matrix = cell_matrix
    dc.objects.append(cell)
    
    # Setup particles
    particles = Particles()
    particles.create_property('Position', data=positions)
    # Hide particles by setting their radius to 0
    particles.create_property('Radius', data=np.zeros(len(positions)))
    
    # Setup bonds
    if len(topology) > 0:
        bonds = Bonds()
        bonds.create_property('Topology', data=topology)
        bonds.create_property('Color', data=colors)
        bonds.create_property('Width', data=widths)
        particles.bonds = bonds
        
    dc.objects.append(particles)
    
    # Create Pipeline
    pipeline = Pipeline(source=StaticSource(data=dc))
    
    # Setup styling for visual representation
    # By default, OVITO binds Bond width to the Width property
    # Let's explicitly configure the ParticlesVis and BondsVis
    particles_vis = pipeline.compute().particles.vis
    particles_vis.shape = particles_vis.Shape.Sphere
    particles_vis.radius = 0.0 # Make nodes invisible safely

    bonds_vis = pipeline.compute().particles.bonds.vis
    bonds_vis.width = 0.1 # This is fallback, it should use the Width property if present
    bonds_vis.flat_shading = False
    
    # Set up viewport
    vp = Viewport()
    vp.type = Viewport.Type.Perspective
    vp.camera_pos = (15, 15, 15)
    vp.camera_dir = (-15, -15, -15)
    vp.fov = np.radians(35)
    
    png_path = filename_format.format(frame_dir, frame)
    
    # Try different rendering engines depending on stability
    try:
        vp.render_image(filename=png_path, size=(2400, 2000), 
                        renderer=TachyonRenderer(), 
                        background=(1.0, 1.0, 1.0))
    except Exception as e:
        print(f"Frame {frame} Tahyon render failed: {e}. Trying OpenGLRenderer.")
        vp.render_image(filename=png_path, size=(2400, 2000), 
                        renderer=OpenGLRenderer(), 
                        background=(1.0, 1.0, 1.0))


def plot(ind, nworkers, frame_indices, frame_dir, filename_format, **kwargs):
    total = len(frame_indices)
    if total == 0:
        return
    start_idx = int(ind * total / nworkers)
    end_idx = int((ind + 1) * total / nworkers)
    end_idx = min(end_idx, total)
    for frame in frame_indices[start_idx:end_idx]:
        plot_frame_ovito(frame=frame, frame_dir=frame_dir, filename_format=filename_format, **kwargs)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--filename", type=str, required=True,
                        help="HDF5 file with 3D trajectory data")
    parser.add_argument("--frame_dir", type=str, default="frames")
    parser.add_argument("--Lx", type=float, default=12)
    parser.add_argument("--Ly", type=float, default=2)
    parser.add_argument("--Lz", type=float, default=2)
    parser.add_argument("--myosin_radius", type=float, default=0.05)
    parser.add_argument("--actin_length", type=float, default=1.0)
    parser.add_argument("--myosin_length", type=float, default=1.5)
    parser.add_argument("--sample_every", type=int, default=1)
    parser.add_argument("--start_frame", type=int, default=0)
    parser.add_argument("--end_frame", type=int, default=None)
    parser.add_argument("--myosin_display", choices=["all", "bonded", "cb_attached"], default="all")
    parser.add_argument("--actin_display", choices=["all", "cb"], default="cb")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    traj = h5py.File(args.filename, 'r')
    data = hdf5_to_dict(traj)
    
    nframes = data["/actin/center"].shape[0]
    print(f"Number of frames: {nframes}")
    
    start_frame = max(0, args.start_frame)
    end_frame = args.end_frame if args.end_frame is not None else nframes
    end_frame = min(end_frame, nframes)
    sample_every = max(1, args.sample_every)
    frame_indices = list(range(start_frame, end_frame, sample_every))

    cpu_workers = joblib.cpu_count()
    print(f"Using {cpu_workers} CPU workers for parallel processing.")
    
    if not os.path.exists(args.frame_dir):
        os.makedirs(args.frame_dir, exist_ok=True)
        
    n_digits = len(str(nframes))
    filename_format = "{}/frame_{:0" + str(n_digits) + "d}.png"

    Parallel(n_jobs=cpu_workers)(
        delayed(plot)(
            i, cpu_workers,
            frame_indices=frame_indices,
            frame_dir=args.frame_dir,
            filename_format=filename_format,
            data=data,
            myosin_radius=args.myosin_radius,
            actin_length=args.actin_length,
            myosin_length=args.myosin_length,
            Lx=args.Lx, Ly=args.Ly, Lz=args.Lz,
            myosin_display=args.myosin_display,
            actin_display=args.actin_display
        )
        for i in range(cpu_workers)
    )
    print("Done rendering.")
