import h5py
import numpy as np
from pathlib import Path
from functools import lru_cache

@lru_cache(maxsize=1)
def get_valid_paths_and_skipped():
    """Return a list of valid (label, Path) tuples and a list of skipped files."""
    base = Path('/scratch/users/xsli/sarcomere/data0313')  # Adjust as needed
    files = sorted(base.glob('**/*.h5'))
    
    valid_paths = []
    skipped_files = []
    
    for h5_path in files:
        label = h5_path.stem
        # Skip files that have already been downsampled to avoid infinite loops
        if "_stride" in label:
            continue
            
        try:
            # Attempt to open the file to verify it is not corrupted
            with h5py.File(h5_path, 'r') as fh:
                # Quick integrity check
                _ = fh.keys()
            valid_paths.append((label, h5_path))
            
        except (OSError, KeyError, RuntimeError, ValueError) as e:
            skipped_files.append((label, str(h5_path), str(e)))
            print(f"Skipping broken/unreadable file: {label} ({h5_path})")
            print(f"  reason: {e}")
            continue
            
    return valid_paths, skipped_files


def downsample_h5_file(src_path: Path, stride: int = 10):
    """
    Creates a copy of the HDF5 file, saving only every `stride` frames 
    for time-series datasets.
    """
    # Create the output filename (e.g., filename_stride10.h5)
    out_path = src_path.with_name(f"{src_path.stem}_stride{stride}{src_path.suffix}")
    
    with h5py.File(src_path, 'r') as src, h5py.File(out_path, 'w') as dst:
        
        # 1. Figure out the total number of frames to safely identify time-series datasets
        n_frames = None
        if "/actin/bonds" in src:
            n_frames = src["/actin/bonds"].shape[0]
        elif "/actin/direction" in src:
            n_frames = src["/actin/direction"].shape[0]
            
        if n_frames is None:
            print(f"  [Warning] Could not determine n_frames for {src_path.name}. Skipping.")
            return False

        # 2. Define the visitor function to process all groups and datasets
        def visit_node(name, node):
            if isinstance(node, h5py.Dataset):
                # If the dataset's first dimension matches n_frames, it's a trajectory array.
                # We slice it to keep only every Nth frame.
                if node.shape and node.shape[0] == n_frames:
                    data = node[::stride]
                else:
                    # If it's not frame-dependent (e.g., /resume snapshots, /catch_bond event logs),
                    # copy it exactly as-is so the file doesn't break.
                    data = node[...]
                    
                # Create the dataset in the new file, preserving compression settings
                dst.create_dataset(
                    name, 
                    data=data, 
                    compression=node.compression, 
                    compression_opts=node.compression_opts
                )
            elif isinstance(node, h5py.Group):
                dst.require_group(name)

        # 3. Copy structure and datasets
        src.visititems(visit_node)
        
        # 4. Copy attributes (metadata) for groups and datasets
        for k, v in src.attrs.items():
            dst.attrs[k] = v
            
        def copy_attrs(name, node):
            if name in dst:
                for k, v in node.attrs.items():
                    dst[name].attrs[k] = v
                    
        src.visititems(copy_attrs)
        
    print(f"  -> Successfully saved downsampled file: {out_path.name}")
    return True


# ==========================================
# Main Execution
# ==========================================
stride_value = 10
paths, skipped = get_valid_paths_and_skipped()

print("-" * 60)
print(f"Found {len(paths)} valid files to process.")
print(f"Skipped {len(skipped)} corrupted/unreadable files.")
print("-" * 60)

for label, h5_path in paths:
    print(f"Processing: {label}")
    downsample_h5_file(h5_path, stride=stride_value)

print("\nAll processing complete!")