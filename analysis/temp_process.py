#!/usr/bin/env python3
import os
import shutil
import time
import h5py
from pathlib import Path

# Set the base directory
BASE_DIR = Path("/scratch/users/xsli/sarcomere")

# The 3 directories to scan for files that might need resuming
INPUT_DIRS = [
    BASE_DIR / "data_2dsweep",
    BASE_DIR / "data_2dsweep_new",
    BASE_DIR / "data_2dsweep_3",
    BASE_DIR / "need_resume"
]

# Destination directory for the files that need to be resumed
OUTPUT_DIR = BASE_DIR / "need_resume2"

def main():
    # Calculate the timestamp for 8 minutes ago
    cutoff_time = time.time() - (8 * 60)

    # List to keep track of files that need to be moved
    files_to_move = []
    
    print("Scanning directories for incomplete .h5 files...\n")
    
    for dir_name in INPUT_DIRS:
        dir_path = BASE_DIR / dir_name
        
        if not dir_path.exists():
            print(f"  [Warning] Directory not found, skipping: {dir_path}")
            continue
            
        print(f"Checking directory: {dir_name}")
        
        # Find all .h5 files in this directory
        for h5_file in dir_path.glob("*.h5"):
            # 1. Check modification time first (it's faster than opening the H5)
            file_mtime = h5_file.stat().st_mtime
            
            # If modified within the last 5 minutes, skip it
            if file_mtime >= cutoff_time:
                continue
                
            # 2. Check the number of frames
            n_frames = 0
            try:
                with h5py.File(h5_file, 'r') as fh:
                    if "/actin/centers" in fh:
                        n_frames = fh["/actin/centers"].shape[0]
                    elif "/actin/center" in fh:
                        n_frames = fh["/actin/center"].shape[0]
                    else:
                        print(f"  [Warning] No actin centers found in {h5_file.name}")
                        continue
            except Exception as e:
                print(f"  [Error] Could not read {h5_file.name}: {e}")
                # We'll skip corrupted/unreadable files to avoid crashing the resume script later
                continue
                
            # If less than 2994 frames, mark for moving
            if n_frames < 2994:
                files_to_move.append((h5_file, n_frames))

    total_to_move = len(files_to_move)
    if total_to_move == 0:
        print("\nNo incomplete files found. Exiting.")
        return

    print(f"\nFound {total_to_move} files with < 2994 frames that are safe to move.")
    
    # Create the output directory if it doesn't exist
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Moving files to: {OUTPUT_DIR}")
    
    # Move the files
    for count, (file_path, frames) in enumerate(files_to_move, 1):
        dest_path = OUTPUT_DIR / file_path.name
        
        # # Use shutil.move to physically move the file
        shutil.move(str(file_path), str(dest_path))
        
        if count % 20 == 0 or count == total_to_move:
            print(f"  Moved {count}/{total_to_move} files...")

    print("\nDone! All incomplete files successfully moved to need_resume.")

if __name__ == "__main__":
    main()