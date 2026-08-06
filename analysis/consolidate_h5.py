#!/usr/bin/env python3
import os
import shutil
from pathlib import Path

# Set the base directory according to your SLURM script
BASE_DIR = Path("/scratch/users/xsli/sarcomere")

# Directories to search
INPUT_DIRS = [
    "data0330", 
    "data0331", 
    "data_resume0331", 
    "data0401"
]

# Destination directory
OUTPUT_DIR = BASE_DIR / "data_2dsweep"

def main():
    # Dictionary to store the path to the most recent version of each file
    # Key: filename (e.g., "A600_M1000_lt0_r01_L12x2x2.h5")
    # Value: Path object of the most recently modified version
    latest_files = {}

    print("Scanning directories for .h5 files...")
    
    for dir_name in INPUT_DIRS:
        dir_path = BASE_DIR / dir_name
        
        if not dir_path.exists():
            print(f"  [Warning] Directory not found, skipping: {dir_path}")
            continue
            
        # Find all .h5 files in this directory
        for h5_file in dir_path.glob("*.h5"):
            filename = h5_file.name
            file_mtime = h5_file.stat().st_mtime # Get modification time
            
            # If we haven't seen this file yet, add it
            if filename not in latest_files:
                latest_files[filename] = h5_file
            else:
                # If we have seen it, compare modification times
                existing_file = latest_files[filename]
                existing_mtime = existing_file.stat().st_mtime
                
                # If the current file is newer, replace the one in the dictionary
                if file_mtime > existing_mtime:
                    latest_files[filename] = h5_file

    total_unique = len(latest_files)
    if total_unique == 0:
        print("No .h5 files found in the specified directories. Exiting.")
        return

    print(f"Found {total_unique} unique .h5 files.")
    
    # Create the output directory if it doesn't exist
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Copying files to: {OUTPUT_DIR}")
    
    # Copy the files
    for count, (filename, file_path) in enumerate(latest_files.items(), 1):
        dest_path = OUTPUT_DIR / filename
        
        # shutil.copy2 preserves original file metadata (like exact modification times)
        shutil.copy2(file_path, dest_path)
        
        # Optional progress printout for large numbers of files
        if count % 20 == 0 or count == total_unique:
            print(f"  Copied {count}/{total_unique} files...")

    print("Done! All unique files successfully consolidated.")

if __name__ == "__main__":
    main()