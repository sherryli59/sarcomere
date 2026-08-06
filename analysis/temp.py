#!/usr/bin/env python3
import time
import shutil
from pathlib import Path

# Set the base directory 
BASE_DIR = Path("/scratch/users/xsli/sarcomere")

# Source and Destination directories
SOURCE_DIR = BASE_DIR / "data_2dsweep"
DEST_DIR = BASE_DIR / "data_2dsweep_3"

def main():
    # Sanity check: Ensure source exists
    if not SOURCE_DIR.exists():
        print(f"Source directory {SOURCE_DIR} does not exist. Exiting.")
        return

    # Create destination if it doesn't exist
    if not DEST_DIR.exists():
        print(f"Creating destination directory: {DEST_DIR}")
        DEST_DIR.mkdir(parents=True, exist_ok=True)

    # Calculate the cutoff (5 minutes ago)
    current_time = time.time()
    cutoff_time = current_time - (5 * 60)
    
    moved_count = 0

    print(f"Scanning {SOURCE_DIR.name} for stagnant files (not modified in >5 mins)...")

    # Iterate through all files in the source directory
    # Using rglob("*") if you want subdirectories, or glob("*") for just top-level
    for file_path in SOURCE_DIR.glob("*.h5"):
        if not file_path.is_file():
            continue
            
        file_mtime = file_path.stat().st_mtime
        
        # LOGIC: If modification time is OLDER than the cutoff
        if file_mtime < cutoff_time:
            target_path = DEST_DIR / file_path.name
            
            try:
                # shutil.move is safer than Path.rename across different filesystems/mounts
                shutil.move(str(file_path), str(target_path))
                moved_count += 1
                print(f"  [->] Moved: {file_path.name} to {DEST_DIR.name}")
            except Exception as e:
                print(f"  [!] Error moving {file_path.name}: {e}")

    print(f"\nDone! Moved {moved_count} stagnant files to {DEST_DIR.name}.")

if __name__ == "__main__":
    main()