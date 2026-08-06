#!/usr/bin/env python3
import os
import shutil
import time
import h5py
import argparse
from pathlib import Path
from collections import defaultdict

# Set the base directory
BASE_DIR = Path("/scratch/users/xsli/sarcomere")

# The directories to scan
INPUT_DIRS = [
    BASE_DIR / "2d_sweep_need_resume",
]

# Destination directories
DIR_FINISHED = BASE_DIR / "2d_sweep_finished"
DIR_RESUME = BASE_DIR / "2d_sweep_need_resume2"

# Create directories if they don't exist
DIR_FINISHED.mkdir(parents=True, exist_ok=True)
DIR_RESUME.mkdir(parents=True, exist_ok=True)

def main():
    parser = argparse.ArgumentParser(description="Identify and migrate Sarcomere H5 files.")
    parser.add_argument("--dry-run", action="store_true", 
                        help="Identify files and print actions without actually moving anything.")
    args = parser.parse_args()

    if args.dry_run:
        print("\n" + "="*50)
        print("🚨 DRY RUN MODE ENABLED 🚨")
        print("Files will be identified but NOT moved.")
        print("="*50 + "\n")

    # Calculate the timestamp for 8 minutes ago (to protect actively running jobs)
    cutoff_time = time.time() - (8 * 60)
    
    # Dictionary to group files by their filename
    file_registry = defaultdict(list)

    print("Scanning directories and reading frame counts...\n")
    
    for dir_path in INPUT_DIRS:
        if not dir_path.exists():
            print(f"  [Warning] Directory not found, skipping: {dir_path.name}")
            continue
            
        print(f"Checking directory: {dir_path.name}")
        
        for h5_file in dir_path.glob("*.h5"):
            mtime = h5_file.stat().st_mtime
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
                n_frames = -1 # Mark as corrupted
                
            file_registry[h5_file.name].append({
                "path": h5_file, 
                "frames": n_frames, 
                "mtime": mtime
            })

    if not file_registry:
        print("\nNo .h5 files found. Exiting.")
        return

    # Create the output directories if this isn't a dry run
    if not args.dry_run:
        DIR_FINISHED.mkdir(parents=True, exist_ok=True)
        DIR_RESUME.mkdir(parents=True, exist_ok=True)
    
    print(f"\nProcessing {len(file_registry)} unique filenames...")
    
    stats = {"finished": 0, "resumed": 0, "skipped_active": 0, "ignored_dupes": 0}
    action_word = "Would move" if args.dry_run else "Moved"

    for fname, records in file_registry.items():
        # Sort the list of files by frame count, descending (highest frames first)
        records.sort(key=lambda x: x["frames"], reverse=True)
        
        best_file = records[0]
        duplicates = records[1:]
        
        # 1. Handle the "Best" File
        if best_file["frames"] >= 2994:
            # It's finished
            dest_path = DIR_FINISHED / fname
            if best_file["path"] != dest_path:
                print(f"  [{action_word}] {fname} -> finished ({best_file['frames']} frames)")
                if not args.dry_run:
                    shutil.move(str(best_file["path"]), str(dest_path))
            stats["finished"] += 1
            
        elif best_file["frames"] >= 0:
            # It's incomplete. Protect it if it's actively running.
            if best_file["mtime"] >= cutoff_time:
                print(f"  [Skipped] {fname} appears to be actively running (< 8 mins old).")
                stats["skipped_active"] += 1
            else:
                dest_path = DIR_RESUME / fname
                if best_file["path"] != dest_path:
                    print(f"  [{action_word}] {fname} -> need_resume ({best_file['frames']} frames)")
                    if not args.dry_run:
                        shutil.move(str(best_file["path"]), str(dest_path))
                stats["resumed"] += 1

        # 2. Handle the Duplicates (No longer deleting)
        for dup in duplicates:
            print(f"  [Duplicate Ignored] Found shorter duplicate of {fname} at: {dup['path']} ({dup['frames']} frames)")
            stats["ignored_dupes"] += 1

    # Final Report
    print("\n" + "="*45)
    print(f"MIGRATION REPORT {'(DRY RUN)' if args.dry_run else '(COMPLETED)'}")
    print("="*45)
    print(f"Files {'to move' if args.dry_run else 'moved'} to finished (>2994): {stats['finished']}")
    print(f"Files {'to move' if args.dry_run else 'moved'} to need_resume:      {stats['resumed']}")
    print(f"Skipped (actively running):           {stats['skipped_active']}")
    print(f"Ignored shorter duplicates:           {stats['ignored_dupes']}")
    print("="*45)

if __name__ == "__main__":
    main()
