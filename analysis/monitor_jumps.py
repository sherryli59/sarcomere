import os
import glob
import time
import h5py
import numpy as np

# ==========================================
# CONFIGURATION & THRESHOLDS
# ==========================================
DATA_DIR = "/scratch/users/xsli/sarcomere/data0416" 
CHECK_INTERVAL = 3

# Simulation Box for Minimum Image Convention
BOX = np.array([12.0, 2.0, 2.0])

# SKIP SETTING
SKIP_FRAMES = 100  # Analysis will start at frame 101

# Thresholds for AVERAGE jump
THRESH_AVG_POS = 0.3   # Max AVERAGE distance (microns) actins should move
THRESH_RATIO   = 0.1   # Max change (10%) in the global strong catch-bond ratio
# ==========================================

def minimal_image(dx, box):
    """Applies minimum image convention to distance vectors."""
    return dx - box * np.round(dx / box)

def check_resume_discontinuity(filepath, last_checked_frame):
    warnings = []
    
    try:
        # Open in SWMR mode to read while the simulation is writing
        with h5py.File(filepath, 'r', swmr=True) as f:
            if 'actin/center' not in f:
                return last_checked_frame, []
                
            a_cent = f['actin/center']
            a_cb   = f['actin/cb_status'] if 'actin/cb_status' in f else None
            
            # Refresh to see newly appended data from the C++ process
            a_cent.refresh()
            if a_cb: a_cb.refresh()

            n_frames = a_cent.shape[0]
            
            # If the file hasn't reached the skip threshold, do nothing yet
            if n_frames <= SKIP_FRAMES + 1 or last_checked_frame >= n_frames - 1:
                return max(last_checked_frame, n_frames), []

            # Ensure we start checking from frame 101 or the last verified frame
            start = max(SKIP_FRAMES + 1, last_checked_frame)
            n_actins = a_cent.shape[1]

            for i in range(start, n_frames):
                # 1. AVERAGE POSITION INTEGRITY CHECK
                # Minimal image used for periodic boundary wrapping
                dx = minimal_image(a_cent[i] - a_cent[i-1], BOX)
                avg_dist = np.mean(np.linalg.norm(dx, axis=-1))
                
                if avg_dist > THRESH_AVG_POS:
                    warnings.append(f"  [FRAME {i}] DISCONTINUITY: Avg position jump of {avg_dist:.4f} μm detected.")

                # 2. STATE INTEGRITY CHECK (Catch-Bond Ratio)
                # Monitors the ratio of status=2 (strong bonds)
                if a_cb is not None:
                    ratio_curr = np.sum(a_cb[i] == 2) / n_actins
                    ratio_prev = np.sum(a_cb[i-1] == 2) / n_actins
                    ratio_diff = abs(ratio_curr - ratio_prev)
                    
                    if ratio_diff > THRESH_RATIO:
                        warnings.append(f"  [FRAME {i}] STATE SHIFT: CB Ratio jumped from {ratio_prev:.3f} to {ratio_curr:.3f}")

            return n_frames, warnings
            
    except (OSError, KeyError):
        return last_checked_frame, []
    except Exception as e:
        return last_checked_frame, [f"  -> Error: {e}"]

def main():
    print(f"Monitoring resume integrity in: {DATA_DIR}")
    print(f"Skipping first {SKIP_FRAMES} frames of each file.")
    print(f"Tracking: Avg Position (Threshold={THRESH_AVG_POS}) and CB Ratio (Threshold={THRESH_RATIO})\n")
    
    file_states = {} # filepath: last_frame_verified
    
    while True:
        h5_files = glob.glob(os.path.join(DATA_DIR, "*.h5"))
        print(f"[{time.strftime('%H:%M:%S')}] Starting scan of {len(h5_files)} files...")
        
        for filepath in h5_files:
            fname = os.path.basename(filepath)
            print(f"  -> Checking: {fname}")
            
            last_frame = file_states.get(filepath, 0)
            new_frame_count, issues = check_resume_discontinuity(filepath, last_frame)
            
            if issues:
                print(f"    *** ALARM in {fname} ***")
                for issue in issues:
                    print(f"    {issue}")
            
            file_states[filepath] = new_frame_count
            
        print(f"[{time.strftime('%H:%M:%S')}] Scan complete. Sleeping for {CHECK_INTERVAL}s...\n")
        time.sleep(CHECK_INTERVAL)

if __name__ == "__main__":
    main()