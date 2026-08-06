#!/bin/bash
# One scan of the wider-box output dir: for each stalled, incomplete .h5 with no
# live job (worker OR resume) of the same name, submit a dedicated resume.
# Guard is by JOB NAME: both the main array tasks and the resume jobs rename
# themselves to their filename, so `squeue -n FILENAME` catches either -> we
# never double-run a file. Run from cpp/ with the sarcomere env active (h5py).
set -uo pipefail

DATA_DIR=/scratch/users/xsli/sarcomere/data_widebox
RESUME_SBATCH=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp/wider_box_resume.sbatch
TARGET=400000
SAVE_EVERY=200
STALLED_MINUTES=${STALLED_MINUTES:-30}

echo "wb-babysitter scan $(date)"
shopt -s nullglob
for FILE in "${DATA_DIR}"/*.h5; do
    case "$FILE" in *.lock) continue;; esac
    FNAME=$(basename "$FILE" .h5)

    # GUARD: a job named FNAME (live worker or in-flight resume) exists -> skip
    if squeue --me -h -n "$FNAME" -o "%i" 2>/dev/null | grep -q .; then
        echo "  [$FNAME] job in queue -> skip"; continue
    fi
    # skip if recently written (a slow-but-alive run whose job already left queue)
    if [ -n "$(find "$FILE" -maxdepth 0 -mmin -"${STALLED_MINUTES}" 2>/dev/null)" ]; then
        echo "  [$FNAME] modified <${STALLED_MINUTES}min ago -> skip"; continue
    fi
    # completeness: frames * save_every >= target
    STEP=$(python - "$FILE" "$SAVE_EVERY" <<'PY'
import sys, h5py
fn, se = sys.argv[1], int(sys.argv[2]); s = -1
try:
    with h5py.File(fn, "r") as f:
        for c in ("/actin/center", "/myosin/center"):
            if c in f and len(f[c].shape) >= 1:
                s = int(f[c].shape[0]) * se; break
except Exception:
    pass
print(s)
PY
)
    STEP=$(echo "$STEP" | tr -d '[:space:]')
    if [[ "$STEP" =~ ^[0-9]+$ ]] && [ "$STEP" -ge "$TARGET" ]; then
        echo "  [$FNAME] complete ($STEP) -> skip"; continue
    fi
    echo "  [$FNAME] stalled at step ${STEP} -> submitting resume"
    sbatch "$RESUME_SBATCH" "$FILE"
done
