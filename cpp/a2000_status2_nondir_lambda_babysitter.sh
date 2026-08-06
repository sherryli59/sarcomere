#!/bin/bash
# One scan for the A2000/M1000 nondirectional, status-2-only lambda sweep.

set -euo pipefail

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
SWEEP_SBATCH="${CPP_DIR}/a2000_status2_nondir_lambda.sbatch"
DATA_DIR="${A2K_ND_OUTPUT_DIR:-/scratch/users/xsli/sarcomere/data_a2000_status2_nondir_lambda}"
TARGET_TOTAL_STEPS="${A2K_ND_NSTEPS:-300000}"
STALLED_MINUTES="${A2K_ND_STALLED_MINUTES:-30}"
LOCK_FILE=/scratch/users/xsli/sarcomere/.a2000_status2_nondir_lambda_babysitter.lock

mkdir -p "${DATA_DIR}"
exec 8>"${LOCK_FILE}"
if ! flock -n 8; then
    echo "another A2000 nondirectional watchdog scan is active; exiting"
    exit 0
fi

LIVE_JOB_NAMES=$(squeue --me -h -o '%j' 2>/dev/null || true)
if grep -Fxq 'a2k_s2_nd_lambda' <<< "${LIVE_JOB_NAMES}"; then
    echo "a nondirectional array has pending/unrenamed elements; deferring retries"
    exit 0
fi

mapfile -t RETRY_IDS < <(python - "${DATA_DIR}" "${TARGET_TOTAL_STEPS}" "${STALLED_MINUTES}" "${LIVE_JOB_NAMES}" <<'PY'
import sys
import time
from pathlib import Path

import h5py

directory = Path(sys.argv[1])
target = int(sys.argv[2])
stalled_seconds = int(sys.argv[3]) * 60
live_names = set(sys.argv[4].splitlines())
labels = (
    "0", "1e-5.5", "1e-5", "1e-4.5", "1e-4", "1e-3.5", "1e-3",
    "1e-2.5", "1e-2", "1e-1.5", "1e-1", "1e-0.5", "1e0",
)
now = time.time()

for rep in range(3):
    for lambda_idx, label in enumerate(labels):
        task_id = rep * len(labels) + lambda_idx
        stem = f"A2000_M1000_s2only_lt{label}_nodir_r{rep:02d}_L12x2x2"
        path = directory / f"{stem}.h5"
        if stem in live_names:
            continue
        if not path.exists() or path.stat().st_size == 0:
            print(task_id)
            continue
        if now - path.stat().st_mtime < stalled_seconds:
            continue
        try:
            with h5py.File(path, "r") as handle:
                if "/resume/current_step" in handle:
                    step = int(handle["/resume/current_step"][...].reshape(-1)[-1])
                elif "/state/current_step" in handle:
                    step = int(handle["/state/current_step"][...].reshape(-1)[-1])
                else:
                    step = -1
        except Exception:
            step = -1
        if step < target:
            print(task_id)
PY
)

if [ "${#RETRY_IDS[@]}" -eq 0 ]; then
    echo "no missing or stalled nondirectional tasks need submission"
    exit 0
fi

ARRAY_SPEC=$(IFS=,; echo "${RETRY_IDS[*]}")
echo "submitting task IDs ${ARRAY_SPEC} without an array concurrency throttle"
if [ "${A2K_ND_WATCHDOG_DRY_RUN:-0}" = 1 ]; then
    echo "dry run: sbatch --array=${ARRAY_SPEC} ${SWEEP_SBATCH}"
    exit 0
fi
sbatch --array="${ARRAY_SPEC}" "${SWEEP_SBATCH}"
