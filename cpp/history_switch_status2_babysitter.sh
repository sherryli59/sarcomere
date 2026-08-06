#!/bin/bash
# One retry scan for the A2000 status-2 history-switch experiment.

set -euo pipefail

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
SWEEP_SBATCH="${CPP_DIR}/history_switch_status2.sbatch"
DATA_DIR="${HISTORY_OUTPUT_DIR:-/scratch/users/xsli/sarcomere/data_history_status2}"
TARGET_TOTAL_STEPS="${HISTORY_TOTAL_STEPS:-400000}"
N_ACTINS="${HISTORY_N_ACTINS:-2000}"
LOCK_FILE=/scratch/users/xsli/sarcomere/.history_switch_status2_babysitter.lock

mkdir -p "${DATA_DIR}"
exec 8>"${LOCK_FILE}"
if ! flock -n 8; then
    echo "another status-2 history watchdog scan is active; exiting"
    exit 0
fi

# Workers retain the generic name, so HDF5 inspection is deferred until the
# entire current array has left the queue.
if squeue --me -h -n cb_history_s2 -o '%i' 2>/dev/null | grep -q .; then
    echo "status-2 history workers remain in the queue; deferring retries"
    exit 0
fi

mapfile -t RETRY_IDS < <(python - "${DATA_DIR}" "${TARGET_TOTAL_STEPS}" "${N_ACTINS}" <<'PY'
import sys
from pathlib import Path

import h5py

directory = Path(sys.argv[1])
target = int(sys.argv[2])
n_actins = int(sys.argv[3])
protocols = ("medium_control", "high_control", "medium_high", "high_medium_high")

for protocol_idx, protocol in enumerate(protocols):
    for rep in range(3):
        task_id = protocol_idx * 3 + rep
        path = directory / f"A{n_actins}_M1000_history_{protocol}_r{rep:02d}_L12x2x2.h5"
        if not path.exists() or path.stat().st_size == 0:
            print(task_id)
            continue
        try:
            with h5py.File(path, "r") as handle:
                if "/resume/current_step" in handle:
                    step = int(handle["/resume/current_step"][...].reshape(-1)[-1])
                elif "/resume_prev/current_step" in handle:
                    step = int(handle["/resume_prev/current_step"][...].reshape(-1)[-1])
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
    echo "no missing or incomplete status-2 history tasks need submission"
    exit 0
fi

ARRAY_SPEC=$(IFS=,; echo "${RETRY_IDS[*]}")
echo "submitting incomplete history task IDs ${ARRAY_SPEC} without a concurrency throttle"
if [ "${HISTORY_WATCHDOG_DRY_RUN:-0}" = 1 ]; then
    echo "dry run: sbatch --array=${ARRAY_SPEC} ${SWEEP_SBATCH}"
    exit 0
fi
sbatch --array="${ARRAY_SPEC}" "${SWEEP_SBATCH}"
