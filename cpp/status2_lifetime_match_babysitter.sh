#!/bin/bash
# One retry scan for status2_lifetime_match.sbatch.  Because workers retain the
# generic Slurm job name, wait until the current array has completely left the
# queue before inspecting HDF5 files or resubmitting exact incomplete task IDs.

set -euo pipefail

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
SWEEP_SBATCH="${CPP_DIR}/status2_lifetime_match.sbatch"
DATA_DIR="${CONTROL_OUTPUT_DIR:-/scratch/users/xsli/sarcomere/data_status2_lifetime_match}"
TARGET_TOTAL_STEPS="${CONTROL_NSTEPS:-150000}"
N_ACTINS="${CONTROL_N_ACTINS:-2000}"
LOCK_FILE=/scratch/users/xsli/sarcomere/.status2_lifetime_match_babysitter.lock

mkdir -p "${DATA_DIR}"
exec 8>"${LOCK_FILE}"
if ! flock -n 8; then
    echo "another status-2 lifetime-match watchdog scan is active; exiting"
    exit 0
fi

if squeue --me -h -n s2_lifetime_match -o '%i' 2>/dev/null | grep -q .; then
    echo "status-2 lifetime-match workers remain in the queue; deferring retries"
    exit 0
fi

mapfile -t RETRY_IDS < <(python - "${DATA_DIR}" "${TARGET_TOTAL_STEPS}" "${N_ACTINS}" <<'PY'
import sys
from pathlib import Path

import h5py

directory = Path(sys.argv[1])
target = int(sys.argv[2])
n_actins = int(sys.argv[3])
tau_values = ("0.001", "0.002", "0.003", "0.004", "0.006")
kon_values = (500, 1000, 2000)
arms = [("catch_lt0p003162_kon1000",)]
arms += [(f"consts2_tau{tau.replace('.', 'p')}_kon{kon}",)
         for tau in tau_values for kon in kon_values]

for arm_idx, (tag,) in enumerate(arms):
    for rep in range(2):
        task_id = arm_idx * 2 + rep
        path = directory / f"A{n_actins}_M1000_{tag}_r{rep:02d}_L12x2x2.h5"
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
    echo "no missing or incomplete status-2 lifetime-match tasks need submission"
    exit 0
fi

ARRAY_SPEC=$(IFS=,; echo "${RETRY_IDS[*]}")
echo "submitting incomplete task IDs ${ARRAY_SPEC} (maximum 8 concurrent)"
if [ "${S2_MATCH_WATCHDOG_DRY_RUN:-0}" = 1 ]; then
    echo "dry run: sbatch --array=${ARRAY_SPEC}%8 ${SWEEP_SBATCH}"
    exit 0
fi
sbatch --array="${ARRAY_SPEC}%8" "${SWEEP_SBATCH}"
