#!/bin/bash
# Scan the 12x4x4 lambda sweep once and resume stale, incomplete trajectories
# that do not already have a live main or resume job.

set -uo pipefail

DATA_DIR=/scratch/users/xsli/sarcomere/data_wb4x4_lambda
RESUME_SBATCH=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp/wb4_lambda_resume.sbatch
TARGET_TOTAL_STEPS=300000
SAVE_EVERY=200
STALLED_MINUTES="${STALLED_MINUTES:-30}"

echo "wb4-lambda babysitter scan $(date)"
shopt -s nullglob
for FILE in "${DATA_DIR}"/A*_M*_lt*_r*_L12x4x4.h5; do
    FNAME=$(basename "${FILE}" .h5)

    # Both original workers and resume workers rename themselves to FNAME.
    if squeue --me -h -n "${FNAME}" -o '%i' 2>/dev/null | grep -q .; then
        echo "  [${FNAME}] job in queue -> skip"
        continue
    fi
    if [ -n "$(find "${FILE}" -maxdepth 0 -mmin -"${STALLED_MINUTES}" 2>/dev/null)" ]; then
        echo "  [${FNAME}] modified <${STALLED_MINUTES} min ago -> skip"
        continue
    fi

    CURRENT_STEP=$(python - "${FILE}" "${SAVE_EVERY}" <<'PY'
import sys
import h5py

path, save_every = sys.argv[1], int(sys.argv[2])
try:
    with h5py.File(path, "r") as handle:
        if "/resume/current_step" in handle:
            step = int(handle["/resume/current_step"][...].reshape(-1)[-1])
        elif "/state/current_step" in handle:
            step = int(handle["/state/current_step"][...].reshape(-1)[-1])
        else:
            step = (int(handle["/actin/center"].shape[0]) - 1) * save_every
    print(step)
except Exception:
    print(-1)
PY
)
    CURRENT_STEP=$(echo "${CURRENT_STEP}" | tr -d '[:space:]')
    if [ "${CURRENT_STEP}" = "-1" ]; then
        echo "  [${FNAME}] unreadable or locked -> retry next scan"
        continue
    fi
    if [[ "${CURRENT_STEP}" =~ ^[0-9]+$ ]] && [ "${CURRENT_STEP}" -ge "${TARGET_TOTAL_STEPS}" ]; then
        echo "  [${FNAME}] complete (${CURRENT_STEP}) -> skip"
        continue
    fi

    echo "  [${FNAME}] stalled at ${CURRENT_STEP} -> submitting resume"
    sbatch --job-name="${FNAME}" "${RESUME_SBATCH}" "${FILE}"
done
