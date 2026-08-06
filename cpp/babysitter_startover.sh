#!/bin/bash

TARGET_TOTAL_STEPS=${TARGET_TOTAL_STEPS:-400000}
PARTITIONS=${PARTITIONS:-normal,hns,owners}
CONSTRAINT=${CONSTRAINT:-CPU_GEN:GEN|CPU_GEN:MLN}
CPUS_PER_TASK=${CPUS_PER_TASK:-24}
STALLED_MINUTES=${STALLED_MINUTES:-20}
SAVE_EVERY=100
SCRATCH_DIR="/scratch/users/xsli/sarcomere"

DATA_DIR=${1:-"${SCRATCH_DIR}/data_2d_sweep_additional"}
DATA_DIR=${DATA_DIR%/}

DIR_NAME=$(basename "${DATA_DIR}")
if [[ "$DIR_NAME" == data_* ]]; then
    LOG_DIR_NAME="logs_${DIR_NAME#data_}"
    BACKUP_DIR_NAME="backups_${DIR_NAME#data_}"
else
    LOG_DIR_NAME="logs_${DIR_NAME}"
    BACKUP_DIR_NAME="backups_${DIR_NAME}"
fi

LOG_DIR="${SCRATCH_DIR}/${LOG_DIR_NAME}"
BACKUP_DIR="${SCRATCH_DIR}/${BACKUP_DIR_NAME}"
mkdir -p "${LOG_DIR}"
mkdir -p "${BACKUP_DIR}"

TIMESTAMP=$(date +%s)
RESUME_LIST="${SCRATCH_DIR}/resume_list_${TIMESTAMP}.txt"
> "${RESUME_LIST}"

echo "Babysitter running at $(date)."
echo "Data Directory: ${DATA_DIR}"
echo "Log Directory: ${LOG_DIR}"
echo "Backup Directory: ${BACKUP_DIR}"
echo "Target Steps: ${TARGET_TOTAL_STEPS}"
echo "Stalled Threshold Minutes: ${STALLED_MINUTES}"

STALLED_FILES=$(find "${DATA_DIR}" -maxdepth 1 -type f -name "*.h5" -mmin +"${STALLED_MINUTES}")

count=0
for FILE in $STALLED_FILES; do
    CURRENT_STEP_RAW=$(python3 - "${FILE}" "${SAVE_EVERY}" << 'PY'
import sys, h5py
try:
    with h5py.File(sys.argv[1], "r") as f:
        for cand in ("/resume/current_step", "resume/current_step", "/state/current_step", "state/current_step"):
            if cand in f:
                d = f[cand][...]
                if getattr(d, "size", 0) > 0:
                    print(int(d.reshape(-1)[-1]))
                    sys.exit(0)
        for cand in ("/actin/centers", "actin/centers", "/actin/center", "actin/center", "/myosin/center"):
            if cand in f:
                d = f[cand]
                if hasattr(d, "shape") and len(d.shape) >= 1:
                    print(int(d.shape[0]) * int(sys.argv[2]))
                    sys.exit(0)
    print(0)
except Exception:
    print(-1)
PY
)
    CURRENT_STEP=$(echo "$CURRENT_STEP_RAW" | tr -d '[:space:]')

    if [[ "${CURRENT_STEP}" =~ ^[0-9]+$ ]] && [ "${CURRENT_STEP}" -ge 0 ] && [ "${CURRENT_STEP}" -lt "${TARGET_TOTAL_STEPS}" ]; then
        # Just append to the list; the sbatch script will handle the copying
        echo "${FILE}" >> "${RESUME_LIST}"
        ((count++))
        echo "Found stalled/incomplete file: $(basename "$FILE") (Steps: $CURRENT_STEP)"
    fi
done

if [ "$count" -gt 0 ]; then
    echo "Submitting resume array for ${count} jobs on partitions '${PARTITIONS}'..."
    # Passing the BACKUP_DIR as the 4th argument
    sbatch \
        --partition="${PARTITIONS}" \
        --constraint="${CONSTRAINT}" \
        --cpus-per-task="${CPUS_PER_TASK}" \
        --output="${LOG_DIR}/%j_%a.out" \
        --error="${LOG_DIR}/%j_%a.err" \
        --array=0-$((count - 1)) \
        sarcomere_startover.sbatch "${RESUME_LIST}" "${LOG_DIR}" "${TARGET_TOTAL_STEPS}" "${BACKUP_DIR}"
else
    echo "No stalled/incomplete jobs found. Everything is running or finished."
    rm -f "${RESUME_LIST}"
fi