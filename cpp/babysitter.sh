#!/bin/bash

TARGET_TOTAL_STEPS=${TARGET_TOTAL_STEPS:-400000}
PARTITIONS=${PARTITIONS:-normal,hns,owners}
CONSTRAINT=${CONSTRAINT:-CPU_GEN:GEN|CPU_GEN:MLN}
CPUS_PER_TASK=${CPUS_PER_TASK:-24}
STALLED_MINUTES=${STALLED_MINUTES:-20}
# Must exceed the babycron interval (30 min) so a still-queuing resume's marker survives one
# cron cycle and we don't pile up duplicate resume arrays for files awaiting their node slot.
RESUME_COOLDOWN_MINUTES=${RESUME_COOLDOWN_MINUTES:-45}
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
mkdir -p "${BACKUP_DIR}/.inflight"

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
    FNAME=$(basename "${FILE}" .h5)

    # GUARD 1 (preemption-requeue conflict): owners is PreemptMode=REQUEUE, so a preempted job
    # returns to the queue under the SAME job name (the dispatcher names each worker after its
    # output file). The worker does NOT take the resume flock, so a babysitter resume running
    # alongside a requeued worker would double-write the HDF5. If the worker job is still in the
    # queue, let Slurm finish it and do NOT resume.
    if squeue --me -h -n "${FNAME}" -o "%i" | grep -q .; then
        echo "Skipping (worker job still in queue, likely requeued after preemption): ${FNAME}"
        continue
    fi
    # GUARD 2 (resume already in flight): skip if a resume job for this file is STILL IN THE QUEUE
    # (PD or R). The marker stores that resume's array-task job id (written at submit time); we test
    # LIVENESS via squeue, not an age timer. Under cluster priority a resume can sit PD for hours, so
    # an age-based marker would expire and trigger duplicate arrays (the PD pileup). If the stored
    # job is gone (finished/died/cancelled), the marker is stale: drop it and allow a fresh resume.
    MARKER="${BACKUP_DIR}/.inflight/${FNAME}"
    if [ -f "${MARKER}" ]; then
        INFLIGHT_ID=$(cat "${MARKER}" 2>/dev/null)
        if [ -n "${INFLIGHT_ID}" ] && squeue --me -h -j "${INFLIGHT_ID}" -o "%i" 2>/dev/null | grep -q .; then
            echo "Skipping (resume ${INFLIGHT_ID} still in queue): ${FNAME}"
            continue
        fi
        rm -f "${MARKER}"   # stale marker: the recorded resume job no longer exists
    fi

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
    # Capture the array job id so we can record per-file in-flight markers (BACKUP_DIR is arg 4).
    ARRAY_ID=$(sbatch --parsable \
        --partition="${PARTITIONS}" \
        --constraint="${CONSTRAINT}" \
        --cpus-per-task="${CPUS_PER_TASK}" \
        --output="${LOG_DIR}/%j_%a.out" \
        --error="${LOG_DIR}/%j_%a.err" \
        --array=0-$((count - 1)) \
        sarcomere_resume_new.sbatch "${RESUME_LIST}" "${LOG_DIR}" "${TARGET_TOTAL_STEPS}" "${BACKUP_DIR}")
    echo "Submitted resume array ${ARRAY_ID}."
    # Record each file's resume task id so the next scan's GUARD 2 can check queue liveness.
    # Array task index == RESUME_LIST line order (the resume sbatch selects line SLURM_ARRAY_TASK_ID+1).
    if [ -n "${ARRAY_ID}" ]; then
        idx=0
        while read -r f; do
            echo "${ARRAY_ID}_${idx}" > "${BACKUP_DIR}/.inflight/$(basename "$f" .h5)"
            idx=$((idx + 1))
        done < "${RESUME_LIST}"
    fi
else
    echo "No stalled/incomplete jobs found. Everything is running or finished."
    rm -f "${RESUME_LIST}"
fi