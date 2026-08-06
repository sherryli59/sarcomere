#!/bin/bash
# One-shot dependent watcher. It is chained after each production attempt and
# submits a fresh --no-requeue job only for explicitly recoverable scheduler
# states. Scientific/configuration failures and user cancellations stop.
#SBATCH --partition=owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=00:10:00
#SBATCH --job-name=sarc_exact_watch
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/prod_exact_setup_20260804/logs/%x_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/prod_exact_setup_20260804/logs/%x_%j.err

set -euo pipefail

ROOT=/scratch/users/xsli/sarcomere/prod_exact_setup_20260804
RUN_SCRIPT=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp/production_exact_minimal_mechanics_A1280_M640_20260804.sbatch
WATCH_SCRIPT=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp/production_exact_minimal_mechanics_A1280_M640_babysitter_20260804.sbatch
OUTPUT=${ROOT}/data/A1280_M640_exact_seed20260804.h5
TARGET_STEP=100000
MAX_RETRIES=4
PREVIOUS_JOB=${WATCH_JOB_ID:?WATCH_JOB_ID must name the production attempt to inspect}
RETRY_FILE=${ROOT}/monitor/retry_count
STATUS_FILE=${ROOT}/monitor/babysitter_status.txt
mkdir -p "${ROOT}/monitor" "${ROOT}/logs"

exec 8>"${ROOT}/monitor/babysitter.lock"
flock -n 8 || { echo "another babysitter owns the monitor lock"; exit 0; }

# Accounting can lag briefly behind an afterany dependency becoming eligible.
state=""
for _ in $(seq 1 24); do
    state=$(sacct -X -j "${PREVIOUS_JOB}" -n -P -o State 2>/dev/null |
        head -n 1 | cut -d'|' -f1 | sed 's/+.*//')
    [ -n "${state}" ] && break
    sleep 5
done
if [ -z "${state}" ]; then
    echo "could not resolve terminal state for ${PREVIOUS_JOB}" >&2
    exit 1
fi

checkpoint_step=$(HDF5_USE_FILE_LOCKING=FALSE python - "${OUTPUT}" <<'PY'
import sys
from pathlib import Path
import h5py

path = Path(sys.argv[1])
if not path.exists() or path.stat().st_size == 0:
    print(-1)
    raise SystemExit
try:
    with h5py.File(path, "r", locking=False) as h:
        steps = []
        for group in ("/resume", "/resume_prev"):
            key = f"{group}/current_step"
            if key in h:
                steps.append(int(h[key][...].reshape(-1)[-1]))
        print(max(steps) if steps else -1)
except Exception:
    print(-1)
PY
)

timestamp=$(date -u +%Y-%m-%dT%H:%M:%SZ)
printf '%s previous_job=%s state=%s checkpoint=%s/%s\n' \
    "${timestamp}" "${PREVIOUS_JOB}" "${state}" \
    "${checkpoint_step}" "${TARGET_STEP}" | tee "${STATUS_FILE}"

if [ "${checkpoint_step}" -ge "${TARGET_STEP}" ]; then
    printf '%s COMPLETE job=%s state=%s checkpoint=%s/%s\n' \
        "${timestamp}" "${PREVIOUS_JOB}" "${state}" "${checkpoint_step}" \
        "${TARGET_STEP}" | tee "${ROOT}/monitor/gate_complete.txt"
    exit 0
fi

if [ "${state}" = COMPLETED ]; then
    printf '%s STOP completed_job_missing_target job=%s checkpoint=%s/%s\n' \
        "${timestamp}" "${PREVIOUS_JOB}" "${checkpoint_step}" \
        "${TARGET_STEP}" | tee "${ROOT}/monitor/babysitter_stopped.txt"
    exit 1
fi

case "${state}" in
    PREEMPTED|NODE_FAIL|TIMEOUT|BOOT_FAIL|REVOKED)
        ;;
    *)
        printf '%s STOP nonrecoverable_state=%s job=%s checkpoint=%s/%s\n' \
            "${timestamp}" "${state}" "${PREVIOUS_JOB}" \
            "${checkpoint_step}" "${TARGET_STEP}" |
            tee "${ROOT}/monitor/babysitter_stopped.txt"
        exit 1
        ;;
esac

retry=0
if [ -s "${RETRY_FILE}" ]; then
    read -r retry < "${RETRY_FILE}"
fi
if ! [[ "${retry}" =~ ^[0-9]+$ ]]; then
    echo "invalid retry counter: ${retry}" >&2
    exit 1
fi
if [ "${retry}" -ge "${MAX_RETRIES}" ]; then
    printf '%s STOP retry_limit=%s state=%s job=%s checkpoint=%s/%s\n' \
        "${timestamp}" "${MAX_RETRIES}" "${state}" "${PREVIOUS_JOB}" \
        "${checkpoint_step}" "${TARGET_STEP}" |
        tee "${ROOT}/monitor/babysitter_stopped.txt"
    exit 1
fi

retry=$((retry + 1))
printf '%s\n' "${retry}" > "${RETRY_FILE}"
new_job=$(sbatch --parsable "${RUN_SCRIPT}")
new_watch=$(sbatch --parsable --dependency="afterany:${new_job}" \
    --export="ALL,WATCH_JOB_ID=${new_job}" "${WATCH_SCRIPT}")
printf '%s RETRY=%s prior=%s state=%s checkpoint=%s/%s new_job=%s watcher=%s\n' \
    "${timestamp}" "${retry}" "${PREVIOUS_JOB}" "${state}" \
    "${checkpoint_step}" "${TARGET_STEP}" "${new_job}" "${new_watch}" |
    tee "${STATUS_FILE}"
