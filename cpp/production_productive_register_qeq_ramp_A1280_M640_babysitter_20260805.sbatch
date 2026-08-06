#!/bin/bash
# One-shot watcher. Slurm auto-requeue remains disabled; this only submits a
# validated checkpoint continuation after recoverable scheduler failures.
#SBATCH --partition=owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=00:10:00
#SBATCH --job-name=sarc_pr_qeq_watch
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/prod_productive_register_qeq_ramp_20260805/logs/%x_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/prod_productive_register_qeq_ramp_20260805/logs/%x_%j.err

set -euo pipefail

ROOT=/scratch/users/xsli/sarcomere/prod_productive_register_qeq_ramp_20260805
RUN_SCRIPT=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp/production_productive_register_qeq_ramp_A1280_M640_100ms_20260805.sbatch
WATCH_SCRIPT=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp/production_productive_register_qeq_ramp_A1280_M640_babysitter_20260805.sbatch
OUTPUT=${ROOT}/data/A1280_M640_productive_register_qeq_ramp_seed20260805_100ms.h5
TARGET_STEP=200000
MAX_RETRIES=4
PREVIOUS_JOB=${WATCH_JOB_ID:?WATCH_JOB_ID must name the production attempt}
MONITOR_ROOT=${ROOT}/monitor/default
RETRY_FILE=${MONITOR_ROOT}/retry_count
STATUS_FILE=${MONITOR_ROOT}/babysitter_status.txt
mkdir -p "${MONITOR_ROOT}" "${ROOT}/logs"

exec 8>"${MONITOR_ROOT}/babysitter.lock"
flock -n 8 || { echo "another babysitter owns the monitor lock"; exit 0; }

state=""
for _ in $(seq 1 24); do
    state=$(sacct -X -j "${PREVIOUS_JOB}" -n -P -o State 2>/dev/null |
        head -n 1 | cut -d'|' -f1 | sed 's/+.*//')
    [ -n "${state}" ] && break
    sleep 5
done
[ -n "${state}" ] || { echo "could not resolve ${PREVIOUS_JOB}" >&2; exit 1; }

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
        meta = h["/meta"].attrs
        if (int(meta.get("am_q_initialized_at_calcium_equilibrium", -1)) != 1 or
                int(meta.get("interaction_ramp_steps", -1)) != 40000 or
                int(meta.get("aa_status_mode", -1)) != 7 or
                abs(float(meta.get("aa_status2_registry_on", -1)) - 0.040) > 1e-12 or
                abs(float(meta.get("myosin_radius", -1)) - 0.025) > 1e-12):
            print(-1)
            raise SystemExit
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
        "${timestamp}" "${PREVIOUS_JOB}" "${state}" \
        "${checkpoint_step}" "${TARGET_STEP}" |
        tee "${MONITOR_ROOT}/production_complete.txt"
    exit 0
fi

case "${state}" in
    PREEMPTED|NODE_FAIL|TIMEOUT|BOOT_FAIL|REVOKED)
        ;;
    *)
        printf '%s STOP state=%s job=%s checkpoint=%s/%s\n' \
            "${timestamp}" "${state}" "${PREVIOUS_JOB}" \
            "${checkpoint_step}" "${TARGET_STEP}" |
            tee "${MONITOR_ROOT}/babysitter_stopped.txt"
        exit 1
        ;;
esac

retry=0
if [ -s "${RETRY_FILE}" ]; then read -r retry <"${RETRY_FILE}"; fi
[[ "${retry}" =~ ^[0-9]+$ ]] || { echo "invalid retry count" >&2; exit 1; }
if [ "${retry}" -ge "${MAX_RETRIES}" ]; then
    echo "retry limit ${MAX_RETRIES} reached" | tee "${MONITOR_ROOT}/babysitter_stopped.txt"
    exit 1
fi
retry=$((retry + 1))
printf '%s\n' "${retry}" >"${RETRY_FILE}"
new_job=$(sbatch --parsable "${RUN_SCRIPT}")
new_watch=$(sbatch --parsable --dependency="afterany:${new_job}" \
    --export="ALL,WATCH_JOB_ID=${new_job}" "${WATCH_SCRIPT}")
printf '%s RETRY=%s prior=%s state=%s checkpoint=%s/%s new_job=%s watcher=%s\n' \
    "${timestamp}" "${retry}" "${PREVIOUS_JOB}" "${state}" \
    "${checkpoint_step}" "${TARGET_STEP}" "${new_job}" "${new_watch}" |
    tee "${STATUS_FILE}"
