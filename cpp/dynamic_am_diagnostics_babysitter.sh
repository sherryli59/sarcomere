#!/bin/bash
# One watchdog scan for a dynamic-AM minimal diagnostic array.

set -euo pipefail

MODE=${DIAG_MODE:?DIAG_MODE must be single, crosslink, or memory}
CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
SCRATCH_DIR=/scratch/users/xsli/sarcomere
TARGET_TOTAL_STEPS=${DIAG_NSTEPS:-200000}
STALLED_MINUTES=${DIAG_STALLED_MINUTES:-15}

case "${MODE}" in
    single)
        SWEEP_SBATCH=${CPP_DIR}/dynamic_am_single_motor.sbatch
        DATA_DIR=${AM_SINGLE_OUTPUT_DIR:-${SCRATCH_DIR}/data_am_single_motor}
        GENERIC_JOB_NAME=am_single_motor
        ;;
    crosslink)
        SWEEP_SBATCH=${CPP_DIR}/dynamic_am_minimal_crosslink.sbatch
        DATA_DIR=${AM_XLINK_OUTPUT_DIR:-${SCRATCH_DIR}/data_am_minimal_crosslink}
        GENERIC_JOB_NAME=am_minimal_xlink
        ;;
    memory)
        SWEEP_SBATCH=${CPP_DIR}/dynamic_am_status2_memory.sbatch
        DATA_DIR=${AM_MEMORY_OUTPUT_DIR:-${SCRATCH_DIR}/data_am_status2_memory}
        GENERIC_JOB_NAME=am_status2_memory
        ;;
    *)
        echo "unknown DIAG_MODE=${MODE}" >&2
        exit 1
        ;;
esac

mkdir -p "${DATA_DIR}"
exec 8>"${SCRATCH_DIR}/.${MODE}_dynamic_am_diagnostics.lock"
flock -n 8 || exit 0

LIVE_JOB_NAMES=$(squeue --me -h -o '%j' 2>/dev/null || true)
if grep -Fxq "${GENERIC_JOB_NAME}" <<<"${LIVE_JOB_NAMES}"; then
    echo "${GENERIC_JOB_NAME} still has live array elements; deferring retries"
    exit 0
fi

mapfile -t RETRY_IDS < <(HDF5_USE_FILE_LOCKING=FALSE python - \
    "${MODE}" "${DATA_DIR}" "${TARGET_TOTAL_STEPS}" "${STALLED_MINUTES}" <<'PY'
import sys
import time
from pathlib import Path

import h5py

mode, directory, target, stalled_minutes = sys.argv[1:]
directory = Path(directory)
target = int(target)
stalled_seconds = int(stalled_minutes) * 60
now = time.time()

total_tasks = 160 if mode == "memory" else 80
for task_id in range(total_tasks):
    arm, rep = divmod(task_id, 20)
    if mode == "single":
        labels = ("0", "0p5", "1", "2")
        stem = f"am_single_F{labels[arm]}_r{rep:02d}"
    elif mode == "crosslink":
        labels = ("free_noaa", "free_persistentaa", "fixed_noaa", "fixed_persistentaa")
        stem = f"am_xlink_{labels[arm]}_r{rep:02d}"
    else:
        geometry = ("free", "fixed")[arm // 4]
        memory = ("0", "1ms", "5ms", "20ms")[arm % 4]
        stem = f"am_xlink_{geometry}_memory{memory}_r{rep:02d}"
    path = directory / f"{stem}.h5"
    if not path.exists() or path.stat().st_size == 0:
        print(task_id)
        continue
    if now - path.stat().st_mtime < stalled_seconds:
        continue
    try:
        with h5py.File(path, "r", locking=False) as handle:
            if "/resume/current_step" in handle:
                step = int(handle["/resume/current_step"][...].reshape(-1)[-1])
            elif "/resume_prev/current_step" in handle:
                step = int(handle["/resume_prev/current_step"][...].reshape(-1)[-1])
            else:
                step = -1
    except Exception:
        step = -1
    if step < target:
        print(task_id)
PY
)

if [ "${#RETRY_IDS[@]}" -eq 0 ]; then
    echo "no ${MODE} diagnostic tasks need submission"
    exit 0
fi
ARRAY_SPEC=$(IFS=,; echo "${RETRY_IDS[*]}")
echo "resubmitting ${MODE} diagnostic task IDs ${ARRAY_SPEC}"
sbatch --partition=owners --account=rotskoff --array="${ARRAY_SPEC}" "${SWEEP_SBATCH}"
