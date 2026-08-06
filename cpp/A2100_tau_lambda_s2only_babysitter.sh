#!/bin/bash
# One watchdog scan for the A2100 tau-lambda status-2-only sweep.

set -euo pipefail

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
SWEEP_SBATCH="${TAULAM_SWEEP_SBATCH:-${CPP_DIR}/A2100_tau_lambda_s2only.sbatch}"
DATA_DIR="${TAULAM_OUTPUT_DIR:-/scratch/users/xsli/sarcomere/data_A2100_tau_lambda_s2only}"
TARGET_TOTAL_STEPS="${TAULAM_NSTEPS:-300000}"
STALLED_MINUTES="${TAULAM_STALLED_MINUTES:-30}"
LOCK_FILE="${TAULAM_WATCHDOG_LOCK_FILE:-/scratch/users/xsli/sarcomere/.A2100_tau_lambda_s2only_babysitter.lock}"
KINETICS_TAG="${TAULAM_KINETICS_TAG:-s2tau}"
GENERIC_JOB_NAME="${TAULAM_GENERIC_JOB_NAME:-A2100_tau_lambda}"

mkdir -p "${DATA_DIR}"
exec 8>"${LOCK_FILE}"
if ! flock -n 8; then
    echo "another A2100 tau-lambda watchdog scan is active; exiting"
    exit 0
fi

LIVE_JOB_NAMES=$(squeue --me -h -o '%j' 2>/dev/null || true)
if grep -Fxq "${GENERIC_JOB_NAME}" <<< "${LIVE_JOB_NAMES}"; then
    echo "the A2100 tau-lambda array still has pending/unrenamed elements; deferring retries"
    exit 0
fi

mapfile -t RETRY_IDS < <(HDF5_USE_FILE_LOCKING=FALSE python - "${DATA_DIR}" "${TARGET_TOTAL_STEPS}" "${STALLED_MINUTES}" "${LIVE_JOB_NAMES}" "${KINETICS_TAG}" <<'PY'
import sys
import time
from pathlib import Path

import h5py

directory = Path(sys.argv[1])
target = int(sys.argv[2])
stalled_seconds = int(sys.argv[3]) * 60
live_names = set(sys.argv[4].splitlines())
tag = sys.argv[5]
taus = ("0", "1e-5", "1e-4")
lambdas = ("0", "1e-5", "1e-4", "1e-3", "1e-2", "1e-1")
now = time.time()

for rep in range(3):
    for tau_idx, tau in enumerate(taus):
        for lambda_idx, lam in enumerate(lambdas):
            task_id = rep * len(taus) * len(lambdas) + tau_idx * len(lambdas) + lambda_idx
            stem = f"A2100_M1000_{tag}_tau{tau}_lt{lam}_r{rep:02d}_L12x2x2"
            path = directory / f"{stem}.h5"
            if stem in live_names:
                continue
            if not path.exists() or path.stat().st_size == 0:
                print(task_id)
                continue
            if now - path.stat().st_mtime < stalled_seconds:
                continue
            try:
                with h5py.File(path, "r", locking=False) as handle:
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
    echo "no missing or stalled A2100 tau-lambda tasks need submission"
    exit 0
fi

ARRAY_SPEC=$(IFS=,; echo "${RETRY_IDS[*]}")
echo "submitting task IDs ${ARRAY_SPEC} without an array concurrency throttle"
if [ "${TAULAM_WATCHDOG_DRY_RUN:-0}" = 1 ]; then
    echo "dry run: sbatch --array=${ARRAY_SPEC} ${SWEEP_SBATCH}"
    exit 0
fi
sbatch --partition=owners --account=rotskoff --array="${ARRAY_SPEC}" "${SWEEP_SBATCH}"
