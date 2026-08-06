#!/bin/bash
# One watchdog scan for the A2000/M1000 constant-tau1 sweep.
# Submit only exact missing/stalled array IDs, never a task with a live writer.

set -euo pipefail

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
SWEEP_SBATCH="${TAU1_SWEEP_SBATCH:-${CPP_DIR}/a2000_tau1_sweep.sbatch}"
DATA_DIR="${TAU1_OUTPUT_DIR:-/scratch/users/xsli/sarcomere/data_a2000_tau1_sweep}"
TARGET_TOTAL_STEPS="${TAU1_NSTEPS:-300000}"
STALLED_MINUTES="${TAU1_STALLED_MINUTES:-30}"
LOCK_FILE="${TAU1_WATCHDOG_LOCK_FILE:-/scratch/users/xsli/sarcomere/.a2000_tau1_sweep_babysitter.lock}"
KINETICS_TAG="${TAU1_KINETICS_TAG:-s2catch}"
GENERIC_JOB_NAME="${TAU1_GENERIC_JOB_NAME:-a2k_tau1}"
WATCHDOG_LABEL="${TAU1_WATCHDOG_LABEL:-A2000 constant-tau1}"

mkdir -p "${DATA_DIR}"
exec 8>"${LOCK_FILE}"
if ! flock -n 8; then
    echo "another ${WATCHDOG_LABEL} watchdog scan is active; exiting"
    exit 0
fi

LIVE_JOB_NAMES=$(squeue --me -h -o '%j' 2>/dev/null || true)

# Pending array elements still carry the generic sbatch job name because the
# worker has not run its per-file rename yet. Defer the whole scan in that case;
# otherwise missing files from pending elements would look retryable.
if grep -Fxq "${GENERIC_JOB_NAME}" <<< "${LIVE_JOB_NAMES}"; then
    echo "a ${WATCHDOG_LABEL} array still has pending/unrenamed elements; deferring retries"
    exit 0
fi

mapfile -t RETRY_IDS < <(python - "${DATA_DIR}" "${TARGET_TOTAL_STEPS}" "${STALLED_MINUTES}" "${LIVE_JOB_NAMES}" "${KINETICS_TAG}" <<'PY'
import sys
import time
from pathlib import Path

import h5py

directory = Path(sys.argv[1])
target = int(sys.argv[2])
stalled_seconds = int(sys.argv[3]) * 60
live_names = set(sys.argv[4].splitlines())
kinetics_tag = sys.argv[5]

# Must mirror a2000_tau1_sweep.sbatch exactly.
tau_labels = ("0", "1e-4.5", "1e-4", "1e-3.5", "1e-3", "1e-2.5")
lambda_labels = ("1e-3", "1e-2")
now = time.time()

per_rep = len(tau_labels) * len(lambda_labels)
for rep in range(2):
    for lambda_idx, lam in enumerate(lambda_labels):
        for tau_idx, tau in enumerate(tau_labels):
            task_id = rep * per_rep + lambda_idx * len(tau_labels) + tau_idx
            stem = (f"A2000_M1000_{kinetics_tag}_tau{tau}_lt{lam}"
                    f"_r{rep:02d}_L12x2x2")
            path = directory / f"{stem}.h5"

            # Main and retry workers rename themselves to this exact stem.
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
                # A stale unreadable file is retried; the worker preserves it
                # before starting fresh, so no trajectory is silently discarded.
                step = -1

            if step < target:
                print(task_id)
PY
)

if [ "${#RETRY_IDS[@]}" -eq 0 ]; then
    echo "no missing or stalled ${WATCHDOG_LABEL} tasks need submission"
    exit 0
fi

ARRAY_SPEC=$(IFS=,; echo "${RETRY_IDS[*]}")
echo "submitting task IDs ${ARRAY_SPEC} without an array concurrency throttle"
if [ "${TAU1_WATCHDOG_DRY_RUN:-0}" = 1 ]; then
    echo "dry run: sbatch --array=${ARRAY_SPEC} ${SWEEP_SBATCH}"
    exit 0
fi
sbatch --array="${ARRAY_SPEC}" "${SWEEP_SBATCH}"
