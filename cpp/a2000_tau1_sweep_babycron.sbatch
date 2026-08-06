#!/bin/bash
#SBATCH --partition=owners,rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=00:10:00
#SBATCH --job-name=a2k_tau1_babycron
#SBATCH --output=/scratch/users/xsli/sarcomere/logs_a2000_tau1_sweep/babycron_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs_a2000_tau1_sweep/babycron_%j.err
#SBATCH --no-requeue

# Sherlock does not permit this account to use Unix crontab. This is the local
# cron convention: scan once, schedule the next scan with --begin, and release
# the allocation immediately. Stop only after all 24 expected files complete.

set -euo pipefail

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
DATA_DIR="${TAU1_OUTPUT_DIR:-/scratch/users/xsli/sarcomere/data_a2000_tau1_sweep}"
TARGET_TOTAL_STEPS="${TAU1_NSTEPS:-300000}"
INTERVAL_MIN="${TAU1_WATCHDOG_INTERVAL_MIN:-30}"
KINETICS_TAG="${TAU1_KINETICS_TAG:-s2catch}"
WATCHDOG_LABEL="${TAU1_WATCHDOG_LABEL:-A2000 constant-tau1}"
WATCHDOG_SCRIPT="${TAU1_WATCHDOG_SCRIPT:-${CPP_DIR}/a2000_tau1_sweep_babysitter.sh}"
BABYCRON_SBATCH="${TAU1_BABYCRON_SBATCH:-${CPP_DIR}/a2000_tau1_sweep_babycron.sbatch}"
LOG_DIR="${TAU1_LOG_DIR:-/scratch/users/xsli/sarcomere/logs_a2000_tau1_sweep}"

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
mkdir -p "${LOG_DIR}"
cd "${CPP_DIR}"

echo "${WATCHDOG_LABEL} babycron scan $(date)"
bash "${WATCHDOG_SCRIPT}" || echo "watchdog scan returned nonzero; chain will continue"

LIVE_JOB_NAMES=$(squeue --me -h -o '%j' 2>/dev/null || true)
read -r COMPLETE EXPECTED <<< "$(python - "${DATA_DIR}" "${TARGET_TOTAL_STEPS}" "${LIVE_JOB_NAMES}" "${KINETICS_TAG}" <<'PY'
import sys
from pathlib import Path

import h5py

directory = Path(sys.argv[1])
target = int(sys.argv[2])
live_names = set(sys.argv[3].splitlines())
kinetics_tag = sys.argv[4]

# Must mirror a2000_tau1_sweep.sbatch exactly.
tau_labels = ("0", "1e-4.5", "1e-4", "1e-3.5", "1e-3", "1e-2.5")
lambda_labels = ("1e-3", "1e-2")

complete = 0
for rep in range(2):
    for lam in lambda_labels:
        for tau in tau_labels:
            stem = (f"A2000_M1000_{kinetics_tag}_tau{tau}_lt{lam}"
                    f"_r{rep:02d}_L12x2x2")
            path = directory / f"{stem}.h5"
            if stem in live_names or not path.exists():
                continue
            try:
                with h5py.File(path, "r") as handle:
                    if "/resume/current_step" in handle:
                        step = int(handle["/resume/current_step"][...].reshape(-1)[-1])
                    elif "/state/current_step" in handle:
                        step = int(handle["/state/current_step"][...].reshape(-1)[-1])
                    else:
                        step = -1
                complete += step >= target
            except Exception:
                pass
print(complete, 2 * len(tau_labels) * len(lambda_labels))
PY
)"

if [ "${COMPLETE:-0}" -ge "${EXPECTED:-24}" ]; then
    echo "all ${EXPECTED} trajectories complete; watchdog chain stopping"
else
    echo "${COMPLETE:-0}/${EXPECTED:-24} trajectories complete; next scan in ${INTERVAL_MIN} min"
    sbatch --begin="now+${INTERVAL_MIN}minutes" "${BABYCRON_SBATCH}"
fi
