#!/bin/bash
#SBATCH --partition=owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=00:10:00
#SBATCH --job-name=A2100_tau_lambda_cron
#SBATCH --output=/scratch/users/xsli/sarcomere/logs_A2100_tau_lambda_s2only/babycron_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs_A2100_tau_lambda_s2only/babycron_%j.err
#SBATCH --no-requeue

set -euo pipefail

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
DATA_DIR="${TAULAM_OUTPUT_DIR:-/scratch/users/xsli/sarcomere/data_A2100_tau_lambda_s2only}"
TARGET_TOTAL_STEPS="${TAULAM_NSTEPS:-300000}"
INTERVAL_MIN="${TAULAM_WATCHDOG_INTERVAL_MIN:-30}"
KINETICS_TAG="${TAULAM_KINETICS_TAG:-s2tau}"
WATCHDOG_SCRIPT="${TAULAM_WATCHDOG_SCRIPT:-${CPP_DIR}/A2100_tau_lambda_s2only_babysitter.sh}"
BABYCRON_SBATCH="${TAULAM_BABYCRON_SBATCH:-${CPP_DIR}/A2100_tau_lambda_s2only_babycron.sbatch}"
LOG_DIR="${TAULAM_LOG_DIR:-/scratch/users/xsli/sarcomere/logs_A2100_tau_lambda_s2only}"

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
mkdir -p "${LOG_DIR}"

echo "A2100 tau-lambda babycron scan $(date)"
bash "${WATCHDOG_SCRIPT}" || echo "watchdog scan returned nonzero; chain will continue"

LIVE_JOB_NAMES=$(squeue --me -h -o '%j' 2>/dev/null || true)
read -r COMPLETE EXPECTED <<< "$(HDF5_USE_FILE_LOCKING=FALSE python - "${DATA_DIR}" "${TARGET_TOTAL_STEPS}" "${LIVE_JOB_NAMES}" "${KINETICS_TAG}" <<'PY'
import sys
from pathlib import Path

import h5py

directory = Path(sys.argv[1])
target = int(sys.argv[2])
live_names = set(sys.argv[3].splitlines())
tag = sys.argv[4]
taus = ("0", "1e-5", "1e-4")
lambdas = ("0", "1e-5", "1e-4", "1e-3", "1e-2", "1e-1")
complete = 0
expected = 0
for rep in range(3):
    for tau in taus:
        for lam in lambdas:
            expected += 1
            stem = f"A2100_M1000_{tag}_tau{tau}_lt{lam}_r{rep:02d}_L12x2x2"
            path = directory / f"{stem}.h5"
            if stem in live_names or not path.exists():
                continue
            try:
                with h5py.File(path, "r", locking=False) as handle:
                    if "/resume/current_step" in handle:
                        step = int(handle["/resume/current_step"][...].reshape(-1)[-1])
                    elif "/state/current_step" in handle:
                        step = int(handle["/state/current_step"][...].reshape(-1)[-1])
                    else:
                        step = -1
                complete += step >= target
            except Exception:
                pass
print(complete, expected)
PY
)"

if [ "${COMPLETE:-0}" -ge "${EXPECTED:-54}" ]; then
    echo "all ${EXPECTED} trajectories complete; watchdog chain stopping"
else
    echo "${COMPLETE:-0}/${EXPECTED:-54} trajectories complete; next scan in ${INTERVAL_MIN} min"
    sbatch --partition=owners --account=rotskoff --begin="now+${INTERVAL_MIN}minutes" "${BABYCRON_SBATCH}"
fi
