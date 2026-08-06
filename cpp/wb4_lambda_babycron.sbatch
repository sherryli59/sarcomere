#!/bin/bash
#SBATCH -p owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=00:10:00
#SBATCH --job-name=wb4_lambda_cron
#SBATCH --output=/scratch/users/xsli/sarcomere/logs/wb4_lambda_cron_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs/wb4_lambda_cron_%j.err

# Self-chaining, cron-like watchdog: scan once, then schedule the next scan only
# while at least one existing 12x4x4 lambda trajectory remains incomplete.

set -uo pipefail

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
DATA_DIR=/scratch/users/xsli/sarcomere/data_wb4x4_lambda
TARGET_TOTAL_STEPS=300000
SAVE_EVERY=200
INTERVAL_MIN="${INTERVAL_MIN:-30}"

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
cd "${CPP_DIR}"

echo "wb4-lambda cron scan $(date)"
bash wb4_lambda_babysitter.sh || echo "  babysitter returned nonzero; chain will continue"

LIVE_JOB_NAMES=$(squeue --me -h -o "%j" 2>/dev/null || true)
read -r TOTAL INCOMPLETE <<< "$(python - "${DATA_DIR}" "${TARGET_TOTAL_STEPS}" "${SAVE_EVERY}" "${LIVE_JOB_NAMES}" <<'PY'
import glob
import sys
import h5py
from pathlib import Path

directory, target, save_every = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
live_names = set(sys.argv[4].splitlines())
paths = glob.glob(f"{directory}/A*_M*_lt*_r*_L12x4x4.h5")
incomplete = 0
for path in paths:
    # Opening a live writer here can take a read lock between its periodic
    # close/reopen cycle and reproduce the errno-11 abort we are recovering.
    if Path(path).stem in live_names:
        incomplete += 1
        continue
    try:
        with h5py.File(path, "r") as handle:
            if "/resume/current_step" in handle:
                step = int(handle["/resume/current_step"][...].reshape(-1)[-1])
            elif "/state/current_step" in handle:
                step = int(handle["/state/current_step"][...].reshape(-1)[-1])
            else:
                step = (int(handle["/actin/center"].shape[0]) - 1) * save_every
        incomplete += step < target
    except Exception:
        incomplete += 1
print(len(paths), incomplete)
PY
)"

if [ "${TOTAL:-0}" -eq 0 ]; then
    echo "no 12x4x4 lambda trajectories found; next scan in ${INTERVAL_MIN} min"
    sbatch --begin="now+${INTERVAL_MIN}minutes" wb4_lambda_babycron.sbatch
elif [ "${INCOMPLETE:-1}" -gt 0 ]; then
    echo "${INCOMPLETE}/${TOTAL} trajectories incomplete; next scan in ${INTERVAL_MIN} min"
    sbatch --begin="now+${INTERVAL_MIN}minutes" wb4_lambda_babycron.sbatch
else
    echo "all ${TOTAL} trajectories complete; watchdog chain stopping"
fi
