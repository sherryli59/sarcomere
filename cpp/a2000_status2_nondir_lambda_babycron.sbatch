#!/bin/bash
#SBATCH --partition=owners,rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=00:10:00
#SBATCH --job-name=a2k_s2_nd_babycron
#SBATCH --output=/scratch/users/xsli/sarcomere/logs_a2000_status2_nondir_lambda/babycron_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs_a2000_status2_nondir_lambda/babycron_%j.err
#SBATCH --no-requeue

# Self-chaining Slurm watchdog used because Unix crontab is unavailable.

set -euo pipefail

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
DATA_DIR="${A2K_ND_OUTPUT_DIR:-/scratch/users/xsli/sarcomere/data_a2000_status2_nondir_lambda}"
TARGET_TOTAL_STEPS="${A2K_ND_NSTEPS:-300000}"
INTERVAL_MIN="${A2K_ND_WATCHDOG_INTERVAL_MIN:-30}"

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
mkdir -p /scratch/users/xsli/sarcomere/logs_a2000_status2_nondir_lambda
cd "${CPP_DIR}"

echo "A2000 nondirectional babycron scan $(date)"
bash a2000_status2_nondir_lambda_babysitter.sh || echo "watchdog scan returned nonzero; chain will continue"

LIVE_JOB_NAMES=$(squeue --me -h -o '%j' 2>/dev/null || true)
if grep -Eq '^(a2k_s2_nd_lambda|A2000_M1000_s2only_lt.*_nodir_r[0-9][0-9]_L12x2x2)$' <<< "${LIVE_JOB_NAMES}"; then
    echo "nondirectional workers remain live; next scan in ${INTERVAL_MIN} min"
    sbatch --begin="now+${INTERVAL_MIN}minutes" a2000_status2_nondir_lambda_babycron.sbatch
    exit 0
fi

# No live writer remains, so it is safe to inspect completion checkpoints.
read -r COMPLETE EXPECTED <<< "$(python - "${DATA_DIR}" "${TARGET_TOTAL_STEPS}" <<'PY'
import sys
from pathlib import Path

import h5py

directory = Path(sys.argv[1])
target = int(sys.argv[2])
labels = (
    "0", "1e-5.5", "1e-5", "1e-4.5", "1e-4", "1e-3.5", "1e-3",
    "1e-2.5", "1e-2", "1e-1.5", "1e-1", "1e-0.5", "1e0",
)
complete = 0
for rep in range(3):
    for label in labels:
        path = directory / f"A2000_M1000_s2only_lt{label}_nodir_r{rep:02d}_L12x2x2.h5"
        if not path.exists():
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
print(complete, 3 * len(labels))
PY
)"

if [ "${COMPLETE:-0}" -ge "${EXPECTED:-39}" ]; then
    echo "all ${EXPECTED} trajectories complete; watchdog chain stopping"
else
    echo "${COMPLETE:-0}/${EXPECTED:-39} trajectories complete; next scan in ${INTERVAL_MIN} min"
    sbatch --begin="now+${INTERVAL_MIN}minutes" a2000_status2_nondir_lambda_babycron.sbatch
fi
