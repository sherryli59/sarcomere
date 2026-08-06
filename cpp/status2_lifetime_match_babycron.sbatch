#!/bin/bash
#SBATCH --partition=owners,rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=00:10:00
#SBATCH --job-name=s2_match_babycron
#SBATCH --output=/scratch/users/xsli/sarcomere/logs/s2_match_babycron_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs/s2_match_babycron_%j.err
#SBATCH --no-requeue

# Self-chaining Slurm watchdog used because Unix crontab is unavailable.

set -euo pipefail

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
DATA_DIR="${CONTROL_OUTPUT_DIR:-/scratch/users/xsli/sarcomere/data_status2_lifetime_match}"
TARGET_TOTAL_STEPS="${CONTROL_NSTEPS:-150000}"
N_ACTINS="${CONTROL_N_ACTINS:-2000}"
INTERVAL_MIN="${S2_MATCH_WATCHDOG_INTERVAL_MIN:-30}"

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
cd "${CPP_DIR}"

echo "status-2 lifetime-match babycron scan $(date)"
bash status2_lifetime_match_babysitter.sh || echo "watchdog scan returned nonzero; chain will continue"

if squeue --me -h -n s2_lifetime_match -o '%i' 2>/dev/null | grep -q .; then
    echo "lifetime-match workers remain live; next scan in ${INTERVAL_MIN} min"
    sbatch --begin="now+${INTERVAL_MIN}minutes" status2_lifetime_match_babycron.sbatch
    exit 0
fi

read -r COMPLETE EXPECTED <<< "$(python - "${DATA_DIR}" "${TARGET_TOTAL_STEPS}" "${N_ACTINS}" <<'PY'
import sys
from pathlib import Path

import h5py

directory = Path(sys.argv[1])
target = int(sys.argv[2])
n_actins = int(sys.argv[3])
tau_values = ("0.001", "0.002", "0.003", "0.004", "0.006")
kon_values = (500, 1000, 2000)
tags = ["catch_lt0p003162_kon1000"]
tags += [f"consts2_tau{tau.replace('.', 'p')}_kon{kon}"
         for tau in tau_values for kon in kon_values]

complete = 0
for tag in tags:
    for rep in range(2):
        path = directory / f"A{n_actins}_M1000_{tag}_r{rep:02d}_L12x2x2.h5"
        if not path.exists():
            continue
        try:
            with h5py.File(path, "r") as handle:
                if "/resume/current_step" in handle:
                    step = int(handle["/resume/current_step"][...].reshape(-1)[-1])
                elif "/resume_prev/current_step" in handle:
                    step = int(handle["/resume_prev/current_step"][...].reshape(-1)[-1])
                elif "/state/current_step" in handle:
                    step = int(handle["/state/current_step"][...].reshape(-1)[-1])
                else:
                    step = -1
            complete += step >= target
        except Exception:
            pass
print(complete, 2 * len(tags))
PY
)"

if [ "${COMPLETE:-0}" -ge "${EXPECTED:-32}" ]; then
    echo "all ${EXPECTED} trajectories complete; watchdog chain stopping"
else
    echo "${COMPLETE:-0}/${EXPECTED:-32} trajectories complete; next scan in ${INTERVAL_MIN} min"
    sbatch --begin="now+${INTERVAL_MIN}minutes" status2_lifetime_match_babycron.sbatch
fi
