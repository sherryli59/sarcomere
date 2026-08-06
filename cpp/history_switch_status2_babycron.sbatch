#!/bin/bash
#SBATCH --partition=owners,rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=00:10:00
#SBATCH --job-name=s2_history_babycron
#SBATCH --output=/scratch/users/xsli/sarcomere/logs/s2_history_babycron_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs/s2_history_babycron_%j.err
#SBATCH --no-requeue

set -euo pipefail

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
DATA_DIR="${HISTORY_OUTPUT_DIR:-/scratch/users/xsli/sarcomere/data_history_status2}"
TARGET_TOTAL_STEPS="${HISTORY_TOTAL_STEPS:-400000}"
N_ACTINS="${HISTORY_N_ACTINS:-2000}"
INTERVAL_MIN="${HISTORY_WATCHDOG_INTERVAL_MIN:-30}"

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
cd "${CPP_DIR}"

echo "status-2 history babycron scan $(date)"
bash history_switch_status2_babysitter.sh || echo "watchdog scan returned nonzero; chain will continue"

if squeue --me -h -n cb_history_s2 -o '%i' 2>/dev/null | grep -q .; then
    echo "history workers remain live; next scan in ${INTERVAL_MIN} min"
    sbatch --begin="now+${INTERVAL_MIN}minutes" history_switch_status2_babycron.sbatch
    exit 0
fi

read -r COMPLETE EXPECTED <<< "$(python - "${DATA_DIR}" "${TARGET_TOTAL_STEPS}" "${N_ACTINS}" <<'PY'
import sys
from pathlib import Path

import h5py

directory = Path(sys.argv[1])
target = int(sys.argv[2])
n_actins = int(sys.argv[3])
protocols = ("medium_control", "high_control", "medium_high", "high_medium_high")
complete = 0

for protocol in protocols:
    for rep in range(3):
        path = directory / f"A{n_actins}_M1000_history_{protocol}_r{rep:02d}_L12x2x2.h5"
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
print(complete, 3 * len(protocols))
PY
)"

if [ "${COMPLETE:-0}" -ge "${EXPECTED:-12}" ]; then
    echo "all ${EXPECTED} histories complete; watchdog chain stopping"
else
    echo "${COMPLETE:-0}/${EXPECTED:-12} histories complete; next scan in ${INTERVAL_MIN} min"
    sbatch --begin="now+${INTERVAL_MIN}minutes" history_switch_status2_babycron.sbatch
fi
