#!/bin/bash
#SBATCH --partition=owners,normal,hns
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=00:10:00
#SBATCH --job-name=nalam_s2fstab_babycron
#SBATCH --output=/scratch/users/xsli/sarcomere/logs_nalam_s2fstab_base1e-3/babycron_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs_nalam_s2fstab_base1e-3/babycron_%j.err
#SBATCH --no-requeue

# Sherlock does not permit this account to use Unix crontab.  Local cron
# convention: scan once, schedule the next scan with --begin, release the
# allocation immediately.  Stop only after all 234 expected files complete.

set -euo pipefail

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
DATA_DIR="${NAL_OUTPUT_DIR:-/scratch/users/xsli/sarcomere/data_nalam_s2fstab_base1e-3}"
TARGET_TOTAL_STEPS="${NAL_NSTEPS:-300000}"
INTERVAL_MIN="${NAL_WATCHDOG_INTERVAL_MIN:-30}"
KINETICS_TAG="${NAL_KINETICS_TAG:-s2fstab}"
WATCHDOG_LABEL="${NAL_WATCHDOG_LABEL:-nalam s2fstab}"
WATCHDOG_SCRIPT="${NAL_WATCHDOG_SCRIPT:-${CPP_DIR}/nalam_s2fstab_base1e-3_babysitter.sh}"
BABYCRON_SBATCH="${NAL_BABYCRON_SBATCH:-${CPP_DIR}/nalam_s2fstab_base1e-3_babycron.sbatch}"
LOG_DIR="${NAL_LOG_DIR:-/scratch/users/xsli/sarcomere/logs_nalam_s2fstab_base1e-3}"

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

n_actins_values = (600, 900, 1200, 1500, 1800, 2100)
labels = (
    "0", "1e-5.5", "1e-5", "1e-4.5", "1e-4", "1e-3.5", "1e-3",
    "1e-2.5", "1e-2", "1e-1.5", "1e-1", "1e-0.5", "1e0",
)

complete = 0
expected = 0
for na in n_actins_values:
    for rep in range(3):
        for label in labels:
            expected += 1
            stem = f"A{na}_M1000_{kinetics_tag}_lt{label}_r{rep:02d}_L12x2x2"
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
print(complete, expected)
PY
)"

if [ "${COMPLETE:-0}" -ge "${EXPECTED:-234}" ]; then
    echo "all ${EXPECTED} trajectories complete; watchdog chain stopping"
else
    echo "${COMPLETE:-0}/${EXPECTED:-234} trajectories complete; next scan in ${INTERVAL_MIN} min"
    sbatch --begin="now+${INTERVAL_MIN}minutes" "${BABYCRON_SBATCH}"
fi
