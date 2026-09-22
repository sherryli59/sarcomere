#!/bin/bash
#SBATCH --partition=owners,rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=00:10:00
#SBATCH --job-name=noalign_babycron
#SBATCH --output=/scratch/users/xsli/sarcomere/production_noalign_5s_20260921/logs/babycron_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/production_noalign_5s_20260921/logs/babycron_%j.err
#SBATCH --no-requeue
#
# Sherlock does not permit this account to use Unix crontab, and rejects
# sleeper jobs.  This is the local cron convention already used across this
# repo: scan once, schedule the next scan with --begin, release the allocation
# immediately.  Stops when all 12 trajectories are complete.
set -uo pipefail
ROOT=/scratch/users/xsli/sarcomere/production_noalign_5s_20260921
SCAN=${ROOT}/source_snapshot/noalign_babysitter.sh
SELF=${ROOT}/source_snapshot/noalign_babycron.sbatch
INTERVAL_MIN=${NOALIGN_INTERVAL_MIN:-20}

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH

echo "noalign babycron scan $(date)"
bash "${SCAN}" || echo "scan returned nonzero; chain continues"

COMPLETE=$(python3 - "${ROOT}" <<'PY'
import os, glob, sys
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
import h5py
n = 0
for p in glob.glob(os.path.join(sys.argv[1], "data", "*.h5")):
    try:
        with h5py.File(p, "r", locking=False) as f:
            if int(f["meta/committed_frames"][0]) >= 1001:
                n += 1
    except Exception:
        pass
print(n)
PY
)
if [ "${COMPLETE:-0}" -ge 12 ]; then
    echo "all 12 trajectories complete; babycron chain stopping"
else
    echo "${COMPLETE:-0}/12 complete; next scan in ${INTERVAL_MIN} min"
    sbatch --begin="now+${INTERVAL_MIN}minutes" "${SELF}"
fi
