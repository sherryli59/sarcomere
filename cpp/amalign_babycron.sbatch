#!/bin/bash
#SBATCH --partition=owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=00:10:00
#SBATCH --job-name=amalign_cron
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/amalign_20261005/logs/%j_cron.out
#SBATCH --error=/scratch/users/xsli/sarcomere/amalign_20261005/logs/%j_cron.err
set -euo pipefail
export PATH=/home/groups/rotskoff/sherry/mamba_root/envs/sarcomere/bin:${PATH}
unset PYTHONPATH
export HDF5_USE_FILE_LOCKING=FALSE

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
ROOT=/scratch/users/xsli/sarcomere/amalign_20261005
mkdir -p "${ROOT}/logs"
date -u +"cron tick %Y-%m-%dT%H:%M:%SZ"
bash "${CPP_DIR}/amalign_babysitter.sh"

# Re-arm unless the campaign is finished.  --no-requeue above so a preempted
# cron tick does not double-fire; the next tick picks the work up anyway.
if [[ ! -e "${ROOT}/.campaign_complete" ]]; then
    sbatch --begin=now+20minutes "${CPP_DIR}/amalign_babycron.sbatch"
else
    echo "campaign complete; not re-arming"
fi
