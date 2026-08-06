#!/bin/bash
#SBATCH --partition=owners,rotskoff
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=00:10:00
#SBATCH --job-name=a2k_explicit_cron
#SBATCH --output=/scratch/users/xsli/sarcomere/logs_a2000_explicit_crossbridge/%j_cron.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs_a2000_explicit_crossbridge/%j_cron.err
#SBATCH --no-requeue

set -euo pipefail
CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
mkdir -p /scratch/users/xsli/sarcomere/logs_a2000_explicit_crossbridge
eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export HDF5_USE_FILE_LOCKING=FALSE
bash "${CPP_DIR}/a2000_explicit_crossbridge_babysitter.sh"

ANALYSIS_DIR="${EXPLICIT_CB_ANALYSIS_DIR:-/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/a2000_explicit_crossbridge_v3}"
if [ ! -e "${ANALYSIS_DIR}/.analysis_complete" ]; then
    sbatch --begin=now+20minutes "${CPP_DIR}/a2000_explicit_crossbridge_babycron.sbatch"
fi
