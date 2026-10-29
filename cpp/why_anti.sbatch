#!/bin/bash
#SBATCH --job-name=whyanti
#SBATCH --partition=owners,normal,hns
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=00:40:00
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/amalign_20261005/logs/early_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/amalign_20261005/logs/early_%j.err
set -euo pipefail
export PATH=/home/groups/rotskoff/sherry/mamba_root/envs/sarcomere/bin:${PATH}
unset PYTHONPATH
echo "host=$(hostname) start=$(date '+%F %T')"
python3 /home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/why_anti.py
echo "end=$(date '+%F %T')"
