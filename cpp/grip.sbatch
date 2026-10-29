#!/bin/bash
#SBATCH --job-name=gripspan
#SBATCH --partition=owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=00:40:00
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/ksoff27_prod_20261003/logs/grip_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/ksoff27_prod_20261003/logs/grip_%j.err
set -euo pipefail
echo "host=$(hostname) start=$(date '+%F %T')"
python3 /home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/grip.py
echo "end=$(date '+%F %T')"
