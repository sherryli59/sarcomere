#!/bin/bash
#SBATCH --job-name=seedscat
#SBATCH --partition=owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/ksoff27_prod_20261003/logs/seedscat_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/ksoff27_prod_20261003/logs/seedscat_%j.err
set -euo pipefail
echo "host=$(hostname) start=$(date '+%F %T')"
python3 /home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/seed_scatter.py
echo "end=$(date '+%F %T')"
