#!/bin/bash
#SBATCH --job-name=probescore
#SBATCH --partition=rotskoff,owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/probe_hph_20261005/logs/score_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/probe_hph_20261005/logs/score_%j.err
set -euo pipefail
V=/scratch/users/xsli/sarcomere/confined_probe_v2_20260910/data
N=/scratch/users/xsli/sarcomere/probe_hph_20261005/data
python3 /home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/probe_score.py \
  "[[\"hph9_tilt60 (ctrl)\",\"$V/span3_tilt60_seed*.h5\"],
    [\"hph27_tilt60 (NEW)\",\"$N/hph27_tilt60_seed*.h5\"],
    [\"hph9_tilt30 (ctrl)\",\"$V/span3_tilt30_seed*.h5\"],
    [\"hph27_tilt30 (NEW)\",\"$N/hph27_tilt30_seed*.h5\"]]"
