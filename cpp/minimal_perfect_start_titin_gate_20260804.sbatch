#!/bin/bash
#SBATCH --job-name=sarc_min_perfect
#SBATCH --partition=owners
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --output=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp/diagnostics/self_loaded_balance_20260804/distance_capture_perfect_start/slurm_%j.out
#SBATCH --error=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp/diagnostics/self_loaded_balance_20260804/distance_capture_perfect_start/slurm_%j.err

set -euo pipefail
cd /home/groups/rotskoff/sherry/sarcomere/sarcomere
exec cpp/run_minimal_perfect_start_titin_gate_20260804.sh
