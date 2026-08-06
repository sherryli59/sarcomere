#!/bin/bash
#SBATCH --job-name=sarc_inline_ann
#SBATCH --partition=dev
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --no-requeue
#SBATCH --output=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp/diagnostics/self_loaded_balance_20260804/axial_inline_annealing/slurm_%j.out
#SBATCH --error=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp/diagnostics/self_loaded_balance_20260804/axial_inline_annealing/slurm_%j.err

set -euo pipefail
cd /home/groups/rotskoff/sherry/sarcomere/sarcomere
exec cpp/run_minimal_axial_inline_annealing_gate_20260804.sh
