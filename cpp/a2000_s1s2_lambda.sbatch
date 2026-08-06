#!/bin/bash
#SBATCH --partition=owners,rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=48GB
#SBATCH --time=2-00:00:00
#SBATCH --job-name=a2k_s1s2_lambda
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/logs_a2000_s1s2_lambda/%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs_a2000_s1s2_lambda/%A_%a.err
#SBATCH --array=0-25

# Parameter-matched A2000/M1000 sweep using the original catch-bond law:
# lambda*f_load stabilizes status 1 and status 2. All other settings are
# delegated to the same worker used by the status-2-only comparison.

set -euo pipefail

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
SCRATCH_DIR=/scratch/users/xsli/sarcomere
export A2K_STABILIZE_STATUS1=true
export A2K_KINETICS_TAG=s1s2
export A2K_OUTPUT_DIR="${SCRATCH_DIR}/data_a2000_s1s2_lambda"
export A2K_LOG_DIR="${SCRATCH_DIR}/logs_a2000_s1s2_lambda"
export A2K_PARAM_DIR="${SCRATCH_DIR}/params_a2000_s1s2_lambda"

exec bash "${CPP_DIR}/a2000_status2_lambda.sbatch"
