#!/bin/bash
# Watchdog configuration for the original status-1 + status-2 catch-bond law.

set -euo pipefail

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
SCRATCH_DIR=/scratch/users/xsli/sarcomere
export A2K_STABILIZE_STATUS1=true
export A2K_KINETICS_TAG=s1s2
export A2K_OUTPUT_DIR="${SCRATCH_DIR}/data_a2000_s1s2_lambda"
export A2K_LOG_DIR="${SCRATCH_DIR}/logs_a2000_s1s2_lambda"
export A2K_PARAM_DIR="${SCRATCH_DIR}/params_a2000_s1s2_lambda"
export A2K_SWEEP_SBATCH="${CPP_DIR}/a2000_s1s2_lambda.sbatch"
export A2K_GENERIC_JOB_NAME=a2k_s1s2_lambda
export A2K_WATCHDOG_LABEL="A2000 status-1+2"
export A2K_WATCHDOG_LOCK_FILE="${SCRATCH_DIR}/.a2000_s1s2_lambda_babysitter.lock"

exec bash "${CPP_DIR}/a2000_status2_lambda_babysitter.sh"
