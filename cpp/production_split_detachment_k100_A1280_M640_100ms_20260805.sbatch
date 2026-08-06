#!/bin/bash
# Matched production comparison arm for the model-21 split detachment law.
# Identical to the k_perp=50 arm except for the stronger transverse bound-head
# spring.  Transverse energy remains excluded from stochastic detachment.
#SBATCH --partition=owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=8GB
#SBATCH --time=08:00:00
#SBATCH --job-name=sarc_splitk100_A1280_100ms
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/prod_split_detachment_k50_20260805/logs/%x_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/prod_split_detachment_k50_20260805/logs/%x_%j.err

set -euo pipefail

export PROD_ROOT=/scratch/users/xsli/sarcomere/prod_split_detachment_k50_20260805
export PROD_EXE=${PROD_ROOT}/bin/sarcomere_split_detachment_k50_model21_20260805
export PROD_EXPECTED_EXE_SHA256=a26d4943fa88d0e7dfbfd0919a73a837ac98f40fbdc2e02a888f802cf1fb2930
export PROD_EXPECTED_AM_MODEL_VERSION=21
export PROD_NSTEPS=200000
export PROD_MAX_NSTEPS=200000
export PROD_SAVE_EVERY=2000
export PROD_AM_HEAD_TRANSVERSE_K=100
export PROD_AM_TRANSVERSE_STRAIN_DETACHMENT_THETA=0
export PROD_TAG_SUFFIX=_split_k100_thetaperp0_kalign120_100ms

exec bash /home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp/production_qeq_coupled_A1280_M640_20260804.sbatch
