#!/bin/bash
#SBATCH -p owners,normal,hns
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=02:00:00
#SBATCH --job-name=s1_vs_lambda
#SBATCH --open-mode=append
#SBATCH --error=/scratch/users/xsli/sarcomere/logs/%j_s1_lambda.err
#SBATCH --output=/scratch/users/xsli/sarcomere/logs/%j_s1_lambda.out

# Count status-1 catch-bond engagement vs lambda (lt) in the OLD (data0430)
# sweep, which ran with the original kinetics where lambda*f_load stabilized
# status-1 bonds too. Only /actin/cb_status carries status-1 information --
# /actin/bonds is status-2-filtered at write time.

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH

# The sweep directory can contain active writers. Never take an HDF5 lock.
export HDF5_USE_FILE_LOCKING=FALSE
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

cd /home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis

python status1_vs_lambda.py \
    --data_dir /scratch/users/xsli/sarcomere/data0430 \
    --glob '*_L12x2x2.h5' \
    --stride 1 \
    --tail_frac 0.25 \
    --workers "${SLURM_CPUS_PER_TASK:-8}" \
    --prefix status1_vs_lambda_old0430
