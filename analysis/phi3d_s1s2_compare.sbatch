#!/bin/bash
#SBATCH -p owners,normal,hns
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --time=01:00:00
#SBATCH --job-name=phi3d_s1s2
#SBATCH --open-mode=append
#SBATCH --output=/scratch/users/xsli/sarcomere/logs/%j_phi3d_s1s2.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs/%j_phi3d_s1s2.err

# Legacy myosin Phi_3D for the in-flight s1s2 sweep vs the completed s2only
# sweep, evaluated over a MATCHED absolute frame window so the comparison is
# not contaminated by the age gap between the two campaigns.

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH

# Both sweeps have live writers. Never take an HDF5 lock.
export HDF5_USE_FILE_LOCKING=FALSE
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

cd /home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis

python compare_phi3d_s1s2_vs_s2only.py \
    --output phi3d_s1s2_vs_s2only_matched.csv \
    --window 300 \
    --stride 10 \
    --min-frames 120 \
    --device numpy
