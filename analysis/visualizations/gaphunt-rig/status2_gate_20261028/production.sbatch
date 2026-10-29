#!/bin/bash
#SBATCH --job-name=status2g130
#SBATCH --partition=hns,normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
set -euo pipefail
export PATH="$GROUP_HOME/sherry/mamba_root/envs/sarcomere/bin:$PATH"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export HDF5_USE_FILE_LOCKING=FALSE
exec "$GROUP_HOME/sherry/mamba_root/envs/sarcomere/bin/python" -u \
  "$GROUP_HOME/sherry/sarcomere/sarcomere/analysis/visualizations/gaphunt-rig/status2_gate_20261028/run_resume.py"
