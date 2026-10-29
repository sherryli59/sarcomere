#!/bin/bash
#SBATCH --job-name=status2readout
#SBATCH --partition=hns,normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=02:00:00
set -euo pipefail
export PATH="$GROUP_HOME/sherry/mamba_root/envs/sarcomere/bin:$PATH"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
exec "$GROUP_HOME/sherry/mamba_root/envs/sarcomere/bin/python" -u \
  "$GROUP_HOME/sherry/sarcomere/sarcomere/analysis/visualizations/gaphunt-rig/status2_gate_20261028/analyze_resume.py"
