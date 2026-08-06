#!/bin/bash
#SBATCH --partition=normal
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=4GB
#SBATCH --time=01:00:00
#SBATCH --job-name=sarc_omp_diag
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/perf_openmp_diagnosis_20260723/slurm_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/perf_openmp_diagnosis_20260723/slurm_%j.err

set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
bash /home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/benchmark_openmp_scaling_diagnosis.sh
