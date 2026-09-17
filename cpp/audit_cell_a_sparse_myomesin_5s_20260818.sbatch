#!/bin/bash
# Submit with an afterany dependency on the two-task Cell A comparison array.
#SBATCH --partition=rotskoff
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --time=00:30:00
#SBATCH --job-name=sarc_cellA_mm_audit
#SBATCH --output=/scratch/users/xsli/sarcomere/cellA_sparse_myomesin_A800_M320_5s_dt10us_20260818/slurm/%j_audit.out
#SBATCH --error=/scratch/users/xsli/sarcomere/cellA_sparse_myomesin_A800_M320_5s_dt10us_20260818/slurm/%j_audit.err

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export HDF5_USE_FILE_LOCKING=FALSE

root=/scratch/users/xsli/sarcomere/cellA_sparse_myomesin_A800_M320_5s_dt10us_20260818/full5s
auditor=/scratch/users/xsli/sarcomere/cellA_sparse_myomesin_A800_M320_5s_dt10us_20260818/source_snapshot/audit_cell_a_sparse_myomesin.py
test -r "${auditor}" || { echo "missing auditor: ${auditor}" >&2; exit 1; }
python "${auditor}" "${root}" --dt 1e-5 --nsteps 500000 --save-every 500
