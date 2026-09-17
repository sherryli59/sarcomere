#!/bin/bash
# Submit with an afterany dependency on the active six-task production array.
#SBATCH --partition=owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --time=00:30:00
#SBATCH --job-name=sarc_mm_300k_audit
#SBATCH --output=/scratch/users/xsli/sarcomere/prod_sparse_myomesin_module_A800_M320_750ms_dt2p5us_300k_20260818/logs/%j_audit.out
#SBATCH --error=/scratch/users/xsli/sarcomere/prod_sparse_myomesin_module_A800_M320_750ms_dt2p5us_300k_20260818/logs/%j_audit.err

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export HDF5_USE_FILE_LOCKING=FALSE

root=/scratch/users/xsli/sarcomere/prod_sparse_myomesin_module_A800_M320_750ms_dt2p5us_300k_20260818
auditor=${root}/source_snapshot/audit_production_300k.py
test -r "${auditor}" || { echo "missing auditor: ${auditor}" >&2; exit 1; }
python "${auditor}" "${root}"
