#!/bin/bash
# Recompute completed paired summaries with the persistent-Z replication gate.
#SBATCH --partition=owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH --time=01:00:00
#SBATCH --job-name=audit_cellA_mm_v2
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/cellA_sparse_myomesin_A800_M320_5s_dt10us_20260818/slurm/%x_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/cellA_sparse_myomesin_A800_M320_5s_dt10us_20260818/slurm/%x_%j.err

set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export HDF5_USE_FILE_LOCKING=FALSE

base=/scratch/users/xsli/sarcomere/cellA_sparse_myomesin_A800_M320_5s_dt10us_20260818
root=${base}/full5s
snapshot=${base}/source_snapshot_nobarrier
analyzer=${snapshot}/check_production_myomesin_pilot.py
auditor=${snapshot}/audit_cell_a_sparse_myomesin.py

for arm in control module40; do
    tag=A800_M320_cellA_${arm}_seed20260805_5000ms
    python "${analyzer}" "${root}/data/${tag}.h5" \
        --output "${root}/analysis/${tag}.json"
done
python "${auditor}" "${root}" --dt 1e-5 --nsteps 500000 \
    --save-every 500 --output "${root}/analysis/cellA_sparse_myomesin_audit_v2.json"
