#!/bin/bash
# Final historical Cell A phenotype gate for the corrected k_mm=0 control.
#SBATCH --partition=owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH --time=01:00:00
#SBATCH --job-name=validate_cellA
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/cellA_sparse_myomesin_A800_M320_5s_dt10us_20260818/slurm/%x_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/cellA_sparse_myomesin_A800_M320_5s_dt10us_20260818/slurm/%x_%j.err

set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH

root=/scratch/users/xsli/sarcomere/cellA_sparse_myomesin_A800_M320_5s_dt10us_20260818
reference=/scratch/users/xsli/sarcomere/prod_brake_A800_M320_5s_20260807/data/A800_M320_brake_seed20260805_5s_fresh.h5
candidate=${root}/full5s_nobarrier/data/A800_M320_cellA_control_seed20260805_5000ms.h5
validator=${root}/source_snapshot_nobarrier/validate_cell_a_replication.py
output=${root}/full5s_nobarrier/analysis/cell_a_replication.json
csv=${root}/full5s_nobarrier/analysis/cell_a_replication_timeseries.csv

python "${validator}" "${reference}" "${candidate}" \
    --output "${output}" --csv "${csv}"
