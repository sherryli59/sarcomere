#!/bin/bash
#SBATCH --partition=rotskoff
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --job-name=assembly_control_analysis
#SBATCH --output=/scratch/users/xsli/sarcomere/assembly_controls_A800_M160_20260828/logs/%x_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/assembly_controls_A800_M160_20260828/logs/%x_%j.err

set -euo pipefail
export HDF5_USE_FILE_LOCKING=FALSE
export OPENBLAS_NUM_THREADS=1
export MPLCONFIGDIR=/tmp/sarcomere-mpl-${SLURM_JOB_ID:-local}
mkdir -p "${MPLCONFIGDIR}"
repo=/home/groups/rotskoff/sherry/sarcomere/sarcomere
root=/scratch/users/xsli/sarcomere/assembly_controls_A800_M160_20260828
mapfile -t trajectories < <(find "${root}/data" -maxdepth 1 -name 'A800_M160_*_1000ms.h5' -type f | sort)
if [[ ${#trajectories[@]} -ne 21 ]]; then
    echo "expected 21 trajectories, found ${#trajectories[@]}" >&2
    exit 1
fi
python "${repo}/analysis/analyze_assembly_control_pilot.py" \
    "${trajectories[@]}" --output="${root}/analysis/assembly_summary.csv"
python "${repo}/analysis/analyze_actin_turnover.py" \
    "${trajectories[@]}" --output="${root}/analysis/turnover_audit.csv"
python "${repo}/analysis/plot_assembly_control_pilot.py" \
    "${root}/analysis/assembly_summary.csv" \
    --output="${root}/analysis/assembly_response.png"

