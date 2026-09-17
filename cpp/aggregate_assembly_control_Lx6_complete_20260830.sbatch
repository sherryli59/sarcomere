#!/bin/bash
#SBATCH --partition=rotskoff,owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64GB
#SBATCH --time=08:00:00
#SBATCH --job-name=assembly_Lx6_complete_analysis
#SBATCH --output=/scratch/users/xsli/sarcomere/assembly_controls_A1600_M320_Lx6_20260828/logs/%x_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/assembly_controls_A1600_M320_Lx6_20260828/logs/%x_%j.err
set -euo pipefail
export HDF5_USE_FILE_LOCKING=FALSE
export OPENBLAS_NUM_THREADS=1
export MPLCONFIGDIR=/tmp/sarcomere-mpl-${SLURM_JOB_ID:-local}
mkdir -p "${MPLCONFIGDIR}"
repo=/home/groups/rotskoff/sherry/sarcomere/sarcomere
root=/scratch/users/xsli/sarcomere/assembly_controls_A1600_M320_Lx6_20260828

# Analyse whatever reached a full 5 s.  Deliberately NOT an equality check on a
# task count: the previous aggregate hard-required 21 files and was cancelled
# along with the turnover tasks, leaving the motor arms with no analysis.
mapfile -t complete < <(
  for f in "${root}"/data/A1600_M320_Lx6_*_5000ms.h5; do
    [ -e "$f" ] || continue
    python -c "
import h5py, sys, numpy as np
try:
    h = h5py.File(sys.argv[1], r, locking=False)
    n = int(h[meta/committed_frames][0])
    s = np.asarray(h[state/current_step][:n]).ravel()
    if s[-1] >= 500000: print(sys.argv[1])
    h.close()
except Exception: pass" "$f"
  done | sort
)
echo "complete 5 s trajectories: ${#complete[@]}"
printf  %sn "${complete[@]}"
if [[ ${#complete[@]} -lt 2 ]]; then
    echo "need at least 2 complete trajectories, found ${#complete[@]}" >&2
    exit 1
fi
python "${repo}/analysis/analyze_assembly_control_pilot.py" \
    "${complete[@]}" --output="${root}/analysis/assembly_summary_complete.csv"
python "${repo}/analysis/analyze_actin_turnover.py" \
    "${complete[@]}" --output="${root}/analysis/turnover_audit_complete.csv"
python "${repo}/analysis/plot_assembly_control_pilot.py" \
    "${root}/analysis/assembly_summary_complete.csv" \
    --output="${root}/analysis/assembly_response_complete.png" \
    --title="A1600/M320 Lx=6, 5 s - all complete trajectories"
