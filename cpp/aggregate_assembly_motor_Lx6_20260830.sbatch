#!/bin/bash
#SBATCH --partition=rotskoff,owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --job-name=assembly_motor_Lx6_analysis
#SBATCH --output=/scratch/users/xsli/sarcomere/assembly_controls_A1600_M320_Lx6_20260828/logs/%x_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/assembly_controls_A1600_M320_Lx6_20260828/logs/%x_%j.err
set -euo pipefail
export HDF5_USE_FILE_LOCKING=FALSE
export OPENBLAS_NUM_THREADS=1
export MPLCONFIGDIR=/tmp/sarcomere-mpl-${SLURM_JOB_ID:-local}
mkdir -p "${MPLCONFIGDIR}"
repo=/home/groups/rotskoff/sherry/sarcomere/sarcomere
root=/scratch/users/xsli/sarcomere/assembly_controls_A1600_M320_Lx6_20260828

# The cancelled turnover tasks left partial trajectories on disk.  Only accept
# files that actually reached the full 500k steps.
complete=()
for f in "${root}"/data/A1600_M320_Lx6_*_5000ms.h5; do
    [[ -e "${f}" ]] || continue
    ok=$(python - "${f}" <<PY
import h5py, sys, numpy as np
try:
    with h5py.File(sys.argv[1], "r", locking=False) as h:
        n = int(h["meta/committed_frames"][0])
        s = np.asarray(h["state/current_step"][:n]).ravel()
        print(1 if len(s) and s[-1] >= 500000 else 0)
except Exception:
    print(0)
PY
)
    [[ "${ok}" == "1" ]] && complete+=("${f}")
done
echo "complete trajectories: ${#complete[@]}"
printf  %sn "${complete[@]##*/}"
if [[ ${#complete[@]} -lt 9 ]]; then
    echo "expected at least 9 complete trajectories, found ${#complete[@]}" >&2
    exit 1
fi
python "${repo}/analysis/analyze_assembly_control_pilot.py" \
    "${complete[@]}" --output="${root}/analysis/motor_summary.csv"
python "${repo}/analysis/analyze_actin_turnover.py" \
    "${complete[@]}" --output="${root}/analysis/motor_turnover_audit.csv"
python "${repo}/analysis/plot_assembly_control_pilot.py" \
    "${root}/analysis/motor_summary.csv" \
    --output="${root}/analysis/motor_response.png" \
    --title="A1600/M320 Lx=6 motor-availability campaign (5 s)"
