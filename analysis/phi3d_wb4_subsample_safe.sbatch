#!/bin/bash
#SBATCH --partition=owners,rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=phi3d_wb4_s10
#SBATCH --output=/scratch/users/xsli/sarcomere/logs/phi3d_wb4_s10_%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs/phi3d_wb4_s10_%A_%a.err
#SBATCH --array=0-9
#SBATCH --no-requeue

# Compute a stride-10 phi_3D trajectory without racing a live HDF5 writer.
# The four inactive A6000 files are analyzed at their current terminal step.
# A8000 lambda files are accepted only after their checkpoint reaches 300000.

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export HDF5_USE_FILE_LOCKING=FALSE
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

REPO_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere
MANIFEST="${WB4_PHI3D_MANIFEST:-${REPO_DIR}/analysis/wb4_phi3d_manifest_20260717.txt}"
OUTPUT_ROOT="${WB4_PHI3D_OUTPUT_ROOT:-${REPO_DIR}/analysis/wb4_phi3d_subsample_20260717}"
PART_DIR="${OUTPUT_ROOT}/parts"
ANALYSIS_SCRIPT="${REPO_DIR}/analysis/extract_sarcomeric_order.py"
mkdir -p "${PART_DIR}"

task_id=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
input=$(sed -n "$((task_id + 1))p" "${MANIFEST}")
if [ -z "${input}" ] || [ ! -f "${input}" ]; then
    echo "manifest task ${task_id} has no readable HDF5 input" >&2
    exit 1
fi

stem=$(basename "${input}" .h5)
output="${PART_DIR}/${stem}.csv"
if [ -s "${output}" ]; then
    echo "[${stem}] output already exists: ${output}"
    exit 0
fi

# This is the same sidecar lock held for the full lifetime of every simulation
# worker.  Never inspect the HDF5 file unless the lock is ours.
exec 9>"${input}.lock"
if ! flock -n 9; then
    echo "[${stem}] simulation writer is live; deferring analysis"
    exit 0
fi

if [[ "${input}" == */data_wb4x4_lambda/* ]]; then
    if ! python - "${input}" <<'PY'
import sys
import h5py

with h5py.File(sys.argv[1], "r") as handle:
    if "/resume/current_step" in handle:
        step = int(handle["/resume/current_step"][...].reshape(-1)[-1])
    elif "/state/current_step" in handle:
        step = int(handle["/state/current_step"][...].reshape(-1)[-1])
    else:
        step = -1
print(f"checkpoint={step}")
raise SystemExit(0 if step >= 300000 else 1)
PY
    then
        echo "[${stem}] lambda trajectory is incomplete; deferring analysis"
        exit 0
    fi
fi

tmp_output="${output}.tmp.${SLURM_JOB_ID}"
trap 'rm -f "${tmp_output}"' EXIT
echo "[${stem}] computing stride-10 phi_3D trajectory"
python -u "${ANALYSIS_SCRIPT}" "${input}" \
    --output "${tmp_output}" \
    --stride 10 \
    --nx 5 \
    --box 12 4 4 \
    --yz-bins 4 4 \
    --save-every 200 \
    --phi3d-n-max-x 30 \
    --phi3d-n-max-perp 5 \
    --shared-all-nx
mv "${tmp_output}" "${output}"
trap - EXIT
echo "[${stem}] wrote ${output}"
