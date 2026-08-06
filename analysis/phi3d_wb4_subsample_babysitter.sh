#!/bin/bash
# Submit only manifest entries that are safe and ready for phi_3D analysis.

set -euo pipefail

REPO_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere
MANIFEST="${WB4_PHI3D_MANIFEST:-${REPO_DIR}/analysis/wb4_phi3d_manifest_20260717.txt}"
OUTPUT_ROOT="${WB4_PHI3D_OUTPUT_ROOT:-${REPO_DIR}/analysis/wb4_phi3d_subsample_20260717}"
PART_DIR="${OUTPUT_ROOT}/parts"
WORKER="${REPO_DIR}/analysis/phi3d_wb4_subsample_safe.sbatch"
LOCK_FILE=/scratch/users/xsli/sarcomere/.phi3d_wb4_subsample_babysitter.lock
mkdir -p "${PART_DIR}"

exec 8>"${LOCK_FILE}"
if ! flock -n 8; then
    echo "another 4x4 phi_3D scan is active; exiting"
    exit 0
fi
if squeue --me -h -n phi3d_wb4_s10 -o '%i' 2>/dev/null | grep -q .; then
    echo "4x4 phi_3D workers remain in the queue; deferring submission"
    exit 0
fi

ready=()
task_id=0
while IFS= read -r input; do
    stem=$(basename "${input}" .h5)
    output="${PART_DIR}/${stem}.csv"
    if [ -s "${output}" ]; then
        task_id=$((task_id + 1))
        continue
    fi

    is_ready=0
    (
        exec 9>"${input}.lock"
        flock -n 9 || exit 1
        if [[ "${input}" == */data_wb4x4_lambda/* ]]; then
            python - "${input}" <<'PY'
import sys
import h5py
with h5py.File(sys.argv[1], "r") as handle:
    if "/resume/current_step" in handle:
        step = int(handle["/resume/current_step"][...].reshape(-1)[-1])
    elif "/state/current_step" in handle:
        step = int(handle["/state/current_step"][...].reshape(-1)[-1])
    else:
        step = -1
raise SystemExit(0 if step >= 300000 else 1)
PY
        fi
    ) && is_ready=1
    if [ "${is_ready}" -eq 1 ]; then
        ready+=("${task_id}")
    fi
    task_id=$((task_id + 1))
done < "${MANIFEST}"

if [ "${#ready[@]}" -eq 0 ]; then
    echo "no new 4x4 trajectories are ready for phi_3D analysis"
    exit 0
fi

array_spec=$(IFS=,; echo "${ready[*]}")
echo "submitting ready phi_3D task IDs ${array_spec}"
if [ "${WB4_PHI3D_DRY_RUN:-0}" = 1 ]; then
    echo "dry run: sbatch --array=${array_spec} ${WORKER}"
    exit 0
fi
sbatch --array="${array_spec}" "${WORKER}"
