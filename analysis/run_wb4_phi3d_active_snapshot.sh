#!/bin/bash
# Run inside an existing allocation.  The user requested an active HDF5
# snapshot, so native HDF5 file locking is explicitly disabled.  Snapshot
# outputs are kept separate from final-at-300k trajectory parts.

set -euo pipefail

REPO_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere
MANIFEST="${REPO_DIR}/analysis/wb4_phi3d_manifest_20260717.txt"
OUTPUT_ROOT="${REPO_DIR}/analysis/wb4_phi3d_subsample_20260717"
OLD_PART_DIR="${OUTPUT_ROOT}/parts"
SNAPSHOT_PART_DIR="${OUTPUT_ROOT}/active_snapshot_parts"
FINAL_CSV="${OUTPUT_ROOT}/wb4_phi3d_stride10_active_snapshot.csv"
PLOT="${OUTPUT_ROOT}/wb4_active_phi3d_trajectories.png"
RUN_TAG="${WB4_ACTIVE_RUN_TAG:-$(hostname -s)}"
LOG="${OUTPUT_ROOT}/active_snapshot_runner_${RUN_TAG}.log"
TASK_IDS="${WB4_ACTIVE_TASK_IDS:-4 5 6 7 8 9}"

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export HDF5_USE_FILE_LOCKING=FALSE
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
mkdir -p "${SNAPSHOT_PART_DIR}"

echo "active snapshot started $(date) on $(hostname)" | tee -a "${LOG}"
echo "HDF5_USE_FILE_LOCKING=${HDF5_USE_FILE_LOCKING}" | tee -a "${LOG}"

for task_id in ${TASK_IDS}; do
    input=$(sed -n "$((task_id + 1))p" "${MANIFEST}")
    stem=$(basename "${input}" .h5)
    output="${SNAPSHOT_PART_DIR}/${stem}.csv"
    tmp_output="${output}.tmp.$$"
    echo "[$(date)] ${stem}" | tee -a "${LOG}"
    completed=0
    for attempt in 1 2 3; do
        if python -u "${REPO_DIR}/analysis/extract_sarcomeric_order.py" "${input}" \
            --output "${tmp_output}" \
            --stride 10 \
            --nx 5 \
            --box 12 4 4 \
            --yz-bins 4 4 \
            --save-every 200 \
            --phi3d-n-max-x 30 \
            --phi3d-n-max-perp 5 \
            --shared-all-nx \
            >> "${LOG}" 2>&1
        then
            mv "${tmp_output}" "${output}"
            completed=1
            break
        fi
        rm -f "${tmp_output}"
        echo "attempt ${attempt} failed for ${stem}; retrying in 30 s" | tee -a "${LOG}"
        sleep 30
    done
    if [ "${completed}" -ne 1 ]; then
        echo "analysis failed after three attempts: ${stem}" | tee -a "${LOG}"
        exit 1
    fi
done

# Multiple GPU-test allocations may split the six inputs.  The last runner to
# finish owns the merge lock and creates the combined table and plot.
snapshot_count=$(find "${SNAPSHOT_PART_DIR}" -maxdepth 1 -type f -name 'A8000_M4000_*.csv' | wc -l)
if [ "${snapshot_count}" -lt 6 ]; then
    echo "${snapshot_count}/6 active snapshot parts complete; another runner will merge" | tee -a "${LOG}"
    exit 0
fi

exec 8>"${OUTPUT_ROOT}/.active_snapshot_merge.lock"
flock 8

python - "${MANIFEST}" "${OLD_PART_DIR}" "${SNAPSHOT_PART_DIR}" "${FINAL_CSV}" <<'PY'
import csv
import sys
from pathlib import Path

manifest = Path(sys.argv[1])
old_parts = Path(sys.argv[2])
snapshot_parts = Path(sys.argv[3])
output = Path(sys.argv[4])
fieldnames = None
rows = []
for index, line in enumerate(manifest.read_text().splitlines()):
    stem = Path(line).stem
    part = (old_parts if index < 4 else snapshot_parts) / f"{stem}.csv"
    with part.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if fieldnames is None:
            fieldnames = reader.fieldnames
        elif reader.fieldnames != fieldnames:
            raise SystemExit(f"CSV header mismatch: {part}")
        rows.extend(reader)
with output.open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)
print(f"merged {len(rows)} rows into {output}")
PY

python "${REPO_DIR}/analysis/plot_wb4_phi3d_active.py" \
    "${FINAL_CSV}" --output "${PLOT}" | tee -a "${LOG}"
echo "active snapshot finished $(date)" | tee -a "${LOG}"
