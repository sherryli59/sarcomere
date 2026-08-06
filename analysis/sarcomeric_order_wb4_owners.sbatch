#!/bin/bash
#SBATCH --partition=owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=sarc_wb4_3d
#SBATCH --output=/scratch/users/xsli/sarcomere/logs/sarc_wb4_3d_%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs/sarc_wb4_3d_%A_%a.err

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH

# The simulations are active writers. Never request or create an HDF5 file lock.
export HDF5_USE_FILE_LOCKING=FALSE
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

REPO_DIR="/home/groups/rotskoff/sherry/sarcomere/sarcomere"
ANALYSIS_SCRIPT="${REPO_DIR}/analysis/extract_sarcomeric_order.py"
OUTPUT_ROOT="${OUTPUT_ROOT:?OUTPUT_ROOT must be exported by sbatch}"
PART_DIR="${OUTPUT_ROOT}/parts"
FINAL_CSV="${OUTPUT_ROOT}/wb4_sarcomeric_order_shared3d_stride10.csv"

mkdir -p "${PART_DIR}"

if [ "${MERGE_ONLY:-0}" = "1" ]; then
    python - "${PART_DIR}" "${FINAL_CSV}" <<'PY'
import csv
import sys
from pathlib import Path

part_dir = Path(sys.argv[1])
output = Path(sys.argv[2])
parts = sorted(part_dir.glob("*.csv"))
if not parts:
    raise SystemExit(f"No CSV parts found in {part_dir}")

fieldnames = None
rows = []
for part in parts:
    with part.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if fieldnames is None:
            fieldnames = reader.fieldnames
        elif reader.fieldnames != fieldnames:
            raise SystemExit(f"CSV header mismatch in {part}")
        rows.extend(reader)

with output.open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print(f"Merged {len(rows)} rows from {len(parts)} files into {output}")
PY
    exit 0
fi

DATA_DIRS=(
    "/scratch/users/xsli/sarcomere/data_wb4x4"
    "/scratch/users/xsli/sarcomere/data_wb4x4_lambda"
)

files=()
for data_dir in "${DATA_DIRS[@]}"; do
    while IFS= read -r path; do
        files+=("${path}")
    done < <(find "${data_dir}" -maxdepth 1 -type f -name '*.h5' | sort)
done

task_id=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
if [ "${task_id}" -ge "${#files[@]}" ]; then
    echo "Array index ${task_id} exceeds ${#files[@]} discovered HDF5 files" >&2
    exit 1
fi

input=${files[$task_id]}
dataset=$(basename "$(dirname "${input}")")
stem=$(basename "${input}" .h5)
output="${PART_DIR}/${dataset}__${stem}.csv"

echo "Input: ${input}"
echo "Output: ${output}"
echo "HDF5_USE_FILE_LOCKING=${HDF5_USE_FILE_LOCKING}"

for attempt in 1 2 3; do
    if python -u "${ANALYSIS_SCRIPT}" "${input}" \
        --output "${output}" \
        --stride 10 \
        --nx 5 \
        --box 12 4 4 \
        --yz-bins 4 4 \
        --save-every 200 \
        --phi3d-n-max-x 30 \
        --phi3d-n-max-perp 5
    then
        exit 0
    fi
    echo "Attempt ${attempt} failed while the simulation was writing; retrying" >&2
    sleep 30
done

echo "Analysis failed after three attempts: ${input}" >&2
exit 1
