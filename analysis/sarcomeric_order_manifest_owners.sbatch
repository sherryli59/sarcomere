#!/bin/bash
#SBATCH --partition=owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --job-name=sarc_2x2_3d
#SBATCH --output=/scratch/users/xsli/sarcomere/logs/sarc_2x2_3d_%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs/sarc_2x2_3d_%A_%a.err

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH

export HDF5_USE_FILE_LOCKING=FALSE
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"

REPO_DIR="/home/groups/rotskoff/sherry/sarcomere/sarcomere"
ANALYSIS_SCRIPT="${REPO_DIR}/analysis/extract_sarcomeric_order.py"
PLOT_SCRIPT="${REPO_DIR}/analysis/plot_sarcomeric_phase_diagrams.py"
MANIFEST="${MANIFEST:?MANIFEST must be exported by sbatch}"
OUTPUT_ROOT="${OUTPUT_ROOT:?OUTPUT_ROOT must be exported by sbatch}"
PART_DIR="${OUTPUT_ROOT}/parts"
FINAL_CSV="${OUTPUT_ROOT}/sarcomeric_order_shared3d_stride10.csv"
PHASE_DIR="${OUTPUT_ROOT}/phase_diagrams"

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

print(f"Merged {len(rows)} rows from {len(parts)} trajectories into {output}")
PY

    python -u "${PLOT_SCRIPT}" "${FINAL_CSV}" \
        --output-dir "${PHASE_DIR}" \
        --tail-samples 20 \
        --default-mmb 6
    exit 0
fi

task_id=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
mapfile -t inputs < <(
    python - "${MANIFEST}" "${task_id}" <<'PY'
import re
import sys
from collections import defaultdict
from pathlib import Path

manifest = Path(sys.argv[1])
task_id = int(sys.argv[2])
groups = defaultdict(list)
pattern = re.compile(
    r"A(\d+)_M(\d+)_lt([^_]+)_r\d+_L([^_]+)\.h5$"
)

def lambda_value(label):
    if "e" in label:
        mantissa, exponent = label.split("e", 1)
        return float(mantissa) * 10.0 ** float(exponent)
    return float(label)

for line in manifest.read_text().splitlines():
    path = Path(line.strip())
    if not str(path):
        continue
    match = pattern.search(path.name)
    if match is None:
        raise SystemExit(f"Cannot parse parameters from {path}")
    A, M, label, box = match.groups()
    key = (int(A), int(M), lambda_value(label), label, box)
    groups[key].append(path)

keys = sorted(groups)
if task_id >= len(keys):
    raise SystemExit(f"Task {task_id} exceeds {len(keys)} parameter groups")
for path in sorted(groups[keys[task_id]]):
    print(path)
PY
)

if [ "${#inputs[@]}" -eq 0 ]; then
    echo "No trajectories found for parameter-group task ${task_id}" >&2
    exit 1
fi

echo "Parameter-group task ${task_id}: ${#inputs[@]} replicas"
echo "HDF5_USE_FILE_LOCKING=${HDF5_USE_FILE_LOCKING}"

for input in "${inputs[@]}"; do
    stem=$(basename "${input}" .h5)
    output="${PART_DIR}/${stem}.csv"
    echo "Input: ${input}"
    echo "Output: ${output}"

    completed=0
    for attempt in 1 2 3; do
        if python -u "${ANALYSIS_SCRIPT}" "${input}" \
            --output "${output}" \
            --stride 10 \
            --nx 5 \
            --box 12 2 2 \
            --yz-bins 2 2 \
            --save-every 200 \
            --phi3d-n-max-x 30 \
            --phi3d-n-max-perp 5
        then
            completed=1
            break
        fi
        echo "Attempt ${attempt} failed for ${input}; retrying" >&2
        sleep 30
    done
    if [ "${completed}" -ne 1 ]; then
        echo "Analysis failed after three attempts: ${input}" >&2
        exit 1
    fi
done
