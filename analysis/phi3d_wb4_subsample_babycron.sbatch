#!/bin/bash
#SBATCH --partition=owners,rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=00:10:00
#SBATCH --job-name=phi3d_wb4_cron
#SBATCH --output=/scratch/users/xsli/sarcomere/logs/phi3d_wb4_cron_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs/phi3d_wb4_cron_%j.err
#SBATCH --no-requeue

set -euo pipefail

REPO_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere
MANIFEST="${WB4_PHI3D_MANIFEST:-${REPO_DIR}/analysis/wb4_phi3d_manifest_20260717.txt}"
OUTPUT_ROOT="${WB4_PHI3D_OUTPUT_ROOT:-${REPO_DIR}/analysis/wb4_phi3d_subsample_20260717}"
PART_DIR="${OUTPUT_ROOT}/parts"
FINAL_CSV="${OUTPUT_ROOT}/wb4_phi3d_stride10.csv"
INTERVAL_MIN="${WB4_PHI3D_INTERVAL_MIN:-30}"

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
cd "${REPO_DIR}"

echo "4x4 phi_3D babycron scan $(date)"
bash analysis/phi3d_wb4_subsample_babysitter.sh || echo "analysis scan returned nonzero; chain will continue"

expected=$(wc -l < "${MANIFEST}")
complete=$(find "${PART_DIR}" -maxdepth 1 -type f -name '*.csv' 2>/dev/null | wc -l)
if [ "${complete}" -lt "${expected}" ]; then
    echo "${complete}/${expected} trajectory parts complete; next scan in ${INTERVAL_MIN} min"
    sbatch --begin="now+${INTERVAL_MIN}minutes" analysis/phi3d_wb4_subsample_babycron.sbatch
    exit 0
fi

python - "${MANIFEST}" "${PART_DIR}" "${FINAL_CSV}" <<'PY'
import csv
import sys
from pathlib import Path

manifest, part_dir, output = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
fieldnames = None
rows = []
for line in manifest.read_text().splitlines():
    stem = Path(line).stem
    part = part_dir / f"{stem}.csv"
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
print(f"merged {len(rows)} rows from {len(manifest.read_text().splitlines())} trajectories into {output}")
PY
echo "all ${expected} phi_3D trajectories complete; analysis chain stopping"
