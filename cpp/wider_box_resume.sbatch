#!/bin/bash
#SBATCH -p rotskoff,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96GB
#SBATCH --time=2-00:00:00
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/logs/%j_wbresume.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs/%j_wbresume.err

# Dedicated resume worker for the wider-box registration run. Unlike the standard
# sarcomere_resume.sbatch, this parses the _mmb token and uses save_every=200 and
# the correct max_myosin_bonds -- the standard one would skip these files (regex)
# or mis-parameterize them (mmb, save_every). Submitted per-file by
# wider_box_babysitter.sh. Continues from the HDF5 /resume checkpoint.
#
# Arg: $1 = full path to the .h5 to resume.

set -uo pipefail
OUTPUT_FILENAME="${1:?missing H5 file}"
LOG_DIR="/scratch/users/xsli/sarcomere/logs"
TARGET_TOTAL_STEPS=400000
SAVE_EVERY=200

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
cd /home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
mkdir -p "${LOG_DIR}"

FILENAME=$(basename "${OUTPUT_FILENAME}" .h5)
scontrol update jobid="${SLURM_JOB_ID}" name="${FILENAME}" 2>/dev/null || true

if [ ! -f "${OUTPUT_FILENAME}" ]; then echo "[${FILENAME}] missing; abort"; exit 1; fi

# single-writer lock (held for the whole run) -- prevents double-writing the HDF5
exec 9>"${OUTPUT_FILENAME}.lock"
if ! flock -n 9; then echo "[${FILENAME}] another resume holds the lock; exit"; exit 0; fi

# parse params incl. the _mmb token (standard regex lacks this)
if [[ "${FILENAME}" =~ ^A([0-9]+)_M([0-9]+)_lt([^_]+)_mmb([0-9]+)_r([0-9]+)_L([0-9]+)x([0-9]+)x([0-9]+)$ ]]; then
    n_actins="${BASH_REMATCH[1]}"; n_myosins="${BASH_REMATCH[2]}"
    lifetime_label="${BASH_REMATCH[3]}"; MAX_MYOSIN_BONDS="${BASH_REMATCH[4]}"
    Lx="${BASH_REMATCH[6]}"; Ly="${BASH_REMATCH[7]}"; Lz="${BASH_REMATCH[8]}"
else
    echo "[${FILENAME}] does not match wider-box pattern; abort"; exit 1
fi
log10_val=${lifetime_label#1e}
lifetime_coeff=$(python - <<EOF
print(f"{10**float('${log10_val}'):.6g}")
EOF
)

# drop any stale temp checkpoint group from a crashed prior resume
python - "${OUTPUT_FILENAME}" <<'PY'
import sys, h5py
with h5py.File(sys.argv[1], "a") as f:
    if "/resume_new" in f:
        del f["/resume_new"]; print("deleted stale /resume_new")
PY

# current step = frames * save_every
CURRENT_STEP=$(python - "${OUTPUT_FILENAME}" "${SAVE_EVERY}" <<'PY'
import sys, h5py
fn, se = sys.argv[1], int(sys.argv[2]); step = -1
try:
    with h5py.File(fn, "r") as f:
        for c in ("/actin/center", "/myosin/center"):
            if c in f and len(f[c].shape) >= 1:
                step = int(f[c].shape[0]) * se; break
except Exception:
    pass
print(step)
PY
)
CURRENT_STEP=$(echo "${CURRENT_STEP}" | tr -d '[:space:]')
if ! [[ "${CURRENT_STEP}" =~ ^[0-9]+$ ]]; then echo "[${FILENAME}] bad current_step; abort"; exit 0; fi
if [ "${CURRENT_STEP}" -ge "${TARGET_TOTAL_STEPS}" ]; then echo "[${FILENAME}] complete (${CURRENT_STEP}); exit"; exit 0; fi

echo "[${FILENAME}] resume from ${CURRENT_STEP}/${TARGET_TOTAL_STEPS}, mmb=${MAX_MYOSIN_BONDS}"
srun build/sarcomere \
    --myosin_radius=0.058 --am_cutoff=0.08 --am_optimal=0.06 \
    --n_actins=${n_actins} --n_myosins=${n_myosins} \
    --filename="${OUTPUT_FILENAME}" --periodic_mask=111 \
    --Lx=${Lx} --Ly=${Ly} --Lz=${Lz} \
    --nsteps=${TARGET_TOTAL_STEPS} --lifetime_coeff=${lifetime_coeff} \
    --save_every=${SAVE_EVERY} --base_lifetime=0.0 --k_mm=0 \
    --k_bundle_max=0.1 --max_myosin_bonds=${MAX_MYOSIN_BONDS} \
    --bundle_ramp_steps=0 --resume=true \
    >> "${LOG_DIR}/${FILENAME}.sim.out" 2>&1
