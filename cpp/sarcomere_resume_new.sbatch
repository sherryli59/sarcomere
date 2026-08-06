#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
# Match the workers' footprint (24 cpus, default DefMemPerCPU*24 ~= 96GB) so resume jobs PACK
# many-per-node instead of grabbing a whole node each. The old `--exclusive` capped concurrency
# at ~3 (149-task arrays sat PD "Priority" forever) and starved backlog recovery.
#SBATCH --cpus-per-task=24
#SBATCH --job-name=sarco_resume
#SBATCH --open-mode=append
#SBATCH --time=48:00:00

set -euo pipefail

# Arguments passed by the babysitter
RESUME_LIST=${1:?missing RESUME_LIST}
LOG_DIR=${2:?missing LOG_DIR}
TARGET_TOTAL_STEPS=${3:?missing TARGET_TOTAL_STEPS}
BACKUP_DIR=${4:?missing BACKUP_DIR}

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH

if [ -n "${SLURM_CPUS_PER_TASK:-}" ]; then
    export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
else
    export OMP_NUM_THREADS=16
fi

# The babysitter submits this job from INSIDE the babycron Slurm job (which sets --mem),
# so SLURM_MEM_PER_NODE leaks in via --export=ALL. This job's own --exclusive allocation
# then also sets SLURM_MEM_PER_CPU, and `srun` aborts: "SLURM_MEM_PER_CPU, SLURM_MEM_PER_GPU,
# and SLURM_MEM_PER_NODE are mutually exclusive." Drop all three; --exclusive already grants
# the whole node's memory, so srun needs no memory constraint.
unset SLURM_MEM_PER_NODE SLURM_MEM_PER_CPU SLURM_MEM_PER_GPU

mkdir -p "${LOG_DIR}" "${BACKUP_DIR}"

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "SLURM_ARRAY_TASK_ID is not set; this script expects a Slurm array task."
    exit 1
fi

# 1. Select the specific file for this array task from the text file.
OUTPUT_FILENAME=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${RESUME_LIST}")
if [ -z "${OUTPUT_FILENAME}" ]; then
    echo "No file listed for array task ${SLURM_ARRAY_TASK_ID}. Exiting."
    exit 0
fi
FILENAME=$(basename "${OUTPUT_FILENAME}" .h5)

# Clear this file's in-flight marker whenever the resume job exits (completion, lock-exit, or
# SIGTERM-driven preemption). Without this, an age-only marker keeps the babysitter from
# re-resuming a file for the full RESUME_COOLDOWN even though nothing is running on it.
trap 'rm -f "${BACKUP_DIR}/.inflight/${FILENAME}"' EXIT

echo "[${FILENAME}] Selected for resume task ${SLURM_ARRAY_TASK_ID}: ${OUTPUT_FILENAME}"

if [ ! -f "${OUTPUT_FILENAME}" ]; then
    echo "[${FILENAME}] File does not exist: ${OUTPUT_FILENAME}"
    exit 1
fi

# Prevent two array jobs or babysitter retries from writing the same HDF5 file.
LOCK_FILE="${OUTPUT_FILENAME}.lock"
exec 9>"${LOCK_FILE}"
if ! flock -n 9; then
    echo "[${FILENAME}] Another resume job already holds ${LOCK_FILE}. Exiting."
    exit 0
fi

# Make a copy before touching the HDF5 file.
TIMESTAMP=$(date +%s)
BACKUP_FILE="${BACKUP_DIR}/${FILENAME}_${TIMESTAMP}_job${SLURM_JOB_ID:-nojob}_array${SLURM_ARRAY_TASK_ID}.bak.h5"
echo "[${FILENAME}] Creating safety backup at: ${BACKUP_FILE}"
cp "${OUTPUT_FILENAME}" "${BACKUP_FILE}"

# If a previous crashed resume left the temporary checkpoint group around,
# remove it after backup. The simulator recreates /resume_new on each save.
python - "${OUTPUT_FILENAME}" <<'PY'
import sys
import h5py

fname = sys.argv[1]
with h5py.File(fname, "a") as f:
    if "/resume_new" in f:
        del f["/resume_new"]
        print("deleted stale /resume_new")
PY

# 2. Dynamically parse the parameters from the filename. Handles _nodir.
if [[ "${FILENAME}" =~ ^A([0-9]+)_M([0-9]+)_lt([^_]+)(_nodir)?_r([0-9]+)_L([0-9]+)x([0-9]+)x([0-9]+)$ ]]; then
    n_actins="${BASH_REMATCH[1]}"
    n_myosins="${BASH_REMATCH[2]}"
    lifetime_label="${BASH_REMATCH[3]}"
    nodir_flag="${BASH_REMATCH[4]}"
    rep_label="${BASH_REMATCH[5]}"
    Lx="${BASH_REMATCH[6]}"
    Ly="${BASH_REMATCH[7]}"
    Lz="${BASH_REMATCH[8]}"
else
    echo "[${FILENAME}] Does not match expected pattern. Skipping."
    exit 1
fi

# 3. Determine directional flag based on filename.
if [ "${nodir_flag}" = "_nodir" ]; then
    DIRECTIONAL_FLAG="--directional=false"
else
    DIRECTIONAL_FLAG="--directional=true"
fi

# 4. Convert lifetime_label back into lifetime_coeff.
if [ "${lifetime_label}" = "0" ]; then
    lifetime_coeff=0
else
    log10_val=${lifetime_label#1e}
    lifetime_coeff=$(python - <<EOF
import math
x = float("${log10_val}")
print(f"{10**x:.6g}")
EOF
)
fi

SAVE_EVERY=$(python - "${OUTPUT_FILENAME}" <<'PY'
import sys
import h5py

fname = sys.argv[1]
save_every = 100
try:
    with h5py.File(fname, "r") as f:
        if "/state/current_step" in f:
            steps = f["/state/current_step"][...].reshape(-1)
            diffs = [int(steps[i] - steps[i - 1]) for i in range(1, len(steps)) if steps[i] > steps[i - 1]]
            if diffs:
                save_every = diffs[-1]
except Exception:
    pass
print(save_every)
PY
)

if ! [[ "${SAVE_EVERY}" =~ ^[0-9]+$ ]] || [ "${SAVE_EVERY}" -le 0 ]; then
    echo "[${FILENAME}] Could not infer save_every. Falling back to 100."
    SAVE_EVERY=100
fi

# 5. Extract current step. Prefer /resume/current_step, then /state/current_step,
# then derive from saved frame count as (n_frames - 1) * save_every.
CURRENT_STEP_RAW=$(python - "${OUTPUT_FILENAME}" "${SAVE_EVERY}" <<'PY'
import sys
fname = sys.argv[1]
save_every = int(sys.argv[2])
step = -1
try:
    import h5py
    with h5py.File(fname, "r") as f:
        for p in ("/resume/current_step", "/resume_prev/current_step", "/state/current_step"):
            if p in f:
                data = f[p][...]
                if data.size > 0:
                    step = max(step, int(data.reshape(-1)[-1]))
        if step < 0:
            for cand in ("/actin/center", "actin/center", "/myosin/center", "myosin/center",
                         "/actin/centers", "actin/centers"):
                if cand in f:
                    d = f[cand]
                    if hasattr(d, "shape") and len(d.shape) >= 1 and d.shape[0] > 0:
                        step = (int(d.shape[0]) - 1) * save_every
                        break
except Exception:
    pass
print(step)
PY
)

CURRENT_STEP=$(echo "${CURRENT_STEP_RAW}" | tr -d '[:space:]')

if ! [[ "${CURRENT_STEP}" =~ ^-?[0-9]+$ ]] || [ "${CURRENT_STEP}" -lt 0 ]; then
    echo "[${FILENAME}] Could not parse current_step. Skipping."
    exit 0
fi

remaining=$(( TARGET_TOTAL_STEPS - CURRENT_STEP ))
if [ "${remaining}" -le 0 ]; then
    echo "[${FILENAME}] Complete at step ${CURRENT_STEP}. Exiting."
    exit 0
fi

run_steps=${remaining}
SEED=$(( n_actins + 1000 * (10#${rep_label}) ))

MYOSIN_RADIUS=0.058
AM_CUTOFF=0.08
AM_OPTIMAL=0.06
BASE_LIFETIME=0.0
DIRECTIONAL_BASE_LIFETIME=0
K_MM=0
K_BUNDLE_MAX=0.1
MAX_MYOSIN_BONDS=$(python - "${OUTPUT_FILENAME}" "${n_myosins}" <<'PY'
import sys
import h5py

fname = sys.argv[1]
n_myosins = int(sys.argv[2])
max_myosin_bonds = 6
try:
    with h5py.File(fname, "r") as f:
        if n_myosins > 0 and "/actin_myo/bonds" in f:
            shape = f["/actin_myo/bonds"].shape
            if len(shape) >= 2 and shape[1] % n_myosins == 0:
                inferred = shape[1] // n_myosins
                if inferred > 0:
                    max_myosin_bonds = inferred
except Exception:
    pass
print(max_myosin_bonds)
PY
)

if ! [[ "${MAX_MYOSIN_BONDS}" =~ ^[0-9]+$ ]] || [ "${MAX_MYOSIN_BONDS}" -le 0 ]; then
    echo "[${FILENAME}] Could not infer max_myosin_bonds. Falling back to 6."
    MAX_MYOSIN_BONDS=6
fi
PERIODIC_MASK=111

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -x "${SCRIPT_DIR}/build/sarcomere" ]; then
    SARCOMERE_EXE="${SCRIPT_DIR}/build/sarcomere"
else
    SARCOMERE_EXE="build/sarcomere"
fi

echo "[${FILENAME}] Resuming from step ${CURRENT_STEP} for ${run_steps} steps to target ${TARGET_TOTAL_STEPS}."
echo "[${FILENAME}] Using save_every=${SAVE_EVERY}, max_myosin_bonds=${MAX_MYOSIN_BONDS}."

cmd=(
    srun "${SARCOMERE_EXE}"
    --myosin_radius="${MYOSIN_RADIUS}"
    --am_cutoff="${AM_CUTOFF}"
    --am_optimal="${AM_OPTIMAL}"
    --n_actins="${n_actins}"
    --n_myosins="${n_myosins}"
    --filename="${OUTPUT_FILENAME}"
    --periodic_mask="${PERIODIC_MASK}"
    --Lx="${Lx}" --Ly="${Ly}" --Lz="${Lz}"
    --nsteps="${TARGET_TOTAL_STEPS}"
    --lifetime_coeff="${lifetime_coeff}"
    --save_every="${SAVE_EVERY}"
    --base_lifetime="${BASE_LIFETIME}"
    --directional_base_lifetime="${DIRECTIONAL_BASE_LIFETIME}"
    "${DIRECTIONAL_FLAG}"
    --k_mm="${K_MM}"
    --k_bundle_max="${K_BUNDLE_MAX}"
    --max_myosin_bonds="${MAX_MYOSIN_BONDS}"
    --bundle_ramp_steps=0
    --seed="${SEED}"
    --resume=true
)

printf '[%s] Command:' "${FILENAME}" >> "${LOG_DIR}/${FILENAME}.sim.out"
printf ' %q' "${cmd[@]}" >> "${LOG_DIR}/${FILENAME}.sim.out"
printf '\n' >> "${LOG_DIR}/${FILENAME}.sim.out"

"${cmd[@]}" >> "${LOG_DIR}/${FILENAME}.sim.out" 2>&1
