#!/bin/bash
#SBATCH -p owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32GB
#SBATCH --time=2-00:00:00
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/logs/%j_wb4lambda_resume.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs/%j_wb4lambda_resume.err

# Resume one A*_M*_lt*_r*_L12x4x4 trajectory from its HDF5 checkpoint.
# Submitted by wb4_lambda_babysitter.sh.

set -euo pipefail

OUTPUT_FILENAME="${1:?missing HDF5 path}"
CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
SARCOMERE_EXE="${CPP_DIR}/build_release/sarcomere"
LOG_DIR=/scratch/users/xsli/sarcomere/logs
TARGET_TOTAL_STEPS=300000
SAVE_EVERY=200

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}"
mkdir -p "${LOG_DIR}"

FILENAME=$(basename "${OUTPUT_FILENAME}" .h5)
scontrol update jobid="${SLURM_JOB_ID}" name="${FILENAME}" 2>/dev/null || true

if [ ! -x "${SARCOMERE_EXE}" ]; then
    echo "[${FILENAME}] missing executable: ${SARCOMERE_EXE}" >&2
    exit 1
fi
if [ ! -f "${OUTPUT_FILENAME}" ]; then
    echo "[${FILENAME}] missing trajectory" >&2
    exit 1
fi

# Hold the same sidecar lock used by the original array task.
exec 9>"${OUTPUT_FILENAME}.lock"
if ! flock -n 9; then
    echo "[${FILENAME}] another writer holds the lock; exiting"
    exit 0
fi

if [[ "${FILENAME}" =~ ^A([0-9]+)_M([0-9]+)_lt([^_]+)_r([0-9]+)_L([0-9]+)x([0-9]+)x([0-9]+)$ ]]; then
    n_actins="${BASH_REMATCH[1]}"
    n_myosins="${BASH_REMATCH[2]}"
    lifetime_label="${BASH_REMATCH[3]}"
    Lx="${BASH_REMATCH[5]}"
    Ly="${BASH_REMATCH[6]}"
    Lz="${BASH_REMATCH[7]}"
else
    echo "[${FILENAME}] does not match the lambda-sweep filename pattern" >&2
    exit 1
fi

lifetime_coeff=$(python - "${lifetime_label}" <<'PY'
import sys
label = sys.argv[1]
print("0" if label == "0" else f"{10 ** float(label.removeprefix('1e')):.8g}")
PY
)

# A killed checkpoint transaction can leave this temporary group behind.
python - "${OUTPUT_FILENAME}" <<'PY'
import sys
import h5py

with h5py.File(sys.argv[1], "a") as handle:
    if "/resume_new" in handle:
        del handle["/resume_new"]
        print("deleted stale /resume_new")
PY

CURRENT_STEP=$(python - "${OUTPUT_FILENAME}" "${SAVE_EVERY}" <<'PY'
import sys
import h5py

path, save_every = sys.argv[1], int(sys.argv[2])
with h5py.File(path, "r") as handle:
    if "/resume/current_step" in handle:
        step = int(handle["/resume/current_step"][...].reshape(-1)[-1])
    elif "/state/current_step" in handle:
        step = int(handle["/state/current_step"][...].reshape(-1)[-1])
    else:
        step = (int(handle["/actin/center"].shape[0]) - 1) * save_every
print(step)
PY
)
CURRENT_STEP=$(echo "${CURRENT_STEP}" | tr -d '[:space:]')
if ! [[ "${CURRENT_STEP}" =~ ^[0-9]+$ ]]; then
    echo "[${FILENAME}] invalid checkpoint step: ${CURRENT_STEP}" >&2
    exit 1
fi
if [ "${CURRENT_STEP}" -ge "${TARGET_TOTAL_STEPS}" ]; then
    echo "[${FILENAME}] already complete (${CURRENT_STEP})"
    exit 0
fi

echo "[${FILENAME}] resuming from ${CURRENT_STEP}/${TARGET_TOTAL_STEPS}"
srun --kill-on-bad-exit=1 "${SARCOMERE_EXE}" \
    --myosin_radius=0.058 \
    --am_cutoff=0.08 \
    --am_optimal=0.06 \
    --n_actins="${n_actins}" \
    --n_myosins="${n_myosins}" \
    --filename="${OUTPUT_FILENAME}" \
    --periodic_mask=111 \
    --Lx="${Lx}" --Ly="${Ly}" --Lz="${Lz}" \
    --nsteps="${TARGET_TOTAL_STEPS}" \
    --lifetime_coeff="${lifetime_coeff}" \
    --base_lifetime=0.0 \
    --directional_base_lifetime=0.0 \
    --directional=true \
    --save_every="${SAVE_EVERY}" \
    --k_mm=0 \
    --k_bundle_max=0.1 \
    --max_myosin_bonds=6 \
    --bundle_ramp_steps=5000 \
    --resume=true \
    >> "${LOG_DIR}/${FILENAME}.sim.out" 2>&1
