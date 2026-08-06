#!/bin/bash
#SBATCH -p rotskoff,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=48GB
#SBATCH --time=2-00:00:00
#SBATCH --job-name=wb8_lambda
#SBATCH --open-mode=append
#SBATCH --requeue
#SBATCH --error=/scratch/users/xsli/sarcomere/logs/%A_%a.err
#SBATCH --output=/scratch/users/xsli/sarcomere/logs/%A_%a.out
#SBATCH --array=0-12

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}"

CPP_DIR="/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp"
SARCOMERE_EXE="${CPP_DIR}/build_release/sarcomere"
SCRATCH_DIR="/scratch/users/xsli/sarcomere"
OUTPUT_DIR="${WB8_OUTPUT_DIR:-${SCRATCH_DIR}/data_wb8x8_lambda}"
LOG_DIR="${WB8_LOG_DIR:-${SCRATCH_DIR}/logs}"
PARAM_DIR="${WB8_PARAM_DIR:-${SCRATCH_DIR}/params_wb8x8_lambda}"
mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}" "${PARAM_DIR}"

if [ ! -x "${SARCOMERE_EXE}" ]; then
    echo "Missing optimized executable: ${SARCOMERE_EXE}" >&2
    exit 1
fi

# Match the established axial-x convention. The transverse 8x8 cross-section
# is 16 times the area of the 2x2 reference box.
Lx=12
Ly=8
Lz=8
n_actins=32000
n_myosins=16000
nsteps="${WB8_NSTEPS:-300000}"

# lambda = lifetime_coeff. Match the existing 13-point sweep: zero, then
# half-decade spacing from 10^-5.5 through 10^0.
LIFETIME_VALUES=(0 3.16228e-6 1e-5 3.16228e-5 1e-4 3.16228e-4 1e-3 3.16228e-3 1e-2 3.16228e-2 1e-1 3.16228e-1 1)
LIFETIME_LABELS=(0 1e-5.5 1e-5 1e-4.5 1e-4 1e-3.5 1e-3 1e-2.5 1e-2 1e-1.5 1e-1 1e-0.5 1e0)

task_id=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
if [ "${task_id}" -lt 0 ] || [ "${task_id}" -ge "${#LIFETIME_VALUES[@]}" ]; then
    echo "Invalid array task ${task_id}" >&2
    exit 1
fi
lifetime_coeff=${LIFETIME_VALUES[$task_id]}
lifetime_label=${LIFETIME_LABELS[$task_id]}

rep_label="r00"
SEED=$(( 98000 + task_id ))
MAX_MYOSIN_BONDS=6
MYOSIN_RADIUS=0.058
AM_CUTOFF=0.08
AM_OPTIMAL=0.06
SAVE_EVERY=200
BASE_LIFETIME=0.0
DIRECTIONAL_BASE_LIFETIME=0.0
K_MM=0
K_BUNDLE_MAX=0.1
BUNDLE_RAMP_STEPS=5000
PERIODIC_MASK=111
PRE_PHASE_STEPS="${WB8_PRE_PHASE_STEPS:-500}"

FILENAME="A${n_actins}_M${n_myosins}_lt${lifetime_label}_${rep_label}_L${Lx}x${Ly}x${Lz}"
OUTPUT_FILENAME="${OUTPUT_DIR}/${FILENAME}.h5"
SIM_LOG="${LOG_DIR}/${FILENAME}.sim.out"
PARAM_FILE="${PARAM_DIR}/${FILENAME}.params"

if [ "${WB8_SKIP_JOB_RENAME:-0}" != "1" ]; then
    scontrol update jobid="${SLURM_JOB_ID}" name="${FILENAME}" 2>/dev/null || true
fi

# Prevent a requeued task or manual retry from writing the same HDF5 file.
exec 9>"${OUTPUT_FILENAME}.lock"
if ! flock -n 9; then
    echo "[${FILENAME}] another process holds the output lock; exiting"
    exit 0
fi

cat > "${PARAM_FILE}" <<EOF
n_actins=${n_actins}
n_myosins=${n_myosins}
lifetime_coeff=${lifetime_coeff}
lifetime_label=${lifetime_label}
Lx=${Lx}
Ly=${Ly}
Lz=${Lz}
nsteps=${nsteps}
save_every=${SAVE_EVERY}
seed=${SEED}
max_myosin_bonds=${MAX_MYOSIN_BONDS}
executable=${SARCOMERE_EXE}
slurm_array_job_id=${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}
slurm_array_task_id=${task_id}
EOF

cmd=(
    srun --kill-on-bad-exit=1 "${SARCOMERE_EXE}"
    --myosin_radius="${MYOSIN_RADIUS}"
    --am_cutoff="${AM_CUTOFF}"
    --am_optimal="${AM_OPTIMAL}"
    --n_actins="${n_actins}"
    --n_myosins="${n_myosins}"
    --filename="${OUTPUT_FILENAME}"
    --periodic_mask="${PERIODIC_MASK}"
    --Lx="${Lx}" --Ly="${Ly}" --Lz="${Lz}"
    --nsteps="${nsteps}"
    --lifetime_coeff="${lifetime_coeff}"
    --base_lifetime="${BASE_LIFETIME}"
    --directional_base_lifetime="${DIRECTIONAL_BASE_LIFETIME}"
    --directional=true
    --initial_structure=myosin_x_noise
    --pre_phase_steps="${PRE_PHASE_STEPS}"
    --save_every="${SAVE_EVERY}"
    --k_mm="${K_MM}"
    --k_bundle_max="${K_BUNDLE_MAX}"
    --max_myosin_bonds="${MAX_MYOSIN_BONDS}"
    --bundle_ramp_steps="${BUNDLE_RAMP_STEPS}"
    --seed="${SEED}"
)

resume_ready=0
if [ -s "${OUTPUT_FILENAME}" ] && python - "${OUTPUT_FILENAME}" "${PRE_PHASE_STEPS}" <<'PY'
import sys
import h5py

minimum_step = int(sys.argv[2])
with h5py.File(sys.argv[1], "r") as handle:
    checkpoint_step = (
        int(handle["/resume/current_step"][...].reshape(-1)[-1])
        if "/resume/current_step" in handle
        else -1
    )
    valid = (
        "/resume/current_step" in handle
        and "/state/current_step" in handle
        and "/actin/center" in handle
        and handle["/actin/center"].shape[0] > 0
        and checkpoint_step >= minimum_step
    )
raise SystemExit(0 if valid else 1)
PY
then
    resume_ready=1
elif [ -s "${OUTPUT_FILENAME}" ]; then
    INCOMPLETE_FILE="${OUTPUT_FILENAME}.incomplete_job${SLURM_JOB_ID}_task${task_id}"
    echo "[${FILENAME}] preserving incomplete non-resumable output as ${INCOMPLETE_FILE}" | tee -a "${SIM_LOG}"
    mv "${OUTPUT_FILENAME}" "${INCOMPLETE_FILE}"
fi

if [ "${resume_ready}" -eq 1 ]; then
    echo "[${FILENAME}] existing output found; resuming to ${nsteps} steps" | tee -a "${SIM_LOG}"
    cmd+=(--resume=true)
else
    echo "[${FILENAME}] starting fresh ${nsteps}-step simulation" | tee -a "${SIM_LOG}"
fi

printf '[%s] Command:' "${FILENAME}" >> "${SIM_LOG}"
printf ' %q' "${cmd[@]}" >> "${SIM_LOG}"
printf '\n' >> "${SIM_LOG}"
"${cmd[@]}" >> "${SIM_LOG}" 2>&1
