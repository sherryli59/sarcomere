#!/bin/bash
#SBATCH --partition=owners,rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=48GB
#SBATCH --time=2-00:00:00
#SBATCH --job-name=a2k_s2_lambda
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/logs_a2000_status2_lambda/%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs_a2000_status2_lambda/%A_%a.err
#SBATCH --array=0-25

# Catch-bond sweep in the 12x2x2 box. The default is status-2-only load
# stabilization; A2K_STABILIZE_STATUS1=true restores the original status-1 +
# status-2 load stabilization law.
#
# Array layout: 13 lambda values x 2 paired replicas.  Seeds depend only on the
# replica, so r00 (and separately r01) starts from the same stochastic initial
# condition at every lambda. Status 1 has the base lifetime (zero in this
# experiment) unless A2K_STABILIZE_STATUS1 is enabled.
#
# This worker is deliberately --no-requeue.  The matching self-chaining Slurm
# watchdog resubmits only missing/incomplete task IDs and the HDF5 checkpoint is
# used to continue an interrupted trajectory.

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}"

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
SARCOMERE_EXE="${CPP_DIR}/build_release/sarcomere"
SCRATCH_DIR=/scratch/users/xsli/sarcomere
OUTPUT_DIR="${A2K_OUTPUT_DIR:-${SCRATCH_DIR}/data_a2000_status2_lambda}"
LOG_DIR="${A2K_LOG_DIR:-${SCRATCH_DIR}/logs_a2000_status2_lambda}"
PARAM_DIR="${A2K_PARAM_DIR:-${SCRATCH_DIR}/params_a2000_status2_lambda}"
KINETICS_TAG="${A2K_KINETICS_TAG:-s2only}"
STABILIZE_STATUS1="${A2K_STABILIZE_STATUS1:-false}"
mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}" "${PARAM_DIR}"

case "${STABILIZE_STATUS1}" in
    true|false) ;;
    *)
        echo "A2K_STABILIZE_STATUS1 must be true or false" >&2
        exit 1
        ;;
esac

if [ ! -x "${SARCOMERE_EXE}" ]; then
    echo "missing optimized executable: ${SARCOMERE_EXE}" >&2
    exit 1
fi

Lx=12
Ly=2
Lz=2
n_actins=2000
n_myosins=1000
nsteps="${A2K_NSTEPS:-300000}"

# Zero plus half-decades from 10^-5.5 through 10^0.
LIFETIME_VALUES=(0 3.16228e-6 1e-5 3.16228e-5 1e-4 3.16228e-4 1e-3 3.16228e-3 1e-2 3.16228e-2 1e-1 3.16228e-1 1)
LIFETIME_LABELS=(0 1e-5.5 1e-5 1e-4.5 1e-4 1e-3.5 1e-3 1e-2.5 1e-2 1e-1.5 1e-1 1e-0.5 1e0)
N_LAMBDAS=${#LIFETIME_VALUES[@]}
N_REPLICAS=2

task_id=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
if [ "${task_id}" -lt 0 ] || [ "${task_id}" -ge "$((N_LAMBDAS * N_REPLICAS))" ]; then
    echo "invalid array task ${task_id}" >&2
    exit 1
fi

lambda_idx=$((task_id % N_LAMBDAS))
rep_idx=$((task_id / N_LAMBDAS))
lifetime_coeff=${LIFETIME_VALUES[$lambda_idx]}
lifetime_label=${LIFETIME_LABELS[$lambda_idx]}
printf -v rep_label 'r%02d' "${rep_idx}"

# Pair initial conditions across lambda within each replica.
SEED=$((97000 + rep_idx))
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
PRE_PHASE_STEPS="${A2K_PRE_PHASE_STEPS:-500}"

FILENAME="A${n_actins}_M${n_myosins}_${KINETICS_TAG}_lt${lifetime_label}_${rep_label}_L${Lx}x${Ly}x${Lz}"
OUTPUT_FILENAME="${OUTPUT_DIR}/${FILENAME}.h5"
SIM_LOG="${LOG_DIR}/${FILENAME}.sim.out"
PARAM_FILE="${PARAM_DIR}/${FILENAME}.params"

if [ "${A2K_SKIP_JOB_RENAME:-0}" != 1 ]; then
    scontrol update jobid="${SLURM_JOB_ID}" name="${FILENAME}" 2>/dev/null || true
fi

# Prevent a watchdog retry and an original task from writing the same HDF5 file.
exec 9>"${OUTPUT_FILENAME}.lock"
if ! flock -n 9; then
    echo "[${FILENAME}] another process holds the output lock; exiting"
    exit 0
fi

if [ "${STABILIZE_STATUS1}" = true ]; then
    KINETICS_DESCRIPTION=status1_and_status2_load_stabilization
    STATUS1_DESCRIPTION="${BASE_LIFETIME}+lambda_force"
else
    KINETICS_DESCRIPTION=status2_only_load_stabilization
    STATUS1_DESCRIPTION="${BASE_LIFETIME}"
fi

cat > "${PARAM_FILE}" <<EOF
n_actins=${n_actins}
n_myosins=${n_myosins}
box=${Lx}x${Ly}x${Lz}
kinetics=${KINETICS_DESCRIPTION}
stabilize_status1=${STABILIZE_STATUS1}
status1_lifetime=${STATUS1_DESCRIPTION}
status2_lifetime=${BASE_LIFETIME}+${DIRECTIONAL_BASE_LIFETIME}+lambda_force
lifetime_coeff=${lifetime_coeff}
lifetime_label=${lifetime_label}
replica=${rep_idx}
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
    --stabilize_status1="${STABILIZE_STATUS1}"
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
current_step=-1
if [ -s "${OUTPUT_FILENAME}" ]; then
    # A killed checkpoint transaction can leave this temporary group behind.
    python - "${OUTPUT_FILENAME}" <<'PY' || true
import sys
import h5py

with h5py.File(sys.argv[1], "a") as handle:
    if "/resume_new" in handle:
        del handle["/resume_new"]
        print("deleted stale /resume_new")
PY

    current_step=$(python - "${OUTPUT_FILENAME}" "${PRE_PHASE_STEPS}" <<'PY' || echo -1
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
        "/state/current_step" in handle
        and "/actin/center" in handle
        and handle["/actin/center"].shape[0] > 0
        and checkpoint_step >= minimum_step
    )
print(checkpoint_step if valid else -1)
PY
    )
    current_step=$(echo "${current_step}" | tail -1 | tr -d '[:space:]')
    if [[ "${current_step}" =~ ^[0-9]+$ ]]; then
        resume_ready=1
    fi
fi

if [ "${resume_ready}" -eq 1 ] && [ "${current_step}" -ge "${nsteps}" ]; then
    echo "[${FILENAME}] already complete (${current_step}/${nsteps})"
    exit 0
elif [ "${resume_ready}" -eq 1 ]; then
    echo "[${FILENAME}] resuming from ${current_step} to ${nsteps}" | tee -a "${SIM_LOG}"
    cmd+=(--resume=true)
elif [ -s "${OUTPUT_FILENAME}" ]; then
    INCOMPLETE_FILE="${OUTPUT_FILENAME}.nonresumable_job${SLURM_JOB_ID}_task${task_id}"
    echo "[${FILENAME}] preserving non-resumable output as ${INCOMPLETE_FILE}" | tee -a "${SIM_LOG}"
    mv "${OUTPUT_FILENAME}" "${INCOMPLETE_FILE}"
    echo "[${FILENAME}] starting fresh ${nsteps}-step simulation" | tee -a "${SIM_LOG}"
else
    echo "[${FILENAME}] starting fresh ${nsteps}-step simulation" | tee -a "${SIM_LOG}"
fi

printf '[%s] Command:' "${FILENAME}" >> "${SIM_LOG}"
printf ' %q' "${cmd[@]}" >> "${SIM_LOG}"
printf '\n' >> "${SIM_LOG}"
"${cmd[@]}" >> "${SIM_LOG}" 2>&1
