#!/bin/bash
#SBATCH -p rotskoff,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32GB
#SBATCH --time=2-00:00:00
#SBATCH --job-name=cb_history_s2
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --array=0-11
#SBATCH --error=/scratch/users/xsli/sarcomere/logs/%A_%a.err
#SBATCH --output=/scratch/users/xsli/sarcomere/logs/%A_%a.out

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}"
export HDF5_USE_FILE_LOCKING=FALSE

CPP_DIR="/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp"
SARCOMERE_EXE="${CPP_DIR}/build_release/sarcomere"
SCRATCH_DIR="/scratch/users/xsli/sarcomere"
OUTPUT_DIR="${HISTORY_OUTPUT_DIR:-${SCRATCH_DIR}/data_history_status2}"
PARAM_DIR="${HISTORY_PARAM_DIR:-${SCRATCH_DIR}/params_history_status2}"
LOG_DIR="${HISTORY_LOG_DIR:-${SCRATCH_DIR}/logs}"
mkdir -p "${OUTPUT_DIR}" "${PARAM_DIR}" "${LOG_DIR}"

if [ ! -x "${SARCOMERE_EXE}" ]; then
    echo "Missing executable: ${SARCOMERE_EXE}" >&2
    exit 1
fi

# Four paired histories x three replicas.  Protocols share the same seed for a
# given replica, so their initial condition and stochastic stream are paired.
PROTOCOLS=(medium_control high_control medium_high high_medium_high)
NREPS=3
task_id=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
protocol=${PROTOCOLS[$((task_id / NREPS))]}
rep_idx=$((task_id % NREPS))
rep_label=$(printf 'r%02d' "${rep_idx}")
SEED=$((92000 + rep_idx))

MEDIUM_LAMBDA="${HISTORY_MEDIUM_LAMBDA:-0.0031622776601683794}"
HIGH_LAMBDA="${HISTORY_HIGH_LAMBDA:-0.1}"
TOTAL_STEPS="${HISTORY_TOTAL_STEPS:-400000}"
MEDIUM_HIGH_SWITCH="${HISTORY_MEDIUM_HIGH_SWITCH:-200000}"
PULSE_START="${HISTORY_PULSE_START:-100000}"
PULSE_END="${HISTORY_PULSE_END:-300000}"

if ! (( PULSE_START < PULSE_END && PULSE_END < TOTAL_STEPS )); then
    echo "Require PULSE_START < PULSE_END < TOTAL_STEPS" >&2
    exit 1
fi
if ! (( MEDIUM_HIGH_SWITCH > 0 && MEDIUM_HIGH_SWITCH < TOTAL_STEPS )); then
    echo "Require 0 < MEDIUM_HIGH_SWITCH < TOTAL_STEPS" >&2
    exit 1
fi

Lx=12
Ly=2
Lz=2
n_actins="${HISTORY_N_ACTINS:-2000}"
n_myosins="${HISTORY_N_MYOSINS:-1000}"
SAVE_EVERY="${HISTORY_SAVE_EVERY:-100}"
MYOSIN_RADIUS=0.058
AM_CUTOFF=0.08
AM_OPTIMAL=0.06
MAX_MYOSIN_BONDS=6
BASE_LIFETIME=0.0
DIRECTIONAL_BASE_LIFETIME=0.0
TAU_REC=0.05
K_MM=0
K_BUNDLE_MAX=0.1
BUNDLE_RAMP_STEPS=5000
PBC_MASK=111

FILENAME="A${n_actins}_M${n_myosins}_history_${protocol}_${rep_label}_L${Lx}x${Ly}x${Lz}"
OUTPUT_FILENAME="${OUTPUT_DIR}/${FILENAME}.h5"
PARAM_FILE="${PARAM_DIR}/${FILENAME}.params"
SIM_LOG="${LOG_DIR}/${FILENAME}.sim.out"

exec 9>"${OUTPUT_FILENAME}.lock"
if ! flock -n 9; then
    echo "[${FILENAME}] another process holds the output lock; exiting"
    exit 0
fi

current_step() {
    if [ ! -s "${OUTPUT_FILENAME}" ]; then
        echo 0
        return
    fi
    python - "${OUTPUT_FILENAME}" <<'PY'
import sys
import h5py

with h5py.File(sys.argv[1], "r") as handle:
    for path in ("/resume/current_step", "/resume_prev/current_step"):
        if path in handle:
            print(int(handle[path][...].reshape(-1)[-1]))
            break
    else:
        raise SystemExit("existing output has no resumable checkpoint")
PY
}

common_args=(
    --myosin_radius="${MYOSIN_RADIUS}"
    --am_cutoff="${AM_CUTOFF}"
    --am_optimal="${AM_OPTIMAL}"
    --n_actins="${n_actins}"
    --n_myosins="${n_myosins}"
    --filename="${OUTPUT_FILENAME}"
    --pbc_mask="${PBC_MASK}"
    --Lx="${Lx}" --Ly="${Ly}" --Lz="${Lz}"
    --directional=true
    --initial_structure=myosin_x_noise
    --pre_phase_steps=500
    --save_every="${SAVE_EVERY}"
    --base_lifetime="${BASE_LIFETIME}"
    --directional_base_lifetime="${DIRECTIONAL_BASE_LIFETIME}"
    --tau_rec="${TAU_REC}"
    --k_mm="${K_MM}"
    --k_bundle_max="${K_BUNDLE_MAX}"
    --max_myosin_bonds="${MAX_MYOSIN_BONDS}"
    --bundle_ramp_steps="${BUNDLE_RAMP_STEPS}"
    --seed="${SEED}"
)

advance_to() {
    local target=$1
    local lambda=$2
    local step
    step=$(current_step)
    if (( step >= target )); then
        echo "[${FILENAME}] checkpoint ${step} already reached target ${target}"
        return
    fi

    local cmd=(
        srun --kill-on-bad-exit=1 "${SARCOMERE_EXE}"
        "${common_args[@]}"
        --nsteps="${target}"
        --lifetime_coeff="${lambda}"
    )
    if [ -s "${OUTPUT_FILENAME}" ]; then
        cmd+=(--resume=true)
    fi

    printf '[%s] segment %d -> %d, lambda=%s:' "${FILENAME}" "${step}" "${target}" "${lambda}" | tee -a "${SIM_LOG}"
    printf ' %q' "${cmd[@]}" | tee -a "${SIM_LOG}"
    printf '\n' | tee -a "${SIM_LOG}"
    "${cmd[@]}" >>"${SIM_LOG}" 2>&1
}

cat >"${PARAM_FILE}" <<EOF
protocol=${protocol}
replica=${rep_label}
seed=${SEED}
medium_lambda=${MEDIUM_LAMBDA}
high_lambda=${HIGH_LAMBDA}
total_steps=${TOTAL_STEPS}
medium_high_switch=${MEDIUM_HIGH_SWITCH}
pulse_start=${PULSE_START}
pulse_end=${PULSE_END}
status1_lifetime=${BASE_LIFETIME}
status2_directional_base_lifetime=${DIRECTIONAL_BASE_LIFETIME}
n_actins=${n_actins}
n_myosins=${n_myosins}
box=${Lx}x${Ly}x${Lz}
save_every=${SAVE_EVERY}
EOF

case "${protocol}" in
    medium_control)
        advance_to "${TOTAL_STEPS}" "${MEDIUM_LAMBDA}"
        ;;
    high_control)
        advance_to "${TOTAL_STEPS}" "${HIGH_LAMBDA}"
        ;;
    medium_high)
        advance_to "${MEDIUM_HIGH_SWITCH}" "${MEDIUM_LAMBDA}"
        advance_to "${TOTAL_STEPS}" "${HIGH_LAMBDA}"
        ;;
    high_medium_high)
        advance_to "${PULSE_START}" "${HIGH_LAMBDA}"
        advance_to "${PULSE_END}" "${MEDIUM_LAMBDA}"
        advance_to "${TOTAL_STEPS}" "${HIGH_LAMBDA}"
        ;;
    *)
        echo "Unknown protocol: ${protocol}" >&2
        exit 1
        ;;
esac

echo "[${FILENAME}] protocol complete at step $(current_step)"
