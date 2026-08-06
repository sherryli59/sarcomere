#!/bin/bash
#SBATCH -p rotskoff,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32GB
#SBATCH --time=1-00:00:00
#SBATCH --job-name=s2_lifetime_match
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --array=0-31%8
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
OUTPUT_DIR="${CONTROL_OUTPUT_DIR:-${SCRATCH_DIR}/data_status2_lifetime_match}"
PARAM_DIR="${CONTROL_PARAM_DIR:-${SCRATCH_DIR}/params_status2_lifetime_match}"
LOG_DIR="${CONTROL_LOG_DIR:-${SCRATCH_DIR}/logs}"
mkdir -p "${OUTPUT_DIR}" "${PARAM_DIR}" "${LOG_DIR}"

if [ ! -x "${SARCOMERE_EXE}" ]; then
    echo "Missing executable: ${SARCOMERE_EXE}" >&2
    exit 1
fi

# Arm 0 is the status-2-only catch reference.  The other 15 arms span a 5x3
# grid of constant status-2 lifetimes and k_on values.  Two kinetic knobs are
# required because a single constant lifetime generally cannot match both
# occupancy and turnover independently.
TAU_VALUES=(0.001 0.002 0.003 0.004 0.006)
KON_VALUES=(500 1000 2000)
N_CONST_ARMS=$((${#TAU_VALUES[@]} * ${#KON_VALUES[@]}))
NREPS=2

task_id=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
arm_idx=$((task_id / NREPS))
rep_idx=$((task_id % NREPS))
rep_label=$(printf 'r%02d' "${rep_idx}")
SEED=$((93000 + rep_idx))

CATCH_LAMBDA="${CONTROL_CATCH_LAMBDA:-0.0031622776601683794}"
REFERENCE_KON="${CONTROL_REFERENCE_KON:-1000}"
if (( arm_idx == 0 )); then
    tag="catch_lt0p003162_kon${REFERENCE_KON}"
    lifetime_coeff="${CATCH_LAMBDA}"
    directional_base_lifetime=0.0
    k_on="${REFERENCE_KON}"
else
    const_idx=$((arm_idx - 1))
    if (( const_idx >= N_CONST_ARMS )); then
        echo "Invalid arm index ${arm_idx}" >&2
        exit 1
    fi
    tau_idx=$((const_idx / ${#KON_VALUES[@]}))
    kon_idx=$((const_idx % ${#KON_VALUES[@]}))
    directional_base_lifetime=${TAU_VALUES[$tau_idx]}
    k_on=${KON_VALUES[$kon_idx]}
    lifetime_coeff=0.0
    tau_label=${directional_base_lifetime/./p}
    tag="consts2_tau${tau_label}_kon${k_on}"
fi

Lx=12
Ly=2
Lz=2
n_actins="${CONTROL_N_ACTINS:-2000}"
n_myosins="${CONTROL_N_MYOSINS:-1000}"
nsteps="${CONTROL_NSTEPS:-150000}"
SAVE_EVERY="${CONTROL_SAVE_EVERY:-100}"
BASE_LIFETIME=0.0
MYOSIN_RADIUS=0.058
AM_CUTOFF=0.08
AM_OPTIMAL=0.06
MAX_MYOSIN_BONDS=6
TAU_REC=0.05
K_MM=0
K_BUNDLE_MAX=0.1
BUNDLE_RAMP_STEPS=5000
PBC_MASK=111

FILENAME="A${n_actins}_M${n_myosins}_${tag}_${rep_label}_L${Lx}x${Ly}x${Lz}"
OUTPUT_FILENAME="${OUTPUT_DIR}/${FILENAME}.h5"
PARAM_FILE="${PARAM_DIR}/${FILENAME}.params"
SIM_LOG="${LOG_DIR}/${FILENAME}.sim.out"

exec 9>"${OUTPUT_FILENAME}.lock"
if ! flock -n 9; then
    echo "[${FILENAME}] another process holds the output lock; exiting"
    exit 0
fi

resume_ready=0
if [ -s "${OUTPUT_FILENAME}" ]; then
    if python - "${OUTPUT_FILENAME}" <<'PY'
import sys
import h5py

with h5py.File(sys.argv[1], "r") as handle:
    valid = any(path in handle for path in ("/resume/current_step", "/resume_prev/current_step"))
raise SystemExit(0 if valid else 1)
PY
    then
        resume_ready=1
    else
        echo "Existing output is not resumable: ${OUTPUT_FILENAME}" >&2
        exit 1
    fi
fi

cat >"${PARAM_FILE}" <<EOF
tag=${tag}
replica=${rep_label}
seed=${SEED}
lifetime_coeff=${lifetime_coeff}
status1_base_lifetime=${BASE_LIFETIME}
status2_directional_base_lifetime=${directional_base_lifetime}
k_on=${k_on}
nsteps=${nsteps}
n_actins=${n_actins}
n_myosins=${n_myosins}
box=${Lx}x${Ly}x${Lz}
save_every=${SAVE_EVERY}
EOF

cmd=(
    srun --kill-on-bad-exit=1 "${SARCOMERE_EXE}"
    --myosin_radius="${MYOSIN_RADIUS}"
    --am_cutoff="${AM_CUTOFF}"
    --am_optimal="${AM_OPTIMAL}"
    --n_actins="${n_actins}"
    --n_myosins="${n_myosins}"
    --filename="${OUTPUT_FILENAME}"
    --pbc_mask="${PBC_MASK}"
    --Lx="${Lx}" --Ly="${Ly}" --Lz="${Lz}"
    --nsteps="${nsteps}"
    --k_on="${k_on}"
    --lifetime_coeff="${lifetime_coeff}"
    --base_lifetime="${BASE_LIFETIME}"
    --directional_base_lifetime="${directional_base_lifetime}"
    --tau_rec="${TAU_REC}"
    --directional=true
    --initial_structure=myosin_x_noise
    --pre_phase_steps=500
    --save_every="${SAVE_EVERY}"
    --k_mm="${K_MM}"
    --k_bundle_max="${K_BUNDLE_MAX}"
    --max_myosin_bonds="${MAX_MYOSIN_BONDS}"
    --bundle_ramp_steps="${BUNDLE_RAMP_STEPS}"
    --seed="${SEED}"
)

if [ "${resume_ready}" -eq 1 ]; then
    cmd+=(--resume=true)
fi

printf '[%s] Command:' "${FILENAME}" | tee -a "${SIM_LOG}"
printf ' %q' "${cmd[@]}" | tee -a "${SIM_LOG}"
printf '\n' | tee -a "${SIM_LOG}"
"${cmd[@]}" >>"${SIM_LOG}" 2>&1
