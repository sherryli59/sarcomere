#!/bin/bash
#SBATCH --partition=owners,rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=48GB
#SBATCH --time=2-00:00:00
#SBATCH --job-name=a2k_tau1x
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/logs_a2000_tau1_sweep/%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs_a2000_tau1_sweep/%A_%a.err
#SBATCH --array=0-3
#
# EXTENSION of a2000_tau1_sweep.sbatch, lambda=1e-2 arm only.
#
# The base sweep topped out at tau1 = 3.16e-3 = 0.316*lambda at lt1e-2, so its
# recovered_frac only reached ~0.48 -- the curve was still climbing where the
# data stopped, and we could not tell whether a constant tau1 reaches the s1s2
# ceiling at lambda=1e-2 the way it does at lambda=1e-3 (tau1~lambda -> ~100%).
#
# This adds the two rungs that bracket the peak, mirroring lambda=1e-3 where
# tau1=lambda hit 103% and tau1=3.16*lambda OVERSHOT and fell below the floor:
#   tau1 = 1e-2     = 1.00 * lambda   (matched-to-old: expect ~ceiling)
#   tau1 = 3.16e-2  = 3.16 * lambda   (overshoot: expect degradation)
#
# Everything else -- executable, kinetics flags, seeds (97000+rep), naming, the
# flock/resume machinery -- is identical to the base sweep, so these files drop
# straight into data_a2000_tau1_sweep and phi3d_tau1_answer.py picks them up.

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}"

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
SARCOMERE_EXE="${CPP_DIR}/build_release/sarcomere"
SCRATCH_DIR=/scratch/users/xsli/sarcomere
OUTPUT_DIR="${TAU1_OUTPUT_DIR:-${SCRATCH_DIR}/data_a2000_tau1_sweep}"
LOG_DIR="${TAU1_LOG_DIR:-${SCRATCH_DIR}/logs_a2000_tau1_sweep}"
PARAM_DIR="${TAU1_PARAM_DIR:-${SCRATCH_DIR}/params_a2000_tau1_sweep}"
KINETICS_TAG="${TAU1_KINETICS_TAG:-s2catch}"
mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}" "${PARAM_DIR}"

if [ ! -x "${SARCOMERE_EXE}" ]; then
    echo "missing optimized executable: ${SARCOMERE_EXE}" >&2
    exit 1
fi

Lx=12
Ly=2
Lz=2
n_actins=2000
n_myosins=1000
nsteps="${TAU1_NSTEPS:-300000}"

# The two new rungs that bracket the peak at lambda=1e-2.
TAU1_VALUES=(1e-2 3.16228e-2)
TAU1_LABELS=(1e-2 1e-1.5)
LAMBDA_VALUES=(1e-2)
LAMBDA_LABELS=(1e-2)
N_TAUS=${#TAU1_VALUES[@]}
N_LAMBDAS=${#LAMBDA_VALUES[@]}
N_REPLICAS=2

task_id=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
PER_REP=$((N_TAUS * N_LAMBDAS))
if [ "${task_id}" -lt 0 ] || [ "${task_id}" -ge "$((PER_REP * N_REPLICAS))" ]; then
    echo "invalid array task ${task_id}" >&2
    exit 1
fi

rep_idx=$((task_id / PER_REP))
within=$((task_id % PER_REP))
lambda_idx=$((within / N_TAUS))
tau_idx=$((within % N_TAUS))

tau1=${TAU1_VALUES[$tau_idx]}
tau1_label=${TAU1_LABELS[$tau_idx]}
lifetime_coeff=${LAMBDA_VALUES[$lambda_idx]}
lambda_label=${LAMBDA_LABELS[$lambda_idx]}
printf -v rep_label 'r%02d' "${rep_idx}"

# Pair initial conditions with the s2only / s1s2 sweeps.
SEED=$((97000 + rep_idx))
MAX_MYOSIN_BONDS=6
MYOSIN_RADIUS=0.058
AM_CUTOFF=0.08
AM_OPTIMAL=0.06
SAVE_EVERY=200
DIRECTIONAL_BASE_LIFETIME=0.0
STABILIZE_STATUS1=false
K_MM=0
K_BUNDLE_MAX=0.1
BUNDLE_RAMP_STEPS=5000
PERIODIC_MASK=111
PRE_PHASE_STEPS="${TAU1_PRE_PHASE_STEPS:-500}"

FILENAME="A${n_actins}_M${n_myosins}_${KINETICS_TAG}_tau${tau1_label}_lt${lambda_label}_${rep_label}_L${Lx}x${Ly}x${Lz}"
OUTPUT_FILENAME="${OUTPUT_DIR}/${FILENAME}.h5"
SIM_LOG="${LOG_DIR}/${FILENAME}.sim.out"
PARAM_FILE="${PARAM_DIR}/${FILENAME}.params"

if [ "${TAU1_SKIP_JOB_RENAME:-0}" != 1 ]; then
    scontrol update jobid="${SLURM_JOB_ID}" name="${FILENAME}" 2>/dev/null || true
fi

# Prevent a watchdog retry and an original task from writing the same HDF5 file.
exec 9>"${OUTPUT_FILENAME}.lock"
if ! flock -n 9; then
    echo "[${FILENAME}] another process holds the output lock; exiting"
    exit 0
fi

cat > "${PARAM_FILE}" <<EOF
n_actins=${n_actins}
n_myosins=${n_myosins}
box=${Lx}x${Ly}x${Lz}
kinetics=constant_tau1_status2_only_load_stabilization
stabilize_status1=${STABILIZE_STATUS1}
status1_lifetime=${tau1}
status2_lifetime=${tau1}+${DIRECTIONAL_BASE_LIFETIME}+lambda_force
base_lifetime=${tau1}
tau1_label=${tau1_label}
lifetime_coeff=${lifetime_coeff}
lifetime_label=${lambda_label}
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
    --base_lifetime="${tau1}"
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
