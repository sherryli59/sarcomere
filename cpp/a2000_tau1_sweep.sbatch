#!/bin/bash
#SBATCH --partition=owners,rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=48GB
#SBATCH --time=2-00:00:00
#SBATCH --job-name=a2k_tau1
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/logs_a2000_tau1_sweep/%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs_a2000_tau1_sweep/%A_%a.err
#SBATCH --array=0-23
#
# CONSTANT-tau1 EXPERIMENT: is status 1 a physical short-lived precursor, or was
# the old assembly an artifact of load-stabilizing status 1?
#
# The old sweep (data0430) passed --base_lifetime=0 and never passed
# --directional_base_lifetime (default 0), so the OLD law collapsed to
#     lifetime(status 1) = lifetime(status 2) = lambda * f_load
# i.e. status 1 and status 2 were kinetically IDENTICAL, and since f_load is
# ~1.0 exactly whenever a pair is loaded (measured: median 1, max 1), the old
# status-1 lifetime was effectively the constant lambda.
#
# Here we instead give status 1 a constant, LOAD-INDEPENDENT lifetime tau1 and
# keep catch stabilization for status 2 only:
#     --base_lifetime=tau1  --stabilize_status1=false  --directional_base_lifetime=0
# which, per cpp/include/catch_bond_kinetics.h, yields
#     lifetime(status 1) = tau1
#     lifetime(status 2) = tau1 + lambda * f_load
# (base_lifetime necessarily floors both statuses -- directional_base_lifetime is
# clamped >= 0 so it cannot subtract. That floor is physical: a doubly engaged
# bond should never be less stable than a single one. Keep tau1 <~ lambda so the
# catch term still dominates status 2.)
#
# Reading:
#   modest tau1 restores assembly            -> status 1 is a real intermediate
#   only tau1 >~ lambda (i.e. ~ old law)     -> old assembly was a model artifact
#
# dt = 1e-5 and break_probability = dt/lifetime, so mean status-1 lifetime in
# steps is exactly tau1/1e-5:
#   tau1 = 0        ->    1 step  (identical to the existing s2only sweep)
#   tau1 = 3.16e-5  ->   ~3 steps
#   tau1 = 1e-4     ->   10 steps
#   tau1 = 3.16e-4  ->   32 steps
#   tau1 = 1e-3     ->  100 steps  (== lambda at lt1e-3: matched to the old law)
#   tau1 = 3.16e-3  ->  316 steps  (deliberately longer than lambda at lt1e-3)
# tau1 = 1e-5 is omitted on purpose: dt/tau1 = 1, identical to tau1 = 0.
#
# Array layout: 6 tau1 x 2 lambda x 2 replicas = 24 tasks. Seeds depend ONLY on
# the replica and match a2000_status2_lambda.sbatch / a2000_s1s2_lambda.sbatch
# (97000 + rep), so every trajectory here starts from the same stochastic initial
# condition as its s2only and s1s2 counterparts.
#
# Deliberately --no-requeue; the matching self-chaining watchdog
# (a2000_tau1_sweep_babycron.sbatch) resubmits only missing/stalled task IDs and
# the HDF5 checkpoint continues an interrupted trajectory.

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

# Constant, load-independent status-1 lifetimes (--base_lifetime).
TAU1_VALUES=(0 3.16228e-5 1e-4 3.16228e-4 1e-3 3.16228e-3)
TAU1_LABELS=(0 1e-4.5 1e-4 1e-3.5 1e-3 1e-2.5)
# Catch-bond strengths: the two regimes where the OLD runs carried the largest
# status-1 populations (f1 = 0.275 at 1e-3, 0.063 at 1e-2).
LAMBDA_VALUES=(1e-3 1e-2)
LAMBDA_LABELS=(1e-3 1e-2)
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
