#!/bin/bash
#SBATCH --partition=owners,normal,hns
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB
#SBATCH --time=2-00:00:00
#SBATCH --job-name=nalam_s2fstab
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/logs_nalam_s2fstab_base1e-3/%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs_nalam_s2fstab_base1e-3/%A_%a.err
#SBATCH --array=0-233
#
# CANONICAL (n_actins, lambda) sweep, new baseline kinetics.
#
# Model, fixed for this campaign:
#   --base_lifetime=1e-3         constant lifetime floor for BOTH status 1 and 2
#   --stabilize_status1=false    f_load (catch) stabilizes STATUS 2 ONLY
#   --directional_base_lifetime=0
#   --directional=true           antiparallel (cos<0) gate = the registration mechanism
# so, per cpp/include/catch_bond_kinetics.h:
#   lifetime(status 1) = 1e-3
#   lifetime(status 2) = 1e-3 + lambda * f_load
#
# This is the "constant tau1 = 1e-3, catch on status 2, directional" model that
# the tau1 sweep showed recovers full sarcomeric assembly without the unphysical
# load-dependent status-1 stabilization of the old (s1s2) runs.
#
# Grid (matches the old data0430 2D sweep): 6 n_actins x 13 lambda x 3 reps = 234.
#   n_actins in {600,900,1200,1500,1800,2100}, n_myosins = 1000 fixed
#   lambda   = 0 plus half-decades 10^-5.5 .. 10^0  (13 values)
# Seeds depend only on (n_actins, rep), so replicas are paired across lambda:
# every lambda at a fixed (n_actins, rep) starts from the same initial condition.
#
# Deliberately --no-requeue; the self-chaining watchdog
# (nalam_s2fstab_base1e-3_babycron.sbatch) resubmits only missing/stalled task
# IDs and the HDF5 checkpoint continues an interrupted trajectory.

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
SARCOMERE_EXE="${CPP_DIR}/build_release/sarcomere"
SCRATCH_DIR=/scratch/users/xsli/sarcomere
OUTPUT_DIR="${NAL_OUTPUT_DIR:-${SCRATCH_DIR}/data_nalam_s2fstab_base1e-3}"
LOG_DIR="${NAL_LOG_DIR:-${SCRATCH_DIR}/logs_nalam_s2fstab_base1e-3}"
PARAM_DIR="${NAL_PARAM_DIR:-${SCRATCH_DIR}/params_nalam_s2fstab_base1e-3}"
KINETICS_TAG="${NAL_KINETICS_TAG:-s2fstab}"
mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}" "${PARAM_DIR}"

if [ ! -x "${SARCOMERE_EXE}" ]; then
    echo "missing optimized executable: ${SARCOMERE_EXE}" >&2
    exit 1
fi

Lx=12
Ly=2
Lz=2
n_myosins=1000
nsteps="${NAL_NSTEPS:-300000}"

# ---- the sweep grid ----
N_ACTINS_VALUES=(600 900 1200 1500 1800 2100)
LIFETIME_VALUES=(0 3.16228e-6 1e-5 3.16228e-5 1e-4 3.16228e-4 1e-3 3.16228e-3 1e-2 3.16228e-2 1e-1 3.16228e-1 1)
LIFETIME_LABELS=(0 1e-5.5 1e-5 1e-4.5 1e-4 1e-3.5 1e-3 1e-2.5 1e-2 1e-1.5 1e-1 1e-0.5 1e0)
N_ACTINS=${#N_ACTINS_VALUES[@]}
N_LAMBDAS=${#LIFETIME_VALUES[@]}
N_REPLICAS=3

task_id=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
PER_ACTIN=$((N_LAMBDAS * N_REPLICAS))
if [ "${task_id}" -lt 0 ] || [ "${task_id}" -ge "$((N_ACTINS * PER_ACTIN))" ]; then
    echo "invalid array task ${task_id}" >&2
    exit 1
fi

actin_idx=$((task_id / PER_ACTIN))
within=$((task_id % PER_ACTIN))
rep_idx=$((within / N_LAMBDAS))
lambda_idx=$((within % N_LAMBDAS))

n_actins=${N_ACTINS_VALUES[$actin_idx]}
lifetime_coeff=${LIFETIME_VALUES[$lambda_idx]}
lifetime_label=${LIFETIME_LABELS[$lambda_idx]}
printf -v rep_label 'r%02d' "${rep_idx}"

# Distinct IC per (n_actins, rep); paired across lambda.
SEED=$((97000 + actin_idx * 100 + rep_idx))

# ---- fixed baseline kinetics ----
BASE_LIFETIME=1e-3
DIRECTIONAL_BASE_LIFETIME=0.0
STABILIZE_STATUS1=false

# ---- other physics params (match the a2000 catch-bond campaigns) ----
MAX_MYOSIN_BONDS=6
MYOSIN_RADIUS=0.058
AM_CUTOFF=0.08
AM_OPTIMAL=0.06
SAVE_EVERY=200
K_MM=0
K_BUNDLE_MAX=0.1
BUNDLE_RAMP_STEPS=5000
PERIODIC_MASK=111
PRE_PHASE_STEPS="${NAL_PRE_PHASE_STEPS:-500}"

FILENAME="A${n_actins}_M${n_myosins}_${KINETICS_TAG}_lt${lifetime_label}_${rep_label}_L${Lx}x${Ly}x${Lz}"
OUTPUT_FILENAME="${OUTPUT_DIR}/${FILENAME}.h5"
SIM_LOG="${LOG_DIR}/${FILENAME}.sim.out"
PARAM_FILE="${PARAM_DIR}/${FILENAME}.params"

if [ "${NAL_SKIP_JOB_RENAME:-0}" != 1 ]; then
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
kinetics=status2_only_load_stabilization_const_base_lifetime
stabilize_status1=${STABILIZE_STATUS1}
status1_lifetime=${BASE_LIFETIME}
status2_lifetime=${BASE_LIFETIME}+${DIRECTIONAL_BASE_LIFETIME}+lambda_force
base_lifetime=${BASE_LIFETIME}
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
