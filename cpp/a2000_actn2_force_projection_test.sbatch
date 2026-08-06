#!/bin/bash
#SBATCH --partition=owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=48GB
#SBATCH --time=12:00:00
#SBATCH --job-name=a2k_actn2_force_test
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/logs_a2000_actn2_force_test/%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs_a2000_actn2_force_test/%j.err

# Single integration test of the motor-limited, signed-AA-force ACTN2 status
# classifier in the established A=2000, M=1000, 12x2x2 geometry.

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}"
export HDF5_USE_FILE_LOCKING=FALSE

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
SARCOMERE_EXE="${CPP_DIR}/build_release/sarcomere"
SCRATCH_DIR=/scratch/users/xsli/sarcomere
OUTPUT_DIR="${ACTN2_FORCE_TEST_OUTPUT_DIR:-${SCRATCH_DIR}/data_a2000_actn2_force_test}"
LOG_DIR="${ACTN2_FORCE_TEST_LOG_DIR:-${SCRATCH_DIR}/logs_a2000_actn2_force_test}"
PARAM_DIR="${ACTN2_FORCE_TEST_PARAM_DIR:-${SCRATCH_DIR}/params_a2000_actn2_force_test}"
mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}" "${PARAM_DIR}"

if [ ! -x "${SARCOMERE_EXE}" ]; then
    echo "missing executable: ${SARCOMERE_EXE}" >&2
    exit 1
fi

N_ACTINS=2000
N_MYOSINS=1000
LX=12
LY=2
LZ=2
DT=1e-5
NSTEPS="${ACTN2_FORCE_TEST_NSTEPS:-20000}"
SAVE_EVERY="${ACTN2_FORCE_TEST_SAVE_EVERY:-200}"
SEED="${ACTN2_FORCE_TEST_SEED:-177001}"
PRE_PHASE_STEPS="${ACTN2_FORCE_TEST_PRE_PHASE_STEPS:-500}"

# Provisional model-unit operating point. The thresholds bracket the ~0.05
# initial strong-head spring load measured in fixed-geometry diagnostics.
FORCE_ON="${ACTN2_FORCE_TEST_FORCE_ON:-0.01}"
FORCE_OFF="${ACTN2_FORCE_TEST_FORCE_OFF:-0.005}"
BASE_LIFETIME="${ACTN2_FORCE_TEST_BASE_LIFETIME:-0.001}"
LIFETIME_COEFF="${ACTN2_FORCE_TEST_LAMBDA:-0.01}"

TAG="A${N_ACTINS}_M${N_MYOSINS}_actn2_forceproj_test_L${LX}x${LY}x${LZ}_seed${SEED}"
OUTPUT_FILE="${OUTPUT_DIR}/${TAG}.h5"
PARAM_FILE="${PARAM_DIR}/${TAG}.params"

exec 9>"${OUTPUT_FILE}.lock"
if ! flock -n 9; then
    echo "[${TAG}] output is already locked; exiting"
    exit 0
fi

resume_flag=""
if [ -s "${OUTPUT_FILE}" ]; then
    if python - "${OUTPUT_FILE}" <<'PY'
import h5py
import sys

with h5py.File(sys.argv[1], "r", locking=False) as handle:
    valid = "/resume/current_step" in handle or "/resume_prev/current_step" in handle
raise SystemExit(0 if valid else 1)
PY
    then
        resume_flag="--resume"
    else
        echo "existing output is not resumable: ${OUTPUT_FILE}" >&2
        exit 1
    fi
fi

exe_sha256=$(sha256sum "${SARCOMERE_EXE}" | awk '{print $1}')
git_head=$(git -C "${CPP_DIR}" rev-parse HEAD 2>/dev/null || echo unknown)
cat >"${PARAM_FILE}" <<EOF
n_actins=${N_ACTINS}
n_myosins=${N_MYOSINS}
box=${LX}x${LY}x${LZ}
dt=${DT}
nsteps=${NSTEPS}
save_every=${SAVE_EVERY}
seed=${SEED}
pre_phase_steps=${PRE_PHASE_STEPS}
dynamic_am=true
aa_status_mode=aa_force_projection
aa_status2_force_on=${FORCE_ON}
aa_status2_force_off=${FORCE_OFF}
aa_load_definition=min(favorable_raw_AA_arm_load,strong_motor_spring_load_per_actin)
base_lifetime=${BASE_LIFETIME}
lifetime_coeff=${LIFETIME_COEFF}
directional_base_lifetime=0
stabilize_status1=false
am_k_on=1000
am_k_weak_off=200
am_k_powerstroke=200
am_k_strong_off=80
am_force_scale=1
am_motor_k=10
am_power_stroke=0.005
am_cycle_period=0
max_myosin_bonds=6
executable=${SARCOMERE_EXE}
executable_sha256=${exe_sha256}
git_head=${git_head}
slurm_job_id=${SLURM_JOB_ID}
EOF

srun --kill-on-bad-exit=1 "${SARCOMERE_EXE}" \
    --n_actins="${N_ACTINS}" --n_myosins="${N_MYOSINS}" \
    --Lx="${LX}" --Ly="${LY}" --Lz="${LZ}" --dimension=3 \
    --periodic_mask=111 \
    --myosin_radius=0.058 --am_cutoff=0.08 --am_optimal=0.06 \
    --initial_structure=myosin_x_noise --pre_phase_steps="${PRE_PHASE_STEPS}" \
    --dynamic_am \
    --am_k_on=1000 --am_k_weak_off=200 --am_k_powerstroke=200 \
    --am_k_strong_off=80 --am_force_scale=1 \
    --am_motor_k=10 --am_power_stroke=0.005 --am_cycle_period=0 \
    --aa_status_mode=aa_force_projection \
    --aa_status2_force_on="${FORCE_ON}" \
    --aa_status2_force_off="${FORCE_OFF}" \
    --base_lifetime="${BASE_LIFETIME}" \
    --directional_base_lifetime=0 \
    --lifetime_coeff="${LIFETIME_COEFF}" \
    --stabilize_status1=false --directional=true \
    --k_mm=0 --k_bundle_max=0.1 --bundle_ramp_steps=5000 \
    --max_myosin_bonds=6 \
    --dt="${DT}" --nsteps="${NSTEPS}" --save_every="${SAVE_EVERY}" \
    --seed="${SEED}" --filename="${OUTPUT_FILE}" \
    ${resume_flag}
