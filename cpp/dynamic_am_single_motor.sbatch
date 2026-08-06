#!/bin/bash
#SBATCH --partition=owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --time=02:00:00
#SBATCH --job-name=am_single_motor
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --array=0-79
#SBATCH --output=/scratch/users/xsli/sarcomere/logs_am_single_motor/%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs_am_single_motor/%A_%a.err

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export OMP_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
SARCOMERE_EXE="${CPP_DIR}/build_release/sarcomere"
SCRATCH_DIR=/scratch/users/xsli/sarcomere
OUTPUT_DIR="${AM_SINGLE_OUTPUT_DIR:-${SCRATCH_DIR}/data_am_single_motor}"
LOG_DIR="${AM_SINGLE_LOG_DIR:-${SCRATCH_DIR}/logs_am_single_motor}"
PARAM_DIR="${AM_SINGLE_PARAM_DIR:-${SCRATCH_DIR}/params_am_single_motor}"
mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}" "${PARAM_DIR}"

if [ ! -x "${SARCOMERE_EXE}" ]; then
    echo "missing executable: ${SARCOMERE_EXE}" >&2
    exit 1
fi

FORCE_RATIOS=(0 0.5 1 2)
FORCE_LABELS=(0 0p5 1 2)
N_REPLICAS=20
task_id=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
ratio_idx=$((task_id / N_REPLICAS))
rep_idx=$((task_id % N_REPLICAS))
force_ratio=${FORCE_RATIOS[$ratio_idx]}
force_label=${FORCE_LABELS[$ratio_idx]}
printf -v rep_label 'r%02d' "${rep_idx}"

DT="${AM_SINGLE_DT:-1e-5}"
NSTEPS="${AM_SINGLE_NSTEPS:-200000}"
SAVE_EVERY="${AM_SINGLE_SAVE_EVERY:-20000}"
MOTOR_K="${AM_SINGLE_MOTOR_K:-10}"
FORCE_SCALE="${AM_SINGLE_FORCE_SCALE:-1}"
STALL_FORCE="${AM_SINGLE_STALL_FORCE:-10}"
POWER_STROKE=$(awk -v ratio="${force_ratio}" -v scale="${FORCE_SCALE}" -v k="${MOTOR_K}" 'BEGIN {printf "%.17g", ratio*scale/k}')
K_ON="${AM_SINGLE_K_ON:-1000}"
K_WEAK_OFF="${AM_SINGLE_K_WEAK_OFF:-200}"
K_POWERSTROKE="${AM_SINGLE_K_POWERSTROKE:-200}"
K_STRONG_OFF="${AM_SINGLE_K_STRONG_OFF:-80}"
SEED=$((120000 + task_id))

tag="am_single_F${force_label}_${rep_label}"
output="${OUTPUT_DIR}/${tag}.h5"
params="${PARAM_DIR}/${tag}.params"
exec 9>"${output}.lock"
flock -n 9 || exit 0
resume_flag=
if [ -s "${output}" ]; then
    if HDF5_USE_FILE_LOCKING=FALSE python - "${output}" <<'PY'
import sys
import h5py
with h5py.File(sys.argv[1], "r", locking=False) as handle:
    raise SystemExit(0 if "/resume/current_step" in handle or "/resume_prev/current_step" in handle else 1)
PY
    then
        resume_flag=--resume
    else
        echo "existing output is not resumable: ${output}" >&2
        exit 1
    fi
fi

cat >"${params}" <<EOF
force_ratio=${force_ratio}
force_scale=${FORCE_SCALE}
stall_force=${STALL_FORCE}
motor_k=${MOTOR_K}
power_stroke=${POWER_STROKE}
k_on=${K_ON}
k_weak_off=${K_WEAK_OFF}
k_powerstroke=${K_POWERSTROKE}
k_strong_off=${K_STRONG_OFF}
dt=${DT}
nsteps=${NSTEPS}
seed=${SEED}
EOF

srun --kill-on-bad-exit=1 "${SARCOMERE_EXE}" \
    --dynamic_am \
    --initial_structure=single_am \
    --pre_phase_steps=0 \
    --n_actins=1 --n_myosins=1 \
    --n_fixed_myosins=1 \
    --actin_diff_coeff_trans=0 --actin_diff_coeff_rot=0 \
    --myosin_diff_coeff_trans=0.05 --myosin_diff_coeff_rot=0.05 \
    --Lx=4 --Ly=1 --Lz=1 --dimension=2 \
    --am_cutoff=0.08 --am_optimal=0.06 \
    --am_k_on="${K_ON}" \
    --am_k_weak_off="${K_WEAK_OFF}" \
    --am_k_powerstroke="${K_POWERSTROKE}" \
    --am_k_strong_off="${K_STRONG_OFF}" \
    --am_force_scale="${FORCE_SCALE}" \
    --am_motor_k="${MOTOR_K}" \
    --am_power_stroke="${POWER_STROKE}" \
    --am_stall_force="${STALL_FORCE}" \
    --am_cycle_period=0 \
    --base_lifetime=0 --directional_base_lifetime=0 --lifetime_coeff=0 \
    --max_myosin_bonds=1 \
    --dt="${DT}" --nsteps="${NSTEPS}" --save_every="${SAVE_EVERY}" \
    --seed="${SEED}" --deterministic --filename="${output}" \
    ${resume_flag:+"${resume_flag}"}
