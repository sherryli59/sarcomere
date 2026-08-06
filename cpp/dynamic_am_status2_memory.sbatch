#!/bin/bash
#SBATCH --partition=owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --time=02:00:00
#SBATCH --job-name=am_status2_memory
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --array=0-159
#SBATCH --output=/scratch/users/xsli/sarcomere/logs_am_status2_memory/%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs_am_status2_memory/%A_%a.err

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export OMP_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
SARCOMERE_EXE="${CPP_DIR}/build_release/sarcomere"
SCRATCH_DIR=/scratch/users/xsli/sarcomere
OUTPUT_DIR="${AM_MEMORY_OUTPUT_DIR:-${SCRATCH_DIR}/data_am_status2_memory}"
LOG_DIR="${AM_MEMORY_LOG_DIR:-${SCRATCH_DIR}/logs_am_status2_memory}"
PARAM_DIR="${AM_MEMORY_PARAM_DIR:-${SCRATCH_DIR}/params_am_status2_memory}"
mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}" "${PARAM_DIR}"

if [ ! -x "${SARCOMERE_EXE}" ]; then
    echo "missing executable: ${SARCOMERE_EXE}" >&2
    exit 1
fi

GEOMETRIES=(free fixed)
MEMORIES=(0 0.001 0.005 0.02)
MEMORY_LABELS=(0 1ms 5ms 20ms)
N_REPLICAS=20
task_id=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
condition_idx=$((task_id / N_REPLICAS))
rep_idx=$((task_id % N_REPLICAS))
geometry_idx=$((condition_idx / 4))
memory_idx=$((condition_idx % 4))
geometry=${GEOMETRIES[$geometry_idx]}
memory=${MEMORIES[$memory_idx]}
memory_label=${MEMORY_LABELS[$memory_idx]}
printf -v rep_label 'r%02d' "${rep_idx}"

if [ "${geometry}" = free ]; then
    actin_diff=0.2
else
    actin_diff=0
fi

DT="${AM_MEMORY_DT:-1e-5}"
NSTEPS="${AM_MEMORY_NSTEPS:-200000}"
SAVE_EVERY="${AM_MEMORY_SAVE_EVERY:-20000}"
SEED=$((160000 + task_id))
condition="${geometry}_memory${memory_label}"
tag="am_xlink_${condition}_${rep_label}"
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
geometry=${geometry}
am_status2_memory=${memory}
base_lifetime=0.001
directional_base_lifetime=0.1
replica=${rep_label}
seed=${SEED}
dt=${DT}
nsteps=${NSTEPS}
EOF

srun --kill-on-bad-exit=1 "${SARCOMERE_EXE}" \
    --dynamic_am \
    --initial_structure=cb \
    --pre_phase_steps=0 \
    --n_actins=2 --n_myosins=4 \
    --n_fixed_myosins=4 \
    --actin_diff_coeff_trans="${actin_diff}" --actin_diff_coeff_rot=0 \
    --myosin_diff_coeff_trans=0.05 --myosin_diff_coeff_rot=0.05 \
    --Lx=4 --Ly=1 --Lz=1 --dimension=2 \
    --am_cutoff=0.08 --am_optimal=0.06 \
    --am_k_on=1000 --am_k_weak_off=200 --am_k_powerstroke=200 \
    --am_k_strong_off=80 --am_force_scale=1 \
    --am_motor_k=10 --am_power_stroke=0.005 --am_cycle_period=0 \
    --am_status2_memory="${memory}" \
    --k_on=1000 --base_lifetime=0.001 \
    --directional_base_lifetime=0.1 --lifetime_coeff=0 \
    --max_myosin_bonds=2 \
    --dt="${DT}" --nsteps="${NSTEPS}" --save_every="${SAVE_EVERY}" \
    --seed="${SEED}" --deterministic --filename="${output}" \
    ${resume_flag:+"${resume_flag}"}
