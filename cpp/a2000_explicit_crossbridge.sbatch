#!/bin/bash
#SBATCH --partition=owners,rotskoff
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=48GB
#SBATCH --time=2-00:00:00
#SBATCH --job-name=a2k_explicit_cb
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --array=0-2
#SBATCH --output=/scratch/users/xsli/sarcomere/logs_a2000_explicit_crossbridge/%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs_a2000_explicit_crossbridge/%A_%a.err

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}"
export HDF5_USE_FILE_LOCKING=FALSE

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
EXE="${CPP_DIR}/build_release/sarcomere"
SCRATCH=/scratch/users/xsli/sarcomere
OUT="${EXPLICIT_CB_OUTPUT_DIR:-${SCRATCH}/data_a2000_explicit_crossbridge_v3}"
LOG="${EXPLICIT_CB_LOG_DIR:-${SCRATCH}/logs_a2000_explicit_crossbridge}"
PARAM="${EXPLICIT_CB_PARAM_DIR:-${SCRATCH}/params_a2000_explicit_crossbridge_v3}"
mkdir -p "${OUT}" "${LOG}" "${PARAM}"

rep=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
printf -v rep_label 'r%02d' "${rep}"
seed=$((181000 + rep))
nsteps="${EXPLICIT_CB_NSTEPS:-300000}"
save_every="${EXPLICIT_CB_SAVE_EVERY:-200}"
lambda="${EXPLICIT_CB_LAMBDA:-0.01}"
tag="A2000_M1000_explicit_cb_lt1e-2_${rep_label}_L12x2x2"
output="${OUT}/${tag}.h5"
params="${PARAM}/${tag}.params"

if [ "${EXPLICIT_CB_SKIP_JOB_RENAME:-0}" != 1 ]; then
    job_element="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}"
    scontrol update jobid="${job_element}" name="${tag}" 2>/dev/null || true
fi
exec 9>"${output}.lock"
flock -n 9 || exit 0

resume_flag=
if [ -s "${output}" ]; then
    if python - "${output}" <<'PY'
import h5py, sys
with h5py.File(sys.argv[1], "r", locking=False) as h:
    ok = "/resume/current_step" in h or "/resume_prev/current_step" in h
raise SystemExit(0 if ok else 1)
PY
    then
        resume_flag=--resume
    else
        echo "existing output is not resumable: ${output}" >&2
        exit 1
    fi
fi

exe_sha256=$(sha256sum "${EXE}" | awk '{print $1}')
git_head=$(git -C "${CPP_DIR}" rev-parse HEAD 2>/dev/null || echo unknown)
cat >"${params}" <<EOF
model=explicit_crossbridge_v3
n_actins=2000
n_myosins=1000
box=12x2x2
nsteps=${nsteps}
save_every=${save_every}
dt=1e-5
seed=${seed}
aa_status_mode=opposed_and_transmitted
aa_load=min(0.5*norm(F_AM_i-F_AM_j),norm(F_AA))
aa_status2_force_on=0.01
aa_status2_force_off=0.005
base_lifetime=0.001
lifetime_coeff=${lambda}
stabilize_status1=false
am_k_on=1000
am_k_weak_off=200
am_k_powerstroke=200
am_k_strong_off=80
am_force_scale=1
am_motor_k=10
am_power_stroke=0.005
am_stall_force=1
executable_sha256=${exe_sha256}
git_head=${git_head}
slurm_job_id=${SLURM_JOB_ID}
EOF

srun --kill-on-bad-exit=1 "${EXE}" \
    --n_actins=2000 --n_myosins=1000 \
    --Lx=12 --Ly=2 --Lz=2 --dimension=3 --periodic_mask=111 \
    --myosin_radius=0.058 --am_cutoff=0.08 --am_optimal=0.06 \
    --initial_structure=myosin_x_noise --pre_phase_steps=500 \
    --dynamic_am=true \
    --am_k_on=1000 --am_k_weak_off=200 --am_k_powerstroke=200 \
    --am_k_strong_off=80 --am_force_scale=1 \
    --am_motor_k=10 --am_power_stroke=0.005 --am_stall_force=1 \
    --am_cycle_period=0 \
    --aa_status_mode=opposed_and_transmitted \
    --aa_status2_force_on=0.01 --aa_status2_force_off=0.005 \
    --base_lifetime=0.001 --directional_base_lifetime=0 \
    --lifetime_coeff="${lambda}" --stabilize_status1=false \
    --directional=true --k_mm=0 --k_bundle_max=0.1 \
    --bundle_ramp_steps=5000 --max_myosin_bonds=6 \
    --dt=1e-5 --nsteps="${nsteps}" --save_every="${save_every}" \
    --seed="${seed}" --filename="${output}" \
    ${resume_flag:+"${resume_flag}"}
