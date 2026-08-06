#!/bin/bash
#SBATCH --partition=owners,hns,normal
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=48GB
#SBATCH --time=12:00:00
#SBATCH --job-name=a2k_theta01_phys
#SBATCH --array=0-1
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/logs_a2000_theta01_physical_pilot/%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs_a2000_theta01_physical_pilot/%A_%a.err

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}"
export HDF5_USE_FILE_LOCKING=FALSE

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
EXE="${CPP_DIR}/build_release/sarcomere"
SCRATCH=/scratch/users/xsli/sarcomere
OUT="${THETA_PILOT_OUTPUT_DIR:-${SCRATCH}/data_a2000_theta01_physical_pilot_overlap_v3}"
LOG="${THETA_PILOT_LOG_DIR:-${SCRATCH}/logs_a2000_theta01_physical_pilot_v2}"
PARAM="${THETA_PILOT_PARAM_DIR:-${SCRATCH}/params_a2000_theta01_physical_pilot_overlap_v3}"
mkdir -p "${OUT}" "${LOG}" "${PARAM}"

case "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}" in
    0) condition=ungated; energy_gate=false ;;
    1) condition=theta0p1; energy_gate=true ;;
    *) echo "invalid array index ${SLURM_ARRAY_TASK_ID}" >&2; exit 2 ;;
esac

N_ACTINS=2000
N_MYOSINS=1000
SEED="${THETA_PILOT_SEED:-203001}"
NSTEPS="${THETA_PILOT_NSTEPS:-10000}"
SAVE_EVERY="${THETA_PILOT_SAVE_EVERY:-200}"
PRE_PHASE_STEPS="${THETA_PILOT_PRE_PHASE_STEPS:-20000}"
INTERACTION_RAMP_STEPS="${THETA_PILOT_INTERACTION_RAMP_STEPS:-10000}"
DT=1e-6
THETA=0.1
TAG="A${N_ACTINS}_M${N_MYOSINS}_physical_v2_${condition}_seed${SEED}_L12x2x2"
OUTPUT="${OUT}/${TAG}.h5"
PARAMS="${PARAM}/${TAG}.params"

if [ ! -x "${EXE}" ]; then
    echo "missing executable: ${EXE}" >&2
    exit 1
fi

job_element="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID}"
scontrol update jobid="${job_element}" name="${TAG}" 2>/dev/null || true

exec 9>"${OUTPUT}.lock"
flock -n 9 || exit 0

# The checkpoint-aware babysitter resubmits a task after preemption. On that
# fresh batch invocation, continue from the last committed HDF5 checkpoint
# instead of reopening the trajectory as a new file. A killed checkpoint
# transaction is repaired by the resume loader; a non-resumable partial file is
# preserved for diagnosis.
resume_flag=()
if [ -s "${OUTPUT}" ] && python - "${OUTPUT}" <<'PY'
import h5py
import sys

with h5py.File(sys.argv[1], "r", locking=False) as handle:
    valid = (
        ("/resume/current_step" in handle or "/resume_prev/current_step" in handle)
        and "/actin/center" in handle
        and handle["/actin/center"].shape[0] > 0
    )
raise SystemExit(0 if valid else 1)
PY
then
    echo "[${TAG}] requeued/retried task found a checkpoint; resuming"
    resume_flag=(--resume=true)
elif [ -e "${OUTPUT}" ]; then
    preserved="${OUTPUT}.incomplete_job${SLURM_JOB_ID}_task${SLURM_ARRAY_TASK_ID}"
    echo "[${TAG}] preserving non-resumable partial output as ${preserved}" >&2
    mv "${OUTPUT}" "${preserved}"
fi

exe_sha256=$(sha256sum "${EXE}" | awk '{print $1}')
git_head=$(git -C "${CPP_DIR}" rev-parse HEAD 2>/dev/null || echo unknown)
cat >"${PARAMS}" <<EOF
model=explicit_crossbridge_physical_theta_pilot
condition=${condition}
aa_energy_gated_association=${energy_gate}
aa_transition_state_theta=${THETA}
aa_passive_detailed_balance=false
aa_encounter_k=10
aa_encounter_kappa=0
n_actins=${N_ACTINS}
n_myosins=${N_MYOSINS}
box=12x2x2
seed=${SEED}
dt=${DT}
nsteps=${NSTEPS}
save_every=${SAVE_EVERY}
physical_mobility=true
temperature_kelvin=310
pre_phase_steps=${PRE_PHASE_STEPS}
interaction_ramp_steps=${INTERACTION_RAMP_STEPS}
actin_steric_k=100
myosin_steric_k=1000
k_am=3000
k_aa=3000
kappa_am=50
kappa_aa=50
aa_spacing_ratio=0.92
aa_antiparallel_overlap=0.1
max_aa_bonds_per_actin=2
base_lifetime=0.001
lifetime_coeff=0.01
stabilize_status1=false
dynamic_am=true
bundle_scaffold=false
executable=${EXE}
executable_sha256=${exe_sha256}
git_head=${git_head}
slurm_job_id=${SLURM_JOB_ID}
EOF

cmd=(srun --kill-on-bad-exit=1 "${EXE}" \
    --n_actins="${N_ACTINS}" --n_myosins="${N_MYOSINS}" \
    --Lx=12 --Ly=2 --Lz=2 --dimension=3 --pbc_mask=111 \
    --initial_structure=myosin_x_noise --pre_phase_steps="${PRE_PHASE_STEPS}" \
    --interaction_ramp_steps="${INTERACTION_RAMP_STEPS}" \
    --physical_mobility=true --temperature_kelvin=310 \
    --myosin_radius=0.025 --actin_radius=0.015 --am_surface_gap=0.015 \
    --actin_steric_k=100 --myosin_steric_k=1000 --am_steric_k=1000 \
    --aa_spacing_ratio=0.92 --aa_capture_margin=0.02 \
    --aa_antiparallel_overlap=0.1 \
    --k_am=3000 --kappa_am=50 --k_aa=3000 --kappa_aa=50 \
    --dynamic_am=true \
    --am_k_on=1000 --am_k_weak_off=200 --am_k_powerstroke=200 \
    --am_k_strong_off=80 --am_force_scale=1 \
    --am_motor_k=1000 --am_power_stroke=0.005 --am_stall_force=4.5 \
    --am_unloaded_velocity=5 \
    --am_cycle_period=0 \
    --aa_status_mode=opposed_and_transmitted \
    --aa_status2_force_on=0.01 --aa_status2_force_off=0.005 \
    --base_lifetime=0.001 --directional_base_lifetime=0 \
    --lifetime_coeff=0.01 --stabilize_status1=false \
    --aa_energy_gated_association="${energy_gate}" \
    --aa_transition_state_theta="${THETA}" \
    --aa_encounter_k=10 --aa_encounter_kappa=0 \
    --aa_passive_detailed_balance=false \
    --tau_rec=0.05 --max_aa_bonds_per_actin=2 \
    --directional=true --sterics=true \
    --k_mm=0 --k_bundle_max=0 --bundle_ramp_steps=0 \
    --max_myosin_bonds=6 \
    --dt="${DT}" --nsteps="${NSTEPS}" --save_every="${SAVE_EVERY}" \
    --seed="${SEED}" --filename="${OUTPUT}")
if [ "${#resume_flag[@]}" -gt 0 ]; then
    cmd+=("${resume_flag[@]}")
fi
"${cmd[@]}"
