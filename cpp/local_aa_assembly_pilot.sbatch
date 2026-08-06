#!/bin/bash
#SBATCH --partition=owners,hns,normal
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64GB
#SBATCH --time=2-00:00:00
#SBATCH --job-name=local_aa_pilot
#SBATCH --array=0-14%6
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/logs_local_aa_assembly_v1/%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs_local_aa_assembly_v1/%A_%a.err

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}"
export HDF5_USE_FILE_LOCKING=FALSE

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
EXE="${LOCAL_AA_EXE:-${CPP_DIR}/build_release/sarcomere}"
ROOT="${LOCAL_AA_ROOT:-/scratch/users/xsli/sarcomere/data_local_aa_assembly_v1}"
LOG="${LOCAL_AA_LOG:-/scratch/users/xsli/sarcomere/logs_local_aa_assembly_v1}"
PARAM="${LOCAL_AA_PARAM:-/scratch/users/xsli/sarcomere/params_local_aa_assembly_v1}"
mkdir -p "${ROOT}" "${LOG}" "${PARAM}"

task=${SLURM_ARRAY_TASK_ID:?}
condition_index=$((task / 3))
replicate=$((task % 3))
case "${condition_index}" in
    0) condition=aa_off;         aa_k_on=0;    links=2; lambda=0; stabilize1=false ;;
    1) condition=one_nocatch;    aa_k_on=1000; links=1; lambda=0; stabilize1=false ;;
    2) condition=two_nocatch;    aa_k_on=1000; links=2; lambda=0; stabilize1=false ;;
    3) condition=two_status2;    aa_k_on=1000; links=2; lambda=1; stabilize1=false ;;
    4) condition=two_status1;    aa_k_on=1000; links=2; lambda=1; stabilize1=true ;;
    *) echo "invalid task ${task}" >&2; exit 2 ;;
esac

N_ACTINS="${LOCAL_AA_N_ACTINS:-2000}"
N_MYOSINS="${LOCAL_AA_N_MYOSINS:-1000}"
LX="${LOCAL_AA_LX:-12}"
LY="${LOCAL_AA_LY:-2}"
LZ="${LOCAL_AA_LZ:-2}"
PRE_STEPS="${LOCAL_AA_PRE_STEPS:-50000}"
RAMP_STEPS="${LOCAL_AA_RAMP_STEPS:-25000}"
NSTEPS="${LOCAL_AA_NSTEPS:-200000}"
SAVE_EVERY="${LOCAL_AA_SAVE_EVERY:-1000}"
DT="${LOCAL_AA_DT:-1e-6}"
SEED_BASE="${LOCAL_AA_SEED_BASE:-204000}"
SEED=$((SEED_BASE + replicate))
POOL="${LOCAL_AA_LINK_POOL:-$((2 * N_ACTINS))}"
TAG="A${N_ACTINS}_M${N_MYOSINS}_${condition}_r${replicate}_seed${SEED}_L${LX}x${LY}x${LZ}"
OUTPUT="${ROOT}/${TAG}.h5"
PARAMS="${PARAM}/${TAG}.params"

test -x "${EXE}" || { echo "missing executable ${EXE}" >&2; exit 1; }
scontrol update jobid="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${task}" name="${TAG}" 2>/dev/null || true
exec 9>"${OUTPUT}.lock"
flock -n 9 || exit 0

resume_flag=()
if [ -s "${OUTPUT}" ] && python - "${OUTPUT}" <<'PY'
import h5py, sys
with h5py.File(sys.argv[1], "r", locking=False) as h:
    ok = "/resume/current_step" in h and "/actin/center" in h and h["/actin/center"].shape[0] > 0
raise SystemExit(0 if ok else 1)
PY
then
    resume_flag=(--resume=true)
elif [ -e "${OUTPUT}" ]; then
    mv "${OUTPUT}" "${OUTPUT}.incomplete_job${SLURM_JOB_ID}"
fi

exe_sha256=$(sha256sum "${EXE}" | awk '{print $1}')
cat >"${PARAMS}" <<EOF
model=sparse_local_actn2_explicit_crossbridge
condition=${condition}
replicate=${replicate}
n_actins=${N_ACTINS}
n_myosins=${N_MYOSINS}
box=${LX}x${LY}x${LZ}
seed=${SEED}
dt=${DT}
nsteps=${NSTEPS}
pre_phase_steps=${PRE_STEPS}
interaction_ramp_steps=${RAMP_STEPS}
aa_k_on=${aa_k_on}
aa_prebundle_links_per_pair=${links}
lifetime_coeff=${lambda}
stabilize_status1=${stabilize1}
aa_total_link_pool=${POOL}
executable=${EXE}
executable_sha256=${exe_sha256}
slurm_job_id=${SLURM_JOB_ID}
EOF

cmd=(srun --kill-on-bad-exit=1 "${EXE}"
    --n_actins="${N_ACTINS}" --n_myosins="${N_MYOSINS}"
    --Lx="${LX}" --Ly="${LY}" --Lz="${LZ}" --dimension=3 --pbc_mask=111
    --initial_structure=myosin_x_noise --pre_phase_steps="${PRE_STEPS}"
    --interaction_ramp_steps="${RAMP_STEPS}"
    --physical_mobility=true --temperature_kelvin=310
    --myosin_radius=0.025 --actin_radius=0.015 --am_surface_gap=0.015
    --actin_steric_k=100 --myosin_steric_k=1000 --am_steric_k=1000
    --aa_spacing_ratio=0.92 --aa_capture_margin=0.02 --aa_antiparallel_overlap=0.1
    --k_am=3000 --kappa_am=50 --k_aa=3000 --kappa_aa=50
    --dynamic_am=true --am_heads_per_half=35
    --am_k_on=1000 --am_k_weak_off=200 --am_k_powerstroke=200
    --am_k_strong_off=80 --am_force_scale=1 --am_motor_k=1000
    --am_power_stroke=0.005 --am_stall_force=4.5 --am_unloaded_velocity=5
    --am_cycle_period=0 --am_thick_feedback=true
    --aa_local_multilink=true --aa_prebundle_links_per_pair="${links}"
    --aa_site_exclusion=0.0192 --aa_overlap_capture=0.1
    --aa_perpendicular_k=50 --aa_total_link_pool="${POOL}"
    --aa_status_mode=productive_opposed
    --aa_status2_force_on=0.01 --aa_status2_force_off=0.005
    --k_on="${aa_k_on}" --base_lifetime=0.001 --directional_base_lifetime=0
    --lifetime_coeff="${lambda}" --stabilize_status1="${stabilize1}"
    --aa_energy_gated_association=true --aa_transition_state_theta=0.1
    --aa_encounter_k=10 --aa_encounter_kappa=0 --aa_passive_detailed_balance=false
    --tau_rec=0.05 --max_aa_bonds_per_actin=2
    --directional=true --sterics=true --k_mm=0 --k_bundle_max=0
    --bundle_ramp_steps=0 --max_myosin_bonds=6
    --dt="${DT}" --nsteps="${NSTEPS}" --save_every="${SAVE_EVERY}"
    --seed="${SEED}" --filename="${OUTPUT}")
if [ "${#resume_flag[@]}" -gt 0 ]; then cmd+=("${resume_flag[@]}"); fi
"${cmd[@]}"
