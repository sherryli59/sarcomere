#!/bin/bash
#SBATCH --partition=hns
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=2-00:00:00
#SBATCH --job-name=prod_status2_reg
#SBATCH --open-mode=append
#SBATCH --requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/prod_status2_register_20260726/logs/worker_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/prod_status2_register_20260726/logs/worker_%j.err

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OMP_DYNAMIC=FALSE
export OMP_SCHEDULE=static
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export MKL_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
JEMALLOC=/home/groups/rotskoff/sherry/mamba_root/envs/sarcomere/lib/libjemalloc.so
test -r "${JEMALLOC}" || { echo "missing jemalloc: ${JEMALLOC}" >&2; exit 1; }
export LD_PRELOAD="${JEMALLOC}"

ROOT="${PROD_ROOT:-/scratch/users/xsli/sarcomere/prod_status2_register_20260726}"
EXE="${PROD_EXE:-${ROOT}/bin/sarcomere_status2_register}"
OUTPUT="${PROD_OUTPUT:-${ROOT}/data/A1280_M640_L7p72_status2_register_seed2072601.h5}"
PARAMS="${PROD_PARAMS:-${ROOT}/params/A1280_M640_L7p72_status2_register_seed2072601.params}"
NSTEPS="${PROD_NSTEPS:-1000000}"
SAVE_EVERY="${PROD_SAVE_EVERY:-1000}"
SEED="${PROD_SEED:-2072601}"
mkdir -p "${ROOT}/data" "${ROOT}/logs" "${ROOT}/params"

test -x "${EXE}" || { echo "missing executable ${EXE}" >&2; exit 1; }

exec 9>"${OUTPUT}.lock"
flock -n 9 || { echo "another writer owns ${OUTPUT}" >&2; exit 0; }

resume_flag=()
if [ -s "${OUTPUT}" ] && python - "${OUTPUT}" <<'PY'
import h5py
import sys
with h5py.File(sys.argv[1], "r", locking=False) as h:
    ok = (
        ("/resume/current_step" in h or "/resume_prev/current_step" in h)
        and "/actin/center" in h
        and h["/actin/center"].shape[0] > 0
    )
raise SystemExit(0 if ok else 1)
PY
then
    echo "Found a valid checkpoint; resuming toward ${NSTEPS} total steps."
    resume_flag=(--resume=true)
elif [ -e "${OUTPUT}" ]; then
    preserved="${OUTPUT}.incomplete_job${SLURM_JOB_ID}"
    echo "Preserving non-resumable output as ${preserved}" >&2
    mv "${OUTPUT}" "${preserved}"
fi

exe_sha256=$(sha256sum "${EXE}" | awk '{print $1}')
git_head=$(git -C /home/groups/rotskoff/sherry/sarcomere/sarcomere rev-parse HEAD 2>/dev/null || echo unknown)
cat >"${PARAMS}" <<EOF
model=explicit_crossbridge_local_actn2_status2_register_production
geometry=sarcomere_length_1.93_um_x_4_periods
n_actins=1280
n_myosins=640
box=7.72x2x2
seed=${SEED}
dt=1e-6
nsteps=${NSTEPS}
save_every=${SAVE_EVERY}
pre_phase_steps=5000
interaction_ramp_steps=20000
aa_status1_register=disabled
aa_overlap_association=linear_extensive
aa_status2_registry_target_um=0.1
aa_status2_registry_on_um=0.1
aa_status2_registry_off_um=0.15
aa_register_saturation_error_um=0.05
executable=${EXE}
executable_sha256=${exe_sha256}
git_head=${git_head}
slurm_job_id=${SLURM_JOB_ID}
omp_num_threads=${OMP_NUM_THREADS}
omp_schedule=${OMP_SCHEDULE}
openblas_num_threads=${OPENBLAS_NUM_THREADS}
ld_preload=${LD_PRELOAD}
EOF

cmd=(srun --kill-on-bad-exit=1 "${EXE}"
    --n_actins=1280 --n_myosins=640
    --Lx=7.72 --Ly=2 --Lz=2 --dimension=3 --pbc_mask=111
    --actin_length=1 --myosin_length=1.5
    --actin_radius=0.015 --myosin_radius=0.025
    --initial_structure=myosin_x_noise
    --actin_direction_noise=-1
    --myosin_direction_noise=1 --myosin_direction_tilt_deg=0
    --pre_phase_steps=5000 --interaction_ramp_steps=20000
    --physical_mobility=true --temperature_kelvin=310
    --thermal_noise_scale=1 --solvent_viscosity=0.001
    --actin_hydrodynamic_diameter=0.010
    --myosin_hydrodynamic_diameter=0.035
    --actin_crowding_parallel=10 --actin_crowding_perpendicular=10
    --actin_crowding_tumble=160 --actin_crowding_roll=160
    --myosin_crowding_parallel=20 --myosin_crowding_perpendicular=20
    --myosin_crowding_tumble=64 --myosin_crowding_roll=64
    --max_actin_displacement=0.005 --max_myosin_displacement=0.005
    --max_actin_rotation=0.005 --max_myosin_rotation=0.005
    --sterics=true
    --actin_steric_k=10000 --myosin_steric_k=1000 --am_steric_k=1000
    --am_surface_gap=0.015
    --k_am=10 --kappa_am=100
    --dynamic_am=true --am_offcenter_torque=true --am_spatial_heads=true
    --am_heads_per_half=141 --am_max_heads_per_actin_half=24
    --am_k_on=100 --am_weak_binding_floor=0.1
    --am_k_weak_off=200 --am_k_powerstroke=200
    --am_k_strong_off=160 --am_force_scale=1
    --am_motor_k=1000 --am_power_stroke=0.005
    --am_stall_force=4.5 --am_unloaded_velocity=5
    --am_weak_stiffness_ratio=0.1
    --am_tether_k=3000 --am_tether_slack=0.002
    --am_alignment_k=30 --am_alignment_slack_deg=2
    --am_tether_saturation_heads=6
    --am_binding_azimuth_selection=true
    --am_binding_azimuth_half_width_deg=30
    --am_rupture_distance=0.1
    --am_strain_dependent_detachment=true
    --am_strain_detachment_theta=0.1
    --am_strain_detachment_max_exponent=20
    --am_strong_off_floor_fraction=1
    --am_cycle_period=0.05
    --am_activation_rise=0.02 --am_activation_decay=0.2
    --am_activation_floor=0 --am_activation_active_fraction=0.5
    --am_thin_activation=true --am_thin_feedback=true
    --am_q_calcium_open_rate=20 --am_q_calcium_half=0.5
    --am_q_calcium_hill=3
    --am_q_strong_feedback_rate=100 --am_q_strong_half=0.02
    --am_q_strong_hill=2 --am_q_close_rate=20
    --am_thick_feedback=false --am_status2_memory=0
    --aa_local_multilink=true
    --aa_spacing_ratio=0.92 --aa_capture_margin=0.02
    --aa_antiparallel_overlap=0.1
    --aa_axial_register_ratio=1.0
    --aa_status2_registry_on=0.1 --aa_status2_registry_off=0.15
    --aa_register_saturation_error=0.05
    --aa_prebundle_links_per_pair=2
    --aa_site_exclusion=0.0192 --aa_overlap_capture=0.1
    --k_aa=3000 --kappa_aa=50 --aa_perpendicular_k=50
    --aa_stiffness_reference_links=3 --aa_total_link_pool=0
    --aa_status_mode=arm_load_productive
    --aa_status2_force_on=1 --aa_status2_force_off=0.5
    --k_on=1000 --base_lifetime=0.1
    --directional_base_lifetime=0.05 --lifetime_coeff=0.4
    --stabilize_status1=false
    --aa_energy_gated_association=true
    --aa_transition_state_theta=0.1
    --aa_encounter_k=10 --aa_encounter_kappa=0
    --aa_passive_detailed_balance=false
    --tau_rec=0.05 --max_aa_bonds_per_actin=3
    --aa_neighbor_identity=true
    --aa_neighbor_slot_memory=0.02 --aa_status2_slot_memory=0.20
    --directional=true --max_myosin_bonds=10
    --k_mm=0 --titin_k=0 --k_bundle_max=0 --bundle_ramp_steps=0
    --dt=1e-6 --nsteps="${NSTEPS}" --save_every="${SAVE_EVERY}"
    --seed="${SEED}" --filename="${OUTPUT}")
if [ "${#resume_flag[@]}" -gt 0 ]; then
    cmd+=("${resume_flag[@]}")
fi

echo "started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf 'command='
printf '%q ' "${cmd[@]}"
printf '\n'
/usr/bin/time -v "${cmd[@]}"
echo "finished_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
