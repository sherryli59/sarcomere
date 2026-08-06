#!/bin/bash
# Fresh 1,280-actin production run of the configuration that passed the
# one-second minimal balanced gate. No external load, fixed lattice, Z node,
# centroid force, automatic extension, or scheduler auto-requeue is used.
# A sterics-only placement relaxation precedes time zero; chemistry and motion
# then start together with q at calcium equilibrium and a 20 ms force ramp.
#SBATCH --partition=owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64GB
#SBATCH --time=08:00:00
#SBATCH --job-name=sarc_pr_qeq_A1280_100ms
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/prod_productive_register_qeq_ramp_20260805/logs/%x_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/prod_productive_register_qeq_ramp_20260805/logs/%x_%j.err

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}"
export OMP_DYNAMIC=FALSE
export OMP_SCHEDULE=static
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export MKL_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

JEMALLOC=/home/groups/rotskoff/sherry/mamba_root/envs/sarcomere/lib/libjemalloc.so
test -r "${JEMALLOC}" || { echo "missing jemalloc: ${JEMALLOC}" >&2; exit 1; }
export LD_PRELOAD="${JEMALLOC}"

ROOT="${PROD_ROOT:-/scratch/users/xsli/sarcomere/prod_productive_register_qeq_ramp_20260805}"
EXE=/scratch/users/xsli/sarcomere/prod_productive_register_qeq_ramp_20260805/bin/sarcomere_productive_register_qeq_ramp_20260805
EXPECTED_EXE_SHA256=bcfaee84bfc9eaf137ce18a314faea142533f2c3660422f6c11580627fa61369
N_ACTINS=1280
N_MYOSINS=640
LX=7.72
SEED="${PROD_SEED:-20260805}"
NSTEPS="${PROD_NSTEPS:-200000}"
SAVE_EVERY=2000
PRE_PHASE_STEPS=10000
RUN_LABEL="${PROD_RUN_LABEL:-100ms}"
TAG="${PROD_TAG:-A1280_M640_productive_register_qeq_ramp_seed${SEED}_${RUN_LABEL}}"
OUTPUT=${ROOT}/data/${TAG}.h5
PARAMS=${ROOT}/params/${TAG}.params
AM_TRANSVERSE_THETA="${PROD_AM_TRANSVERSE_THETA:-0}"
TITIN_K="${PROD_TITIN_K:-100}"
TITIN_REST_LENGTH="${PROD_TITIN_REST_LENGTH:-0.64}"
TITIN_CAPTURE_DISTANCE="${PROD_TITIN_CAPTURE_DISTANCE:-0.80}"
TITIN_K_OFF="${PROD_TITIN_K_OFF:-0.05}"
TITIN_TRANSITION_THETA="${PROD_TITIN_TRANSITION_THETA:-0.95}"
TITIN_K_OFF_UNSUPPORTED="${PROD_TITIN_K_OFF_UNSUPPORTED:-0}"
TITIN_REBIND_REQUIRES_ACTN2="${PROD_TITIN_REBIND_REQUIRES_ACTN2:-false}"

mkdir -p "${ROOT}"/{data,logs,params,monitor}
test -x "${EXE}" || { echo "missing executable ${EXE}" >&2; exit 1; }
exe_sha256=$(sha256sum "${EXE}" | awk '{print $1}')
if [ "${exe_sha256}" != "${EXPECTED_EXE_SHA256}" ]; then
    echo "frozen executable hash mismatch: ${exe_sha256}" >&2
    exit 1
fi

python - "${N_ACTINS}" "${N_MYOSINS}" "${LX}" <<'PY'
import sys
n_actin, n_myosin, lx = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
if n_actin != 2 * n_myosin:
    raise SystemExit(f"actin:myosin must be 2:1, got {n_actin}:{n_myosin}")
density = n_actin / (lx * 2.0 * 2.0)
if abs(density - 41.4507772) > 1e-3:
    raise SystemExit(f"unexpected actin density {density:.6f}/um^3")
if lx <= 3.68:
    raise SystemExit(f"Lx={lx} violates the minimum-image neighbour floor")
# The validated 55 nm A-M centerline target is exactly the sum of the
# 25 nm myosin radius, 15 nm actin radius, and 15 nm surface gap.
if abs((0.025 + 0.015 + 0.015) - 0.055) > 1e-12:
    raise SystemExit("inconsistent steric and preferred A-M spacing")
PY

exec 9>"${OUTPUT}.lock"
flock -n 9 || { echo "another writer owns ${OUTPUT}" >&2; exit 1; }

resume_flag=()
checkpoint_step=-1
if [ -s "${OUTPUT}" ]; then
    checkpoint_step=$(python - "${OUTPUT}" "${AM_TRANSVERSE_THETA}" \
        "${TITIN_K}" "${TITIN_REST_LENGTH}" "${TITIN_CAPTURE_DISTANCE}" \
        "${TITIN_K_OFF}" "${TITIN_TRANSITION_THETA}" \
        "${TITIN_K_OFF_UNSUPPORTED}" "${TITIN_REBIND_REQUIRES_ACTN2}" <<'PY'
import h5py
import math
import sys

path = sys.argv[1]
am_transverse_theta = float(sys.argv[2])
titin_k = float(sys.argv[3])
titin_rest = float(sys.argv[4])
titin_capture = float(sys.argv[5])
titin_k_off = float(sys.argv[6])
titin_transition_theta = float(sys.argv[7])
titin_k_off_unsupported = float(sys.argv[8])
titin_rebind_requires_actn2 = sys.argv[9].lower() == "true"
with h5py.File(path, "r", locking=False) as h:
    meta = h["/meta"].attrs
    def text(name):
        value = meta.get(name, "")
        return value.decode() if isinstance(value, bytes) else str(value)
    checks = (
        h["/actin/center"].shape[1] == 1280,
        h["/myosin/center"].shape[1] == 640,
        int(meta.get("am_model_version", -1)) == 21,
        int(meta.get("am_q_initialized_at_calcium_equilibrium", -1)) == 1,
        math.isclose(float(meta.get("am_initial_q", -1)), 8.0 / 17.0),
        int(meta.get("interaction_ramp_steps", -1)) == 40000,
        int(meta.get("aa_status_mode", -1)) == 7,
        int(meta.get("aa_status2_locked_until_detach", -1)) == 1,
        math.isclose(float(meta.get("aa_status2_registry_on", -1)), 0.040),
        math.isclose(float(meta.get("aa_status2_registry_off", -1)), 0.040),
        math.isclose(float(meta.get("aa_status2_force_on", -1)), 0.25),
        math.isclose(float(meta.get("aa_status1_site_diffusion", -1)), 0.010),
        math.isclose(float(meta.get("myosin_radius", -1)), 0.025),
        math.isclose(float(meta.get("actin_radius", -1)), 0.015),
        math.isclose(float(meta.get("am_surface_gap", -1)), 0.015),
        math.isclose(float(meta.get("am_head_transverse_k", -1)), 100.0),
        math.isclose(float(meta.get("am_transverse_strain_detachment_theta", -1)),
                     am_transverse_theta),
        math.isclose(float(meta.get("am_alignment_k", -1)), 120.0),
        int(meta.get("titin_mode", -1)) == 2,
        int(meta.get("titin_require_live_actn2", -1)) == 1,
        math.isclose(float(meta.get("titin_k", -1)), titin_k),
        math.isclose(float(meta.get("titin_rest_length", -1)), titin_rest),
        math.isclose(float(meta.get("titin_capture_distance", -1)), titin_capture),
        math.isclose(float(meta.get("titin_k_off", -1)), titin_k_off),
        math.isclose(float(meta.get("titin_transition_theta", -1)),
                     titin_transition_theta),
        math.isclose(float(meta.get("titin_k_off_unsupported", -1)),
                     titin_k_off_unsupported),
        int(meta.get("titin_rebind_requires_actn2", -1)) ==
            int(titin_rebind_requires_actn2),
        math.isclose(float(meta.get("dt", -1)), 5e-7),
        text("am_head_association_axial_gate") ==
            "material_head_projection_must_lie_inside_finite_actin_segment",
        text("titin_association_law") ==
            "distance_limited_constant_rate_then_tension_only_spring",
    )
    if not all(checks):
        raise SystemExit("checkpoint metadata does not match the validated production setup")
    steps = []
    for group in ("/resume", "/resume_prev"):
        key = f"{group}/current_step"
        if key in h:
            steps.append(int(h[key][...].reshape(-1)[-1]))
    if not steps:
        raise SystemExit("no promoted resume checkpoint")
    print(max(steps))
PY
    )
    if [ "${checkpoint_step}" -ge "${NSTEPS}" ]; then
        echo "validated checkpoint already reached ${checkpoint_step}/${NSTEPS}; nothing to do"
        exit 0
    fi
    echo "resuming validated checkpoint ${checkpoint_step}/${NSTEPS}"
    resume_flag=(--resume=true)
elif [ -e "${OUTPUT}" ]; then
    preserved="${OUTPUT}.empty_job${SLURM_JOB_ID}_$(date -u +%Y%m%dT%H%M%SZ)"
    mv "${OUTPUT}" "${preserved}"
fi

git_head=$(git -C /home/groups/rotskoff/sherry/sarcomere/sarcomere rev-parse HEAD 2>/dev/null || echo unknown)
cat >"${PARAMS}" <<EOF
model=productive_register_qeq_ramp_model21
stage=production_1280_actin_100ms_no_auto_extension
n_actins=${N_ACTINS}
n_myosins=${N_MYOSINS}
box_um=${LX}x2x2_periodic
seed=${SEED}
dt_s=5e-7
nsteps=${NSTEPS}
save_every=${SAVE_EVERY}
sterics_only_prephase_steps=${PRE_PHASE_STEPS}
fixed_lattice_chemistry_steps=0
chemistry_and_motion_start_together=true
interaction_ramp_steps=40000
interaction_ramp_ms=20
q_initialized_at_calcium_equilibrium=true
initial_q=8/17
actin_radius_um=0.015
myosin_radius_um=0.025
am_surface_gap_um=0.015
am_preferred_centerline_um=0.055
am_spacing_identity=0.025+0.015+0.015=0.055
am_heads_per_half=35
am_head_transverse_k=100
am_transverse_strain_detachment_theta=${AM_TRANSVERSE_THETA}
am_alignment_k=120
actn2_status_mode=productive_register
actn2_formation_registry_um=0.040
actn2_productive_support_threshold_pN=0.25
actn2_k_on_per_s=2500
actn2_status1_site_diffusion_um2_per_s=0.010
actn2_target_overlap_um=0.10
actn2_register_all_statuses=false
titin_k_pN_per_um=${TITIN_K}
titin_slack_length_um=${TITIN_REST_LENGTH}
titin_capture_distance_um=${TITIN_CAPTURE_DISTANCE}
titin_k_off_per_s=${TITIN_K_OFF}
titin_transition_theta=${TITIN_TRANSITION_THETA}
titin_k_off_unsupported_per_s=${TITIN_K_OFF_UNSUPPORTED}
titin_rebind_requires_actn2=${TITIN_REBIND_REQUIRES_ACTN2}
external_load=false
preexisting_z_disc=false
centroid_or_z_node_forces=false
automatic_extension=false
executable=${EXE}
executable_sha256=${exe_sha256}
git_head=${git_head}
slurm_job_id=${SLURM_JOB_ID}
resume_checkpoint_step=${checkpoint_step}
EOF

cmd=(srun --kill-on-bad-exit=1 "${EXE}"
    --n_actins=${N_ACTINS} --n_myosins=${N_MYOSINS}
    --Lx=${LX} --Ly=2 --Lz=2 --dimension=3 --pbc_mask=111
    --actin_length=1 --myosin_length=1.5
    --actin_radius=0.015 --myosin_radius=0.025
    --initial_structure=myosin_x_noise
    --actin_direction_noise=-1
    --myosin_direction_noise=1 --myosin_direction_tilt_deg=0
    --pre_phase_steps=${PRE_PHASE_STEPS} --interaction_ramp_steps=40000
    --fixed_lattice_preequil_steps=0
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
    --actin_steric_k=2000 --myosin_steric_k=1000 --am_steric_k=1000
    --am_cutoff=0.10 --am_optimal=0.055 --am_surface_gap=0.015
    --k_am=0 --kappa_am=0
    --dynamic_am=true --am_offcenter_torque=true --am_spatial_heads=true
    --am_heads_per_half=35 --am_max_heads_per_actin_half=24
    --am_k_on=100 --am_weak_binding_floor=0.1
    --am_k_weak_off=200 --am_k_powerstroke=200
    --am_k_strong_off=160 --am_force_scale=1
    --am_motor_k=1000 --am_power_stroke=0.005
    --am_stall_force=4.5 --am_unloaded_velocity=5
    --am_weak_stiffness_ratio=0.1
    --am_tether_k=3000 --am_tether_slack=0.002
    --am_head_transverse_k=100
    --am_transverse_strain_detachment_theta=${AM_TRANSVERSE_THETA}
    --am_alignment_k=120 --am_alignment_slack_deg=2
    --am_tether_saturation_heads=6
    --am_binding_azimuth_selection=true
    --am_binding_azimuth_half_width_deg=30
    --am_rupture_distance=0.1
    --am_strain_dependent_detachment=true
    --am_strain_detachment_theta=0.1
    --am_strain_detachment_max_exponent=20
    --am_strong_off_floor_fraction=1
    --am_cycle_period=0 --am_activation_floor=0
    --am_activation_active_fraction=0.5
    --am_thin_activation=true --am_thin_feedback=true
    --am_q_calcium_open_rate=20 --am_q_calcium_half=0.5
    --am_q_calcium_hill=3
    --am_q_strong_feedback_rate=100 --am_q_strong_half=0.02
    --am_q_strong_hill=2 --am_q_close_rate=20
    --am_initialize_q_calcium_equilibrium=true
    --am_thick_feedback=false --am_status2_memory=0
    --am_bidirectional_crossbridge=false --am_overlap_target_ratio=0
    --aa_local_multilink=true
    --aa_spacing_ratio=1 --aa_preferred_centerline=0.055
    --aa_capture_margin=0.045
    --aa_antiparallel_overlap=0.10 --aa_barbed_zone=0.20
    --aa_axial_register_ratio=3 --aa_axial_register_force_max=36
    --aa_register_all_statuses=false
    --aa_status2_registry_on=0.040 --aa_status2_registry_off=0.040
    --aa_status2_support_window=0.010
    --aa_register_saturation_error=0.05
    --aa_prebundle_links_per_pair=3
    --aa_site_exclusion=0.0192 --aa_overlap_capture=0.10
    --k_aa=1000 --kappa_aa=100 --aa_perpendicular_k=50
    --aa_status1_radial_k=75 --aa_status1_site_diffusion=0.010
    --aa_status2_locked_until_detach=true
    --aa_stiffness_reference_links=3 --aa_total_link_pool=0
    --aa_status_mode=productive_register
    --aa_status2_force_on=0.25 --aa_status2_force_off=0.10
    --k_on=2500 --k_off=1 --base_lifetime=0.1
    --directional_base_lifetime=0.05 --lifetime_coeff=0.4
    --stabilize_status1=false
    --aa_energy_gated_association=false
    --aa_transition_state_theta=0.1
    --aa_encounter_k=10 --aa_encounter_kappa=0
    --aa_passive_detailed_balance=false
    --tau_rec=0 --max_aa_bonds_per_actin=2
    --aa_neighbor_identity=true
    --aa_neighbor_slot_memory=0.02 --aa_status2_slot_memory=0.20
    --aa_association_max_direction_dot=-0.5
    --aa_link_rupture_length=0.12
    --directional=true --max_myosin_bonds=24
    --titin_mode=explicit_dynamic
    --titin_k=${TITIN_K} --titin_rest_length=${TITIN_REST_LENGTH}
    --titin_links_per_half=3 --titin_capture_distance=${TITIN_CAPTURE_DISTANCE}
    --titin_k_on=5000 --titin_k_off=${TITIN_K_OFF}
    --titin_strain_gated_kinetics=true --titin_transition_theta=${TITIN_TRANSITION_THETA}
    --titin_rebind_memory=0.020 --titin_rebind_rate_multiplier=100
    --titin_require_live_actn2=true
    --titin_k_off_unsupported=${TITIN_K_OFF_UNSUPPORTED}
    --titin_rebind_requires_actn2=${TITIN_REBIND_REQUIRES_ACTN2}
    --k_mm=0 --k_bundle_max=0 --bundle_ramp_steps=0
    --dt=5e-7 --nsteps=${NSTEPS} --save_every=${SAVE_EVERY}
    --seed=${SEED} --filename=${OUTPUT})
if [ "${#resume_flag[@]}" -gt 0 ]; then
    cmd+=("${resume_flag[@]}")
fi

echo "protocol=productive_register_qeq_20ms_ramp"
echo "frozen_executable=${EXE}"
echo "executable_sha256=${exe_sha256}"
echo "resume_checkpoint_step=${checkpoint_step}"
echo "started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf 'command='
printf '%q ' "${cmd[@]}"
printf '\n'
/usr/bin/time -v "${cmd[@]}"
echo "finished_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
