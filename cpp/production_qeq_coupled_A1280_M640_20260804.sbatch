#!/bin/bash
# Production-size 20 ms A/B gate with calcium-equilibrium q initialization,
# finite-actin axial-inline head association, restored k_align=120, and an
# optional explicit passive transverse spring on every bound A-M head.
# Molecular bonds start empty; chemistry and filament motion start together.
# No fixed-geometry chemical pre-equilibration and no automatic extension.
#SBATCH --partition=owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=8GB
#SBATCH --time=08:00:00
#SBATCH --job-name=sarc_headtether_ab_A1280
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/prod_head_tether_ab_20260804/logs/%x_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/prod_head_tether_ab_20260804/logs/%x_%j.err

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

ROOT="${PROD_ROOT:-/scratch/users/xsli/sarcomere/prod_head_tether_ab_20260804}"
EXE="${PROD_EXE:-${ROOT}/bin/sarcomere_head_tether_ab_20260804}"
EXPECTED_EXE_SHA256="${PROD_EXPECTED_EXE_SHA256:-f0cbaf635e38991601df6ae34736b518085382c99b8e46e002608efc6cd9b9ff}"
EXPECTED_AM_MODEL_VERSION="${PROD_EXPECTED_AM_MODEL_VERSION:-20}"
N_ACTINS=1280
N_MYOSINS=640
LX=7.72
SEED="${PROD_SEED:-20260804}"
NSTEPS="${PROD_NSTEPS:-40000}"
SAVE_EVERY="${PROD_SAVE_EVERY:-2000}"
MYOSIN_DIRECTION_NOISE="${PROD_MYOSIN_DIRECTION_NOISE:-1}"
AM_HEAD_TRANSVERSE_K="${PROD_AM_HEAD_TRANSVERSE_K:-0}"
AM_TRANSVERSE_STRAIN_DETACHMENT_THETA="${PROD_AM_TRANSVERSE_STRAIN_DETACHMENT_THETA:-legacy}"
MAX_NSTEPS="${PROD_MAX_NSTEPS:-40000}"
PRE_PHASE_STEPS=10000
if [ "${PROD_LOCAL_SMOKE:-0}" = 1 ]; then
    PRE_PHASE_STEPS="${PROD_SMOKE_PRE_PHASE_STEPS:-0}"
fi
TAG_SUFFIX="${PROD_TAG_SUFFIX:-}"
TAG="A1280_M640_qeq_coupled_seed${SEED}${TAG_SUFFIX}"
OUTPUT="${PROD_OUTPUT:-${ROOT}/data/${TAG}.h5}"

mkdir -p "${ROOT}/data" "${ROOT}/logs" "${ROOT}/monitor"
test -x "${EXE}" || { echo "missing executable ${EXE}" >&2; exit 1; }
exe_sha256=$(sha256sum "${EXE}" | awk '{print $1}')
if [ "${exe_sha256}" != "${EXPECTED_EXE_SHA256}" ]; then
    echo "frozen executable hash mismatch: ${exe_sha256}" >&2
    exit 1
fi
if [ "${NSTEPS}" -gt "${MAX_NSTEPS}" ]; then
    echo "gate refuses PROD_NSTEPS=${NSTEPS}; maximum is ${MAX_NSTEPS}" >&2
    exit 1
fi

python - "${N_ACTINS}" "${N_MYOSINS}" "${LX}" "${AM_HEAD_TRANSVERSE_K}" <<'PY'
import sys
n_actin, n_myosin, lx = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
head_k = float(sys.argv[4])
if n_actin != 2 * n_myosin:
    raise SystemExit(f"actin:myosin must be 2:1, got {n_actin}:{n_myosin}")
density = n_actin / (lx * 2.0 * 2.0)
if abs(density - 41.4507772) > 1e-3:
    raise SystemExit(f"unexpected actin density {density:.6f}/um^3")
if lx <= 3.68:
    raise SystemExit(f"Lx={lx} violates the minimum-image neighbour floor")
if head_k < 0.0:
    raise SystemExit(f"negative AM head transverse stiffness {head_k}")
PY

exec 9>"${OUTPUT}.lock"
flock -n 9 || { echo "another writer owns ${OUTPUT}" >&2; exit 1; }

resume_flag=()
run_nsteps="${NSTEPS}"
checkpoint_step=-1
if [ -s "${OUTPUT}" ]; then
    if checkpoint_step=$(python - "${OUTPUT}" "${AM_HEAD_TRANSVERSE_K}" \
        "${AM_TRANSVERSE_STRAIN_DETACHMENT_THETA}" \
        "${EXPECTED_AM_MODEL_VERSION}" <<'PY'
import h5py
import math
import sys

path = sys.argv[1]
expected_head_k = float(sys.argv[2])
expected_transverse_theta = sys.argv[3]
expected_model_version = int(sys.argv[4])
with h5py.File(path, "r", locking=False) as h:
    meta = h["/meta"].attrs
    def attr_text(name):
        value = meta.get(name, "")
        return value.decode() if isinstance(value, bytes) else str(value)

    checks = (
        h["/actin/center"].shape[1] == 1280,
        h["/myosin/center"].shape[1] == 640,
        int(meta.get("am_model_version", -1)) == expected_model_version,
        int(meta.get("titin_mode", -1)) == 2,
        int(meta.get("titin_schema_version", -1)) == 4,
        int(meta.get("titin_require_live_actn2", -1)) == 1,
        int(meta.get("aa_status_mode", -1)) == 6,
        int(meta.get("aa_status2_locked_until_detach", -1)) == 1,
        int(meta.get("am_heads_per_half", -1)) == 35,
        int(meta.get("am_q_initialized_at_calcium_equilibrium", -1)) == 1,
        int(meta.get("release_step_origin", -1)) == 0,
        math.isclose(float(meta.get("am_initial_q", -1)), 8.0 / 17.0),
        math.isclose(float(meta.get("aa_status1_site_diffusion", -1)), 0.010),
        math.isclose(float(meta.get("titin_k", -1)), 100.0),
        math.isclose(float(meta.get("titin_rest_length", -1)), 0.64),
        math.isclose(float(meta.get("am_alignment_k", -1)), 120.0),
        math.isclose(float(meta.get("am_head_transverse_k", -1)),
                     expected_head_k),
        math.isclose(float(meta.get("dt", -1)), 5e-7),
        attr_text("am_head_association_axial_gate") ==
            "material_head_projection_must_lie_inside_finite_actin_segment",
        attr_text("titin_association_law") ==
            "distance_limited_constant_rate_then_tension_only_spring",
    )
    if expected_transverse_theta != "legacy":
        checks += (
            math.isclose(
                float(meta.get("am_transverse_strain_detachment_theta", -1)),
                float(expected_transverse_theta)),
            attr_text("am_head_transverse_law") ==
                "per_bound_material_head_flat_bottom_transverse_spring_with_independent_escape_coupling",
        )
    if not all(checks):
        raise SystemExit("checkpoint metadata does not match q-equilibrium setup")
    candidates = []
    for group in ("/resume", "/resume_prev"):
        key = f"{group}/current_step"
        if key in h:
            candidates.append(int(h[key][...].reshape(-1)[-1]))
    if not candidates:
        raise SystemExit("no promoted resume checkpoint")
    print(max(candidates))
PY
    ); then
        if [ "${checkpoint_step}" -ge "${NSTEPS}" ]; then
            echo "validated checkpoint already reached ${checkpoint_step}/${NSTEPS}; nothing to do"
            exit 0
        fi
        echo "resuming validated checkpoint ${checkpoint_step}/${NSTEPS}"
        resume_flag=(--resume=true)
        run_nsteps="${NSTEPS}"
    else
        preserved="${OUTPUT}.nonresumable_job${SLURM_JOB_ID}_$(date -u +%Y%m%dT%H%M%SZ)"
        echo "preserving non-resumable partial output as ${preserved}" >&2
        mv "${OUTPUT}" "${preserved}"
        checkpoint_step=-1
    fi
elif [ -e "${OUTPUT}" ]; then
    preserved="${OUTPUT}.empty_job${SLURM_JOB_ID}_$(date -u +%Y%m%dT%H%M%SZ)"
    echo "preserving empty partial output as ${preserved}" >&2
    mv "${OUTPUT}" "${preserved}"
fi

launch_prefix=(srun --kill-on-bad-exit=1)
if [ "${PROD_LOCAL_SMOKE:-0}" = 1 ]; then
    launch_prefix=(env)
fi
cmd=("${launch_prefix[@]}" "${EXE}"
    --n_actins="${N_ACTINS}" --n_myosins="${N_MYOSINS}"
    --Lx="${LX}" --Ly=2 --Lz=2 --dimension=3 --pbc_mask=111
    --actin_length=1 --myosin_length=1.5
    --actin_radius=0.015 --myosin_radius=0.025
    --initial_structure=myosin_x_noise
    --actin_direction_noise=-1
    --myosin_direction_noise="${MYOSIN_DIRECTION_NOISE}" --myosin_direction_tilt_deg=0
    --pre_phase_steps="${PRE_PHASE_STEPS}" --interaction_ramp_steps=0
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
    --am_head_transverse_k="${AM_HEAD_TRANSVERSE_K}"
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
    --aa_status2_registry_on=0.020 --aa_status2_registry_off=0.040
    --aa_status2_support_window=0.010
    --aa_register_saturation_error=0.05
    --aa_prebundle_links_per_pair=3
    --aa_site_exclusion=0.0192 --aa_overlap_capture=0.10
    --k_aa=1000 --kappa_aa=100 --aa_perpendicular_k=50
    --aa_status1_radial_k=75 --aa_status1_site_diffusion=0.010
    --aa_status2_locked_until_detach=true
    --aa_stiffness_reference_links=3 --aa_total_link_pool=0
    --aa_status_mode=arm_load_productive
    --aa_status2_force_on=0.25 --aa_status2_force_off=0.10
    --k_on=1000 --k_off=1 --base_lifetime=0.1
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
    --titin_k=100 --titin_rest_length=0.64
    --titin_links_per_half=3 --titin_capture_distance=0.80
    --titin_k_on=5000 --titin_k_off=0.05
    --titin_strain_gated_kinetics=true --titin_transition_theta=0.95
    --titin_rebind_memory=0.020 --titin_rebind_rate_multiplier=100
    --titin_require_live_actn2=true
    --titin_k_off_unsupported=0 --titin_rebind_requires_actn2=false
    --k_mm=0 --k_bundle_max=0 --bundle_ramp_steps=0
    --dt=5e-7 --nsteps="${run_nsteps}" --save_every="${SAVE_EVERY}"
    --seed="${SEED}" --filename="${OUTPUT}")
if [ "${AM_TRANSVERSE_STRAIN_DETACHMENT_THETA}" != legacy ]; then
    cmd+=(--am_transverse_strain_detachment_theta="${AM_TRANSVERSE_STRAIN_DETACHMENT_THETA}")
fi
if [ "${#resume_flag[@]}" -gt 0 ]; then
    cmd+=("${resume_flag[@]}")
fi

echo "protocol=q_equilibrium_coupled_chemistry_motion"
echo "initial_q=8/17"
echo "myosin_direction_noise=${MYOSIN_DIRECTION_NOISE}"
echo "am_alignment_k=120"
echo "am_head_transverse_k=${AM_HEAD_TRANSVERSE_K}"
echo "am_transverse_strain_detachment_theta=${AM_TRANSVERSE_STRAIN_DETACHMENT_THETA}"
echo "pre_phase_steps=${PRE_PHASE_STEPS}"
echo "fixed_lattice_preequil_steps=0"
echo "am_head_association_axial_gate=material_head_projection_must_lie_inside_finite_actin_segment"
echo "titin_association_law=distance_limited_constant_rate_then_tension_only_spring"
echo "frozen_executable=${EXE}"
echo "executable_sha256=${exe_sha256}"
echo "resume_checkpoint_step=${checkpoint_step}"
echo "started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf 'command='
printf '%q ' "${cmd[@]}"
printf '\n'
/usr/bin/time -v "${cmd[@]}"
echo "finished_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
