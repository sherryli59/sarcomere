#!/bin/bash
# Production-size gate for the exact molecular mechanics validated in the
# mobile 2-myosin/4-actin periodic cell.  This stage releases a randomly
# initialized 1,280-actin assembly for 20 ms after 30 ms of fixed-geometry
# chemistry. It never auto-extends; the longer continuation is submitted only
# after this checkpoint passes numerical, occupancy, and structural gates.
#SBATCH --partition=owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64GB
#SBATCH --time=08:00:00
#SBATCH --job-name=sarc_exact_A1280_gate
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/prod_exact_setup_20260804/logs/%x_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/prod_exact_setup_20260804/logs/%x_%j.err

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

ROOT="${PROD_ROOT:-/scratch/users/xsli/sarcomere/prod_exact_setup_20260804}"
EXE="${PROD_EXE:-${ROOT}/bin/sarcomere_exact_20260804}"
EXPECTED_EXE_SHA256=6a4dda5eb1858fc304dba09560871a02686db0b0214d3cc4f6d0c4e88dd04395
N_ACTINS=1280
N_MYOSINS=640
LX=7.72
SEED="${PROD_SEED:-20260804}"
NSTEPS="${PROD_NSTEPS:-40000}"
SAVE_EVERY="${PROD_SAVE_EVERY:-2000}"
FIXED_PREEQ_STEPS=60000
RESUME_TARGET=$((FIXED_PREEQ_STEPS + NSTEPS))
TAG="A1280_M640_exact_seed${SEED}"
OUTPUT="${PROD_OUTPUT:-${ROOT}/data/${TAG}.h5}"
PARAMS="${PROD_PARAMS:-${ROOT}/params/${TAG}.params}"
mkdir -p "${ROOT}/data" "${ROOT}/logs" "${ROOT}/params" "${ROOT}/monitor"

test -x "${EXE}" || { echo "missing executable ${EXE}" >&2; exit 1; }
exe_sha256=$(sha256sum "${EXE}" | awk '{print $1}')
if [ "${exe_sha256}" != "${EXPECTED_EXE_SHA256}" ]; then
    echo "frozen executable hash mismatch: ${exe_sha256}" >&2
    exit 1
fi
if [ "${NSTEPS}" -gt 40000 ]; then
    echo "gate refuses PROD_NSTEPS=${NSTEPS}; maximum is 40000" >&2
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
PY

exec 9>"${OUTPUT}.lock"
flock -n 9 || { echo "another writer owns ${OUTPUT}" >&2; exit 1; }

# Auto-requeue is deliberately disabled at the Slurm level. A replacement job
# may resume only from a validated, fully promoted simulator checkpoint. The
# global KMC counter includes the 60,000 fixed-lattice chemistry updates, so a
# resumed Langevin run targets 60,000 + NSTEPS rather than NSTEPS alone.
resume_flag=()
run_nsteps="${NSTEPS}"
checkpoint_step=-1
if [ -s "${OUTPUT}" ]; then
    if checkpoint_step=$(python - "${OUTPUT}" <<'PY'
import h5py
import math
import sys

path = sys.argv[1]
with h5py.File(path, "r", locking=False) as h:
    meta = h["/meta"].attrs
    checks = (
        h["/actin/center"].shape[1] == 1280,
        h["/myosin/center"].shape[1] == 640,
        int(meta.get("titin_mode", -1)) == 2,
        int(meta.get("titin_schema_version", -1)) == 4,
        int(meta.get("titin_require_live_actn2", -1)) == 1,
        int(meta.get("aa_status_mode", -1)) == 6,
        int(meta.get("aa_status2_locked_until_detach", -1)) == 1,
        int(meta.get("am_heads_per_half", -1)) == 35,
        math.isclose(float(meta.get("aa_status1_site_diffusion", -1)), 0.010),
        math.isclose(float(meta.get("titin_k", -1)), 100.0),
        math.isclose(float(meta.get("titin_rest_length", -1)), 0.64),
        math.isclose(float(meta.get("dt", -1)), 5e-7),
    )
    if not all(checks):
        raise SystemExit("checkpoint metadata does not match exact setup")
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
        if [ "${checkpoint_step}" -ge "${RESUME_TARGET}" ]; then
            echo "validated checkpoint already reached ${checkpoint_step}/${RESUME_TARGET}; nothing to do"
            exit 0
        fi
        echo "resuming validated checkpoint ${checkpoint_step}/${RESUME_TARGET}"
        resume_flag=(--resume=true)
        run_nsteps="${RESUME_TARGET}"
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

git_head=$(git -C /home/groups/rotskoff/sherry/sarcomere/sarcomere rev-parse HEAD 2>/dev/null || echo unknown)
cat >"${PARAMS}" <<EOF
model=exact_mobile_balanced_cell_mechanics_scaled_to_random_assembly
stage=production_size_20ms_gate_no_auto_extension
n_actins=${N_ACTINS}
n_myosins=${N_MYOSINS}
box_um=${LX}x2x2_periodic
seed=${SEED}
dt_s=5e-7
fixed_lattice_chemistry_steps=${FIXED_PREEQ_STEPS}
fixed_lattice_chemistry_ms=30
released_steps=${NSTEPS}
released_ms=$(python -c "print(${NSTEPS}*5e-4)")
status1_actn2_site_diffusion_um2_per_s=0.010
status2_actn2_locked_until_molecular_detachment=true
actn2_molecules_per_pair=3
actn2_target_overlap_um=0.10
actn2_barbed_zone_um=0.20
actn2_status1_capture_centerline_um=0.10
actn2_status1_radial_k_pN_per_um=75
actn2_status2_preferred_centerline_um=0.055
actn2_contour_length_um=0.12
actn2_status_mode=arm_load_productive
titin_topology=explicit_thick_tip_domains_to_individual_actin_barbed_domains
titin_domains_per_half=3
titin_slack_length_um=0.64
titin_capture_distance_um=0.80
titin_k_pN_per_um=100
titin_k_on_per_s=5000
titin_k_off_per_s=0.05
titin_transition_theta=0.95
titin_new_capture_requires_live_actn2=true
titin_unsupported_forced_release=false
titin_force_law=tension_only_zero_under_compression
am_heads_per_half=35
am_activation=tonic
am_alignment_k_pN_um_per_rad2=30
thermal_translation_and_rotation=true
temperature_kelvin=310
preexisting_z_disc=false
centroid_or_z_node_forces=false
external_load=false
myomesin_or_bundle_spring=false
executable=${EXE}
executable_sha256=${exe_sha256}
git_head=${git_head}
slurm_job_id=${SLURM_JOB_ID}
resume_checkpoint_step=${checkpoint_step}
simulator_nsteps_argument=${run_nsteps}
EOF

cmd=(srun --kill-on-bad-exit=1 "${EXE}"
    --n_actins="${N_ACTINS}" --n_myosins="${N_MYOSINS}"
    --Lx="${LX}" --Ly=2 --Lz=2 --dimension=3 --pbc_mask=111
    --actin_length=1 --myosin_length=1.5
    --actin_radius=0.015 --myosin_radius=0.025
    --initial_structure=myosin_x_noise
    --actin_direction_noise=-1
    --myosin_direction_noise=1 --myosin_direction_tilt_deg=0
    --pre_phase_steps=10000 --interaction_ramp_steps=0
    --fixed_lattice_preequil_steps="${FIXED_PREEQ_STEPS}"
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
    --am_alignment_k=30 --am_alignment_slack_deg=2
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
if [ "${#resume_flag[@]}" -gt 0 ]; then
    cmd+=("${resume_flag[@]}")
fi

echo "started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf 'command='
printf '%q ' "${cmd[@]}"
printf '\n'
/usr/bin/time -v "${cmd[@]}"
echo "finished_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
