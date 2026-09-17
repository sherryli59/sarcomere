#!/bin/bash
# Guarded half-size, production-like assembly-from-scratch pilot.
#
# This preserves the A1280/M640 density and both 2-um transverse dimensions,
# while shortening only Lx.  The executable is frozen and hash-checked.  The
# output namespace is fresh, resume is forbidden, and Slurm requeue/automatic
# extension are disabled.
#SBATCH --partition=rotskoff
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64GB
#SBATCH --time=16:00:00
#SBATCH --job-name=sarc_titin457_A800_M320_5s
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/prod_titin457_A800_M320_5s_20260807/logs/%x_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/prod_titin457_A800_M320_5s_20260807/logs/%x_%j.err

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

ROOT=/scratch/users/xsli/sarcomere/prod_titin457_A800_M320_5s_20260807
EXE=${ROOT}/bin/sarcomere_dt1e5_20260807
EXPECTED_EXE_SHA256=08750ac02736f43448f35d7e34891c11b466796784826207ad1ced487ef06810
N_ACTINS=800
N_MYOSINS=320
LX=4.825
LY=2
LZ=2
SEED="${PROD_SEED:-20260805}"
NSTEPS=500000
SAVE_EVERY=500
PRE_PHASE_STEPS=500
TAG="A800_M320_titin457_seed${SEED}_5s_fresh"
OUTPUT=${ROOT}/data/${TAG}.h5
PARAMS=${ROOT}/params/${TAG}.params

# Latest explicit-titin sandbox values from the validated bilateral rig.
TITIN_K=100
TITIN_REST_LENGTH=0.457
TITIN_CAPTURE_DISTANCE=0.72
TITIN_K_ON=5000
TITIN_K_OFF=0.01
TITIN_TRANSITION_THETA=1.0

mkdir -p "${ROOT}"/{data,logs,params,monitor}
test -x "${EXE}" || { echo "missing executable: ${EXE}" >&2; exit 1; }
exe_sha256=$(sha256sum "${EXE}" | awk '{print $1}')
if [ "${exe_sha256}" != "${EXPECTED_EXE_SHA256}" ]; then
    echo "frozen executable hash mismatch: ${exe_sha256}" >&2
    exit 1
fi

python - "${N_ACTINS}" "${N_MYOSINS}" "${LX}" "${LY}" "${LZ}" <<'PY'
import sys
n_actin, n_myosin = map(int, sys.argv[1:3])
lx, ly, lz = map(float, sys.argv[3:6])
if (n_actin, n_myosin) != (800, 320):
    raise SystemExit(f"expected 800:320 for this low-drive pilot, got {n_actin}:{n_myosin}")
density = n_actin / (lx * ly * lz)
if abs(density - 41.4507772) > 1e-3:
    raise SystemExit(f"unexpected actin density {density:.7f}/um^3")
if lx <= 3.68:
    raise SystemExit(f"Lx={lx} violates the conservative minimum-image floor")
if ly != 2.0 or lz != 2.0:
    raise SystemExit("transverse box must remain 2 x 2 um for this pilot")
PY

if [ -e "${OUTPUT}" ]; then
    echo "refusing to resume or overwrite existing output: ${OUTPUT}" >&2
    exit 1
fi
exec 9>"${OUTPUT}.lock"
flock -n 9 || { echo "another writer owns ${OUTPUT}" >&2; exit 1; }

git_head=$(git -C /home/groups/rotskoff/sherry/sarcomere/sarcomere rev-parse HEAD 2>/dev/null || echo unknown)
cat >"${PARAMS}" <<EOF
model=current_sandbox_model22_explicit_bilateral_titin_global_actin_sites
stage=titin_rest_length_commensurate_sarcomere_A800_M320_5s
titin_rest_length_um=0.457
rationale=myosin_half_0.750_plus_titin_rest_sets_the_half_sarcomere
old_rest_0.665_implies_2.83um_sarcomere_which_Lx_4.825_cannot_host
new_rest_0.457_implies_2.413um_equals_Lx_over_2
baseline_defect=75.6pct_of_live_titin_slack__8.75pct_of_myosins_bilaterally_taut
bilaterally_taut_myosins_reach_Mband_order_0.818_vs_0.177_overall
control_arm=job_38011216_titin_rest_0.665
single_variable_change=titin_rest_length_only
n_actins=${N_ACTINS}
n_myosins=${N_MYOSINS}
box_um=${LX}x${LY}x${LZ}_periodic
actin_density_per_um3=41.4507772
seed=${SEED}
dt_s=1e-5
nsteps=${NSTEPS}
duration_ms=5000
save_every_steps=${SAVE_EVERY}
sterics_only_prephase_steps=${PRE_PHASE_STEPS}
sterics_only_prephase_ms=5
fixed_lattice_chemistry_steps=0
chemistry_and_motion_start_together=true
interaction_ramp_steps=2000
interaction_ramp_ms=20
q_initialized_at_calcium_equilibrium=true
initial_q=8/17
am_head_association_axial_gate=material_head_projection_must_lie_inside_finite_actin_segment
am_sites_per_actin=182
am_site_spacing_nm=5.49
am_global_head_cap_per_actin=24
am_neighbor_pair_cap=overlap_scaled_max_24
am_duplicate_site_ownership=forbidden_globally_across_myosins
am_heads_per_half=9
am_crossover_pull_through=true
am_head_transverse_k=100
am_transverse_strain_detachment_theta=0
am_alignment_k=120
actn2_status_mode=productive_register
actn2_k_on_per_s=2500
actn2_target_overlap_um=0.10
actn2_formation_register_um=0.040
actn2_status1_site_diffusion_um2_per_s=0.010
titin_mode=explicit_dynamic_per_half
titin_links_per_myosin_half=3
titin_k_pN_per_um=${TITIN_K}
titin_slack_length_um=${TITIN_REST_LENGTH}
titin_capture_distance_um=${TITIN_CAPTURE_DISTANCE}
titin_k_on_per_s=${TITIN_K_ON}
titin_k_off_per_s=${TITIN_K_OFF}
titin_transition_theta=${TITIN_TRANSITION_THETA}
titin_new_partner_requires_live_actn2=true
titin_same_partner_rebind_memory_ms=20
titin_same_partner_rebind_multiplier=100
compression_bearing_titin=false
short_length_penalty=false
pointed_end_contact=false
external_load=false
preexisting_z_disc=false
centroid_or_z_node_forces=false
automatic_extension=false
resume=false
scheduler_requeue=false
executable=${EXE}
executable_sha256=${exe_sha256}
git_head=${git_head}
source_tracked_diff_sha256=1b7253c6e6fe3f96c900be8a694b5b142956d1277544eff8152507e54a60e026
slurm_job_id=${SLURM_JOB_ID}
EOF

cmd=(srun --kill-on-bad-exit=1 "${EXE}"
    --n_actins=${N_ACTINS} --n_myosins=${N_MYOSINS}
    --Lx=${LX} --Ly=${LY} --Lz=${LZ} --dimension=3 --pbc_mask=111
    --actin_length=1 --myosin_length=1.5
    --actin_radius=0.015 --myosin_radius=0.025
    --initial_structure=myosin_x_noise
    --actin_direction_noise=-1
    --myosin_direction_noise=1 --myosin_direction_tilt_deg=0
    --pre_phase_steps=${PRE_PHASE_STEPS} --interaction_ramp_steps=2000
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
    --am_cutoff=0.05 --am_optimal=0.055 --am_surface_gap=0.015
    --k_am=0 --kappa_am=0
    --dynamic_am=true --am_offcenter_torque=true --am_spatial_heads=true
    --am_heads_per_half=9 --am_crossover_pull_through=true
    --am_max_heads_per_actin_half=24
    --am_sites_per_actin=182
    --am_max_heads_per_actin=24
    --am_k_on=100 --am_weak_binding_floor=0.1
    --am_k_weak_off=200 --am_k_powerstroke=200
    --am_k_strong_off=160 --am_force_scale=1
    --am_motor_k=1000 --am_power_stroke=0.005
    --am_stall_force=4.5 --am_unloaded_velocity=5
    --am_weak_stiffness_ratio=0.1
    --am_tether_k=3000 --am_tether_slack=0.002
    --am_head_transverse_k=100
    --am_transverse_strain_detachment_theta=0
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
    --aa_double_overlap_contact=true --aa_double_overlap_contact_diameter=0.012 --aa_double_overlap_contact_epsilon=0.00428 --aa_double_overlap_contact_force_cap=0
    --directional=true --max_myosin_bonds=24
    --titin_mode=explicit_dynamic
    --titin_k=${TITIN_K} --titin_rest_length=${TITIN_REST_LENGTH}
    --titin_links_per_half=3
    --titin_capture_distance=${TITIN_CAPTURE_DISTANCE}
    --titin_k_on=${TITIN_K_ON} --titin_k_off=${TITIN_K_OFF}
    --titin_strain_gated_kinetics=true
    --titin_transition_theta=${TITIN_TRANSITION_THETA}
    --titin_rebind_memory=0.020 --titin_rebind_rate_multiplier=100
    --titin_require_live_actn2=true
    --titin_k_off_unsupported=0
    --titin_rebind_requires_actn2=false
    --k_mm=0 --k_bundle_max=0 --bundle_ramp_steps=0
    --dt=1e-5 --nsteps=${NSTEPS} --save_every=${SAVE_EVERY}
    --seed=${SEED} --filename=${OUTPUT})

echo "protocol=small_bilateral_sites_assembly_from_scratch_20ms"
echo "titin_rest_length=${TITIN_REST_LENGTH} (baseline 0.665)"
echo "frozen_executable=${EXE}"
echo "executable_sha256=${exe_sha256}"
echo "started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf 'command='
printf '%q ' "${cmd[@]}"
printf '\n'
/usr/bin/time -v "${cmd[@]}"

# Fail loudly if the completed file does not encode the intended model.
python - "${OUTPUT}" <<'PY'
import h5py
import math
import sys

path = sys.argv[1]
with h5py.File(path, "r", locking=False) as h:
    m = h["/meta"].attrs
    def txt(key):
        value = m.get(key, "")
        return value.decode() if isinstance(value, bytes) else str(value)
    checks = {
        "actin_count": h["/actin/center"].shape[1] == 800,
        "myosin_count": h["/myosin/center"].shape[1] == 320,
        "model22": int(m.get("am_model_version", -1)) == 22,
        "global_sites": int(m.get("am_sites_per_actin", -1)) == 182,
        "global_cap": int(m.get("am_max_heads_per_actin", -1)) == 24,
        "q_equilibrium": int(m.get("am_q_initialized_at_calcium_equilibrium", -1)) == 1,
        "ramp": int(m.get("interaction_ramp_steps", -1)) == 2000,
        "explicit_titin": int(m.get("titin_mode", -1)) == 2,
        "titin_slack": math.isclose(float(m.get("titin_rest_length", -1)), 0.457),
        "titin_capture": math.isclose(float(m.get("titin_capture_distance", -1)), 0.72),
        "titin_off": math.isclose(float(m.get("titin_k_off", -1)), 0.01),
        "titin_theta": math.isclose(float(m.get("titin_transition_theta", -1)), 1.0),
        "titin_actn2_gate": int(m.get("titin_require_live_actn2", -1)) == 1,
        "axial_inline": txt("am_head_association_axial_gate") ==
            "material_head_projection_must_lie_inside_finite_actin_segment",
        "contact_on": int(m.get("aa_double_overlap_contact_enabled", -1)) == 1,
    }
    failed = [name for name, ok in checks.items() if not ok]
    if failed:
        raise SystemExit("metadata validation failed: " + ", ".join(failed))
print("completed_metadata_validation=PASS")
PY

echo "finished_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
