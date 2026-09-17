#!/bin/bash
#SBATCH --partition=rotskoff,owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64GB
# owners caps wall time at 2 days; ~15 h expected for 500k steps at A1600.
#SBATCH --time=2-00:00:00
#SBATCH --open-mode=append
#SBATCH --requeue
#SBATCH --array=0-5
#SBATCH --job-name=actn2_rest74
#SBATCH --output=/scratch/users/xsli/sarcomere/actn2_rest74_lateral_20260831/logs/%x_%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/actn2_rest74_lateral_20260831/logs/%x_%A_%a.err

set -euo pipefail
# Lx=6 is mandatory.  The model selects a 2.00 um sarcomere period (dominant
# harmonic k=3 in Lx=6).  Lx=3 would need k=1.5, which PBC cannot host, so the
# 3x2x2 box shows a flat band spectrum and no axial order.  Density is held at
# 66.7 actin/um^3 and 13.3 myosin/um^3, so Lx=6 implies A1600/M320.
export HDF5_USE_FILE_LOCKING=FALSE
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-32}

repo=/home/groups/rotskoff/sherry/sarcomere/sarcomere
root=${ASSEMBLY_PILOT_ROOT:-/scratch/users/xsli/sarcomere/actn2_rest74_lateral_20260831}
exe=${root}/bin/sarcomere_assembly_controls_20260828
monitor=${root}/source_snapshot/monitor_live_production_hdf5.py
mkdir -p "${root}"/{data,params,analysis,logs}

# ACTN2 rest length raised 0.055 -> 0.0736 um = 0.92 * am_preferred_centerline
# (lattice_spacing.h physiological_aa_am_ratio).  aa_capture_margin stays
# 0.045, so aa_cutoff grows 0.100 -> 0.1186 as an intended consequence.
# Rationale: the lateral repulsion stalls at ~46 nm because ACTN2 pulls back
# to its 55 nm rest length with ~40 pN, exactly matching the lateral cap.
arms=(rest74_off rest74_on)
seeds=(20260805 20260815 20260825)
task=${SLURM_ARRAY_TASK_ID}
arm=${arms[$((task / 3))]}
seed=${seeds[$((task % 3))]}

# Full-ON motor, turnover off: isolate the lateral term against the arm whose
# collapse we characterised (turn0 seed 20260805).
turnover=false
turnover_rate=0
thick_feedback=false
recruit=2
fold=18
initial_on=0.1
clearance_ratio=0.80

# Sized against the measured collapse: the capsule steric delivers only ~8 pN
# at the 26 nm the ends currently settle to, so at 26 nm this term must exceed
# that to open the spacing.  range 0.075 sits between the 55 nm ACTN2 rest
# length and the 81 nm that would tile a Z-disc face; k=300 gives 14.7 pN at
# 26 nm.  The cap matches the other contact caps in the model (20-40 pN).
lateral_range=0
lateral_k=0
lateral_fmax=0
case "${arm}" in
    rest74_off) ;;
    rest74_on) lateral_range=0.080; lateral_k=2000; lateral_fmax=40 ;;
    *) exit 2 ;;
esac

dt=1e-5
# Production duration per the plan: 5 s at dt=1e-5, 1000 saved frames,
# matching the A2400 Lx=6 reference cadence (~12 GB per trajectory).
nsteps=${ASSEMBLY_PILOT_NSTEPS:-100000}
save_every=${ASSEMBLY_PILOT_SAVE_EVERY:-500}
turnover_start=${ASSEMBLY_TURNOVER_START:-0.2}
duration_ms=$((nsteps / 100))
tag=A1600_M320_Lx6_${arm}_seed${seed}_${duration_ms}ms
output=${root}/data/${tag}.h5
params=${root}/params/${tag}.params
live_log=${root}/analysis/${tag}.live.jsonl

test -x "${exe}"
test -r "${monitor}"
# A task preempted before its first checkpoint leaves a file with zero
# committed frames.  Resuming from that aborts the binary with
#   "No frames available in actin center dataset for resume."
# so decide on committed frames, not on mere file existence.
committed=0
if [[ -e "${output}" ]]; then
    committed=$(python - "${output}" <<'PYFRAMES'
import sys
try:
    import h5py
    with h5py.File(sys.argv[1], "r", locking=False) as handle:
        n = int(handle["meta/committed_frames"][0]) if "meta/committed_frames" in handle else 0
        m = handle["actin/center"].shape[0] if "actin/center" in handle else 0
        print(min(n, m))
except Exception:
    print(0)
PYFRAMES
)
fi
if [[ -e "${output}" && "${committed}" -gt 0 && "${SLURM_RESTART_COUNT:-0}" -eq 0 && "${ASSEMBLY_PILOT_RESUME_EXISTING:-false}" != true ]]; then
    echo "refusing to overwrite ${output}" >&2
    exit 1
fi
resume=()
if [[ "${committed}" -gt 0 ]]; then
    resume=(--resume=true)
elif [[ -e "${output}" ]]; then
    echo "discarding zero-frame ${output} (preempted before first checkpoint)" >&2
    rm -f "${output}"
fi
exec 9>"${output}.lock"
flock -n 9

printf '%s\n' \
    "arm=${arm}" "seed=${seed}" "turnover_rate=${turnover_rate}" \
    "thick_feedback=${thick_feedback}" "recruit=${recruit}" "fold=${fold}" \
    "initial_on=${initial_on}" "tension_scale=3" \
    "clearance_ratio=${clearance_ratio}" \
    "lateral_range=${lateral_range}" "lateral_k=${lateral_k}" "lateral_fmax=${lateral_fmax}" \
    "turnover_start=${turnover_start}" >"${params}"

cmd=(srun --kill-on-bad-exit=1 "${exe}"
    --n_actins=1600 --n_myosins=320 --Lx=6 --Ly=2 --Lz=2
    --dimension=3 --pbc_mask=111 --actin_length=1 --myosin_length=1.5
    --actin_radius=0.015 --myosin_radius=0.05
    --initial_structure=myosin_x_noise --actin_direction_noise=-1
    --myosin_direction_noise=1 --myosin_direction_tilt_deg=0
    --pre_phase_steps=500 --interaction_ramp_steps=2000
    --fixed_lattice_preequil_steps=0 --physical_mobility=true
    --temperature_kelvin=310 --thermal_noise_scale=1 --solvent_viscosity=0.001
    --actin_hydrodynamic_diameter=0.010 --myosin_hydrodynamic_diameter=0.035
    --actin_crowding_parallel=10 --actin_crowding_perpendicular=10
    --actin_crowding_tumble=160 --actin_crowding_roll=160
    --myosin_crowding_parallel=20 --myosin_crowding_perpendicular=20
    --myosin_crowding_tumble=64 --myosin_crowding_roll=64
    --max_actin_displacement=0.005 --max_myosin_displacement=0.005
    --max_actin_rotation=0.005 --max_myosin_rotation=0.005
    --sterics=true --actin_steric_k=2000 --myosin_steric_k=1000
    --am_steric_k=1000 --myosin_exclusion_distance=0
    --am_cutoff=0.05 --am_optimal=0.055 --am_surface_gap=0.015
    --k_am=0 --kappa_am=0 --dynamic_am=true --am_offcenter_torque=true
    --am_spatial_heads=true --am_heads_per_half=9 --am_crossover_pull_through=true
    --am_max_heads_per_actin_half=24 --am_sites_per_actin=182
    --am_max_heads_per_actin=24 --am_k_on=100 --am_weak_binding_floor=0.1
    --am_k_weak_off=200 --am_k_powerstroke=200 --am_k_strong_off=160
    --am_force_scale=1 --am_motor_k=1000 --am_power_stroke=0.005
    --am_stall_force=4.5 --am_unloaded_velocity=5 --am_weak_stiffness_ratio=0.1
    --am_tether_k=3000 --am_tether_slack=0.002 --am_head_transverse_k=100
    --am_transverse_strain_detachment_theta=0 --am_alignment_k=120
    --am_alignment_slack_deg=2 --am_tether_saturation_heads=6
    --am_binding_azimuth_selection=true --am_binding_azimuth_half_width_deg=30
    --am_rupture_distance=0.1 --am_strain_dependent_detachment=true
    --am_strain_detachment_theta=0.1 --am_strain_detachment_max_exponent=20
    --am_strong_off_floor_fraction=1 --am_cycle_period=0 --am_activation_floor=0
    --am_activation_active_fraction=0.5 --am_thin_activation=true
    --am_thin_feedback=true --am_q_calcium_open_rate=20 --am_q_calcium_half=0.5
    --am_q_calcium_hill=3 --am_q_strong_feedback_rate=100
    --am_q_strong_half=0.02 --am_q_strong_hill=2 --am_q_close_rate=20
    --am_initialize_q_calcium_equilibrium=true --am_thick_feedback="${thick_feedback}"
    --am_thick_recruit_rate="${recruit}" --am_thick_fold_rate="${fold}"
    --am_thick_tension_scale=3 --am_thick_initial_on_fraction="${initial_on}"
    --am_status2_memory=0 --am_bidirectional_crossbridge=false
    --am_overlap_target_ratio=0 --actin_turnover="${turnover}"
    --actin_turnover_rate="${turnover_rate}" --actin_turnover_start_time="${turnover_start}"
    --actin_turnover_max_placement_attempts=5000
    --actin_turnover_clearance_ratio="${clearance_ratio}"
    --aa_barbed_lateral_range="${lateral_range}"
    --aa_barbed_lateral_k="${lateral_k}"
    --aa_barbed_lateral_force_max="${lateral_fmax}"
    --aa_barbed_lateral_axial_window=0.250
    --aa_local_multilink=true --aa_spacing_ratio=1 --aa_preferred_centerline=0.0736
    --aa_capture_margin=0.045 --aa_antiparallel_overlap=0.10 --aa_barbed_zone=0.20
    --aa_axial_register_ratio=3 --aa_axial_register_force_max=36
    --aa_register_all_statuses=false --aa_status2_registry_on=0.040
    --aa_status2_registry_off=0.040 --aa_status2_support_window=0.010
    --aa_register_saturation_error=0.05 --aa_prebundle_links_per_pair=3
    --aa_site_exclusion=0.0192 --aa_overlap_capture=0.10 --k_aa=1000
    --kappa_aa=100 --aa_perpendicular_k=50 --aa_status1_radial_k=75
    --aa_status1_site_diffusion=0.010 --aa_arm_state_model=false
    --aa_actn2_contact_radius=0.005 --aa_myosin_contact_epsilon=0.00428
    --aa_myosin_contact_force_cap=40 --aa_myosin_contact_shell_width=0.015
    --aa_myosin_contact_shell_energy=0.0214 --aa_occlusion_rupture_penetration=0.010
    --aa_zdisc_tip_contact=true --aa_zdisc_tip_onset=0.050
    --aa_zdisc_tip_stiffness=1000 --aa_zdisc_tip_force_max=20
    --aa_zdisc_tip_transverse_flat=0.100 --aa_zdisc_tip_transverse_cutoff=0.130
    --aa_status2_locked_until_detach=true --aa_stiffness_reference_links=3
    --aa_total_link_pool=0 --aa_status_mode=productive_register
    --aa_status2_force_on=0.25 --aa_status2_force_off=0.10 --k_on=2500
    --k_off=1 --base_lifetime=0.1 --directional_base_lifetime=0.05
    --lifetime_coeff=0.4 --stabilize_status1=false
    --aa_energy_gated_association=false --aa_transition_state_theta=0.1
    --aa_encounter_k=10 --aa_encounter_kappa=0 --aa_passive_detailed_balance=false
    --tau_rec=0 --max_aa_bonds_per_actin=4 --aa_neighbor_identity=true
    --aa_neighbor_slot_memory=0.02 --aa_status2_slot_memory=0.20
    --aa_association_max_direction_dot=-0.5 --aa_link_rupture_length=0.12
    --aa_double_overlap_contact=true --aa_double_overlap_contact_diameter=0.012
    --aa_double_overlap_contact_epsilon=0.00428 --aa_double_overlap_contact_force_cap=0
    --directional=true --max_myosin_bonds=24 --titin_mode=explicit_dynamic
    --titin_k=100 --titin_rest_length=0.665 --titin_links_per_half=3
    --titin_capture_distance=0.72 --titin_k_on=5000 --titin_k_off=0.01
    --titin_strain_gated_kinetics=true --titin_transition_theta=1
    --titin_rebind_memory=0.020 --titin_rebind_rate_multiplier=100
    --titin_require_live_actn2=true --titin_k_off_unsupported=0
    --titin_rebind_requires_actn2=false --k_mm=0 --k_bundle_max=0
    --bundle_ramp_steps=0 --dt="${dt}" --nsteps="${nsteps}"
    --save_every="${save_every}" --seed="${seed}" --filename="${output}"
    ${resume[@]+"${resume[@]}"})

"${cmd[@]}" &
simulation_pid=$!
python "${monitor}" "${output}" "${live_log}" --pid "${simulation_pid}" \
    --initial-seconds 60 --interval-seconds 300 \
    --actin-force-cap-fraction-max 0.75 --myosin-force-cap-fraction-max 0.50 &
monitor_pid=$!
set +e
wait "${simulation_pid}"
status=$?
wait "${monitor_pid}"
set -e
exit "${status}"

