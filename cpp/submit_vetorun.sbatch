#!/bin/bash
#SBATCH --partition=owners,normal,hns
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4GB
#SBATCH --time=04:00:00
#SBATCH --open-mode=append
#SBATCH --array=0-7
#SBATCH --job-name=vetorun
#SBATCH --output=/scratch/users/xsli/sarcomere/vetorun_20260928/logs/%x_%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/vetorun_20260928/logs/%x_%A_%a.err
set -euo pipefail
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-2}

# ---------------------------------------------------------------------------
# Lx=4 EXCEPTION, FLAGGED.  Lx=6 is mandatory for any run whose result touches
# lattice period or registration.  This is the minimal-sarcomere rig -- 4
# actins, 2 myosins, one Z-disc pair per side -- built solely to test the
# ACTN2 am_overlap_occlusion veto radius.  These runs are INADMISSIBLE as
# evidence about lattice period, registration, or nematic order at scale.
#
# ARMS.  Exactly one knob differs between them:
#   veto130  --aa_occlusion_capture_radius unset  -> 0.130 um (myosin HEAD reach)
#            = current behaviour, bit-identical to HEAD (gate 1 of job 45791760)
#   veto080  --aa_occlusion_capture_radius=0.080  -> steric footprint only
#            (myosin_radius 0.05 + actin_radius 0.015 + am_surface_gap 0.015)
#
# PREDICTION being tested (from 80 logged occlusion events in sarc4_20260925):
#   the veto currently fires at d_tip = +0.078 um, OUTSIDE the 0.050 um
#   aa_zdisc_tip_onset, so the end-stop is engaged at only 15% of deletions.
#   Shrinking the capsule removes a median 0.0747 um of axial overhang, which
#   should move the veto to d_tip ~ +0.004 um, i.e. 46 nm INSIDE the onset.
# ---------------------------------------------------------------------------

root=/scratch/users/xsli/sarcomere/vetorun_20260928
exe=${root}/bin/sarcomere_veto
monitor=${root}/source_snapshot/monitor_live_production_hdf5.py

arms=(veto130 veto080)
# Same 4 seeds as sarc4_20260925, so the veto130 control arm is directly
# comparable to the runs the prediction was measured on.
seeds=(20260805 20260815 20260825 20260835)
task=${SLURM_ARRAY_TASK_ID}
armtag=${arms[$(( task / 4 ))]}
seed=${seeds[$(( task % 4 ))]}

vetoflag=()
case "${armtag}" in
    veto130) ;;                                              # unset: 0.130 head reach
    veto080) vetoflag=(--aa_occlusion_capture_radius=0.080) ;;
    *) echo "unknown arm ${armtag}" >&2; exit 1 ;;
esac

nsteps=500000
save_every=100
tag=vetorun_${armtag}_seed${seed}
output=${root}/data/${tag}.h5
live_log=${root}/analysis/${tag}.live.jsonl

test -x "${exe}"
echo "task ${task}: arm=${armtag} seed=${seed} flag='${vetoflag[*]-none}' -> ${tag}"

# Fresh run only.  aa_occlusion_capture_radius is resume-guarded, and a resume
# would silently keep the old value.
if [[ -e "${output}" ]]; then
    echo "refusing to overwrite ${output}" >&2; exit 1
fi

exec 9>"${output}.lock"
flock -n 9

cmd=(srun --kill-on-bad-exit=1 "${exe}"
  --n_actins=4 --n_myosins=2 --Lx=4 --Ly=1 --Lz=1
  --dimension=3 --pbc_mask=111 --actin_length=1 --myosin_length=1.5
  --actin_radius=0.015 --myosin_radius=0.05
  --initial_structure=minimal_catch_balanced
  --actin_direction_polar_fraction=0.5 --myosin_direction_polar_fraction=0.5
  --pre_phase_steps=500 --interaction_ramp_steps=2000
  --fixed_lattice_preequil_steps=0 --physical_mobility=true
  --temperature_kelvin=310 --thermal_noise_scale=1 --solvent_viscosity=0.001
  --actin_hydrodynamic_diameter=0.010 --myosin_hydrodynamic_diameter=0.035
  --actin_crowding_parallel=10 --actin_crowding_perpendicular=1000
  --actin_crowding_tumble=160 --actin_crowding_roll=160
  --myosin_crowding_parallel=20 --myosin_crowding_perpendicular=1000
  --myosin_crowding_tumble=64 --myosin_crowding_roll=64
  --max_actin_displacement=0.005 --max_myosin_displacement=0.005
  --max_actin_rotation=0.005 --max_myosin_rotation=0.005
  --sterics=true --actin_steric_k=2000 --myosin_steric_k=1000
  --am_steric_k=1000 --myosin_exclusion_distance=0
  --am_cutoff=0.05 --am_optimal=0.055 --am_surface_gap=0.015
  --k_am=0 --kappa_am=0 --dynamic_am=true --am_offcenter_torque=true
  --am_head_azimuth=false --torque_breakdown=true
  --am_spatial_heads=true --am_heads_per_half=9 --am_crossover_pull_through=true
  --am_max_heads_per_actin_half=24 --am_sites_per_actin=182
  --am_max_heads_per_actin=24 --am_k_on=100 --am_weak_binding_floor=0.1
  --am_k_weak_off=200 --am_k_powerstroke=200 --am_k_strong_off=160
  --am_force_scale=1 --am_motor_k=1000 --am_power_stroke=0.005
  --am_stall_force=4.5 --am_unloaded_velocity=5 --am_weak_stiffness_ratio=0.1
  --am_tether_k=3000 --am_tether_slack=0.002 --am_head_transverse_k=100
  --am_transverse_strain_detachment_theta=0 --am_alignment_k=0
  --am_alignment_slack_deg=2 --am_tether_saturation_heads=6
  --am_binding_azimuth_selection=true --am_binding_azimuth_half_width_deg=30
  --am_rupture_distance=0.1 --am_strain_dependent_detachment=true
  --am_strain_detachment_theta=0.1 --am_strain_detachment_max_exponent=20
  --am_strong_off_floor_fraction=1 --am_cycle_period=0 --am_activation_floor=0
  --am_activation_active_fraction=0.5 --am_thin_activation=true
  --am_thin_feedback=true --am_q_calcium_open_rate=20 --am_q_calcium_half=0.5
  --am_q_calcium_hill=3 --am_q_strong_feedback_rate=100
  --am_q_strong_half=0.02 --am_q_strong_hill=2 --am_q_close_rate=20
  --am_initialize_q_calcium_equilibrium=true --am_thick_feedback=true
  --am_thick_recruit_rate=15 --am_thick_fold_rate=5
  --am_thick_tension_scale=3 --am_thick_initial_on_fraction=0.75
  --am_status2_memory=0 --am_bidirectional_crossbridge=false
  --am_overlap_target_ratio=0 --actin_turnover=false
  --actin_turnover_rate=0 --actin_turnover_start_time=0.2
  --actin_turnover_max_placement_attempts=5000
  --actin_turnover_clearance_ratio=0.80
  --aa_barbed_lateral_range=0.080 --aa_barbed_lateral_k=2000
  --aa_barbed_lateral_force_max=40 --aa_barbed_lateral_axial_window=0.250
  --aa_local_multilink=true --aa_spacing_ratio=0.9875 --aa_preferred_centerline=0.079
  --aa_capture_margin=0.0396 --aa_antiparallel_overlap=0.10 --aa_barbed_zone=0.20
  --aa_axial_register_ratio=3 --aa_axial_register_force_max=36
  --aa_register_all_statuses=false --aa_status2_registry_on=0.040
  --aa_status2_registry_off=0.040 --aa_status2_support_window=0.010
  --aa_register_saturation_error=0.05 --aa_prebundle_links_per_pair=3
  --aa_site_exclusion=0.0192 --aa_overlap_capture=0.10 --k_aa=1000
  --kappa_aa=3 --aa_perpendicular_k=50 --aa_status1_radial_k=75
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
  --titin_require_live_actn2=true --titin_k_off_unsupported=10
  --titin_rebind_requires_actn2=true --k_mm=0 --k_bundle_max=0
  --bundle_ramp_steps=0 --dt=1e-5 --nsteps="${nsteps}" --save_every="${save_every}" --seed="${seed}" --filename="${output}"
  ${vetoflag[@]+"${vetoflag[@]}"})

"${cmd[@]}" &
simulation_pid=$!
python "${monitor}" "${output}" "${live_log}" --pid "${simulation_pid}" || true
wait "${simulation_pid}"
echo "task ${task} done"
