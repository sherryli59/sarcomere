#!/bin/bash
# PHASE 1 -- torque attribution at am_alignment_k = 0, pre-aligned.
#
# Phase 0 established that order here is ACTIVELY MAINTAINED: S holds at
# 0.887-0.899 over 3.03 P2 rotational relaxation times against a free-diffusion
# prediction of 0.036-0.107, with zero turnover so it can only be rotation.
# The equilibrium angular width gives kappa_eff = 23.5-25.0 kT/rad^2.  This run
# asks WHICH of the 12 instrumented torque sources supplies it.
#
# These are RESUMED COPIES.  The originals in prealigned_k0_5s_20260922/ are
# left untouched: resuming appends frames and adds ~37 per-actin datasets, and
# those trajectories are campaign results.
#
# normal,hns -- neither preempts.  An owners requeue would restart from step 0
# and overwrite the incrementally written trajectory.
#SBATCH --partition=normal,hns
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=16GB
#SBATCH --time=06:00:00
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --array=0-2
#SBATCH --job-name=torqattr
#SBATCH --output=/scratch/users/xsli/sarcomere/torqueattr_20260923/logs/%x_%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/torqueattr_20260923/logs/%x_%A_%a.err
set -euo pipefail

root=/scratch/users/xsli/sarcomere/torqueattr_20260923
# The SAME staged binary the resumed trajectories were produced with.  The
# repo's cpp/build/sarcomere is a Sep 2 build and would be a different model.
exe=${root}/bin/sarcomere_prealign_20260922
monitor=${root}/source_snapshot/monitor_live_production_hdf5.py

# Lx=6 is mandatory (standing instruction) and is unchanged here: these are
# resumes of Lx=6 trajectories.
arm=intermediate
sigmas=(0.30 0.30 0.125)
seeds=(20260825 20260845 20260855)
task=${SLURM_ARRAY_TASK_ID}
sigma=${sigmas[${task}]}
seed=${seeds[${task}]}

save_every=500
tag=A1600_M320_Lx6_${arm}_sig${sigma}_seed${seed}_5000ms
output=${root}/data/${tag}.h5
live_log=${root}/analysis/${tag}.live.jsonl

test -x "${exe}"
test -r "${output}"

# Resume must ADVANCE.  nsteps is absolute, so asking for the step count the
# file already reached runs zero steps and passes vacuously.  Derive it from
# the committed frames actually present and add 30 frames.
EXTRA_FRAMES=30
committed=$(python - "${output}" <<'PYFRAMES'
import sys
import h5py
with h5py.File(sys.argv[1], "r", locking=False) as h:
    n = int(h["meta/committed_frames"][0])
    m = h["actin/center"].shape[0]
    print(min(n, m))
PYFRAMES
)
if [[ "${committed}" -lt 1 ]]; then
    echo "FATAL: ${output} has no committed frames to resume from" >&2
    exit 2
fi
nsteps=$(( (committed + EXTRA_FRAMES) * save_every ))
echo "task ${task}: sigma=${sigma} seed=${seed} committed=${committed} -> nsteps=${nsteps}"
resume=(--resume=true)

exec 9>"${output}.lock"
flock -n 9

cmd=(srun --kill-on-bad-exit=1 "${exe}"
    --n_actins=1600 --n_myosins=320 --Lx=6 --Ly=2 --Lz=2
    --dimension=3 --pbc_mask=111 --actin_length=1 --myosin_length=1.5
    --actin_radius=0.015 --myosin_radius=0.05
    --initial_structure=myosin_x_noise --actin_direction_noise="${sigma}"
    --myosin_direction_noise="${sigma}" --myosin_direction_tilt_deg=0
    --actin_direction_polar_fraction=0.5
    --myosin_direction_polar_fraction=0.5
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
    --am_head_azimuth=false
    --torque_breakdown=true
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
    --aa_barbed_lateral_range=0.080
    --aa_barbed_lateral_k=2000
    --aa_barbed_lateral_force_max=40
    --aa_barbed_lateral_axial_window=0.250
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
    --bundle_ramp_steps=0 --dt=1e-5 --nsteps="${nsteps}"
    --save_every="${save_every}" --seed="${seed}" --filename="${output}"
    ${resume[@]+"${resume[@]}"})

"${cmd[@]}" &
simulation_pid=$!
python "${monitor}" "${output}" "${live_log}" --pid "${simulation_pid}" || true
wait "${simulation_pid}"
echo "task ${task} done"
