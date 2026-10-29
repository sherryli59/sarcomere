#!/bin/bash
# Matched four-seed/two-arm 5 s production campaign.
# Arms differ only in kappa_aa (3 versus 0), matched at each seed.
# Launched only after the expanded minimal gate and 200 ms production gate passed.
# Requeue-safe resume handling is retained from the matched production script.
#SBATCH --partition=owners,normal,hns
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16GB
#SBATCH --time=2-00:00:00
#SBATCH --open-mode=append
#SBATCH --requeue
#SBATCH --array=0-3
#SBATCH --job-name=stallf
#SBATCH --output=/scratch/users/xsli/sarcomere/stallforce_20261005/logs/%x_%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/stallforce_20261005/logs/%x_%A_%a.err
set -euo pipefail

# ---------------------------------------------------------------------------
# PRODUCTION, new binary.  Lx=6 as mandated; A1600/M320 unchanged.
#
# The physical argument list matches submit_lowS_v3_20260925 except nsteps
# (500000, i.e. 5 s), plus the explicit steric radius below.
# --aa_zdisc_tip_contact=true and --aa_occlusion_rupture_penetration=0.010 are
# unchanged; --aa_occlusion_capture_radius=0.080 separates sterics from AM reach.
#
#   1. DetachmentReason 3 REDEFINED.  Was decide_am_overlap_occlusion(), firing
#      when an actin had no crosslinkable barbed contour left -- MEASURED to fire
#      86 nm past the thick-filament end with a 12 nm surface gap, i.e. with
#      nothing touching, and accounting for 23% of all ACTN2 detachments (100%
#      once the bond was made permanent).  Now
#      decide_barbed_tip_capsule_penetration(), firing only when a barbed tip
#      has actually penetrated the capture capsule by
#      aa_occlusion_rupture_penetration = 0.010 um.
#   2. compute_barbed_tip_myosin_endstop(): capsule-normal barbed-tip repulsion
#      from true point-to-finite-segment distance, with a pre-contact onset,
#      stiffness and force cap borrowed from aa_zdisc_tip_*.  Momentum
#      conserving and independent of any live ACTN2 link.
#   3. NEW flags aa_overlap_occlusion_veto (default true) and
#      aa_occlusion_capture_radius, explicitly set to 0.080 here.
#      AM capture remains 0.130 and is not changed by this steric radius.
#
# So this campaign pairs seed-for-seed against lowS_20260924 with exactly one
# variable changed.  Because the arguments cannot distinguish the two binaries,
# the output TAG carries "_endstop" and source_snapshot/ holds BINARY_SHA256.txt,
# GIT_STATE.txt and WORKING_TREE.patch.  Nothing else marks the generation.
#
# COST: eight 5 s production trajectories at 16 CPUs.
# --requeue is kept: the resume guard below
# resumes from committed_frames on a preemption restart instead of clobbering.
# save_every=500 yields 1001 frames per trajectory.
# ---------------------------------------------------------------------------
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

# ---------------------------------------------------------------------------
# MORE BINDING, NOT MORE PULLING FORCE -- am_stall_force 4.5 -> 1.5 at hph=27.
#
# Arm C of the head-geometry ladder (hph=27) tripled the crowns AND tripled the
# bound heads, and it failed: nematic S 0.738 -> 0.508, ACTN2 links per actin
# 0.939 -> 0.629 (-33%). The leading mechanism is MECHANICAL RUPTURE of the
# ACTN2 crosslinks. Reason 1 in local_link_events_v3 is named
# "geometric_contour" but catch_bond_kinetics.h:102 defines it purely as
# `connector_length > rupture_length`, i.e. overextension past
# aa_link_rupture_length = 0.12 um. It is 55% of detachments in the baseline and
# 65% in arm C, while stochastic death FALLS 39% -> 27%. Connector means barely
# move (83.2 -> 83.7 nm) but the near-rupture TAIL grows ~50% (fraction above
# 0.110 um: 0.66% -> 1.03%), and rupture is a threshold process.
#
# THIS ARM keeps every geometric gain of hph=27 and removes the extra pull.
# Per-head force is STALL-LIMITED here: am_motor_k 1000 x am_power_stroke 0.005
# = 5 pN, capped by am_stall_force = 4.5. Setting it to 1.5 cuts per-head force
# ~3x while the head COUNT stays tripled, so total drive is roughly preserved.
# One knob differs from arm C; seeds match A/B/C seed-for-seed.
#
# PREDICTION (sharp, so it can fail): ACTN2 links per actin recover toward
# 0.94 and S toward 0.74 WHILE multi-crown stays near 29%. If links recover but
# S does not, the rupture account is incomplete. If bound heads collapse, the
# arm is confounded -- lowering stall force also weakens head RETENTION, so
# check bound-head count before reading anything else.
#
# NOTE the principled fix is not this knob: _am_per_head_tether_stiffness
# (declared sarcomere.h:1345, defined sarcomere.cpp:3485) takes the bound-head
# count and has ZERO call sites -- per-head force is never renormalised by
# occupancy. This arm tests the idea; wiring that function up is the real fix.
# ---------------------------------------------------------------------------
root=/scratch/users/xsli/sarcomere/stallforce_20261005
exe=/scratch/users/xsli/sarcomere/barbed_capsule_normal_20260930/bin/sarcomere_barbed_capsule_normal

# Lx=6 is mandatory and is unchanged here.
arm=intermediate
# sigma is the Cartesian perturbation in set_direction_axis_noise, NOT an
# angle.  Calibrated against the initializer: 0.40 -> S0 ~ 0.62,
# 0.65 -> S0 ~ 0.35 (existing runs: 0.125 -> 0.95, 0.30 -> 0.76).
kappas=(3 3 3 3)
seeds=(20260805 20260815 20260825 20260835)
labels=(stallf15 stallf15 stallf15 stallf15)
hphs=(27 27 27 27)
spokess=(3 3 3 3)
sigmas=(0.40 0.40 0.40 0.40)
task=${SLURM_ARRAY_TASK_ID}
sigma=${sigmas[${task}]}
kappa_aa=${kappas[${task}]}
seed=${seeds[${task}]}
label=${labels[${task}]}
hph=${hphs[${task}]}
spokes=${spokess[${task}]}
armtag=sig${sigma}_kaa${kappa_aa}

nsteps=500000
save_every=500
tag=A1600_M320_Lx6_${arm}_${armtag}_seed${seed}_5000ms_${label}
output=${root}/data/${tag}.h5

test -x "${exe}"
sha256sum "${exe}" || true
echo "task ${task}: arm=${label} hph=${hph} spokes=${spokes} kappa_aa=${kappa_aa} seed=${seed} -> ${tag}"

# Fresh run, not a resume.  Refuse to clobber an existing trajectory.
committed=0
if [[ -e "${output}" ]]; then
    committed=$(python - "${output}" <<'PYFRAMES'
import sys
try:
    import h5py
    with h5py.File(sys.argv[1], "r", locking=False) as h:
        print(min(int(h["meta/committed_frames"][0]), h["actin/center"].shape[0]))
except Exception:
    print(0)
PYFRAMES
)
fi
resume=()
if [[ "${committed}" -gt 0 ]]; then
    if [[ "${SLURM_RESTART_COUNT:-0}" -eq 0 && "${PROD_ALLOW_RESUME:-0}" -ne 1 ]]; then
        echo "refusing to overwrite ${output} (${committed} frames); set" >&2
        echo "PROD_ALLOW_RESUME=1 to continue it deliberately" >&2; exit 1
    fi
    echo "resuming ${output} from ${committed} committed frames"
    resume=(--resume=true)
elif [[ -e "${output}" ]]; then
    rm -f "${output}"
fi

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
    --am_spatial_heads=true --am_heads_per_half=${hph}
    --am_spokes_per_crown=${spokes} --am_crossover_pull_through=true
    --am_max_heads_per_actin_half=24 --am_sites_per_actin=182
    --am_max_heads_per_actin=24 --am_k_on=100 --am_weak_binding_floor=0.1
    --am_k_weak_off=200 --am_k_powerstroke=200 --am_k_strong_off=160
    --am_force_scale=1 --am_motor_k=1000 --am_power_stroke=0.005
    --am_stall_force=1.5 --am_unloaded_velocity=5 --am_weak_stiffness_ratio=0.1
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
    --kappa_aa="${kappa_aa}" --aa_perpendicular_k=50 --aa_status1_radial_k=75
    --aa_status1_site_diffusion=0.010 --aa_arm_state_model=false
    --aa_actn2_contact_radius=0.005 --aa_myosin_contact_epsilon=0.00428
    --aa_myosin_contact_force_cap=40 --aa_myosin_contact_shell_width=0.015
    --aa_myosin_contact_shell_energy=0.0214 --aa_occlusion_rupture_penetration=0.010
    --aa_occlusion_capture_radius=0.080
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

"${cmd[@]}"
echo "task ${task} done"
