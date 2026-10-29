# Shared production argv.  Copied verbatim from
# azgate_20261015/submit_hdbasekam.sbatch (the hph27 5 s production campaign),
# with ONLY the burial block parameterized.  Anything else differing from that
# file is a bug.
#
# Callers set, before sourcing:  EXE OUT SEED HPH NSTEPS SAVE_EVERY RESUME_ARGS
#
# THE BURIAL BLOCK NOW DEFAULTS TO RULE-ONLY (standing instruction 2026-10-23).
# A caller that sets none of KILL/RATE/SLAB/CUTOFF gets:
#   kill line 50 nm, hazard rate 1/s, slab OFF.
# MEASURED across a 2 head-count x 2 rule x 3 alignment-level factorial at 200 ms:
# rule-only adds +1.2 to +2.1 ACTN2 links per actin in ALL SIX cells -- with the
# AM aligner off, head-ramped, or normalised -- and cuts burial deaths ~27000 to
# ~1000.  The slab geometry is NOT included: it adds nothing at hph27 (t=1.4),
# costs nematic order, and runs ~3x slower.
# TITIN AND MYOMESIN OFF (standing instruction 2026-10-27): titin_k=0 with the default
# legacy_crosslinked topology skips every titin code path (all guarded by titin_k > 0);
# k_mm=0 disables myomesin. Neither term is what its name suggests in this model.
# Z-DISC STOP GATED ON ACTN2 (standing instruction 2026-10-25):
#   --aa_zdisc_tip_contact_requires_link=true.  The myosin capsule rejects the
#   Z-disc (an actin holding a live ACTN2 link), not a bare barbed tip.  Before
#   this every production run shoved uncrosslinked actins at the 20 pN cap.
build_prod_cmd() {
  KILL=${KILL:-0.050}
  RATE=${RATE:-1}
  SLAB=${SLAB:-false}
  CUTOFF=${CUTOFF:-0}
  PROD_CMD=(srun --kill-on-bad-exit=1 "${EXE}"
    --n_actins=1600 --n_myosins=320 --Lx=6 --Ly=2 --Lz=2
    --dimension=3 --pbc_mask=111 --actin_length=1 --myosin_length=1.5
    --actin_radius=0.015 --myosin_radius=0.05
    --initial_structure=myosin_x_noise --actin_direction_noise=0.40
    --myosin_direction_noise=0.40 --myosin_direction_tilt_deg=0
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
    --am_spatial_heads=true --am_heads_per_half="${HPH}" --am_spokes_per_crown=3
    --am_crossover_pull_through=true
    --am_max_heads_per_actin_half=24 --am_sites_per_actin=182
    --am_max_heads_per_actin=24 --am_k_on=100 --am_weak_binding_floor=0.1
    --am_k_weak_off=200 --am_k_powerstroke=200 --am_k_strong_off=160
    --am_force_scale=1 --am_motor_k=1000 --am_power_stroke=0.005
    --am_stall_force=4.5 --am_unloaded_velocity=5 --am_weak_stiffness_ratio=0.1
    --am_tether_k=3000 --am_tether_slack=0.002 --am_head_transverse_k=100
    --am_transverse_strain_detachment_theta=0 --am_alignment_k=3
    # STANDING INSTRUCTION 2026-10-22: sliding force must NOT scale with
    # head count.  msat6 caps each head's axial stiffness and stall cap by
    # 1/min(bound,6) so a pair's TOTAL sliding force saturates.  MEASURED at
    # hph27/soff160/kam3, 2 seeds: bound heads unchanged (0.95x) but overlap
    # 0.87x, free contour 1.19x and ACTN2 links 1.42x.
    --am_motor_saturation_heads=6
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
    --kappa_aa=3 --aa_perpendicular_k=50 --aa_status1_radial_k=1000
    --aa_status1_site_diffusion=0.010 --aa_arm_state_model=false
    --aa_actn2_contact_radius=0.005 --aa_myosin_contact_epsilon=0.00428
    --aa_myosin_contact_force_cap=40 --aa_myosin_contact_shell_width=0.015
    --aa_myosin_contact_shell_energy=0.0214
    ### ---- THE BURIAL BLOCK: the only thing that varies between arms ----
    --aa_occlusion_rupture_penetration="${KILL}"
    --aa_occlusion_rupture_rate="${RATE}"
    --aa_barbed_tip_endstop_slab="${SLAB}"
    --aa_barbed_tip_endstop_transverse_cutoff="${CUTOFF}"
    ### -------------------------------------------------------------------
    --aa_occlusion_capture_radius=0.080
    --aa_zdisc_tip_contact=true --aa_zdisc_tip_contact_requires_link=true --aa_zdisc_tip_onset=0.050
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
    --directional=true --max_myosin_bonds=24 --titin_mode=legacy_crosslinked
    --titin_k=0 --k_mm=0 --k_bundle_max=0
    --bundle_ramp_steps=0 --dt=1e-5 --nsteps="${NSTEPS}"
    --save_every="${SAVE_EVERY}" --seed="${SEED}" --filename="${OUT}"
    ${RESUME_ARGS[@]+"${RESUME_ARGS[@]}"})
}
