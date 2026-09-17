// HDF5 metadata writer, split out of sarcomere.cpp on 2026-09-17.
//
// Sarcomere::new_file() is the single place every /meta attribute is
// emitted.  It lived 8600 lines into a 13,109-line file, which made the
// attribute set hard to find and hard to audit -- a survey that day found
// 84 of 259 CLI options had no attribute at all.  The body below is a
// byte-for-byte move: no logic was changed in the split.

#include "sarcomere.h"
#include "aa_refractory.h"
#include "aa_association_kinetics.h"
#include "catch_bond_kinetics.h"
#include "counter_rng.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <thread>
#include <tuple>
#include <unordered_set>
#include <gsl/gsl_randist.h>

void Sarcomere::new_file(){
    create_file(filename, actin, myosin, max_myosin_bonds);
    H5::H5File file(filename, H5F_ACC_RDWR);
    H5::Group meta = file.openGroup("/meta");
    H5::DataSpace scalar(H5S_SCALAR);
    const int mode_value = static_cast<int>(aa_status_mode);
    H5::Attribute mode_attr = meta.createAttribute(
        "aa_status_mode", H5::PredType::NATIVE_INT, scalar);
    mode_attr.write(H5::PredType::NATIVE_INT, &mode_value);
    H5::Attribute on_attr = meta.createAttribute(
        "aa_status2_force_on", H5::PredType::NATIVE_DOUBLE, scalar);
    on_attr.write(H5::PredType::NATIVE_DOUBLE, &aa_status2_force_on);
    H5::Attribute off_attr = meta.createAttribute(
        "aa_status2_force_off", H5::PredType::NATIVE_DOUBLE, scalar);
    off_attr.write(H5::PredType::NATIVE_DOUBLE, &aa_status2_force_off);
    const int aa_load_schema = 2;
    H5::Attribute load_schema_attr = meta.createAttribute(
        "aa_load_transition_schema", H5::PredType::NATIVE_INT, scalar);
    load_schema_attr.write(H5::PredType::NATIVE_INT, &aa_load_schema);
    const int am_model_version = dynamic_am ? 22 : 0;
    H5::Attribute am_model_attr = meta.createAttribute(
        "am_model_version", H5::PredType::NATIVE_INT, scalar);
    am_model_attr.write(H5::PredType::NATIVE_INT, &am_model_version);
    H5::Attribute stall_attr = meta.createAttribute(
        "am_stall_force", H5::PredType::NATIVE_DOUBLE, scalar);
    stall_attr.write(H5::PredType::NATIVE_DOUBLE, &am_stall_force);
    H5::Attribute motor_k_attr = meta.createAttribute(
        "am_motor_k", H5::PredType::NATIVE_DOUBLE, scalar);
    motor_k_attr.write(H5::PredType::NATIVE_DOUBLE, &am_motor_k);
    H5::Attribute stroke_attr = meta.createAttribute(
        "am_power_stroke", H5::PredType::NATIVE_DOUBLE, scalar);
    stroke_attr.write(H5::PredType::NATIVE_DOUBLE, &am_power_stroke);
    H5::Attribute heads_attr = meta.createAttribute(
        "am_heads_per_half", H5::PredType::NATIVE_INT, scalar);
    heads_attr.write(H5::PredType::NATIVE_INT, &am_heads_per_half);
    auto write_int_attr = [&](const char* name, int value) {
        H5::Attribute attr = meta.createAttribute(
            name, H5::PredType::NATIVE_INT, scalar);
        attr.write(H5::PredType::NATIVE_INT, &value);
    };
    auto write_double_attr = [&](const char* name, double value) {
        H5::Attribute attr = meta.createAttribute(
            name, H5::PredType::NATIVE_DOUBLE, scalar);
        attr.write(H5::PredType::NATIVE_DOUBLE, &value);
    };
    auto write_string_attr = [&](const char* name, const std::string& value) {
        H5::StrType type(H5::PredType::C_S1, value.size() + 1);
        type.setStrpad(H5T_STR_NULLTERM);
        H5::Attribute attr = meta.createAttribute(name, type, scalar);
        attr.write(type, value);
    };
    write_string_attr("cli_invocation", cli_invocation);
    write_string_attr("cli_defaults", cli_defaults);
    write_double_attr("box_x", box.size() > 0 ? box[0] : 0.0);
    write_double_attr("box_y", box.size() > 1 ? box[1] : 0.0);
    write_double_attr("box_z", box.size() > 2 ? box[2] : 0.0);
    write_int_attr("periodic_x", is_periodic[0] ? 1 : 0);
    write_int_attr("periodic_y", is_periodic[1] ? 1 : 0);
    write_int_attr("periodic_z", is_periodic[2] ? 1 : 0);
    write_int_attr("actin_turnover_enabled",
                   actin_turnover_enabled ? 1 : 0);
    write_double_attr("actin_turnover_rate", actin_turnover_rate);
    write_double_attr("actin_turnover_start_time",
                      actin_turnover_start_time);
    write_int_attr("actin_turnover_max_placement_attempts",
                   actin_turnover_max_placement_attempts);
    write_double_attr("actin_turnover_clearance_ratio",
                      actin_turnover_clearance_ratio);
    write_string_attr(
        "actin_turnover_semantics",
        "fixed_population_whole_filament_reseed_uniform_hazard_"
        "steric_rejection_preserve_x_polarity");
    write_string_attr(
        "actin_turnover_events_columns",
        "actin,step,age_steps,success,placement_attempts,released_am_heads,"
        "released_actn2,released_titin,old_center_x,old_center_y,"
        "old_center_z,new_center_x,new_center_y,new_center_z,old_direction_x,"
        "old_direction_y,old_direction_z,new_direction_x,new_direction_y,"
        "new_direction_z");
    if (!initial_geometry_source.empty()) {
        write_string_attr("initial_geometry_source", initial_geometry_source);
        write_int_attr("initial_geometry_source_frame_requested",
                       initial_geometry_source_frame_requested);
        write_int_attr("initial_geometry_source_frame_loaded",
                       initial_geometry_source_frame_loaded);
        write_string_attr("initial_geometry_import_scope",
                          "filament_centers_and_directions_only");
    }
    write_double_attr("am_tether_k", am_tether_k);
    write_double_attr("am_tether_occupancy_floor", am_tether_occupancy_floor);
    write_double_attr("am_tether_occupancy_floor_reach",
                      am_tether_occupancy_floor_reach);
    write_double_attr("probe_tilt_deg", probe_tilt_deg);
    write_double_attr("probe_confine_k", probe_confine_k);
    write_double_attr("am_tether_slack", am_tether_slack);
    write_double_attr("am_head_transverse_k", am_head_transverse_k);
    write_double_attr("am_alignment_k", am_alignment_k);
    write_double_attr("am_alignment_slack_rad", am_alignment_slack);
    write_int_attr("am_tether_saturation_heads", am_tether_saturation_heads);
    write_int_attr("am_max_heads_per_actin_half",
                   am_max_heads_per_actin_half);
    write_int_attr("am_sites_per_actin", am_sites_per_actin);
    write_int_attr("am_max_heads_per_actin", am_max_heads_per_actin);
    write_string_attr(
        "am_site_ownership_law",
        "global_material_actin_sites_shared_across_all_myosins_strong_then_low_strain_priority");
    write_string_attr(
        "am_double_overlap_law",
        "wrong_half_heads_may_bind_weakly_but_cannot_powerstroke");
    write_double_attr("am_rupture_distance", am_rupture_distance);
    write_int_attr("am_strain_dependent_detachment",
                   am_strain_dependent_detachment ? 1 : 0);
    write_double_attr("am_strain_detachment_theta",
                      am_strain_detachment_theta);
    write_double_attr("am_transverse_strain_detachment_theta",
                      am_transverse_strain_detachment_theta);
    write_double_attr("am_strain_detachment_max_exponent",
                      am_strain_detachment_max_exponent);
    write_double_attr("am_strong_off_floor_fraction",
                      am_strong_off_floor_fraction);
    write_double_attr("am_cycle_period", am_cycle_period);
    write_double_attr("am_activation_rise", am_activation_rise);
    write_double_attr("am_activation_decay", am_activation_decay);
    write_double_attr("am_activation_floor", am_activation_floor);
    write_double_attr("am_activation_active_fraction",
                      am_activation_active_fraction);
    write_double_attr("am_k_on", am_k_on);
    write_double_attr("am_weak_binding_floor", am_weak_binding_floor);
    write_double_attr("am_k_weak_off", am_k_weak_off);
    write_double_attr("am_k_powerstroke", am_k_powerstroke);
    write_double_attr("am_k_strong_off", am_k_strong_off);
    write_double_attr("am_binding_azimuth_half_width_rad",
                      am_binding_azimuth_half_width);
    write_int_attr("am_binding_azimuth_selection",
                   am_binding_azimuth_selection ? 1 : 0);
    write_double_attr("am_unloaded_velocity", am_unloaded_velocity);
    write_string_attr("am_head_target_semantics",
                      "fixed_actin_material_coordinate_um");
    write_string_attr("am_head_mechanics",
                      "objective_two_spring_anisotropic_strain_kmc_v4");
    write_string_attr(
        "am_head_transverse_law",
        "per_bound_material_head_flat_bottom_transverse_spring_with_independent_escape_coupling");
    write_string_attr(
        "am_head_emergency_rupture_law",
        "connector_excess_beyond_relaxed_reach_only_after_stochastic_kmc");
    write_string_attr("am_head_crown_layout",
                      "axial_crowns_with_three_120deg_sites_and_40deg_twist");
    write_string_attr(
        "am_head_association_axial_gate",
        "material_head_projection_must_lie_inside_finite_actin_segment");
    write_double_attr("am_bare_zone_fraction",
                      geometry::am_bare_zone_fraction);
    write_string_attr("unit_length", "um");
    write_string_attr("unit_time", "s");
    write_string_attr("unit_force", "pN");
    write_string_attr("unit_energy", "pN um");
    write_string_attr("unit_stiffness", "pN/um");
    write_string_attr("unit_torque", "pN um");
    write_int_attr("am_thin_activation_enabled",
                   am_thin_activation_enabled ? 1 : 0);
    write_int_attr("am_spatial_heads", am_spatial_heads ? 1 : 0);
    write_int_attr("sterics_enabled", sterics_enabled ? 1 : 0);
    write_double_attr("myosin_radius", myosin.radius);
    write_double_attr("actin_radius", actin_radius);
    write_double_attr("am_surface_gap", am_surface_gap);
    write_double_attr("am_steric_k", am_steric_k);
    write_double_attr("actin_steric_k", actin_steric_k);
    write_double_attr("myosin_steric_k", myosin_steric_k);
    write_double_attr("myosin_exclusion_distance", myosin_exclusion_distance);
    write_int_attr("myomesin_schema_version", k_mm > 0.0 ? 1 : 0);
    write_double_attr("k_mm", k_mm);
    write_int_attr("myomesin_sites_per_filament",
                   myomesin_sites_per_filament);
    write_double_attr("myomesin_k_on", myomesin_k_on);
    write_double_attr("myomesin_k_off", myomesin_k_off);
    write_double_attr("myomesin_axial_capture", myomesin_axial_capture);
    write_double_attr("myomesin_radial_rest", myomesin_radial_rest);
    write_double_attr("myomesin_radial_capture", myomesin_radial_capture);
    write_double_attr("myomesin_radial_k", myomesin_radial_k);
    write_double_attr("myomesin_alignment_cos", myomesin_alignment_cos);
    write_double_attr("myomesin_bell_distance", myomesin_bell_distance);
    write_double_attr("myomesin_force_max", myomesin_force_max);
    write_string_attr(
        "myomesin_model",
        k_mm > 0.0
            ? "sparse_local_dynamic_finite_valency_bare_zone_modules"
            : "disabled");
    write_string_attr(
        "myomesin_force_law",
        "coarse_module_pairwise_axial_zero_register_plus_transverse_nonzero_rest_spacing");
    write_string_attr(
        "myomesin_bridge_state_columns",
        "id_hi,id_lo,myosin_i,myosin_j,attach_step,axial_offset,"
        "radial_distance,force");
    write_string_attr(
        "myomesin_events_columns",
        "step,event_code,id,myosin_i,myosin_j,axial_offset,"
        "radial_distance,force,age");
    write_string_attr("myomesin_event_codes",
                      "1=association,2=stochastic_release,"
                      "3=invalid_geometry_release");
    write_int_attr("aa_double_overlap_contact_enabled",
                   aa_double_overlap_contact_enabled ? 1 : 0);
    write_double_attr("aa_double_overlap_contact_diameter",
                      aa_double_overlap_contact_diameter);
    write_double_attr("aa_double_overlap_contact_epsilon",
                      aa_double_overlap_contact_epsilon);
    write_double_attr("aa_double_overlap_contact_force_cap",
                      aa_double_overlap_contact_force_cap);
    write_string_attr(
        "aa_double_overlap_contact_law",
        "local_antiparallel_pointed_end_unilateral_WCA_equal_opposite_force_and_torque");
    write_double_attr("am_preferred_centerline", am_preferred_centerline);
    write_double_attr("aa_preferred_centerline", aa_optimal);
    write_double_attr("aa_spacing_ratio", aa_spacing_ratio);
    write_double_attr("aa_capture_margin", aa_capture_margin);
    write_double_attr("aa_antiparallel_overlap", aa_antiparallel_overlap);
    write_double_attr("aa_axial_register_ratio", aa_axial_register_ratio);
    write_int_attr("aa_register_all_statuses",
                   aa_register_all_statuses ? 1 : 0);
    write_double_attr("aa_axial_register_force_max",
                      aa_axial_register_force_max);
    write_double_attr("aa_barbed_zone", aa_barbed_zone);
    write_double_attr("aa_status2_registry_on", aa_status2_registry_on);
    write_double_attr("aa_status2_registry_off", aa_status2_registry_off);
    write_double_attr(
        "aa_register_saturation_error", aa_register_saturation_error);
    write_int_attr("aa_local_multilink", aa_local_multilink ? 1 : 0);
    write_int_attr("aa_prebundle_links_per_pair", aa_prebundle_links_per_pair);
    write_double_attr("aa_site_exclusion", aa_site_exclusion);
    write_double_attr("aa_overlap_capture", aa_overlap_capture);
    write_double_attr("aa_perpendicular_k", aa_perpendicular_k);
    write_double_attr("aa_status1_radial_k", aa_status1_radial_k);
    write_double_attr("aa_status1_site_diffusion", aa_status1_site_diffusion);
    write_int_attr("aa_status2_locked_until_detach",
                   aa_status2_locked_until_detach ? 1 : 0);
    write_double_attr("aa_status2_unlock_rate", aa_status2_unlock_rate);
    write_double_attr("aa_status2_unlock_force_scale",
                      aa_status2_unlock_force_scale);
    write_double_attr("aa_status2_unlock_registry_scale",
                      aa_status2_unlock_registry_scale);
    write_int_attr("aa_arm_state_model", aa_arm_state_model ? 1 : 0);
    write_double_attr("aa_arm_open_rate_weak", aa_arm_open_rate_weak);
    write_double_attr("aa_arm_open_rate_mature", aa_arm_open_rate_mature);
    write_double_attr("aa_arm_bell_distance", aa_arm_bell_distance);
    write_double_attr("aa_tether_off_rate", aa_tether_off_rate);
    write_double_attr("aa_tether_rebind_rate", aa_tether_rebind_rate);
    write_double_attr("aa_tether_capture_distance",
                      aa_tether_capture_distance);
    write_double_attr("aa_arm_rebind_theta", aa_arm_rebind_theta);
    write_double_attr("aa_actn2_contact_radius", aa_actn2_contact_radius);
    write_double_attr("aa_myosin_contact_epsilon",
                      aa_myosin_contact_epsilon);
    write_double_attr("aa_myosin_contact_force_cap",
                      aa_myosin_contact_force_cap);
    write_double_attr("aa_myosin_contact_shell_width",
                      aa_myosin_contact_shell_width);
    write_double_attr("aa_myosin_contact_shell_energy",
                      aa_myosin_contact_shell_energy);
    write_int_attr("aa_zdisc_tip_contact", aa_zdisc_tip_contact ? 1 : 0);
    write_double_attr("aa_zdisc_tip_onset", aa_zdisc_tip_onset);
    write_double_attr("aa_zdisc_tip_stiffness", aa_zdisc_tip_stiffness);
    write_double_attr("aa_zdisc_tip_force_max", aa_zdisc_tip_force_max);
    write_double_attr("aa_barbed_lateral_range", aa_barbed_lateral_range);
    write_double_attr("aa_barbed_lateral_k", aa_barbed_lateral_k);
    write_double_attr("aa_barbed_lateral_force_max",
                      aa_barbed_lateral_force_max);
    write_double_attr("aa_barbed_lateral_axial_window",
                      aa_barbed_lateral_axial_window);
    write_double_attr("aa_zdisc_tip_transverse_flat",
                      aa_zdisc_tip_transverse_flat);
    write_double_attr("aa_zdisc_tip_transverse_cutoff",
                      aa_zdisc_tip_transverse_cutoff);
    write_string_attr("aa_tip_accessibility_scope",
                      "all_nearby_myosin_tips_head_independent");
    write_string_attr("aa_zdisc_tip_nonpenetration",
                      "mobility_weighted_projected_gauss_seidel_g_ge_0");
    write_string_attr("aa_zdisc_tip_projection_traversal",
                      "spatial_incidence_exact_conservative_broadphase");
    write_double_attr("aa_occlusion_rupture_penetration",
                      aa_occlusion_rupture_penetration);
    write_double_attr(
        "aa_status1_radial_k_effective",
        aa_link_radial_stiffness(static_cast<BondFlag>(1)));
    write_double_attr("aa_association_max_direction_dot",
                      aa_association_max_direction_dot);
    write_double_attr("aa_link_rupture_length", aa_link_rupture_length);
    write_int_attr("aa_stiffness_reference_links",
                   aa_stiffness_reference_links);
    write_double_attr(
        "aa_single_link_k",
        actn2_load_kinetics::fixed_single_link_stiffness(
            k_aa, aa_stiffness_reference_links));
    write_double_attr(
        "aa_single_link_kappa",
        actn2_load_kinetics::fixed_single_link_stiffness(
            kappa_aa, aa_stiffness_reference_links));
    write_int_attr("aa_total_link_pool", aa_total_link_pool);
    write_int_attr("aa_neighbor_identity",
                   aa_neighbor_identity_enabled ? 1 : 0);
    write_double_attr("aa_neighbor_slot_memory",
                      aa_neighbor_slot_memory_steps * dt);
    write_double_attr("aa_status2_slot_memory",
                      aa_status2_slot_memory_steps * dt);
    write_int_attr("aa_local_link_schema",
                   aa_local_multilink ? (aa_arm_state_model ? 5 : 4) : 0);
    write_double_attr("titin_k", titin_k);
    write_double_attr("titin_rest_length", titin_rest_length);
    write_int_attr("titin_mode", static_cast<int>(titin_mode));
    // Schema 3 adds the ACTN2-support capture gate and the per-frame
    // association-provenance counters. The bump is deliberate: it makes the
    // resume guard reject every schema-2 trajectory, including the cancelled
    // saturated v4 pilot, rather than silently continuing it under new rules.
    write_int_attr("titin_schema_version",
                   titin_mode == TitinMode::explicit_dynamic ? 4 : 0);
    write_double_attr("titin_k_off_unsupported",
                      explicit_titin_k_off_unsupported);
    write_int_attr("titin_rebind_requires_actn2",
                   explicit_titin_rebind_requires_actn2 ? 1 : 0);
    write_string_attr(
        "titin_anchor_topology",
        titin_mode != TitinMode::explicit_dynamic
            ? "not_applicable"
            : (explicit_titin_k_off_unsupported > 0.0
                ? "serial_actin_actn2_titin_composite_lifetime"
                : "parallel_single_k_off"));
    write_int_attr("titin_require_live_actn2",
                   explicit_titin_require_live_actn2 ? 1 : 0);
    write_string_attr(
        "titin_capture_gate",
        titin_mode != TitinMode::explicit_dynamic
            ? "not_applicable"
            : (explicit_titin_require_live_actn2
                ? "formation_requires_live_actn2_on_barbed_end_any_status"
                : "unrestricted_within_capture_shell"));
    write_string_attr(
        "titin_association_counts_columns",
        "actn2_supported,actn2_unsupported,same_partner_rebind,"
        "unsupported_releases");
    write_string_attr(
        "titin_event_codes",
        "1=association,2=release_with_live_actn2,3=release_without_live_actn2");
    write_int_attr("titin_links_per_half", explicit_titin_per_half);
    write_double_attr("titin_capture_distance",
                      explicit_titin_capture_distance);
    write_double_attr("titin_k_on", explicit_titin_k_on);
    write_double_attr("titin_k_off", explicit_titin_k_off);
    write_int_attr("titin_strain_gated_kinetics",
                   explicit_titin_strain_gated_kinetics ? 1 : 0);
    write_string_attr("titin_association_law",
                      "distance_limited_constant_rate_then_tension_only_spring");
    write_double_attr("titin_transition_theta",
                      explicit_titin_transition_theta);
    write_double_attr("titin_rebind_memory",
                      explicit_titin_rebind_memory);
    write_double_attr("titin_rebind_rate_multiplier",
                      explicit_titin_rebind_rate_multiplier);
    write_string_attr("titin_force_law",
                      titin_mode == TitinMode::explicit_dynamic
                          ? "tension_only_flat_bottom"
                          : "legacy_bidirectional_hookean");
    write_string_attr("titin_capture_semantics",
                      titin_mode == TitinMode::explicit_dynamic
                          ? "material_tip_domain_to_individual_actin_barbed_end"
                          : "legacy_effective_centering");
    write_string_attr(
        "titin_topology",
        titin_mode == TitinMode::explicit_dynamic
            ? "explicit_thick_tip_domains_to_dynamic_actin_barbed_ends"
            : (titin_mode == TitinMode::balanced_owning
                ? "balanced_barbed_end_to_owning_myosin_center"
                : "legacy_status2_crosslinked_partner"));
    write_int_attr("release_step_origin", static_cast<int>(
        std::min(release_step_origin,
                 static_cast<size_t>(std::numeric_limits<int>::max()))));
    // The release origin is exactly the number of fixed-geometry chemistry
    // steps requested by main. Store the protocol name explicitly so a
    // production checkpoint can prove that chemistry and motion were coupled
    // from step zero without relying on that implementation detail.
    write_int_attr("fixed_lattice_preequil_steps", static_cast<int>(
        std::min(release_step_origin,
                 static_cast<size_t>(std::numeric_limits<int>::max()))));
    write_string_attr(
        "aa_local_link_mechanics_columns",
        "connector_x,connector_y,connector_z,length,extension,arm_i,arm_j,"
        "bilateral_load,productive_support,status2_gate,lifetime,single_link_tension,"
        "raw_spring_arm_i,raw_spring_arm_j,allocated_motor_i,allocated_motor_j");
    write_string_attr(
        "aa_local_link_mechanics_v2_columns",
        "connector_x,connector_y,connector_z,length,extension,arm_i,arm_j,"
        "bilateral_load,productive_support,status2_gate,lifetime,single_link_tension,"
        "raw_spring_arm_i,raw_spring_arm_j,allocated_motor_i,allocated_motor_j,"
        "direction_dot,registry_error,registry_threshold,registry_ok,arm_load_ok,"
        "productive_support_ok,actual_status2_gate");
    write_string_attr(
        "aa_local_link_mechanics_v2_frame_offsets_columns",
        "step,identity_start,mechanics_v2_start,count");
    write_string_attr(
        "aa_local_link_arm_state_columns",
        "state,state_since_step,last_partner,arm_i_bound,arm_j_bound,"
        "arm_force_i,arm_force_j,actn2_myosin_contact_force,"
        "actn2_myosin_contact_penetration");
    write_string_attr(
        "aa_arm_state_codes",
        "1=weak_two_arm_bridge,2=mature_two_arm_bridge,"
        "3=tethered_to_actin_i,4=tethered_to_actin_j");
    write_string_attr(
        "aa_arm_events_columns",
        "link_id,step,event_code,old_state,new_state,bound_actin,old_partner,"
        "new_or_opened_actin,bound_site,free_or_opened_site,arm_force,"
        "contact_force,contact_penetration,event_rate");
    write_string_attr(
        "aa_arm_event_codes",
        "1=one_arm_opened,2=free_arm_rebound,3=tether_detached");
    write_string_attr(
        "aa_local_link_motor_support_columns",
        "candidate_myosin_i,candidate_myosin_j,candidate_half_i,"
        "candidate_half_j,support_ema,support_instantaneous,candidate_age,"
        "identity_locked,candidate_since_step,candidate_last_seen_step,"
        "raw_source_support_i,raw_source_support_j,"
        "allocated_source_support_i,allocated_source_support_j,"
        "challenger_myosin_i,challenger_myosin_j,challenger_half_i,"
        "challenger_half_j,challenger_support_ema,challenger_age");
    write_string_attr(
        "aa_motor_pair_events_columns",
        "link_id_hi,link_id_lo,step,event_code,old_myosin_i,old_myosin_j,"
        "old_half_i,old_half_j,new_myosin_i,new_myosin_j,new_half_i,"
        "new_half_j,current_support,best_support,support_ema,candidate_age,"
        "raw_source_support_i,raw_source_support_j,"
        "allocated_source_support_i,allocated_source_support_j");
    write_string_attr(
        "aa_motor_pair_event_codes",
        "1=selected,2=switched,3=dropped_after_support_timeout,"
        "4=locked_at_maturation,5=released_by_actn2_detachment,"
        "6=unlocked_to_weak_state");
    write_string_attr(
        "aa_motor_profile_semantics",
        "local_update_is_wall_time;source_build_and_candidate_scan_are_sums_"
        "of_worker_thread_time;printed_at_process_exit_for_motor_transmitted_register");
    write_string_attr(
        "aa_local_link_events_v2_columns",
        "link_id,actin_i,actin_j,step,event_code,old_status,new_status,"
        "effective_load,lifetime,s_i,s_j,arm_i,arm_j,bilateral_load,"
        "productive_support,status2_gate,single_link_tension,raw_spring_arm_i,"
        "raw_spring_arm_j,allocated_motor_i,allocated_motor_j");
    write_string_attr(
        "aa_local_link_events_v3_columns",
        "link_id,actin_i,actin_j,step,event_code,old_status,new_status,"
        "effective_load,lifetime,s_i,s_j,arm_i,arm_j,bilateral_load,"
        "productive_support,status2_gate,single_link_tension,raw_spring_arm_i,"
        "raw_spring_arm_j,allocated_motor_i,allocated_motor_j,connector_length,"
        "rupture_length,direction_dot,registry_error,registry_threshold,registry_ok,"
        "arm_load_ok,productive_support_ok,actual_status2_gate,detachment_reason,"
        "break_probability,stochastic_draw,geometric_detachment,"
        "stochastic_detachment");
    write_string_attr(
        "aa_local_link_detachment_reason_codes",
        "0=not_a_detachment,1=geometric_contour,2=stochastic,"
        "3=am_overlap_occlusion,4=actn2_myosin_collision");
    write_string_attr(
        "aa_tip_accessibility_events_columns",
        "link_id,step,event_code,actin_i,actin_j,s_i,s_j,nearest_myosin,"
        "tip_sign,axial_gap,transverse_distance,tip_zone_occupied,"
        "motor_engaged");
    write_string_attr(
        "aa_tip_accessibility_event_codes",
        "1=accepted_association,2=accepted_tether_rebinding");
    write_string_attr(
        "aa_catch_load_semantics",
        "per_actin_productive_motor_load_allocated_across_transmitting_arms_"
        "and_bounded_by_raw_arm_force");
    write_string_attr(
        "am_forced_ruptures_v2_columns",
        "actin,myosin,step,head_id,old_state,reason,failure_mask,"
        "segment_distance,rupture_distance,actin_myosin_dot,expected_half,"
        "head_half,productive_axial_overlap,binding_ratio,partial_ratio,"
        "head_connector_distance,head_surface_distance,capture_reach,"
        "retention_reach");
    write_int_attr("physical_rod_mobility", physical_rod_mobility ? 1 : 0);
    write_int_attr("isometric_actin_axial", isometric_actin_axial ? 1 : 0);
    write_double_attr("temperature_kelvin", temperature_kelvin);
    write_double_attr("thermal_beta", thermal_beta);
    write_double_attr("legacy_actin_D_trans", legacy_actin_D_trans);
    write_double_attr("legacy_actin_D_rot", legacy_actin_D_rot);
    write_double_attr("legacy_myosin_D_trans", legacy_myosin_D_trans);
    write_double_attr("legacy_myosin_D_rot", legacy_myosin_D_rot);
    write_double_attr("thermal_noise_scale", thermal_noise_scale_metadata);
    write_double_attr("actin_length", actin.length);
    write_double_attr("myosin_length", myosin.length);
    const double expected_beta = rod_mobility::beta(temperature_kelvin);
    write_int_attr("thermal_beta_overridden",
                   std::abs(thermal_beta - expected_beta) >
                       1e-10 * expected_beta ? 1 : 0);
    write_double_attr("aa_refractory_time", bond_recovery_steps * dt);
    write_double_attr("dt", dt);
    write_double_attr("k_am", k_am);
    write_double_attr("k_aa", k_aa);
    write_double_attr("kappa_am", kappa_am);
    write_double_attr("kappa_aa", kappa_aa);
    write_double_attr("aa_base_lifetime", base_lifetime);
    write_double_attr("aa_directional_base_lifetime", directional_base_lifetime);
    write_double_attr("aa_lifetime_coefficient", lifetime_coeff);
    write_double_attr("aa_motor_pair_support_window",
                      aa_motor_pair_support_window);
    write_double_attr("aa_motor_pair_commit_time",
                      aa_motor_pair_commit_time);
    write_double_attr("aa_motor_pair_switch_ratio",
                      aa_motor_pair_switch_ratio);
    write_double_attr("aa_motor_pair_exchange_time",
                      aa_motor_pair_exchange_time);
    write_double_attr("aa_motor_pair_dropout_time",
                      aa_motor_pair_dropout_time);
    write_string_attr(
        "aa_motor_transmitted_register_semantics",
        "distinct_myosin_pair_opposes_actual_local_actn2_connector_on_both_arms;"
        "weak_state_exchange_then_identity_lock_at_maturation");
    write_int_attr("max_aa_bonds_per_actin", max_aa_bonds_per_actin);
    write_int_attr("aa_energy_gated_association",
                   aa_energy_gated_association ? 1 : 0);
    write_double_attr("aa_transition_state_theta", aa_transition_state_theta);
    write_int_attr("aa_passive_detailed_balance",
                   aa_passive_detailed_balance ? 1 : 0);
    write_double_attr("aa_encounter_k", aa_encounter_k);
    write_double_attr("aa_encounter_kappa", aa_encounter_kappa);
    write_int_attr("interaction_ramp_start_step",
                   static_cast<int>(interaction_ramp_start_step));
    write_int_attr("interaction_ramp_steps",
                   static_cast<int>(interaction_ramp_steps));
    auto write_mobility = [&](const char* prefix,
                              const rod_mobility::Diffusion& mobility) {
        const std::string base(prefix);
        write_double_attr((base + "_hydrodynamic_length").c_str(), mobility.length);
        write_double_attr((base + "_hydrodynamic_diameter").c_str(), mobility.diameter);
        write_double_attr((base + "_viscosity_pas").c_str(), mobility.viscosity);
        write_double_attr((base + "_crowding_parallel").c_str(), mobility.crowding_parallel);
        write_double_attr((base + "_crowding_perpendicular").c_str(), mobility.crowding_perpendicular);
        write_double_attr((base + "_crowding_tumble").c_str(), mobility.crowding_tumble);
        write_double_attr((base + "_crowding_roll").c_str(), mobility.crowding_roll);
        write_double_attr((base + "_D_parallel").c_str(), mobility.parallel);
        write_double_attr((base + "_D_perpendicular").c_str(), mobility.perpendicular);
        write_double_attr((base + "_D_tumble").c_str(), mobility.tumble);
        write_double_attr((base + "_D_roll").c_str(), mobility.roll);
    };
    write_mobility("actin", actin_rod_mobility);
    write_mobility("myosin", myosin_rod_mobility);
    write_int_attr("am_thin_feedback_enabled",
                   am_thin_feedback_enabled ? 1 : 0);
    write_double_attr("am_q_calcium_open_rate", am_q_calcium_open_rate);
    write_double_attr("am_q_calcium_half", am_q_calcium_half);
    write_double_attr("am_q_calcium_hill", am_q_calcium_hill);
    write_double_attr("am_q_strong_feedback_rate", am_q_strong_feedback_rate);
    write_double_attr("am_q_strong_half", am_q_strong_half);
    write_double_attr("am_q_strong_hill", am_q_strong_hill);
    write_double_attr("am_q_close_rate", am_q_close_rate);
    write_int_attr("am_q_initialized_at_calcium_equilibrium",
                   am_q_initialized_at_calcium_equilibrium ? 1 : 0);
    write_double_attr("am_initial_q", am_initial_q);
    write_int_attr("am_thick_feedback_enabled",
                   am_thick_feedback_enabled ? 1 : 0);
    write_double_attr("am_thick_recruit_rate", am_thick_recruit_rate);
    write_double_attr("am_thick_fold_rate", am_thick_fold_rate);
    write_double_attr("am_thick_tension_scale", am_thick_tension_scale);
    write_double_attr("am_thick_initial_on_fraction",
                      am_thick_initial_on_fraction);
    write_double_attr("am_thick_force_clamp", am_thick_force_clamp);
    file.flush(H5F_SCOPE_GLOBAL);
}
