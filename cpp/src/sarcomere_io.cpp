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

namespace {
// Moved from sarcomere.cpp with the persistence functions on 2026-09-17.
// Verified before the move that none of these is used elsewhere: of the 13
// names in that anonymous namespace, only boltzmann_exclusion_factor is
// referenced outside the moved block, and it stayed behind.
void delete_link_if_exists(H5::H5File& file, const std::string& full_path) {
    if (file.nameExists(full_path)) {
        H5Ldelete(file.getId(), full_path.c_str(), H5P_DEFAULT);
    }
}

void move_link(H5::H5File& file, const std::string& src_path, const std::string& dst_path) {
    H5Lmove(file.getId(), src_path.c_str(), file.getId(), dst_path.c_str(), H5P_DEFAULT, H5P_DEFAULT);
}

void promote_temp_group(H5::H5File& file,
                        const std::string& final_path,
                        const std::string& temp_path,
                        const std::string& backup_path) {
    delete_link_if_exists(file, backup_path);
    if (file.nameExists(final_path)) {
        move_link(file, final_path, backup_path);
    }
    move_link(file, temp_path, final_path);
    delete_link_if_exists(file, backup_path);
    file.flush(H5F_SCOPE_GLOBAL);
}

void maybe_test_sleep(const char* env_name, const char* label) {
    const char* raw = std::getenv(env_name);
    if (raw == nullptr || raw[0] == '\0') {
        return;
    }
    const int sleep_seconds = std::max(0, std::atoi(raw));
    if (sleep_seconds <= 0) {
        return;
    }
    fflush(stdout);
    std::this_thread::sleep_for(std::chrono::seconds(sleep_seconds));
}

void truncate_dataset_first_dim(H5::Group& group, const std::string& dataset_name, hsize_t keep_frames) {
    if (!group.nameExists(dataset_name)) {
        return;
    }
    H5::DataSet dataset = group.openDataSet(dataset_name);
    H5::DataSpace filespace = dataset.getSpace();
    const int rank = filespace.getSimpleExtentNdims();
    if (rank <= 0) {
        return;
    }
    std::vector<hsize_t> dims(rank, 0);
    filespace.getSimpleExtentDims(dims.data(), nullptr);
    if (dims.empty() || dims[0] <= keep_frames) {
        return;
    }
    dims[0] = keep_frames;
    dataset.extend(dims.data());
}

struct PackedAABondState {
    std::vector<int> pairs;
    std::vector<int> status;
    std::vector<int> lifetime;
};

struct PackedRecoveryState {
    std::vector<int> pairs;
    std::vector<int> until;
};

template <typename BondMatrix, typename StatusMatrix>
PackedAABondState pack_aa_bond_state(
    const std::vector<std::pair<int, int>>& active_pairs,
    const BondMatrix& bonds,
    const StatusMatrix& status,
    const std::vector<std::vector<int>>& lifetime) {
    PackedAABondState packed;
    packed.pairs.reserve(active_pairs.size() * 2);
    packed.status.reserve(active_pairs.size());
    packed.lifetime.reserve(active_pairs.size());
    for (const auto& [i, j] : active_pairs) {
        if (bonds[i][j] != 1) {
            continue;
        }
        packed.pairs.push_back(i);
        packed.pairs.push_back(j);
        packed.status.push_back(status[i][j]);
        packed.lifetime.push_back(lifetime[i][j]);
    }
    return packed;
}

template <typename BondMatrix>
std::vector<int> pack_am_bond_pairs(
    const std::vector<std::pair<int, int>>& active_pairs,
    const BondMatrix& bonds) {
    std::vector<int> packed;
    packed.reserve(active_pairs.size() * 2);
    for (const auto& [i, j] : active_pairs) {
        if (bonds[i][j] != 1) {
            continue;
        }
        packed.push_back(i);
        packed.push_back(j);
    }
    return packed;
}

template <typename RecoveryMatrix>
PackedRecoveryState pack_recovery_state(
    const std::vector<std::pair<int, int>>& active_pairs,
    const RecoveryMatrix& recovery_until,
    size_t current_step) {
    PackedRecoveryState packed;
    packed.pairs.reserve(active_pairs.size() * 2);
    packed.until.reserve(active_pairs.size());
    for (const auto& [i, j] : active_pairs) {
        const size_t until = static_cast<size_t>(recovery_until[i][j]);
        if (until <= current_step) {
            continue;
        }
        packed.pairs.push_back(i);
        packed.pairs.push_back(j);
        packed.until.push_back(
            static_cast<int>(std::min(until, static_cast<size_t>(std::numeric_limits<int>::max()))));
    }
    return packed;
}

void pad_int_vector(std::vector<int>& values, size_t target_size, int pad_value) {
    // Fit to exactly target_size: pad with pad_value if shorter, TRUNCATE if longer.
    // Truncation lets the fixed-width compact /state bond snapshot degrade gracefully when the
    // bond count exceeds capacity (dense, low-myosin networks) instead of overflowing the dataset.
    // Trajectory data lives in separate datasets and is unaffected.
    values.resize(target_size, pad_value);
}

std::vector<double> flatten_vec_values(const std::vector<vec>& values) {
    std::vector<double> flat;
    flat.reserve(3 * values.size());
    for (const auto& value : values) {
        flat.push_back(value.x);
        flat.push_back(value.y);
        flat.push_back(value.z);
    }
    return flat;
}

}  // namespace


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
    write_double_attr("minimal_catch_tilt_deg", minimal_catch_tilt_deg);
    write_double_attr("minimal_catch_myosin_tilt_deg",
                      minimal_catch_myosin_tilt_deg);
    write_double_attr("am_spokes_per_crown",
                      static_cast<double>(am_spokes_per_crown));
    write_int_attr("am_head_azimuth", am_head_azimuth ? 1 : 0);
    write_int_attr("legacy_actin_torque_cap", legacy_actin_torque_cap ? 1 : 0);
    write_double_attr("actin_direction_polar_fraction",
                      actin_direction_polar_fraction);
    write_double_attr("myosin_direction_polar_fraction",
                      myosin_direction_polar_fraction);
    // The rotation caps set the torque ceiling
    // (max_actin_torque = max_rotation / (beta * D_rot * dt)) and were not
    // recorded, so analyses had to back-derive them.  Record them.
    // The rotation caps themselves live in Langevin, but the ceilings they
    // derive are what analysis needs: max_rotation = torque * beta*D_rot*dt.
    write_double_attr("max_actin_torque", max_actin_torque);
    write_double_attr("max_myosin_torque", max_myosin_torque);
    write_double_attr("max_actin_force", max_actin_force);
    write_double_attr("max_myosin_force", max_myosin_force);
    write_int_attr("torque_breakdown", torque_breakdown ? 1 : 0);
    write_double_attr("am_head_radial_force_max", am_head_radial_force_max);
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


// save_state, save_resume_snapshot and load_state moved here from
// sarcomere.cpp on 2026-09-17.  Byte-for-byte move, no logic change.

void Sarcomere::save_state(){
    auto bondData =
        _extract_bonded_pairs(actin_actin_bonds, actin_actin_status, myosinIndicesPerActin);
    std::vector<double> flatActinBonds = std::get<0>(bondData);
    std::vector<double> flatMyosinBonds = std::get<1>(bondData);
    std::vector<double> flatActinMyosinBonds = std::get<2>(bondData);
    append_to_file(filename, actin, myosin, flatActinBonds,
                   flatMyosinBonds, flatActinMyosinBonds, max_myosin_bonds);
    
    // Save compact bond state for resume using pair lists instead of dense matrices.
    {
        H5::H5File file(filename, H5F_ACC_RDWR);
        H5::Group group_state;
        try {
            group_state = file.openGroup("/state");
        } catch (H5::Exception&) {
            group_state = file.createGroup("/state");
        }
        auto ensure_state_dataset_int = [&](const std::string& dataset_name, hsize_t width) {
            const std::string full_path = "/state/" + dataset_name;
            if (!file.nameExists(full_path)) {
                std::vector<hsize_t> initDims = {0, width};
                std::vector<hsize_t> maxDims = {H5S_UNLIMITED, width};
                std::vector<hsize_t> chunkDims = {10, width};
                create_empty_dataset_int(file, "/state", dataset_name, initDims, maxDims, chunkDims);
            }
        };
        auto ensure_state_dataset_double = [&](const std::string& dataset_name, hsize_t width) {
            const std::string full_path = "/state/" + dataset_name;
            if (!file.nameExists(full_path)) {
                std::vector<hsize_t> initDims = {0, width};
                std::vector<hsize_t> maxDims = {H5S_UNLIMITED, width};
                std::vector<hsize_t> chunkDims = {10, width};
                create_empty_dataset(file, "/state", dataset_name, initDims, maxDims, chunkDims);
            }
        };

        const hsize_t actin_width = static_cast<hsize_t>(actin.n);
        const hsize_t myosin_width = static_cast<hsize_t>(myosin.n);
        const hsize_t max_aa_pairs = static_cast<hsize_t>(actin.n) *
                                     static_cast<hsize_t>(std::max(10, 2 * max_myosin_bonds));
        const hsize_t max_am_pairs = static_cast<hsize_t>(myosin.n) *
                                     static_cast<hsize_t>(max_myosin_bonds);
        const hsize_t max_recovery_pairs = max_aa_pairs;

        auto ensure_state_pair_dataset = [&](const std::string& dataset_name, hsize_t rows, hsize_t cols) {
            const std::string full_path = "/state/" + dataset_name;
            if (!file.nameExists(full_path)) {
                std::vector<hsize_t> initDims = {0, rows, cols};
                std::vector<hsize_t> maxDims = {H5S_UNLIMITED, rows, cols};
                std::vector<hsize_t> chunkDims = {10, rows, cols};
                create_empty_dataset_int(file, "/state", dataset_name, initDims, maxDims, chunkDims);
            }
        };
        auto ensure_state_pair_dataset_double = [&](const std::string& dataset_name, hsize_t rows, hsize_t cols) {
            const std::string full_path = "/state/" + dataset_name;
            if (!file.nameExists(full_path)) {
                std::vector<hsize_t> initDims = {0, rows, cols};
                std::vector<hsize_t> maxDims = {H5S_UNLIMITED, rows, cols};
                std::vector<hsize_t> chunkDims = {10, rows, cols};
                create_empty_dataset(file, "/state", dataset_name, initDims, maxDims, chunkDims);
            }
        };

        ensure_state_dataset_int("current_step", 1);
        ensure_state_dataset_int("state_commit_count", 1);
        ensure_state_dataset_int("aa_pairs_prev_count", 1);
        ensure_state_dataset_int("aa_pairs_current_count", 1);
        ensure_state_dataset_int("am_pairs_prev_count", 1);
        ensure_state_dataset_int("am_pairs_current_count", 1);
        ensure_state_pair_dataset("aa_pairs_prev", max_aa_pairs, 2);
        ensure_state_pair_dataset("aa_status_prev", max_aa_pairs, 1);
        ensure_state_pair_dataset("aa_lifetime_prev", max_aa_pairs, 1);
        ensure_state_pair_dataset("aa_pairs_current", max_aa_pairs, 2);
        ensure_state_pair_dataset("aa_status_current", max_aa_pairs, 1);
        ensure_state_pair_dataset("aa_lifetime_current", max_aa_pairs, 1);
        ensure_state_pair_dataset("aa_attach_step_current", max_aa_pairs, 1);
        ensure_state_pair_dataset("aa_status2_memory_until_current", max_aa_pairs, 1);
        ensure_state_pair_dataset("am_pairs_prev", max_am_pairs, 2);
        ensure_state_pair_dataset("am_pairs_current", max_am_pairs, 2);
        ensure_state_dataset_int("am_state_current", max_am_pairs);
        ensure_state_dataset_double("am_target_axial_current", max_am_pairs);
        const hsize_t heads_per_myosin =
            static_cast<hsize_t>(2 * am_heads_per_half);
        ensure_state_pair_dataset("am_head_state_current", myosin_width,
                                  heads_per_myosin);
        ensure_state_pair_dataset("am_head_bound_actin_current", myosin_width,
                                  heads_per_myosin);
        ensure_state_pair_dataset_double("am_head_target_axial_current", myosin_width,
                                         heads_per_myosin);
        ensure_state_pair_dataset_double("am_head_binding_azimuth_current",
                                         myosin_width, heads_per_myosin);
        ensure_state_pair_dataset_double("am_head_stroke_progress_current",
                                         myosin_width, heads_per_myosin);
        ensure_state_pair_dataset_double("am_head_strain_energy_current",
                                         myosin_width, heads_per_myosin);
        ensure_state_pair_dataset_double("am_head_detachment_rate_current",
                                         myosin_width, heads_per_myosin);
        ensure_state_pair_dataset_double("myosin_material_e1_current",
                                         myosin_width, 3);
        ensure_state_pair_dataset_double("myosin_material_e2_current",
                                         myosin_width, 3);
        ensure_state_dataset_int("aa_recovery_count", 1);
        ensure_state_pair_dataset("aa_recovery_pairs", max_recovery_pairs, 2);
        ensure_state_pair_dataset("aa_recovery_until_values", max_recovery_pairs, 1);
        ensure_state_dataset_double("neighbor_last_actin_x", actin_width);
        ensure_state_dataset_double("neighbor_last_actin_y", actin_width);
        ensure_state_dataset_double("neighbor_last_actin_z", actin_width);
        ensure_state_dataset_double("neighbor_last_myosin_x", myosin_width);
        ensure_state_dataset_double("neighbor_last_myosin_y", myosin_width);
        ensure_state_dataset_double("neighbor_last_myosin_z", myosin_width);
        ensure_state_dataset_double("am_q_current", actin_width);
        ensure_state_dataset_int("actin_last_turnover_step_current",
                                 actin_width);
        const hsize_t titin_width =
            static_cast<hsize_t>(explicit_titins.size());
        if (titin_mode == TitinMode::explicit_dynamic) {
            ensure_state_dataset_int("titin_bound_actin_current", titin_width);
            ensure_state_dataset_int("titin_last_bound_actin_current", titin_width);
            ensure_state_dataset_int("titin_attach_step_current", titin_width);
            ensure_state_dataset_int("titin_last_detach_step_current", titin_width);
            ensure_state_dataset_double("titin_distance_current", titin_width);
            ensure_state_dataset_double("titin_extension_current", titin_width);
            ensure_state_dataset_double("titin_force_current", titin_width);
            ensure_state_dataset_double(
                "titin_attachment_distance_current", titin_width);
            ensure_state_dataset_double(
                "titin_attachment_energy_current", titin_width);
            ensure_state_dataset_double(
                "titin_attachment_force_current", titin_width);
            ensure_state_dataset_int("titin_association_counts", 4);
        }

        // Save current_step (scalar value)
        std::vector<int> step_vec = {static_cast<int>(current_step)};
        append_to_dataset_int(group_state, "current_step", step_vec, {1, 1});

        auto aa_prev = pack_aa_bond_state(
            previous_aa_pairs, actin_actin_bonds_prev,
            actin_actin_status_prev, actin_actin_lifetime_prev);
        auto aa_current = pack_aa_bond_state(
            active_aa_pairs, actin_actin_bonds,
            actin_actin_status, actin_actin_lifetime);
        auto am_prev = pack_am_bond_pairs(previous_am_pairs, am_bonds_prev);
        auto am_current = pack_am_bond_pairs(active_am_pairs, am_bonds);
        auto aa_recovery = pack_recovery_state(
            active_recovery_pairs, actin_recovery_until, current_step);

        const hsize_t aa_prev_count = static_cast<hsize_t>(aa_prev.status.size());
        const hsize_t aa_current_count = static_cast<hsize_t>(aa_current.status.size());
        const hsize_t am_prev_count = static_cast<hsize_t>(am_prev.size() / 2);
        const hsize_t am_current_count = static_cast<hsize_t>(am_current.size() / 2);
        const hsize_t aa_recovery_count = static_cast<hsize_t>(aa_recovery.until.size());
        if (aa_prev_count > max_aa_pairs || aa_current_count > max_aa_pairs ||
            am_prev_count > max_am_pairs || am_current_count > max_am_pairs ||
            aa_recovery_count > max_recovery_pairs) {
            // Degrade gracefully: the compact /state bond snapshot is truncated to capacity
            // (see pad_int_vector) and its stored counts are clamped below. Trajectory and
            // /actin per-frame fields are unaffected; only a later resume would see the
            // over-capacity bonds dropped from this checkpoint.
            static bool warned_state_capacity = false;
            if (!warned_state_capacity) {
                fprintf(stderr,
                        "Warning: bond pair list exceeds compact /state capacity at step %llu "
                        "(aa<=%llu/%llu, am<=%llu/%llu, rec=%llu/%llu); truncating resume snapshot, "
                        "trajectory data unaffected.\n",
                        static_cast<unsigned long long>(current_step),
                        static_cast<unsigned long long>(std::max(aa_prev_count, aa_current_count)),
                        static_cast<unsigned long long>(max_aa_pairs),
                        static_cast<unsigned long long>(std::max(am_prev_count, am_current_count)),
                        static_cast<unsigned long long>(max_am_pairs),
                        static_cast<unsigned long long>(aa_recovery_count),
                        static_cast<unsigned long long>(max_recovery_pairs));
                warned_state_capacity = true;
            }
        }

        std::vector<int> aa_prev_pairs = aa_prev.pairs;
        std::vector<int> aa_prev_status = aa_prev.status;
        std::vector<int> aa_prev_lifetime = aa_prev.lifetime;
        std::vector<int> aa_current_pairs = aa_current.pairs;
        std::vector<int> aa_current_status = aa_current.status;
        std::vector<int> aa_current_lifetime = aa_current.lifetime;
        std::vector<int> aa_current_attach_step;
        std::vector<int> aa_current_status2_memory_until;
        aa_current_attach_step.reserve(aa_current.status.size());
        aa_current_status2_memory_until.reserve(aa_current.status.size());
        for (size_t idx = 0; idx + 1 < aa_current.pairs.size(); idx += 2) {
            const int a = aa_current.pairs[idx];
            const int b = aa_current.pairs[idx + 1];
            aa_current_attach_step.push_back(aa_attach_step[a][b]);
            aa_current_status2_memory_until.push_back(
                static_cast<int>(aa_status2_memory_until[a][b]));
        }
        std::vector<int> am_prev_pairs = am_prev;
        std::vector<int> am_current_pairs = am_current;
        std::vector<int> am_current_state;
        std::vector<double> am_current_target;
        std::vector<int> am_head_state_current;
        std::vector<int> am_head_bound_actin_current;
        std::vector<double> am_head_target_current;
        std::vector<double> am_head_binding_azimuth_current;
        std::vector<double> am_head_stroke_progress_current;
        std::vector<double> am_head_strain_energy_current;
        std::vector<double> am_head_detachment_rate_current;
        am_current_state.reserve(am_current_count);
        am_current_target.reserve(am_current_count);
        for (size_t idx = 0; idx + 1 < am_current.size(); idx += 2) {
            const int a = am_current[idx];
            const int m = am_current[idx + 1];
            am_current_state.push_back(dynamic_am ? static_cast<int>(am_state[a][m]) : 2);
            am_current_target.push_back(dynamic_am ? am_target_axial[a][m]
                                                   : _am_axial_coordinate(a, m));
        }
        const size_t total_head_count = static_cast<size_t>(myosin.n) *
                                        static_cast<size_t>(2 * am_heads_per_half);
        am_head_state_current.reserve(total_head_count);
        am_head_bound_actin_current.reserve(total_head_count);
        am_head_target_current.reserve(total_head_count);
        am_head_binding_azimuth_current.reserve(total_head_count);
        am_head_stroke_progress_current.reserve(total_head_count);
        am_head_strain_energy_current.reserve(total_head_count);
        am_head_detachment_rate_current.reserve(total_head_count);
        for (int m = 0; m < myosin.n; ++m) {
            const auto& myosin_halves =
                am_heads_by_myosin[static_cast<size_t>(m)];
            for (int half_index = 0; half_index < 2; ++half_index) {
                const auto& half_heads =
                    myosin_halves[static_cast<size_t>(half_index)];
                for (int h = 0; h < static_cast<int>(half_heads.size()); ++h) {
                    const auto& head = half_heads[static_cast<size_t>(h)];
                    am_head_state_current.push_back(static_cast<int>(head.state));
                    am_head_bound_actin_current.push_back(head.bound_actin);
                    am_head_target_current.push_back(head.target_axial);
                    am_head_binding_azimuth_current.push_back(
                        head.binding_azimuth);
                    am_head_stroke_progress_current.push_back(
                        head.stroke_progress);
                    const auto state = static_cast<am_motor_kinetics::State>(
                        head.state);
                    double strain_energy = 0.0;
                    double detachment_rate = 0.0;
                    if (head.bound_actin >= 0 &&
                        (state == am_motor_kinetics::State::weak ||
                         state == am_motor_kinetics::State::strong)) {
                        const int head_id = _am_head_id(half_index, h);
                        strain_energy = _am_head_strain_energy(
                            head.bound_actin, m, head_id, state);
                        double chemical_rate = am_k_weak_off;
                        if (state == am_motor_kinetics::State::strong) {
                            chemical_rate =
                                am_motor_kinetics::load_slowed_rate_with_floor(
                                am_k_strong_off,
                                _am_head_active_force(
                                    head.bound_actin, m, head_id),
                                am_force_scale,
                                am_strong_off_floor_fraction);
                        }
                        detachment_rate = _am_head_detachment_rate(
                            head.bound_actin, m, head_id, state,
                            chemical_rate);
                    }
                    am_head_strain_energy_current.push_back(strain_energy);
                    am_head_detachment_rate_current.push_back(detachment_rate);
                }
            }
        }
        std::vector<int> aa_recovery_pairs = aa_recovery.pairs;
        std::vector<int> aa_recovery_until_values = aa_recovery.until;
        std::vector<int> titin_bound_actin_current;
        std::vector<int> titin_last_bound_actin_current;
        std::vector<int> titin_attach_step_current;
        std::vector<int> titin_last_detach_step_current;
        std::vector<double> titin_distance_current;
        std::vector<double> titin_extension_current;
        std::vector<double> titin_force_current;
        std::vector<double> titin_attachment_distance_current;
        std::vector<double> titin_attachment_energy_current;
        std::vector<double> titin_attachment_force_current;
        if (titin_mode == TitinMode::explicit_dynamic) {
            for (const auto& titin : explicit_titins) {
                titin_bound_actin_current.push_back(
                    titin.bound ? titin.bound_actin : -1);
                titin_last_bound_actin_current.push_back(
                    titin.last_bound_actin);
                titin_attach_step_current.push_back(
                    static_cast<int>(titin.attach_step));
                titin_last_detach_step_current.push_back(
                    static_cast<int>(titin.last_detach_step));
                titin_distance_current.push_back(titin.distance);
                titin_extension_current.push_back(titin.extension);
                titin_force_current.push_back(titin.force);
                titin_attachment_distance_current.push_back(
                    titin.attachment_distance);
                titin_attachment_energy_current.push_back(
                    titin.attachment_energy);
                titin_attachment_force_current.push_back(
                    titin.attachment_force);
            }
        }

        pad_int_vector(aa_prev_pairs, static_cast<size_t>(max_aa_pairs * 2), -1);
        pad_int_vector(aa_prev_status, static_cast<size_t>(max_aa_pairs), 0);
        pad_int_vector(aa_prev_lifetime, static_cast<size_t>(max_aa_pairs), 0);
        pad_int_vector(aa_current_pairs, static_cast<size_t>(max_aa_pairs * 2), -1);
        pad_int_vector(aa_current_status, static_cast<size_t>(max_aa_pairs), 0);
        pad_int_vector(aa_current_lifetime, static_cast<size_t>(max_aa_pairs), 0);
        pad_int_vector(aa_current_attach_step, static_cast<size_t>(max_aa_pairs), -1);
        pad_int_vector(aa_current_status2_memory_until,
                       static_cast<size_t>(max_aa_pairs), 0);
        pad_int_vector(am_prev_pairs, static_cast<size_t>(max_am_pairs * 2), -1);
        pad_int_vector(am_current_pairs, static_cast<size_t>(max_am_pairs * 2), -1);
        pad_int_vector(am_current_state, static_cast<size_t>(max_am_pairs), 0);
        am_current_target.resize(static_cast<size_t>(max_am_pairs), 0.0);
        pad_int_vector(am_head_state_current, total_head_count, 0);
        pad_int_vector(am_head_bound_actin_current, total_head_count, -1);
        am_head_target_current.resize(total_head_count, 0.0);
        am_head_binding_azimuth_current.resize(total_head_count, 0.0);
        am_head_stroke_progress_current.resize(total_head_count, 0.0);
        am_head_strain_energy_current.resize(total_head_count, 0.0);
        am_head_detachment_rate_current.resize(total_head_count, 0.0);
        pad_int_vector(aa_recovery_pairs, static_cast<size_t>(max_recovery_pairs * 2), -1);
        pad_int_vector(aa_recovery_until_values, static_cast<size_t>(max_recovery_pairs), 0);

        append_to_dataset_int(group_state, "aa_pairs_prev_count", {static_cast<int>(std::min(aa_prev_count, max_aa_pairs))}, {1, 1});
        append_to_dataset_int(group_state, "aa_pairs_current_count", {static_cast<int>(std::min(aa_current_count, max_aa_pairs))}, {1, 1});
        append_to_dataset_int(group_state, "am_pairs_prev_count", {static_cast<int>(std::min(am_prev_count, max_am_pairs))}, {1, 1});
        append_to_dataset_int(group_state, "am_pairs_current_count", {static_cast<int>(std::min(am_current_count, max_am_pairs))}, {1, 1});
        append_to_dataset_int(group_state, "aa_recovery_count", {static_cast<int>(std::min(aa_recovery_count, max_recovery_pairs))}, {1, 1});
        append_to_dataset_int(group_state, "aa_pairs_prev", aa_prev_pairs, {1, max_aa_pairs, 2});
        append_to_dataset_int(group_state, "aa_status_prev", aa_prev_status, {1, max_aa_pairs, 1});
        append_to_dataset_int(group_state, "aa_lifetime_prev", aa_prev_lifetime, {1, max_aa_pairs, 1});
        append_to_dataset_int(group_state, "aa_pairs_current", aa_current_pairs, {1, max_aa_pairs, 2});
        append_to_dataset_int(group_state, "aa_status_current", aa_current_status, {1, max_aa_pairs, 1});
        append_to_dataset_int(group_state, "aa_lifetime_current", aa_current_lifetime, {1, max_aa_pairs, 1});
        append_to_dataset_int(group_state, "aa_attach_step_current", aa_current_attach_step, {1, max_aa_pairs, 1});
        append_to_dataset_int(group_state, "aa_status2_memory_until_current",
                              aa_current_status2_memory_until,
                              {1, max_aa_pairs, 1});
        append_to_dataset_int(group_state, "am_pairs_prev", am_prev_pairs, {1, max_am_pairs, 2});
        append_to_dataset_int(group_state, "am_pairs_current", am_current_pairs, {1, max_am_pairs, 2});
        append_to_dataset_int(group_state, "am_state_current", am_current_state, {1, max_am_pairs});
        append_to_dataset(group_state, "am_target_axial_current", am_current_target, {1, max_am_pairs});
        append_to_dataset_int(group_state, "am_head_state_current", am_head_state_current,
                              {1, myosin_width, heads_per_myosin});
        append_to_dataset_int(group_state, "am_head_bound_actin_current",
                              am_head_bound_actin_current,
                              {1, myosin_width, heads_per_myosin});
        append_to_dataset(group_state, "am_head_target_axial_current", am_head_target_current,
                          {1, myosin_width, heads_per_myosin});
        append_to_dataset(group_state, "am_head_binding_azimuth_current",
                          am_head_binding_azimuth_current,
                          {1, myosin_width, heads_per_myosin});
        append_to_dataset(group_state, "am_head_stroke_progress_current",
                          am_head_stroke_progress_current,
                          {1, myosin_width, heads_per_myosin});
        append_to_dataset(group_state, "am_head_strain_energy_current",
                          am_head_strain_energy_current,
                          {1, myosin_width, heads_per_myosin});
        append_to_dataset(group_state, "am_head_detachment_rate_current",
                          am_head_detachment_rate_current,
                          {1, myosin_width, heads_per_myosin});
        append_to_dataset(group_state, "myosin_material_e1_current",
                          flatten_vec_values(myosin_material_e1),
                          {1, myosin_width, 3});
        append_to_dataset(group_state, "myosin_material_e2_current",
                          flatten_vec_values(myosin_material_e2),
                          {1, myosin_width, 3});
        append_to_dataset_int(group_state, "aa_recovery_pairs", aa_recovery_pairs,
                              {1, max_recovery_pairs, 2});
        append_to_dataset_int(group_state, "aa_recovery_until_values", aa_recovery_until_values,
                              {1, max_recovery_pairs, 1});

        std::vector<double> neighbor_last_actin_x;
        std::vector<double> neighbor_last_actin_y;
        std::vector<double> neighbor_last_actin_z;
        std::vector<double> neighbor_last_myosin_x;
        std::vector<double> neighbor_last_myosin_y;
        std::vector<double> neighbor_last_myosin_z;
        neighbor_list.get_last_species_positions(
            neighbor_last_actin_x, neighbor_last_actin_y, neighbor_last_actin_z,
            neighbor_last_myosin_x, neighbor_last_myosin_y, neighbor_last_myosin_z);
        append_to_dataset(group_state, "neighbor_last_actin_x", neighbor_last_actin_x, {1, actin_width});
        append_to_dataset(group_state, "neighbor_last_actin_y", neighbor_last_actin_y, {1, actin_width});
        append_to_dataset(group_state, "neighbor_last_actin_z", neighbor_last_actin_z, {1, actin_width});
        append_to_dataset(group_state, "neighbor_last_myosin_x", neighbor_last_myosin_x, {1, myosin_width});
        append_to_dataset(group_state, "neighbor_last_myosin_y", neighbor_last_myosin_y, {1, myosin_width});
        append_to_dataset(group_state, "neighbor_last_myosin_z", neighbor_last_myosin_z, {1, myosin_width});
        append_to_dataset(group_state, "am_q_current", am_thin_activation,
                          {1, actin_width});
        std::vector<int> turnover_steps(actin_last_turnover_step.begin(),
                                        actin_last_turnover_step.end());
        append_to_dataset_int(
            group_state, "actin_last_turnover_step_current", turnover_steps,
            {1, actin_width});
        if (titin_mode == TitinMode::explicit_dynamic) {
            append_to_dataset_int(group_state, "titin_bound_actin_current",
                                  titin_bound_actin_current,
                                  {1, titin_width});
            append_to_dataset_int(group_state,
                                  "titin_last_bound_actin_current",
                                  titin_last_bound_actin_current,
                                  {1, titin_width});
            append_to_dataset_int(group_state, "titin_attach_step_current",
                                  titin_attach_step_current,
                                  {1, titin_width});
            append_to_dataset_int(group_state,
                                  "titin_last_detach_step_current",
                                  titin_last_detach_step_current,
                                  {1, titin_width});
            append_to_dataset(group_state, "titin_distance_current",
                              titin_distance_current, {1, titin_width});
            append_to_dataset(group_state, "titin_extension_current",
                              titin_extension_current, {1, titin_width});
            append_to_dataset(group_state, "titin_force_current",
                              titin_force_current, {1, titin_width});
            append_to_dataset(group_state,
                              "titin_attachment_distance_current",
                              titin_attachment_distance_current,
                              {1, titin_width});
            append_to_dataset(group_state,
                              "titin_attachment_energy_current",
                              titin_attachment_energy_current,
                              {1, titin_width});
            append_to_dataset(group_state,
                              "titin_attachment_force_current",
                              titin_attachment_force_current,
                              {1, titin_width});
            append_to_dataset_int(
                group_state, "titin_association_counts",
                {static_cast<int>(explicit_titin_supported_associations),
                 static_cast<int>(explicit_titin_unsupported_associations),
                 static_cast<int>(explicit_titin_rebind_associations),
                 static_cast<int>(explicit_titin_unsupported_releases)},
                {1, 4});
        }

        // Save RNG states so resumed trajectories can be bitwise reproducible.
        if (rng != nullptr) {
            const size_t main_state_size = gsl_rng_size(rng);
            if (main_state_size > 0) {
                ensure_state_dataset_int("rng_main_state", static_cast<hsize_t>(main_state_size));
                const auto* state_ptr = static_cast<const unsigned char*>(gsl_rng_state(rng));
                std::vector<int> main_state(main_state_size, 0);
                for (size_t idx = 0; idx < main_state_size; ++idx) {
                    main_state[idx] = static_cast<int>(state_ptr[idx]);
                }
                append_to_dataset_int(group_state, "rng_main_state", main_state,
                                      {1, static_cast<hsize_t>(main_state_size)});
            }
        }
        if (!rng_engines.empty() && rng_engines[0] != nullptr) {
            const int thread_count = static_cast<int>(rng_engines.size());
            const size_t thread_state_size = gsl_rng_size(rng_engines[0]);
            const hsize_t flat_width = static_cast<hsize_t>(thread_count) *
                                       static_cast<hsize_t>(thread_state_size);
            if (thread_state_size > 0) {
                ensure_state_dataset_int("rng_thread_state", flat_width);
                ensure_state_dataset_int("rng_thread_count", 1);
                ensure_state_dataset_int("rng_thread_state_size", 1);

                std::vector<int> thread_state(flat_width, 0);
                for (int t = 0; t < thread_count; ++t) {
                    if (rng_engines[t] == nullptr) {
                        continue;
                    }
                    const size_t local_size = gsl_rng_size(rng_engines[t]);
                    if (local_size != thread_state_size) {
                        continue;
                    }
                    const auto* local_ptr =
                        static_cast<const unsigned char*>(gsl_rng_state(rng_engines[t]));
                    size_t base = static_cast<size_t>(t) * thread_state_size;
                    for (size_t idx = 0; idx < thread_state_size; ++idx) {
                        thread_state[base + idx] = static_cast<int>(local_ptr[idx]);
                    }
                }
                append_to_dataset_int(group_state, "rng_thread_state", thread_state, {1, flat_width});
                append_to_dataset_int(group_state, "rng_thread_count", {thread_count}, {1, 1});
                append_to_dataset_int(group_state, "rng_thread_state_size",
                                      {static_cast<int>(thread_state_size)}, {1, 1});
            }
        }
        maybe_test_sleep("SARCOMERE_TEST_SLEEP_BEFORE_STATE_COMMIT", "state_commit");
        append_to_dataset_int(group_state, "state_commit_count", {1}, {1, 1});
        file.flush(H5F_SCOPE_GLOBAL);
    }

    // Sparse, frame-ragged local ACTN2 state.  The offset table maps each
    // saved step to a contiguous slice of identity and coordinate records,
    // avoiding an O(N_actin^2) trajectory representation.
    if (aa_local_multilink) {
        H5::H5File file(filename, H5F_ACC_RDWR);
        H5::Group group_cb(file.openGroup("/catch_bond"));
        if (!file.nameExists("/catch_bond/local_link_identity")) {
            create_empty_dataset_int(
                file, "/catch_bond", "local_link_identity",
                {0, 7}, {H5S_UNLIMITED, 7}, {256, 7});
            create_empty_dataset(
                file, "/catch_bond", "local_link_coordinates",
                {0, 2}, {H5S_UNLIMITED, 2}, {256, 2});
            create_empty_dataset_int(
                file, "/catch_bond", "local_link_frame_offsets",
                {0, 3}, {H5S_UNLIMITED, 3}, {64, 3});
        }
        // One row per identity record. The primary arm loads are the
        // motor-allocated, mechanically transmitted loads used by KMC. Raw
        // spring-arm forces and the allocated motor budgets are appended for
        // auditing and must never be confused with the lifetime input.
        if (!file.nameExists("/catch_bond/local_link_mechanics")) {
            create_empty_dataset(
                file, "/catch_bond", "local_link_mechanics",
                {0, 16}, {H5S_UNLIMITED, 16}, {256, 16});
        }
        if (!file.nameExists("/catch_bond/local_link_mechanics_v2")) {
            create_empty_dataset(
                file, "/catch_bond", "local_link_mechanics_v2",
                {0, 23}, {H5S_UNLIMITED, 23}, {256, 23});
        }
        if (!file.nameExists(
                "/catch_bond/local_link_mechanics_v2_frame_offsets")) {
            create_empty_dataset_int(
                file, "/catch_bond", "local_link_mechanics_v2_frame_offsets",
                {0, 4}, {H5S_UNLIMITED, 4}, {64, 4});
        }
        if (!file.nameExists("/catch_bond/local_link_arm_state")) {
            create_empty_dataset(
                file, "/catch_bond", "local_link_arm_state",
                {0, 9}, {H5S_UNLIMITED, 9}, {256, 9});
        }
        if (aa_status_mode ==
                actn2_load_kinetics::StatusMode::motor_transmitted_register &&
            !file.nameExists("/catch_bond/local_link_motor_support")) {
            // Rows align one-for-one with local_link_identity.
            create_empty_dataset(
                file, "/catch_bond", "local_link_motor_support",
                {0, 20}, {H5S_UNLIMITED, 20}, {256, 20});
        }
        H5::DataSet identity_dataset = group_cb.openDataSet("local_link_identity");
        H5::DataSpace identity_space = identity_dataset.getSpace();
        hsize_t identity_dims[2] = {0, 0};
        identity_space.getSimpleExtentDims(identity_dims);
        H5::DataSet mechanics_v2_dataset =
            group_cb.openDataSet("local_link_mechanics_v2");
        H5::DataSpace mechanics_v2_space = mechanics_v2_dataset.getSpace();
        hsize_t mechanics_v2_dims[2] = {0, 0};
        mechanics_v2_space.getSimpleExtentDims(mechanics_v2_dims);
        const size_t count = aa_active_link_slots.size();
        std::vector<int> identity;
        std::vector<double> coordinates;
        std::vector<double> mechanics_rows;
        std::vector<double> mechanics_rows_v2;
        std::vector<double> motor_support_rows;
        std::vector<double> arm_state_rows;
        identity.reserve(count * 7);
        coordinates.reserve(count * 2);
        mechanics_rows.reserve(count * 16);
        mechanics_rows_v2.reserve(count * 23);
        motor_support_rows.reserve(count * 20);
        arm_state_rows.reserve(count * 9);
        const double ramp = _current_interaction_ramp();
        const double single_link_alignment_k = ramp *
            actn2_load_kinetics::fixed_single_link_stiffness(
                kappa_aa, aa_stiffness_reference_links);
        std::vector<LocalAALinkResult> saved_mechanics(aa_links.size());
        std::vector<actn2_load_kinetics::ArmLoads> raw_arm_loads(aa_links.size());
        std::vector<double> total_raw_favorable_arm_load(
            static_cast<size_t>(actin.n), 0.0);
        for (size_t slot : aa_active_link_slots) {
            const auto& link = aa_links[slot];
            if (aa_status_mode ==
                actn2_load_kinetics::StatusMode::motor_transmitted_register) {
                const double age = link.candidate_myosin_i >= 0 &&
                        current_step >=
                            static_cast<size_t>(link.candidate_since)
                    ? (current_step -
                       static_cast<size_t>(link.candidate_since)) * dt
                    : 0.0;
                const double challenger_age = link.challenger_myosin_i >= 0 &&
                        current_step >=
                            static_cast<size_t>(link.challenger_since)
                    ? (current_step -
                       static_cast<size_t>(link.challenger_since)) * dt
                    : 0.0;
                motor_support_rows.insert(motor_support_rows.end(), {
                    static_cast<double>(link.candidate_myosin_i),
                    static_cast<double>(link.candidate_myosin_j),
                    static_cast<double>(link.candidate_half_i),
                    static_cast<double>(link.candidate_half_j),
                    link.candidate_support_ema,
                    link.candidate_support_instantaneous, age,
                    link.candidate_locked ? 1.0 : 0.0,
                    static_cast<double>(link.candidate_since),
                    static_cast<double>(link.candidate_last_seen),
                    link.candidate_raw_support_i,
                    link.candidate_raw_support_j,
                    link.candidate_allocated_support_i,
                    link.candidate_allocated_support_j,
                    static_cast<double>(link.challenger_myosin_i),
                    static_cast<double>(link.challenger_myosin_j),
                    static_cast<double>(link.challenger_half_i),
                    static_cast<double>(link.challenger_half_j),
                    link.challenger_support_ema, challenger_age});
            }
            if (aa_arm_state_model &&
                actn2_arm_kinetics::is_tether(link.arm_state)) {
                continue;
            }
            const double link_radial_k =
                ramp * aa_link_radial_stiffness(link.status);
            saved_mechanics[slot] = compute_local_aa_link(
                actin, link.actin_i, link.actin_j, link.s_i, link.s_j, box,
                link_radial_k, ramp * aa_link_perpendicular_stiffness(link.status),
                single_link_alignment_k, aa_optimal, aa_antiparallel_overlap,
                aa_link_register_ratio(link.status),
                aa_register_saturation_error, aa_axial_register_force_max,
            ramp * aa_radial_force_max);
            const auto& mech = saved_mechanics[slot];
            const std::array<double, 3> link_force{
                mech.force_on_first.x, mech.force_on_first.y,
                mech.force_on_first.z};
            const std::array<double, 3> direction_i{
                actin.direction[link.actin_i].x,
                actin.direction[link.actin_i].y,
                actin.direction[link.actin_i].z};
            const std::array<double, 3> direction_j{
                actin.direction[link.actin_j].x,
                actin.direction[link.actin_j].y,
                actin.direction[link.actin_j].z};
            raw_arm_loads[slot] = actn2_load_kinetics::signed_arm_loads(
                link_force, direction_i, direction_j);
            total_raw_favorable_arm_load[static_cast<size_t>(link.actin_i)] +=
                std::max(0.0, raw_arm_loads[slot].first);
            total_raw_favorable_arm_load[static_cast<size_t>(link.actin_j)] +=
                std::max(0.0, raw_arm_loads[slot].second);
        }
        std::vector<double> productive_motor_load(
            static_cast<size_t>(actin.n), 0.0);
        for (int i = 0; i < actin.n; ++i) {
            productive_motor_load[static_cast<size_t>(i)] =
                actn2_load_kinetics::norm(
                    _actin_productive_motor_force_vector(i));
        }
        for (size_t slot : aa_active_link_slots) {
            const auto& link = aa_links[slot];
            identity.insert(identity.end(), {
                static_cast<int>(link.id >> 32),
                static_cast<int>(link.id & 0xffffffffULL),
                link.actin_i, link.actin_j, static_cast<int>(link.status),
                static_cast<int>(link.attach_step),
                static_cast<int>(link.status2_memory_until)});
            coordinates.insert(coordinates.end(), {link.s_i, link.s_j});
            arm_state_rows.insert(arm_state_rows.end(), {
                static_cast<double>(static_cast<int>(link.arm_state)),
                static_cast<double>(link.arm_state_since),
                static_cast<double>(link.last_partner),
                actn2_arm_kinetics::arm_i_bound(link.arm_state) ? 1.0 : 0.0,
                actn2_arm_kinetics::arm_j_bound(link.arm_state) ? 1.0 : 0.0,
                link.arm_force_i, link.arm_force_j, link.contact_force,
                link.contact_penetration});
            if (aa_arm_state_model &&
                actn2_arm_kinetics::is_tether(link.arm_state)) {
                mechanics_rows.insert(mechanics_rows.end(), 16, 0.0);
                mechanics_rows_v2.insert(mechanics_rows_v2.end(), 23, 0.0);
                continue;
            }
            const auto& mech = saved_mechanics[slot];
            const auto raw_arms = raw_arm_loads[slot];
            const bool motor_allocated_mode = aa_status_mode ==
                actn2_load_kinetics::StatusMode::arm_load_productive;
            const bool productive_register_mode = aa_status_mode ==
                actn2_load_kinetics::StatusMode::productive_register;
            const double allocated_motor_i = motor_allocated_mode
                ? actn2_load_kinetics::allocated_motor_load(
                    raw_arms.first,
                    total_raw_favorable_arm_load[
                        static_cast<size_t>(link.actin_i)],
                    productive_motor_load[static_cast<size_t>(link.actin_i)])
                : 0.0;
            const double allocated_motor_j = motor_allocated_mode
                ? actn2_load_kinetics::allocated_motor_load(
                    raw_arms.second,
                    total_raw_favorable_arm_load[
                        static_cast<size_t>(link.actin_j)],
                    productive_motor_load[static_cast<size_t>(link.actin_j)])
                : 0.0;
            const actn2_load_kinetics::ArmLoads arms = motor_allocated_mode
                ? actn2_load_kinetics::ArmLoads{
                      std::min(std::max(0.0, raw_arms.first), allocated_motor_i),
                      std::min(std::max(0.0, raw_arms.second), allocated_motor_j)}
                : raw_arms;
            const double productive_support = _aa_productive_opposed_load(
                link.actin_i, link.actin_j);
            const double bilateral_load = arms.favorable_bilateral();
            const double registry_threshold = link.status == 2
                ? aa_status2_registry_off : aa_status2_registry_on;
            const double arm_load_threshold = link.status == 2
                ? aa_status2_force_off : aa_status2_force_on;
            const bool registry_ok = actn2_load_kinetics::registry_gate(
                static_cast<int>(link.status), mech.direction_dot,
                mech.registry_error, aa_status2_registry_on,
                aa_status2_registry_off);
            const bool arm_load_ok =
                arms.first > arm_load_threshold &&
                arms.second > arm_load_threshold;
            const bool productive_support_ok = productive_support > 0.0;
            const bool status2_gate = productive_register_mode
                ? actn2_load_kinetics::classify_productive_register(
                      static_cast<int>(link.status), productive_support,
                      aa_status2_force_on, mech.direction_dot,
                      mech.registry_error, aa_status2_registry_on,
                      aa_status2_registry_off) == 2
                : motor_allocated_mode
                ? arm_load_ok && productive_support_ok && registry_ok
                : actn2_load_kinetics::classify_arm_loaded_productive(
                      static_cast<int>(link.status), arms, productive_support,
                      aa_status2_force_on, aa_status2_force_off) == 2;
            const double lifetime_load = link.status == 2
                ? (productive_register_mode
                       ? productive_support : bilateral_load)
                : 0.0;
            const double lifetime = catch_bond_kinetics::lifetime_for_status(
                link.status, base_lifetime, directional_base_lifetime,
                lifetime_coeff, lifetime_load, stabilize_status1);
            const double legacy_tension = std::max(
                0.0, ramp * aa_link_radial_stiffness(link.status) *
                    (mech.length - aa_optimal));
            mechanics_rows.insert(mechanics_rows.end(), {
                mech.connector.x, mech.connector.y, mech.connector.z,
                mech.length, mech.length - aa_optimal,
                arms.first, arms.second, bilateral_load,
                productive_support, status2_gate ? 1.0 : 0.0,
                lifetime, legacy_tension,
                raw_arms.first, raw_arms.second,
                allocated_motor_i, allocated_motor_j});
            mechanics_rows_v2.insert(mechanics_rows_v2.end(), {
                mech.connector.x, mech.connector.y, mech.connector.z,
                mech.length, mech.length - aa_optimal,
                arms.first, arms.second, bilateral_load,
                productive_support, status2_gate ? 1.0 : 0.0,
                lifetime, legacy_tension,
                raw_arms.first, raw_arms.second,
                allocated_motor_i, allocated_motor_j,
                mech.direction_dot, mech.registry_error, registry_threshold,
                registry_ok ? 1.0 : 0.0,
                arm_load_ok ? 1.0 : 0.0,
                productive_support_ok ? 1.0 : 0.0,
                status2_gate ? 1.0 : 0.0});
        }
        if (count > 0) {
            append_to_dataset_int(group_cb, "local_link_identity", identity,
                                  {static_cast<hsize_t>(count), 7});
            append_to_dataset(group_cb, "local_link_coordinates", coordinates,
                              {static_cast<hsize_t>(count), 2});
            append_to_dataset(group_cb, "local_link_mechanics", mechanics_rows,
                              {static_cast<hsize_t>(count), 16});
            append_to_dataset(group_cb, "local_link_mechanics_v2",
                              mechanics_rows_v2,
                              {static_cast<hsize_t>(count), 23});
            append_to_dataset(group_cb, "local_link_arm_state",
                              arm_state_rows,
                              {static_cast<hsize_t>(count), 9});
            if (aa_status_mode ==
                actn2_load_kinetics::StatusMode::motor_transmitted_register) {
                append_to_dataset(group_cb, "local_link_motor_support",
                                  motor_support_rows,
                                  {static_cast<hsize_t>(count), 20});
            }
        }
        append_to_dataset_int(
            group_cb, "local_link_frame_offsets",
            {static_cast<int>(current_step), static_cast<int>(identity_dims[0]),
             static_cast<int>(count)}, {1, 3});
        append_to_dataset_int(
            group_cb, "local_link_mechanics_v2_frame_offsets",
            {static_cast<int>(current_step), static_cast<int>(identity_dims[0]),
             static_cast<int>(mechanics_v2_dims[0]), static_cast<int>(count)},
            {1, 4});
        file.flush(H5F_SCOPE_GLOBAL);
    }

    // Sparse frame-ragged myomesin occupancy.  Empty frames still receive an
    // offset row, so occupancy is unambiguous and arbitrary-frame resume can
    // restore the exact molecular graph.
    if (k_mm > 0.0) {
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/myomesin")) file.createGroup("/myomesin");
        if (!file.nameExists("/myomesin/bridge_state")) {
            create_empty_dataset(
                file, "/myomesin", "bridge_state",
                {0, 8}, {H5S_UNLIMITED, 8}, {128, 8});
            create_empty_dataset_int(
                file, "/myomesin", "frame_offsets",
                {0, 3}, {H5S_UNLIMITED, 3}, {128, 3});
        }
        if (!file.nameExists("/myomesin/next_id_by_frame")) {
            create_empty_dataset_int(
                file, "/myomesin", "next_id_by_frame",
                {0, 2}, {H5S_UNLIMITED, 2}, {128, 2});
        }
        H5::Group group(file.openGroup("/myomesin"));
        H5::DataSet state_dataset = group.openDataSet("bridge_state");
        H5::DataSpace state_space = state_dataset.getSpace();
        hsize_t state_dims[2] = {0, 0};
        state_space.getSimpleExtentDims(state_dims, nullptr);
        std::vector<double> rows;
        rows.reserve(myomesin_bridges.size() * 8);
        for (const auto& bridge : myomesin_bridges) {
            rows.insert(rows.end(), {
                static_cast<double>(bridge.id >> 32),
                static_cast<double>(bridge.id & 0xffffffffULL),
                static_cast<double>(bridge.myosin_i),
                static_cast<double>(bridge.myosin_j),
                static_cast<double>(bridge.attach_step),
                bridge.axial_offset, bridge.radial_distance, bridge.force});
        }
        if (!myomesin_bridges.empty()) {
            append_to_dataset(
                group, "bridge_state", rows,
                {static_cast<hsize_t>(myomesin_bridges.size()), 8});
        }
        append_to_dataset_int(
            group, "frame_offsets",
            {static_cast<int>(current_step), static_cast<int>(state_dims[0]),
             static_cast<int>(myomesin_bridges.size())},
            {1, 3});
        append_to_dataset_int(
            group, "next_id_by_frame",
            {static_cast<int>(myomesin_next_id >> 32),
             static_cast<int>(myomesin_next_id & 0xffffffffULL)},
            {1, 2});
    }

    // Flush any recorded catch-bond events to the HDF5 file
    if (!cb_breakage_events.empty() || !cb_limit_events.empty() || !aa_completed_lifetimes.empty()) {
        H5::H5File file(filename, H5F_ACC_RDWR);
        H5::Group group_cb(file.openGroup("/catch_bond"));

        if (!cb_breakage_events.empty()) {
            hsize_t event_width = 11 + 2 * max_myosin_bonds;
            hsize_t n_events = cb_breakage_events.size() / event_width;
            append_to_dataset(group_cb, "breakage", cb_breakage_events,
                               {n_events, event_width});
            cb_breakage_events.clear();
        }

        if (!cb_limit_events.empty()) {
            hsize_t limit_width = 5;
            hsize_t n_limit = cb_limit_events.size() / limit_width;
            append_to_dataset(group_cb, "limit_removal", cb_limit_events,
                               {n_limit, limit_width});
            cb_limit_events.clear();
        }

        if (!aa_completed_lifetimes.empty()) {
            hsize_t n = aa_completed_lifetimes.size();
            append_to_dataset(group_cb, "completed_lifetimes", aa_completed_lifetimes, {n});
            aa_completed_lifetimes.clear();
        }
    }

    if (!am_transition_events.empty()) {
        constexpr hsize_t transition_width = 11;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/actin_myo/head_transitions")) {
            create_empty_dataset(
                file, "/actin_myo", "head_transitions",
                {0, transition_width}, {H5S_UNLIMITED, transition_width},
                {100, transition_width});
        }
        H5::Group group_am(file.openGroup("/actin_myo"));
        const hsize_t n_events = am_transition_events.size() / transition_width;
        append_to_dataset(group_am, "head_transitions", am_transition_events,
                          {n_events, transition_width});
        am_transition_events.clear();
    }

    if (!am_forced_rupture_events.empty()) {
        constexpr hsize_t rupture_width = 15;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/actin_myo/forced_ruptures")) {
            create_empty_dataset(
                file, "/actin_myo", "forced_ruptures",
                {0, rupture_width}, {H5S_UNLIMITED, rupture_width},
                {100, rupture_width});
        }
        H5::Group group_am(file.openGroup("/actin_myo"));
        const hsize_t n_events =
            am_forced_rupture_events.size() / rupture_width;
        append_to_dataset(group_am, "forced_ruptures",
                          am_forced_rupture_events,
                          {n_events, rupture_width});
        am_forced_rupture_events.clear();
    }
    if (!am_forced_rupture_events_v2.empty()) {
        constexpr hsize_t rupture_width = 19;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/actin_myo/forced_ruptures_v2")) {
            create_empty_dataset(
                file, "/actin_myo", "forced_ruptures_v2",
                {0, rupture_width}, {H5S_UNLIMITED, rupture_width},
                {100, rupture_width});
        }
        H5::Group group_am(file.openGroup("/actin_myo"));
        const hsize_t n_events =
            am_forced_rupture_events_v2.size() / rupture_width;
        append_to_dataset(group_am, "forced_ruptures_v2",
                          am_forced_rupture_events_v2,
                          {n_events, rupture_width});
        am_forced_rupture_events_v2.clear();
    }

    if (!am_thick_transition_events.empty()) {
        constexpr hsize_t transition_width = 7;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/actin_myo/thick_transitions")) {
            create_empty_dataset(
                file, "/actin_myo", "thick_transitions",
                {0, transition_width}, {H5S_UNLIMITED, transition_width},
                {100, transition_width});
        }
        H5::Group group_am(file.openGroup("/actin_myo"));
        const hsize_t n_events =
            am_thick_transition_events.size() / transition_width;
        append_to_dataset(group_am, "thick_transitions",
                          am_thick_transition_events,
                          {n_events, transition_width});
        am_thick_transition_events.clear();
    }

    if (!actin_turnover_events.empty()) {
        constexpr hsize_t event_width = 20;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/actin/turnover_events")) {
            create_empty_dataset(
                file, "/actin", "turnover_events",
                {0, event_width}, {H5S_UNLIMITED, event_width},
                {100, event_width});
        }
        H5::Group group_actin(file.openGroup("/actin"));
        const hsize_t n_events = actin_turnover_events.size() / event_width;
        append_to_dataset(group_actin, "turnover_events",
                          actin_turnover_events, {n_events, event_width});
        actin_turnover_events.clear();
    }

    if (!aa_state_events.empty()) {
        constexpr hsize_t state_width = 13;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/catch_bond/state_transitions")) {
            create_empty_dataset(
                file, "/catch_bond", "state_transitions",
                {0, state_width}, {H5S_UNLIMITED, state_width},
                {100, state_width});
        }
        H5::Group group_cb(file.openGroup("/catch_bond"));
        const hsize_t n_events = aa_state_events.size() / state_width;
        append_to_dataset(group_cb, "state_transitions", aa_state_events,
                          {n_events, state_width});
        aa_state_events.clear();
    }

    if (!aa_load_events.empty()) {
        constexpr hsize_t load_width = 11;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/catch_bond/load_transitions")) {
            create_empty_dataset(
                file, "/catch_bond", "load_transitions",
                {0, load_width}, {H5S_UNLIMITED, load_width},
                {100, load_width});
        }
        H5::Group group_cb(file.openGroup("/catch_bond"));
        const hsize_t n_events = aa_load_events.size() / load_width;
        append_to_dataset(group_cb, "load_transitions", aa_load_events,
                          {n_events, load_width});
        aa_load_events.clear();
    }

    if (!aa_local_link_events.empty()) {
        constexpr hsize_t event_width = 11;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/catch_bond/local_link_events")) {
            create_empty_dataset(
                file, "/catch_bond", "local_link_events",
                {0, event_width}, {H5S_UNLIMITED, event_width},
                {100, event_width});
        }
        H5::Group group_cb(file.openGroup("/catch_bond"));
        const hsize_t n_events = aa_local_link_events.size() / event_width;
        append_to_dataset(group_cb, "local_link_events",
                          aa_local_link_events, {n_events, event_width});
        aa_local_link_events.clear();
    }
    if (!aa_local_link_events_v2.empty()) {
        constexpr hsize_t event_width = 21;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/catch_bond/local_link_events_v2")) {
            create_empty_dataset(
                file, "/catch_bond", "local_link_events_v2",
                {0, event_width}, {H5S_UNLIMITED, event_width},
                {100, event_width});
        }
        H5::Group group_cb(file.openGroup("/catch_bond"));
        const hsize_t n_events = aa_local_link_events_v2.size() / event_width;
        append_to_dataset(group_cb, "local_link_events_v2",
                          aa_local_link_events_v2, {n_events, event_width});
        aa_local_link_events_v2.clear();
    }
    if (!aa_local_link_events_v3.empty()) {
        constexpr hsize_t event_width = 35;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/catch_bond/local_link_events_v3")) {
            create_empty_dataset(
                file, "/catch_bond", "local_link_events_v3",
                {0, event_width}, {H5S_UNLIMITED, event_width},
                {100, event_width});
        }
        H5::Group group_cb(file.openGroup("/catch_bond"));
        const hsize_t n_events = aa_local_link_events_v3.size() / event_width;
        append_to_dataset(group_cb, "local_link_events_v3",
                          aa_local_link_events_v3, {n_events, event_width});
        aa_local_link_events_v3.clear();
    }
    if (!aa_tip_accessibility_events.empty()) {
        constexpr hsize_t event_width = 13;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/catch_bond/tip_accessibility_events")) {
            create_empty_dataset(
                file, "/catch_bond", "tip_accessibility_events",
                {0, event_width}, {H5S_UNLIMITED, event_width},
                {100, event_width});
        }
        H5::Group group_cb(file.openGroup("/catch_bond"));
        const hsize_t n_events =
            aa_tip_accessibility_events.size() / event_width;
        append_to_dataset(group_cb, "tip_accessibility_events",
                          aa_tip_accessibility_events,
                          {n_events, event_width});
        aa_tip_accessibility_events.clear();
    }
    if (!aa_arm_events.empty()) {
        constexpr hsize_t event_width = 14;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/catch_bond/arm_events")) {
            create_empty_dataset(
                file, "/catch_bond", "arm_events",
                {0, event_width}, {H5S_UNLIMITED, event_width},
                {100, event_width});
        }
        H5::Group group_cb(file.openGroup("/catch_bond"));
        const hsize_t n_events = aa_arm_events.size() / event_width;
        append_to_dataset(group_cb, "arm_events", aa_arm_events,
                          {n_events, event_width});
        aa_arm_events.clear();
    }
    if (!aa_motor_pair_events.empty()) {
        constexpr hsize_t event_width = 20;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/catch_bond/motor_pair_events")) {
            create_empty_dataset(
                file, "/catch_bond", "motor_pair_events",
                {0, event_width}, {H5S_UNLIMITED, event_width},
                {100, event_width});
        }
        H5::Group group_cb(file.openGroup("/catch_bond"));
        const hsize_t n_events = aa_motor_pair_events.size() / event_width;
        append_to_dataset(group_cb, "motor_pair_events",
                          aa_motor_pair_events, {n_events, event_width});
        aa_motor_pair_events.clear();
    }
    if (!explicit_titin_events.empty()) {
        constexpr hsize_t event_width = 10;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/titin")) file.createGroup("/titin");
        if (!file.nameExists("/titin/events")) {
            create_empty_dataset(
                file, "/titin", "events",
                {0, event_width}, {H5S_UNLIMITED, event_width},
                {100, event_width});
        }
        H5::Group group_titin(file.openGroup("/titin"));
        const hsize_t n_events = explicit_titin_events.size() / event_width;
        append_to_dataset(group_titin, "events", explicit_titin_events,
                          {n_events, event_width});
        explicit_titin_events.clear();
    }
    if (!myomesin_events.empty()) {
        constexpr hsize_t event_width = 9;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/myomesin")) file.createGroup("/myomesin");
        if (!file.nameExists("/myomesin/events")) {
            create_empty_dataset(
                file, "/myomesin", "events",
                {0, event_width}, {H5S_UNLIMITED, event_width},
                {100, event_width});
        }
        H5::Group group(file.openGroup("/myomesin"));
        append_to_dataset(
            group, "events", myomesin_events,
            {static_cast<hsize_t>(myomesin_events.size() / event_width),
             event_width});
        myomesin_events.clear();
    }
}

void Sarcomere::save_resume_snapshot() {
    H5::H5File file(filename, H5F_ACC_RDWR);
    const std::string final_resume_path = "/resume";
    const std::string temp_resume_path = "/resume_new";
    const std::string backup_resume_path = "/resume_prev";
    delete_link_if_exists(file, temp_resume_path);
    H5::Group group_resume = file.createGroup(temp_resume_path);

    const hsize_t n_actins = static_cast<hsize_t>(actin.n);
    const hsize_t n_myosins = static_cast<hsize_t>(myosin.n);
    const hsize_t max_aa_pairs = static_cast<hsize_t>(actin.n) *
                                 static_cast<hsize_t>(std::max(10, 2 * max_myosin_bonds));
    const hsize_t max_am_pairs = static_cast<hsize_t>(myosin.n) *
                                 static_cast<hsize_t>(max_myosin_bonds);
    const hsize_t max_recovery_pairs = max_aa_pairs;

    auto replace_dataset_int = [&](const std::string& name,
                                   const std::vector<hsize_t>& dims,
                                   const std::vector<int>& data) {
        const std::string full_path = temp_resume_path + "/" + name;
        if (file.nameExists(full_path)) {
            H5Ldelete(file.getId(), full_path.c_str(), H5P_DEFAULT);
        }
        H5::DataSpace dataspace(dims.size(), dims.data());
        H5::IntType datatype(H5::PredType::STD_I32LE);
        H5::DataSet dataset = group_resume.createDataSet(name, datatype, dataspace);
        dataset.write(data.data(), H5::PredType::STD_I32LE);
    };

    auto replace_dataset_double = [&](const std::string& name,
                                      const std::vector<hsize_t>& dims,
                                      const std::vector<double>& data) {
        const std::string full_path = temp_resume_path + "/" + name;
        if (file.nameExists(full_path)) {
            H5Ldelete(file.getId(), full_path.c_str(), H5P_DEFAULT);
        }
        H5::DataSpace dataspace(dims.size(), dims.data());
        H5::FloatType datatype(H5::PredType::IEEE_F64LE);
        H5::DataSet dataset = group_resume.createDataSet(name, datatype, dataspace);
        dataset.write(data.data(), H5::PredType::IEEE_F64LE);
    };

    auto aa_prev = pack_aa_bond_state(
        previous_aa_pairs, actin_actin_bonds_prev,
        actin_actin_status_prev, actin_actin_lifetime_prev);
    auto aa_current = pack_aa_bond_state(
        active_aa_pairs, actin_actin_bonds,
        actin_actin_status, actin_actin_lifetime);
    auto am_prev = pack_am_bond_pairs(previous_am_pairs, am_bonds_prev);
    auto am_current = pack_am_bond_pairs(active_am_pairs, am_bonds);
    auto aa_recovery = pack_recovery_state(
        active_recovery_pairs, actin_recovery_until, current_step);

    if (aa_recovery.until.size() > static_cast<size_t>(max_recovery_pairs)) {
        // Degrade gracefully (truncate) instead of aborting; see pad_int_vector / save_state.
        static bool warned_resume_capacity = false;
        if (!warned_resume_capacity) {
            fprintf(stderr, "Warning: recovery pair list (%zu) exceeds compact resume capacity (%llu); "
                    "truncating resume snapshot.\n",
                    aa_recovery.until.size(), static_cast<unsigned long long>(max_recovery_pairs));
            warned_resume_capacity = true;
        }
    }

    std::vector<int> aa_prev_pairs = aa_prev.pairs;
    std::vector<int> aa_prev_status = aa_prev.status;
    std::vector<int> aa_prev_lifetime = aa_prev.lifetime;
    std::vector<int> aa_current_pairs = aa_current.pairs;
    std::vector<int> aa_current_status = aa_current.status;
    std::vector<int> aa_current_lifetime = aa_current.lifetime;
    std::vector<int> aa_current_attach_step;
    std::vector<int> aa_current_status2_memory_until;
    aa_current_attach_step.reserve(aa_current.status.size());
    aa_current_status2_memory_until.reserve(aa_current.status.size());
    for (size_t idx = 0; idx + 1 < aa_current.pairs.size(); idx += 2) {
        const int a = aa_current.pairs[idx];
        const int b = aa_current.pairs[idx + 1];
        aa_current_attach_step.push_back(aa_attach_step[a][b]);
        aa_current_status2_memory_until.push_back(
            static_cast<int>(aa_status2_memory_until[a][b]));
    }
    std::vector<int> am_prev_pairs = am_prev;
    std::vector<int> am_current_pairs = am_current;
    std::vector<int> am_current_state;
    std::vector<double> am_current_target;
    std::vector<int> am_head_state_current;
    std::vector<int> am_head_bound_actin_current;
    std::vector<double> am_head_target_current;
    std::vector<double> am_head_binding_azimuth_current;
    std::vector<double> am_head_stroke_progress_current;
    am_current_state.reserve(am_current.size() / 2);
    am_current_target.reserve(am_current.size() / 2);
    for (size_t idx = 0; idx + 1 < am_current.size(); idx += 2) {
        const int a = am_current[idx];
        const int m = am_current[idx + 1];
        am_current_state.push_back(dynamic_am ? static_cast<int>(am_state[a][m]) : 2);
        am_current_target.push_back(dynamic_am ? am_target_axial[a][m]
                                               : _am_axial_coordinate(a, m));
    }
    const size_t total_head_count = static_cast<size_t>(myosin.n) *
                                    static_cast<size_t>(2 * am_heads_per_half);
    am_head_state_current.reserve(total_head_count);
    am_head_bound_actin_current.reserve(total_head_count);
    am_head_target_current.reserve(total_head_count);
    am_head_binding_azimuth_current.reserve(total_head_count);
    am_head_stroke_progress_current.reserve(total_head_count);
    for (const auto& myosin_halves : am_heads_by_myosin) {
        for (const auto& half_heads : myosin_halves) {
            for (const auto& head : half_heads) {
                am_head_state_current.push_back(static_cast<int>(head.state));
                am_head_bound_actin_current.push_back(head.bound_actin);
                am_head_target_current.push_back(head.target_axial);
                am_head_binding_azimuth_current.push_back(
                    head.binding_azimuth);
                am_head_stroke_progress_current.push_back(
                    head.stroke_progress);
            }
        }
    }
    std::vector<int> aa_recovery_pairs = aa_recovery.pairs;
    std::vector<int> aa_recovery_until_values = aa_recovery.until;
    const hsize_t local_link_rows = static_cast<hsize_t>(
        std::max<size_t>(1, aa_active_link_slots.size()));
    std::vector<int> local_link_identity(
        static_cast<size_t>(local_link_rows) * 7, -1);
    std::vector<double> local_link_coordinates(
        static_cast<size_t>(local_link_rows) * 2, 0.0);
    std::vector<int> local_link_motor_identity(
        static_cast<size_t>(local_link_rows) * 12, -1);
    std::vector<double> local_link_motor_evidence(
        static_cast<size_t>(local_link_rows) * 7, 0.0);
    std::vector<int> local_link_arm_state(
        static_cast<size_t>(local_link_rows) * 3, -1);
    std::vector<double> local_link_arm_loads(
        static_cast<size_t>(local_link_rows) * 4, 0.0);
    for (size_t record = 0; record < aa_active_link_slots.size(); ++record) {
        const auto& link = aa_links[aa_active_link_slots[record]];
        local_link_identity[7 * record] = static_cast<int>(link.id >> 32);
        local_link_identity[7 * record + 1] =
            static_cast<int>(link.id & 0xffffffffULL);
        local_link_identity[7 * record + 2] = link.actin_i;
        local_link_identity[7 * record + 3] = link.actin_j;
        local_link_identity[7 * record + 4] = static_cast<int>(link.status);
        local_link_identity[7 * record + 5] = static_cast<int>(link.attach_step);
        local_link_identity[7 * record + 6] =
            static_cast<int>(link.status2_memory_until);
        local_link_coordinates[2 * record] = link.s_i;
        local_link_coordinates[2 * record + 1] = link.s_j;
        local_link_motor_identity[12 * record] = link.candidate_myosin_i;
        local_link_motor_identity[12 * record + 1] = link.candidate_myosin_j;
        local_link_motor_identity[12 * record + 2] = link.candidate_half_i;
        local_link_motor_identity[12 * record + 3] = link.candidate_half_j;
        local_link_motor_identity[12 * record + 4] =
            static_cast<int>(link.candidate_since);
        local_link_motor_identity[12 * record + 5] =
            static_cast<int>(link.candidate_last_seen);
        local_link_motor_identity[12 * record + 6] =
            link.candidate_locked ? 1 : 0;
        local_link_motor_identity[12 * record + 7] =
            link.challenger_myosin_i;
        local_link_motor_identity[12 * record + 8] =
            link.challenger_myosin_j;
        local_link_motor_identity[12 * record + 9] =
            link.challenger_half_i;
        local_link_motor_identity[12 * record + 10] =
            link.challenger_half_j;
        local_link_motor_identity[12 * record + 11] =
            static_cast<int>(link.challenger_since);
        local_link_motor_evidence[7 * record] = link.candidate_support_ema;
        local_link_motor_evidence[7 * record + 1] =
            link.candidate_support_instantaneous;
        local_link_motor_evidence[7 * record + 2] =
            link.candidate_raw_support_i;
        local_link_motor_evidence[7 * record + 3] =
            link.candidate_raw_support_j;
        local_link_motor_evidence[7 * record + 4] =
            link.candidate_allocated_support_i;
        local_link_motor_evidence[7 * record + 5] =
            link.candidate_allocated_support_j;
        local_link_motor_evidence[7 * record + 6] =
            link.challenger_support_ema;
        local_link_arm_state[3 * record] =
            static_cast<int>(link.arm_state);
        local_link_arm_state[3 * record + 1] =
            static_cast<int>(link.arm_state_since);
        local_link_arm_state[3 * record + 2] = link.last_partner;
        local_link_arm_loads[4 * record] = link.arm_force_i;
        local_link_arm_loads[4 * record + 1] = link.arm_force_j;
        local_link_arm_loads[4 * record + 2] = link.contact_force;
        local_link_arm_loads[4 * record + 3] = link.contact_penetration;
    }
    const hsize_t titin_rows = static_cast<hsize_t>(explicit_titins.size());
    std::vector<int> titin_state_int;
    std::vector<double> titin_state_double;
    if (titin_mode == TitinMode::explicit_dynamic) {
        titin_state_int.reserve(explicit_titins.size() * 5);
        titin_state_double.reserve(explicit_titins.size() * 6);
        for (const auto& titin : explicit_titins) {
            titin_state_int.insert(titin_state_int.end(), {
                titin.bound ? 1 : 0,
                titin.bound ? titin.bound_actin : -1,
                titin.last_bound_actin,
                static_cast<int>(titin.attach_step),
                static_cast<int>(titin.last_detach_step)});
            titin_state_double.insert(titin_state_double.end(), {
                titin.distance, titin.extension, titin.force,
                titin.attachment_distance, titin.attachment_energy,
                titin.attachment_force});
        }
    }
    const hsize_t myomesin_rows = static_cast<hsize_t>(
        std::max<size_t>(1, static_cast<size_t>(myosin.n) *
                               static_cast<size_t>(myomesin_sites_per_filament) /
                               2));
    std::vector<int> myomesin_state_int(
        static_cast<size_t>(myomesin_rows) * 5, -1);
    std::vector<double> myomesin_state_double(
        static_cast<size_t>(myomesin_rows) * 3, 0.0);
    for (size_t index = 0; index < myomesin_bridges.size(); ++index) {
        const auto& bridge = myomesin_bridges[index];
        myomesin_state_int[5 * index] =
            static_cast<int>(bridge.id >> 32);
        myomesin_state_int[5 * index + 1] =
            static_cast<int>(bridge.id & 0xffffffffULL);
        myomesin_state_int[5 * index + 2] = bridge.myosin_i;
        myomesin_state_int[5 * index + 3] = bridge.myosin_j;
        myomesin_state_int[5 * index + 4] =
            static_cast<int>(bridge.attach_step);
        myomesin_state_double[3 * index] = bridge.axial_offset;
        myomesin_state_double[3 * index + 1] = bridge.radial_distance;
        myomesin_state_double[3 * index + 2] = bridge.force;
    }

    pad_int_vector(aa_prev_pairs, static_cast<size_t>(max_aa_pairs * 2), -1);
    pad_int_vector(aa_prev_status, static_cast<size_t>(max_aa_pairs), 0);
    pad_int_vector(aa_prev_lifetime, static_cast<size_t>(max_aa_pairs), 0);
    pad_int_vector(aa_current_pairs, static_cast<size_t>(max_aa_pairs * 2), -1);
    pad_int_vector(aa_current_status, static_cast<size_t>(max_aa_pairs), 0);
    pad_int_vector(aa_current_lifetime, static_cast<size_t>(max_aa_pairs), 0);
    pad_int_vector(aa_current_attach_step, static_cast<size_t>(max_aa_pairs), -1);
    pad_int_vector(aa_current_status2_memory_until,
                   static_cast<size_t>(max_aa_pairs), 0);
    pad_int_vector(am_prev_pairs, static_cast<size_t>(max_am_pairs * 2), -1);
    pad_int_vector(am_current_pairs, static_cast<size_t>(max_am_pairs * 2), -1);
    pad_int_vector(am_current_state, static_cast<size_t>(max_am_pairs), 0);
    am_current_target.resize(static_cast<size_t>(max_am_pairs), 0.0);
    pad_int_vector(am_head_state_current, total_head_count, 0);
    pad_int_vector(am_head_bound_actin_current, total_head_count, -1);
    am_head_target_current.resize(total_head_count, 0.0);
    am_head_binding_azimuth_current.resize(total_head_count, 0.0);
    am_head_stroke_progress_current.resize(total_head_count, 0.0);
    pad_int_vector(aa_recovery_pairs, static_cast<size_t>(max_recovery_pairs * 2), -1);
    pad_int_vector(aa_recovery_until_values, static_cast<size_t>(max_recovery_pairs), 0);

    replace_dataset_int("current_step", {1}, {static_cast<int>(current_step)});
    replace_dataset_double("actin_center", {n_actins, 3}, flatten_3d_array(actin.center));
    replace_dataset_double("actin_direction", {n_actins, 3}, flatten_3d_array(actin.direction));
    replace_dataset_double("myosin_center", {n_myosins, 3}, flatten_3d_array(myosin.center));
    replace_dataset_double("myosin_direction", {n_myosins, 3}, flatten_3d_array(myosin.direction));
    replace_dataset_double("myosin_material_e1_current", {n_myosins, 3},
                           flatten_vec_values(myosin_material_e1));
    replace_dataset_double("myosin_material_e2_current", {n_myosins, 3},
                           flatten_vec_values(myosin_material_e2));
    replace_dataset_double("am_q_current", {n_actins}, am_thin_activation);
    replace_dataset_int(
        "actin_last_turnover_step_current", {n_actins},
        std::vector<int>(actin_last_turnover_step.begin(),
                         actin_last_turnover_step.end()));
    replace_dataset_int("aa_pairs_prev_count", {1}, {static_cast<int>(std::min(aa_prev.status.size(), static_cast<size_t>(max_aa_pairs)))});
    replace_dataset_int("aa_pairs_current_count", {1}, {static_cast<int>(std::min(aa_current.status.size(), static_cast<size_t>(max_aa_pairs)))});
    replace_dataset_int("am_pairs_prev_count", {1}, {static_cast<int>(std::min(am_prev.size() / 2, static_cast<size_t>(max_am_pairs)))});
    replace_dataset_int("am_pairs_current_count", {1}, {static_cast<int>(std::min(am_current.size() / 2, static_cast<size_t>(max_am_pairs)))});
    replace_dataset_int("aa_pairs_prev", {max_aa_pairs, 2}, aa_prev_pairs);
    replace_dataset_int("aa_status_prev", {max_aa_pairs, 1}, aa_prev_status);
    replace_dataset_int("aa_lifetime_prev", {max_aa_pairs, 1}, aa_prev_lifetime);
    replace_dataset_int("aa_pairs_current", {max_aa_pairs, 2}, aa_current_pairs);
    replace_dataset_int("aa_status_current", {max_aa_pairs, 1}, aa_current_status);
    replace_dataset_int("aa_lifetime_current", {max_aa_pairs, 1}, aa_current_lifetime);
    replace_dataset_int("aa_attach_step_current", {max_aa_pairs, 1}, aa_current_attach_step);
    replace_dataset_int("aa_status2_memory_until_current", {max_aa_pairs, 1},
                        aa_current_status2_memory_until);
    replace_dataset_int("am_pairs_prev", {max_am_pairs, 2}, am_prev_pairs);
    replace_dataset_int("am_pairs_current", {max_am_pairs, 2}, am_current_pairs);
    replace_dataset_int("am_state_current", {max_am_pairs}, am_current_state);
    replace_dataset_double("am_target_axial_current", {max_am_pairs}, am_current_target);
    const hsize_t heads_per_myosin = static_cast<hsize_t>(2 * am_heads_per_half);
    replace_dataset_int("am_head_state_current",
                        {n_myosins, heads_per_myosin},
                        am_head_state_current);
    replace_dataset_int("am_head_bound_actin_current",
                        {n_myosins, heads_per_myosin},
                        am_head_bound_actin_current);
    replace_dataset_double("am_head_target_axial_current",
                           {n_myosins, heads_per_myosin},
                           am_head_target_current);
    replace_dataset_double("am_head_binding_azimuth_current",
                           {n_myosins, heads_per_myosin},
                           am_head_binding_azimuth_current);
    replace_dataset_double("am_head_stroke_progress_current",
                           {n_myosins, heads_per_myosin},
                           am_head_stroke_progress_current);
    replace_dataset_int("aa_recovery_count", {1}, {static_cast<int>(std::min(aa_recovery.until.size(), static_cast<size_t>(max_recovery_pairs)))});
    replace_dataset_int("aa_recovery_pairs", {max_recovery_pairs, 2}, aa_recovery_pairs);
    replace_dataset_int("aa_recovery_until_values", {max_recovery_pairs, 1}, aa_recovery_until_values);
    replace_dataset_int("aa_local_link_count", {1},
                        {static_cast<int>(aa_active_link_slots.size())});
    replace_dataset_int("aa_local_link_identity", {local_link_rows, 7},
                        local_link_identity);
    replace_dataset_double("aa_local_link_coordinates", {local_link_rows, 2},
                           local_link_coordinates);
    replace_dataset_int("aa_local_link_motor_identity", {local_link_rows, 12},
                        local_link_motor_identity);
    replace_dataset_double("aa_local_link_motor_evidence", {local_link_rows, 7},
                           local_link_motor_evidence);
    replace_dataset_int("aa_local_link_arm_state", {local_link_rows, 3},
                        local_link_arm_state);
    replace_dataset_double("aa_local_link_arm_loads", {local_link_rows, 4},
                           local_link_arm_loads);
    replace_dataset_int("aa_local_next_id", {2},
                        {static_cast<int>(aa_next_link_id >> 32),
                         static_cast<int>(aa_next_link_id & 0xffffffffULL)});
    if (titin_mode == TitinMode::explicit_dynamic) {
        replace_dataset_int("titin_state_int", {titin_rows, 5},
                            titin_state_int);
        replace_dataset_double("titin_state_double", {titin_rows, 6},
                               titin_state_double);
        // Carried across resume so the per-frame provenance series stays a
        // monotone cumulative count instead of restarting at zero.
        replace_dataset_int(
            "titin_association_counts", {1, 4},
            {static_cast<int>(explicit_titin_supported_associations),
             static_cast<int>(explicit_titin_unsupported_associations),
             static_cast<int>(explicit_titin_rebind_associations),
             static_cast<int>(explicit_titin_unsupported_releases)});
    }
    replace_dataset_int(
        "myomesin_bridge_count", {1},
        {static_cast<int>(myomesin_bridges.size())});
    replace_dataset_int("myomesin_state_int", {myomesin_rows, 5},
                        myomesin_state_int);
    replace_dataset_double("myomesin_state_double", {myomesin_rows, 3},
                           myomesin_state_double);
    replace_dataset_int(
        "myomesin_next_id", {2},
        {static_cast<int>(myomesin_next_id >> 32),
         static_cast<int>(myomesin_next_id & 0xffffffffULL)});

    std::vector<double> neighbor_last_actin_x;
    std::vector<double> neighbor_last_actin_y;
    std::vector<double> neighbor_last_actin_z;
    std::vector<double> neighbor_last_myosin_x;
    std::vector<double> neighbor_last_myosin_y;
    std::vector<double> neighbor_last_myosin_z;
    neighbor_list.get_last_species_positions(
        neighbor_last_actin_x, neighbor_last_actin_y, neighbor_last_actin_z,
        neighbor_last_myosin_x, neighbor_last_myosin_y, neighbor_last_myosin_z);
    replace_dataset_double("neighbor_last_actin_x", {n_actins}, neighbor_last_actin_x);
    replace_dataset_double("neighbor_last_actin_y", {n_actins}, neighbor_last_actin_y);
    replace_dataset_double("neighbor_last_actin_z", {n_actins}, neighbor_last_actin_z);
    replace_dataset_double("neighbor_last_myosin_x", {n_myosins}, neighbor_last_myosin_x);
    replace_dataset_double("neighbor_last_myosin_y", {n_myosins}, neighbor_last_myosin_y);
    replace_dataset_double("neighbor_last_myosin_z", {n_myosins}, neighbor_last_myosin_z);
    replace_dataset_double("cb_breakage_pending",
                           {static_cast<hsize_t>(cb_breakage_events.size())},
                           cb_breakage_events);
    replace_dataset_double("cb_limit_pending",
                           {static_cast<hsize_t>(cb_limit_events.size())},
                           cb_limit_events);
    replace_dataset_double("aa_completed_lifetimes_pending",
                           {static_cast<hsize_t>(aa_completed_lifetimes.size())},
                           aa_completed_lifetimes);
    replace_dataset_double("am_transition_events_pending",
                           {static_cast<hsize_t>(am_transition_events.size())},
                           am_transition_events);
    replace_dataset_double("am_forced_rupture_events_pending",
                           {static_cast<hsize_t>(am_forced_rupture_events.size())},
                           am_forced_rupture_events);
    replace_dataset_double("am_forced_rupture_events_v2_pending",
                           {static_cast<hsize_t>(am_forced_rupture_events_v2.size())},
                           am_forced_rupture_events_v2);
    replace_dataset_double(
        "am_thick_transition_events_pending",
        {static_cast<hsize_t>(am_thick_transition_events.size())},
        am_thick_transition_events);
    replace_dataset_double(
        "actin_turnover_events_pending",
        {static_cast<hsize_t>(actin_turnover_events.size())},
        actin_turnover_events);
    replace_dataset_double("aa_state_events_pending",
                           {static_cast<hsize_t>(aa_state_events.size())},
                           aa_state_events);
    replace_dataset_double("aa_load_events_pending",
                           {static_cast<hsize_t>(aa_load_events.size())},
                           aa_load_events);
    replace_dataset_double("aa_local_link_events_pending",
                           {static_cast<hsize_t>(aa_local_link_events.size())},
                           aa_local_link_events);
    replace_dataset_double("aa_local_link_events_v2_pending",
                           {static_cast<hsize_t>(aa_local_link_events_v2.size())},
                           aa_local_link_events_v2);
    replace_dataset_double("aa_local_link_events_v3_pending",
                           {static_cast<hsize_t>(aa_local_link_events_v3.size())},
                           aa_local_link_events_v3);
    replace_dataset_double(
        "aa_tip_accessibility_events_pending",
        {static_cast<hsize_t>(aa_tip_accessibility_events.size())},
        aa_tip_accessibility_events);
    replace_dataset_double("aa_motor_pair_events_pending",
                           {static_cast<hsize_t>(aa_motor_pair_events.size())},
                           aa_motor_pair_events);
    replace_dataset_double("aa_arm_events_pending",
                           {static_cast<hsize_t>(aa_arm_events.size())},
                           aa_arm_events);
    replace_dataset_double("explicit_titin_events_pending",
                           {static_cast<hsize_t>(explicit_titin_events.size())},
                           explicit_titin_events);
    replace_dataset_double("myomesin_events_pending",
                           {static_cast<hsize_t>(myomesin_events.size())},
                           myomesin_events);

    if (rng != nullptr) {
        const size_t main_state_size = gsl_rng_size(rng);
        std::vector<int> main_state(main_state_size, 0);
        const auto* state_ptr = static_cast<const unsigned char*>(gsl_rng_state(rng));
        for (size_t idx = 0; idx < main_state_size; ++idx) {
            main_state[idx] = static_cast<int>(state_ptr[idx]);
        }
        replace_dataset_int("rng_main_state", {static_cast<hsize_t>(main_state_size)}, main_state);
    }
    if (!rng_engines.empty() && rng_engines[0] != nullptr) {
        const int thread_count = static_cast<int>(rng_engines.size());
        const size_t thread_state_size = gsl_rng_size(rng_engines[0]);
        const hsize_t flat_width = static_cast<hsize_t>(thread_count) *
                                   static_cast<hsize_t>(thread_state_size);
        std::vector<int> thread_state(flat_width, 0);
        for (int t = 0; t < thread_count; ++t) {
            if (rng_engines[t] == nullptr) {
                continue;
            }
            const auto* local_ptr =
                static_cast<const unsigned char*>(gsl_rng_state(rng_engines[t]));
            const size_t base = static_cast<size_t>(t) * thread_state_size;
            for (size_t idx = 0; idx < thread_state_size; ++idx) {
                thread_state[base + idx] = static_cast<int>(local_ptr[idx]);
            }
        }
        replace_dataset_int("rng_thread_state", {flat_width}, thread_state);
        replace_dataset_int("rng_thread_count", {1}, {thread_count});
        replace_dataset_int("rng_thread_state_size", {1}, {static_cast<int>(thread_state_size)});
    }
    file.flush(H5F_SCOPE_GLOBAL);
    maybe_test_sleep("SARCOMERE_TEST_SLEEP_BEFORE_RESUME_PROMOTE", "resume_promote");
    promote_temp_group(file, final_resume_path, temp_resume_path, backup_resume_path);
}

int Sarcomere::load_state(int& n_frames, int frame_index){
    {
        H5::H5File file(filename, H5F_ACC_RDONLY);
        H5::Group meta = file.openGroup("/meta");
        if (aa_zdisc_tip_contact) {
            constexpr const char* expected_scope =
                "all_nearby_myosin_tips_head_independent";
            if (!meta.attrExists("aa_tip_accessibility_scope")) {
                throw std::runtime_error(
                    "Cannot resume ACTN2 tip contact from a checkpoint "
                    "without head-independent accessibility metadata");
            }
            H5::Attribute scope_attr =
                meta.openAttribute("aa_tip_accessibility_scope");
            H5::StrType scope_type = scope_attr.getStrType();
            std::string stored_scope;
            scope_attr.read(scope_type, stored_scope);
            if (stored_scope != expected_scope) {
                throw std::runtime_error(
                    "Cannot resume ACTN2 tip contact with legacy "
                    "motor-scoped accessibility semantics");
            }
        }
        if (dynamic_am) {
            int stored_model_version = 0;
            if (!meta.attrExists("am_model_version")) {
                throw std::runtime_error(
                    "Cannot resume dynamic AM: trajectory predates AA/AM spacing schema v9");
            }
            meta.openAttribute("am_model_version").read(
                H5::PredType::NATIVE_INT, &stored_model_version);
            if (stored_model_version < 22 || !meta.attrExists("am_heads_per_half")) {
                throw std::runtime_error(
                    "Cannot resume dynamic AM checkpoint predating continuous "
                    "strain-dependent retention, crown mechanics, and the ATP-limited "
                    "strong-head off-rate, motor-allocated ACTN2 catch load, and "
                    "overlap-scaled neighbor head capacity, and fixed per-link "
                    "ACTN2 mechanics, finite-systole heartbeat gating, and "
                    "occupied-pair angular compliance and optional explicit "
                    "per-head transverse compliance with independently coupled "
                    "axial/transverse escape barriers, global actin-site "
                    "ownership, and nonproductive double-overlap binding v22");
            }
            int stored_heads = 1;
            meta.openAttribute("am_heads_per_half").read(
                H5::PredType::NATIVE_INT, &stored_heads);
            if (stored_heads != am_heads_per_half) {
                throw std::runtime_error(
                    "Cannot resume with a different am_heads_per_half than trajectory metadata");
            }
            auto require_int_attr = [&](const char* name, int expected) {
                if (!meta.attrExists(name)) {
                    throw std::runtime_error(
                        std::string("Cannot resume: missing AM metadata ") + name);
                }
                int stored = 0;
                meta.openAttribute(name).read(H5::PredType::NATIVE_INT, &stored);
                if (stored != expected) {
                    throw std::runtime_error(
                        std::string("Cannot resume with a different ") + name);
                }
            };
            auto require_double_attr = [&](const char* name, double expected) {
                if (!meta.attrExists(name)) {
                    throw std::runtime_error(
                        std::string("Cannot resume: missing AM metadata ") + name);
                }
                double stored = 0.0;
                meta.openAttribute(name).read(
                    H5::PredType::NATIVE_DOUBLE, &stored);
                if (std::abs(stored - expected) >
                    1e-12 * std::max({1.0, std::abs(stored), std::abs(expected)})) {
                    throw std::runtime_error(
                        std::string("Cannot resume with a different ") + name);
                }
            };
            require_int_attr("am_max_heads_per_actin_half",
                             am_max_heads_per_actin_half);
            require_int_attr("am_sites_per_actin", am_sites_per_actin);
            require_int_attr("am_max_heads_per_actin",
                             am_max_heads_per_actin);
            require_int_attr("am_thin_activation_enabled",
                             am_thin_activation_enabled ? 1 : 0);
            require_int_attr("am_spatial_heads", am_spatial_heads ? 1 : 0);
            require_int_attr("sterics_enabled", sterics_enabled ? 1 : 0);
            require_int_attr("am_tether_saturation_heads",
                             am_tether_saturation_heads);
            require_double_attr("am_head_transverse_k",
                                am_head_transverse_k);
            require_double_attr("am_alignment_k", am_alignment_k);
            require_double_attr("am_alignment_slack_rad",
                                am_alignment_slack);
            require_int_attr("am_strain_dependent_detachment",
                             am_strain_dependent_detachment ? 1 : 0);
            require_double_attr("am_strain_detachment_theta",
                                am_strain_detachment_theta);
            require_double_attr(
                "am_transverse_strain_detachment_theta",
                am_transverse_strain_detachment_theta);
            require_double_attr("am_strain_detachment_max_exponent",
                                am_strain_detachment_max_exponent);
            require_double_attr("am_strong_off_floor_fraction",
                                am_strong_off_floor_fraction);
            require_double_attr("am_cycle_period", am_cycle_period);
            require_double_attr("am_activation_rise", am_activation_rise);
            require_double_attr("am_activation_decay", am_activation_decay);
            require_double_attr("am_activation_floor", am_activation_floor);
            require_double_attr("am_activation_active_fraction",
                                am_activation_active_fraction);
            require_double_attr("actin_radius", actin_radius);
            require_double_attr("am_surface_gap", am_surface_gap);
            require_double_attr("am_steric_k", am_steric_k);
            require_double_attr("am_preferred_centerline", am_preferred_centerline);
            require_double_attr("aa_preferred_centerline", aa_optimal);
            require_double_attr("aa_spacing_ratio", aa_spacing_ratio);
            require_double_attr("aa_capture_margin", aa_capture_margin);
            require_int_attr("am_thin_feedback_enabled",
                             am_thin_feedback_enabled ? 1 : 0);
            require_double_attr("am_q_calcium_open_rate", am_q_calcium_open_rate);
            require_double_attr("am_q_calcium_half", am_q_calcium_half);
            require_double_attr("am_q_calcium_hill", am_q_calcium_hill);
            require_double_attr("am_q_strong_feedback_rate",
                                am_q_strong_feedback_rate);
            require_double_attr("am_q_strong_half", am_q_strong_half);
            require_double_attr("am_q_strong_hill", am_q_strong_hill);
            require_double_attr("am_q_close_rate", am_q_close_rate);
            if (meta.attrExists("am_q_initialized_at_calcium_equilibrium")) {
                require_int_attr(
                    "am_q_initialized_at_calcium_equilibrium",
                    am_q_initialized_at_calcium_equilibrium ? 1 : 0);
                require_double_attr("am_initial_q", am_initial_q);
            } else if (am_q_initialized_at_calcium_equilibrium) {
                throw std::runtime_error(
                    "Cannot resume calcium-equilibrium q initialization from "
                    "a trajectory without q-initialization metadata");
            }
            require_int_attr("am_thick_feedback_enabled",
                             am_thick_feedback_enabled ? 1 : 0);
            require_double_attr("am_thick_recruit_rate", am_thick_recruit_rate);
            require_double_attr("am_thick_fold_rate", am_thick_fold_rate);
            require_double_attr("am_thick_tension_scale", am_thick_tension_scale);
            require_double_attr("am_thick_initial_on_fraction",
                                am_thick_initial_on_fraction);
            require_double_attr("am_thick_force_clamp", am_thick_force_clamp);
        }
        auto require_promoted_int = [&](const char* name, int expected,
                                        int legacy_default) {
            if (!meta.attrExists(name)) {
                if (expected != legacy_default) {
                    throw std::runtime_error(
                        std::string("Cannot resume: trajectory lacks metadata for ") + name);
                }
                return;
            }
            int stored = 0;
            meta.openAttribute(name).read(H5::PredType::NATIVE_INT, &stored);
            if (stored != expected) {
                throw std::runtime_error(
                    std::string("Cannot resume with a different ") + name);
            }
        };
        auto require_promoted_double = [&](const char* name, double expected,
                                           double legacy_default) {
            if (!meta.attrExists(name)) {
                if (std::abs(expected - legacy_default) > 1e-12) {
                    throw std::runtime_error(
                        std::string("Cannot resume: trajectory lacks metadata for ") + name);
                }
                return;
            }
            double stored = 0.0;
            meta.openAttribute(name).read(
                H5::PredType::NATIVE_DOUBLE, &stored);
            if (std::abs(stored - expected) >
                1e-12 * std::max({1.0, std::abs(stored), std::abs(expected)})) {
                throw std::runtime_error(
                    std::string("Cannot resume with a different ") + name);
            }
        };
        require_promoted_int(
            "actin_turnover_enabled", actin_turnover_enabled ? 1 : 0, 0);
        require_promoted_double(
            "actin_turnover_rate", actin_turnover_rate, 0.0);
        require_promoted_double(
            "actin_turnover_start_time", actin_turnover_start_time, 0.0);
        require_promoted_int(
            "actin_turnover_max_placement_attempts",
            actin_turnover_max_placement_attempts, 1000);
        require_promoted_double(
            "actin_turnover_clearance_ratio",
            actin_turnover_clearance_ratio, 1.0);
        require_promoted_double(
            "aa_barbed_lateral_range", aa_barbed_lateral_range, 0.0);
        require_promoted_double(
            "aa_barbed_lateral_k", aa_barbed_lateral_k, 0.0);
        require_promoted_double(
            "aa_barbed_lateral_force_max",
            aa_barbed_lateral_force_max, 0.0);
        require_promoted_double(
            "aa_barbed_lateral_axial_window",
            aa_barbed_lateral_axial_window, 0.25);
        require_promoted_int(
            "aa_register_all_statuses", aa_register_all_statuses ? 1 : 0, 0);
        require_promoted_double(
            "aa_axial_register_force_max", aa_axial_register_force_max, 0.0);
        require_promoted_double("aa_barbed_zone", aa_barbed_zone, -1.0);
        require_promoted_int(
            "titin_mode", static_cast<int>(titin_mode),
            static_cast<int>(TitinMode::legacy_crosslinked));
        require_promoted_double("titin_k", titin_k, 0.0);
        require_promoted_double("k_mm", k_mm, 0.0);
        if (k_mm > 0.0) {
            require_promoted_int("myomesin_schema_version", 1, 0);
            require_promoted_int("myomesin_sites_per_filament",
                                 myomesin_sites_per_filament, 2);
            require_promoted_double("myomesin_k_on", myomesin_k_on, 5000.0);
            require_promoted_double("myomesin_k_off", myomesin_k_off, 1.0);
            require_promoted_double("myomesin_axial_capture",
                                    myomesin_axial_capture, 0.25);
            require_promoted_double("myomesin_radial_rest",
                                    myomesin_radial_rest, 0.11);
            require_promoted_double("myomesin_radial_capture",
                                    myomesin_radial_capture, 0.04);
            require_promoted_double("myomesin_radial_k",
                                    myomesin_radial_k, 0.0);
            require_promoted_double("myomesin_alignment_cos",
                                    myomesin_alignment_cos, 0.94);
            require_promoted_double("myomesin_bell_distance",
                                    myomesin_bell_distance, 0.0005);
            require_promoted_double("myomesin_force_max",
                                    myomesin_force_max, 10.0);
        }
        if (titin_k > 0.0 || meta.attrExists("titin_rest_length")) {
            require_promoted_double(
                "titin_rest_length", titin_rest_length, 0.0);
        }
        if (titin_mode == TitinMode::explicit_dynamic) {
            require_promoted_int("titin_schema_version", 4, 0);
            require_promoted_double(
                "titin_k_off_unsupported",
                explicit_titin_k_off_unsupported, 0.0);
            require_promoted_int(
                "titin_rebind_requires_actn2",
                explicit_titin_rebind_requires_actn2 ? 1 : 0, 0);
            require_promoted_int(
                "titin_require_live_actn2",
                explicit_titin_require_live_actn2 ? 1 : 0, 0);
            require_promoted_int(
                "titin_links_per_half", explicit_titin_per_half, 3);
            require_promoted_double(
                "titin_capture_distance",
                explicit_titin_capture_distance, 0.70);
            require_promoted_double("titin_k_on", explicit_titin_k_on, 5000.0);
            require_promoted_double("titin_k_off", explicit_titin_k_off, 0.5);
            require_promoted_int(
                "titin_strain_gated_kinetics",
                explicit_titin_strain_gated_kinetics ? 1 : 0, 1);
            require_promoted_double(
                "titin_transition_theta",
                explicit_titin_transition_theta, 0.5);
            require_promoted_double(
                "titin_rebind_memory",
                explicit_titin_rebind_memory, 0.0);
            require_promoted_double(
                "titin_rebind_rate_multiplier",
                explicit_titin_rebind_rate_multiplier, 1.0);
        }
        if (aa_local_multilink || meta.attrExists("aa_local_multilink")) {
            auto require_local_int = [&](const char* name, int expected) {
                if (!meta.attrExists(name)) {
                    throw std::runtime_error(
                        std::string("Cannot resume: missing local ACTN2 metadata ") + name);
                }
                int stored = 0;
                meta.openAttribute(name).read(H5::PredType::NATIVE_INT, &stored);
                if (stored != expected) {
                    throw std::runtime_error(
                        std::string("Cannot resume with a different ") + name);
                }
            };
            auto require_local_double = [&](const char* name, double expected) {
                if (!meta.attrExists(name)) {
                    throw std::runtime_error(
                        std::string("Cannot resume: missing local ACTN2 metadata ") + name);
                }
                double stored = 0.0;
                meta.openAttribute(name).read(
                    H5::PredType::NATIVE_DOUBLE, &stored);
                if (std::abs(stored - expected) >
                    1e-12 * std::max({1.0, std::abs(stored), std::abs(expected)})) {
                    throw std::runtime_error(
                        std::string("Cannot resume with a different ") + name);
                }
            };
            require_local_int("aa_local_multilink", aa_local_multilink ? 1 : 0);
            if (aa_local_multilink) {
                require_local_int("aa_local_link_schema",
                                  aa_arm_state_model ? 5 : 4);
                require_local_int("aa_prebundle_links_per_pair",
                                  aa_prebundle_links_per_pair);
                require_local_double("aa_site_exclusion", aa_site_exclusion);
                require_local_double("aa_overlap_capture", aa_overlap_capture);
                require_local_double(
                    "aa_axial_register_ratio", aa_axial_register_ratio);
                require_local_double(
                    "aa_status2_registry_on", aa_status2_registry_on);
                require_local_double(
                    "aa_status2_registry_off", aa_status2_registry_off);
                require_local_double(
                    "aa_register_saturation_error",
                    aa_register_saturation_error);
                require_local_double("aa_perpendicular_k", aa_perpendicular_k);
                if (meta.attrExists("aa_status1_radial_k")) {
                    require_local_double(
                        "aa_status1_radial_k", aa_status1_radial_k);
                } else if (aa_status1_radial_k > 0.0) {
                    throw std::runtime_error(
                        "Cannot resume an older trajectory with a non-legacy "
                        "aa_status1_radial_k");
                }
                require_promoted_double(
                    "aa_status1_site_diffusion",
                    aa_status1_site_diffusion, 0.0);
                require_promoted_int(
                    "aa_status2_locked_until_detach",
                    aa_status2_locked_until_detach ? 1 : 0, 0);
                require_promoted_double(
                    "aa_status2_unlock_rate", aa_status2_unlock_rate, 0.0);
                require_promoted_double(
                    "aa_status2_unlock_force_scale",
                    aa_status2_unlock_force_scale, 0.1);
                require_promoted_double(
                    "aa_status2_unlock_registry_scale",
                    aa_status2_unlock_registry_scale, 0.02);
                require_promoted_int(
                    "aa_arm_state_model", aa_arm_state_model ? 1 : 0, 0);
                require_promoted_double(
                    "aa_arm_open_rate_weak", aa_arm_open_rate_weak, 20.0);
                require_promoted_double(
                    "aa_arm_open_rate_mature", aa_arm_open_rate_mature, 1.0);
                require_promoted_double(
                    "aa_arm_bell_distance", aa_arm_bell_distance, 0.001);
                require_promoted_double(
                    "aa_tether_off_rate", aa_tether_off_rate, 5.0);
                require_promoted_double(
                    "aa_tether_rebind_rate", aa_tether_rebind_rate, 100.0);
                require_promoted_double(
                    "aa_tether_capture_distance",
                    aa_tether_capture_distance, 0.10);
                require_promoted_double(
                    "aa_arm_rebind_theta", aa_arm_rebind_theta, 0.5);
                require_promoted_double(
                    "aa_actn2_contact_radius", aa_actn2_contact_radius, 0.005);
                require_promoted_double(
                    "aa_myosin_contact_epsilon",
                    aa_myosin_contact_epsilon, 0.00428);
                require_promoted_double(
                    "aa_myosin_contact_force_cap",
                    aa_myosin_contact_force_cap, 40.0);
                require_promoted_double(
                    "aa_myosin_contact_shell_width",
                    aa_myosin_contact_shell_width, 0.0);
                require_promoted_double(
                    "aa_myosin_contact_shell_energy",
                    aa_myosin_contact_shell_energy, 0.0);
                require_promoted_int(
                    "aa_zdisc_tip_contact", aa_zdisc_tip_contact ? 1 : 0, 0);
                require_promoted_double(
                    "aa_zdisc_tip_onset", aa_zdisc_tip_onset, 0.050);
                require_promoted_double(
                    "aa_zdisc_tip_stiffness", aa_zdisc_tip_stiffness, 1000.0);
                require_promoted_double(
                    "aa_zdisc_tip_force_max", aa_zdisc_tip_force_max, 20.0);
                require_promoted_double(
                    "aa_zdisc_tip_transverse_flat",
                    aa_zdisc_tip_transverse_flat, 0.100);
                require_promoted_double(
                    "aa_zdisc_tip_transverse_cutoff",
                    aa_zdisc_tip_transverse_cutoff, 0.130);
                require_promoted_double(
                    "aa_occlusion_rupture_penetration",
                    aa_occlusion_rupture_penetration, 0.010);
                require_promoted_double(
                    "aa_association_max_direction_dot",
                    aa_association_max_direction_dot, 1.0);
                require_promoted_double(
                    "aa_link_rupture_length", aa_link_rupture_length, 0.0);
                require_local_int("aa_stiffness_reference_links",
                                  aa_stiffness_reference_links);
                require_local_int("aa_total_link_pool", aa_total_link_pool);
            }
        }
        if (meta.attrExists("aa_status_mode")) {
            int stored_mode = 0;
            meta.openAttribute("aa_status_mode").read(
                H5::PredType::NATIVE_INT, &stored_mode);
            if (stored_mode != static_cast<int>(aa_status_mode)) {
                throw std::runtime_error(
                    "Cannot resume with a different aa_status_mode than the trajectory metadata");
            }
            if (meta.attrExists("aa_status2_force_on") &&
                meta.attrExists("aa_status2_force_off")) {
                double stored_on = 0.0;
                double stored_off = 0.0;
                meta.openAttribute("aa_status2_force_on").read(
                    H5::PredType::NATIVE_DOUBLE, &stored_on);
                meta.openAttribute("aa_status2_force_off").read(
                    H5::PredType::NATIVE_DOUBLE, &stored_off);
                if (std::abs(stored_on - aa_status2_force_on) > 1e-12 ||
                    std::abs(stored_off - aa_status2_force_off) > 1e-12) {
                    throw std::runtime_error(
                        "Cannot resume with different AA status force thresholds than the trajectory metadata");
                }
            } else if (aa_status_mode !=
                       actn2_load_kinetics::StatusMode::legacy_motor_identity) {
                throw std::runtime_error(
                    "Cannot resume load-classified AA kinetics without stored force thresholds");
            }
        } else if (aa_status_mode !=
                   actn2_load_kinetics::StatusMode::legacy_motor_identity) {
            throw std::runtime_error(
                "Cannot resume load-classified AA kinetics from a trajectory without AA status-mode metadata");
        }
    }
    int target_frame = load_from_file(filename, actin, myosin, actin_actin_bonds, n_frames, frame_index);

    const size_t aa_stride = static_cast<size_t>(actin.n) * static_cast<size_t>(actin.n);
    const size_t am_stride = static_cast<size_t>(actin.n) * static_cast<size_t>(myosin.n);

    auto load_state_vector = [&](H5::Group& group, const std::string& dataset_name,
                                 size_t expected_width, std::vector<int>& out) -> bool {
        try {
            if (!group.nameExists(dataset_name)) {
                return false;
            }
            std::vector<hsize_t> dims;
            std::vector<double> raw = load_from_dataset(group, dataset_name, dims);
            if (dims.size() < 2 || target_frame < 0 || target_frame >= static_cast<int>(dims[0])) {
                return false;
            }
            if (static_cast<size_t>(dims[1]) != expected_width) {
                return false;
            }
            size_t offset = static_cast<size_t>(target_frame) * expected_width;
            out.resize(expected_width);
            for (size_t i = 0; i < expected_width; ++i) {
                out[i] = static_cast<int>(std::llround(raw[offset + i]));
            }
            return true;
        } catch (H5::Exception&) {
            return false;
        }
    };

    auto load_state_vector_double = [&](H5::Group& group, const std::string& dataset_name,
                                        size_t expected_width, std::vector<double>& out) -> bool {
        try {
            if (!group.nameExists(dataset_name)) {
                return false;
            }
            std::vector<hsize_t> dims;
            std::vector<double> raw = load_from_dataset(group, dataset_name, dims);
            if (dims.size() < 2 || target_frame < 0 || target_frame >= static_cast<int>(dims[0])) {
                return false;
            }
            if (static_cast<size_t>(dims[1]) != expected_width) {
                return false;
            }
            size_t offset = static_cast<size_t>(target_frame) * expected_width;
            out.resize(expected_width);
            std::copy_n(raw.begin() + static_cast<std::ptrdiff_t>(offset),
                        static_cast<std::ptrdiff_t>(expected_width),
                        out.begin());
            return true;
        } catch (H5::Exception&) {
            return false;
        }
    };

    auto assign_aa_matrix = [&](const std::vector<int>& flat, auto& matrix) {
        if (flat.size() != aa_stride) {
            return;
        }
        for (int i = 0; i < actin.n; ++i) {
            for (int j = 0; j < actin.n; ++j) {
                matrix[i][j] = flat[static_cast<size_t>(i) * actin.n + j];
            }
        }
    };

    auto assign_am_matrix = [&](const std::vector<int>& flat, auto& matrix) {
        if (flat.size() != am_stride) {
            return;
        }
        for (int i = 0; i < actin.n; ++i) {
            for (int j = 0; j < myosin.n; ++j) {
                matrix[i][j] = flat[static_cast<size_t>(i) * myosin.n + j];
            }
        }
    };

    auto load_state_tensor_int = [&](H5::Group& group, const std::string& dataset_name,
                                     size_t expected_rows, size_t expected_cols,
                                     std::vector<int>& out) -> bool {
        try {
            if (!group.nameExists(dataset_name)) {
                return false;
            }
            std::vector<hsize_t> dims;
            std::vector<double> raw = load_from_dataset(group, dataset_name, dims);
            if (dims.size() < 3 || target_frame < 0 || target_frame >= static_cast<int>(dims[0])) {
                return false;
            }
            if (static_cast<size_t>(dims[1]) != expected_rows ||
                static_cast<size_t>(dims[2]) != expected_cols) {
                return false;
            }
            const size_t frame_width = expected_rows * expected_cols;
            const size_t offset = static_cast<size_t>(target_frame) * frame_width;
            out.resize(frame_width);
            for (size_t i = 0; i < frame_width; ++i) {
                out[i] = static_cast<int>(std::llround(raw[offset + i]));
            }
            return true;
        } catch (H5::Exception&) {
            return false;
        }
    };

    auto load_state_tensor_double = [&](H5::Group& group, const std::string& dataset_name,
                                        size_t expected_rows, size_t expected_cols,
                                        std::vector<double>& out) -> bool {
        try {
            if (!group.nameExists(dataset_name)) {
                return false;
            }
            std::vector<hsize_t> dims;
            std::vector<double> raw = load_from_dataset(group, dataset_name, dims);
            if (dims.size() < 3 || target_frame < 0 ||
                target_frame >= static_cast<int>(dims[0])) {
                return false;
            }
            if (static_cast<size_t>(dims[1]) != expected_rows ||
                static_cast<size_t>(dims[2]) != expected_cols) {
                return false;
            }
            const size_t frame_width = expected_rows * expected_cols;
            const size_t offset = static_cast<size_t>(target_frame) * frame_width;
            out.resize(frame_width);
            std::copy_n(raw.begin() + static_cast<std::ptrdiff_t>(offset),
                        static_cast<std::ptrdiff_t>(frame_width), out.begin());
            return true;
        } catch (H5::Exception&) {
            return false;
        }
    };

    auto load_fixed_vector_int = [&](H5::Group& group, const std::string& dataset_name,
                                     size_t expected_width, std::vector<int>& out) -> bool {
        try {
            if (!group.nameExists(dataset_name)) {
                return false;
            }
            std::vector<hsize_t> dims;
            std::vector<double> raw = load_from_dataset(group, dataset_name, dims);
            if (dims.size() != 1 || static_cast<size_t>(dims[0]) != expected_width) {
                return false;
            }
            out.resize(expected_width);
            for (size_t i = 0; i < expected_width; ++i) {
                out[i] = static_cast<int>(std::llround(raw[i]));
            }
            return true;
        } catch (H5::Exception&) {
            return false;
        }
    };

    auto load_fixed_vector_double = [&](H5::Group& group, const std::string& dataset_name,
                                        size_t expected_width, std::vector<double>& out) -> bool {
        try {
            if (!group.nameExists(dataset_name)) {
                return false;
            }
            std::vector<hsize_t> dims;
            out = load_from_dataset(group, dataset_name, dims);
            return dims.size() == 1 && static_cast<size_t>(dims[0]) == expected_width;
        } catch (H5::Exception&) {
            return false;
        }
    };

    auto load_any_vector_double = [&](H5::Group& group, const std::string& dataset_name,
                                      std::vector<double>& out) -> bool {
        try {
            if (!group.nameExists(dataset_name)) {
                return false;
            }
            std::vector<hsize_t> dims;
            out = load_from_dataset(group, dataset_name, dims);
            return dims.size() == 1;
        } catch (H5::Exception&) {
            return false;
        }
    };

    auto load_fixed_matrix_int = [&](H5::Group& group, const std::string& dataset_name,
                                     size_t expected_rows, size_t expected_cols,
                                     std::vector<int>& out) -> bool {
        try {
            if (!group.nameExists(dataset_name)) {
                return false;
            }
            std::vector<hsize_t> dims;
            std::vector<double> raw = load_from_dataset(group, dataset_name, dims);
            if (dims.size() != 2 || static_cast<size_t>(dims[0]) != expected_rows ||
                static_cast<size_t>(dims[1]) != expected_cols) {
                return false;
            }
            out.resize(expected_rows * expected_cols);
            for (size_t i = 0; i < out.size(); ++i) {
                out[i] = static_cast<int>(std::llround(raw[i]));
            }
            return true;
        } catch (H5::Exception&) {
            return false;
        }
    };

    auto load_fixed_matrix_double = [&](H5::Group& group, const std::string& dataset_name,
                                        size_t expected_rows, size_t expected_cols,
                                        std::vector<double>& out) -> bool {
        try {
            if (!group.nameExists(dataset_name)) {
                return false;
            }
            std::vector<hsize_t> dims;
            out = load_from_dataset(group, dataset_name, dims);
            return dims.size() == 2 && static_cast<size_t>(dims[0]) == expected_rows &&
                   static_cast<size_t>(dims[1]) == expected_cols;
        } catch (H5::Exception&) {
            return false;
        }
    };

    auto clear_aa_state = [&](auto& bonds, auto& status, auto& lifetime) {
        for (int i = 0; i < actin.n; ++i) {
            std::fill(bonds[i].begin(), bonds[i].end(), 0);
            std::fill(status[i].begin(), status[i].end(), 0);
            std::fill(lifetime[i].begin(), lifetime[i].end(), 0);
        }
    };

    auto clear_am_state = [&](auto& bonds) {
        for (int i = 0; i < actin.n; ++i) {
            std::fill(bonds[i].begin(), bonds[i].end(), 0);
        }
    };

    auto clear_recovery_state = [&]() {
        for (int i = 0; i < actin.n; ++i) {
            std::fill(actin_recovery_until[i].begin(), actin_recovery_until[i].end(), 0);
        }
    };

    auto rebuild_explicit_titin_ownership = [&]() {
        explicit_titin_owner_by_actin.assign(static_cast<size_t>(actin.n), -1);
        for (size_t index = 0; index < explicit_titins.size(); ++index) {
            auto& titin = explicit_titins[index];
            if (!titin.bound) {
                titin.bound_actin = -1;
                continue;
            }
            if (titin.bound_actin < 0 || titin.bound_actin >= actin.n) {
                throw std::runtime_error(
                    "Explicit titin checkpoint contains an invalid bound-actin index");
            }
            int& owner = explicit_titin_owner_by_actin[
                static_cast<size_t>(titin.bound_actin)];
            if (owner >= 0) {
                throw std::runtime_error(
                    "Explicit titin checkpoint assigns one actin to multiple titins");
            }
            owner = static_cast<int>(index);
        }
    };

    auto restore_aa_state_from_pairs = [&](H5::Group& group,
                                           const std::string& count_name,
                                           const std::string& pair_name,
                                           const std::string& status_name,
                                           const std::string& lifetime_name,
                                           auto& bonds,
                                           auto& status,
                                           auto& lifetime) -> bool {
        std::vector<int> counts;
        std::vector<int> pairs;
        std::vector<int> flat_status;
        std::vector<int> flat_lifetime;
        const size_t max_aa_pairs = static_cast<size_t>(actin.n) *
                                    static_cast<size_t>(std::max(10, 2 * max_myosin_bonds));
        if (!load_state_vector(group, count_name, 1, counts) ||
            !load_state_tensor_int(group, pair_name, max_aa_pairs, 2, pairs) ||
            !load_state_tensor_int(group, status_name, max_aa_pairs, 1, flat_status) ||
            !load_state_tensor_int(group, lifetime_name, max_aa_pairs, 1, flat_lifetime)) {
            return false;
        }
        const size_t count = static_cast<size_t>(std::max(0, counts[0]));
        clear_aa_state(bonds, status, lifetime);
        for (size_t idx = 0; idx < count; ++idx) {
            const int a = pairs[2 * idx];
            const int b = pairs[2 * idx + 1];
            if (a < 0 || a >= actin.n || b < 0 || b >= actin.n || a == b) {
                continue;
            }
            bonds[a][b] = 1;
            bonds[b][a] = 1;
            status[a][b] = flat_status[idx];
            status[b][a] = flat_status[idx];
            lifetime[a][b] = flat_lifetime[idx];
            lifetime[b][a] = flat_lifetime[idx];
        }
        return true;
    };

    auto restore_am_state_from_pairs = [&](H5::Group& group,
                                           const std::string& count_name,
                                           const std::string& pair_name,
                                           auto& bonds) -> bool {
        std::vector<int> counts;
        std::vector<int> pairs;
        const size_t max_am_pairs = static_cast<size_t>(myosin.n) *
                                    static_cast<size_t>(max_myosin_bonds);
        if (!load_state_vector(group, count_name, 1, counts) ||
            !load_state_tensor_int(group, pair_name, max_am_pairs, 2, pairs)) {
            return false;
        }
        const size_t count = static_cast<size_t>(std::max(0, counts[0]));
        clear_am_state(bonds);
        for (size_t idx = 0; idx < count; ++idx) {
            const int a = pairs[2 * idx];
            const int m = pairs[2 * idx + 1];
            if (a < 0 || a >= actin.n || m < 0 || m >= myosin.n) {
                continue;
            }
            bonds[a][m] = 1;
        }
        return true;
    };

    auto restore_recovery_state_from_pairs = [&](H5::Group& group,
                                                 const std::string& count_name,
                                                 const std::string& pair_name,
                                                 const std::string& value_name) -> bool {
        std::vector<int> counts;
        std::vector<int> pairs;
        std::vector<int> until_values;
        const size_t max_recovery_pairs = static_cast<size_t>(actin.n) *
                                          static_cast<size_t>(std::max(10, 2 * max_myosin_bonds));
        if (!load_state_vector(group, count_name, 1, counts) ||
            !load_state_tensor_int(group, pair_name, max_recovery_pairs, 2, pairs) ||
            !load_state_tensor_int(group, value_name, max_recovery_pairs, 1, until_values)) {
            return false;
        }
        const size_t count = static_cast<size_t>(std::max(0, counts[0]));
        clear_recovery_state();
        for (size_t idx = 0; idx < count; ++idx) {
            const int a = pairs[2 * idx];
            const int b = pairs[2 * idx + 1];
            if (a < 0 || a >= actin.n || b < 0 || b >= actin.n || a == b) {
                continue;
            }
            const RecoveryStep until =
                static_cast<RecoveryStep>(std::max(0, until_values[idx]));
            actin_recovery_until[a][b] = until;
            actin_recovery_until[b][a] = until;
        }
        return true;
    };

    auto restore_aa_state_from_fixed_pairs = [&](H5::Group& group,
                                                 const std::string& count_name,
                                                 const std::string& pair_name,
                                                 const std::string& status_name,
                                                 const std::string& lifetime_name,
                                                 auto& bonds,
                                                 auto& status,
                                                 auto& lifetime) -> bool {
        std::vector<int> counts;
        std::vector<int> pairs;
        std::vector<int> flat_status;
        std::vector<int> flat_lifetime;
        const size_t max_aa_pairs = static_cast<size_t>(actin.n) *
                                    static_cast<size_t>(std::max(10, 2 * max_myosin_bonds));
        if (!load_fixed_vector_int(group, count_name, 1, counts) ||
            !load_fixed_matrix_int(group, pair_name, max_aa_pairs, 2, pairs) ||
            !load_fixed_matrix_int(group, status_name, max_aa_pairs, 1, flat_status) ||
            !load_fixed_matrix_int(group, lifetime_name, max_aa_pairs, 1, flat_lifetime)) {
            return false;
        }
        const size_t count = static_cast<size_t>(std::max(0, counts[0]));
        clear_aa_state(bonds, status, lifetime);
        for (size_t idx = 0; idx < count; ++idx) {
            const int a = pairs[2 * idx];
            const int b = pairs[2 * idx + 1];
            if (a < 0 || a >= actin.n || b < 0 || b >= actin.n || a == b) {
                continue;
            }
            bonds[a][b] = 1;
            bonds[b][a] = 1;
            status[a][b] = flat_status[idx];
            status[b][a] = flat_status[idx];
            lifetime[a][b] = flat_lifetime[idx];
            lifetime[b][a] = flat_lifetime[idx];
        }
        return true;
    };

    auto restore_am_state_from_fixed_pairs = [&](H5::Group& group,
                                                 const std::string& count_name,
                                                 const std::string& pair_name,
                                                 auto& bonds) -> bool {
        std::vector<int> counts;
        std::vector<int> pairs;
        const size_t max_am_pairs = static_cast<size_t>(myosin.n) *
                                    static_cast<size_t>(max_myosin_bonds);
        if (!load_fixed_vector_int(group, count_name, 1, counts) ||
            !load_fixed_matrix_int(group, pair_name, max_am_pairs, 2, pairs)) {
            return false;
        }
        const size_t count = static_cast<size_t>(std::max(0, counts[0]));
        clear_am_state(bonds);
        for (size_t idx = 0; idx < count; ++idx) {
            const int a = pairs[2 * idx];
            const int m = pairs[2 * idx + 1];
            if (a < 0 || a >= actin.n || m < 0 || m >= myosin.n) {
                continue;
            }
            bonds[a][m] = 1;
        }
        return true;
    };

    auto restore_recovery_state_from_fixed_pairs = [&](H5::Group& group,
                                                       const std::string& count_name,
                                                       const std::string& pair_name,
                                                       const std::string& value_name) -> bool {
        std::vector<int> counts;
        std::vector<int> pairs;
        std::vector<int> until_values;
        const size_t max_recovery_pairs = static_cast<size_t>(actin.n) *
                                          static_cast<size_t>(std::max(10, 2 * max_myosin_bonds));
        if (!load_fixed_vector_int(group, count_name, 1, counts) ||
            !load_fixed_matrix_int(group, pair_name, max_recovery_pairs, 2, pairs) ||
            !load_fixed_matrix_int(group, value_name, max_recovery_pairs, 1, until_values)) {
            return false;
        }
        const size_t count = static_cast<size_t>(std::max(0, counts[0]));
        clear_recovery_state();
        for (size_t idx = 0; idx < count; ++idx) {
            const int a = pairs[2 * idx];
            const int b = pairs[2 * idx + 1];
            if (a < 0 || a >= actin.n || b < 0 || b >= actin.n || a == b) {
                continue;
            }
            const RecoveryStep until =
                static_cast<RecoveryStep>(std::max(0, until_values[idx]));
            actin_recovery_until[a][b] = until;
            actin_recovery_until[b][a] = until;
        }
        return true;
    };

    bool loaded_state_group = false;
    bool restored_neighbor_cache = false;
    bool loaded_resume_snapshot = false;
    try {
        H5::H5File file(filename, H5F_ACC_RDWR);
        H5::Group group_state(file.openGroup("/state"));
        loaded_state_group = true;

        std::vector<int> flat;
        std::vector<hsize_t> dims;

        // Detect interrupted /state appends and rewind to the last complete checkpoint.
        // state_commit_count is appended last in save_state(), so its frame count marks
        // the number of fully committed state snapshots.
        const std::string state_commit_dataset =
            group_state.nameExists("state_commit_count") ? "state_commit_count" : "rng_thread_state_size";
        if (group_state.nameExists(state_commit_dataset)) {
            std::vector<double> dummy = load_from_dataset(group_state, state_commit_dataset, dims);
            (void)dummy;
            if (dims.size() >= 2) {
                const int safe_frames = static_cast<int>(dims[0]);
                const std::vector<std::string> state_datasets = {
                    "current_step",
                    "state_commit_count",
                    "aa_pairs_prev_count",
                    "aa_pairs_current_count",
                    "am_pairs_prev_count",
                    "am_pairs_current_count",
                    "aa_pairs_prev",
                    "aa_status_prev",
                    "aa_lifetime_prev",
                    "aa_pairs_current",
                    "aa_status_current",
                    "aa_lifetime_current",
                    "aa_attach_step_current",
                    "aa_status2_memory_until_current",
                    "am_pairs_prev",
                    "am_pairs_current",
                    "aa_recovery_count",
                    "aa_recovery_pairs",
                    "aa_recovery_until_values",
                    "neighbor_last_actin_x",
                    "neighbor_last_actin_y",
                    "neighbor_last_actin_z",
                    "neighbor_last_myosin_x",
                    "neighbor_last_myosin_y",
                    "neighbor_last_myosin_z",
                    "titin_bound_actin_current",
                    "titin_last_bound_actin_current",
                    "titin_attach_step_current",
                    "titin_last_detach_step_current",
                    "titin_distance_current",
                    "titin_extension_current",
                    "titin_force_current",
                    "titin_attachment_distance_current",
                    "titin_attachment_energy_current",
                    "titin_attachment_force_current",
                    "rng_main_state",
                    "rng_thread_state",
                    "rng_thread_count",
                    "rng_thread_state_size"
                };
                if (safe_frames >= 0 && dims[0] > static_cast<hsize_t>(safe_frames)) {
                    for (const auto& name : state_datasets) {
                        truncate_dataset_first_dim(group_state, name, static_cast<hsize_t>(safe_frames));
                    }
                    file.flush(H5F_SCOPE_GLOBAL);
                }
                if (target_frame >= safe_frames) {
                    target_frame = safe_frames - 1;
                    if (target_frame >= 0) {
                        load_from_file(filename, actin, myosin, actin_actin_bonds, n_frames, target_frame);
                    }
                }
            }
        }

        // Load current_step (if missing, fallback handled below).
        if (group_state.nameExists("current_step")) {
            std::vector<double> step_data = load_from_dataset(group_state, "current_step", dims);
            if (dims.size() >= 2 && target_frame >= 0 && target_frame < static_cast<int>(dims[0])) {
                current_step = static_cast<size_t>(std::llround(step_data[target_frame]));
            } else {
                current_step = 0;
            }
        } else {
            current_step = 0;
        }

        bool have_prev_aa_pairs = restore_aa_state_from_pairs(
            group_state,
            "aa_pairs_prev_count",
            "aa_pairs_prev",
            "aa_status_prev",
            "aa_lifetime_prev",
            actin_actin_bonds_prev,
            actin_actin_status_prev,
            actin_actin_lifetime_prev);
        if (!have_prev_aa_pairs) {
            if (load_state_vector(group_state, "actin_actin_bonds_prev", aa_stride, flat)) {
                assign_aa_matrix(flat, actin_actin_bonds_prev);
            }
            if (load_state_vector(group_state, "actin_actin_status_prev", aa_stride, flat)) {
                assign_aa_matrix(flat, actin_actin_status_prev);
            }
            if (load_state_vector(group_state, "actin_actin_lifetime_prev", aa_stride, flat)) {
                assign_aa_matrix(flat, actin_actin_lifetime_prev);
            }
        }

        bool have_prev_am_pairs = restore_am_state_from_pairs(
            group_state, "am_pairs_prev_count", "am_pairs_prev", am_bonds_prev);
        if (!have_prev_am_pairs && load_state_vector(group_state, "am_bonds_prev", am_stride, flat)) {
            assign_am_matrix(flat, am_bonds_prev);
        }

        if (restore_recovery_state_from_pairs(
                group_state, "aa_recovery_count", "aa_recovery_pairs", "aa_recovery_until_values")) {
            // loaded sparse recovery state
        } else if (load_state_tensor_int(group_state, "actin_recovery_until", aa_stride, 1, flat) ||
                   load_state_vector(group_state, "actin_recovery_until", aa_stride, flat)) {
            for (int i = 0; i < actin.n; ++i) {
                for (int j = 0; j < actin.n; ++j) {
                    actin_recovery_until[i][j] = static_cast<RecoveryStep>(std::max(
                        0, flat[static_cast<size_t>(i) * actin.n + j]));
                }
            }
        } else {
            for (int i = 0; i < actin.n; ++i) {
                std::fill(actin_recovery_until[i].begin(), actin_recovery_until[i].end(), 0);
            }
        }

        bool have_current_aa_pairs = restore_aa_state_from_pairs(
            group_state,
            "aa_pairs_current_count",
            "aa_pairs_current",
            "aa_status_current",
            "aa_lifetime_current",
            actin_actin_bonds,
            actin_actin_status,
            actin_actin_lifetime);
        bool have_current_bonds = false;
        bool have_current_status = false;
        bool have_current_lifetime = false;
        if (!have_current_aa_pairs) {
            have_current_bonds = load_state_vector(group_state, "actin_actin_bonds_current", aa_stride, flat);
            if (have_current_bonds) {
                assign_aa_matrix(flat, actin_actin_bonds);
            }

            have_current_status = load_state_vector(group_state, "actin_actin_status_current", aa_stride, flat);
            if (have_current_status) {
                assign_aa_matrix(flat, actin_actin_status);
            } else {
                // Legacy fallback: /actin/bonds stores only strong bonds.
                for (int i = 0; i < actin.n; ++i) {
                    for (int j = 0; j < actin.n; ++j) {
                        actin_actin_status[i][j] = (actin_actin_bonds[i][j] == 1) ? 2 : 0;
                    }
                }
            }

            have_current_lifetime = load_state_vector(group_state, "actin_actin_lifetime_current", aa_stride, flat);
            if (have_current_lifetime) {
                assign_aa_matrix(flat, actin_actin_lifetime);
            } else {
                for (int i = 0; i < actin.n; ++i) {
                    for (int j = 0; j < actin.n; ++j) {
                        if (actin_actin_bonds[i][j] == 1) {
                            actin_actin_lifetime[i][j] = std::max(1, actin_actin_lifetime_prev[i][j]);
                        } else {
                            actin_actin_lifetime[i][j] = 0;
                        }
                    }
                }
            }
        }
        for (int i = 0; i < actin.n; ++i) {
            std::fill(aa_attach_step[i].begin(), aa_attach_step[i].end(), -1);
        }
        {
            std::vector<int> counts_attach;
            std::vector<int> pairs_attach;
            std::vector<int> attach_values;
            const size_t max_aa_pairs = static_cast<size_t>(actin.n) *
                                        static_cast<size_t>(std::max(10, 2 * max_myosin_bonds));
            if (load_state_vector(group_state, "aa_pairs_current_count", 1, counts_attach) &&
                load_state_tensor_int(group_state, "aa_pairs_current", max_aa_pairs, 2, pairs_attach) &&
                load_state_tensor_int(group_state, "aa_attach_step_current", max_aa_pairs, 1, attach_values)) {
                const size_t count = static_cast<size_t>(std::max(0, counts_attach[0]));
                for (size_t idx = 0; idx < count; ++idx) {
                    const int a = pairs_attach[2 * idx];
                    const int b = pairs_attach[2 * idx + 1];
                    if (a < 0 || a >= actin.n || b < 0 || b >= actin.n || a == b) {
                        continue;
                    }
                    aa_attach_step[a][b] = attach_values[idx];
                    aa_attach_step[b][a] = attach_values[idx];
                }
            }
        }
        for (int i = 0; i < actin.n; ++i) {
            std::fill(aa_status2_memory_until[i].begin(),
                      aa_status2_memory_until[i].end(), 0);
        }
        {
            std::vector<int> counts_memory;
            std::vector<int> pairs_memory;
            std::vector<int> memory_values;
            const size_t max_aa_pairs = static_cast<size_t>(actin.n) *
                                        static_cast<size_t>(std::max(10, 2 * max_myosin_bonds));
            if (load_state_vector(group_state, "aa_pairs_current_count", 1, counts_memory) &&
                load_state_tensor_int(group_state, "aa_pairs_current", max_aa_pairs, 2, pairs_memory) &&
                load_state_tensor_int(group_state, "aa_status2_memory_until_current",
                                      max_aa_pairs, 1, memory_values)) {
                const size_t count = static_cast<size_t>(std::max(0, counts_memory[0]));
                for (size_t idx = 0; idx < count; ++idx) {
                    const int a = pairs_memory[2 * idx];
                    const int b = pairs_memory[2 * idx + 1];
                    if (a < 0 || a >= actin.n || b < 0 || b >= actin.n || a == b) continue;
                    const RecoveryStep until = static_cast<RecoveryStep>(
                        std::max(0, memory_values[idx]));
                    aa_status2_memory_until[a][b] = until;
                    aa_status2_memory_until[b][a] = until;
                }
            }
        }

        bool loaded_current_am_from_pairs = restore_am_state_from_pairs(
            group_state, "am_pairs_current_count", "am_pairs_current", am_bonds);
        bool loaded_current_am_from_dense = false;
        if (!loaded_current_am_from_pairs) {
            loaded_current_am_from_dense =
                load_state_vector(group_state, "am_bonds_current", am_stride, flat);
        }
        if (loaded_current_am_from_dense) {
            assign_am_matrix(flat, am_bonds);
        } else if (!loaded_current_am_from_pairs) {
            for (int i = 0; i < actin.n; ++i) {
                std::fill(am_bonds[i].begin(), am_bonds[i].end(), 0);
            }
            // Fallback for legacy files: reconstruct current am_bonds from /actin_myo/bonds.
            try {
                H5::Group group_am(file.openGroup("/actin_myo"));
                std::vector<double> am_bonds_all = load_from_dataset(group_am, "bonds", dims);
                if (dims.size() >= 3 && target_frame >= 0 && target_frame < static_cast<int>(dims[0])) {
                    size_t bonds_per_frame = static_cast<size_t>(dims[1]) * static_cast<size_t>(dims[2]);
                    size_t start = static_cast<size_t>(target_frame) * bonds_per_frame;
                    for (size_t idx = 0; idx + 1 < bonds_per_frame; idx += 2) {
                        int a = static_cast<int>(std::llround(am_bonds_all[start + idx]));
                        int m = static_cast<int>(std::llround(am_bonds_all[start + idx + 1]));
                        if (a >= 0 && a < actin.n && m >= 0 && m < myosin.n) {
                            am_bonds[a][m] = 1;
                        }
                    }
                }
            } catch (H5::Exception&) {
                // Keep zeros if legacy AM bonds are unavailable.
            }
        }

        for (int i = 0; i < actin.n; ++i) {
            std::fill(am_state[i].begin(), am_state[i].end(), 0);
            std::fill(am_seen[i].begin(), am_seen[i].end(), 0);
            std::fill(am_target_axial[i].begin(), am_target_axial[i].end(), 0.0);
        }
        if (dynamic_am) {
            _initialize_am_head_pools();
            std::vector<double> material_e1_values;
            std::vector<double> material_e2_values;
            if (!load_state_tensor_double(
                    group_state, "myosin_material_e1_current",
                    static_cast<size_t>(myosin.n), 3, material_e1_values) ||
                !load_state_tensor_double(
                    group_state, "myosin_material_e2_current",
                    static_cast<size_t>(myosin.n), 3, material_e2_values)) {
                throw std::runtime_error(
                    "Dynamic AM schema v8 checkpoint is missing material frames");
            }
            for (int m = 0; m < myosin.n; ++m) {
                myosin_material_e1[static_cast<size_t>(m)] = {
                    material_e1_values[3 * m], material_e1_values[3 * m + 1],
                    material_e1_values[3 * m + 2]};
                myosin_material_e2[static_cast<size_t>(m)] = {
                    material_e2_values[3 * m], material_e2_values[3 * m + 1],
                    material_e2_values[3 * m + 2]};
                const vec axis = myosin.direction[m];
                transport_myosin_material_frame(m, axis, axis, 0.0);
            }
            std::vector<int> am_head_states;
            std::vector<int> am_head_bound_actins;
            std::vector<double> am_head_targets;
            std::vector<double> am_head_binding_azimuths;
            std::vector<double> am_head_stroke_progress;
            const size_t heads_per_myosin = static_cast<size_t>(2 * am_heads_per_half);
            const bool have_head_state =
                load_state_tensor_int(group_state, "am_head_state_current",
                                      static_cast<size_t>(myosin.n), heads_per_myosin,
                                      am_head_states) &&
                load_state_tensor_int(group_state, "am_head_bound_actin_current",
                                      static_cast<size_t>(myosin.n), heads_per_myosin,
                                      am_head_bound_actins) &&
                load_state_tensor_double(group_state, "am_head_target_axial_current",
                                         static_cast<size_t>(myosin.n), heads_per_myosin,
                                         am_head_targets);
            const bool have_binding_azimuths = load_state_tensor_double(
                group_state, "am_head_binding_azimuth_current",
                static_cast<size_t>(myosin.n), heads_per_myosin,
                am_head_binding_azimuths);
            const bool have_stroke_progress = load_state_tensor_double(
                group_state, "am_head_stroke_progress_current",
                static_cast<size_t>(myosin.n), heads_per_myosin,
                am_head_stroke_progress);
            if (!have_head_state) {
                throw std::runtime_error(
                    "Dynamic AM schema v5 checkpoint is missing myosin-owned head state");
            }
            for (int a = 0; a < actin.n; ++a) {
                std::fill(am_bonds[a].begin(), am_bonds[a].end(), 0);
            }
            for (int m = 0; m < myosin.n; ++m) {
                for (int head_id = 0; head_id < static_cast<int>(heads_per_myosin); ++head_id) {
                    const size_t flat_index = static_cast<size_t>(m) * heads_per_myosin +
                                              static_cast<size_t>(head_id);
                    auto* head = _am_head(m, head_id);
                    if (head == nullptr) continue;
                    head->state = static_cast<BondFlag>(
                        std::clamp(am_head_states[flat_index], 0, 3));
                    head->bound_actin = am_head_bound_actins[flat_index];
                    head->target_axial = am_head_targets[flat_index];
                    if (have_binding_azimuths) {
                        head->binding_azimuth =
                            am_head_binding_azimuths[flat_index];
                    }
                    head->stroke_progress = have_stroke_progress
                        ? std::clamp(am_head_stroke_progress[flat_index],
                                     0.0, am_power_stroke)
                        : (head->state == static_cast<BondFlag>(
                               am_motor_kinetics::State::strong)
                               ? am_power_stroke : 0.0);
                    if (head->state == static_cast<BondFlag>(
                            am_motor_kinetics::State::detached) ||
                        head->state == static_cast<BondFlag>(
                            am_motor_kinetics::State::off_srx)) {
                        head->bound_actin = -1;
                    } else if (head->bound_actin < 0 || head->bound_actin >= actin.n) {
                        throw std::runtime_error(
                            "Dynamic AM checkpoint contains an invalid bound-actin index");
                    } else {
                        am_bonds[head->bound_actin][m] = 1;
                    }
                }
            }
            for (int a = 0; a < actin.n; ++a) {
                for (int m = 0; m < myosin.n; ++m) {
                    if (am_bonds[a][m] != 0) _sync_am_pair_summary(a, m);
                }
            }
            std::vector<double> restored_q;
            if (!load_state_vector_double(
                    group_state, "am_q_current",
                    static_cast<size_t>(actin.n), restored_q)) {
                throw std::runtime_error(
                    "Dynamic AM schema v4 checkpoint is missing thin activation q_i");
            }
            am_thin_activation = std::move(restored_q);
        }
        {
            std::vector<int> restored_turnover_steps;
            if (load_state_vector(
                    group_state, "actin_last_turnover_step_current",
                    static_cast<size_t>(actin.n), restored_turnover_steps)) {
                actin_last_turnover_step.resize(static_cast<size_t>(actin.n));
                for (int i = 0; i < actin.n; ++i) {
                    actin_last_turnover_step[static_cast<size_t>(i)] =
                        static_cast<RecoveryStep>(std::max(
                            0, restored_turnover_steps[static_cast<size_t>(i)]));
                }
            } else if (actin_turnover_enabled) {
                throw std::runtime_error(
                    "Turnover checkpoint is missing actin renewal ages");
            } else {
                actin_last_turnover_step.assign(
                    static_cast<size_t>(actin.n), 0);
            }
        }

        if (titin_mode == TitinMode::explicit_dynamic) {
            const size_t titin_count = explicit_titins.size();
            std::vector<int> bound_actin;
            std::vector<int> last_bound_actin;
            std::vector<int> attach_step;
            std::vector<int> last_detach_step;
            std::vector<double> distance;
            std::vector<double> extension;
            std::vector<double> force;
            std::vector<double> attachment_distance;
            std::vector<double> attachment_energy;
            std::vector<double> attachment_force;
            const bool have_titin_state =
                load_state_vector(group_state, "titin_bound_actin_current",
                                  titin_count, bound_actin) &&
                load_state_vector(group_state,
                                  "titin_last_bound_actin_current",
                                  titin_count, last_bound_actin) &&
                load_state_vector(group_state, "titin_attach_step_current",
                                  titin_count, attach_step) &&
                load_state_vector(group_state,
                                  "titin_last_detach_step_current",
                                  titin_count, last_detach_step) &&
                load_state_vector_double(group_state, "titin_distance_current",
                                         titin_count, distance) &&
                load_state_vector_double(group_state, "titin_extension_current",
                                         titin_count, extension) &&
                load_state_vector_double(group_state, "titin_force_current",
                                         titin_count, force) &&
                load_state_vector_double(
                    group_state, "titin_attachment_distance_current",
                    titin_count, attachment_distance) &&
                load_state_vector_double(
                    group_state, "titin_attachment_energy_current",
                    titin_count, attachment_energy) &&
                load_state_vector_double(
                    group_state, "titin_attachment_force_current",
                    titin_count, attachment_force);
            if (!have_titin_state) {
                throw std::runtime_error(
                    "Explicit titin checkpoint is missing molecular state");
            }
            for (size_t index = 0; index < titin_count; ++index) {
                auto& titin = explicit_titins[index];
                titin.bound_actin = bound_actin[index];
                titin.bound = titin.bound_actin >= 0;
                titin.last_bound_actin = last_bound_actin[index];
                titin.attach_step = static_cast<RecoveryStep>(
                    std::max(0, attach_step[index]));
                titin.last_detach_step = static_cast<RecoveryStep>(
                    std::max(0, last_detach_step[index]));
                titin.distance = distance[index];
                titin.extension = extension[index];
                titin.force = force[index];
                titin.attachment_distance = attachment_distance[index];
                titin.attachment_energy = attachment_energy[index];
                titin.attachment_force = attachment_force[index];
                titin.candidate_actin = -1;
                titin.candidate_distance = 0.0;
                titin.candidate_on_rate = 0.0;
            }
            rebuild_explicit_titin_ownership();
        }

        myomesin_bridges.clear();
        if (k_mm > 0.0) {
            if (!file.nameExists("/myomesin/frame_offsets") ||
                !file.nameExists("/myomesin/bridge_state")) {
                throw std::runtime_error(
                    "Sparse myomesin checkpoint is missing molecular state");
            }
            H5::Group myomesin_group(file.openGroup("/myomesin"));
            std::vector<hsize_t> offset_dims;
            const std::vector<double> offsets = load_from_dataset(
                myomesin_group, "frame_offsets", offset_dims);
            if (offset_dims.size() != 2 || offset_dims[1] != 3 ||
                target_frame < 0 ||
                target_frame >= static_cast<int>(offset_dims[0])) {
                throw std::runtime_error(
                    "Sparse myomesin frame offsets do not match trajectory");
            }
            const size_t offset_base =
                static_cast<size_t>(target_frame) * 3;
            const size_t start = static_cast<size_t>(std::max(
                0.0, offsets[offset_base + 1]));
            const size_t count = static_cast<size_t>(std::max(
                0.0, offsets[offset_base + 2]));
            std::vector<hsize_t> state_dims;
            const std::vector<double> states = load_from_dataset(
                myomesin_group, "bridge_state", state_dims);
            if (state_dims.size() != 2 || state_dims[1] != 8 ||
                start + count > static_cast<size_t>(state_dims[0])) {
                throw std::runtime_error(
                    "Sparse myomesin frame contains invalid state bounds");
            }
            std::vector<int> degree(static_cast<size_t>(myosin.n), 0);
            std::unordered_set<std::uint64_t> pairs;
            std::uint64_t largest_id = 0;
            for (size_t record = 0; record < count; ++record) {
                const size_t base = 8 * (start + record);
                MyomesinBridge bridge;
                bridge.id =
                    (static_cast<std::uint64_t>(static_cast<std::uint32_t>(
                         std::llround(states[base]))) << 32) |
                    static_cast<std::uint32_t>(
                        std::llround(states[base + 1]));
                bridge.myosin_i = static_cast<int>(
                    std::llround(states[base + 2]));
                bridge.myosin_j = static_cast<int>(
                    std::llround(states[base + 3]));
                bridge.attach_step = static_cast<RecoveryStep>(std::max(
                    0.0, states[base + 4]));
                bridge.axial_offset = states[base + 5];
                bridge.radial_distance = states[base + 6];
                bridge.force = states[base + 7];
                if (bridge.id == 0 || bridge.myosin_i < 0 ||
                    bridge.myosin_i >= myosin.n ||
                    bridge.myosin_j <= bridge.myosin_i ||
                    bridge.myosin_j >= myosin.n) {
                    throw std::runtime_error(
                        "Sparse myomesin checkpoint contains invalid bridge");
                }
                const auto key = counter_rng::canonical_pair_key(
                    bridge.myosin_i, bridge.myosin_j);
                if (!pairs.insert(key).second ||
                    ++degree[static_cast<size_t>(bridge.myosin_i)] >
                        myomesin_sites_per_filament ||
                    ++degree[static_cast<size_t>(bridge.myosin_j)] >
                        myomesin_sites_per_filament) {
                    throw std::runtime_error(
                        "Sparse myomesin checkpoint violates valency");
                }
                largest_id = std::max(largest_id, bridge.id);
                myomesin_bridges.push_back(bridge);
            }
            myomesin_next_id = largest_id + 1;
            if (myomesin_group.nameExists("next_id_by_frame")) {
                std::vector<hsize_t> next_dims;
                const std::vector<double> next_ids = load_from_dataset(
                    myomesin_group, "next_id_by_frame", next_dims);
                if (next_dims.size() != 2 || next_dims[1] != 2 ||
                    target_frame >= static_cast<int>(next_dims[0])) {
                    throw std::runtime_error(
                        "Sparse myomesin next-ID state does not match trajectory");
                }
                const size_t next_base =
                    static_cast<size_t>(target_frame) * 2;
                myomesin_next_id =
                    (static_cast<std::uint64_t>(
                         static_cast<std::uint32_t>(std::llround(
                             next_ids[next_base]))) << 32) |
                    static_cast<std::uint32_t>(std::llround(
                        next_ids[next_base + 1]));
            } else if (myomesin_group.nameExists("events")) {
                // Compatibility with schema-1 files written before the next
                // ID was stored per frame.  Events retain detached molecules
                // that are absent from the occupied-state rows.
                std::vector<hsize_t> event_dims;
                const std::vector<double> events = load_from_dataset(
                    myomesin_group, "events", event_dims);
                if (event_dims.size() == 2 && event_dims[1] == 9) {
                    for (size_t row = 0;
                         row < static_cast<size_t>(event_dims[0]); ++row) {
                        if (events[9 * row] >
                            static_cast<double>(current_step)) {
                            continue;
                        }
                        const auto event_id = static_cast<std::uint64_t>(
                            std::llround(events[9 * row + 2]));
                        myomesin_next_id = std::max(
                            myomesin_next_id, event_id + 1);
                    }
                }
            }
        }

        // Restore per-actin cb_status from the matrix so diagnostic output remains consistent.
        for (int i = 0; i < actin.n; ++i) {
            int max_status = 0;
            for (int j = 0; j < actin.n; ++j) {
                max_status = std::max(max_status, static_cast<int>(actin_actin_status[i][j]));
            }
            actin.cb_status[i] = max_status;
        }

        // Restore RNG states for exact reproducibility.
        if (rng != nullptr) {
            if (group_state.nameExists("rng_main_state")) {
                std::vector<double> rng_main_data = load_from_dataset(group_state, "rng_main_state", dims);
                const size_t expected = gsl_rng_size(rng);
                if (dims.size() >= 2 && target_frame >= 0 && target_frame < static_cast<int>(dims[0]) &&
                    static_cast<size_t>(dims[1]) == expected) {
                    size_t start = static_cast<size_t>(target_frame) * expected;
                    auto* state_ptr = static_cast<unsigned char*>(gsl_rng_state(rng));
                    for (size_t idx = 0; idx < expected; ++idx) {
                        int byte_value = static_cast<int>(std::llround(rng_main_data[start + idx]));
                        byte_value = std::clamp(byte_value, 0, 255);
                        state_ptr[idx] = static_cast<unsigned char>(byte_value);
                    }
                }
            }
        }
        if (!rng_engines.empty() && rng_engines[0] != nullptr &&
            group_state.nameExists("rng_thread_state")) {
                std::vector<double> thread_data = load_from_dataset(group_state, "rng_thread_state", dims);
                if (dims.size() >= 2 && target_frame >= 0 && target_frame < static_cast<int>(dims[0])) {
                    int file_thread_count = 0;
                    int file_state_size = 0;
                    if (group_state.nameExists("rng_thread_count") && group_state.nameExists("rng_thread_state_size")) {
                        std::vector<hsize_t> scalar_dims;
                        std::vector<double> tc = load_from_dataset(group_state, "rng_thread_count", scalar_dims);
                        if (scalar_dims.size() >= 2 && target_frame < static_cast<int>(scalar_dims[0])) {
                            file_thread_count = static_cast<int>(std::llround(tc[target_frame]));
                        }
                        std::vector<double> ts = load_from_dataset(group_state, "rng_thread_state_size", scalar_dims);
                        if (scalar_dims.size() >= 2 && target_frame < static_cast<int>(scalar_dims[0])) {
                            file_state_size = static_cast<int>(std::llround(ts[target_frame]));
                        }
                    }

                    if (file_thread_count <= 0) {
                        file_thread_count = static_cast<int>(rng_engines.size());
                    }
                    if (file_state_size <= 0) {
                        if (file_thread_count > 0) {
                            file_state_size = static_cast<int>(dims[1] / static_cast<hsize_t>(file_thread_count));
                        } else {
                            file_state_size = 0;
                        }
                    }

                    const size_t frame_width = static_cast<size_t>(dims[1]);
                    const size_t frame_start = static_cast<size_t>(target_frame) * frame_width;
                    const int local_thread_count = static_cast<int>(rng_engines.size());
                    const int restore_threads = std::min(local_thread_count, file_thread_count);
                    for (int t = 0; t < restore_threads; ++t) {
                        if (rng_engines[t] == nullptr) {
                            continue;
                        }
                        const size_t local_state_size = gsl_rng_size(rng_engines[t]);
                        const size_t copy_size = std::min(local_state_size, static_cast<size_t>(file_state_size));
                        auto* local_ptr = static_cast<unsigned char*>(gsl_rng_state(rng_engines[t]));
                        size_t thread_offset = frame_start + static_cast<size_t>(t) * static_cast<size_t>(file_state_size);
                        if (thread_offset + copy_size > thread_data.size()) {
                            break;
                        }
                        for (size_t idx = 0; idx < copy_size; ++idx) {
                            int byte_value = static_cast<int>(std::llround(thread_data[thread_offset + idx]));
                            byte_value = std::clamp(byte_value, 0, 255);
                            local_ptr[idx] = static_cast<unsigned char>(byte_value);
                        }
                    }
                    if (local_thread_count > file_thread_count) {
                        for (int t = file_thread_count; t < local_thread_count; ++t) {
                            if (rng_engines[t] == nullptr) {
                                continue;
                            }
                            gsl_rng_set(
                                rng_engines[t],
                                static_cast<unsigned long>(initial_seed) +
                                    static_cast<unsigned long>(t) +
                                    static_cast<unsigned long>(current_step));
                        }
                    }
                }
        }

        std::vector<double> neighbor_last_actin_x;
        std::vector<double> neighbor_last_actin_y;
        std::vector<double> neighbor_last_actin_z;
        std::vector<double> neighbor_last_myosin_x;
        std::vector<double> neighbor_last_myosin_y;
        std::vector<double> neighbor_last_myosin_z;
        bool have_neighbor_cache =
            load_state_vector_double(group_state, "neighbor_last_actin_x", static_cast<size_t>(actin.n), neighbor_last_actin_x) &&
            load_state_vector_double(group_state, "neighbor_last_actin_y", static_cast<size_t>(actin.n), neighbor_last_actin_y) &&
            load_state_vector_double(group_state, "neighbor_last_actin_z", static_cast<size_t>(actin.n), neighbor_last_actin_z) &&
            load_state_vector_double(group_state, "neighbor_last_myosin_x", static_cast<size_t>(myosin.n), neighbor_last_myosin_x) &&
            load_state_vector_double(group_state, "neighbor_last_myosin_y", static_cast<size_t>(myosin.n), neighbor_last_myosin_y) &&
            load_state_vector_double(group_state, "neighbor_last_myosin_z", static_cast<size_t>(myosin.n), neighbor_last_myosin_z);
        if (have_neighbor_cache) {
            // Reconstruct cached neighbor pairs from the original rebuild reference positions.
            neighbor_list.set_species_positions(neighbor_last_actin_x, neighbor_last_actin_y, neighbor_last_actin_z,
                                                neighbor_last_myosin_x, neighbor_last_myosin_y, neighbor_last_myosin_z);
            neighbor_list.rebuild_neighbor_list();
            // Then restore current particle positions while preserving the cached "last" positions.
            neighbor_list.set_species_positions(actin.center_x, actin.center_y, actin.center_z,
                                                myosin.center_x, myosin.center_y, myosin.center_z);
            restored_neighbor_cache = true;
        }

        if (frame_index < 0) {
            const std::string resume_group_path =
                file.nameExists("/resume") ? "/resume" :
                (file.nameExists("/resume_prev") ? "/resume_prev" : "");
            if (!resume_group_path.empty()) {
                try {
                    if (resume_group_path == "/resume_prev") {
                    }
                    H5::Group group_resume(file.openGroup(resume_group_path));
                    std::vector<int> ints;
                    std::vector<double> doubles;
                    size_t resume_step_value = current_step;

                    if (load_fixed_vector_int(group_resume, "current_step", 1, ints)) {
                        const size_t resume_step = static_cast<size_t>(std::max(0, ints[0]));
                        if (resume_step < current_step) {
                            throw H5::Exception("Sarcomere::load_state", "Stale resume snapshot");
                        }
                        current_step = resume_step;
                        resume_step_value = resume_step;
                    }
                    if (load_fixed_matrix_double(group_resume, "actin_center", static_cast<size_t>(actin.n), 3, doubles)) {
                        for (int i = 0; i < actin.n; ++i) {
                            actin.center[i].x = doubles[3 * i];
                            actin.center[i].y = doubles[3 * i + 1];
                            actin.center[i].z = doubles[3 * i + 2];
                        }
                    }
                    if (load_fixed_matrix_double(group_resume, "actin_direction", static_cast<size_t>(actin.n), 3, doubles)) {
                        for (int i = 0; i < actin.n; ++i) {
                            actin.direction[i].x = doubles[3 * i];
                            actin.direction[i].y = doubles[3 * i + 1];
                            actin.direction[i].z = doubles[3 * i + 2];
                        }
                    }
                    if (load_fixed_matrix_double(group_resume, "myosin_center", static_cast<size_t>(myosin.n), 3, doubles)) {
                        for (int i = 0; i < myosin.n; ++i) {
                            myosin.center[i].x = doubles[3 * i];
                            myosin.center[i].y = doubles[3 * i + 1];
                            myosin.center[i].z = doubles[3 * i + 2];
                        }
                    }
                    if (load_fixed_matrix_double(group_resume, "myosin_direction", static_cast<size_t>(myosin.n), 3, doubles)) {
                        for (int i = 0; i < myosin.n; ++i) {
                            myosin.direction[i].x = doubles[3 * i];
                            myosin.direction[i].y = doubles[3 * i + 1];
                            myosin.direction[i].z = doubles[3 * i + 2];
                        }
                    }
                    std::vector<double> material_e1_values;
                    std::vector<double> material_e2_values;
                    if (!load_fixed_matrix_double(
                            group_resume, "myosin_material_e1_current",
                            static_cast<size_t>(myosin.n), 3,
                            material_e1_values) ||
                        !load_fixed_matrix_double(
                            group_resume, "myosin_material_e2_current",
                            static_cast<size_t>(myosin.n), 3,
                            material_e2_values)) {
                        throw std::runtime_error(
                            "Resume snapshot is missing transported myosin frames");
                    }
                    for (int m = 0; m < myosin.n; ++m) {
                        myosin_material_e1[static_cast<size_t>(m)] = {
                            material_e1_values[3 * m],
                            material_e1_values[3 * m + 1],
                            material_e1_values[3 * m + 2]};
                        myosin_material_e2[static_cast<size_t>(m)] = {
                            material_e2_values[3 * m],
                            material_e2_values[3 * m + 1],
                            material_e2_values[3 * m + 2]};
                        const vec axis = myosin.direction[m];
                        transport_myosin_material_frame(m, axis, axis, 0.0);
                    }
                    actin.update_endpoints();
                    myosin.update_endpoints();

                    restore_aa_state_from_fixed_pairs(
                        group_resume,
                        "aa_pairs_prev_count",
                        "aa_pairs_prev",
                        "aa_status_prev",
                        "aa_lifetime_prev",
                        actin_actin_bonds_prev,
                        actin_actin_status_prev,
                        actin_actin_lifetime_prev);
                    restore_aa_state_from_fixed_pairs(
                        group_resume,
                        "aa_pairs_current_count",
                        "aa_pairs_current",
                        "aa_status_current",
                        "aa_lifetime_current",
                        actin_actin_bonds,
                        actin_actin_status,
                        actin_actin_lifetime);
                    for (int i = 0; i < actin.n; ++i) {
                        std::fill(aa_attach_step[i].begin(), aa_attach_step[i].end(), -1);
                    }
                    {
                        std::vector<int> counts_attach;
                        std::vector<int> pairs_attach;
                        std::vector<int> attach_values;
                        const size_t max_aa_pairs = static_cast<size_t>(actin.n) *
                                                    static_cast<size_t>(std::max(10, 2 * max_myosin_bonds));
                        if (load_fixed_vector_int(group_resume, "aa_pairs_current_count", 1, counts_attach) &&
                            load_fixed_matrix_int(group_resume, "aa_pairs_current", max_aa_pairs, 2, pairs_attach) &&
                            load_fixed_matrix_int(group_resume, "aa_attach_step_current", max_aa_pairs, 1, attach_values)) {
                            const size_t count = static_cast<size_t>(std::max(0, counts_attach[0]));
                            for (size_t idx = 0; idx < count; ++idx) {
                                const int a = pairs_attach[2 * idx];
                                const int b = pairs_attach[2 * idx + 1];
                                if (a < 0 || a >= actin.n || b < 0 || b >= actin.n || a == b) {
                                    continue;
                                }
                                aa_attach_step[a][b] = attach_values[idx];
                                aa_attach_step[b][a] = attach_values[idx];
                            }
                        }
                    }
                    for (int i = 0; i < actin.n; ++i) {
                        std::fill(aa_status2_memory_until[i].begin(),
                                  aa_status2_memory_until[i].end(), 0);
                    }
                    {
                        std::vector<int> counts_memory;
                        std::vector<int> pairs_memory;
                        std::vector<int> memory_values;
                        const size_t max_aa_pairs = static_cast<size_t>(actin.n) *
                                                    static_cast<size_t>(std::max(10, 2 * max_myosin_bonds));
                        if (load_fixed_vector_int(group_resume, "aa_pairs_current_count", 1, counts_memory) &&
                            load_fixed_matrix_int(group_resume, "aa_pairs_current", max_aa_pairs, 2, pairs_memory) &&
                            load_fixed_matrix_int(group_resume, "aa_status2_memory_until_current",
                                                  max_aa_pairs, 1, memory_values)) {
                            const size_t count = static_cast<size_t>(std::max(0, counts_memory[0]));
                            for (size_t idx = 0; idx < count; ++idx) {
                                const int a = pairs_memory[2 * idx];
                                const int b = pairs_memory[2 * idx + 1];
                                if (a < 0 || a >= actin.n || b < 0 || b >= actin.n || a == b) continue;
                                const RecoveryStep until = static_cast<RecoveryStep>(
                                    std::max(0, memory_values[idx]));
                                aa_status2_memory_until[a][b] = until;
                                aa_status2_memory_until[b][a] = until;
                            }
                        }
                    }

                    aa_links.clear();
                    aa_free_link_slots.clear();
                    aa_active_link_slots.clear();
                    aa_links_by_pair.clear();
                    aa_links_by_actin.assign(static_cast<size_t>(actin.n), {});
                    aa_next_link_id = 1;
                    if (aa_local_multilink) {
                        std::vector<int> local_counts;
                        std::vector<int> local_identity;
                        std::vector<double> local_coordinates;
                        std::vector<int> local_motor_identity;
                        std::vector<double> local_motor_evidence;
                        std::vector<int> local_arm_state;
                        std::vector<double> local_arm_loads;
                        std::vector<int> local_next_id;
                        if (!load_fixed_vector_int(
                                group_resume, "aa_local_link_count", 1,
                                local_counts)) {
                            throw std::runtime_error(
                                "Resume snapshot is missing local ACTN2 link count");
                        }
                        const size_t local_count = static_cast<size_t>(
                            std::max(0, local_counts[0]));
                        const size_t local_rows = std::max<size_t>(1, local_count);
                        if (!load_fixed_matrix_int(
                                group_resume, "aa_local_link_identity",
                                local_rows, 7, local_identity) ||
                            !load_fixed_matrix_double(
                                group_resume, "aa_local_link_coordinates",
                                local_rows, 2, local_coordinates) ||
                            !load_fixed_vector_int(
                                group_resume, "aa_local_next_id", 2,
                                local_next_id)) {
                            throw std::runtime_error(
                                "Resume snapshot is missing local ACTN2 link state");
                        }
                        const bool needs_motor_identity =
                            aa_status_mode == actn2_load_kinetics::StatusMode::
                                motor_transmitted_register;
                        const bool has_motor_identity =
                            load_fixed_matrix_int(
                                group_resume,
                                "aa_local_link_motor_identity",
                                local_rows, 12, local_motor_identity) &&
                            load_fixed_matrix_double(
                                group_resume,
                                "aa_local_link_motor_evidence",
                                local_rows, 7, local_motor_evidence);
                        if (needs_motor_identity && !has_motor_identity) {
                            throw std::runtime_error(
                                "Resume snapshot is missing identity-resolved "
                                "ACTN2 motor-pair evidence");
                        }
                        const bool has_arm_state =
                            load_fixed_matrix_int(
                                group_resume, "aa_local_link_arm_state",
                                local_rows, 3, local_arm_state) &&
                            load_fixed_matrix_double(
                                group_resume, "aa_local_link_arm_loads",
                                local_rows, 4, local_arm_loads);
                        if (aa_arm_state_model && !has_arm_state) {
                            throw std::runtime_error(
                                "Resume snapshot is missing ACTN2 two-arm state");
                        }
                        aa_links.reserve(local_count);
                        for (size_t record = 0; record < local_count; ++record) {
                            AACrosslink link;
                            link.id =
                                (static_cast<std::uint64_t>(static_cast<std::uint32_t>(
                                     local_identity[7 * record])) << 32) |
                                static_cast<std::uint32_t>(
                                    local_identity[7 * record + 1]);
                            link.actin_i = local_identity[7 * record + 2];
                            link.actin_j = local_identity[7 * record + 3];
                            const int status = local_identity[7 * record + 4];
                            if (link.actin_i < 0 || link.actin_i >= actin.n ||
                                link.actin_j < 0 || link.actin_j >= actin.n ||
                                link.actin_i == link.actin_j ||
                                (status != 1 && status != 2)) {
                                throw std::runtime_error(
                                    "Invalid local ACTN2 link in resume snapshot");
                            }
                            link.status = static_cast<BondFlag>(status);
                            link.attach_step = static_cast<RecoveryStep>(
                                static_cast<std::uint32_t>(
                                    local_identity[7 * record + 5]));
                            link.status2_memory_until = static_cast<RecoveryStep>(
                                static_cast<std::uint32_t>(
                                    local_identity[7 * record + 6]));
                            link.s_i = local_coordinates[2 * record];
                            link.s_j = local_coordinates[2 * record + 1];
                            if (has_arm_state) {
                                const int arm_state_value =
                                    local_arm_state[3 * record];
                                if (arm_state_value < static_cast<int>(
                                        actn2_arm_kinetics::State::weak_bridge) ||
                                    arm_state_value > static_cast<int>(
                                        actn2_arm_kinetics::State::tethered_to_j)) {
                                    throw std::runtime_error(
                                        "Resume snapshot contains invalid ACTN2 arm state");
                                }
                                link.arm_state = static_cast<
                                    actn2_arm_kinetics::State>(arm_state_value);
                                link.arm_state_since = static_cast<RecoveryStep>(
                                    std::max(0, local_arm_state[3 * record + 1]));
                                link.last_partner =
                                    local_arm_state[3 * record + 2];
                                link.arm_force_i = local_arm_loads[4 * record];
                                link.arm_force_j =
                                    local_arm_loads[4 * record + 1];
                                link.contact_force =
                                    local_arm_loads[4 * record + 2];
                                link.contact_penetration =
                                    local_arm_loads[4 * record + 3];
                            } else {
                                link.arm_state =
                                    actn2_arm_kinetics::bridge_state_for_status(
                                        status);
                                link.arm_state_since = link.attach_step;
                            }
                            if (has_motor_identity) {
                                link.candidate_myosin_i =
                                    local_motor_identity[12 * record];
                                link.candidate_myosin_j =
                                    local_motor_identity[12 * record + 1];
                                link.candidate_half_i = static_cast<int8_t>(
                                    local_motor_identity[12 * record + 2]);
                                link.candidate_half_j = static_cast<int8_t>(
                                    local_motor_identity[12 * record + 3]);
                                link.candidate_since =
                                    static_cast<RecoveryStep>(
                                        static_cast<std::uint32_t>(
                                            local_motor_identity[
                                                12 * record + 4]));
                                link.candidate_last_seen =
                                    static_cast<RecoveryStep>(
                                        static_cast<std::uint32_t>(
                                            local_motor_identity[
                                                12 * record + 5]));
                                link.candidate_locked =
                                    local_motor_identity[12 * record + 6] != 0;
                                link.candidate_support_ema =
                                    local_motor_evidence[7 * record];
                                link.candidate_support_instantaneous =
                                    local_motor_evidence[7 * record + 1];
                                link.candidate_raw_support_i =
                                    local_motor_evidence[7 * record + 2];
                                link.candidate_raw_support_j =
                                    local_motor_evidence[7 * record + 3];
                                link.candidate_allocated_support_i =
                                    local_motor_evidence[7 * record + 4];
                                link.candidate_allocated_support_j =
                                    local_motor_evidence[7 * record + 5];
                                link.challenger_myosin_i =
                                    local_motor_identity[12 * record + 7];
                                link.challenger_myosin_j =
                                    local_motor_identity[12 * record + 8];
                                link.challenger_half_i = static_cast<int8_t>(
                                    local_motor_identity[12 * record + 9]);
                                link.challenger_half_j = static_cast<int8_t>(
                                    local_motor_identity[12 * record + 10]);
                                link.challenger_since =
                                    static_cast<RecoveryStep>(
                                        static_cast<std::uint32_t>(
                                            local_motor_identity[
                                                12 * record + 11]));
                                link.challenger_support_ema =
                                    local_motor_evidence[7 * record + 6];
                            }
                            link.active = true;
                            aa_links.push_back(link);
                        }
                        aa_next_link_id =
                            (static_cast<std::uint64_t>(static_cast<std::uint32_t>(
                                 local_next_id[0])) << 32) |
                            static_cast<std::uint32_t>(local_next_id[1]);
                        _rebuild_local_aa_indices();
                    }

                    restore_am_state_from_fixed_pairs(
                        group_resume, "am_pairs_prev_count", "am_pairs_prev", am_bonds_prev);
                    restore_am_state_from_fixed_pairs(
                        group_resume, "am_pairs_current_count", "am_pairs_current", am_bonds);

                    if (dynamic_am) {
                        for (int a = 0; a < actin.n; ++a) {
                            std::fill(am_state[a].begin(), am_state[a].end(), 0);
                            std::fill(am_seen[a].begin(), am_seen[a].end(), 0);
                            std::fill(am_target_axial[a].begin(), am_target_axial[a].end(), 0.0);
                        }
                        _initialize_am_head_pools();
                        std::vector<int> am_head_states;
                        std::vector<int> am_head_bound_actins;
                        std::vector<double> am_head_targets;
                        std::vector<double> am_head_binding_azimuths;
                        std::vector<double> am_head_stroke_progress;
                        const size_t heads_per_myosin =
                            static_cast<size_t>(2 * am_heads_per_half);
                        const bool have_head_state =
                            load_fixed_matrix_int(
                                group_resume, "am_head_state_current",
                                static_cast<size_t>(myosin.n), heads_per_myosin,
                                am_head_states) &&
                            load_fixed_matrix_int(
                                group_resume, "am_head_bound_actin_current",
                                static_cast<size_t>(myosin.n), heads_per_myosin,
                                am_head_bound_actins) &&
                            load_fixed_matrix_double(
                                group_resume, "am_head_target_axial_current",
                                static_cast<size_t>(myosin.n), heads_per_myosin,
                                am_head_targets);
                        const bool have_binding_azimuths =
                            load_fixed_matrix_double(
                                group_resume,
                                "am_head_binding_azimuth_current",
                                static_cast<size_t>(myosin.n),
                                heads_per_myosin,
                                am_head_binding_azimuths);
                        const bool have_stroke_progress =
                            load_fixed_matrix_double(
                                group_resume,
                                "am_head_stroke_progress_current",
                                static_cast<size_t>(myosin.n),
                                heads_per_myosin,
                                am_head_stroke_progress);
                        if (!have_head_state) {
                            throw std::runtime_error(
                                "Resume snapshot is missing myosin-owned AM head state");
                        }
                        for (int a = 0; a < actin.n; ++a) {
                            std::fill(am_bonds[a].begin(), am_bonds[a].end(), 0);
                        }
                        for (int m = 0; m < myosin.n; ++m) {
                            for (int head_id = 0;
                                 head_id < static_cast<int>(heads_per_myosin); ++head_id) {
                                const size_t flat_index =
                                    static_cast<size_t>(m) * heads_per_myosin +
                                    static_cast<size_t>(head_id);
                                auto* head = _am_head(m, head_id);
                                if (head == nullptr) continue;
                                head->state = static_cast<BondFlag>(
                                    std::clamp(am_head_states[flat_index], 0, 3));
                                head->bound_actin = am_head_bound_actins[flat_index];
                                head->target_axial = am_head_targets[flat_index];
                                if (have_binding_azimuths) {
                                    head->binding_azimuth =
                                        am_head_binding_azimuths[flat_index];
                                }
                                head->stroke_progress = have_stroke_progress
                                    ? std::clamp(
                                          am_head_stroke_progress[flat_index],
                                          0.0, am_power_stroke)
                                    : (head->state == static_cast<BondFlag>(
                                           am_motor_kinetics::State::strong)
                                           ? am_power_stroke : 0.0);
                                if (head->state == static_cast<BondFlag>(
                                        am_motor_kinetics::State::detached) ||
                                    head->state == static_cast<BondFlag>(
                                        am_motor_kinetics::State::off_srx)) {
                                    head->bound_actin = -1;
                                } else if (head->bound_actin < 0 ||
                                           head->bound_actin >= actin.n) {
                                    throw std::runtime_error(
                                        "Resume snapshot contains invalid AM head ownership");
                                } else {
                                    am_bonds[head->bound_actin][m] = 1;
                                }
                            }
                        }
                        for (int a = 0; a < actin.n; ++a) {
                            for (int m = 0; m < myosin.n; ++m) {
                                if (am_bonds[a][m] != 0) {
                                    _sync_am_pair_summary(a, m);
                                }
                            }
                        }
                        std::vector<double> restored_q;
                        if (!load_fixed_vector_double(
                                group_resume, "am_q_current",
                                static_cast<size_t>(actin.n), restored_q)) {
                            throw std::runtime_error(
                                "Resume snapshot is missing thin activation q_i");
                        }
                        am_thin_activation = std::move(restored_q);
                    }
                    {
                        std::vector<int> restored_turnover_steps;
                        if (load_fixed_vector_int(
                                group_resume,
                                "actin_last_turnover_step_current",
                                static_cast<size_t>(actin.n),
                                restored_turnover_steps)) {
                            actin_last_turnover_step.resize(
                                static_cast<size_t>(actin.n));
                            for (int i = 0; i < actin.n; ++i) {
                                actin_last_turnover_step[static_cast<size_t>(i)] =
                                    static_cast<RecoveryStep>(std::max(
                                        0, restored_turnover_steps[
                                               static_cast<size_t>(i)]));
                            }
                        } else if (actin_turnover_enabled) {
                            throw std::runtime_error(
                                "Resume snapshot is missing actin renewal ages");
                        } else {
                            actin_last_turnover_step.assign(
                                static_cast<size_t>(actin.n), 0);
                        }
                    }

                    if (titin_mode == TitinMode::explicit_dynamic) {
                        const size_t titin_count = explicit_titins.size();
                        std::vector<int> titin_int;
                        std::vector<double> titin_double;
                        if (!load_fixed_matrix_int(
                                group_resume, "titin_state_int", titin_count,
                                5, titin_int) ||
                            !load_fixed_matrix_double(
                                group_resume, "titin_state_double", titin_count,
                                6, titin_double)) {
                            throw std::runtime_error(
                                "Resume snapshot is missing explicit titin molecular state");
                        }
                        for (size_t index = 0; index < titin_count; ++index) {
                            auto& titin = explicit_titins[index];
                            titin.bound = titin_int[5 * index] != 0;
                            titin.bound_actin = titin.bound
                                ? titin_int[5 * index + 1] : -1;
                            titin.last_bound_actin =
                                titin_int[5 * index + 2];
                            titin.attach_step = static_cast<RecoveryStep>(
                                std::max(0, titin_int[5 * index + 3]));
                            titin.last_detach_step = static_cast<RecoveryStep>(
                                std::max(0, titin_int[5 * index + 4]));
                            titin.distance = titin_double[6 * index];
                            titin.extension = titin_double[6 * index + 1];
                            titin.force = titin_double[6 * index + 2];
                            titin.attachment_distance =
                                titin_double[6 * index + 3];
                            titin.attachment_energy =
                                titin_double[6 * index + 4];
                            titin.attachment_force =
                                titin_double[6 * index + 5];
                            titin.candidate_actin = -1;
                            titin.candidate_distance = 0.0;
                            titin.candidate_on_rate = 0.0;
                        }
                        rebuild_explicit_titin_ownership();
                        // Association provenance is optional on resume: a
                        // snapshot written before schema 3 simply restarts the
                        // counters, which is visible in the series rather than
                        // silently wrong.
                        std::vector<int> association_counts;
                        if (load_fixed_matrix_int(
                                group_resume, "titin_association_counts", 1, 4,
                                association_counts) &&
                            association_counts.size() == 4) {
                            explicit_titin_supported_associations =
                                static_cast<std::uint64_t>(
                                    std::max(0, association_counts[0]));
                            explicit_titin_unsupported_associations =
                                static_cast<std::uint64_t>(
                                    std::max(0, association_counts[1]));
                            explicit_titin_rebind_associations =
                                static_cast<std::uint64_t>(
                                    std::max(0, association_counts[2]));
                            explicit_titin_unsupported_releases =
                                static_cast<std::uint64_t>(
                                    std::max(0, association_counts[3]));
                        }
                    }

                    if (k_mm > 0.0) {
                        const size_t bridge_rows = std::max<size_t>(
                            1, static_cast<size_t>(myosin.n) *
                                   static_cast<size_t>(
                                       myomesin_sites_per_filament) /
                                   2);
                        std::vector<int> bridge_count;
                        std::vector<int> bridge_int;
                        std::vector<double> bridge_double;
                        std::vector<int> next_id;
                        if (!load_fixed_vector_int(
                                group_resume, "myomesin_bridge_count", 1,
                                bridge_count) ||
                            !load_fixed_matrix_int(
                                group_resume, "myomesin_state_int",
                                bridge_rows, 5, bridge_int) ||
                            !load_fixed_matrix_double(
                                group_resume, "myomesin_state_double",
                                bridge_rows, 3, bridge_double) ||
                            !load_fixed_vector_int(
                                group_resume, "myomesin_next_id", 2,
                                next_id)) {
                            throw std::runtime_error(
                                "Resume snapshot is missing sparse myomesin state");
                        }
                        const size_t count = static_cast<size_t>(
                            std::max(0, bridge_count[0]));
                        if (count > bridge_rows) {
                            throw std::runtime_error(
                                "Resume snapshot has invalid myomesin count");
                        }
                        myomesin_bridges.clear();
                        std::vector<int> degree(
                            static_cast<size_t>(myosin.n), 0);
                        std::unordered_set<std::uint64_t> pairs;
                        for (size_t index = 0; index < count; ++index) {
                            MyomesinBridge bridge;
                            bridge.id =
                                (static_cast<std::uint64_t>(
                                     static_cast<std::uint32_t>(
                                         bridge_int[5 * index])) << 32) |
                                static_cast<std::uint32_t>(
                                    bridge_int[5 * index + 1]);
                            bridge.myosin_i = bridge_int[5 * index + 2];
                            bridge.myosin_j = bridge_int[5 * index + 3];
                            bridge.attach_step = static_cast<RecoveryStep>(
                                std::max(0, bridge_int[5 * index + 4]));
                            bridge.axial_offset = bridge_double[3 * index];
                            bridge.radial_distance =
                                bridge_double[3 * index + 1];
                            bridge.force = bridge_double[3 * index + 2];
                            if (bridge.id == 0 || bridge.myosin_i < 0 ||
                                bridge.myosin_i >= myosin.n ||
                                bridge.myosin_j <= bridge.myosin_i ||
                                bridge.myosin_j >= myosin.n) {
                                throw std::runtime_error(
                                    "Resume snapshot contains invalid myomesin bridge");
                            }
                            const auto key = counter_rng::canonical_pair_key(
                                bridge.myosin_i, bridge.myosin_j);
                            if (!pairs.insert(key).second ||
                                ++degree[static_cast<size_t>(bridge.myosin_i)] >
                                    myomesin_sites_per_filament ||
                                ++degree[static_cast<size_t>(bridge.myosin_j)] >
                                    myomesin_sites_per_filament) {
                                throw std::runtime_error(
                                    "Resume snapshot violates myomesin valency");
                            }
                            myomesin_bridges.push_back(bridge);
                        }
                        myomesin_next_id =
                            (static_cast<std::uint64_t>(
                                 static_cast<std::uint32_t>(next_id[0])) << 32) |
                            static_cast<std::uint32_t>(next_id[1]);
                    } else {
                        myomesin_bridges.clear();
                    }

                    if (restore_recovery_state_from_fixed_pairs(
                            group_resume,
                            "aa_recovery_count",
                            "aa_recovery_pairs",
                            "aa_recovery_until_values")) {
                        // loaded sparse recovery state
                    } else if (load_fixed_matrix_int(group_resume, "actin_recovery_until",
                                                     static_cast<size_t>(actin.n),
                                                     static_cast<size_t>(actin.n), ints)) {
                        for (int i = 0; i < actin.n; ++i) {
                            for (int j = 0; j < actin.n; ++j) {
                                actin_recovery_until[i][j] = static_cast<RecoveryStep>(std::max(
                                    0, ints[static_cast<size_t>(i) * actin.n + j]));
                            }
                        }
                    }

                    if (rng != nullptr &&
                        load_fixed_vector_int(group_resume, "rng_main_state", gsl_rng_size(rng), ints)) {
                        auto* state_ptr = static_cast<unsigned char*>(gsl_rng_state(rng));
                        for (size_t idx = 0; idx < ints.size(); ++idx) {
                            state_ptr[idx] = static_cast<unsigned char>(std::clamp(ints[idx], 0, 255));
                        }
                    }
                    if (!rng_engines.empty() && rng_engines[0] != nullptr &&
                        load_fixed_vector_int(group_resume, "rng_thread_state_size", 1, ints)) {
                        const int file_state_size = ints[0];
                        if (load_fixed_vector_int(group_resume, "rng_thread_count", 1, ints)) {
                            const int file_thread_count = ints[0];
                            const size_t flat_width =
                                static_cast<size_t>(std::max(0, file_thread_count)) *
                                static_cast<size_t>(std::max(0, file_state_size));
                            if (load_fixed_vector_int(group_resume, "rng_thread_state", flat_width, ints)) {
                                const int restore_threads =
                                    std::min(static_cast<int>(rng_engines.size()), file_thread_count);
                                for (int t = 0; t < restore_threads; ++t) {
                                    if (rng_engines[t] == nullptr) {
                                        continue;
                                    }
                                    auto* local_ptr =
                                        static_cast<unsigned char*>(gsl_rng_state(rng_engines[t]));
                                    const size_t local_size = gsl_rng_size(rng_engines[t]);
                                    const size_t copy_size =
                                        std::min(local_size, static_cast<size_t>(std::max(0, file_state_size)));
                                    const size_t base =
                                        static_cast<size_t>(t) * static_cast<size_t>(std::max(0, file_state_size));
                                    for (size_t idx = 0; idx < copy_size; ++idx) {
                                        local_ptr[idx] =
                                            static_cast<unsigned char>(std::clamp(ints[base + idx], 0, 255));
                                    }
                                }
                                const int local_thread_count = static_cast<int>(rng_engines.size());
                                if (local_thread_count > file_thread_count) {
                                    for (int t = file_thread_count; t < local_thread_count; ++t) {
                                        if (rng_engines[t] == nullptr) {
                                            continue;
                                        }
                                        gsl_rng_set(
                                            rng_engines[t],
                                            static_cast<unsigned long>(initial_seed) +
                                                static_cast<unsigned long>(t) +
                                                static_cast<unsigned long>(current_step));
                                    }
                                }
                            }
                        }
                    }
                    
                    if (load_fixed_vector_double(group_resume, "neighbor_last_actin_x", static_cast<size_t>(actin.n), neighbor_last_actin_x) &&
                        load_fixed_vector_double(group_resume, "neighbor_last_actin_y", static_cast<size_t>(actin.n), neighbor_last_actin_y) &&
                        load_fixed_vector_double(group_resume, "neighbor_last_actin_z", static_cast<size_t>(actin.n), neighbor_last_actin_z) &&
                        load_fixed_vector_double(group_resume, "neighbor_last_myosin_x", static_cast<size_t>(myosin.n), neighbor_last_myosin_x) &&
                        load_fixed_vector_double(group_resume, "neighbor_last_myosin_y", static_cast<size_t>(myosin.n), neighbor_last_myosin_y) &&
                        load_fixed_vector_double(group_resume, "neighbor_last_myosin_z", static_cast<size_t>(myosin.n), neighbor_last_myosin_z)) {
                        neighbor_list.set_species_positions(neighbor_last_actin_x, neighbor_last_actin_y, neighbor_last_actin_z,
                                                            neighbor_last_myosin_x, neighbor_last_myosin_y, neighbor_last_myosin_z);
                        neighbor_list.rebuild_neighbor_list();
                        neighbor_list.set_species_positions(actin.center_x, actin.center_y, actin.center_z,
                                                            myosin.center_x, myosin.center_y, myosin.center_z);
                        restored_neighbor_cache = true;
                    } else {
                        restored_neighbor_cache = false;
                    }

                    if (!load_any_vector_double(group_resume, "cb_breakage_pending", cb_breakage_events)) {
                        cb_breakage_events.clear();
                    }
                    if (!load_any_vector_double(group_resume, "cb_limit_pending", cb_limit_events)) {
                        cb_limit_events.clear();
                    }
                    if (!load_any_vector_double(group_resume, "aa_completed_lifetimes_pending",
                                                aa_completed_lifetimes)) {
                        aa_completed_lifetimes.clear();
                    }
                    if (!load_any_vector_double(group_resume, "am_transition_events_pending",
                                                am_transition_events)) {
                        am_transition_events.clear();
                    }
                    if (!load_any_vector_double(
                            group_resume, "am_forced_rupture_events_pending",
                            am_forced_rupture_events)) {
                        am_forced_rupture_events.clear();
                    }
                    if (!load_any_vector_double(
                            group_resume, "am_forced_rupture_events_v2_pending",
                            am_forced_rupture_events_v2)) {
                        am_forced_rupture_events_v2.clear();
                    }
                    if (!load_any_vector_double(
                            group_resume, "am_thick_transition_events_pending",
                            am_thick_transition_events)) {
                        am_thick_transition_events.clear();
                    }
                    if (!load_any_vector_double(
                            group_resume, "actin_turnover_events_pending",
                            actin_turnover_events)) {
                        actin_turnover_events.clear();
                    }
                    if (!load_any_vector_double(group_resume, "aa_state_events_pending",
                                                aa_state_events)) {
                        aa_state_events.clear();
                    }
                    if (!load_any_vector_double(group_resume, "aa_load_events_pending",
                                                aa_load_events)) {
                        aa_load_events.clear();
                    }
                    if (!load_any_vector_double(
                            group_resume, "aa_local_link_events_pending",
                            aa_local_link_events)) {
                        aa_local_link_events.clear();
                    }
                    if (!load_any_vector_double(
                            group_resume, "aa_local_link_events_v2_pending",
                            aa_local_link_events_v2)) {
                        aa_local_link_events_v2.clear();
                    }
                    if (!load_any_vector_double(
                            group_resume, "aa_local_link_events_v3_pending",
                            aa_local_link_events_v3)) {
                        aa_local_link_events_v3.clear();
                    }
                    if (!load_any_vector_double(
                            group_resume,
                            "aa_tip_accessibility_events_pending",
                            aa_tip_accessibility_events)) {
                        aa_tip_accessibility_events.clear();
                    }
                    if (!load_any_vector_double(
                            group_resume, "aa_motor_pair_events_pending",
                            aa_motor_pair_events)) {
                        aa_motor_pair_events.clear();
                    }
                    if (!load_any_vector_double(
                            group_resume, "aa_arm_events_pending",
                            aa_arm_events)) {
                        aa_arm_events.clear();
                    }
                    if (!load_any_vector_double(
                            group_resume, "explicit_titin_events_pending",
                            explicit_titin_events)) {
                        explicit_titin_events.clear();
                    }
                    if (!load_any_vector_double(
                            group_resume, "myomesin_events_pending",
                            myomesin_events)) {
                        myomesin_events.clear();
                    }
                    if (resume_step_value < current_step) {
                        throw H5::Exception("Sarcomere::load_state", "Resume snapshot older than committed /state");
                    }
                    loaded_resume_snapshot = true;
                } catch (H5::Exception&) {
                    loaded_resume_snapshot = loaded_state_group;
                }
            } else {
                loaded_resume_snapshot = loaded_state_group;
            }
        }

    } catch (H5::Exception&) {
        loaded_state_group = false;
    }

    if (!loaded_state_group) {
        throw std::runtime_error("Could not load /state group for resume. Refusing approximate fallback.");
    }

    if (!restored_neighbor_cache) {
        // Legacy checkpoint fallback: rebuild from current positions.
        neighbor_list.set_species_positions(actin.center_x, actin.center_y, actin.center_z,
                                            myosin.center_x, myosin.center_y, myosin.center_z);
        neighbor_list.rebuild_neighbor_list();
    }

    for (int i = 0; i < actin.n; ++i) {
        myosinIndicesPerActin.deleteAllConnections(i);
    }
    for (int m = 0; m < myosin.n; ++m) {
        actinIndicesPerMyosin.deleteAllConnections(m);
        for (auto& temp_conn : actinIndicesPerMyosin_temp) {
            temp_conn.deleteAllConnections(m);
        }
    }
    for (int i = 0; i < actin.n; ++i) {
        for (int m = 0; m < myosin.n; ++m) {
            if (am_bonds[i][m] != 1) {
                continue;
            }
            myosinIndicesPerActin.addConnection(i, m);
            actinIndicesPerMyosin.addConnection(m, i);
        }
    }
    has_myosin_bond_pairs = false;

    _rebuild_active_pair_lists_from_dense();
    return target_frame;
}

