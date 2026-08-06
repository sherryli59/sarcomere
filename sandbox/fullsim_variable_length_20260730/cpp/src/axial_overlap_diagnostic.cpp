#include "axial_cell.h"
#include "langevin.h"
#include "rod_mobility.h"

#include <gsl/gsl_rng.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

struct HeadSummary {
    int weak = 0;
    int strong = 0;
    double weak_extension = 0.0;
    double strong_extension = 0.0;
};

struct AAPairGeometry {
    double signed_axial = 0.0;
    double transverse = 0.0;
    double distance = 0.0;
    int status1 = 0;
    int status2 = 0;
};

struct HexGeometry {
    double am_radius_mean = 0.0;
    double am_radius_cv = 0.0;
    double triplet_chord_mean = 0.0;
    double triplet_chord_cv = 0.0;
    double triplet_shape_error = 0.0;
    double z_radius_mean = 0.0;
    double z_radius_cv = 0.0;
    double z_edge_mean = 0.0;
    double z_edge_cv = 0.0;
    double z_angle_max_error_deg = 0.0;
};

std::pair<double, double> mean_cv(const std::vector<double>& values) {
    if (values.empty()) return {0.0, 0.0};
    double mean = 0.0;
    for (double value : values) mean += value;
    mean /= static_cast<double>(values.size());
    double variance = 0.0;
    for (double value : values) {
        const double delta = value - mean;
        variance += delta * delta;
    }
    variance /= static_cast<double>(values.size());
    return {mean, mean > 0.0 ? std::sqrt(variance) / mean : 0.0};
}

HexGeometry hex_geometry(const Sarcomere& model) {
    HexGeometry result;
    if (model.actin.n != 12 || model.myosin.n != 2) return result;
    const std::array<std::array<int, 5>, 4> populations{{
        {{0, 0, 9, 10, 11}}, {{0, 1, 0, 1, 2}},
        {{1, 0, 3, 4, 5}}, {{1, 1, 6, 7, 8}},
    }};
    std::vector<double> radii;
    std::vector<double> chords;
    double shape_error = 0.0;
    for (const auto& population : populations) {
        const int myosin = population[0];
        vec axis = model.myosin.direction[myosin];
        axis.normalize();
        double triplet_radius = 0.0;
        for (int offset = 0; offset < 3; ++offset) {
            const int actin = population[2 + offset];
            vec displacement =
                static_cast<vec>(model.actin.center[actin]) -
                static_cast<vec>(model.myosin.center[myosin]);
            displacement.pbc_wrap(model.box, model.is_periodic);
            const vec radial = displacement - displacement.dot(axis) * axis;
            radii.push_back(radial.norm());
            triplet_radius += radial.norm();
        }
        triplet_radius /= 3.0;
        for (int offset = 0; offset < 3; ++offset) {
            const int first = population[2 + offset];
            const int second = population[2 + (offset + 1) % 3];
            vec displacement =
                static_cast<vec>(model.actin.center[second]) -
                static_cast<vec>(model.actin.center[first]);
            displacement.pbc_wrap(model.box, model.is_periodic);
            const vec transverse =
                displacement - displacement.dot(axis) * axis;
            const double chord = transverse.norm();
            chords.push_back(chord);
            if (triplet_radius > 0.0) {
                shape_error += std::abs(
                    chord / (std::sqrt(3.0) * triplet_radius) - 1.0);
            }
        }
    }
    const auto [radius_mean, radius_cv] = mean_cv(radii);
    const auto [chord_mean, chord_cv] = mean_cv(chords);
    result.am_radius_mean = radius_mean;
    result.am_radius_cv = radius_cv;
    result.triplet_chord_mean = chord_mean;
    result.triplet_chord_cv = chord_cv;
    result.triplet_shape_error = shape_error / 12.0;

    std::vector<double> z_radii;
    std::vector<double> z_edges;
    double max_angle_error = 0.0;
    for (int ring = 0; ring < 2; ++ring) {
        const int first = 6 * ring;
        double centroid_y = 0.0;
        double centroid_z = 0.0;
        for (int offset = 0; offset < 6; ++offset) {
            const vec point = model.actin.left_end[first + offset];
            centroid_y += point.y;
            centroid_z += point.z;
        }
        centroid_y /= 6.0;
        centroid_z /= 6.0;
        std::array<std::pair<double, int>, 6> angles{};
        for (int offset = 0; offset < 6; ++offset) {
            const vec point = model.actin.left_end[first + offset];
            const double dy = point.y - centroid_y;
            const double dz = point.z - centroid_z;
            z_radii.push_back(std::hypot(dy, dz));
            double angle = std::atan2(dz, dy);
            if (angle < 0.0) angle += 2.0 * M_PI;
            angles[static_cast<size_t>(offset)] = {angle, first + offset};
        }
        std::sort(angles.begin(), angles.end());
        for (int offset = 0; offset < 6; ++offset) {
            const auto [angle, actin] = angles[static_cast<size_t>(offset)];
            const auto [next_angle_raw, next_actin] =
                angles[static_cast<size_t>((offset + 1) % 6)];
            const double next_angle = offset == 5
                ? next_angle_raw + 2.0 * M_PI : next_angle_raw;
            max_angle_error = std::max(
                max_angle_error,
                std::abs(next_angle - angle - M_PI / 3.0));
            const vec first_point = model.actin.left_end[actin];
            const vec second_point = model.actin.left_end[next_actin];
            z_edges.push_back(std::hypot(
                second_point.y - first_point.y,
                second_point.z - first_point.z));
        }
    }
    const auto [z_radius_mean, z_radius_cv] = mean_cv(z_radii);
    const auto [z_edge_mean, z_edge_cv] = mean_cv(z_edges);
    result.z_radius_mean = z_radius_mean;
    result.z_radius_cv = z_radius_cv;
    result.z_edge_mean = z_edge_mean;
    result.z_edge_cv = z_edge_cv;
    result.z_angle_max_error_deg = max_angle_error * 180.0 / M_PI;
    return result;
}

AAPairGeometry aa_pair_geometry(
    const Sarcomere& model, int first, int second) {
    AAPairGeometry result;
    vec displacement =
        static_cast<vec>(model.actin.left_end[second]) -
        static_cast<vec>(model.actin.left_end[first]);
    displacement.pbc_wrap(model.box, model.is_periodic);
    result.signed_axial = displacement.x;
    result.transverse = std::hypot(displacement.y, displacement.z);
    result.distance = displacement.norm();
    for (const auto& link : model.aa_links) {
        if (!link.active) continue;
        const bool matches =
            (link.actin_i == first && link.actin_j == second) ||
            (link.actin_i == second && link.actin_j == first);
        if (!matches) continue;
        if (link.status == 2) ++result.status2;
        else ++result.status1;
    }
    return result;
}

HeadSummary summarize_heads(const Sarcomere& model) {
    HeadSummary result;
    for (int m = 0; m < model.myosin.n; ++m) {
        for (int half = 0; half < 2; ++half) {
            for (const auto& head :
                 model.am_heads_by_myosin[static_cast<size_t>(m)]
                                          [static_cast<size_t>(half)]) {
                if (head.bound_actin < 0 ||
                    head.bound_actin >= model.actin.n) continue;
                const auto state =
                    static_cast<am_motor_kinetics::State>(head.state);
                if (state != am_motor_kinetics::State::weak &&
                    state != am_motor_kinetics::State::strong) continue;
                const bool strong =
                    state == am_motor_kinetics::State::strong;
                const auto mechanics = compute_am_two_spring_head(
                    model.actin, model.myosin, head.bound_actin, m,
                    model.box, head.target_axial, head.anchor_axial,
                    head.binding_azimuth,
                    model.am_spatial_heads ? model.myosin.radius : 0.0,
                    model.am_spatial_heads
                        ? &model.myosin_material_e1[static_cast<size_t>(m)]
                        : nullptr,
                    model.am_spatial_heads
                        ? &model.myosin_material_e2[static_cast<size_t>(m)]
                        : nullptr,
                    model.actin_radius + model.am_surface_gap,
                    0.0, model.am_tether_slack,
                    model.am_motor_k *
                        (strong ? 1.0 : model.am_weak_stiffness_ratio),
                    head.stroke_progress, model.am_stall_force, strong,
                    model.am_bidirectional_crossbridge);
                if (strong) {
                    ++result.strong;
                    result.strong_extension +=
                        mechanics.longitudinal_extension;
                } else {
                    ++result.weak;
                    result.weak_extension +=
                        mechanics.longitudinal_extension;
                }
            }
        }
    }
    if (result.weak > 0) result.weak_extension /= result.weak;
    if (result.strong > 0) result.strong_extension /= result.strong;
    return result;
}

double mean_thin_activation(const Sarcomere& model) {
    if (model.am_thin_activation.empty()) return 0.0;
    double total = 0.0;
    for (double value : model.am_thin_activation) total += value;
    return total /
        static_cast<double>(model.am_thin_activation.size());
}

std::array<int, 2> link_counts(const Sarcomere& model) {
    std::array<int, 2> counts{{0, 0}};
    for (const auto& link : model.aa_links) {
        if (link.active) ++counts[link.status == 2 ? 1 : 0];
    }
    return counts;
}

double mean_abs_tilt(const Filament& filament) {
    double total = 0.0;
    for (int i = 0; i < filament.n; ++i) {
        total += std::acos(std::clamp(
            std::abs(filament.direction[i].x), 0.0, 1.0));
    }
    return filament.n > 0 ? total / filament.n : 0.0;
}

double mean_feature(const Filament& filament, const std::string& name) {
    const auto found = filament.custom_features.find(name);
    if (found == filament.custom_features.end() || found->second.empty()) {
        return 0.0;
    }
    double total = 0.0;
    for (double value : found->second) total += value;
    return total / static_cast<double>(found->second.size());
}

void write_header(std::ofstream& out, int n_actins, int n_myosins) {
    out << "step,time,Lx,aa_overlap,aa_linked_overlap,am_overlap,"
           "aa01_axial,aa01_transverse,aa01_distance,aa01_status1,aa01_status2,"
           "aa23_axial,aa23_transverse,aa23_distance,aa23_status1,aa23_status2,"
           "thin_q,q_am,q_aa,q_titin,q_register,q_midpoint,q_external,q_total,"
           "dLx_dt,weak_heads,strong_heads,weak_extension,strong_extension,"
           "titin_bound,titin_candidates,titin_tension,titin_compressed,"
           "titin_released,titin_newly_bound,"
           "titin_mean_distance,titin_total_force,"
           "titin_attachment_energy_total,titin_attachment_energy_max,"
           "titin_attachment_force_max,"
           "aa_status1,aa_status2,actin_tilt,myosin_tilt,"
           "am_radius_mean,am_radius_cv,triplet_chord_mean,"
           "triplet_chord_cv,triplet_shape_error,z_radius_mean,z_radius_cv,"
           "z_edge_mean,z_edge_cv,z_angle_max_error_deg,"
           "actin_force_cap,actin_torque_cap,actin_rotation_cap,"
           "myosin_force_cap,myosin_torque_cap,myosin_rotation_cap";
    for (int i = 0; i < n_actins; ++i) {
        out << ",a" << i << "_x,a" << i << "_y,a" << i << "_z"
            << ",a" << i << "_ux,a" << i << "_uy,a" << i << "_uz"
            << ",a" << i << "_fx,a" << i << "_fy,a" << i << "_fz"
            << ",a" << i << "_tx,a" << i << "_ty,a" << i << "_tz";
    }
    for (int i = 0; i < n_myosins; ++i) {
        out << ",m" << i << "_x,m" << i << "_y,m" << i << "_z"
            << ",m" << i << "_ux,m" << i << "_uy,m" << i << "_uz"
            << ",m" << i << "_fx,m" << i << "_fy,m" << i << "_fz"
            << ",m" << i << "_tx,m" << i << "_ty,m" << i << "_tz";
    }
    out << '\n';
}

void write_frame(std::ofstream& out, const Sarcomere& model,
                 const axial_cell::Diagnostics& diag, int step,
                 double dt) {
    const auto heads = summarize_heads(model);
    const auto links = link_counts(model);
    const auto aa01 = aa_pair_geometry(model, 0, 1);
    const auto aa23 = aa_pair_geometry(model, 2, 3);
    const auto hex = hex_geometry(model);
    out << step << ',' << step * dt << ',' << model.box[0] << ','
        << axial_cell::mean_aa_overlap(model) << ','
        << axial_cell::mean_linked_aa_overlap(model) << ','
        << axial_cell::mean_am_overlap(model) << ','
        << aa01.signed_axial << ',' << aa01.transverse << ','
        << aa01.distance << ',' << aa01.status1 << ',' << aa01.status2 << ','
        << aa23.signed_axial << ',' << aa23.transverse << ','
        << aa23.distance << ',' << aa23.status1 << ',' << aa23.status2 << ','
        << mean_thin_activation(model) << ','
        << diag.q_am << ',' << diag.q_aa << ',' << diag.q_titin << ','
        << diag.q_register << ',' << diag.q_midpoint << ','
        << diag.q_external << ',' << diag.q_total << ','
        << diag.dlength_dt << ',' << heads.weak << ',' << heads.strong << ','
        << heads.weak_extension << ',' << heads.strong_extension << ','
        << diag.titin_bound_links << ','
        << diag.titin_candidate_links << ','
        << diag.titin_tension_links << ','
        << diag.titin_compressed_links << ','
        << diag.titin_released_links << ','
        << diag.titin_newly_bound_links << ','
        << diag.titin_mean_distance << ','
        << diag.titin_total_force << ','
        << diag.titin_attachment_energy_total << ','
        << diag.titin_attachment_energy_max << ','
        << diag.titin_attachment_force_max << ','
        << links[0] << ',' << links[1] << ','
        << mean_abs_tilt(model.actin) << ','
        << mean_abs_tilt(model.myosin) << ','
        << hex.am_radius_mean << ',' << hex.am_radius_cv << ','
        << hex.triplet_chord_mean << ',' << hex.triplet_chord_cv << ','
        << hex.triplet_shape_error << ','
        << hex.z_radius_mean << ',' << hex.z_radius_cv << ','
        << hex.z_edge_mean << ',' << hex.z_edge_cv << ','
        << hex.z_angle_max_error_deg << ','
        << mean_feature(model.actin, "cumulative_force_cap_events") << ','
        << mean_feature(model.actin, "cumulative_torque_cap_events") << ','
        << mean_feature(model.actin, "cumulative_rotation_cap_events") << ','
        << mean_feature(model.myosin, "cumulative_force_cap_events") << ','
        << mean_feature(model.myosin, "cumulative_torque_cap_events") << ','
        << mean_feature(model.myosin, "cumulative_rotation_cap_events");
    for (int i = 0; i < model.actin.n; ++i) {
        const vec c = model.actin.center[i];
        const vec u = model.actin.direction[i];
        const vec f = model.actin.force[i];
        const vec tau = model.actin.torque[i];
        out << ',' << c.x << ',' << c.y << ',' << c.z
            << ',' << u.x << ',' << u.y << ',' << u.z
            << ',' << f.x << ',' << f.y << ',' << f.z
            << ',' << tau.x << ',' << tau.y << ',' << tau.z;
    }
    for (int i = 0; i < model.myosin.n; ++i) {
        const vec c = model.myosin.center[i];
        const vec u = model.myosin.direction[i];
        const vec f = model.myosin.force[i];
        const vec tau = model.myosin.torque[i];
        out << ',' << c.x << ',' << c.y << ',' << c.z
            << ',' << u.x << ',' << u.y << ',' << u.z
            << ',' << f.x << ',' << f.y << ',' << f.z
            << ',' << tau.x << ',' << tau.y << ',' << tau.z;
    }
    out << '\n';
}

void write_link_header(std::ofstream& out) {
    out << "step,time,link_id,actin_i,actin_j,status,attach_step,s_i,s_j,"
           "first_x,first_y,first_z,second_x,second_y,second_z,"
           "connector_x,connector_y,connector_z,length,extension,"
           "force_first_x,force_first_y,force_first_z,"
           "torque_first_x,torque_first_y,torque_first_z,"
           "torque_second_x,torque_second_y,torque_second_z,"
           "direction_dot,registry_error\n";
}

void write_links(std::ofstream& out, const Sarcomere& model,
                 int step, double dt) {
    const double single_link_k =
        actn2_load_kinetics::fixed_single_link_stiffness(
            model.k_aa, model.aa_stiffness_reference_links);
    const double single_link_alignment_k =
        actn2_load_kinetics::fixed_single_link_stiffness(
            model.kappa_aa, model.aa_stiffness_reference_links);
    for (const auto& link : model.aa_links) {
        if (!link.active) continue;
        const auto mechanics = compute_local_aa_link(
            model.actin, link.actin_i, link.actin_j,
            link.s_i, link.s_j, model.box,
            single_link_k, model.aa_perpendicular_k,
            single_link_alignment_k, model.aa_optimal,
            model.aa_antiparallel_overlap,
            model.aa_link_register_ratio(link.status),
            model.aa_register_saturation_error,
            model.aa_axial_register_force_max);
        const vec first =
            static_cast<vec>(model.actin.center[link.actin_i]) +
            link.s_i * static_cast<vec>(
                model.actin.direction[link.actin_i]);
        vec second =
            static_cast<vec>(model.actin.center[link.actin_j]) +
            link.s_j * static_cast<vec>(
                model.actin.direction[link.actin_j]);
        vec connector = second - first;
        connector.pbc_wrap(model.box, model.is_periodic);
        second = first + connector;
        out << step << ',' << step * dt << ',' << link.id << ','
            << link.actin_i << ',' << link.actin_j << ','
            << static_cast<int>(link.status) << ',' << link.attach_step << ','
            << link.s_i << ',' << link.s_j << ','
            << first.x << ',' << first.y << ',' << first.z << ','
            << second.x << ',' << second.y << ',' << second.z << ','
            << mechanics.connector.x << ',' << mechanics.connector.y << ','
            << mechanics.connector.z << ',' << mechanics.length << ','
            << mechanics.length - model.aa_optimal << ','
            << mechanics.force_on_first.x << ','
            << mechanics.force_on_first.y << ','
            << mechanics.force_on_first.z << ','
            << mechanics.torque_on_first.x << ','
            << mechanics.torque_on_first.y << ','
            << mechanics.torque_on_first.z << ','
            << mechanics.torque_on_second.x << ','
            << mechanics.torque_on_second.y << ','
            << mechanics.torque_on_second.z << ','
            << mechanics.direction_dot << ','
            << mechanics.registry_error << '\n';
    }
}

void write_titin_header(std::ofstream& out) {
    out << "step,time,titin_id,myosin,half,track,anchor_azimuth,"
           "bound,bound_actin,candidate_actin,last_bound_actin,"
           "attach_step,last_detach_step,"
           "distance,extension,force,candidate_distance,candidate_energy,"
           "candidate_on_rate,attachment_distance,attachment_energy,"
           "attachment_force,thick_x,thick_y,thick_z,"
           "actin_x,actin_y,actin_z\n";
}

void write_titin(std::ofstream& out, const Sarcomere& model,
                 const axial_cell::Config& config, int step, double dt) {
    for (size_t index = 0; index < config.explicit_titins.size(); ++index) {
        const auto& titin = config.explicit_titins[index];
        const vec thick = axial_cell::explicit_titin_thick_anchor(model, titin);
        const vec thin = titin.bound && titin.bound_actin >= 0
            ? static_cast<vec>(model.actin.left_end[titin.bound_actin])
            : vec{0.0, 0.0, 0.0};
        out << step << ',' << step * dt << ',' << index << ','
            << titin.myosin << ',' << titin.half << ',' << titin.track << ','
            << titin.anchor_azimuth << ',' << titin.bound << ','
            << titin.bound_actin << ',' << titin.candidate_actin << ','
            << titin.last_bound_actin << ','
            << titin.attach_step << ',' << titin.last_detach_step << ','
            << titin.distance << ',' << titin.extension << ',' << titin.force
            << ',' << titin.candidate_distance
            << ',' << titin.candidate_energy
            << ',' << titin.candidate_on_rate
            << ',' << titin.attachment_distance
            << ',' << titin.attachment_energy
            << ',' << titin.attachment_force
            << ',' << thick.x << ',' << thick.y << ',' << thick.z
            << ',' << thin.x << ',' << thin.y << ',' << thin.z
            << '\n';
    }
}

void write_head_header(std::ofstream& out) {
    out << "step,time,myosin,half,head_id,state,bound_actin,"
           "anchor_x,anchor_y,anchor_z,binding_x,binding_y,binding_z,"
           "extension,stroke_progress,active_force\n";
}

void write_heads(std::ofstream& out, const Sarcomere& model,
                 int step, double dt) {
    for (int m = 0; m < model.myosin.n; ++m) {
        vec axis = model.myosin.direction[m];
        axis.normalize();
        const vec e1 =
            model.myosin_material_e1[static_cast<size_t>(m)];
        const vec e2 =
            model.myosin_material_e2[static_cast<size_t>(m)];
        for (int half = 0; half < 2; ++half) {
            const auto& heads =
                model.am_heads_by_myosin[static_cast<size_t>(m)]
                                         [static_cast<size_t>(half)];
            for (size_t head_id = 0; head_id < heads.size(); ++head_id) {
                const auto& head = heads[head_id];
                const auto state =
                    static_cast<am_motor_kinetics::State>(head.state);
                const bool bound =
                    head.bound_actin >= 0 &&
                    head.bound_actin < model.actin.n &&
                    (state == am_motor_kinetics::State::weak ||
                     state == am_motor_kinetics::State::strong);
                const double azimuth =
                    bound ? head.binding_azimuth : head.anchor_azimuth;
                vec anchor = static_cast<vec>(model.myosin.center[m]) +
                    head.anchor_axial * axis +
                    model.myosin.radius *
                        (std::cos(azimuth) * e1 +
                         std::sin(azimuth) * e2);
                vec binding{
                    std::numeric_limits<double>::quiet_NaN(),
                    std::numeric_limits<double>::quiet_NaN(),
                    std::numeric_limits<double>::quiet_NaN()};
                double extension = 0.0;
                double active_force = 0.0;
                if (bound) {
                    const bool strong =
                        state == am_motor_kinetics::State::strong;
                    const auto mechanics = compute_am_two_spring_head(
                        model.actin, model.myosin, head.bound_actin, m,
                        model.box, head.target_axial, head.anchor_axial,
                        head.binding_azimuth,
                        model.am_spatial_heads ? model.myosin.radius : 0.0,
                        &e1, &e2,
                        model.actin_radius + model.am_surface_gap,
                        0.0, model.am_tether_slack,
                        model.am_motor_k *
                            (strong ? 1.0 :
                                      model.am_weak_stiffness_ratio),
                        head.stroke_progress, model.am_stall_force, strong,
                        model.am_bidirectional_crossbridge);
                    anchor = mechanics.myosin_attachment;
                    binding = mechanics.actin_attachment;
                    extension = mechanics.longitudinal_extension;
                    active_force = mechanics.active_force;
                }
                out << step << ',' << step * dt << ',' << m << ','
                    << half << ',' << head_id << ','
                    << static_cast<int>(head.state) << ','
                    << head.bound_actin << ','
                    << anchor.x << ',' << anchor.y << ',' << anchor.z << ','
                    << binding.x << ',' << binding.y << ',' << binding.z << ','
                    << extension << ',' << head.stroke_progress << ','
                    << active_force << '\n';
            }
        }
    }
}

axial_cell::Config mechanism_config(const std::string& mechanism,
                                    Sarcomere& model) {
    axial_cell::Config config;
    config.mode = axial_cell::Mode::stress_clamp;
    config.mobility = 0.20;
    config.min_length = 3.4;
    config.max_length = 5.2;

    if (mechanism == "turnover") return config;
    if (mechanism == "turnover_fixed") {
        config.mode = axial_cell::Mode::fixed;
        return config;
    }
    if (mechanism == "titin_fixed") {
        config.mode = axial_cell::Mode::fixed;
        config.apply_internal_forces_when_fixed = true;
        config.titin_k = 220.0;
        config.titin_rest_length = 1.125;
        return config;
    }
    if (mechanism == "explicit_titin_fixed" ||
        mechanism == "explicit_titin_hex_fixed" ||
        mechanism == "explicit_titin") {
        config.mode = (mechanism == "explicit_titin_fixed" ||
                       mechanism == "explicit_titin_hex_fixed")
            ? axial_cell::Mode::fixed : axial_cell::Mode::stress_clamp;
        config.apply_internal_forces_when_fixed = true;
        config.explicit_titin = true;
        config.explicit_titin_per_half =
            mechanism == "explicit_titin_hex_fixed" ? 3 : 1;
        // The assembled thick-tip--barbed-end distance is approximately
        // 0.48--0.55 um.  The former 0.375-um slack length and 220-pN/um
        // stiffness loaded a newly captured tether by tens of pN.  Use a
        // low-force coarse-grained tether and charge strained capture its
        // elastic free-energy cost.
        config.explicit_titin_k = 40.0;
        config.explicit_titin_slack_length = 0.500;
        config.explicit_titin_capture_distance = 0.70;
        config.explicit_titin_k_on = 5000.0;
        config.explicit_titin_k_off = 0.5;
        config.explicit_titin_strain_gated_kinetics = true;
        config.explicit_titin_transition_theta = 0.5;
        axial_cell::initialize_explicit_titins(model, config);
        return config;
    }
    if (mechanism != "titin" &&
        mechanism != "register_titin" &&
        mechanism != "register_titin_consistent" &&
        mechanism != "register_titin_strong" &&
        mechanism != "opposing_restore_titin" &&
        mechanism != "opposing_assist_titin" &&
        mechanism != "bidirectional_titin") {
        throw std::invalid_argument("unknown mechanism: " + mechanism);
    }
    config.titin_k = 220.0;
    config.titin_rest_length = 1.125;
    if (mechanism == "register_titin" ||
        mechanism == "register_titin_consistent" ||
        mechanism == "register_titin_strong") {
        config.register_k = 180.0;
        config.register_target = 0.10;
        config.register_force_cap =
            mechanism == "register_titin_strong" ? 100.0 : 12.0;
        config.register_endpoint_torque =
            mechanism != "register_titin";
    } else if (mechanism == "opposing_restore_titin" ||
               mechanism == "opposing_assist_titin") {
        config.midpoint_feedback_k = 180.0;
        config.midpoint_target = 0.10;
        config.midpoint_force_cap = 12.0;
        config.midpoint_feedback_sign =
            mechanism == "opposing_restore_titin" ? 1.0 : -1.0;
    } else if (mechanism == "bidirectional_titin") {
        model.am_bidirectional_crossbridge = true;
    }
    return config;
}

void write_event_summary(const Sarcomere& model,
                         const std::string& output_prefix) {
    constexpr size_t transition_width = 11;
    constexpr size_t forced_width = 19;
    std::array<size_t, 7> reasons{{0, 0, 0, 0, 0, 0, 0}};
    for (size_t offset = 0;
         offset + transition_width <= model.am_transition_events.size();
         offset += transition_width) {
        const int reason = static_cast<int>(
            std::llround(model.am_transition_events[offset + 5]));
        if (reason >= 0 && reason < static_cast<int>(reasons.size())) {
            ++reasons[static_cast<size_t>(reason)];
        }
    }
    std::array<size_t, 7> failure_bits{{0, 0, 0, 0, 0, 0, 0}};
    for (size_t offset = 0;
         offset + forced_width <= model.am_forced_rupture_events_v2.size();
         offset += forced_width) {
        const int mask = static_cast<int>(
            std::llround(model.am_forced_rupture_events_v2[offset + 6]));
        for (int bit = 0; bit < static_cast<int>(failure_bits.size()); ++bit) {
            if ((mask & (1 << bit)) != 0) {
                ++failure_bits[static_cast<size_t>(bit)];
            }
        }
    }
    std::ofstream out(output_prefix + "_event_summary.csv");
    out << "attachment,weak_detachment,powerstroke,strong_detachment,"
           "geometric_rupture,neighbor_loss,capacity_limit,"
           "forced_not_neighbor,forced_segment_distance,"
           "forced_wrong_half,forced_no_productive_overlap,"
           "forced_association_gate,forced_unexpected_candidate,"
           "forced_head_reach\n";
    for (size_t index = 0; index < reasons.size(); ++index) {
        if (index > 0) out << ',';
        out << reasons[index];
    }
    for (size_t index = 0; index < failure_bits.size(); ++index) {
        out << ',' << failure_bits[index];
    }
    out << '\n';
}

void write_event_details(const Sarcomere& model,
                         const std::string& output_prefix,
                         int release_step_origin, double dt) {
    constexpr size_t transition_width = 11;
    std::ofstream transitions(output_prefix + "_events.csv");
    transitions
        << "actin,myosin,release_step,release_time,old_state,new_state,"
           "reason,activation,strain,active_force,previous_aa_status,head_id\n";
    for (size_t offset = 0;
         offset + transition_width <= model.am_transition_events.size();
         offset += transition_width) {
        const auto event_step = static_cast<long long>(std::llround(
            model.am_transition_events[offset + 2]));
        const auto release_step =
            event_step - static_cast<long long>(release_step_origin);
        transitions
            << model.am_transition_events[offset] << ','
            << model.am_transition_events[offset + 1] << ','
            << release_step << ',' << release_step * dt;
        for (size_t field = 3; field < transition_width; ++field) {
            transitions << ',' << model.am_transition_events[offset + field];
        }
        transitions << '\n';
    }

    constexpr size_t forced_width = 19;
    std::ofstream ruptures(output_prefix + "_ruptures.csv");
    ruptures
        << "actin,myosin,release_step,release_time,head_id,old_state,reason,"
           "failure_mask,segment_distance,rupture_distance,actin_myosin_dot,"
           "expected_half,head_half,productive_axial_overlap,binding_ratio,"
           "partial_ratio,head_connector_distance,head_surface_distance,"
           "capture_reach,retention_reach\n";
    for (size_t offset = 0;
         offset + forced_width <= model.am_forced_rupture_events_v2.size();
         offset += forced_width) {
        const auto event_step = static_cast<long long>(std::llround(
            model.am_forced_rupture_events_v2[offset + 2]));
        const auto release_step =
            event_step - static_cast<long long>(release_step_origin);
        ruptures
            << model.am_forced_rupture_events_v2[offset] << ','
            << model.am_forced_rupture_events_v2[offset + 1] << ','
            << release_step << ',' << release_step * dt;
        for (size_t field = 3; field < forced_width; ++field) {
            ruptures << ','
                     << model.am_forced_rupture_events_v2[offset + field];
        }
        ruptures << '\n';
    }

    constexpr size_t aa_width = 21;
    std::ofstream aa_events(output_prefix + "_aa_events.csv");
    aa_events
        << "link_id,actin_i,actin_j,release_step,release_time,event_code,"
           "old_status,new_status,effective_load,lifetime,s_i,s_j,"
           "arm_i,arm_j,bilateral_load,productive_support,status2_gate,"
           "single_link_tension,raw_spring_arm_i,raw_spring_arm_j,"
           "allocated_motor_i,allocated_motor_j\n";
    for (size_t offset = 0;
         offset + aa_width <= model.aa_local_link_events_v2.size();
         offset += aa_width) {
        const auto event_step = static_cast<long long>(std::llround(
            model.aa_local_link_events_v2[offset + 3]));
        const auto release_step =
            event_step - static_cast<long long>(release_step_origin);
        aa_events
            << model.aa_local_link_events_v2[offset] << ','
            << model.aa_local_link_events_v2[offset + 1] << ','
            << model.aa_local_link_events_v2[offset + 2] << ','
            << release_step << ',' << release_step * dt;
        for (size_t field = 4; field < aa_width; ++field) {
            aa_events << ','
                      << model.aa_local_link_events_v2[offset + field];
        }
        aa_events << '\n';
    }
}

void reset_release_diagnostics(Sarcomere& model) {
    model.am_transition_events.clear();
    model.am_forced_rupture_events_v2.clear();
    model.aa_local_link_events.clear();
    model.aa_local_link_events_v2.clear();
    for (auto& values : model.am_transition_events_temp) values.clear();
    const auto reset_filament_counters = [](auto& filament) {
        for (const std::string& name : {
                 "cumulative_force_cap_events",
                 "cumulative_torque_cap_events",
                 "cumulative_rotation_cap_events"}) {
            auto found = filament.custom_features.find(name);
            if (found != filament.custom_features.end()) {
                std::fill(found->second.begin(), found->second.end(), 0.0);
            }
        }
    };
    reset_filament_counters(model.actin);
    reset_filament_counters(model.myosin);
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 8 || argc > 12) {
        std::cerr
            << "usage: axial_overlap_diagnostic MECHANISM INITIAL_LX "
               "INITIAL_AA SEED NSTEPS SAVE_EVERY OUTPUT_PREFIX "
               "[full|no_rotation|deterministic] [HEADS_PER_HALF] "
               "[PREEQUIL_STEPS] "
               "[baseline|no_head_torque|no_aa_mechanics|no_am_alignment|"
               "legacy_aa_status|aa_k10x|aa_k100x|"
               "aa_k100x_register_r025|aa_k100x_register_r05|"
               "aa_k100x_register_r1|aa_k100x_register_r1_l3|"
               "aa_k100x_register_all_l2|"
               "aa_k100x_register_all_l3|actn2_unseeded_control|"
               "actn2_relaxed_capture|actn2_relaxed_titin060|"
               "actn2_relaxed_titin070|"
               "actn2_relaxed_titin070_rebind080|"
               "actn2_relaxed_titin070_rebind080_off005|"
               "actn2_relaxed_titin075_rebind085|"
               "actn2_relaxed_titin075_rebind085_kon50k|"
               "actn2_relaxed_titin075_rebind085_kon50k_memory20|"
               "actn2_relaxed_multititin3]\n";
        return 2;
    }
    const std::string mechanism = argv[1];
    const double initial_lx = std::stod(argv[2]);
    const double initial_aa = std::stod(argv[3]);
    const int seed = std::stoi(argv[4]);
    const int nsteps = std::stoi(argv[5]);
    const int diagnostic_every = std::stoi(argv[6]);
    const std::string output_prefix = argv[7];
    const std::string noise_mode = argc >= 9 ? argv[8] : "full";
    const int heads_per_half = argc >= 10 ? std::stoi(argv[9]) : 35;
    const int preequil_steps = argc >= 11 ? std::stoi(argv[10]) : 0;
    const std::string force_variant = argc >= 12 ? argv[11] :
        (mechanism == "explicit_titin_hex_fixed" ? "actn2_tuned"
                                                   : "baseline");
    if (!(initial_lx > 0.0) || !(initial_aa > 0.0) ||
        seed < 0 || nsteps < 0 || diagnostic_every <= 0 ||
        heads_per_half <= 0 || preequil_steps < 0) {
        throw std::invalid_argument("invalid diagnostic argument");
    }
    if (noise_mode != "full" && noise_mode != "no_rotation" &&
        noise_mode != "deterministic") {
        throw std::invalid_argument("invalid noise mode: " + noise_mode);
    }
    if (force_variant != "baseline" &&
        force_variant != "no_head_torque" &&
        force_variant != "no_aa_mechanics" &&
        force_variant != "no_am_alignment" &&
        force_variant != "legacy_aa_status" &&
        force_variant != "aa_k10x" &&
        force_variant != "aa_k100x" &&
        force_variant != "aa_k100x_register_r025" &&
        force_variant != "aa_k100x_register_r05" &&
        force_variant != "aa_k100x_register_r1" &&
        force_variant != "aa_k100x_register_r2" &&
        force_variant != "aa_k100x_register_r4" &&
        force_variant != "aa_k100x_register_r8" &&
        force_variant != "actn2_tuned" &&
        force_variant != "actn2_unseeded_control" &&
        force_variant != "actn2_relaxed_capture" &&
        force_variant != "actn2_relaxed_titin060" &&
        force_variant != "actn2_relaxed_titin070" &&
        force_variant != "actn2_relaxed_titin070_rebind080" &&
        force_variant != "actn2_relaxed_titin070_rebind080_off005" &&
        force_variant != "actn2_relaxed_titin075_rebind085" &&
        force_variant != "actn2_relaxed_titin075_rebind085_kon50k" &&
        force_variant !=
            "actn2_relaxed_titin075_rebind085_kon50k_memory20" &&
        force_variant != "actn2_relaxed_multititin3" &&
        force_variant != "aa_k100x_register_r1_l3" &&
        force_variant != "aa_k100x_register_all_l2" &&
        force_variant != "aa_k100x_register_all_l3") {
        throw std::invalid_argument(
            "invalid force variant: " + force_variant);
    }

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, static_cast<unsigned long>(seed));

    const bool hex_topology = mechanism == "explicit_titin_hex_fixed";
    const bool relaxed_titin060 =
        force_variant == "actn2_relaxed_titin060";
    const bool relaxed_titin_off005 =
        force_variant == "actn2_relaxed_titin070_rebind080_off005";
    const bool relaxed_titin_memory20 = force_variant ==
        "actn2_relaxed_titin075_rebind085_kon50k_memory20";
    const bool relaxed_multititin3 =
        force_variant == "actn2_relaxed_multititin3";
    const bool relaxed_titin_kon50k =
        force_variant == "actn2_relaxed_titin075_rebind085_kon50k" ||
        relaxed_titin_memory20;
    const bool relaxed_titin075_rebind085 =
        force_variant == "actn2_relaxed_titin075_rebind085" ||
        relaxed_titin_kon50k;
    const bool relaxed_titin_rebind080 =
        force_variant == "actn2_relaxed_titin070_rebind080" ||
        relaxed_titin_off005;
    const bool relaxed_titin070 =
        force_variant == "actn2_relaxed_titin070" ||
        relaxed_titin_rebind080 || relaxed_titin075_rebind085;
    const bool relaxed_actn2 =
        force_variant == "actn2_relaxed_capture" || relaxed_titin060 ||
        relaxed_titin070 || relaxed_multititin3;
    const bool unseeded_actn2 = relaxed_actn2 ||
        force_variant == "actn2_unseeded_control";
    int n_actins = hex_topology ? 12 : 4;
    int n_myosins = 2;
    std::vector<double> box{initial_lx, 3.0, 3.0};
    double actin_length = 1.0;
    double myosin_length = 1.5;
    double myosin_radius = 0.025;
    double actin_radius = 0.015;
    double surface_gap = 0.015;
    double am_cutoff = 0.05;
    double am_optimal = 0.055;
    // A regular six-member ring has edge length equal to its circumradius.
    // Match AA and AM preferred spacings in the hexagonal diagnostic so the
    // closed ACTN2 ring is not geometrically frustrated at frame zero.
    double aa_optimal = (hex_topology || unseeded_actn2)
        ? am_optimal : 0.0506;
    // The relaxed status-1 diagnostic captures out to 100 nm while the
    // locked spring still minimizes at the 55 nm preferred centerline.
    double aa_cutoff = relaxed_actn2 ? 0.10 : aa_optimal + 0.02;
    double zero = 0.0;
    double aa_k_on = hex_topology ? 5000.0 : 1000.0;
    double aa_k_off = 1.0;
    double base_lifetime = hex_topology ? 1.0 : 0.1;
    double directional_lifetime = hex_topology ? 0.5 : 0.05;
    double lifetime_coeff = 0.4;
    double diffusion_ratio = 4.0;
    double k_aa = 10.0;
    double kappa_aa = 100.0;
    double dt = 1.0e-6;
    const std::string h5_file = output_prefix + ".h5";
    std::string mutable_h5_file = h5_file;
    int mutable_seed = seed;
    int initial_fixed_myosins = 0;
    bool directional = true;
    const double temperature = 310.0;
    double beta = rod_mobility::beta(temperature);
    const auto actin_mobility = rod_mobility::finite_cylinder(
        actin_length, 0.010, temperature, 0.001,
        10.0, 10.0, 160.0, 160.0);
    const auto myosin_mobility = rod_mobility::finite_cylinder(
        myosin_length, 0.035, temperature, 0.001,
        20.0, 20.0, 64.0, 64.0);
    const double max_actin_force =
        0.005 / (beta * std::max(
            actin_mobility.parallel, actin_mobility.perpendicular) * dt);
    const double max_myosin_force =
        0.005 / (beta * std::max(
            myosin_mobility.parallel, myosin_mobility.perpendicular) * dt);
    const double max_actin_torque =
        0.005 / (beta * actin_mobility.tumble * dt);
    const double max_myosin_torque =
        0.005 / (beta * myosin_mobility.tumble * dt);

    Sarcomere model(
        n_actins, n_myosins, box, actin_length, myosin_length,
        myosin_radius, am_cutoff, am_optimal, aa_cutoff, aa_optimal,
        aa_k_on, aa_k_off, base_lifetime, directional_lifetime,
        lifetime_coeff, diffusion_ratio, k_aa, kappa_aa,
        zero, zero, zero, zero, mutable_h5_file, rng, mutable_seed,
        initial_fixed_myosins, dt, zero, zero, zero, directional, 24,
        max_actin_force, max_myosin_force,
        max_actin_torque, max_myosin_torque);
    model.set_periodicity({true, false, false});
    const int aa_links_per_pair =
        (force_variant == "aa_k100x_register_r1_l3" ||
         force_variant == "aa_k100x_register_all_l3") ? 3 : 2;
    model.set_local_aa_model(
        true, aa_links_per_pair, 0.0192, 0.1, 50.0, 0, 3);
    if (hex_topology || unseeded_actn2) {
        model.set_aa_neighbor_identity(true, 0.02, 0.20);
        model.max_aa_bonds_per_actin = 2;
        model.aa_association_max_direction_dot = -0.5;
        model.aa_link_rupture_length = 0.12;
    }
    // The initializer reads this value to construct the requested displaced
    // starting overlap.  After initialization, the mechanical target is reset
    // to 100 nm below so the unseeded probes can anneal toward it.
    model.set_aa_axial_overlap(initial_aa, 0.0, 0.1, 0.15, 0.05);
    model.aa_barbed_zone = relaxed_actn2 ? 0.20 : 0.10;
    model.aa_status1_radial_k = relaxed_actn2 ? 75.0 : -1.0;
    model.set_am_radial_geometry(actin_radius, surface_gap, 1000.0);
    model.set_steric_stiffness(2000.0, 1000.0);
    // Match the production initialization path. The occupied-pair radial
    // tether uses am_preferred_centerline, which is configured by
    // set_aa_spacing_geometry rather than set_am_radial_geometry.
    const double am_preferred_centerline =
        myosin_radius + actin_radius + surface_gap;
    model.set_aa_spacing_geometry(
        aa_optimal / am_preferred_centerline,
        am_preferred_centerline, aa_cutoff - aa_optimal);
    model.set_am_kinetics(
        true, 100.0, 200.0, 200.0, 160.0, 1.0,
        1000.0, 0.005, 4.5, 5.0, 0.1, 0.1,
        0.0, 0.02, 0.2, 0.0, 0.0, 0.5);
    model.set_am_proximity_tether(
        3000.0, 0.002, true, 30.0, 30.0, 2.0);
    model.set_am_retention_kinetics(true, 0.1, 20.0, 6, 1.0);
    model.am_offcenter_torque = true;
    model.am_spatial_heads = true;
    model.set_am_head_count(heads_per_half);
    model.set_am_neighbor_head_cap(24);
    model.set_am_weak_binding_floor(0.1);
    if (hex_topology) {
        model.init_aa_overlap = initial_aa;
        // Near-contact but unassembled: no ACTN2 bonds, no seeded motor heads,
        // and no perfect transverse ring. All bonds must form dynamically.
        model.three_membered_lattice(false, 0.012, 0.025, 0.02);
    } else {
        model.minimal_catch_variant("balanced");
        if (unseeded_actn2) {
            // Remove the deterministic seed links: this diagnostic is meant
            // to test whether the wider weak gate can nucleate status 1.
            model.set_local_aa_model(
                true, aa_links_per_pair, 0.0192, 0.1, 50.0, 0, 3);
            model.set_aa_neighbor_identity(true, 0.02, 0.20);
            // Begin at a 90 nm A-A centerline separation: outside the 75 nm
            // control gate but inside the relaxed 100 nm status-1 gate. Shift
            // each associated myosin with its actins so A-M radial geometry is
            // unchanged and capture reach is the only intended perturbation.
            constexpr double half_probe_spacing = 0.045;
            model.actin.center[0].y = -half_probe_spacing;
            model.actin.center[2].y = -half_probe_spacing;
            model.actin.center[1].y = half_probe_spacing;
            model.actin.center[3].y = half_probe_spacing;
            model.myosin.center[0].y = -half_probe_spacing;
            model.myosin.center[1].y = half_probe_spacing;
            model.actin.update_endpoints();
            model.myosin.update_endpoints();
        }
    }
    // Match the production driver's default ACTN2 state classifier. Direct
    // Sarcomere construction otherwise retains the legacy library default,
    // which cannot mature a relaxed local link from opposed motor support.
    model.set_aa_status_kinetics(
        actn2_load_kinetics::StatusMode::productive_register, 1.0, 0.5);
    // Keep the preferred register independent of the starting geometry.
    model.aa_antiparallel_overlap = 0.10;
    if (force_variant == "no_head_torque") {
        model.am_offcenter_torque = false;
    } else if (force_variant == "no_aa_mechanics") {
        model.k_aa = 0.0;
        model.kappa_aa = 0.0;
        model.aa_perpendicular_k = 0.0;
    } else if (force_variant == "no_am_alignment") {
        model.am_alignment_k = 0.0;
    } else if (force_variant == "legacy_aa_status") {
        model.set_aa_status_kinetics(
            actn2_load_kinetics::StatusMode::legacy_motor_identity,
            1.0, 0.5);
    } else if (force_variant == "aa_k10x") {
        model.k_aa *= 10.0;
    } else if (force_variant == "aa_k100x") {
        model.k_aa *= 100.0;
    } else if (force_variant == "aa_k100x_register_r025") {
        model.k_aa *= 100.0;
        model.aa_axial_register_ratio = 0.25;
        model.aa_axial_register_force_max = 6.0;
    } else if (force_variant == "aa_k100x_register_r05") {
        model.k_aa *= 100.0;
        model.aa_axial_register_ratio = 0.5;
        model.aa_axial_register_force_max = 9.0;
    } else if (force_variant == "aa_k100x_register_r1" ||
               force_variant == "aa_k100x_register_r1_l3" ||
               force_variant == "aa_k100x_register_all_l2" ||
               force_variant == "aa_k100x_register_all_l3") {
        model.k_aa *= 100.0;
        model.aa_axial_register_ratio = 1.0;
        model.aa_axial_register_force_max = 12.0;
        model.aa_register_all_statuses =
            force_variant == "aa_k100x_register_all_l2" ||
            force_variant == "aa_k100x_register_all_l3";
    } else if (force_variant == "aa_k100x_register_r2" ||
               force_variant == "aa_k100x_register_r4" ||
               force_variant == "aa_k100x_register_r8" ||
               force_variant == "actn2_tuned" ||
               force_variant == "actn2_unseeded_control" ||
               relaxed_actn2) {
        model.k_aa *= 100.0;
        const double register_ratio =
            force_variant == "aa_k100x_register_r2" ? 2.0 :
            (force_variant == "actn2_tuned" ||
             force_variant == "actn2_unseeded_control" ||
             relaxed_actn2) ? 3.0 :
            force_variant == "aa_k100x_register_r4" ? 4.0 : 8.0;
        model.aa_axial_register_ratio = register_ratio;
        model.aa_axial_register_force_max = 12.0 * register_ratio;
    }

    model.set_physical_mobility(
        temperature, beta, actin_mobility, myosin_mobility);

    double actin_trans = 0.2;
    double actin_rot = 0.1;
    double myosin_trans = 0.05;
    double myosin_rot = 0.05;
    int h5_save_every = std::max(nsteps + 1, 1);
    bool resume = false;
    Langevin sim(
        model, beta, dt, actin_trans, actin_rot,
        myosin_trans, myosin_rot, h5_save_every, resume, true,
        0.005, 0.005, 0.005, 0.005);
    sim.set_rod_mobility(actin_mobility, myosin_mobility);
    sim.set_thermal_noise_scale(noise_mode == "deterministic" ? 0.0 : 1.0);
    if (noise_mode == "no_rotation") {
        sim.set_rotational_noise_scale(0.0);
    }
    auto cell_config = mechanism_config(mechanism, model);
    if (relaxed_titin060) {
        // At the deliberately shallow 160-nm A-A starting register, the
        // intended individual titin domains begin 0.578--0.584 um from their
        // actin barbed ends. A 0.500-um slack length assigns this encounter
        // more than 25 kBT of elastic cost and kinetically excludes one half
        // before it can be motor captured. A 0.600-um slack domain permits
        // binding in that assembled range while remaining force-free until it
        // is extended beyond 0.600 um.
        cell_config.explicit_titin_slack_length = 0.600;
    } else if (relaxed_titin070) {
        // Match the force-free binding-domain reach to the existing 700-nm
        // capture radius. A bound domain is still a local, reciprocal,
        // tension-only tether: it generates no force in compression and only
        // resists separation after its end-to-end distance exceeds 700 nm.
        cell_config.explicit_titin_slack_length = 0.700;
        if (relaxed_titin075_rebind085) {
            cell_config.explicit_titin_slack_length = 0.750;
            cell_config.explicit_titin_capture_distance = 0.850;
            if (relaxed_titin_kon50k) {
                cell_config.explicit_titin_k_on = 50000.0;
            }
            if (relaxed_titin_memory20) {
                cell_config.explicit_titin_k = 200.0;
                cell_config.explicit_titin_rebind_memory = 0.020;
                cell_config.explicit_titin_rebind_rate_multiplier = 100.0;
                cell_config.explicit_titin_transition_theta = 0.90;
            }
        }
        if (relaxed_titin_rebind080) {
            // Keep a finite search shell beyond the slack boundary so a
            // strained tether that dissociates at 0.70--0.80 um can rebind.
            // Detailed-balance strain gating still suppresses association by
            // exp(-beta U); this only removes the irreversible hard cutoff.
            cell_config.explicit_titin_capture_distance = 0.800;
        }
        if (relaxed_titin_off005) {
            // Effective avidity for the coarse-grained half-lattice tether.
            // Turnover remains stochastic and strain gated; the zero-load
            // mean lifetime is 20 s rather than 2 s.
            cell_config.explicit_titin_k_off = 0.05;
        }
    }
    if (relaxed_multititin3) {
        cell_config.explicit_titin_per_half = 3;
        cell_config.explicit_titin_max_per_actin = 3;
        cell_config.explicit_titin_k = 40.0;
        cell_config.explicit_titin_slack_length = 0.750;
        cell_config.explicit_titin_capture_distance = 0.850;
        cell_config.explicit_titin_k_on = 5000.0;
        cell_config.explicit_titin_k_off = 0.5;
        cell_config.explicit_titin_transition_theta = 1.0;
        cell_config.explicit_titin_rebind_memory = 0.020;
        cell_config.explicit_titin_rebind_rate_multiplier = 100.0;
        // mechanism_config initialized the default one-track topology. Rebuild
        // the explicit material anchors after changing the molecular count.
        axial_cell::initialize_explicit_titins(model, cell_config);
    }
    sim.set_axial_cell(cell_config);

    // Advance activation and all binding-state kinetics at a completely fixed
    // lattice geometry. No Langevin or axial-cell step is taken here, so
    // neither filament centers nor orientations can move.
    for (int step = 0; step < preequil_steps; ++step) {
        model.update_system();
    }
    const auto release_heads = summarize_heads(model);
    const double release_q = mean_thin_activation(model);
    reset_release_diagnostics(model);

    std::ofstream out(output_prefix + ".csv");
    if (!out) throw std::runtime_error("could not open diagnostic CSV");
    std::ofstream head_out(output_prefix + "_heads.csv");
    if (!head_out) throw std::runtime_error("could not open head CSV");
    std::ofstream link_out(output_prefix + "_links.csv");
    if (!link_out) throw std::runtime_error("could not open link CSV");
    std::ofstream titin_out(output_prefix + "_titin.csv");
    if (!titin_out) throw std::runtime_error("could not open titin CSV");
    out << std::setprecision(10);
    head_out << std::setprecision(10);
    link_out << std::setprecision(10);
    titin_out << std::setprecision(10);
    write_header(out, model.actin.n, model.myosin.n);
    write_head_header(head_out);
    write_link_header(link_out);
    write_titin_header(titin_out);
    axial_cell::Diagnostics initial_diag;
    initial_diag.length_before = model.box[0];
    initial_diag.length_after = model.box[0];
    write_frame(out, model, initial_diag, 0, dt);
    write_heads(head_out, model, 0, dt);
    write_links(link_out, model, 0, dt);
    write_titin(titin_out, model, sim.axial_cell_config, 0, dt);

    int fixed_myosins = 0;
    for (int step = 1; step <= nsteps; ++step) {
        model.update_system();
        sim.sample_step(dt, rng, fixed_myosins);
        if (step % diagnostic_every == 0 || step == nsteps) {
            write_frame(out, model, sim.last_axial_cell, step, dt);
            write_heads(head_out, model, step, dt);
            write_links(link_out, model, step, dt);
            write_titin(
                titin_out, model, sim.axial_cell_config, step, dt);
        }
    }
    out.close();
    head_out.close();
    link_out.close();
    titin_out.close();
    write_event_summary(model, output_prefix);
    write_event_details(
        model, output_prefix, preequil_steps, dt);
    std::cout << mechanism << " seed=" << seed
              << " heads_per_half=" << heads_per_half
              << " preequil_ms=" << preequil_steps * dt * 1.0e3
              << " force_variant=" << force_variant
              << " release_q=" << release_q
              << " release_weak=" << release_heads.weak
              << " release_strong=" << release_heads.strong
              << " initial_Lx=" << initial_lx
              << " initial_AA=" << initial_aa
              << " final_Lx=" << model.box[0]
              << " final_AA=" << axial_cell::mean_aa_overlap(model)
              << " final_linked_AA="
              << axial_cell::mean_linked_aa_overlap(model)
              << " final_AM=" << axial_cell::mean_am_overlap(model)
              << '\n';
    gsl_rng_free(rng);
    return 0;
}
