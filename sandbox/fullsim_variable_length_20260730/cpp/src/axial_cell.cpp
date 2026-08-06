#include "axial_cell.h"

#include "actn2_load_kinetics.h"
#include "counter_rng.h"
#include "geometry.h"
#include "interaction.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
#include <vector>

namespace axial_cell {
namespace {

constexpr std::array<std::array<int, 2>, 2> aa_pairs{{{{0, 1}}, {{2, 3}}}};
constexpr std::array<int, 4> owning_myosin{{0, 1, 0, 1}};
constexpr std::array<std::array<int, 2>, 12> hex_aa_pairs{{
    {{0, 3}}, {{0, 5}}, {{1, 4}}, {{1, 3}}, {{2, 5}}, {{2, 4}},
    {{6, 9}}, {{6, 11}}, {{7, 10}}, {{7, 9}}, {{8, 11}}, {{8, 10}},
}};
constexpr std::array<int, 12> hex_owning_myosin{{
    0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0,
}};

bool supported_topology(const Sarcomere& model) {
    return (model.actin.n == 4 || model.actin.n == 12) &&
           model.myosin.n == 2 &&
           model.box.size() >= 1 && model.box[0] > 0.0 &&
           model.is_periodic[0];
}

double capped(double value, double cap) {
    if (!(cap > 0.0)) return value;
    return std::clamp(value, -cap, cap);
}

double interaction_ramp(const Sarcomere& model) {
    if (model.interaction_ramp_steps == 0) return 1.0;
    if (model.current_step <= model.interaction_ramp_start_step) return 0.0;
    const double progress = static_cast<double>(
        model.current_step - model.interaction_ramp_start_step) /
        static_cast<double>(model.interaction_ramp_steps);
    return std::clamp(progress, 0.0, 1.0);
}

double am_cell_force(const Sarcomere& model) {
    if (!model.dynamic_am) return 0.0;
    const double length = model.box[0];
    const double ramp = interaction_ramp(model);
    double q = 0.0;
    for (int myosin_index = 0; myosin_index < model.myosin.n;
         ++myosin_index) {
        for (int half_index = 0; half_index < 2; ++half_index) {
            const auto& heads =
                model.am_heads_by_myosin[static_cast<size_t>(myosin_index)]
                                         [static_cast<size_t>(half_index)];
            for (const auto& head : heads) {
                const int actin_index = head.bound_actin;
                if (actin_index < 0 || actin_index >= model.actin.n) continue;
                const auto state =
                    static_cast<am_motor_kinetics::State>(head.state);
                if (state == am_motor_kinetics::State::detached ||
                    state == am_motor_kinetics::State::off_srx) {
                    continue;
                }
                const bool strong =
                    state == am_motor_kinetics::State::strong;
                const double axial_k = ramp * model.am_motor_k *
                    (strong ? 1.0 : model.am_weak_stiffness_ratio);
                const auto mechanics = compute_am_two_spring_head(
                    model.actin, model.myosin, actin_index, myosin_index,
                    model.box, head.target_axial, head.anchor_axial,
                    head.binding_azimuth,
                    model.am_spatial_heads ? model.myosin.radius : 0.0,
                    model.am_spatial_heads
                        ? &model.myosin_material_e1[
                              static_cast<size_t>(myosin_index)]
                        : nullptr,
                    model.am_spatial_heads
                        ? &model.myosin_material_e2[
                              static_cast<size_t>(myosin_index)]
                        : nullptr,
                    model.actin_radius + model.am_surface_gap,
                    0.0, model.am_tether_slack, axial_k,
                    head.stroke_progress, model.am_stall_force, strong,
                    model.am_bidirectional_crossbridge);
                q += pair_generalized_force(
                    mechanics.force_on_actin.x,
                    model.actin.center[actin_index].x,
                    model.myosin.center[myosin_index].x, length);
            }
        }

        // Pair-level radial containment is part of the force actually
        // integrated by Sarcomere::_calc_am_force_velocity().
        for (int actin_index = 0; actin_index < model.actin.n;
             ++actin_index) {
            int occupied = 0;
            for (const auto& half_heads :
                 model.am_heads_by_myosin[
                     static_cast<size_t>(myosin_index)]) {
                for (const auto& head : half_heads) {
                    if (head.bound_actin == actin_index &&
                        head.state != static_cast<Sarcomere::BondFlag>(
                            am_motor_kinetics::State::detached) &&
                        head.state != static_cast<Sarcomere::BondFlag>(
                            am_motor_kinetics::State::off_srx)) {
                        ++occupied;
                    }
                }
            }
            if (occupied == 0) continue;
            const double occupancy = std::min(
                1.0, static_cast<double>(occupied) /
                         static_cast<double>(
                             std::max(1, model.am_tether_saturation_heads)));
            const auto tether = compute_am_head_tether(
                model.actin, model.myosin, actin_index, myosin_index,
                model.box, 0.0, ramp * model.am_tether_k * occupancy,
                model.am_preferred_centerline, 0.0, 0.0,
                nullptr, nullptr, model.actin_radius, model.am_tether_slack);
            q += pair_generalized_force(
                tether.force_on_actin.x,
                model.actin.center[actin_index].x,
                model.myosin.center[myosin_index].x, length);
        }
    }
    return q;
}

double aa_cell_force(const Sarcomere& model) {
    if (!model.aa_local_multilink) return 0.0;
    const double length = model.box[0];
    const double ramp = interaction_ramp(model);
    const double link_k = ramp *
        actn2_load_kinetics::fixed_single_link_stiffness(
            model.k_aa, model.aa_stiffness_reference_links);
    const double alignment_k = ramp *
        actn2_load_kinetics::fixed_single_link_stiffness(
            model.kappa_aa, model.aa_stiffness_reference_links);
    double q = 0.0;
    for (const auto& link : model.aa_links) {
        if (!link.active) continue;
        const auto mechanics = compute_local_aa_link(
            model.actin, link.actin_i, link.actin_j,
            link.s_i, link.s_j, model.box,
            link_k, ramp * model.aa_perpendicular_k, alignment_k,
            model.aa_optimal, model.aa_antiparallel_overlap,
            model.aa_link_register_ratio(link.status),
            model.aa_register_saturation_error,
            model.aa_axial_register_force_max);
        q += pair_generalized_force(
            mechanics.force_on_first.x,
            model.actin.center[link.actin_i].x,
            model.actin.center[link.actin_j].x, length);
    }
    return q;
}

double add_titin(Sarcomere& model, const Config& config,
                 double& net_added_force) {
    if (!(config.titin_k > 0.0) ||
        !(config.titin_rest_length > 0.0)) return 0.0;
    const double length = model.box[0];
    double q = 0.0;
    for (int actin_index = 0; actin_index < 4; ++actin_index) {
        const int myosin_index =
            owning_myosin[static_cast<size_t>(actin_index)];
        vec displacement =
            static_cast<vec>(model.myosin.center[myosin_index]) -
            static_cast<vec>(model.actin.left_end[actin_index]);
        displacement.pbc_wrap(model.box, model.is_periodic);
        const double distance = displacement.norm();
        if (!(distance > 1e-12)) continue;
        const vec force_on_actin =
            (config.titin_k *
             (distance - config.titin_rest_length) / distance) *
            displacement;
        model.actin.force[actin_index] += force_on_actin;
        model.myosin.force[myosin_index] -= force_on_actin;
        net_added_force += force_on_actin.x - force_on_actin.x;
        q += pair_generalized_force(
            force_on_actin.x, model.actin.center[actin_index].x,
            model.myosin.center[myosin_index].x, length);
    }
    return q;
}

vec titin_thick_anchor(const Sarcomere& model,
                       const ExplicitTitin& titin) {
    vec axis = model.myosin.direction[titin.myosin];
    axis.normalize();
    vec e1 = model.myosin_material_e1[static_cast<size_t>(titin.myosin)];
    e1 -= axis * e1.dot(axis);
    if (e1.norm() <= 1e-12) {
        const vec reference = std::abs(axis.z) < 0.9
            ? vec{0.0, 0.0, 1.0} : vec{0.0, 1.0, 0.0};
        e1 = axis.cross(reference);
    }
    e1.normalize();
    vec e2 = axis.cross(e1);
    e2.normalize();
    const double half_sign = titin.half == 0 ? -1.0 : 1.0;
    const vec radial = std::cos(titin.anchor_azimuth) * e1 +
                       std::sin(titin.anchor_azimuth) * e2;
    return static_cast<vec>(model.myosin.center[titin.myosin]) +
           half_sign * 0.5 * model.myosin.length * axis +
           model.myosin.radius * radial;
}

bool titin_actin_is_owned(const Config& config, int actin_index,
                          size_t except_titin) {
    int occupancy = 0;
    for (size_t index = 0; index < config.explicit_titins.size(); ++index) {
        if (index == except_titin) continue;
        const auto& titin = config.explicit_titins[index];
        if (titin.bound && titin.bound_actin == actin_index) ++occupancy;
    }
    return occupancy >= std::max(1, config.explicit_titin_max_per_actin);
}

int best_explicit_titin_candidate(const Sarcomere& model,
                                  const Config& config,
                                  size_t titin_index,
                                  const vec& thick_anchor) {
    const auto& titin = config.explicit_titins[titin_index];
    if (config.explicit_titin_rebind_memory > 0.0 &&
        titin.last_bound_actin >= 0 &&
        titin.last_bound_actin < model.actin.n &&
        !titin_actin_is_owned(
            config, titin.last_bound_actin, titin_index)) {
        const double elapsed = model.dt * static_cast<double>(
            model.current_step >= titin.last_detach_step
                ? model.current_step - titin.last_detach_step : 0);
        if (elapsed <= config.explicit_titin_rebind_memory) {
            vec displacement = thick_anchor - static_cast<vec>(
                model.actin.left_end[titin.last_bound_actin]);
            displacement.pbc_wrap(model.box, model.is_periodic);
            if (displacement.norm() <=
                config.explicit_titin_capture_distance) {
                return titin.last_bound_actin;
            }
        }
    }
    vec myosin_axis = model.myosin.direction[titin.myosin];
    myosin_axis.normalize();
    int best_actin = -1;
    double best_distance = config.explicit_titin_capture_distance;
    for (int actin_index = 0; actin_index < model.actin.n; ++actin_index) {
        if (titin_actin_is_owned(config, actin_index, titin_index)) continue;
        vec from_myosin =
            static_cast<vec>(model.actin.left_end[actin_index]) -
            static_cast<vec>(model.myosin.center[titin.myosin]);
        from_myosin.pbc_wrap(model.box, model.is_periodic);
        const double half_projection = from_myosin.dot(myosin_axis);
        if ((titin.half == 0 && half_projection >= 0.0) ||
            (titin.half == 1 && half_projection <= 0.0)) {
            continue;
        }
        vec displacement = thick_anchor -
            static_cast<vec>(model.actin.left_end[actin_index]);
        displacement.pbc_wrap(model.box, model.is_periodic);
        vec actin_axis = model.actin.direction[actin_index];
        actin_axis.normalize();
        // The actin pointed direction must lead from its barbed-end domain
        // toward this thick-filament half.
        if (actin_axis.dot(displacement) <= 0.0) continue;
        const double distance = displacement.norm();
        if (distance > best_distance) continue;
        if (best_actin < 0 || distance < best_distance ||
            (distance == best_distance && actin_index < best_actin)) {
            best_actin = actin_index;
            best_distance = distance;
        }
    }
    return best_actin;
}

double add_explicit_titin(Sarcomere& model, Config& config,
                          Diagnostics& diagnostics) {
    if (!config.explicit_titin) return 0.0;
    if (!(config.explicit_titin_k > 0.0) ||
        !(config.explicit_titin_slack_length > 0.0) ||
        !(config.explicit_titin_capture_distance > 0.0) ||
        config.explicit_titin_max_per_actin <= 0 ||
        config.explicit_titin_k_on < 0.0 ||
        config.explicit_titin_k_off < 0.0 ||
        config.explicit_titin_transition_theta < 0.0 ||
        config.explicit_titin_transition_theta > 1.0 ||
        !std::isfinite(config.explicit_titin_rebind_memory) ||
        config.explicit_titin_rebind_memory < 0.0 ||
        !std::isfinite(config.explicit_titin_rebind_rate_multiplier) ||
        config.explicit_titin_rebind_rate_multiplier < 0.0 ||
        (config.explicit_titin_strain_gated_kinetics &&
         !(model.thermal_beta > 0.0))) {
        throw std::invalid_argument("invalid explicit titin configuration");
    }
    if (config.explicit_titins.empty()) {
        initialize_explicit_titins(model, config);
    }
    double q = 0.0;
    for (size_t index = 0; index < config.explicit_titins.size(); ++index) {
        auto& titin = config.explicit_titins[index];
        const vec thick_anchor = titin_thick_anchor(model, titin);
        if (!titin.bound) {
            titin.distance = 0.0;
            titin.extension = 0.0;
            titin.force = 0.0;
            titin.candidate_distance = 0.0;
            titin.candidate_energy = 0.0;
            titin.candidate_on_rate = 0.0;
            titin.candidate_actin = best_explicit_titin_candidate(
                model, config, index, thick_anchor);
            if (titin.candidate_actin < 0) continue;
            ++diagnostics.titin_candidate_links;
            vec candidate_displacement = thick_anchor -
                static_cast<vec>(model.actin.left_end[titin.candidate_actin]);
            candidate_displacement.pbc_wrap(model.box, model.is_periodic);
            titin.candidate_distance = candidate_displacement.norm();
            titin.candidate_energy = tension_only_energy(
                titin.candidate_distance,
                config.explicit_titin_slack_length,
                config.explicit_titin_k);
            double candidate_base_rate = config.explicit_titin_k_on;
            if (config.explicit_titin_rebind_memory > 0.0 &&
                titin.candidate_actin == titin.last_bound_actin &&
                model.current_step >= titin.last_detach_step &&
                model.dt * static_cast<double>(
                    model.current_step - titin.last_detach_step) <=
                    config.explicit_titin_rebind_memory) {
                candidate_base_rate *=
                    config.explicit_titin_rebind_rate_multiplier;
            }
            titin.candidate_on_rate = candidate_base_rate;
            if (config.explicit_titin_strain_gated_kinetics) {
                titin.candidate_on_rate = strain_gated_association_rate(
                    candidate_base_rate, model.thermal_beta,
                    config.explicit_titin_transition_theta,
                    titin.candidate_energy);
            }
            const double bind_probability =
                -std::expm1(-titin.candidate_on_rate * model.dt);
            counter_rng::CounterRng binding_rng(
                static_cast<std::uint64_t>(model.initial_seed),
                static_cast<std::uint64_t>(model.current_step),
                counter_rng::Stream::explicit_titin_binding,
                static_cast<std::uint64_t>(index),
                static_cast<std::uint64_t>(titin.candidate_actin));
            if (binding_rng.uniform() >= bind_probability) continue;
            titin.bound = true;
            titin.bound_actin = titin.candidate_actin;
            titin.last_bound_actin = titin.bound_actin;
            titin.candidate_actin = -1;
            titin.attach_step = model.current_step;
            titin.attachment_distance = titin.candidate_distance;
            titin.attachment_energy = titin.candidate_energy;
            titin.attachment_force = tension_only_force(
                titin.attachment_distance,
                config.explicit_titin_slack_length,
                config.explicit_titin_k);
            ++diagnostics.titin_newly_bound_links;
            diagnostics.titin_attachment_energy_total +=
                titin.attachment_energy;
            diagnostics.titin_attachment_energy_max = std::max(
                diagnostics.titin_attachment_energy_max,
                titin.attachment_energy);
            diagnostics.titin_attachment_force_max = std::max(
                diagnostics.titin_attachment_force_max,
                titin.attachment_force);
        }
        if (titin.bound_actin < 0 ||
            titin.bound_actin >= model.actin.n) {
            throw std::runtime_error("explicit titin has invalid actin site");
        }
        const int actin_index = titin.bound_actin;
        const vec actin_anchor = model.actin.left_end[actin_index];
        vec displacement = thick_anchor - actin_anchor;
        displacement.pbc_wrap(model.box, model.is_periodic);
        titin.distance = displacement.norm();
        titin.extension = tension_only_extension(
            titin.distance, config.explicit_titin_slack_length);
        titin.force = tension_only_force(
            titin.distance, config.explicit_titin_slack_length,
            config.explicit_titin_k);

        if (config.explicit_titin_k_off > 0.0) {
            double off_rate = config.explicit_titin_k_off;
            const double elastic_energy = tension_only_energy(
                titin.distance, config.explicit_titin_slack_length,
                config.explicit_titin_k);
            if (config.explicit_titin_strain_gated_kinetics) {
                off_rate = strain_gated_dissociation_rate(
                    config.explicit_titin_k_off, model.thermal_beta,
                    config.explicit_titin_transition_theta,
                    elastic_energy);
            } else if (config.explicit_titin_force_scale > 0.0) {
                off_rate *= std::exp(std::min(
                    50.0, titin.force / config.explicit_titin_force_scale));
            }
            const double release_probability =
                -std::expm1(-off_rate * model.dt);
            counter_rng::CounterRng release_rng(
                static_cast<std::uint64_t>(model.initial_seed),
                static_cast<std::uint64_t>(model.current_step),
                counter_rng::Stream::explicit_titin_release,
                static_cast<std::uint64_t>(index),
                static_cast<std::uint64_t>(actin_index));
            if (release_rng.uniform() < release_probability) {
                titin.last_bound_actin = actin_index;
                titin.bound = false;
                titin.bound_actin = -1;
                titin.candidate_actin = -1;
                titin.last_detach_step = model.current_step;
                titin.distance = 0.0;
                titin.extension = 0.0;
                titin.force = 0.0;
                ++diagnostics.titin_released_links;
                continue;
            }
        }

        ++diagnostics.titin_bound_links;
        diagnostics.titin_mean_distance += titin.distance;
        diagnostics.titin_total_force += titin.force;
        if (!(titin.extension > 0.0) || titin.distance <= 1e-12) {
            ++diagnostics.titin_compressed_links;
            continue;
        }
        ++diagnostics.titin_tension_links;
        const vec force_on_actin =
            (titin.force / titin.distance) * displacement;
        const vec force_on_myosin = -force_on_actin;
        model.actin.force[actin_index] += force_on_actin;
        model.myosin.force[titin.myosin] += force_on_myosin;

        const auto actin_torque = compute_am_offcenter_torque(
            model.actin.center[actin_index],
            model.actin.direction[actin_index], actin_anchor,
            force_on_actin, model.box, model.is_periodic);
        model.actin.torque[actin_index] += actin_torque.orientation_drive;
        const auto myosin_torque = compute_am_offcenter_torque(
            model.myosin.center[titin.myosin],
            model.myosin.direction[titin.myosin], thick_anchor,
            force_on_myosin, model.box, model.is_periodic);
        vec myosin_axis = model.myosin.direction[titin.myosin];
        myosin_axis.normalize();
        const vec axial_roll =
            myosin_torque.physical_torque.dot(myosin_axis) * myosin_axis;
        model.myosin.torque[titin.myosin] +=
            myosin_torque.orientation_drive + axial_roll;
        q += pair_generalized_force(
            force_on_actin.x, model.actin.center[actin_index].x,
            model.myosin.center[titin.myosin].x, model.box[0]);
    }
    if (diagnostics.titin_bound_links > 0) {
        diagnostics.titin_mean_distance /=
            static_cast<double>(diagnostics.titin_bound_links);
    }
    return q;
}

double add_register(Sarcomere& model, const Config& config,
                    double& net_added_force) {
    if (!(config.register_k > 0.0)) return 0.0;
    const double length = model.box[0];
    double q = 0.0;
    for (const auto& pair : aa_pairs) {
        const int first = pair[0];
        const int second = pair[1];
        const double delta = wrap_delta(
            model.actin.left_end[first].x -
            model.actin.left_end[second].x, length);
        const double overlap = std::abs(delta);
        const double sign = delta >= 0.0 ? 1.0 : -1.0;
        const double force_first = capped(
            -config.register_k *
                (overlap - config.register_target) * sign,
            config.register_force_cap);
        const vec force{force_first, 0.0, 0.0};
        model.actin.force[first] += force;
        model.actin.force[second] -= force;
        if (config.register_endpoint_torque) {
            // Filament::left_end = center - 0.5 * length * direction.
            // `torque` is the negative orientational energy gradient consumed
            // by the projected director integrator, not a cross-product
            // physical torque.
            const vec direction_drive{
                -0.5 * model.actin.length * force_first, 0.0, 0.0};
            model.actin.torque[first] += direction_drive;
            model.actin.torque[second] -= direction_drive;
        }
        net_added_force += force_first - force_first;
        q += pair_generalized_force(
            force_first, model.actin.center[first].x,
            model.actin.center[second].x, length);
    }
    return q;
}

double add_midpoint_feedback(Sarcomere& model, const Config& config,
                             double& net_added_force) {
    if (!(config.midpoint_feedback_k > 0.0)) return 0.0;
    const double length = model.box[0];
    double q = 0.0;
    for (const auto& pair : aa_pairs) {
        const int first = pair[0];
        const int second = pair[1];
        const double delta = wrap_delta(
            model.actin.left_end[first].x -
            model.actin.left_end[second].x, length);
        const double overlap = std::abs(delta);
        const double excess = overlap - config.midpoint_target;
        if (!(excess > 0.0)) continue;
        const double sign = delta >= 0.0 ? 1.0 : -1.0;
        const double restoring_first =
            -config.midpoint_feedback_k * excess * sign;
        const double force_first = capped(
            config.midpoint_feedback_sign * restoring_first,
            config.midpoint_force_cap);
        const double force_second = -force_first;
        const int myosin_first =
            owning_myosin[static_cast<size_t>(first)];
        const int myosin_second =
            owning_myosin[static_cast<size_t>(second)];
        model.actin.force[first] += {force_first, 0.0, 0.0};
        model.myosin.force[myosin_first] -= {force_first, 0.0, 0.0};
        model.actin.force[second] += {force_second, 0.0, 0.0};
        model.myosin.force[myosin_second] -= {force_second, 0.0, 0.0};
        net_added_force +=
            force_first - force_first + force_second - force_second;
        q += pair_generalized_force(
            force_first, model.actin.center[first].x,
            model.myosin.center[myosin_first].x, length);
        q += pair_generalized_force(
            force_second, model.actin.center[second].x,
            model.myosin.center[myosin_second].x, length);
    }
    return q;
}

}  // namespace

Mode parse_mode(const std::string& name) {
    if (name == "fixed") return Mode::fixed;
    if (name == "stress_clamp") return Mode::stress_clamp;
    if (name == "series_spring") return Mode::series_spring;
    throw std::invalid_argument("unknown axial cell mode: " + name);
}

const char* mode_name(Mode mode) {
    switch (mode) {
        case Mode::fixed: return "fixed";
        case Mode::stress_clamp: return "stress_clamp";
        case Mode::series_spring: return "series_spring";
    }
    return "unknown";
}

double wrap_delta(double value, double length) {
    if (!(length > 0.0)) return value;
    return value - length * std::round(value / length);
}

double pair_generalized_force(double force_on_first_x,
                              double first_x,
                              double second_x,
                              double length) {
    if (!(length > 0.0)) {
        throw std::invalid_argument("pair virial requires positive length");
    }
    return force_on_first_x *
           wrap_delta(first_x - second_x, length) / length;
}

double integrate_length(double length, double generalized_force,
                        double mobility, double dt,
                        double min_length, double max_length) {
    if (!(length > 0.0) || !(dt >= 0.0) || !(mobility >= 0.0) ||
        !(min_length > 0.0) || !(max_length >= min_length)) {
        throw std::invalid_argument("invalid axial cell integration input");
    }
    return std::clamp(
        length + mobility * generalized_force * dt,
        min_length, max_length);
}

double tension_only_extension(double distance, double rest_length) {
    if (!(distance >= 0.0) || !(rest_length >= 0.0)) {
        throw std::invalid_argument("invalid tension-only tether length");
    }
    return std::max(0.0, distance - rest_length);
}

double tension_only_force(double distance, double rest_length,
                          double stiffness) {
    if (!(stiffness >= 0.0)) {
        throw std::invalid_argument("invalid tension-only tether stiffness");
    }
    return stiffness * tension_only_extension(distance, rest_length);
}

double tension_only_energy(double distance, double rest_length,
                           double stiffness) {
    if (!(stiffness >= 0.0)) {
        throw std::invalid_argument("invalid tension-only tether stiffness");
    }
    const double extension = tension_only_extension(distance, rest_length);
    return 0.5 * stiffness * extension * extension;
}

double strain_gated_association_rate(double base_rate, double beta,
                                     double theta, double elastic_energy) {
    if (!(base_rate >= 0.0) || !(beta >= 0.0) || theta < 0.0 || theta > 1.0 ||
        !(elastic_energy >= 0.0)) {
        throw std::invalid_argument("invalid titin association kinetics");
    }
    const double exponent = std::max(-50.0, -beta * theta * elastic_energy);
    return base_rate * std::exp(exponent);
}

double strain_gated_dissociation_rate(double base_rate, double beta,
                                      double theta, double elastic_energy) {
    if (!(base_rate >= 0.0) || !(beta >= 0.0) || theta < 0.0 || theta > 1.0 ||
        !(elastic_energy >= 0.0)) {
        throw std::invalid_argument("invalid titin dissociation kinetics");
    }
    const double exponent = std::min(
        50.0, beta * (1.0 - theta) * elastic_energy);
    return base_rate * std::exp(exponent);
}

void initialize_explicit_titins(const Sarcomere& model, Config& config) {
    if (config.explicit_titin_per_half <= 0 || model.myosin.n <= 0) {
        throw std::invalid_argument(
            "explicit titin requires a positive tether count per myosin half");
    }
    config.explicit_titins.clear();
    config.explicit_titins.reserve(static_cast<size_t>(
        2 * model.myosin.n * config.explicit_titin_per_half));
    const int count = config.explicit_titin_per_half;
    for (int myosin = 0; myosin < model.myosin.n; ++myosin) {
        for (int half = 0; half < 2; ++half) {
            // Opposite halves use staggered material tracks. This is a local
            // thick-filament binding-site pattern, not a Z-ring template.
            const double phase = count > 1 && half == 0
                ? M_PI / static_cast<double>(count) : 0.0;
            for (int track = 0; track < count; ++track) {
                ExplicitTitin titin;
                titin.myosin = myosin;
                titin.half = half;
                titin.track = track;
                titin.anchor_azimuth = phase +
                    2.0 * M_PI * static_cast<double>(track) /
                        static_cast<double>(count);
                config.explicit_titins.push_back(titin);
            }
        }
    }
}

vec explicit_titin_thick_anchor(const Sarcomere& model,
                                const ExplicitTitin& titin) {
    if (titin.myosin < 0 || titin.myosin >= model.myosin.n ||
        (titin.half != 0 && titin.half != 1)) {
        throw std::invalid_argument("invalid explicit titin thick anchor");
    }
    return titin_thick_anchor(model, titin);
}

double aa_overlap(const Sarcomere& model, int first, int second) {
    if (first < 0 || second < 0 ||
        first >= model.actin.n || second >= model.actin.n ||
        model.box.empty()) {
        throw std::out_of_range("invalid AA overlap indices");
    }
    return std::abs(wrap_delta(
        model.actin.left_end[first].x -
        model.actin.left_end[second].x, model.box[0]));
}

double mean_aa_overlap(const Sarcomere& model) {
    if (!supported_topology(model)) {
        throw std::invalid_argument(
            "AA diagnostic requires balanced 4- or 12-actin/2-myosin topology");
    }
    if (model.actin.n == 4) {
        return 0.5 * (
            aa_overlap(model, 0, 1) + aa_overlap(model, 2, 3));
    }
    double total = 0.0;
    for (const auto& pair : hex_aa_pairs) {
        total += aa_overlap(model, pair[0], pair[1]);
    }
    return total / static_cast<double>(hex_aa_pairs.size());
}

double mean_linked_aa_overlap(const Sarcomere& model) {
    if (!supported_topology(model)) {
        throw std::invalid_argument(
            "linked AA diagnostic requires balanced 4- or 12-actin/2-myosin topology");
    }
    std::set<std::pair<int, int>> linked_pairs;
    for (const auto& link : model.aa_links) {
        if (!link.active) continue;
        const int first = std::min(link.actin_i, link.actin_j);
        const int second = std::max(link.actin_i, link.actin_j);
        if (model.actin.direction[first].dot(
                model.actin.direction[second]) >= 0.0) {
            continue;
        }
        linked_pairs.emplace(first, second);
    }
    if (linked_pairs.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    double total = 0.0;
    for (const auto& pair : linked_pairs) {
        total += aa_overlap(model, pair.first, pair.second);
    }
    return total / static_cast<double>(linked_pairs.size());
}

double mean_am_overlap(const Sarcomere& model) {
    if (!supported_topology(model)) {
        throw std::invalid_argument(
            "AM diagnostic requires balanced 4- or 12-actin/2-myosin topology");
    }
    double total = 0.0;
    for (int actin_index = 0; actin_index < model.actin.n; ++actin_index) {
        const int myosin_index = model.actin.n == 4
            ? owning_myosin[static_cast<size_t>(actin_index)]
            : hex_owning_myosin[static_cast<size_t>(actin_index)];
        vec center = static_cast<vec>(model.actin.center[actin_index]) -
                     static_cast<vec>(model.myosin.center[myosin_index]);
        center.pbc_wrap(model.box, model.is_periodic);
        vec axis = model.myosin.direction[myosin_index];
        axis.normalize();
        const double center_projection = center.dot(axis);
        const double half_projection = 0.5 * model.actin.length *
            model.actin.direction[actin_index].dot(axis);
        const double first = center_projection - half_projection;
        const double second = center_projection + half_projection;
        const double actin_min = std::min(first, second);
        const double actin_max = std::max(first, second);
        const double half_myosin = 0.5 * model.myosin.length;
        total += std::max(
            0.0, std::min(actin_max, half_myosin) -
                     std::max(actin_min, -half_myosin));
    }
    return total / static_cast<double>(model.actin.n);
}

Diagnostics prepare_step(Sarcomere& model, Config& config) {
    Diagnostics diagnostics;
    diagnostics.length_before =
        model.box.empty() ? 0.0 : model.box[0];
    diagnostics.length_after = diagnostics.length_before;
    if (config.mode == Mode::fixed &&
        !config.apply_internal_forces_when_fixed) return diagnostics;
    if (!supported_topology(model)) {
        throw std::invalid_argument(
            "variable axial cell requires balanced 4- or 12-actin/2-myosin topology");
    }
    if (model.actin.n != 4 &&
        (config.titin_k > 0.0 || config.register_k > 0.0 ||
         config.midpoint_feedback_k > 0.0)) {
        throw std::invalid_argument(
            "legacy pairwise axial restraints require the 4-actin topology");
    }

    diagnostics.q_am = am_cell_force(model);
    diagnostics.q_aa = aa_cell_force(model);
    diagnostics.q_titin =
        add_titin(model, config, diagnostics.net_added_force);
    diagnostics.q_titin +=
        add_explicit_titin(model, config, diagnostics);
    diagnostics.q_register =
        add_register(model, config, diagnostics.net_added_force);
    diagnostics.q_midpoint =
        add_midpoint_feedback(model, config, diagnostics.net_added_force);
    diagnostics.q_external = config.applied_tension;
    if (config.mode == Mode::series_spring) {
        diagnostics.q_external -= config.series_k *
            (diagnostics.length_before - config.series_rest_length);
    }
    diagnostics.q_total =
        diagnostics.q_am + diagnostics.q_aa + diagnostics.q_titin +
        diagnostics.q_register + diagnostics.q_midpoint +
        diagnostics.q_external;
    diagnostics.dlength_dt = config.mode == Mode::fixed
        ? 0.0 : config.mobility * diagnostics.q_total;
    return diagnostics;
}

void set_axial_length(Sarcomere& model, double new_length) {
    if (model.box.empty() || !(model.box[0] > 0.0) ||
        !(new_length > 0.0)) {
        throw std::invalid_argument("axial cell lengths must be positive");
    }
    const double scale = new_length / model.box[0];
    for (int i = 0; i < model.actin.n; ++i) {
        model.actin.center_x[static_cast<size_t>(i)] *= scale;
    }
    for (int i = 0; i < model.myosin.n; ++i) {
        model.myosin.center_x[static_cast<size_t>(i)] *= scale;
    }
    model.box[0] = new_length;
    model.actin.box[0] = new_length;
    model.myosin.box[0] = new_length;
    model.actin.set_periodic_axes(model.is_periodic);
    model.myosin.set_periodic_axes(model.is_periodic);
    model.actin.update_endpoints();
    model.myosin.update_endpoints();
    const double skin = 0.15 * model.cutoff_radius;
    model.neighbor_list = NeighborList(
        model.cutoff_radius + skin, model.box, skin / 2.0,
        model.is_periodic);
    model.neighbor_list.initialize(
        model.actin.center_x, model.actin.center_y, model.actin.center_z,
        model.myosin.center_x, model.myosin.center_y,
        model.myosin.center_z);
}

void advance_cell(Sarcomere& model, const Config& config,
                  double dt, Diagnostics& diagnostics) {
    if (config.mode == Mode::fixed) return;
    const double new_length = integrate_length(
        diagnostics.length_before, diagnostics.q_total,
        config.mobility, dt, config.min_length, config.max_length);
    if (new_length != diagnostics.length_before) {
        set_axial_length(model, new_length);
    }
    diagnostics.length_after = new_length;
}

}  // namespace axial_cell
