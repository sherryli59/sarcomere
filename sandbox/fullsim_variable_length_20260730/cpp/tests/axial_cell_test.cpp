#include "axial_cell.h"

#include <cassert>
#include <cmath>
#include <memory>
#include <string>

namespace {

std::unique_ptr<Sarcomere> make_model(gsl_rng* rng) {
    int n_actins = 4;
    int n_myosins = 2;
    std::vector<double> box{4.3, 3.0, 3.0};
    double actin_length = 1.0;
    double myosin_length = 1.5;
    double myosin_radius = 0.025;
    double am_cutoff = 0.05;
    double am_optimal = 0.055;
    double aa_cutoff = 0.0706;
    double aa_optimal = 0.0506;
    double k_on = 0.0;
    double k_off = 0.0;
    double base_lifetime = 1.0e9;
    double directional_lifetime = 1.0e9;
    double lifetime_coeff = 0.0;
    double diff_ratio = 1.0;
    double k_aa = 3000.0;
    double kappa_aa = 50.0;
    double zero = 0.0;
    std::string filename = "/tmp/axial_cell_test.h5";
    int seed = 4901;
    int fixed_myosins = 0;
    double dt = 1e-5;
    bool directional = true;
    auto model = std::make_unique<Sarcomere>(
        n_actins, n_myosins, box, actin_length, myosin_length,
        myosin_radius, am_cutoff, am_optimal, aa_cutoff, aa_optimal,
        k_on, k_off, base_lifetime, directional_lifetime,
        lifetime_coeff, diff_ratio, k_aa, kappa_aa,
        zero, zero, zero, zero, filename, rng, seed, fixed_myosins,
        dt, zero, zero, zero, directional, 2,
        zero, zero, zero, zero);
    model->set_periodicity({true, false, false});
    model->set_local_aa_model(true, 2, 0.0192, 0.1, 50.0, 0);
    model->set_aa_axial_overlap(0.1);
    model->aa_barbed_zone = 0.1;
    model->set_am_radial_geometry(0.015, 0.015, 1000.0);
    model->set_am_kinetics(
        true, 1000.0, 200.0, 200.0, 80.0, 1.0,
        1000.0, 0.005, 4.5, 5.0, 0.1, 0.1,
        0.05, 0.005, 0.03, 0.0, 0.0005);
    model->set_am_proximity_tether(1000.0, 0.005, true, 180.0);
    model->set_am_head_count(35);
    model->minimal_catch_variant("balanced");
    return model;
}

vec explicit_titin_anchor(const Sarcomere& model,
                          const axial_cell::ExplicitTitin& titin) {
    vec axis = model.myosin.direction[titin.myosin];
    axis.normalize();
    vec e1 = model.myosin_material_e1[titin.myosin];
    e1 -= axis * e1.dot(axis);
    e1.normalize();
    vec e2 = axis.cross(e1);
    e2.normalize();
    const double half_sign = titin.half == 0 ? -1.0 : 1.0;
    return static_cast<vec>(model.myosin.center[titin.myosin]) +
        half_sign * 0.5 * model.myosin.length * axis +
        model.myosin.radius *
            (std::cos(titin.anchor_azimuth) * e1 +
             std::sin(titin.anchor_azimuth) * e2);
}

}  // namespace

int main() {
    using namespace axial_cell;
    constexpr double tol = 1e-12;

    // Pair virial is invariant when either endpoint is represented by a
    // periodic image.
    const double q0 = pair_generalized_force(3.0, 2.05, -2.05, 4.3);
    const double q1 = pair_generalized_force(3.0, 2.05, 2.25, 4.3);
    assert(std::abs(q0 - q1) < tol);

    // Passive strain gating preserves detailed balance for the elastic
    // contribution: (kon/koff)/(kon0/koff0) = exp(-beta U).
    const double elastic_energy = tension_only_energy(0.503, 0.375, 220.0);
    assert(std::abs(
        elastic_energy - 0.5 * 220.0 * 0.128 * 0.128) < tol);
    const double beta = 1.0 / 0.00428;
    const double gated_on = strain_gated_association_rate(
        5000.0, beta, 0.5, elastic_energy);
    const double gated_off = strain_gated_dissociation_rate(
        0.5, beta, 0.5, elastic_energy);
    assert(gated_on < 1e-15);
    assert(gated_off > 1e15);
    assert(strain_gated_association_rate(7.0, beta, 0.5, 0.0) == 7.0);
    assert(strain_gated_dissociation_rate(3.0, beta, 0.5, 0.0) == 3.0);
    assert(q0 < 0.0);

    // Cell integration signs and clamping.
    assert(integrate_length(4.3, 0.0, 0.1, 0.01, 3.0, 6.0) == 4.3);
    assert(integrate_length(4.3, 2.0, 0.1, 0.01, 3.0, 6.0) > 4.3);
    assert(integrate_length(4.3, -2.0, 0.1, 0.01, 3.0, 6.0) < 4.3);

    // A physical titin-like tether is slack in compression.
    assert(tension_only_extension(1.0, 1.125) == 0.0);
    assert(tension_only_force(1.0, 1.125, 220.0) == 0.0);
    assert(tension_only_force(1.125, 1.125, 220.0) == 0.0);
    assert(std::abs(
        tension_only_force(1.25, 1.125, 220.0) - 27.5) < tol);

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, 4901);
    auto model_ptr = make_model(rng);
    auto& model = *model_ptr;
    const double initial_length = model.box[0];
    const double initial_fraction = model.actin.center[0].x / initial_length;
    assert(std::abs(mean_aa_overlap(model) - 0.1) < 1e-10);
    assert(std::abs(mean_am_overlap(model) - 0.625) < 1e-10);

    // Fixed mode is a strict no-op: it adds no force and does not rescale.
    Config fixed;
    const vec force_before = model.actin.force[0];
    auto fixed_diag = prepare_step(model, fixed);
    advance_cell(model, fixed, 1e-3, fixed_diag);
    assert(model.box[0] == initial_length);
    assert((static_cast<vec>(model.actin.force[0]) - force_before).norm() == 0.0);

    // Affine change preserves fractional centers and rod directions.
    const vec direction_before = model.actin.direction[0];
    set_axial_length(model, 4.7);
    assert(std::abs(model.actin.center[0].x / model.box[0] -
                    initial_fraction) < tol);
    assert((static_cast<vec>(model.actin.direction[0]) -
            direction_before).norm() < tol);

    // A compressed always-connected titin network is expansile.
    set_axial_length(model, 3.9);
    Config titin;
    titin.mode = Mode::stress_clamp;
    titin.mobility = 0.01;
    titin.titin_k = 100.0;
    titin.titin_rest_length = 1.25;
    auto titin_diag = prepare_step(model, titin);
    assert(titin_diag.q_titin > 0.0);
    assert(std::abs(titin_diag.net_added_force) < tol);

    // Register force has the correct restoring sign on either side.
    Config register_config;
    register_config.mode = Mode::stress_clamp;
    register_config.register_k = 100.0;
    register_config.register_target = 0.1;
    model.actin.center_x[0] += 0.02;
    model.actin.update_endpoints();
    auto register_diag = prepare_step(model, register_config);
    assert(register_diag.q_register != 0.0);
    assert(std::abs(register_diag.net_added_force) < tol);

    // When the register coordinate is an endpoint separation, the optional
    // orientational term is the exact matching negative energy gradient.
    model.actin.direction[0] =
        {std::cos(0.2), std::sin(0.2), 0.0};
    model.actin.update_endpoints();
    const double endpoint_delta = wrap_delta(
        model.actin.left_end[0].x - model.actin.left_end[1].x,
        model.box[0]);
    const double endpoint_sign = endpoint_delta >= 0.0 ? 1.0 : -1.0;
    const double expected_force =
        -register_config.register_k *
        (std::abs(endpoint_delta) - register_config.register_target) *
        endpoint_sign;
    const vec torque_first_before = model.actin.torque[0];
    const vec torque_second_before = model.actin.torque[1];
    register_config.register_endpoint_torque = true;
    prepare_step(model, register_config);
    const vec torque_first_added =
        static_cast<vec>(model.actin.torque[0]) - torque_first_before;
    const vec torque_second_added =
        static_cast<vec>(model.actin.torque[1]) - torque_second_before;
    assert(std::abs(
        torque_first_added.x +
        0.5 * model.actin.length * expected_force) < tol);
    assert((torque_first_added + torque_second_added).norm() < tol);

    // Explicit titin discovers four individual barbed-end sites without any
    // myosin-head occupancy or ACTN2-connected component.
    auto dynamic_model_ptr = make_model(rng);
    auto& dynamic_model = *dynamic_model_ptr;
    for (auto& link : dynamic_model.aa_links) link.active = false;
    for (auto& myosin_heads : dynamic_model.am_heads_by_myosin) {
        for (auto& half_heads : myosin_heads) {
            for (auto& head : half_heads) {
                head.bound_actin = -1;
                head.state = static_cast<Sarcomere::BondFlag>(
                    am_motor_kinetics::State::detached);
            }
        }
    }
    Config dynamic;
    dynamic.mode = Mode::fixed;
    dynamic.apply_internal_forces_when_fixed = true;
    dynamic.explicit_titin = true;
    dynamic.explicit_titin_per_half = 1;
    dynamic.explicit_titin_k = 220.0;
    dynamic.explicit_titin_slack_length = 1.125;
    dynamic.explicit_titin_capture_distance = 1.5;
    dynamic.explicit_titin_k_on = 1.0e12;
    initialize_explicit_titins(dynamic_model, dynamic);
    auto capture_diag = prepare_step(dynamic_model, dynamic);
    assert(capture_diag.titin_bound_links == 4);
    assert(capture_diag.titin_candidate_links == 4);
    std::array<bool, 4> captured_actins{};
    for (const auto& titin_link : dynamic.explicit_titins) {
        assert(titin_link.bound);
        assert(titin_link.bound_actin >= 0 && titin_link.bound_actin < 4);
        captured_actins[static_cast<size_t>(titin_link.bound_actin)] = true;
    }
    for (bool captured : captured_actins) assert(captured);

    // Binding identity persists without relying on heads or ACTN2 memory.
    ++dynamic_model.current_step;
    auto memory_diag = prepare_step(dynamic_model, dynamic);
    assert(memory_diag.titin_bound_links == 4);

    // A coarse minimal actin may represent several resolved lattice tracks.
    // Multiple explicit molecules may bind it independently when the local
    // occupancy limit is raised; the default one-per-actin rule is unchanged.
    Config multivalent = dynamic;
    multivalent.explicit_titin_per_half = 3;
    multivalent.explicit_titin_max_per_actin = 3;
    multivalent.explicit_titins.clear();
    initialize_explicit_titins(dynamic_model, multivalent);
    auto multivalent_diag = prepare_step(dynamic_model, multivalent);
    assert(multivalent_diag.titin_bound_links == 12);
    std::array<int, 4> titins_per_actin{};
    for (const auto& titin_link : multivalent.explicit_titins) {
        assert(titin_link.bound);
        ++titins_per_actin[static_cast<size_t>(titin_link.bound_actin)];
    }
    for (int occupancy : titins_per_actin) assert(occupancy == 3);

    // One explicit titin works with one actin. It is slack below its rest
    // length and, above rest, exerts a reciprocal force only on that actin and
    // its own myosin rather than distributing force across an AA component.
    Config single_link;
    single_link.mode = Mode::fixed;
    single_link.apply_internal_forces_when_fixed = true;
    single_link.explicit_titin = true;
    single_link.explicit_titin_per_half = 1;
    single_link.explicit_titin_k = 220.0;
    single_link.explicit_titin_slack_length = 1.125;
    single_link.explicit_titin_capture_distance = 1.5;
    ExplicitTitin one_titin;
    one_titin.myosin = 0;
    one_titin.half = 0;
    one_titin.track = 0;
    one_titin.bound = true;
    one_titin.bound_actin = 0;
    single_link.explicit_titins.push_back(one_titin);
    const auto zero_forces = [&]() {
        for (int i = 0; i < dynamic_model.actin.n; ++i) {
            dynamic_model.actin.force[i] = {0.0, 0.0, 0.0};
            dynamic_model.actin.torque[i] = {0.0, 0.0, 0.0};
        }
        for (int i = 0; i < dynamic_model.myosin.n; ++i) {
            dynamic_model.myosin.force[i] = {0.0, 0.0, 0.0};
            dynamic_model.myosin.torque[i] = {0.0, 0.0, 0.0};
        }
    };
    dynamic_model.actin.direction[0] = {1.0, 0.0, 0.0};
    const vec negative_anchor = explicit_titin_anchor(
        dynamic_model, single_link.explicit_titins[0]);
    vec barbed_site = negative_anchor - vec{1.0, 0.0, 0.0};
    dynamic_model.actin.center[0] = barbed_site + vec{0.5, 0.0, 0.0};
    int actin_zero = 0;
    dynamic_model.actin.update_endpoints(actin_zero);
    zero_forces();
    auto compressed_diag = prepare_step(dynamic_model, single_link);
    assert(compressed_diag.titin_compressed_links == 1);
    assert(compressed_diag.titin_total_force == 0.0);
    assert(static_cast<vec>(dynamic_model.actin.force[0]).norm() == 0.0);
    assert(static_cast<vec>(dynamic_model.actin.torque[0]).norm() == 0.0);

    barbed_site = negative_anchor - vec{1.25, 0.0, 0.0};
    dynamic_model.actin.center[0] = barbed_site + vec{0.5, 0.0, 0.0};
    dynamic_model.actin.update_endpoints(actin_zero);
    zero_forces();
    auto tension_diag = prepare_step(dynamic_model, single_link);
    assert(tension_diag.titin_tension_links == 1);
    assert(std::abs(tension_diag.titin_total_force - 27.5) < 1e-9);
    assert(std::abs(dynamic_model.actin.force[0].x - 27.5) < 1e-9);
    assert(static_cast<vec>(dynamic_model.actin.force[1]).norm() == 0.0);
    assert((static_cast<vec>(dynamic_model.actin.force[0]) +
            static_cast<vec>(dynamic_model.myosin.force[0])).norm() < tol);
    assert(static_cast<vec>(dynamic_model.actin.torque[0]).norm() == 0.0);

    // Turnover is titin-local and independent of motor/ACTN2 state.
    single_link.explicit_titin_k_off = 1.0e12;
    ++dynamic_model.current_step;
    auto release_diag = prepare_step(dynamic_model, single_link);
    assert(release_diag.titin_released_links == 1);
    assert(!single_link.explicit_titins[0].bound);
    assert(single_link.explicit_titins[0].last_detach_step ==
           dynamic_model.current_step);
    assert(single_link.explicit_titins[0].last_bound_actin == 0);

    // A just-detached molecule may rebind the same nearby molecular partner
    // even if a global half/polarity discovery gate changed sign at the
    // detachment boundary. The memory applies no force while unbound and
    // expires after the configured finite interval.
    single_link.explicit_titin_k_off = 0.0;
    single_link.explicit_titin_k_on = 1.0e12;
    single_link.explicit_titin_rebind_memory = 0.020;
    barbed_site = static_cast<vec>(dynamic_model.myosin.center[0]) +
        vec{0.10, 0.0, 0.0};
    dynamic_model.actin.center[0] = barbed_site + vec{0.5, 0.0, 0.0};
    dynamic_model.actin.update_endpoints(actin_zero);
    ++dynamic_model.current_step;
    auto rebind_diag = prepare_step(dynamic_model, single_link);
    assert(rebind_diag.titin_newly_bound_links == 1);
    assert(single_link.explicit_titins[0].bound);
    assert(single_link.explicit_titins[0].bound_actin == 0);

    gsl_rng_free(rng);
    return 0;
}
