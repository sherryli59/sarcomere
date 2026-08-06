#include "sarcomere.h"

#include <array>
#include <cassert>
#include <cmath>
#include <set>
#include <string>
#include <tuple>

namespace {

int active_links(const Sarcomere& model, int first, int second) {
    if (first > second) std::swap(first, second);
    int count = 0;
    for (const auto& link : model.aa_links) {
        if (link.active && link.actin_i == first && link.actin_j == second) {
            ++count;
        }
    }
    return count;
}

int occupied_heads(const Sarcomere& model, int myosin, int half, int actin) {
    int count = 0;
    for (const auto& head : model.am_heads_by_myosin[myosin][half]) {
        if (head.bound_actin == actin &&
            head.state == static_cast<Sarcomere::BondFlag>(
                am_motor_kinetics::State::weak)) {
            ++count;
        }
    }
    return count;
}

double maximum_seeded_axial_strain(const Sarcomere& model, int myosin,
                                   int half, int actin) {
    double maximum = 0.0;
    for (const auto& head : model.am_heads_by_myosin[myosin][half]) {
        if (head.bound_actin != actin ||
            head.state != static_cast<Sarcomere::BondFlag>(
                am_motor_kinetics::State::weak)) {
            continue;
        }
        const auto geometry = compute_am_two_spring_head(
            model.actin, model.myosin, actin, myosin, model.box,
            head.target_axial, head.anchor_axial, head.binding_azimuth,
            model.am_spatial_heads ? model.myosin.radius : 0.0,
            model.am_spatial_heads ? &model.myosin_material_e1[myosin] : nullptr,
            model.am_spatial_heads ? &model.myosin_material_e2[myosin] : nullptr,
            model.actin_radius + model.am_surface_gap, 0.0, 0.0, 0.0,
            0.0, model.am_stall_force, false);
        maximum = std::max(maximum, std::abs(geometry.productive_strain));
    }
    return maximum;
}

}  // namespace

int main() {
    int n_actins = 4;
    int n_myosins = 4;
    std::vector<double> box{5.32, 3.0, 3.0};
    double actin_length = 1.0;
    double myosin_length = 1.5;
    double myosin_radius = 0.025;
    double am_cutoff = 0.05;
    double am_optimal = 0.055;
    double aa_cutoff = 0.0706;
    double aa_optimal = 0.0506;
    double k_on = 500.0;
    double k_off = 1.0;
    double base_lifetime = 0.05;
    double directional_lifetime = 0.05;
    double lifetime_coeff = 0.01;
    double diff_ratio = 1.0;
    double k_aa = 3000.0;
    double kappa_aa = 50.0;
    double zero = 0.0;
    std::string filename =
        "/tmp/minimal_double_myosin_balanced_initialization_test.h5";
    int seed = 3217;
    int fixed_myosins = 0;
    double dt = 1e-6;
    bool directional = true;

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, seed);
    Sarcomere model(
        n_actins, n_myosins, box, actin_length, myosin_length,
        myosin_radius, am_cutoff, am_optimal, aa_cutoff, aa_optimal,
        k_on, k_off, base_lifetime, directional_lifetime,
        lifetime_coeff, diff_ratio, k_aa, kappa_aa,
        zero, zero, zero, zero,
        filename, rng, seed, fixed_myosins, dt, 0.002, zero, zero,
        directional, 2, zero, zero, zero, zero);
    model.set_periodicity({true, false, false});
    model.set_local_aa_model(true, 2, 0.0192, 0.1, 50.0, 0);
    model.set_aa_axial_overlap(0.1);
    model.set_am_radial_geometry(0.015, 0.015, 1000.0);
    model.set_am_kinetics(
        true, 1000.0, 200.0, 200.0, 80.0, 1.0,
        1000.0, 0.005, 4.5, 5.0, 0.1, 0.1,
        0.05, 0.005, 0.03, 0.0, 0.0005);
    model.set_am_proximity_tether(1000.0, 0.005, true, 180.0);
    model.set_am_head_count(35);
    model.init_aa_overlap = 0.16;
    model.minimal_catch_variant("balanced_double_myosin");

    assert(active_links(model, 0, 1) == 2);
    assert(active_links(model, 2, 3) == 2);

    const std::array<std::tuple<int, int, int>, 8> expected_pairs{{
        {0, 0, 0}, {2, 0, 0},
        {0, 1, 2}, {2, 1, 2},
        {1, 1, 1}, {3, 1, 1},
        {1, 0, 3}, {3, 0, 3},
    }};
    const double seeded_strain_limit =
        3.0 * model.actin.length / model.am_sites_per_actin + 1e-12;
    for (const auto& [myosin, half, actin] : expected_pairs) {
        assert(occupied_heads(model, myosin, half, actin) == 6);
        assert(maximum_seeded_axial_strain(model, myosin, half, actin) <
               seeded_strain_limit);
    }

    // Binding-site occupancy is finite at the head level, but is not an
    // artificial one-myosin-per-actin ownership rule.  Every actin must have
    // heads from exactly two distinct thick filaments at frame zero.
    for (int actin = 0; actin < n_actins; ++actin) {
        std::set<int> attached_myosins;
        std::set<double> occupied_material_sites;
        int attached_heads = 0;
        for (int myosin = 0; myosin < n_myosins; ++myosin) {
            for (int half = 0; half < 2; ++half) {
                const int count = occupied_heads(model, myosin, half, actin);
                if (count > 0) attached_myosins.insert(myosin);
                attached_heads += count;
                for (const auto& head :
                     model.am_heads_by_myosin[myosin][half]) {
                    if (head.bound_actin == actin &&
                        head.state == static_cast<Sarcomere::BondFlag>(
                            am_motor_kinetics::State::weak)) {
                        assert(occupied_material_sites.insert(
                                   head.target_axial).second);
                    }
                }
            }
        }
        assert(attached_myosins.size() == 2);
        assert(attached_heads == 12);
        assert(occupied_material_sites.size() == 12);
    }

    // The two thick filaments sharing an actin are axially registered and sit
    // on opposite transverse sides at the preferred centerline separation.
    const double preferred_reach =
        model.myosin.radius + model.actin_radius + model.am_surface_gap;
    const std::array<std::tuple<int, int, int>, 4> transverse_pairs{{
        {0, 2, 0}, {0, 2, 2}, {1, 3, 1}, {1, 3, 3},
    }};
    for (const auto& [upper, lower, actin] : transverse_pairs) {
        assert(std::abs(model.myosin.center[upper].x -
                        model.myosin.center[lower].x) < 1e-12);
        for (const int myosin : {upper, lower}) {
            vec displacement = model.actin.center[actin] -
                               model.myosin.center[myosin];
            displacement.pbc_wrap(box, model.is_periodic);
            const double transverse =
                std::hypot(displacement.y, displacement.z);
            assert(std::abs(transverse - preferred_reach) < 1e-12);
        }
    }

    // The annealing initializer changes only relative thick-filament
    // register: the two pair midpoints and every transverse separation remain
    // unchanged, while each pair begins 150 nm out of axial register.
    const double right_midpoint_before =
        0.5 * (model.myosin.center[0].x + model.myosin.center[2].x);
    const double left_midpoint_before =
        0.5 * (model.myosin.center[1].x + model.myosin.center[3].x);
    model.minimal_catch_variant("balanced_double_myosin_anneal");
    assert(std::abs(model.myosin.center[2].x - model.myosin.center[0].x -
                    0.150) < 1e-12);
    assert(std::abs(model.myosin.center[3].x - model.myosin.center[1].x +
                    0.150) < 1e-12);
    assert(std::abs(0.5 * (model.myosin.center[0].x +
                           model.myosin.center[2].x) -
                    right_midpoint_before) < 1e-12);
    assert(std::abs(0.5 * (model.myosin.center[1].x +
                           model.myosin.center[3].x) -
                    left_midpoint_before) < 1e-12);
    for (const auto& [myosin, half, actin] : expected_pairs) {
        assert(occupied_heads(model, myosin, half, actin) == 6);
        assert(maximum_seeded_axial_strain(model, myosin, half, actin) <
               seeded_strain_limit);
    }

    gsl_rng_free(rng);
    return 0;
}
