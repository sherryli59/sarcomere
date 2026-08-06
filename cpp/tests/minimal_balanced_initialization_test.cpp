#include "sarcomere.h"

#include <cassert>
#include <cmath>
#include <string>

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
    int n_myosins = 2;
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
    std::string filename = "/tmp/minimal_balanced_initialization_test.h5";
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
    model.minimal_catch_variant("balanced");

    assert(model.aa_link_register_ratio(1) == 0.0);
    assert(model.aa_link_register_ratio(2) == 1.0);
    model.aa_register_all_statuses = true;
    assert(model.aa_link_register_ratio(1) == 1.0);

    // Both the central and periodic-image Z-discs must exist at frame zero.
    assert(active_links(model, 0, 1) == 2);
    assert(active_links(model, 2, 3) == 2);

    // Each myosin half has the correct-polarity actin and six seeded heads.
    assert(occupied_heads(model, 0, 0, 0) == 6);
    assert(occupied_heads(model, 0, 1, 2) == 6);
    assert(occupied_heads(model, 1, 1, 1) == 6);
    assert(occupied_heads(model, 1, 0, 3) == 6);
    // Seeding may not manufacture endpoint grabs from heads beyond the finite
    // actin. Explicit monomer-scale sites replace the old continuously
    // sliding target, so a three-head crown can begin within a few site
    // spacings of zero axial strain rather than at machine-zero strain.
    const double seeded_strain_limit =
        3.0 * model.actin.length / model.am_sites_per_actin + 1e-12;
    assert(maximum_seeded_axial_strain(model, 0, 0, 0) < seeded_strain_limit);
    assert(maximum_seeded_axial_strain(model, 0, 1, 2) < seeded_strain_limit);
    assert(maximum_seeded_axial_strain(model, 1, 1, 1) < seeded_strain_limit);
    assert(maximum_seeded_axial_strain(model, 1, 0, 3) < seeded_strain_limit);

    // Minimum-image center separation of each antiparallel pair is L-overlap.
    vec central = model.actin.center[1] - model.actin.center[0];
    vec periodic = model.actin.center[3] - model.actin.center[2];
    central.pbc_wrap(box, model.is_periodic);
    periodic.pbc_wrap(box, model.is_periodic);
    const double axial_separation = actin_length - model.init_aa_overlap;
    assert(std::abs(central.norm() -
                    std::sqrt(axial_separation * axial_separation +
                              aa_optimal * aa_optimal)) < 1e-8);
    assert(std::abs(periodic.norm() -
                    std::sqrt(axial_separation * axial_separation +
                              aa_optimal * aa_optimal)) < 1e-8);

    // Splitting the A-A overlap symmetrically keeps both Z-pair centroids at
    // their nominal lattice positions.  Consequently, all four actins begin
    // with the same axial center separation from their owning myosin.
    const std::array<int, 4> initial_owner{{0, 1, 0, 1}};
    const double expected_am_center_separation =
        0.25 * box[0] - 0.5 * actin_length +
        0.5 * model.init_aa_overlap;
    for (int i = 0; i < 4; ++i) {
        vec displacement = model.actin.center[i] -
                           model.myosin.center[initial_owner[i]];
        displacement.pbc_wrap(box, model.is_periodic);
        assert(std::abs(std::abs(displacement.x) -
                        expected_am_center_separation) < 1e-12);
    }

    // The promoted titin topology is reciprocal and independent of current
    // AM/AA occupancy: each barbed end maps to its owning myosin.
    model.dynamic_am = false;
    model.aa_local_multilink = false;
    model.k_aa = 0.0;
    model.kappa_aa = 0.0;
    model.aa_perpendicular_k = 0.0;
    model.sterics_enabled = false;
    model.titin_k = 220.0;
    model.titin_rest_length = 1.125;
    model.titin_mode = Sarcomere::TitinMode::balanced_owning;
    const std::array<int, 4> owner{{0, 1, 0, 1}};
    std::array<vec, 4> expected_actin;
    std::array<vec, 2> expected_myosin{{vec{0, 0, 0}, vec{0, 0, 0}}};
    for (int i = 0; i < 4; ++i) {
        vec displacement = model.myosin.center[owner[i]] -
                           model.actin.left_end[i];
        displacement.pbc_wrap(box, model.is_periodic);
        const double distance = displacement.norm();
        expected_actin[i] =
            (model.titin_k * (distance - model.titin_rest_length) / distance) *
            displacement;
        expected_myosin[owner[i]] -= expected_actin[i];
    }
    model.update_system();
    for (int i = 0; i < 4; ++i) {
        assert((static_cast<vec>(model.actin.force[i]) -
                expected_actin[i]).norm() < 1e-10);
    }
    for (int i = 0; i < 2; ++i) {
        assert((static_cast<vec>(model.myosin.force[i]) -
                expected_myosin[i]).norm() < 1e-10);
    }

    gsl_rng_free(rng);
    return 0;
}
