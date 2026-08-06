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
    const double locked_k = model.k_aa /
        static_cast<double>(model.aa_stiffness_reference_links);
    assert(std::abs(model.aa_link_radial_stiffness(1) - locked_k) < 1e-12);
    assert(std::abs(model.aa_link_radial_stiffness(2) - locked_k) < 1e-12);
    model.aa_status1_radial_k = 75.0;
    assert(std::abs(model.aa_link_radial_stiffness(1) - 75.0) < 1e-12);
    assert(std::abs(model.aa_link_radial_stiffness(2) - locked_k) < 1e-12);
    model.set_aa_axial_overlap(0.1);
    model.set_am_radial_geometry(0.015, 0.015, 1000.0);
    model.set_am_kinetics(
        true, 1000.0, 200.0, 200.0, 80.0, 1.0,
        1000.0, 0.005, 4.5, 5.0, 0.1, 0.1,
        0.05, 0.005, 0.03, 0.0, 0.0005);
    model.set_am_proximity_tether(1000.0, 0.005, true, 180.0);
    model.set_am_head_count(35);
    model.minimal_catch_variant("balanced");

    // Both the central and periodic-image Z-discs must exist at frame zero.
    assert(active_links(model, 0, 1) == 2);
    assert(active_links(model, 2, 3) == 2);

    // Each myosin half has the correct-polarity actin and six seeded heads.
    assert(occupied_heads(model, 0, 0, 0) == 6);
    assert(occupied_heads(model, 0, 1, 2) == 6);
    assert(occupied_heads(model, 1, 1, 1) == 6);
    assert(occupied_heads(model, 1, 0, 3) == 6);

    // Minimum-image center separation of each antiparallel pair is L-overlap.
    vec central = model.actin.center[1] - model.actin.center[0];
    vec periodic = model.actin.center[3] - model.actin.center[2];
    central.pbc_wrap(box, model.is_periodic);
    periodic.pbc_wrap(box, model.is_periodic);
    assert(std::abs(central.norm() -
                    std::sqrt(0.9 * 0.9 + aa_optimal * aa_optimal)) < 1e-8);
    assert(std::abs(periodic.norm() -
                    std::sqrt(0.9 * 0.9 + aa_optimal * aa_optimal)) < 1e-8);

    // The imposed A-A interdigitation must be centered on both periodic
    // Z-register planes. Otherwise alternate thick filaments begin with
    // systematically different tip-to-barbed-end (titin) gaps.
    assert(std::abs(
        model.actin.center[0].x + model.actin.center[1].x) < 1e-12);
    assert(std::abs(
        model.actin.center[2].x + model.actin.center[3].x) < 1e-12);

    gsl_rng_free(rng);
    return 0;
}
