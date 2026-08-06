#include "sarcomere.h"
#include "interaction.h"

#include <cassert>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <string>

namespace {
int active_link_count(const Sarcomere& model) {
    int count = 0;
    for (const auto& link : model.aa_links) {
        if (link.active) ++count;
    }
    return count;
}

int occupied_heads(const Sarcomere& model, int myosin, int half) {
    int count = 0;
    for (const auto& head : model.am_heads_by_myosin[myosin][half]) {
        if (head.bound_actin >= 0) ++count;
    }
    return count;
}
}

int main() {
    constexpr double d_am = 0.055;
    int n_actins = 16;
    int n_myosins = 4;
    std::vector<double> box{
        2.2, 2.0 * std::sqrt(3.0) * d_am, 3.0 * d_am};
    double actin_length = 1.0;
    double myosin_length = 1.5;
    double myosin_radius = 0.025;
    double am_cutoff = 0.05;
    double am_optimal = d_am;
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
    std::string filename = "/tmp/shared_lattice_initialization_test.h5";
    int seed = 9012;
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
        zero, zero, zero, zero, filename, rng, seed, fixed_myosins, dt,
        0.002, zero, zero, directional, 24, zero, zero, zero, zero);
    model.set_periodicity({true, true, true});
    model.set_local_aa_model(true, 2, 0.0192, 0.1, 50.0, 0);
    model.set_aa_neighbor_identity(true, 0.02, 0.20);
    model.max_aa_bonds_per_actin = 3;
    model.set_aa_axial_overlap(0.1);
    model.set_am_radial_geometry(0.015, 0.015, 1000.0);
    model.set_am_kinetics(
        true, 1000.0, 200.0, 200.0, 80.0, 1.0,
        1000.0, 0.005, 4.5, 5.0, 0.1, 0.1,
        0.05, 0.005, 0.03, 0.0, 0.0005);
    model.set_am_proximity_tether(1000.0, 0.005, true, 30.0);
    model.set_am_head_count(141);
    model.shared_thick_thin_lattice();
    model.prebind_shared_thick_thin_lattice(1);
    model.seed_shared_thick_thin_heads(6);

    assert(active_link_count(model) == 24);
    for (int actin = 0; actin < 16; ++actin) {
        assert(model.aa_neighbor_degree(actin) == 3);
    }
    for (int myosin = 0; myosin < 4; ++myosin) {
        assert(occupied_heads(model, myosin, 0) == 36);
        assert(occupied_heads(model, myosin, 1) == 36);
    }

    // No initialized actin-actin capsules overlap.
    for (int first = 0; first < 16; ++first) {
        for (int second = first + 1; second < 16; ++second) {
            const auto repulsion = compute_actin_repulsion(
                model.actin, first, second, box, 0.030, 100.0);
            assert(!repulsion.applied);
        }
    }

    gsl_rng_free(rng);
    return 0;
}
