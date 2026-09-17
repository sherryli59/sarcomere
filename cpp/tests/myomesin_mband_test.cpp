// Sparse, local, finite-valency myomesin/M-band bridge tests.

#include "sarcomere.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>

namespace {

std::unique_ptr<Sarcomere> make_model(
    const std::string& output, gsl_rng* rng, int& seed, double k_mm,
    int n_myosins = 2) {
    int n_actins = 1;
    std::vector<double> box{8.0, 8.0, 8.0};
    double actin_length = 1.0;
    double myosin_length = 1.5;
    double myosin_radius = 0.025;
    double am_cutoff = 0.10;
    double am_optimal = 0.055;
    double aa_cutoff = 0.08;
    double aa_optimal = 0.0506;
    double zero = 0.0;
    double diff_ratio = 1.0;
    int fixed_myosins = 0;
    double dt = 1e-6;
    bool directional = true;
    int max_myosin_bonds = 6;
    std::string filename = output;
    auto model = std::make_unique<Sarcomere>(
        n_actins, n_myosins, box, actin_length, myosin_length,
        myosin_radius, am_cutoff, am_optimal, aa_cutoff, aa_optimal,
        zero, zero, zero, zero, zero, diff_ratio, zero, zero, zero,
        zero, k_mm, zero, filename, rng, seed, fixed_myosins, dt, zero,
        zero, zero, directional, max_myosin_bonds, zero, zero, zero,
        zero);
    model->set_periodicity({false, false, false});
    model->dynamic_am = false;
    model->aa_local_multilink = false;
    model->sterics_enabled = true;
    model->actin.center[0] = {0.0, 3.5, 3.5};
    model->actin.direction[0] = {1.0, 0.0, 0.0};
    model->actin.update_endpoints();
    for (int i = 0; i < n_myosins; ++i) {
        model->myosin.direction[i] = {1.0, 0.0, 0.0};
    }
    model->reset_myosin_material_frames();
    return model;
}

void configure(Sarcomere& model, int valency, double k_on, double k_off,
               double radial_k = 100.0, double force_max = 100.0) {
    model.configure_myomesin_bridges(
        valency, k_on, k_off, 0.25, 0.11, 0.04, radial_k,
        0.94, 0.0, force_max);
}

void place_pair(Sarcomere& model, double axial, double transverse) {
    model.myosin.center[0] = {0.0, 0.0, 0.0};
    model.myosin.center[1] = {axial, transverse, 0.0};
    model.myosin.update_endpoints();
}

}  // namespace

int main() {
    const std::string output = "/tmp/myomesin_mband_test.h5";
    int seed = 5501;
    constexpr double K = 100.0;

    // An occupied bridge restores axial register, has zero transverse force
    // at its nonzero rest spacing, and is exactly reciprocal.
    {
        gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng, seed);
        auto model = make_model(output, rng, seed, K);
        configure(*model, 1, 0.0, 0.0);
        place_pair(*model, 0.15, 0.11);
        assert(model->prebind_myomesin_bridge(0, 1));
        model->update_system();
        const vec f0 = model->myosin.force[0];
        const vec f1 = model->myosin.force[1];
        assert(f0.x > 0.0 && f1.x < 0.0);
        assert(std::abs(f0.y) < 1e-9 && std::abs(f1.y) < 1e-9);
        assert((f0 + f1).norm() < 1e-9);
        assert(model->myomesin_bridges.size() == 1);
        model->sterics_enabled = false;
        model->update_system();
        assert(model->myosin.force[0].x > 0.0);
        assert(model->myosin.force[1].x < 0.0);
        gsl_rng_free(rng);
    }

    // A coarse-grained module's force capacity scales independently from its
    // axial stiffness.  Raising k_mm without raising the module cap must still
    // clamp the resultant force; raising both permits a larger restoring load.
    {
        gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
        auto capped = make_model(output, rng, seed, 100.0);
        configure(*capped, 1, 0.0, 0.0, 300.0, 10.0);
        place_pair(*capped, 0.15, 0.11);
        assert(capped->prebind_myomesin_bridge(0, 1));
        capped->update_system();
        assert(std::abs(capped->myomesin_bridges[0].force - 10.0) < 1e-9);

        auto scaled = make_model(output, rng, seed, 200.0);
        configure(*scaled, 1, 0.0, 0.0, 300.0, 20.0);
        place_pair(*scaled, 0.15, 0.11);
        assert(scaled->prebind_myomesin_bridge(0, 1));
        scaled->update_system();
        assert(std::abs(scaled->myomesin_bridges[0].force - 20.0) < 1e-9);
        gsl_rng_free(rng);
    }

    // The transverse spring preserves a finite lattice spacing instead of
    // collapsing thick-filament centerlines.
    {
        gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
        auto model = make_model(output, rng, seed, K);
        configure(*model, 1, 0.0, 0.0);
        place_pair(*model, 0.0, 0.14);
        assert(model->prebind_myomesin_bridge(0, 1));
        model->update_system();
        assert(model->myosin.force[0].y > 0.0);
        assert(model->myosin.force[1].y < 0.0);
        place_pair(*model, 0.0, 0.08);
        model->update_system();
        assert(model->myosin.force[0].y < 0.0);
        assert(model->myosin.force[1].y > 0.0);
        gsl_rng_free(rng);
    }

    // Local capture: an aligned near pair binds, while an axially or radially
    // distant pair does not.
    {
        gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
        auto model = make_model(output, rng, seed, K);
        configure(*model, 1, 1e12, 0.0);
        place_pair(*model, 0.15, 0.11);
        model->update_system();
        assert(model->myomesin_bridges.size() == 1);
        model->myomesin_bridges.clear();
        place_pair(*model, 0.40, 0.11);
        model->update_system();
        assert(model->myomesin_bridges.empty());
        place_pair(*model, 0.10, 0.50);
        model->update_system();
        assert(model->myomesin_bridges.empty());
        gsl_rng_free(rng);
    }

    // Finite valency: three mutually local thick filaments with one site each
    // can form only one bridge, never the old three-edge complete graph.
    {
        gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
        auto model = make_model(output, rng, seed, K, 3);
        configure(*model, 1, 1e12, 0.0);
        model->myosin.center[0] = {0.00, 0.00, 0.0};
        model->myosin.center[1] = {0.02, 0.11, 0.0};
        model->myosin.center[2] = {-0.02, 0.055, 0.0952628};
        model->myosin.update_endpoints();
        model->update_system();
        assert(model->myomesin_bridges.size() == 1);
        for (int i = 0; i < 3; ++i) {
            assert(model->myomesin_degree(i) <= 1);
        }
        gsl_rng_free(rng);
    }

    // Turnover changes only bridge occupancy; it does not alter AM ownership.
    {
        gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
        auto model = make_model(output, rng, seed, K);
        configure(*model, 1, 0.0, 1e12);
        place_pair(*model, 0.15, 0.11);
        assert(model->prebind_myomesin_bridge(0, 1));
        model->update_system();
        assert(model->myomesin_bridges.empty());
        gsl_rng_free(rng);
    }

    // k_mm=0 remains an exact opt-out.
    {
        gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
        auto model = make_model(output, rng, seed, 0.0);
        configure(*model, 1, 1e12, 0.0);
        place_pair(*model, 0.15, 0.11);
        model->update_system();
        assert(model->myomesin_bridges.empty());
        assert(model->myosin.force[0].norm() < 1e-12);
        gsl_rng_free(rng);
    }

    std::remove(output.c_str());
    return 0;
}
