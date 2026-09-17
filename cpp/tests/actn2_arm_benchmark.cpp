#include "actn2_arm_kinetics.h"
#include "sarcomere.h"

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>

int main(int argc, char** argv) {
    const std::string mode = argc > 1 ? argv[1] : "legacy";
    const bool tether_model = mode == "tether";
    const bool arm_model = mode == "arm" || tether_model;
    const int steps = argc > 2 ? std::max(1, std::atoi(argv[2])) : 10000;
    int n_actins = 2;
    int n_myosins = 1;
    std::vector<double> box{4.0, 4.0, 4.0};
    double actin_length = 1.0;
    double myosin_length = 1.5;
    double myosin_radius = 0.05;
    double am_cutoff = 0.05;
    double am_optimal = 0.055;
    double aa_cutoff = 0.08;
    double aa_optimal = 0.05;
    double zero = 0.0;
    double base_lifetime = 1e9;
    int seed = 1729;
    int fixed_myosins = 1;
    double dt = 1e-5;
    bool directional = true;
    std::string output = "/tmp/actn2_arm_benchmark_unused.h5";
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, seed);
    auto model = std::make_unique<Sarcomere>(
        n_actins, n_myosins, box, actin_length, myosin_length,
        myosin_radius, am_cutoff, am_optimal, aa_cutoff, aa_optimal,
        zero, zero, base_lifetime, zero, zero, zero, zero, zero, zero,
        zero, zero, zero, output, rng, seed, fixed_myosins, dt, zero, zero,
        zero, directional, 24, zero, zero, zero, zero);
    model->set_periodicity({false, false, false});
    model->dynamic_am = false;
    model->aa_local_multilink = true;
    model->aa_prebundle_links_per_pair = 1;
    model->aa_neighbor_identity_enabled = false;
    model->aa_arm_state_model = arm_model;
    model->aa_arm_open_rate_weak = 0.0;
    model->aa_arm_open_rate_mature = 0.0;
    model->aa_tether_off_rate = 0.0;
    model->aa_tether_rebind_rate = 0.0;
    model->aa_link_rupture_length = 0.0;
    model->aa_occlusion_rupture_penetration = 0.0;
    model->sterics_enabled = false;

    model->actin.center[0] = {0.0, -0.025, 0.0};
    model->actin.center[1] = {0.0, 0.025, 0.0};
    model->actin.direction[0] = {1.0, 0.0, 0.0};
    model->actin.direction[1] = {-1.0, 0.0, 0.0};
    model->actin.update_endpoints();
    // Perpendicular thick filament: it is a geometric ACTN2 contact candidate
    // but not an antiparallel A-M binding partner.
    model->myosin.center[0] = {0.0, 0.0, 0.040};
    model->myosin.direction[0] = {0.0, 0.0, 1.0};
    model->myosin.update_endpoints();
    model->prebind_local_aa_pair(0, 1, 1);
    auto& link = model->aa_links[model->aa_active_link_slots.front()];
    link.status = static_cast<Sarcomere::BondFlag>(2);
    link.arm_state = tether_model
        ? actn2_arm_kinetics::State::tethered_to_i
        : actn2_arm_kinetics::State::mature_bridge;

    for (int warmup = 0; warmup < 100; ++warmup) model->update_system();
    const auto start = std::chrono::steady_clock::now();
    for (int step = 0; step < steps; ++step) model->update_system();
    const auto elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    size_t active_links = 0;
    double contact_force = 0.0;
    for (const auto& current : model->aa_links) {
        if (!current.active) continue;
        ++active_links;
        contact_force = std::max(contact_force, current.contact_force);
    }
    std::cout << mode
              << " steps=" << steps << " seconds=" << elapsed
              << " microseconds_per_step=" << 1e6 * elapsed / steps
              << " active_links=" << active_links
              << " contact_force_pN=" << contact_force
              << '\n';
    gsl_rng_free(rng);
    return 0;
}
