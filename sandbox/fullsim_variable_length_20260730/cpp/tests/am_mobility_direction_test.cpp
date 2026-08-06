#include "langevin.h"
#include "rod_mobility.h"
#include "sarcomere.h"

#include <cassert>
#include <cmath>
#include <string>

namespace {

bool near(double first, double second, double tolerance = 1e-11) {
    return std::abs(first - second) <= tolerance;
}

rod_mobility::Diffusion mobility(double parallel, double perpendicular) {
    rod_mobility::Diffusion result;
    result.enabled = true;
    result.length = 1.0;
    result.diameter = 0.01;
    result.viscosity = 0.001;
    result.parallel = parallel;
    result.perpendicular = perpendicular;
    result.tumble = 0.1;
    result.roll = 0.1;
    return result;
}

}  // namespace

int main() {
    int n_actins = 1;
    int n_myosins = 1;
    std::vector<double> box{4.0, 4.0, 4.0};
    double actin_length = 1.0;
    double myosin_length = 1.5;
    double myosin_radius = 0.025;
    double am_cutoff = 0.1;
    double am_optimal = 0.05;
    double aa_cutoff = 0.1;
    double aa_optimal = 0.05;
    double k_on = 0.0;
    double k_off = 0.0;
    double base_lifetime = 0.0;
    double directional_lifetime = 0.0;
    double lifetime_coeff = 0.0;
    double diff_ratio = 1.0;
    double zero = 0.0;
    std::string filename = "/tmp/am_mobility_direction_test.h5";
    int seed = 731;
    int fix_myosin = 0;
    double dt = 0.01;
    bool directional = true;

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, seed);
    Sarcomere model(
        n_actins, n_myosins, box, actin_length, myosin_length,
        myosin_radius, am_cutoff, am_optimal, aa_cutoff, aa_optimal,
        k_on, k_off, base_lifetime, directional_lifetime, lifetime_coeff,
        diff_ratio, zero, zero, zero, zero, zero, zero,
        filename, rng, seed, fix_myosin, dt, 0.0, 0.0, 0.0,
        directional, 1, zero, zero, zero, zero);
    model.sterics_enabled = false;
    model.actin.direction[0] = {1.0, 0.0, 0.0};
    model.myosin.direction[0] = {1.0, 0.0, 0.0};

    const auto actin_mobility = mobility(0.30, 0.20);
    const auto myosin_mobility = mobility(0.10, 0.07);
    const double temperature = 310.0;
    double beta = 2.0;
    model.set_physical_mobility(
        temperature, beta, actin_mobility, myosin_mobility);

    double scalar_actin_translation = 0.0;
    double scalar_actin_rotation = 0.0;
    double scalar_myosin_translation = 0.0;
    double scalar_myosin_rotation = 0.0;
    int save_every = 1000;
    bool resume = false;
    double maximum_rotation = 0.01;
    Langevin integrator(
        model, beta, dt,
        scalar_actin_translation, scalar_actin_rotation,
        scalar_myosin_translation, scalar_myosin_rotation,
        save_every, resume, true, 0.0, 0.0,
        maximum_rotation, maximum_rotation);
    integrator.set_rod_mobility(actin_mobility, myosin_mobility);
    integrator.set_thermal_noise_scale(0.0);

    auto reset = [&]() {
        model.actin.center[0] = {-0.30, 0.0, 0.0};
        model.myosin.center[0] = {0.0, 0.0, 0.0};
        model.actin.direction[0] = {1.0, 0.0, 0.0};
        model.myosin.direction[0] = {1.0, 0.0, 0.0};
        model.actin.force[0] = {1.0, 0.0, 0.0};
        model.myosin.force[0] = {-1.0, 0.0, 0.0};
        model.actin.velocity[0] = {0.0, 0.0, 0.0};
        model.myosin.velocity[0] = {0.0, 0.0, 0.0};
        model.actin.torque[0] = {0.0, 0.0, 0.0};
        model.myosin.torque[0] = {0.0, 0.0, 0.0};
        model.actin.update_endpoints();
        model.myosin.update_endpoints();
    };

    // Negative-half actin: +x points inward toward the M-line.
    reset();
    fix_myosin = 1;
    integrator.set_isometric_actin_axial(false);
    const double fixed_myo_actin_x = model.actin.center[0].x;
    integrator.sample_step(dt, rng, fix_myosin);
    assert(model.actin.center[0].x > fixed_myo_actin_x);
    assert(near(model.myosin.center[0].x, 0.0));

    reset();
    fix_myosin = 0;
    integrator.set_isometric_actin_axial(true);
    integrator.sample_step(dt, rng, fix_myosin);
    assert(near(model.actin.center[0].x, -0.30));
    assert(model.myosin.center[0].x < 0.0);

    reset();
    integrator.set_isometric_actin_axial(false);
    const double initial_relative =
        model.actin.center[0].x - model.myosin.center[0].x;
    const double initial_drag_center =
        model.actin.center[0].x / actin_mobility.parallel +
        model.myosin.center[0].x / myosin_mobility.parallel;
    integrator.sample_step(dt, rng, fix_myosin);
    const double final_relative =
        model.actin.center[0].x - model.myosin.center[0].x;
    const double final_drag_center =
        model.actin.center[0].x / actin_mobility.parallel +
        model.myosin.center[0].x / myosin_mobility.parallel;
    assert(final_relative > initial_relative);
    assert(model.actin.center[0].x > -0.30);
    assert(model.myosin.center[0].x < 0.0);
    assert(near(final_drag_center, initial_drag_center));

    // The realized director increment includes deterministic torque and noise.
    // A large increment must be rescaled to displacement_slack (1.3) times
    // the configured trust radius, and the event must be observable in HDF5.
    reset();
    model.actin.force[0] = {0.0, 0.0, 0.0};
    model.myosin.force[0] = {0.0, 0.0, 0.0};
    model.actin.torque[0] = {0.0, 100.0, 0.0};
    fix_myosin = 1;
    integrator.sample_step(dt, rng, fix_myosin);
    const double actin_rotation = std::acos(std::clamp(
        static_cast<vec>(model.actin.direction[0]).dot({1.0, 0.0, 0.0}),
        -1.0, 1.0));
    assert(near(actin_rotation, 1.3 * maximum_rotation, 1e-10));
    assert(near(model.actin["cumulative_rotation_cap_events"][0], 1.0));

    reset();
    model.actin.force[0] = {0.0, 0.0, 0.0};
    model.myosin.force[0] = {0.0, 0.0, 0.0};
    model.myosin.torque[0] = {0.0, 100.0, 0.0};
    fix_myosin = 0;
    integrator.sample_step(dt, rng, fix_myosin);
    const double myosin_rotation = std::acos(std::clamp(
        static_cast<vec>(model.myosin.direction[0]).dot({1.0, 0.0, 0.0}),
        -1.0, 1.0));
    assert(near(myosin_rotation, 1.3 * maximum_rotation, 1e-10));
    assert(near(model.myosin["cumulative_rotation_cap_events"][0], 1.0));

    gsl_rng_free(rng);
    return 0;
}
