#include <benchmark/benchmark.h>
#include "langevin.h"
#include "sarcomere.h"
#include "components.h"
#include "cxxopts.hpp"  // if needed, though not used in benchmark
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// Fixture for benchmarking run_langevin.
class RunLangevinBenchmark : public benchmark::Fixture {
public:
    void SetUp(const ::benchmark::State& state) override {
        // Set up simulation parameters (use smaller nsteps for benchmarking).
        nsteps = 10;          // Use a moderate number of steps for timing.
        seed = 0;
        dt = 0.00001;
        beta = 241.0;
        actin_diff_coeff_trans = 1;
        actin_diff_coeff_rot = 2;
        myosin_diff_coeff_trans = 0.05;
        myosin_diff_coeff_rot = 0.05;
        save_every = 200;
        k_on = 100;
        k_off = 1;
        base_lifetime = 0.001;
        lifetime_coeff = 0.4;
        k_aa = 300;
        kappa_aa = 50;
        k_am = 50;
        kappa_am = 50;
        k_mm = 0;
        v_am = 5;
        n_actins = 800;
        n_myosins = 400;
        Lx = 12;
        Ly = 6;
        Lz = 3.5;
        actin_length = 1;
        myosin_length = 1.5;
        myosin_radius = 0.015;
        am_cutoff = 0.05;
        am_optimal = 0.03;
        aa_cutoff = 0.05;
        aa_optimal = 0.03;
        resume = false;
        directional = true;
        n_fixed_myosins = 0;
        filename = "traj.h5";
        init_struc = "random"; // or "sarcomere", "partial", etc.
        double tau_rec = 0.0;
        double titin_k = 0.0;
        double titin_rest_length = 0.5 * myosin_length + (2.0 / 3.0) * actin_length;

        // Create the simulation box.
        std::vector<double> box = {Lx, Ly, Lz};

        // Allocate the GSL RNG.
        rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng, seed);

        // Create the Sarcomere model.
        double diff_coeff_ratio = actin_diff_coeff_trans/myosin_diff_coeff_trans;
        double max_actin_displacement = 0.01;
        double max_myosin_displacement = 0.01;
        auto compute_max_force = [&](double max_disp, double D_coeff) {
            if (max_disp <= 0.0 || D_coeff <= 0.0 || !std::isfinite(beta) || !std::isfinite(dt)) {
                return std::numeric_limits<double>::infinity();
            }
            double val = max_disp / (beta * D_coeff * dt);
            if (!(val > 0.0)) {
                return std::numeric_limits<double>::infinity();
            }
            return val;
        };
        double max_actin_force = compute_max_force(max_actin_displacement, actin_diff_coeff_trans);
        double max_myosin_force = compute_max_force(max_myosin_displacement, myosin_diff_coeff_trans);
        model = new Sarcomere(n_actins, n_myosins, box, actin_length, myosin_length,
            myosin_radius, am_cutoff, am_optimal, aa_cutoff, aa_optimal,
            k_on, k_off,
            base_lifetime, lifetime_coeff, diff_coeff_ratio,
              k_aa, kappa_aa, k_am, kappa_am, k_mm, v_am,
            filename,rng, seed, n_fixed_myosins, dt, tau_rec,
            titin_k, titin_rest_length, directional, 5, max_actin_force, max_myosin_force,
            std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());

        // Create the Langevin simulation instance.
        bool is3D = true;
        const double max_actin_rotation = 0.01;
        const double max_myosin_rotation = 0.01;
        sim = new Langevin(*model, beta, dt, actin_diff_coeff_trans,actin_diff_coeff_rot,
             myosin_diff_coeff_trans, myosin_diff_coeff_rot, save_every, resume, is3D,
             max_actin_displacement, max_myosin_displacement,
             max_actin_rotation, max_myosin_rotation, -1);

        // Set up the initial structure.
        if (!resume) {
            if (init_struc == "sarcomere") {
                sim->model.sarcomeric_structure();
            } else if (init_struc == "partial") {
                sim->model.partial_fix(n_fixed_myosins);
            } else if (init_struc == "cb") {
                sim->model.cb();
            }
        }

        // Perform volume exclusion as in the simulation.
        sim->volume_exclusion(1, rng, n_fixed_myosins);
    }

    void TearDown(const ::benchmark::State& state) override {
        gsl_rng_free(rng);
        delete sim;
        delete model;
    }

protected:
    // Simulation parameters.
    int nsteps, seed, save_every;
    int n_actins, n_myosins, n_fixed_myosins;
    double dt, beta, actin_diff_coeff_trans, actin_diff_coeff_rot, myosin_diff_coeff_trans,
              myosin_diff_coeff_rot;
    double k_on, k_off, base_lifetime, lifetime_coeff;
    double k_aa, kappa_aa, k_am, kappa_am, k_mm, v_am;
    double Lx, Ly, Lz, actin_length, myosin_length, myosin_radius;
    double am_cutoff, am_optimal, aa_cutoff, aa_optimal;
    bool resume, directional;
    std::string filename, init_struc;

    gsl_rng* rng;
    Sarcomere* model;
    Langevin* sim;
};

// Benchmark for the run_langevin function.
BENCHMARK_DEFINE_F(RunLangevinBenchmark, RunLangevin)(benchmark::State& state) {
    for (auto _ : state) {
        // Time the run_langevin call.
        sim->run_langevin(nsteps, rng, n_fixed_myosins);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK_REGISTER_F(RunLangevinBenchmark, RunLangevin)
    ->Iterations(1);  // Adjust the iterations count as needed for your workload.

BENCHMARK_MAIN();
