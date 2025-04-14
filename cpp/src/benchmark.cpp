#include <benchmark/benchmark.h>
#include "langevin.h"
#include "sarcomere.h"
#include "components.h"
#include "cxxopts.hpp"  // if needed, though not used in benchmark
#include <vector>
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
        v_am = 5;
        n_actins = 800;
        n_myosins = 400;
        Lx = 12;
        Ly = 6;
        Lz = 3.5;
        actin_length = 1;
        myosin_length = 1.5;
        myosin_radius = 0.2;
        myosin_radius_ratio = 0.75;
        crosslinker_length = 0.06;
        resume = false;
        directional = true;
        n_fixed_myosins = 0;
        filename = "traj.h5";
        init_struc = "random"; // or "sarcomere", "partial", etc.

        // Create the simulation box.
        std::vector<double> box = {Lx, Ly, Lz};

        // Allocate the GSL RNG.
        rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng, seed);

        // Create the Sarcomere model.
        double diff_coeff_ratio = actin_diff_coeff_trans/myosin_diff_coeff_trans;
        model = new Sarcomere(n_actins, n_myosins, box, actin_length, myosin_length,
            myosin_radius, myosin_radius_ratio, crosslinker_length, k_on, k_off,
            base_lifetime, lifetime_coeff, diff_coeff_ratio,
              k_aa, kappa_aa, k_am, kappa_am, v_am,
            filename,rng, seed, n_fixed_myosins, dt, directional);

        // Create the Langevin simulation instance.
        sim = new Langevin(*model, beta, dt, actin_diff_coeff_trans,actin_diff_coeff_rot,
             myosin_diff_coeff_trans, myosin_diff_coeff_rot, save_every, resume);

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
              myosin_diff_coeff_rot, myosin_radius_ratio;
    double k_on, k_off, base_lifetime, lifetime_coeff;
    double k_aa, kappa_aa, k_am, kappa_am, v_am;
    double Lx, Ly, Lz, actin_length, myosin_length, myosin_radius, crosslinker_length;
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
