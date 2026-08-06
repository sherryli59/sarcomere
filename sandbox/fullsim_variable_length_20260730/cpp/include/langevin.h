#ifndef LANGEVIN_H
#define LANGEVIN_H

#include "sarcomere.h"
#include "utils.h"
#include "h5_utils.h"
#include "rod_mobility.h"
#include "axial_cell.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <omp.h>
#include <mutex>

// Declare the global mutex for saving (defined in the source file)
extern std::mutex save_mutex;

class Langevin {
public:
    // Constructor.
    // Note: Parameters are passed by reference as in the original code.
    Langevin(Sarcomere& model0, double& beta0, double& dt0, double& D0_actin_trans,
        double& D0_actin_rot, double& D0_myosin_trans, double& D0_myosin_rot,
        int& save_every0, bool& resume, bool is3D,
        double max_actin_displacement, double max_myosin_displacement,
        double max_actin_rotation, double max_myosin_rotation,
        int resume_frame = -1);
    // Destructor.
    ~Langevin();

    // Run the Langevin simulation for a given number of steps.
    void run_langevin(int nsteps, gsl_rng* rng, int& fix_myosin);
    void volume_exclusion(int nsteps, gsl_rng* rng, int& fix_myosin,
                          bool save_frames = true);
    // Take a single simulation step.
    void sample_step(double& dt, gsl_rng* rng, int& fix_myosin);
    void set_rod_mobility(const rod_mobility::Diffusion& actin_mobility,
                          const rod_mobility::Diffusion& myosin_mobility);
    void set_isometric_actin_axial(bool enabled) {
        constrain_actin_axial = enabled;
    }
    void set_thermal_noise_scale(double scale) {
        thermal_noise_scale = scale;
        rotational_noise_scale = scale;
    }
    void set_rotational_noise_scale(double scale) {
        rotational_noise_scale = scale;
    }
    void set_axial_cell(const axial_cell::Config& config) {
        axial_cell_config = config;
    }

    // Data members.
    double dt, beta, D_actin_trans, D_actin_rot, D_myosin_trans, D_myosin_rot;
    double max_actin_displacement, max_myosin_displacement;
    double max_actin_rotation, max_myosin_rotation;
    rod_mobility::Diffusion actin_rod_mobility, myosin_rod_mobility;
    bool physical_rod_mobility = false;
    bool constrain_actin_axial = false;
    double thermal_noise_scale = 1.0;
    double rotational_noise_scale = 1.0;
    bool is3D;
    axial_cell::Config axial_cell_config;
    axial_cell::Diagnostics last_axial_cell;
    bool skip_initial_save;
    int save_every, start_step, loaded_frame_index;
    // When true, interpret `nsteps` argument to run_* as the total target
    // number of steps (i.e., inclusive of previously-run steps) while
    // resuming. Set in constructor when resume mode is active.
    bool interpret_nsteps_as_total_when_resuming = false;
    Sarcomere& model;
};

#endif // LANGEVIN_H
