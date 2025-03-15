#include "langevin.h"
#include <cstdio>
#include <cmath>
#include <vector>

// Define the global mutex.
std::mutex save_mutex;

//---------------------------------------------------------------------
// Constructor: Initializes the Langevin object and, if resume is true,
// loads the previous state and prints details. Otherwise, creates a new file.
//---------------------------------------------------------------------
Langevin::Langevin(Sarcomere& model0, double& beta0, double& dt0, double& D0_actin, double& D0_myosin,
                   int& update_myosin_every0, int& update_dt_every0, int& save_every0, bool& resume)
    : model(model0), beta(beta0), dt(dt0), D_actin(D0_actin), D_myosin(D0_myosin),
      update_myosin_every(update_myosin_every0), update_dt_every(update_dt_every0),
      save_every(save_every0)
{
    if (resume) {
        model.load_state();
        printf("Resuming from file %s\n", model.filename.c_str());
        // Print myosin details.
        for (int i = 0; i < model.myosin.n; i++) {
            printf("Myosin %d: %f %f %f\n", i, model.myosin.center[i].x, model.myosin.center[i].y, model.myosin.center[i].z);
            printf("Myosin endpoints: (%f %f %f), (%f %f %f)\n",
                   model.myosin.left_end[i].x, model.myosin.left_end[i].y, model.myosin.left_end[i].z,
                   model.myosin.right_end[i].x, model.myosin.right_end[i].y, model.myosin.right_end[i].z);
        }
        // Print actin details.
        for (int i = 0; i < model.actin.n; i++) {
            printf("Actin %d: %f %f %f\n", i, model.actin.center[i].x, model.actin.center[i].y, model.actin.center[i].z);
            printf("Actin endpoints: (%f %f %f), (%f %f %f)\n",
                   model.actin.left_end[i].x, model.actin.left_end[i].y, model.actin.left_end[i].z,
                   model.actin.right_end[i].x, model.actin.right_end[i].y, model.actin.right_end[i].z);
        }
    } else {
        model.new_file();
    }
}

//---------------------------------------------------------------------
// Destructor
//---------------------------------------------------------------------
Langevin::~Langevin() {
    // No dynamic resources need explicit deallocation.
}

//---------------------------------------------------------------------
// run_langevin: Runs the simulation for nsteps, periodically saving the
// state and taking sample steps.
//---------------------------------------------------------------------
void Langevin::run_langevin(int nsteps, gsl_rng* rng, int& fix_myosin) {
    double start, end;
    for (int i = 0; i < nsteps; i++) {
        if (i % save_every == 0) {
            std::cout << "Step " << i << std::endl;
            // Optionally protect saving with the mutex:
            // std::lock_guard<std::mutex> lock(save_mutex);
            model.save_state();
            start = omp_get_wtime();
        }
        model.update_system();
        sample_step(dt, D_actin, D_myosin, rng, fix_myosin);
        if (i % save_every == 0) {
            end = omp_get_wtime();
            printf("Step %d took %f seconds\n", i, end - start);
        }
    }
}

void Langevin::volume_exclusion(int nsteps, gsl_rng* rng, int& fix_myosin) {
    double start, end;
    for (int i = 0; i < nsteps; i++) {
        if (i % save_every == 0) {
            std::cout << "Step " << i << std::endl;
            // Optionally protect saving with the mutex:
            // std::lock_guard<std::mutex> lock(save_mutex);
            start = omp_get_wtime();
        }
        model.update_system_sterics_only();
        sample_step(dt, D_actin, D_myosin, rng, fix_myosin);
        if (i % save_every == 0) {
            end = omp_get_wtime();
            printf("Step %d took %f seconds\n", i, end - start);
        }
    }
}

//---------------------------------------------------------------------
// sample_step: Performs a single Langevin dynamics step by updating the 
// system, generating noise, and displacing myosin and actin particles.
//---------------------------------------------------------------------
void Langevin::sample_step(double& dt, double& D_actin, double& D_myosin, gsl_rng* rng, int& fix_myosin) {
    model.update_system();

    // Generate noise for both myosin and actin particles.
    int n_randns = (model.myosin.n + model.actin.n) * 5;
    std::vector<double> noise(n_randns);
    for (int i = 0; i < n_randns; i++) {
        noise[i] = gsl_rng_uniform(rng) * 2 - 1;
    }

    int n_acc_randns = model.myosin.n + model.actin.n;
    std::vector<double> acc_rand(n_acc_randns);
    for (int i = 0; i < n_acc_randns; i++) {
        acc_rand[i] = gsl_rng_uniform(rng);
    }

    int offset = model.myosin.n * 5;
    double D = D_myosin;
    // Update myosin particles.
    for (int i = fix_myosin; i < model.myosin.n; i++) {
        double dx = model.myosin.force[i].x * beta * D * dt +
                    model.myosin.velocity[i].x * dt +
                    sqrt(2 * D * dt) * noise[i * 5];
        double dy = model.myosin.force[i].y * beta * D * dt +
                    model.myosin.velocity[i].y * dt +
                    sqrt(2 * D * dt) * noise[i * 5 + 1];
        double dz = model.myosin.force[i].z * beta * D * dt +
                    model.myosin.velocity[i].z * dt +
                    sqrt(2 * D * dt) * noise[i * 5 + 2];
        double dtheta = model.myosin.angular_force[i][0] * beta * D * dt +
                        sqrt(2 * D * dt) * noise[i * 5 + 3] * M_PI / 5;
        double dphi = model.myosin.angular_force[i][1] * beta * D * dt +
                        sqrt(2 * D * dt) * noise[i * 5 + 4] * M_PI / 5;
        model.myosin.displace(i, dx, dy, dz, dtheta, dphi);
    }
    D = D_actin;
    // Update actin particles.
    for (int i = 0; i < model.actin.n; i++) {
        double dx = model.actin.force[i].x * beta * D * dt +
                    model.actin.velocity[i].x * dt +
                    sqrt(2 * D * dt) * noise[offset + i * 5];
        double dy = model.actin.force[i].y * beta * D * dt +
                    model.actin.velocity[i].y * dt +
                    sqrt(2 * D * dt) * noise[offset + i * 5 + 1];
        double dz = model.actin.force[i].z * beta * D * dt +
                    model.actin.velocity[i].z * dt +
                    sqrt(2 * D * dt) * noise[offset + i * 5 + 2];
        double dtheta = model.actin.angular_force[i][0] * beta * D * dt +
                        sqrt(2 * D * dt) * noise[offset + i * 5 + 3] * M_PI / 5;
        double dphi = model.actin.angular_force[i][1] * beta * D * dt +
                        sqrt(2 * D * dt) * noise[offset + i * 5 + 4] * M_PI / 5;

        // (Optional debug printing if the displacement is large.)
        if (dx > 0.04 || dy > 0.04) {
            printf("actin %d displacement too large \n", i);
        }
        model.actin.displace(i, dx, dy, dz, dtheta, dphi);
    }
}
