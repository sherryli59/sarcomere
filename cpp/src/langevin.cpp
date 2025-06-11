#include "langevin.h"
#include <cstdio>
#include <cmath>
#include <vector>

using vec = utils::vec;
// Define the global mutex.
std::mutex save_mutex;

//---------------------------------------------------------------------
// Constructor: Initializes the Langevin object and, if resume is true,
// loads the previous state and prints details. Otherwise, creates a new file.
//---------------------------------------------------------------------
Langevin::Langevin(Sarcomere& model0, double& beta0, double& dt0, double& D0_actin_trans,
    double& D0_actin_rot, double& D0_myosin_trans, double& D0_myosin_rot, int& save_every0, bool& resume)
        : model(model0), beta(beta0), dt(dt0), D_actin_trans(D0_actin_trans),
        D_actin_rot(D0_actin_rot),D_myosin_trans(D0_myosin_trans), D_myosin_rot(D0_myosin_rot),
        save_every(save_every0)
{
    if (resume) {
        int n_frames;
        model.load_state(n_frames);
        start_step = (n_frames-1) * save_every + 1;
        printf("Resuming from file %s at step %d\n", model.filename.c_str(), start_step);
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
        start_step = 0;
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
        sample_step(dt, rng, fix_myosin);
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
        sample_step(dt, rng, fix_myosin);
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
void Langevin::sample_step(double& dt, gsl_rng* rng, int& fix_myosin) {
    model.update_system();

    // Generate noise for both myosin and actin particles.
    int n_randns = (model.myosin.n + model.actin.n) * 6;
    std::vector<double> noise(n_randns);
    for (int i = 0; i < n_randns; i++) {
        noise[i] = gsl_ran_gaussian(rng, 1.0);
    }

    int n_acc_randns = model.myosin.n + model.actin.n;
    std::vector<double> acc_rand(n_acc_randns);
    for (int i = 0; i < n_acc_randns; i++) {
        acc_rand[i] = gsl_rng_uniform(rng);
    }

    int offset = model.myosin.n * 6;
    double D = D_myosin_trans;
    double D_rot = D_myosin_rot;
    // Update myosin particles.
    for (int i = fix_myosin; i < model.myosin.n; i++) {
        double dx = model.myosin.force[i].x * beta * D * dt +
                    model.myosin.velocity[i].x * dt +
                    sqrt(2 * D * dt) * noise[i * 6];
        double dy = model.myosin.force[i].y * beta * D * dt +
                    model.myosin.velocity[i].y * dt +
                    sqrt(2 * D * dt) * noise[i * 6 + 1];
        double dz = model.myosin.force[i].z * beta * D * dt +
                    model.myosin.velocity[i].z * dt +
                    sqrt(2 * D * dt) * noise[i * 6 + 2];
        // double dtheta = model.myosin.angular_force[i][0] * beta * D_rot * dt +
        //                 sqrt(2 * D_rot * dt) * noise[i * 5 + 3] * M_PI;
        // double dphi = model.myosin.angular_force[i][1] * beta * D_rot * dt +
        //                 sqrt(2 * D_rot * dt) * noise[i * 5 + 4] * M_PI;
        model.myosin.displace(i, dx, dy, dz);
        vec rot_noise={noise[i * 6 + 3], noise[i * 6 + 4], noise[i * 6 + 5]};
        vec delta_u = sqrt(2 * D_rot * dt) * rot_noise.cross(model.myosin.direction[i])
                        + dt * model.myosin.torque[i].cross(model.myosin.direction[i]);
        model.myosin.direction[i] += delta_u;
        model.myosin.direction[i].normalize();
    }
    // Update actin particles.
    for (int i = 0; i < model.actin.n; i++) {
        if (model.actin.cb_strength[i] > 1e-3){
            D = D_myosin_trans;
            D_rot = D_myosin_rot;
        }
        else{
            D = D_actin_trans;
            D_rot = D_actin_rot;
        }
        double dx = model.actin.force[i].x * beta * D * dt +
                    model.actin.velocity[i].x * dt +
                    sqrt(2 * D * dt) * noise[offset + i * 6];
        double dy = model.actin.force[i].y * beta * D * dt +
                    model.actin.velocity[i].y * dt +
                    sqrt(2 * D * dt) * noise[offset + i * 6 + 1];
        double dz = model.actin.force[i].z * beta * D * dt +
                    model.actin.velocity[i].z * dt +
                    sqrt(2 * D * dt) * noise[offset + i * 6 + 2];
        // double dtheta = model.actin.angular_force[i][0] * beta * D_rot * dt +
        //                 sqrt(2 * D_rot * dt) * noise[offset + i * 5 + 3] * M_PI;
        // double dphi = model.actin.angular_force[i][1] * beta * D_rot * dt +
        //                 sqrt(2 * D_rot * dt) * noise[offset + i * 5 + 4] * M_PI;

        // (Optional debug printing if the displacement is large.)
        if (dx > 0.04 || dy > 0.04) {
            printf("actin %d displacement too large \n", i);
        }
        model.actin.displace(i, dx, dy, dz);
        vec rot_noise={noise[offset + i * 6 + 3], noise[offset + i * 6 + 4], noise[offset + i * 6 + 5]};
        vec delta_u = sqrt(2 * D_rot * dt) * rot_noise.cross(model.actin.direction[i])
                        + dt * model.actin.torque[i].cross(model.actin.direction[i]);
        model.actin.direction[i] += delta_u;
        model.actin.direction[i].normalize();
    }
}
