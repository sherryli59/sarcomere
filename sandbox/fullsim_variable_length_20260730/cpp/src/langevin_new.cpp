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
    double& D0_actin_rot, double& D0_myosin_trans, double& D0_myosin_rot, int& save_every0, bool& resume, bool is3D,
    double max_actin_displacement0, double max_myosin_displacement0,
    double max_actin_rotation0, double max_myosin_rotation0,
    int resume_frame)
        : model(model0), beta(beta0), dt(dt0), D_actin_trans(D0_actin_trans),
        D_actin_rot(D0_actin_rot),D_myosin_trans(D0_myosin_trans), D_myosin_rot(D0_myosin_rot),
        save_every(save_every0), is3D(is3D),
        max_actin_displacement(max_actin_displacement0), max_myosin_displacement(max_myosin_displacement0),
        max_actin_rotation(max_actin_rotation0), max_myosin_rotation(max_myosin_rotation0),
        loaded_frame_index(-1)
{
    if (resume) {
        int n_frames;
        int requested_frame = resume_frame;
        int frame_idx = model.load_state(n_frames, requested_frame);
        loaded_frame_index = frame_idx;
        start_step = frame_idx * save_every + 1;
        printf("Resuming from file %s at step %d (frame %d of %d)\n",
               model.filename.c_str(), start_step, frame_idx, n_frames);
        if (requested_frame >= 0 && frame_idx != requested_frame) {
            printf("Requested resume frame %d adjusted to %d (available frames: %d)\n",
                   requested_frame, frame_idx, n_frames);
        }
        // Print myosin details.
        // for (int i = 0; i < model.myosin.n; i++) {
        //     printf("Myosin %d: %f %f %f\n", i, model.myosin.center[i].x, model.myosin.center[i].y, model.myosin.center[i].z);
        //     printf("Myosin endpoints: (%f %f %f), (%f %f %f)\n",
        //            model.myosin.left_end[i].x, model.myosin.left_end[i].y, model.myosin.left_end[i].z,
        //            model.myosin.right_end[i].x, model.myosin.right_end[i].y, model.myosin.right_end[i].z);
        // }
        // // Print actin details.
        // for (int i = 0; i < model.actin.n; i++) {
        //     printf("Actin %d: %f %f %f\n", i, model.actin.center[i].x, model.actin.center[i].y, model.actin.center[i].z);
        //     printf("Actin endpoints: (%f %f %f), (%f %f %f)\n",
        //            model.actin.left_end[i].x, model.actin.left_end[i].y, model.actin.left_end[i].z,
        //            model.actin.right_end[i].x, model.actin.right_end[i].y, model.actin.right_end[i].z);
        // }
    } else {
        start_step = 0;
        loaded_frame_index = -1;
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
    int end_step = start_step + nsteps;
    for (int step = start_step; step < end_step; ++step) {
        if (step % save_every == 0) {
            std::cout << "Step " << step << std::endl;
            // Optionally protect saving with the mutex:
            // std::lock_guard<std::mutex> lock(save_mutex);
            model.save_state();
            start = omp_get_wtime();
        }
        model.update_system();
        sample_step(dt, rng, fix_myosin);
        if (step % save_every == 0) {
            end = omp_get_wtime();
            printf("Step %d took %f seconds\n", step, end - start);
            //model.debug_cb_stats();
        } 
    }
    start_step = end_step;
}

void Langevin::volume_exclusion(int nsteps, gsl_rng* rng, int& fix_myosin) {
    double start, end;
    int end_step = start_step + nsteps;
    for (int step = start_step; step < end_step; ++step) {
        if (step % save_every == 0) {
            std::cout << "Step " << step << std::endl;
            // Optionally protect saving with the mutex:
            // std::lock_guard<std::mutex> lock(save_mutex);
            start = omp_get_wtime();
        }
        model.update_system_sterics_only();
        sample_step(dt, rng, fix_myosin);
        if (step % save_every == 0) {
            end = omp_get_wtime();
            printf("Step %d took %f seconds\n", step, end - start);
        }
    }
    start_step = end_step;
}

//---------------------------------------------------------------------
// sample_step: Performs a single Langevin dynamics step by updating the 
// system, generating noise, and displacing myosin and actin particles.
//---------------------------------------------------------------------
void Langevin::sample_step(double& dt, gsl_rng* rng, int& fix_myosin) {
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
    const double displacement_slack = 1.3;
    const std::array<bool,3>& periodic = model.is_periodic;
    for (int i = fix_myosin; i < model.myosin.n; i++) {
        if (!is3D) {
            model.myosin.force[i].z = 0;
            model.myosin.velocity[i].z = 0;
            model.myosin.torque[i].z = 0;
        }
        vec delta_pos{
            model.myosin.force[i].x * beta * D * dt +
                model.myosin.velocity[i].x * dt +
                std::sqrt(2 * D * dt) * noise[i * 6],
            model.myosin.force[i].y * beta * D * dt +
                model.myosin.velocity[i].y * dt +
                std::sqrt(2 * D * dt) * noise[i * 6 + 1],
            is3D ? (model.myosin.force[i].z * beta * D * dt +
                    model.myosin.velocity[i].z * dt +
                    std::sqrt(2 * D * dt) * noise[i * 6 + 2]) : 0.0
        };
        double disp_sq = delta_pos.x * delta_pos.x + delta_pos.y * delta_pos.y +
                         (is3D ? delta_pos.z * delta_pos.z : 0.0);
        double disp_limit = max_myosin_displacement;
        if (disp_limit > 0.0 && std::isfinite(disp_limit)) {
            double disp_mag = std::sqrt(disp_sq);
            double allowed = displacement_slack * disp_limit;
            if (disp_mag > allowed && disp_mag > 1e-12) {
                double scale = allowed / disp_mag;
                delta_pos.x *= scale;
                delta_pos.y *= scale;
                if (is3D) {
                    delta_pos.z *= scale;
                }
            }
        }
        vec rot_noise={noise[i * 6 + 3], noise[i * 6 + 4], is3D ? noise[i * 6 + 5] : 0.0};
        vec delta_u = std::sqrt(2 * D_rot * dt) * rot_noise + dt * model.myosin.torque[i] * D_rot * beta;
        if (!is3D) {
            delta_pos.z = 0.0;
            delta_u.z = 0.0;
        }
        double rot_limit = max_myosin_rotation;
        // if (rot_limit > 0.0 && std::isfinite(rot_limit)) {
        //     double rot_mag = delta_u.norm();
        //     double allowed = displacement_slack * rot_limit;
        // }
        double dx = delta_pos.x;
        double dy = delta_pos.y;
        double dz = delta_pos.z;
        model.myosin.displace(i, dx, dy, dz);
        vec wrapped_center_myo = model.myosin.center[i];
        wrapped_center_myo.pbc_wrap(model.box, periodic);
        model.myosin.center[i] = wrapped_center_myo;
        model.myosin.direction[i] += delta_u;
        if (!is3D) {
            model.myosin.direction[i].z = 0.0;
        }
        model.myosin.direction[i].normalize();
        model.myosin.update_endpoints(i);
    }
    // Update actin particles.
    for (int i = 0; i < model.actin.n; i++) {
        if (model.actin.cb_status[i] > 1){
            D = D_myosin_trans;
            D_rot = D_myosin_rot;
        }
        else{
            D = D_actin_trans;
            D_rot = D_actin_rot;
        }
        if (!is3D) {
            model.actin.force[i].z = 0;
            model.actin.velocity[i].z = 0;
            model.actin.torque[i].z = 0;
        }
        vec delta_pos{
            model.actin.force[i].x * beta * D * dt +
                model.actin.velocity[i].x * dt +
                std::sqrt(2 * D * dt) * noise[offset + i * 6],
            model.actin.force[i].y * beta * D * dt +
                model.actin.velocity[i].y * dt +
                std::sqrt(2 * D * dt) * noise[offset + i * 6 + 1],
            is3D ? (model.actin.force[i].z * beta * D * dt +
                    model.actin.velocity[i].z * dt +
                    std::sqrt(2 * D * dt) * noise[offset + i * 6 + 2]) : 0.0
        };
        double disp_sq = delta_pos.x * delta_pos.x + delta_pos.y * delta_pos.y +
                         (is3D ? delta_pos.z * delta_pos.z : 0.0);
        double disp_limit = max_actin_displacement;
        if (disp_limit > 0.0 && std::isfinite(disp_limit)) {
            double disp_mag = std::sqrt(disp_sq);
            double allowed = displacement_slack * disp_limit;
            if (disp_mag > allowed && disp_mag > 1e-12) {
                double scale = allowed / disp_mag;
                delta_pos.x *= scale;
                delta_pos.y *= scale;
                if (is3D) {
                    delta_pos.z *= scale;
                }
            }
        }
        vec rot_noise={noise[offset + i * 6 + 3], noise[offset + i * 6 + 4], is3D ? noise[offset + i * 6 + 5] : 0.0};
        vec delta_u = std::sqrt(2 * D_rot * dt) * rot_noise + dt * model.actin.torque[i] * D_rot * beta;
        if (!is3D) {
            delta_pos.z = 0.0;
            delta_u.z = 0.0;
        }
        double rot_limit = max_actin_rotation;
        if (rot_limit > 0.0 && std::isfinite(rot_limit)) {
            double rot_mag = delta_u.norm();
            double allowed = displacement_slack * rot_limit;
            if (rot_mag > allowed) {
                double noise_mag = (rot_noise).norm();
            }
        }
        double dx = delta_pos.x;
        double dy = delta_pos.y;
        double dz = delta_pos.z;
        model.actin.displace(i, dx, dy, dz);
        vec wrapped_center_act = model.actin.center[i];
        wrapped_center_act.pbc_wrap(model.box, periodic);
        model.actin.center[i] = wrapped_center_act;
        model.actin.direction[i] += delta_u;
        if (!is3D) {
            model.actin.direction[i].z = 0.0;
        }
        model.actin.direction[i].normalize();
        model.actin.update_endpoints(i);
    }
}
