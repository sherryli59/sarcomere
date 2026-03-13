#include "langevin.h"
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>

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
        skip_initial_save(false),
        loaded_frame_index(-1)
{
    if (resume) {
        int n_frames;
        int requested_frame = resume_frame;
        int frame_idx = model.load_state(n_frames, requested_frame);
        loaded_frame_index = frame_idx;
        // Frames are written at step multiples of save_every before integrating that step.
        // Resume should continue from that same step index.
        start_step = static_cast<int>(model.current_step);
        skip_initial_save = true;
        printf("Resuming from file %s at step %d (frame %d of %d, current_step=%zu)\n",
               model.filename.c_str(), start_step, frame_idx, n_frames, model.current_step);
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
        skip_initial_save = false;
        loaded_frame_index = -1;
        model.new_file();
    }
    // If resume mode active, interpret nsteps argument as total target
    // steps (i.e., inclusive of previously-run steps) when running.
    interpret_nsteps_as_total_when_resuming = resume;
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
    bool skip_first_save = skip_initial_save;
    int end_step;
    if (interpret_nsteps_as_total_when_resuming) {
        // Here `nsteps` is interpreted as the absolute target total number of steps.
        end_step = nsteps;
    } else {
        end_step = start_step + nsteps;
    }
    for (int step = start_step; step < end_step; ++step) {
        bool should_save = (step % save_every == 0);
        if (skip_first_save && step == start_step) {
            should_save = false;
        }
        if (should_save) {
            std::cout << "Step " << step << std::endl;
            // Optionally protect saving with the mutex:
            // std::lock_guard<std::mutex> lock(save_mutex);
            model.save_state();
            start = omp_get_wtime();
        }
        model.update_system();
        sample_step(dt, rng, fix_myosin);
        if (should_save) {
            end = omp_get_wtime();
            printf("Step %d took %f seconds\n", step, end - start);
            //model.debug_cb_stats();
        }
        if (skip_first_save && step == start_step) {
            skip_first_save = false;
        }
    }
    skip_initial_save = false;
    start_step = end_step;
    model.save_resume_snapshot();
}

void Langevin::volume_exclusion(int nsteps, gsl_rng* rng, int& fix_myosin) {
    double start, end;
    bool skip_first_save = skip_initial_save;
    int end_step;
    if (interpret_nsteps_as_total_when_resuming) {
        end_step = nsteps;
    } else {
        end_step = start_step + nsteps;
    }
    for (int step = start_step; step < end_step; ++step) {
        bool should_save = (step % save_every == 0);
        if (skip_first_save && step == start_step) {
            should_save = false;
        }
        if (should_save) {
            std::cout << "Step " << step << std::endl;
            // Optionally protect saving with the mutex:
            // std::lock_guard<std::mutex> lock(save_mutex);
            start = omp_get_wtime();
            model.save_state();
        }
        model.update_system_sterics_only();
        sample_step(dt, rng, fix_myosin);
        if (should_save) {
            end = omp_get_wtime();
            printf("Step %d took %f seconds\n", step, end - start);
        }
        if (skip_first_save && step == start_step) {
            skip_first_save = false;
        }
    }
    skip_initial_save = false;
    start_step = end_step;
    model.save_resume_snapshot();
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

    const double wall_eps = 1e-9;
    auto reflect_segment = [&](Filament& filament, int idx, vec& delta_pos, vec& delta_dir) {
        if (model.is_periodic[0] && model.is_periodic[1] && model.is_periodic[2]) {
            return;
        }
        int guard = 0;
        while (guard < 6) {
            vec center_current = filament.center[idx];
            vec direction_current = filament.direction[idx];
            vec center_prop = center_current + delta_pos;
            vec direction_prop = direction_current + delta_dir;
            if (!is3D) {
                direction_prop.z = 0.0;
            }
            double prop_norm = direction_prop.norm();
            if (prop_norm <= 1e-12) {
                direction_prop = direction_current;
            } else {
                direction_prop = direction_prop / prop_norm;
            }
            vec dir_norm = direction_prop;
            dir_norm.normalize();
            vec left_prop = center_prop - 0.5 * filament.length * dir_norm;
            vec right_prop = center_prop + 0.5 * filament.length * dir_norm;
            bool reflected = false;
            for (int axis = 0; axis < 3; ++axis) {
                if (model.is_periodic[axis]) {
                    continue;
                }
                if (axis >= static_cast<int>(model.box.size()) || model.box[axis] <= 0.0) {
                    continue;
                }
                double half = 0.5 * model.box[axis];
                double lower = -half;
                double upper = half;
                auto component = [&](const vec& v) -> double {
                    return axis == 0 ? v.x : (axis == 1 ? v.y : v.z);
                };
                double left_c = component(left_prop);
                double right_c = component(right_prop);
                if (left_c < lower - wall_eps || right_c < lower - wall_eps ||
                    left_c > upper + wall_eps || right_c > upper + wall_eps) {
                    if (axis == 0) {
                        delta_pos.x *= -1.0;
                        delta_dir.x *= -1.0;
                    } else if (axis == 1) {
                        delta_pos.y *= -1.0;
                        delta_dir.y *= -1.0;
                    } else {
                        delta_pos.z *= -1.0;
                        delta_dir.z *= -1.0;
                    }
                    reflected = true;
                }
            }
            if (!reflected) {
                break;
            }
            ++guard;
        }
    };

    auto clamp_center = [&](Filament& filament, int idx) {
        if (model.is_periodic[0] && model.is_periodic[1] && model.is_periodic[2]) {
            return;
        }
        vec center = filament.center[idx];
        vec dir = filament.direction[idx];
        for (int axis = 0; axis < 3; ++axis) {
            if (model.is_periodic[axis]) {
                continue;
            }
            if (axis >= static_cast<int>(model.box.size()) || model.box[axis] <= 0.0) {
                continue;
            }
            double half = 0.5 * model.box[axis];
            double lower = -half;
            double upper = half;
            double dir_component = axis == 0 ? dir.x : (axis == 1 ? dir.y : dir.z);
            double half_span = 0.5 * filament.length * std::abs(dir_component);
            double min_center = lower + half_span + wall_eps;
            double max_center = upper - half_span - wall_eps;
            double* center_component = (axis == 0) ? &center.x : (axis == 1 ? &center.y : &center.z);
            if (min_center > max_center) {
                *center_component = std::clamp(*center_component, lower + wall_eps, upper - wall_eps);
            } else {
                *center_component = std::clamp(*center_component, min_center, max_center);
            }
        }
        filament.center[idx] = center;
        filament.update_endpoints(idx);
    };

    int offset = model.myosin.n * 6;
    double D = D_myosin_trans;
    double D_rot = D_myosin_rot;
    // Update myosin particles.
    const double displacement_slack = 1.3;
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
        reflect_segment(model.myosin, i, delta_pos, delta_u);

        double dx = delta_pos.x;
        double dy = delta_pos.y;
        double dz = delta_pos.z;
        model.myosin.displace(i, dx, dy, dz);
        //print myosin force and displacement
        // printf("Myosin %d force: (%f, %f, %f), displacement: (%f, %f, %f)\n",
        //        i, model.myosin.force[i].x, model.myosin.force[i].y, model.myosin.force[i].z,
        //        dx, dy, dz);
        vec new_dir = static_cast<vec>(model.myosin.direction[i]) + delta_u;
        if (!is3D) {
            new_dir.z = 0.0;
            new_dir.normalize();
        } else {
            new_dir.normalize();
        }
        model.myosin.direction[i] = new_dir;
        model.myosin.update_endpoints(i);
        clamp_center(model.myosin, i);
    }
    // Update actin particles.
    for (int i = 0; i < model.actin.n; i++) {
        // if (model.actin.cb_status[i] > 1){
        //     D = D_myosin_trans;
        //     D_rot = D_myosin_rot;
        // }
        // else{
        //     D = D_actin_trans;
        //     D_rot = D_actin_rot;
        // }
        D = D_actin_trans;
        D_rot = D_actin_rot;
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
        reflect_segment(model.actin, i, delta_pos, delta_u);

        double dx = delta_pos.x;
        double dy = delta_pos.y;
        double dz = delta_pos.z;
        model.actin.displace(i, dx, dy, dz);
        vec new_dir_act = static_cast<vec>(model.actin.direction[i]) + delta_u;
        if (!is3D) {
            new_dir_act.z = 0.0;
            new_dir_act.normalize();
        } else {
            new_dir_act.normalize();
        }
        model.actin.direction[i] = new_dir_act;
        model.actin.update_endpoints(i);
        clamp_center(model.actin, i);
    }
}
