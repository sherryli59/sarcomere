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
        if (requested_frame >= 0 && frame_idx != requested_frame) {
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

void Langevin::set_rod_mobility(
    const rod_mobility::Diffusion& actin_mobility,
    const rod_mobility::Diffusion& myosin_mobility) {
    if (actin_mobility.enabled != myosin_mobility.enabled) {
        throw std::invalid_argument(
            "actin and myosin tensor mobility must be enabled together");
    }
    actin_rod_mobility = actin_mobility;
    myosin_rod_mobility = myosin_mobility;
    physical_rod_mobility = actin_mobility.enabled;
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
            model.save_resume_snapshot();
            end = omp_get_wtime();
            std::cout << "Step " << step << " took " << (end - start) << " seconds" << std::endl;
            //model.debug_cb_stats();
        }
        if (skip_first_save && step == start_step) {
            skip_first_save = false;
        }
    }
    if (end_step > start_step && (end_step % save_every == 0)) {
        std::cout << "Step " << end_step << std::endl;
        model.save_state();
    }
    skip_initial_save = false;
    start_step = end_step;
    model.save_resume_snapshot();
}

void Langevin::volume_exclusion(int nsteps, gsl_rng* rng, int& fix_myosin,
                                bool save_frames) {
    double start, end;
    bool skip_first_save = skip_initial_save;
    int end_step;
    if (interpret_nsteps_as_total_when_resuming) {
        end_step = nsteps;
    } else {
        end_step = start_step + nsteps;
    }
    for (int step = start_step; step < end_step; ++step) {
        bool should_save = save_frames && (step % save_every == 0);
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
            model.save_resume_snapshot();
            end = omp_get_wtime();
            std::cout << "Step " << step << " took " << (end - start) << " seconds" << std::endl;
        }
        if (skip_first_save && step == start_step) {
            skip_first_save = false;
        }
    }
    if (save_frames && end_step > start_step && (end_step % save_every == 0)) {
        std::cout << "Step " << end_step << std::endl;
        model.save_state();
    }
    skip_initial_save = false;
    start_step = end_step;
    if (save_frames) model.save_resume_snapshot();
}

void Langevin::resolve_myosin_hard_contacts(int fix_myosin, int max_sweeps) {
    if (!model.sterics_enabled || model.myosin.n < 2 || max_sweeps <= 0) {
        return;
    }

    // Refresh the cached positions.  The Verlet skin preserves the pair set
    // between rebuilds; needs_rebuild() triggers before a filament can move
    // across that skin without entering the list.
    model.neighbor_list.set_species_positions(
        model.actin.center_x, model.actin.center_y, model.actin.center_z,
        model.myosin.center_x, model.myosin.center_y, model.myosin.center_z);
    if (model.neighbor_list.needs_rebuild()) {
        model.neighbor_list.rebuild_neighbor_list();
    }

    auto& contact_events =
        model.myosin["cumulative_hard_contact_events"];
    auto& contact_correction =
        model.myosin["cumulative_hard_contact_correction"];
    const double target = 2.0 * model.myosin.radius;

    // Sequential Gauss-Seidel sweeps are deterministic and conserve the pair
    // midpoint whenever both rods are mobile.  Repeated sweeps are needed
    // because resolving one contact in a dense cluster can close another.
    for (int sweep = 0; sweep < max_sweeps; ++sweep) {
        int corrected = 0;
        for (int i = 0; i < model.myosin.n; ++i) {
            for (const auto& [neighbor_index, neighbor_type] :
                 model.neighbor_list.get_neighbors(i + model.actin.n)) {
                if (neighbor_type != ParticleType::Myosin) continue;
                const int j = neighbor_index - model.actin.n;
                if (i >= j) continue;
                const vec old_direction_i = model.myosin.direction[i];
                const vec old_direction_j = model.myosin.direction[j];
                const double translation_weight = physical_rod_mobility
                    ? myosin_rod_mobility.perpendicular
                    : std::max(D_myosin_trans, 1e-12);
                const double rotation_weight = physical_rod_mobility
                    ? myosin_rod_mobility.tumble
                    : std::max(D_myosin_rot, 0.0);
                const auto projection = project_myosin_hard_contact(
                    model.myosin, i, j, model.box, fix_myosin, target,
                    translation_weight, rotation_weight);
                if (!projection.applied) continue;
                if (i >= fix_myosin &&
                    (static_cast<vec>(model.myosin.direction[i]) -
                     old_direction_i).norm() > 0.0) {
                    model.transport_myosin_material_frame(
                        i, old_direction_i, model.myosin.direction[i], 0.0);
                }
                if (j >= fix_myosin &&
                    (static_cast<vec>(model.myosin.direction[j]) -
                     old_direction_j).norm() > 0.0) {
                    model.transport_myosin_material_frame(
                        j, old_direction_j, model.myosin.direction[j], 0.0);
                }
                ++corrected;
                contact_events[i] += 1.0;
                contact_events[j] += 1.0;
                contact_correction[i] += 0.5 * projection.correction;
                contact_correction[j] += 0.5 * projection.correction;
            }
        }
        if (corrected == 0) break;
    }
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
        noise[i] = thermal_noise_scale * gsl_ran_gaussian(rng, 1.0);
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
        vec translation_noise{noise[i * 6], noise[i * 6 + 1],
                              is3D ? noise[i * 6 + 2] : 0.0};
        vec delta_pos;
        vec myosin_axis = model.myosin.direction[i];
        myosin_axis.normalize();
        if (physical_rod_mobility) {
            delta_pos = beta * dt * rod_mobility::apply_diffusion_tensor(
                model.myosin.force[i], myosin_axis, myosin_rod_mobility);
            delta_pos += dt * model.myosin.velocity[i];
            delta_pos += std::sqrt(2.0 * dt) *
                rod_mobility::apply_sqrt_diffusion_tensor(
                    translation_noise, myosin_axis, myosin_rod_mobility);
        } else {
            delta_pos = beta * D * dt * model.myosin.force[i] +
                        dt * model.myosin.velocity[i] +
                        std::sqrt(2.0 * D * dt) * translation_noise;
        }
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
        vec delta_u;
        double axial_roll = 0.0;
        if (physical_rod_mobility) {
            if (is3D) {
                const vec tangent_drive = rod_mobility::project_perpendicular(
                    model.myosin.torque[i], myosin_axis);
                const vec tangent_noise = rod_mobility::project_perpendicular(
                    rot_noise, myosin_axis);
                delta_u = beta * myosin_rod_mobility.tumble * dt * tangent_drive +
                    std::sqrt(2.0 * myosin_rod_mobility.tumble * dt) * tangent_noise;
                axial_roll = beta * myosin_rod_mobility.roll * dt *
                                 model.myosin.torque[i].dot(myosin_axis) +
                             std::sqrt(2.0 * myosin_rod_mobility.roll * dt) *
                                 rot_noise.dot(myosin_axis);
            } else {
                const vec tangent{-myosin_axis.y, myosin_axis.x, 0.0};
                const double scalar = beta * myosin_rod_mobility.tumble * dt *
                                          model.myosin.torque[i].dot(tangent) +
                                      std::sqrt(2.0 * myosin_rod_mobility.tumble * dt) *
                                          rot_noise.dot(tangent);
                delta_u = scalar * tangent;
            }
        } else {
            delta_u = std::sqrt(2 * D_rot * dt) * rot_noise +
                      dt * model.myosin.torque[i] * D_rot * beta;
            axial_roll = delta_u.dot(myosin_axis);
        }
        if (!is3D) {
            delta_pos.z = 0.0;
            delta_u.z = 0.0;
        }
        double rot_limit = max_myosin_rotation;
        const double rotation_scale = rod_mobility::angular_increment_scale(
            delta_u, rot_limit, displacement_slack);
        if (rotation_scale < 1.0) {
            delta_u = rotation_scale * delta_u;
            // A large tangent step and material-frame roll belong to the same
            // realized Brownian/mechanical update. Keep their time increment
            // coherent whenever the director trust radius is activated.
            axial_roll *= rotation_scale;
            model.myosin["cumulative_rotation_cap_events"][i] += 1.0;
        }
        reflect_segment(model.myosin, i, delta_pos, delta_u);

        double dx = delta_pos.x;
        double dy = delta_pos.y;
        double dz = delta_pos.z;
        model.myosin.displace(i, dx, dy, dz);
        //print myosin force and displacement
        // printf("Myosin %d force: (%f, %f, %f), displacement: (%f, %f, %f)\n",
        //        i, model.myosin.force[i].x, model.myosin.force[i].y, model.myosin.force[i].z,
        //        dx, dy, dz);
        const vec old_dir = model.myosin.direction[i];
        vec new_dir = physical_rod_mobility
            ? rod_mobility::sphere_exp(myosin_axis, delta_u)
            : old_dir + delta_u;
        if (!is3D) {
            new_dir.z = 0.0;
            new_dir.normalize();
        } else {
            new_dir.normalize();
        }
        model.myosin.direction[i] = new_dir;
        model.transport_myosin_material_frame(
            i, old_dir, new_dir, axial_roll);
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
        vec translation_noise{noise[offset + i * 6], noise[offset + i * 6 + 1],
                              is3D ? noise[offset + i * 6 + 2] : 0.0};
        vec delta_pos;
        vec actin_axis = model.actin.direction[i];
        actin_axis.normalize();
        if (physical_rod_mobility) {
            delta_pos = beta * dt * rod_mobility::apply_diffusion_tensor(
                model.actin.force[i], actin_axis, actin_rod_mobility);
            delta_pos += dt * model.actin.velocity[i];
            delta_pos += std::sqrt(2.0 * dt) *
                rod_mobility::apply_sqrt_diffusion_tensor(
                    translation_noise, actin_axis, actin_rod_mobility);
        } else {
            delta_pos = beta * D * dt * model.actin.force[i] +
                        dt * model.actin.velocity[i] +
                        std::sqrt(2.0 * D * dt) * translation_noise;
        }
        if (constrain_actin_axial) {
            delta_pos -= actin_axis * delta_pos.dot(actin_axis);
        }
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
        vec delta_u;
        if (physical_rod_mobility) {
            if (is3D) {
                delta_u = beta * actin_rod_mobility.tumble * dt *
                              rod_mobility::project_perpendicular(
                                  model.actin.torque[i], actin_axis) +
                          std::sqrt(2.0 * actin_rod_mobility.tumble * dt) *
                              rod_mobility::project_perpendicular(rot_noise, actin_axis);
            } else {
                const vec tangent{-actin_axis.y, actin_axis.x, 0.0};
                const double scalar = beta * actin_rod_mobility.tumble * dt *
                                          model.actin.torque[i].dot(tangent) +
                                      std::sqrt(2.0 * actin_rod_mobility.tumble * dt) *
                                          rot_noise.dot(tangent);
                delta_u = scalar * tangent;
            }
        } else {
            delta_u = std::sqrt(2 * D_rot * dt) * rot_noise +
                      dt * model.actin.torque[i] * D_rot * beta;
        }
        if (!is3D) {
            delta_pos.z = 0.0;
            delta_u.z = 0.0;
        }
        double rot_limit = max_actin_rotation;
        const double rotation_scale = rod_mobility::angular_increment_scale(
            delta_u, rot_limit, displacement_slack);
        if (rotation_scale < 1.0) {
            delta_u = rotation_scale * delta_u;
            model.actin["cumulative_rotation_cap_events"][i] += 1.0;
        }
        reflect_segment(model.actin, i, delta_pos, delta_u);

        double dx = delta_pos.x;
        double dy = delta_pos.y;
        double dz = delta_pos.z;
        model.actin.displace(i, dx, dy, dz);
        vec new_dir_act = physical_rod_mobility
            ? rod_mobility::sphere_exp(actin_axis, delta_u)
            : static_cast<vec>(model.actin.direction[i]) + delta_u;
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
    resolve_myosin_hard_contacts(fix_myosin);
}
