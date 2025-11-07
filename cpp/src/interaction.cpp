#include "interaction.h"
#include "geometry.h"
#include <cstdio>
#include <cstdlib>  // For exit()
#include <cmath>
#include <limits>
#include <algorithm>

//---------------------------------------------------------------------
// Non-Template Function Definitions
//---------------------------------------------------------------------




real am_energy1(const ArrayXreal& center1, const double& length1, const ArrayXreal& dir1,
    const ArrayXreal& center2, const double& length2, const ArrayXreal& dir2,
    const std::vector<double>& box, const std::array<bool,3>& periodic, const double k_am, const double kappa_am,
    const double cutoff, const double optimal)
{
    // Compute endpoints for both filaments
    ArrayXreal a = center1 - 0.5 * length1 * dir1;
    ArrayXreal b = center1 + 0.5 * length1 * dir1;
    ArrayXreal c = center2 - 0.5 * length2 * dir2;
    ArrayXreal d = center2 + 0.5 * length2 * dir2;

    // Convert to raw arrays for segment-segment distance
    real a_arr[3] = {a[0], a[1], a[2]};
    real b_arr[3] = {b[0], b[1], b[2]};
    real c_arr[3] = {c[0], c[1], c[2]};
    real d_arr[3] = {d[0], d[1], d[2]};
    real dist = segment_segment_distance(a_arr, b_arr, c_arr, d_arr, box, periodic);
    // Angle and strength
    real dot_val = dir1[0]*dir2[0] + dir1[1]*dir2[1] + dir1[2]*dir2[2];
    real strength = abs(dot_val);
    dot_val = std::max(real(-1), std::min(real(1), dot_val));  // clamp
    real angle_energy = 0.5 * kappa_am * (1.0 - dot_val * dot_val); 
    if (dist > cutoff) {
        printf("something's wrong. dist: %f\n", dist.val());
        return angle_energy;
    }
    else {
        real offset = dist - optimal;
        real rel = (optimal != 0.0) ? offset / optimal : offset;
        return 0.5 * k_am * strength * rel * rel + angle_energy;
    }
}

real am_energy(const ArrayXreal& dir1, const ArrayXreal& dir2, const double kappa_am)
{
    real dot_val = dir1[0]*dir2[0] + dir1[1]*dir2[1] + dir1[2]*dir2[2];
    dot_val = std::max(real(-1), std::min(real(1), dot_val));  
    real strength = abs(dot_val);
    real angle_energy = 0.5 * kappa_am * (1.0 - dot_val * dot_val); 
    return angle_energy;
}


real aa_energy(const ArrayXreal& center1, const double& length1, 
    const ArrayXreal& dir1,
    const ArrayXreal& center2, const double& length2, 
    const ArrayXreal& dir2,
    const std::vector<double>& box, const std::array<bool,3>& periodic, 
    const double k_aa, const double kappa_aa,
    const double cutoff, const double optimal)
{
    // Compute endpoints for filament 1 in 3D.
    ArrayXreal a = center1 - 0.5 * length1 * dir1;
    ArrayXreal b = center1 + 0.5 * length1 * dir1;

    // Compute endpoints for filament 2 in 3D.
    ArrayXreal c = center2 - 0.5 * length2 * dir2;
    ArrayXreal d = center2 + 0.5 * length2 * dir2;

    // Convert to raw arrays for segment_segment_distance
    real a_arr[3] = {a[0], a[1], a[2]};
    real b_arr[3] = {b[0], b[1], b[2]};
    real c_arr[3] = {c[0], c[1], c[2]};
    real d_arr[3] = {d[0], d[1], d[2]};

    real raw_dist = segment_segment_distance(a_arr, b_arr, c_arr, d_arr, box, periodic);
    real dot_val = dir1[0]*dir2[0] + dir1[1]*dir2[1] + dir1[2]*dir2[2];
    dot_val = std::max(real(-1), std::min(real(1), dot_val));  // clamp for safety
    real angle_energy = 0.5 * kappa_aa * (1.0 - dot_val * dot_val);
    if (cutoff > 0.0 && raw_dist > cutoff) {
        return angle_energy;
    }
    real dist = raw_dist - optimal;
    real rel = (optimal != 0.0) ? dist / optimal : dist;

    // ===== Axial-gap soft wall between barbed ends a (filament1) and c (filament2) =====
    // Compute minimum-image displacement c - a
    real dx = c_arr[0] - a_arr[0];
    real dy = c_arr[1] - a_arr[1];
    real dz = c_arr[2] - a_arr[2];
    real box_arr[3] = {box.size() > 0 ? box[0] : 0.0,
                       box.size() > 1 ? box[1] : 0.0,
                       box.size() > 2 ? box[2] : 0.0};
    if (periodic[0] && box_arr[0] != real(0.0)) {
        dx -= box_arr[0] * smooth_round(dx / box_arr[0]);
    }
    if (periodic[1] && box_arr[1] != real(0.0)) {
        dy -= box_arr[1] * smooth_round(dy / box_arr[1]);
    }
    if (periodic[2] && box_arr[2] != real(0.0)) {
        dz -= box_arr[2] * smooth_round(dz / box_arr[2]);
    }
    // Axial offset along filament 1's axis
    real ds = dx*dir1[0] + dy*dir1[1] + dz*dir1[2];
    // One-sided quadratic penalty beyond s_tol
    real U_ax = 0.0;
    if (ds < 0.0) {
        U_ax = 0.5 * k_aa * ds * ds;
    }    
    return 0.5 * (k_aa * rel * rel) + angle_energy + U_ax;
}


std::vector<double> compute_aa_force_and_energy(Filament& actin,
                                                int& actin1_index, int& actin2_index,
                                                const std::vector<double>& box,
                                                const double k_aa, const double kappa_aa,
                                                const double cutoff, const double optimal)
{
    const double EPS = 1e-12;
    std::vector<double> forces(9, 0.0);

    vec left1 = actin.left_end[actin1_index];
    vec right1 = actin.right_end[actin1_index];
    vec left2 = actin.left_end[actin2_index];
    vec right2 = actin.right_end[actin2_index];

    vec dir1 = actin.direction[actin1_index];
    vec dir2 = actin.direction[actin2_index];
    dir1.normalize();
    dir2.normalize();

    auto geom = geometry::segment_segment_distance_w_normal(
        left1, right1, left2, right2, box, actin.periodic_axes);

    double distance = geom.first;
    bool within_cutoff = (cutoff <= 0.0) || (distance <= cutoff);
    if (!within_cutoff) {
        return forces;}

    vec shortest = {0.0, 0.0, 0.0};
    auto normal_it = geom.second.find("normal");
    if (normal_it != geom.second.end()) {
        shortest = normal_it->second;
    }

    double dist_norm = shortest.norm();
    vec unit = {0.0, 0.0, 0.0};
    if (dist_norm <= EPS) {
        vec cross_dir = dir1.cross(dir2);
        double cross_norm = cross_dir.norm();
        if (cross_norm > EPS) {
            unit = cross_dir / cross_norm;
        } else {
            unit = utils::pbc_diff_masked(actin.center[actin1_index],
                                              actin.center[actin2_index],
                                              box,
                                              actin.periodic_axes);
            unit = unit / unit.norm();
        }
    }
    else {
        unit = shortest / dist_norm;
    }
    double delta = dist_norm - optimal;
    double relative_extension = (optimal != 0.0) ? delta / optimal : delta;
    vec force_vec = -k_aa * relative_extension * unit;
    double dot = std::clamp(dir1.dot(dir2), -1.0, 1.0);
    vec cross12 = dir1.cross(dir2);
    vec torque1 = kappa_aa * dot * cross12;
    vec torque2 = -torque1;
    forces[0] = force_vec.x;
    forces[1] = force_vec.y;
    forces[2] = force_vec.z;
    forces[3] = torque1.x;
    forces[4] = torque1.y;
    forces[5] = torque1.z;
    forces[6] = torque2.x;
    forces[7] = torque2.y;
    forces[8] = torque2.z;
    return forces;
}

std::vector<double> compute_aa_force_and_energy_autodiff(Filament& actin,
                                                         int& actin1_index, int& actin2_index,
                                                         const std::vector<double>& box,
                                                         const double k_aa, const double kappa_aa,
                                                         const double cutoff, const double optimal)
{
    // Construct 3D centers.
    ArrayXreal center1(3);
    center1 << actin.center[actin1_index].x, actin.center[actin1_index].y, actin.center[actin1_index].z;
    ArrayXreal dir1(3);
    dir1 << actin.direction[actin1_index].x, actin.direction[actin1_index].y, actin.direction[actin1_index].z;

    ArrayXreal center2(3);
    center2 << actin.center[actin2_index].x, actin.center[actin2_index].y, actin.center[actin2_index].z;
    ArrayXreal dir2(3);
    dir2 << actin.direction[actin2_index].x, actin.direction[actin2_index].y, actin.direction[actin2_index].z;

    real u;
    std::vector<double> forces;
    // Compute the gradient of aa_energy with respect to center1, theta1, phi1, theta2, and phi2.
    VectorXd forces_2 = -gradient(aa_energy, wrt(center1, dir1, dir2),
                                  at(center1, actin.length, dir1,
                                     center2, actin.length, dir2, box, actin.periodic_axes,
                                     k_aa, kappa_aa, cutoff, optimal), u);
    forces.resize(forces_2.size());
    VectorXd::Map(&forces[0], forces_2.size()) = forces_2;
    return forces;
}

std::vector<double> compute_am_force_and_energy(Filament& actin, Myosin& myosin,
                                                int& actin_index, int& myosin_index,
                                                const std::vector<double>& box,
                                                const double k_am, const double kappa_am,
                                                const double cutoff, const double optimal)
{
    const double EPS = 1e-12;
    std::vector<double> forces(9, 0.0);

    vec act_left = actin.left_end[actin_index];
    vec act_right = actin.right_end[actin_index];
    vec myo_left = myosin.left_end[myosin_index];
    vec myo_right = myosin.right_end[myosin_index];

    vec act_dir = actin.direction[actin_index];
    vec myo_dir = myosin.direction[myosin_index];
    act_dir.normalize();
    myo_dir.normalize();
    if (k_am > EPS){
        auto geom = geometry::segment_segment_distance_w_normal(
            act_left, act_right, myo_left, myo_right, box, actin.periodic_axes);

        double distance = geom.first;
        bool within_cutoff = (cutoff <= 0.0) || (distance <= cutoff);
        if (!within_cutoff) {
            return forces;
        }
        vec shortest = {0.0, 0.0, 0.0};
        auto normal_it = geom.second.find("normal");
        if (normal_it != geom.second.end()) {
            shortest = normal_it->second;
        }
        double dist_norm = shortest.norm();
        vec unit = {0.0, 0.0, 0.0};
        if (dist_norm <= EPS) {
            vec cross_dir = act_dir.cross(myo_dir);
            double cross_norm = cross_dir.norm();
            if (cross_norm > EPS) {
                unit = cross_dir / cross_norm;
            } else {
                unit = utils::pbc_diff_masked(actin.center[actin_index],
                                                myosin.center[myosin_index],
                                                box,
                                                actin.periodic_axes);
                unit = unit / unit.norm();
            }
        }
        else {
            unit = shortest / dist_norm;
        }
        double delta = distance - optimal;
        double relative_extension = (optimal != 0.0) ? delta / optimal : delta;
        vec force_vec = -k_am * relative_extension * unit;    
        forces[0] = force_vec.x;
        forces[1] = force_vec.y;
        forces[2] = force_vec.z;
        // Apply end-stop forces near myosin termini, preserving previous behaviour.
        vec endstop_force = {0.0, 0.0, 0.0};
        double am_dot = act_dir.dot(myo_dir);
        vec myosin_end = (am_dot > 0.0) ? myosin.left_end[myosin_index]
                                        : myosin.right_end[myosin_index];
        vec actin_tip = actin.right_end[actin_index];
        vec tip_disp = actin_tip - myosin_end;
        tip_disp.pbc_wrap(box, actin.periodic_axes);
        double s = std::fabs(tip_disp.dot(myo_dir));
        double Lm = 0.4 * myosin.length;
        if (s >= Lm) {
            double dist_to_end = s - Lm;
            double mag_mid = k_am * dist_to_end;
            endstop_force = -mag_mid * act_dir;
            forces[0] += endstop_force.x;
            forces[1] += endstop_force.y;
            forces[2] += endstop_force.z;
        }
    }
    double dot = std::clamp(act_dir.dot(myo_dir), -1.0, 1.0);
    vec cross_am = act_dir.cross(myo_dir);
    vec torque_act = kappa_am * dot * cross_am;
    vec torque_myo = -torque_act;
    forces[3] = torque_act.x;
    forces[4] = torque_act.y;
    forces[5] = torque_act.z;
    forces[6] = torque_myo.x;
    forces[7] = torque_myo.y;
    forces[8] = torque_myo.z;
    return forces;
}

std::vector<double> compute_am_force_and_energy_autodiff(Filament& actin, Myosin& myosin,
                                                         int& actin_index, int& myosin_index,
                                                         const std::vector<double>& box,
                                                         const double k_am, const double kappa_am,
                                                         const double cutoff, const double optimal)
{
    // Define 3D center positions
    ArrayXreal center1(3);
    center1[0] = actin.center[actin_index].x;
    center1[1] = actin.center[actin_index].y;
    center1[2] = actin.center[actin_index].z;
    ArrayXreal dir1(3);
    dir1[0] = actin.direction[actin_index].x;
    dir1[1] = actin.direction[actin_index].y;
    dir1[2] = actin.direction[actin_index].z;
    ArrayXreal center2(3);
    center2[0] = myosin.center[myosin_index].x;
    center2[1] = myosin.center[myosin_index].y;
    center2[2] = myosin.center[myosin_index].z;
    ArrayXreal dir2(3);
    dir2[0] = myosin.direction[myosin_index].x;
    dir2[1] = myosin.direction[myosin_index].y;
    dir2[2] = myosin.direction[myosin_index].z;
    real u;
    std::vector<double> forces;
    VectorXd forces_3; // 3D version of force vector
    if (k_am > 1e-6) {
        forces_3 = -gradient(am_energy1, wrt(center1, dir1, dir2),
                             at(center1, actin.length, dir1,
                                center2, myosin.length, dir2, box, actin.periodic_axes,
                                k_am, kappa_am, cutoff, optimal), u);
        forces.resize(forces_3.size());
        VectorXd::Map(&forces[0], forces_3.size()) = forces_3;
        // Now add in the myosin centerstop forces
        vec actin_dir = actin.direction[actin_index];
        double am_dot = actin.direction[actin_index].dot(myosin.direction[myosin_index]);
        vec myosin_end;
        if (am_dot > 0) {
            myosin_end = myosin.left_end[myosin_index];
        }
        else {
            myosin_end = myosin.right_end[myosin_index];
        }
        vec actin_tip = actin.right_end[actin_index];
        vec tip_disp = actin_tip - myosin_end;
        tip_disp.pbc_wrap(box, actin.periodic_axes);
        double s = std::fabs(tip_disp.dot(myosin.direction[myosin_index]));  // +s toward "right" tip, −s toward "left"
        double Lm = 0.4 * myosin.length;                     // half-length of myosin
        vec endstop_force = {0.0, 0.0, 0.0};
        if (s >= Lm) {
            double dist_to_end = s - Lm;
            double mag_mid = k_am * dist_to_end;         // increase toward center
            // choose outward sign; for s==0 pick +1 by convention
            endstop_force = - mag_mid * actin_dir;         // pushes toward the nearer end
        }
        forces[0] += endstop_force.x;
        forces[1] += endstop_force.y;
        forces[2] += endstop_force.z;
    }
    else {
        forces_3 = -gradient(am_energy, wrt(dir1, dir2),
                             at(dir1, dir2, kappa_am), u);
        // Prepend three zeros for the force vector (x, y, z)
        forces.resize(forces_3.size() + 3);
        forces[0] = 0;
        forces[1] = 0;
        forces[2] = 0;
        VectorXd::Map(&forces[3], forces_3.size()) = forces_3;
    }
    return forces;
}

RepulsionResult compute_myosin_repulsion(const Filament& actin,
                                         const Myosin& myosin,
                                         int i,
                                         int j,
                                         const std::vector<double>& box,
                                         int fix_myosin,
                                         const utils::MoleculeConnection& actinIndicesPerMyosin,
                                         double stiffness,
                                         double max_force_cap)
{
    RepulsionResult result{};
    const double cutoff = 2.0 * myosin.radius;
    const double EPS = 1e-9;
    const double max_force_limit = (max_force_cap > 0.0 && std::isfinite(max_force_cap))
                                       ? max_force_cap
                                       : std::numeric_limits<double>::infinity();

    vec center_displacement = myosin.center[i] - myosin.center[j];
    center_displacement.pbc_wrap(box, myosin.periodic_axes);
    const double center_distance = center_displacement.norm();
    if (center_distance > cutoff + myosin.length) {
        return result;
    }

    auto geom = geometry::segment_segment_distance_w_normal(
        myosin.left_end[i], myosin.right_end[i],
        myosin.left_end[j], myosin.right_end[j],
        box, myosin.periodic_axes);

    const double distance = geom.first;
    if (distance >= cutoff) {
        return result;
    }

    auto normal_it = geom.second.find("normal");
    if (normal_it == geom.second.end()) {
        return result;
    }

    vec normal_vector = normal_it->second;
    double norm = normal_vector.norm();
    if (norm <= EPS) {
        const double denom = std::max(center_distance, EPS);
        if (denom <= EPS) {
            return result;
        }
        normal_vector = center_displacement / denom;
    } else {
        normal_vector = normal_vector / norm;
    }

    double overlap = cutoff - distance;
    if (overlap <= 0.0) {
        return result;
    }

    double effective_stiffness = (std::isfinite(stiffness) && stiffness > 0.0) ? stiffness : 0.0;
    if (effective_stiffness <= 0.0) {
        return result;
    }

    double force_mag = effective_stiffness * overlap;
    if (!std::isfinite(force_mag) || force_mag < 0.0) {
        force_mag = 0.0;
    }
    if (force_mag > max_force_limit) {
        force_mag = max_force_limit;
    }

    vec force_vec = force_mag * normal_vector;

    result.applied = true;

    if (i < fix_myosin) {
        result.force_on_second -= force_vec * 2.0;
        return result;
    }
    if (j < fix_myosin) {
        result.force_on_first += force_vec * 2.0;
        return result;
    }

    int status_i = 0;
    const auto& connections_i = actinIndicesPerMyosin.getConnections(i);
    for (int act_idx : connections_i) {
        status_i = std::max(status_i, actin.cb_status[act_idx]);
    }

    int status_j = 0;
    const auto& connections_j = actinIndicesPerMyosin.getConnections(j);
    for (int act_idx : connections_j) {
        status_j = std::max(status_j, actin.cb_status[act_idx]);
    }

    if (status_i > 0 && status_j == 0) {
        result.force_on_second -= force_vec * 2.0;
    } else if (status_i == 0 && status_j > 0) {
        result.force_on_first += force_vec * 2.0;
    } else {
        result.force_on_first += force_vec;
        result.force_on_second -= force_vec;
    }
    return result;
}

RepulsionResult compute_actin_repulsion(const Filament& actin,
                                        int i,
                                        int j,
                                        const std::vector<double>& box,
                                        double crosslinker_length,
                                        double stiffness,
                                        double max_force_cap)
{
    RepulsionResult result{};
    const double max_force_limit = (max_force_cap > 0.0 && std::isfinite(max_force_cap))
                                       ? max_force_cap
                                       : std::numeric_limits<double>::infinity();
    const double EPS = 1e-9;

    vec center_displacement = actin.center[i] - actin.center[j];
    center_displacement.pbc_wrap(box, actin.periodic_axes);
    const double center_distance = center_displacement.norm();
    if (center_distance > crosslinker_length + actin.length) {
        return result;
    }

    auto geom = geometry::segment_segment_distance_w_normal(
        actin.left_end[i], actin.right_end[i],
        actin.left_end[j], actin.right_end[j],
        box, actin.periodic_axes);

    const double distance = geom.first;
    if (distance >= crosslinker_length) {
        return result;
    }

    auto normal_it = geom.second.find("normal");
    if (normal_it == geom.second.end()) {
        return result;
    }

    vec normal_vector = normal_it->second;
    double norm = normal_vector.norm();
    if (norm <= EPS) {
        normal_vector = actin.direction[i].cross(actin.direction[j]);
        double fallback_norm = normal_vector.norm();
        if (fallback_norm <= EPS) {
            return result;
        }
        normal_vector = normal_vector / fallback_norm;
    } else {
        normal_vector = normal_vector / norm;
    }

    double overlap = crosslinker_length - distance;
    if (overlap <= 0.0) {
        return result;
    }

    double force_mag = stiffness * overlap;
    if (!std::isfinite(force_mag) || force_mag < 0.0) {
        force_mag = 0.0;
    }
    // if (force_mag > max_force_limit) {
    //     force_mag = max_force_limit;
    // }
    vec force_vec = force_mag * normal_vector;

    const int status_i = actin.cb_status[i];
    const int status_j = actin.cb_status[j];
    const bool both_status_two = (status_i == 2 && status_j == 2);
    const bool apply_first = (status_i < 2) || both_status_two;
    const bool apply_second = (status_j < 2) || both_status_two;

    if (apply_first) {
        result.force_on_first += force_vec;
    }
    if (apply_second) {
        result.force_on_second -= force_vec;
    }
    result.applied = apply_first || apply_second;
    return result;
}

vec compute_actin_myosin_repulsion(const Filament& actin,
                                   const Myosin& myosin,
                                   int act_idx,
                                   int myo_idx,
                                   const std::vector<double>& box,
                                   double radius,
                                   double stiffness,
                                   double max_force_cap)
{
    // const double EPS = 1e-9;
    // const double max_force_limit = (max_force_cap > 0.0 && std::isfinite(max_force_cap))
    //                                    ? max_force_cap
    //                                    : std::numeric_limits<double>::infinity();

    auto geom = geometry::segment_segment_distance_w_normal(
        actin.left_end[act_idx], actin.right_end[act_idx],
        myosin.left_end[myo_idx], myosin.right_end[myo_idx],
        box, actin.periodic_axes);

    const double distance = geom.first;
    if (distance >= radius) {
        return {0.0, 0.0, 0.0};
    }

    // auto normal_it = geom.second.find("normal");
    // vec normal_vector = normal_it->second;
    // double norm = normal_vector.norm();
    // vec dir;
    // if (norm <= EPS) {
    //     vec center_displacement = actin.center[act_idx] - myosin.center[myo_idx];
    //     center_displacement.pbc_wrap(box, actin.periodic_axes);
    //     double center_norm = center_displacement.norm();
    //     dir = center_displacement / center_norm;
    // } else {
    //     dir = normal_vector / norm;
    // }

    // double overlap = radius - distance;

    // double magnitude = stiffness * overlap;
    
    // if (magnitude > max_force_limit) {
    //     magnitude = max_force_limit;
    // }
    // vec repulsive_force = magnitude * dir;


    vec u = myosin.direction[myo_idx];
    // Signed axial coordinate of the actin *left* tip relative to myosin center
    vec tip_disp = actin.left_end[act_idx] - myosin.center[myo_idx];
    tip_disp.pbc_wrap(box, actin.periodic_axes);
    double s  = tip_disp.dot(u);                         // +s toward "right" tip, −s toward "left"
    double Lm = 0.5 * myosin.length;                     // half-length of myosin
    vec endstop_force = {0.0, 0.0, 0.0};
    if (std::fabs(s) <= Lm) {
        double dist_to_end = Lm - std::fabs(s);       // ∈ [0, Lm]
        double mag_mid = stiffness * dist_to_end;         // increase toward center
        vec direction = (s >= 0) ? u : -u; // +s toward "right" tip, −s toward "left"
        endstop_force = - mag_mid * direction; // pushes outward along myosin direction
        // printf("actin-myosin (%d, %d) repulsion force magnitude: %f\n", act_idx, myo_idx,
        //     endstop_force.norm());
    }
    return endstop_force; //force on myosin
}


bool apply_myosin_repulsion(const Filament& actin,
                            const Myosin& myosin,
                            int i,
                            int j,
                            const std::vector<double>& box,
                            int fix_myosin,
                            const utils::MoleculeConnection& actinIndicesPerMyosin,
                            double stiffness,
                            double max_force_cap,
                            vec& force_on_first,
                            vec& force_on_second)
{
    auto result = compute_myosin_repulsion(
        actin,
        myosin,
        i,
        j,
        box,
        fix_myosin,
        actinIndicesPerMyosin,
        stiffness,
        max_force_cap);
    if (!result.applied) {
        return false;
    }
    force_on_first += result.force_on_first;
    force_on_second += result.force_on_second;
    return true;
}

bool apply_actin_repulsion(const Filament& actin,
                           int i,
                           int j,
                           const std::vector<double>& box,
                           double crosslinker_length,
                           double stiffness,
                           double max_force_cap,
                           vec& force_on_first,
                           vec& force_on_second)
{
    auto result = compute_actin_repulsion(
        actin,
        i,
        j,
        box,
        crosslinker_length,
        stiffness,
        max_force_cap);
    if (!result.applied) {
        return false;
    }
    force_on_first += result.force_on_first;
    force_on_second += result.force_on_second;
    return true;
}

bool apply_actin_myosin_repulsion(const Filament& actin,
                                  const Myosin& myosin,
                                  int act_idx,
                                  int myo_idx,
                                  const std::vector<double>& box,
                                  double radius,
                                  double stiffness,
                                  double max_force_cap,
                                  vec& force_on_actin,
                                  vec& force_on_myosin)
{
    constexpr double EPS_FORCE = 1e-9;
    vec repulsive_force = compute_actin_myosin_repulsion(
        actin,
        myosin,
        act_idx,
        myo_idx,
        box,
        radius,
        stiffness,
        max_force_cap);
    if (repulsive_force.norm() <= EPS_FORCE) {
        return false;
    }
    force_on_myosin += repulsive_force;
    return true;
}
