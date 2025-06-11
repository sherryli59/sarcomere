#include "interaction.h"
#include <cstdio>
#include <cstdlib>  // For exit()

//---------------------------------------------------------------------
// Non-Template Function Definitions
//---------------------------------------------------------------------



real am_energy1(const ArrayXreal& center1, const double& length1, const ArrayXreal& dir1,
    const ArrayXreal& center2, const double& length2, const ArrayXreal& dir2,
    const std::vector<double>& box, const double k_am, const double kappa_am,
    const double myosin_radius)
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

    real dist = segment_segment_distance(a_arr, b_arr, c_arr, d_arr, box);

    // Angle and strength
    real dot_val = dir1[0]*dir2[0] + dir1[1]*dir2[1] + dir1[2]*dir2[2];
    real strength = abs(dot_val);
    dot_val = std::max(real(-1), std::min(real(1), dot_val));  // clamp
    real angle_energy = 0.5 * kappa_am * (1.0 - dot_val * dot_val); 

    if (dist > 0.8 * myosin_radius) {
    if (dist > myosin_radius) {
    printf("something's wrong. dist: %f\n", dist.val());
    exit(1);
    }
    return 0.5 * k_am * strength * dist * dist + angle_energy;
    } else {
    return 0.5 * (k_am / 10) * strength * dist * dist + angle_energy;
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
    const std::vector<double>& box, 
    const double k_aa, const double kappa_aa)
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

    real dist = segment_segment_distance(a_arr, b_arr, c_arr, d_arr, box);
    dist = dist - 0.03;

    // Compute angle between dir1 and dir2
    real dot_val = dir1[0]*dir2[0] + dir1[1]*dir2[1] + dir1[2]*dir2[2];
    dot_val = std::max(real(-1), std::min(real(1), dot_val));  // clamp for safety

    real angle;
    if (dot_val == 1) {
        angle = 0;
    } else if (dot_val == -1) {
        angle = M_PI;
    } else {
        angle = acos(dot_val);
    }
    real offset = abs(angle - M_PI);

    return 0.5 * (k_aa * dist * dist + kappa_aa * offset * offset);
}


std::vector<double> compute_aa_force_and_energy(Filament& actin,
                                                int& actin1_index, int& actin2_index,
                                                const std::vector<double>& box,
                                                const double k_aa, const double kappa_aa)
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
                                       center2, actin.length, dir2, box, k_aa, kappa_aa), u);
    forces.resize(forces_2.size());
    VectorXd::Map(&forces[0], forces_2.size()) = forces_2;
    return forces;
}

std::vector<double> compute_am_force_and_energy(Filament& actin, Myosin& myosin,
                                                int& actin_index, int& myosin_index,
                                                const std::vector<double>& box,
                                                const double k_am, const double kappa_am,
                                                const double myosin_radius)
{
    // Define 3D center positions
    ArrayXreal center1(3);
    center1 << actin.center[actin_index].x, actin.center[actin_index].y, actin.center[actin_index].z;
    ArrayXreal dir1(3);
    dir1 << actin.direction[actin_index].x, actin.direction[actin_index].y, actin.direction[actin_index].z;

    ArrayXreal center2(3);
    center2 << myosin.center[myosin_index].x, myosin.center[myosin_index].y, myosin.center[myosin_index].z;
    ArrayXreal dir2(3);
    dir2 << myosin.direction[myosin_index].x, myosin.direction[myosin_index].y, myosin.direction[myosin_index].z;

    real u;
    std::vector<double> forces;
    VectorXd forces_3; // 3D version of force vector

    if (k_am > 1e-6) {
        forces_3 = -gradient(am_energy1, wrt(center1, dir1, dir2),
                                at(center1, actin.length, dir1,
                                   center2, myosin.length, dir2, box, k_am, kappa_am, myosin_radius), u);
        forces.resize(forces_3.size());
        VectorXd::Map(&forces[0], forces_3.size()) = forces_3;
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



