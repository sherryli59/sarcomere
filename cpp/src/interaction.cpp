#include "interaction.h"
#include <cstdio>
#include <cstdlib>  // For exit()

//---------------------------------------------------------------------
// Non-Template Function Definitions
//---------------------------------------------------------------------



real am_energy1(const ArrayXreal& center1, const double& length1, const real& theta1, const real& phi1,
                const ArrayXreal& center2, const double& length2, const real& theta2, const real& phi2,
                const std::vector<double>& box, const double k_am, const double kappa_am,
                const double myosin_radius)
{
    // Compute endpoints for the first segment in 3D.
    real x1_start = center1[0] - (length1 / 2) * sin(phi1) * cos(theta1);
    real y1_start = center1[1] - (length1 / 2) * sin(phi1) * sin(theta1);
    real z1_start = center1[2] - (length1 / 2) * cos(phi1);
    real x1_end   = center1[0] + (length1 / 2) * sin(phi1) * cos(theta1);
    real y1_end   = center1[1] + (length1 / 2) * sin(phi1) * sin(theta1);
    real z1_end   = center1[2] + (length1 / 2) * cos(phi1);

    // Compute endpoints for the second segment in 3D.
    real x2_start = center2[0] - (length2 / 2) * sin(phi2) * cos(theta2);
    real y2_start = center2[1] - (length2 / 2) * sin(phi2) * sin(theta2);
    real z2_start = center2[2] - (length2 / 2) * cos(phi2);
    real x2_end   = center2[0] + (length2 / 2) * sin(phi2) * cos(theta2);
    real y2_end   = center2[1] + (length2 / 2) * sin(phi2) * sin(theta2);
    real z2_end   = center2[2] + (length2 / 2) * cos(phi2);

    // Create endpoint arrays (3D).
    real a[3] = {x1_start, y1_start, z1_start};
    real b[3] = {x1_end,   y1_end,   z1_end};
    real c[3] = {x2_start, y2_start, z2_start};
    real d[3] = {x2_end,   y2_end,   z2_end};

    // Compute the distance between segments.
    real dist = segment_segment_distance(a, b, c, d, box);

    // Compute strength as the absolute cosine of the angle between the two filaments.
    real dir1[3] = { sin(phi1)*cos(theta1), sin(phi1)*sin(theta1), cos(phi1) };
    real dir2[3] = { sin(phi2)*cos(theta2), sin(phi2)*sin(theta2), cos(phi2) };
    real dot_val = dir1[0]*dir2[0] + dir1[1]*dir2[1] + dir1[2]*dir2[2];
    real strength = abs(cos(acos(dot_val))); // essentially abs(dot_val)

    // Compute the angular difference between the filaments.
    real angle = acos(dot_val);
    // Adjust angle to minimal value.
    angle = (angle < M_PI - angle) ? angle : (M_PI - angle);

    if (dist > 0.8 * myosin_radius) {
        if (dist > myosin_radius) {
            // Optionally print debug messages.
            printf("dist: %f\n", dist.val());
            // Additional debugging info...
        }
        return 0.5 * k_am * strength * dist * dist + 0.5 * kappa_am * angle * angle;
    } else {
        real energy = 0.5 * (k_am / 10) * strength * dist * dist + 0.5 * kappa_am * angle * angle;
        return energy;
    }
}

real am_energy(const real& theta1, const real& phi1,
               const real& theta2, const real& phi2, const double kappa_am)
{
    // Compute unit direction vectors in spherical coordinates.
    real dir1_x = sin(phi1) * cos(theta1);
    real dir1_y = sin(phi1) * sin(theta1);
    real dir1_z = cos(phi1);

    real dir2_x = sin(phi2) * cos(theta2);
    real dir2_y = sin(phi2) * sin(theta2);
    real dir2_z = cos(phi2);

    real dot_val = dir1_x * dir2_x + dir1_y * dir2_y + dir1_z * dir2_z;
    real angle = acos(dot_val);  // angle in [0, π]
    return 0.5 * kappa_am * angle * angle;
}


real aa_energy(const ArrayXreal& center1, const double& length1, 
               const real& theta1, const real& phi1,
               const ArrayXreal& center2, const double& length2, 
               const real& theta2, const real& phi2,
               const std::vector<double>& box, 
               const double k_aa, const double kappa_aa)
{
    // Compute endpoints for filament 1 in 3D.
    real x1_start = center1[0] - (length1/2)*sin(phi1)*cos(theta1);
    real y1_start = center1[1] - (length1/2)*sin(phi1)*sin(theta1);
    real z1_start = center1[2] - (length1/2)*cos(phi1);
    real x1_end   = center1[0] + (length1/2)*sin(phi1)*cos(theta1);
    real y1_end   = center1[1] + (length1/2)*sin(phi1)*sin(theta1);
    real z1_end   = center1[2] + (length1/2)*cos(phi1);

    // Compute endpoints for filament 2 in 3D.
    real x2_start = center2[0] - (length2/2)*sin(phi2)*cos(theta2);
    real y2_start = center2[1] - (length2/2)*sin(phi2)*sin(theta2);
    real z2_start = center2[2] - (length2/2)*cos(phi2);
    real x2_end   = center2[0] + (length2/2)*sin(phi2)*cos(theta2);
    real y2_end   = center2[1] + (length2/2)*sin(phi2)*sin(theta2);
    real z2_end   = center2[2] + (length2/2)*cos(phi2);

    // Set up endpoints as 3D arrays.
    real a[3] = {x1_start, y1_start, z1_start};
    real b[3] = {x1_end,   y1_end,   z1_end};
    real c[3] = {x2_start, y2_start, z2_start};
    real d[3] = {x2_end,   y2_end,   z2_end};

    // Compute distance between segments using a 3D segment-segment function with PBC.
    real dist = segment_segment_distance(a, b, c, d, box);
    dist = dist - 0.03;

    // Compute direction vectors from the angles.
    real dir1[3] = { sin(phi1)*cos(theta1), sin(phi1)*sin(theta1), cos(phi1) };
    real dir2[3] = { sin(phi2)*cos(theta2), sin(phi2)*sin(theta2), cos(phi2) };

    real dot_val = dir1[0]*dir2[0] + dir1[1]*dir2[1] + dir1[2]*dir2[2];
    real angle = acos(dot_val);

    return 0.5 * (k_aa * dist * dist + kappa_aa * angle * angle);
}

std::vector<double> compute_aa_force_and_energy(Filament& actin,
                                                int& actin1_index, int& actin2_index,
                                                const std::vector<double>& box,
                                                const double k_aa, const double kappa_aa)
{
    // Construct 3D centers.
    ArrayXreal center1(3);
    center1 << actin.center[actin1_index].x, actin.center[actin1_index].y, actin.center[actin1_index].z;
    real theta1 = actin.theta[actin1_index];
    real phi1   = actin.phi[actin1_index];

    ArrayXreal center2(3);
    center2 << actin.center[actin2_index].x, actin.center[actin2_index].y, actin.center[actin2_index].z;
    real theta2 = actin.theta[actin2_index];
    real phi2   = actin.phi[actin2_index];

    real u;
    std::vector<double> forces;
    // Compute the gradient of aa_energy with respect to center1, theta1, phi1, theta2, and phi2.
    VectorXd forces_2 = -gradient(aa_energy, wrt(center1, theta1, phi1, theta2, phi2),
                                    at(center1, actin.length, theta1, phi1,
                                       center2, actin.length, theta2, phi2, box, k_aa, kappa_aa), u);
    forces.resize(forces_2.size());
    VectorXd::Map(&forces[0], forces_2.size()) = forces_2;
    return forces;
}

std::vector<double> compute_am_force_and_energy(Filament& actin, Myosin& myosin,
                                                int& actin_index, int& myosin_index,
                                                const std::vector<double>& box,
                                                const double k_am, const double kappa_am,
                                                const double myosin_radius, const bool turn_on_spring)
{
    // Define 3D center positions
    ArrayXreal center1(3);
    center1 << actin.center[actin_index].x, actin.center[actin_index].y, actin.center[actin_index].z;
    real theta1 = actin.theta[actin_index];
    real phi1 = actin.phi[actin_index]; // New angle for 3D

    ArrayXreal center2(3);
    center2 << myosin.center[myosin_index].x, myosin.center[myosin_index].y, myosin.center[myosin_index].z;
    real theta2 = myosin.theta[myosin_index];
    real phi2 = myosin.phi[myosin_index]; // New angle for 3D

    real u;
    std::vector<double> forces;
    VectorXd forces_3; // 3D version of force vector

    if (turn_on_spring) {
        forces_3 = -gradient(am_energy1, wrt(center1, theta1, phi1, theta2, phi2),
                             at(center1, actin.length, theta1, phi1,
                                center2, myosin.length, theta2, phi2, box, k_am, kappa_am, myosin_radius), u);
        forces.resize(forces_3.size());
        VectorXd::Map(&forces[0], forces_3.size()) = forces_3;
    }
    else {
        forces_3 = -gradient(am_energy, wrt(theta1, phi1, theta2, phi2),
                             at(theta1, phi1, theta2, phi2, kappa_am), u);
        // Prepend three zeros for the force vector (x, y, z)
        forces.resize(forces_3.size() + 3);
        forces[0] = 0;
        forces[1] = 0;
        forces[2] = 0;
        VectorXd::Map(&forces[3], forces_3.size()) = forces_3;
    }
    return forces;
}
