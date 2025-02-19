#ifndef INTERACTION_H
#define INTERACTION_H

#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <limits>
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#include "components.h"  // Ensure this header defines Filament, Myosin, etc.

using Eigen::VectorXd;
using autodiff::ArrayXreal;
using autodiff::real;

//---------------------------------------------------------------------
// Template Functions (must be in the header)
//---------------------------------------------------------------------

// Case 1: For autodiff types
template<typename T, typename std::enable_if_t<!std::is_arithmetic_v<T>, int> = 0>
T smooth_round(T x) {
    return round(val(x));  // Use autodiff's round function
}

// Case 2: For standard numeric types
template<typename T, typename std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
T smooth_round(T x) {
    return std::round(x);  // Use standard round function
}
// // Function to calculate the norm of an array
// template<typename T>
// T norm(T* x, int dim) {
//     T norm_val = 0;
//     for (int i = 0; i < dim; i++) {
//         norm_val += x[i] * x[i];
//     }
//     return sqrt(norm_val);
// }

// // Function to calculate the norm of a vector
// template<typename T>
// T norm(const std::vector<T>& x) {
//     T norm_val = 0;
//     for (int i = 0; i < x.size(); i++) {
//         norm_val += x[i] * x[i];
//     }
//     return sqrt(norm_val);
// }


// Template function for computing the shortest distance between two segments with PBC in 3D.
// The segments are defined as:
//   Segment 1: P(s) = a + s*(b - a), with s in [0,1]
//   Segment 2: Q(t) = c + t*(d - c), with t in [0,1]
template<typename T>
T segment_segment_distance(const T* a, const T* b, const T* c, const T* d, const std::vector<double>& box) {
    double EPS = 1e-9;
    // Helper lambda: dot product of two 3D vectors.
    auto dot3 = [](const T* x, const T* y) -> T {
        return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
    };

    // Helper lambda: subtract two 3D vectors: out = x - y.
    auto subtract3 = [&](const T* x, const T* y, T* out) {
        out[0] = x[0] - y[0];
        out[1] = x[1] - y[1];
        out[2] = x[2] - y[2];
    };

    // Helper lambda: add two 3D vectors: out = x + y.
    auto add3 = [&](const T* x, const T* y, T* out) {
        out[0] = x[0] + y[0];
        out[1] = x[1] + y[1];
        out[2] = x[2] + y[2];
    };

    // Helper lambda: multiply vector x by scalar, store in out.
    auto multiply3 = [&](const T* x, T scalar, T* out) {
        out[0] = x[0] * scalar;
        out[1] = x[1] * scalar;
        out[2] = x[2] * scalar;
    };

    // Helper lambda: compute the norm (magnitude) of a 3D vector.
    auto norm3 = [&](const T* x) -> T {
        return sqrt(dot3(x, x));
    };

    // Helper lambda: Adjust a vector difference for periodic boundary conditions.
    auto pbc_adjust = [&](const T* diff, const std::vector<double>& box, T* out) {
        out[0] = diff[0] - box[0] * smooth_round(diff[0] / box[0]);
        out[1] = diff[1] - box[1] * smooth_round(diff[1] / box[1]);
        out[2] = diff[2] - box[2] * smooth_round(diff[2] / box[2]);
    };

    // Compute direction vectors for the segments.
    T u[3], v[3], w[3];
    subtract3(b, a, u);   // u = b - a
    subtract3(d, c, v);   // v = d - c
    subtract3(a, c, w);   // w = a - c

    // Apply PBC adjustments.
    T u_pbc[3], v_pbc[3], w_pbc[3];
    pbc_adjust(u, box, u_pbc);
    pbc_adjust(v, box, v_pbc);
    pbc_adjust(w, box, w_pbc);

    // Compute scalar coefficients.
    T a_val = dot3(u_pbc, u_pbc);
    T b_val = dot3(u_pbc, v_pbc);
    T c_val = dot3(v_pbc, v_pbc);
    T d_val = dot3(u_pbc, w_pbc);
    T e_val = dot3(v_pbc, w_pbc);
    T D_val = a_val * c_val - b_val * b_val;

    T s, t;
    if (abs(D_val) < EPS) {
        s = 0;
        t = (c_val > EPS ? e_val / c_val : 0);
    } else {
        s = (b_val * e_val - c_val * d_val) / D_val;
        t = (a_val * e_val - b_val * d_val) / D_val;
    }

    // Clamp parameters to [0,1].
    s = std::max((T)0, std::min((T)1, s));
    t = std::max((T)0, std::min((T)1, t));

    // Compute the closest points on the segments.
    T P[3], Q[3], s_u[3], t_v[3];
    multiply3(u_pbc, s, s_u);
    add3(a, s_u, P);
    multiply3(v_pbc, t, t_v);
    add3(c, t_v, Q);

    // Compute the difference vector between the closest points, then apply PBC.
    T dP[3];
    subtract3(P, Q, dP);
    T dP_pbc[3];
    pbc_adjust(dP, box, dP_pbc);

    return norm3(dP_pbc);
}


//---------------------------------------------------------------------
// Declarations of Non-Template Functions
//---------------------------------------------------------------------

// Compute the energy between two actin segments in 3D.
real aa_energy(const ArrayXreal& center1, const double& length1, const real& theta1, const real& phi1,
               const ArrayXreal& center2, const double& length2, const real& theta2, const real& phi2,
               const std::vector<double>& box, const double k_aa, const double kappa_aa);

// Compute the energy between an actin and a myosin segment in 3D.
real am_energy1(const ArrayXreal& center1, const double& length1, const real& theta1, const real& phi1,
                const ArrayXreal& center2, const double& length2, const real& theta2, const real& phi2,
                const std::vector<double>& box, const double k_am, const double kappa_am,
                const double myosin_radius);

// Compute the energy between two segments based solely on their angles in 3D.
real am_energy(const real& theta1, const real& phi1,
               const real& theta2, const real& phi2, const double kappa_am);

// Compute forces and energy for actin–actin interaction in 3D.
std::vector<double> compute_aa_force_and_energy(Filament& actin,
                                                int& actin1_index, int& actin2_index,
                                                const std::vector<double>& box,
                                                const double k_aa, const double kappa_aa);

// Compute forces and energy for actin–myosin interaction in 3D.
std::vector<double> compute_am_force_and_energy(Filament& actin, Myosin& myosin,
                                                int& actin_index, int& myosin_index,
                                                const std::vector<double>& box,
                                                const double k_am, const double kappa_am,
                                                const double myosin_radius, const bool turn_on_spring);

#endif  // INTERACTION_H
