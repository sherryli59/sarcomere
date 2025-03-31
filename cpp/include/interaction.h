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

const autodiff::real autodiff_infinity = autodiff::real(std::numeric_limits<double>::infinity());

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


// Helper clamp function.
template<typename T>
T clamp(T x, T lo, T hi) {
    return max(lo, min(x, hi));
}

template<typename T>
T segment_segment_distance(const T* A, const T* B, 
                           const T* C, const T* D, 
                           const std::vector<double>& box) {
    const double EPS = 1e-9;
    // Helper lambda: adjust a point to its minimum image
    auto pbc_diff = [](const T* x,const T* y, const std::vector<double>& box, T* out) {
        out[0] = x[0] - y[0];
        out[1] = x[1] - y[1];
        out[2] = x[2] - y[2];
        out[0] -= T(box[0]) * smooth_round(out[0] / T(box[0]));
        out[1] -= T(box[1]) * smooth_round(out[1] / T(box[1]));
        out[2] -= T(box[2]) * smooth_round(out[2] / T(box[2]));
    };
    // Compute direction vectors with periodic boundary conditions.
    // u = pbc_diff(B, A, box)
    // v = pbc_diff(D, C, box)
    // w = pbc_diff(A, C, box)
    T u[3], v[3], w[3];
    pbc_diff(B, A, box, u);
    pbc_diff(D, C, box, v);
    pbc_diff(A, C, box, w);

    // Lambda for 3D dot product.
    auto dot3 = [](const T* x, const T* y) -> T {
        return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
    };

    // Lambda for computing squared norm.
    auto norm_sq = [&](const T* x) -> T {
        return dot3(x, x);
    };

    // Compute coefficients.
    T a = dot3(u, u);      // |u|²
    T b = dot3(u, v);
    T c = dot3(v, v);      // |v|²
    T d_val = dot3(u, w);
    T e_val = dot3(v, w);
    T denom = a * c - b * b;

    // Unconstrained optimum:
    T t_opt = 0, s_opt = 0;
    if (abs(denom) > EPS) {
        t_opt = (b * e_val - c * d_val) / denom;
        s_opt = (a * e_val - b * d_val) / denom;
    } else {
        // Nearly parallel: set t_opt = 0 and optimize s.
        t_opt = 0;
        s_opt = (c > EPS ? dot3(v, w) / c : 0);
    }

    // Evaluate candidate (t,s) by computing squared distance.
    T best_dist2 = autodiff_infinity;
    T best_t = 0, best_s = 0;
    auto evaluate = [&](T t, T s) {
        T diff[3];
        for (int i = 0; i < 3; ++i)
            diff[i] = w[i] + u[i] * t - v[i] * s;
        T d2 = norm_sq(diff);
        if (d2 < best_dist2) {
            best_dist2 = d2;
            best_t = t;
            best_s = s;
        }
    };

    // If unconstrained optimum lies within [0,1]^2, evaluate it.
    if (t_opt >= 0 && t_opt <= 1 && s_opt >= 0 && s_opt <= 1) {
        evaluate(t_opt, s_opt);
    }
    // Otherwise, check boundaries.
    // For t fixed at 0 and 1, optimize s.
    for (int i = 0; i < 2; ++i) {
        T t_candidate = (i == 0) ? 0 : 1;
        T temp[3];
        for (int j = 0; j < 3; ++j)
            temp[j] = w[j] + u[j] * t_candidate;
        T s_candidate = clamp(dot3(v, temp) / c, T(0), T(1));
        evaluate(t_candidate, s_candidate);
    }

    // For s fixed at 0 and 1, optimize t.
    for (int j = 0; j < 2; ++j) {
        T s_candidate = (j == 0) ? 0 : 1;
        T temp[3];
        for (int i = 0; i < 3; ++i)
            temp[i] = s_candidate * v[i] - w[i];
        T t_candidate = clamp(dot3(u, temp) / a, T(0), T(1));
        evaluate(t_candidate, s_candidate);
    }
    return sqrt(best_dist2);
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
                                                const double myosin_radius);

#endif  // INTERACTION_H
