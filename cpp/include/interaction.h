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
    
    // --- Global PBC Adjustment via Midpoints ---
    // Create local copies for A and B.
    T A_adj[3], B_adj[3];
    for (int i = 0; i < 3; ++i) {
        A_adj[i] = A[i];
        B_adj[i] = B[i];
    }
    
    // Compute midpoints of segment AB and segment CD.
    T midA[3], midCD[3];
    for (int i = 0; i < 3; ++i) {
        midA[i] = (A_adj[i] + B_adj[i]) / T(2);
        midCD[i] = (C[i] + D[i]) / T(2);
    }
    
    // Compute displacement between the midpoints.
    T disp[3];
    for (int i = 0; i < 3; ++i) {
        disp[i] = midA[i] - midCD[i];
    }
    
    // Compute the shift vector using the box dimensions.
    T shift[3];
    for (int i = 0; i < 3; ++i) {
        shift[i] = -T(box[i]) * smooth_round(disp[i] / T(box[i]));
    }
    
    // Apply the shift to the first segment's endpoints.
    for (int i = 0; i < 3; ++i) {
        A_adj[i] += shift[i];
        B_adj[i] += shift[i];
    }
    
    // --- End Global PBC Adjustment ---
    
    // Now compute the direction vectors using simple arithmetic.
    // u = B_adj - A_adj (direction of first segment)
    // v = D - C        (direction of second segment)
    // w = A_adj - C    (vector from C to adjusted A)
    T u[3], v[3], w[3];
    for (int i = 0; i < 3; ++i) {
        u[i] = B_adj[i] - A_adj[i];
        v[i] = D[i] - C[i];
        w[i] = A_adj[i] - C[i];
    }
    
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
    
    // Compute unconstrained optimum.
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
    
    // If the unconstrained optimum lies within [0,1]^2, evaluate it.
    if (t_opt >= 0 && t_opt <= 1 && s_opt >= 0 && s_opt <= 1) {
        evaluate(t_opt, s_opt);
    }
    
    // Evaluate boundaries: fix t = 0 and t = 1, optimizing s.
    for (int i = 0; i < 2; ++i) {
        T t_candidate = (i == 0) ? 0 : 1;
        T temp[3];
        for (int j = 0; j < 3; ++j)
            temp[j] = w[j] + u[j] * t_candidate;
        T s_candidate = clamp(dot3(v, temp) / c, T(0), T(1));
        evaluate(t_candidate, s_candidate);
    }
    
    // Evaluate boundaries: fix s = 0 and s = 1, optimizing t.
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
real aa_energy(const ArrayXreal& center1, const double& length1, 
    const ArrayXreal& dir1,
    const ArrayXreal& center2, const double& length2, 
    const ArrayXreal& dir2,
    const std::vector<double>& box, 
    const double k_aa, const double kappa_aa);

// Compute the energy between an actin and a myosin segment in 3D.
real am_energy1(const ArrayXreal& center1, const double& length1, const ArrayXreal& dir1,
    const ArrayXreal& center2, const double& length2, const ArrayXreal& dir2,
    const std::vector<double>& box, const double k_am, const double kappa_am,
    const double cutoff, const double optimal);

// Compute the energy between two segments based solely on their angles in 3D.
real am_energy(const ArrayXreal& dir1, const ArrayXreal& dir2, const double kappa_am);

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
                                                const double cutoff, const double optimal);

#endif  // INTERACTION_H
