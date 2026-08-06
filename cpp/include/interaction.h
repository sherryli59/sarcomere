#ifndef INTERACTION_H
#define INTERACTION_H

#include <array>
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
                           const std::vector<double>& box,
                           const std::array<bool,3>& periodic) {
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
        double L = (i < static_cast<int>(box.size())) ? box[i] : 0.0;
        if (periodic[i] && L > 0.0) {
            shift[i] = -T(L) * smooth_round(disp[i] / T(L));
        } else {
            shift[i] = T(0);
        }
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
    const std::vector<double>& box, const std::array<bool,3>& periodic, 
    const double k_aa, const double kappa_aa,
    const double cutoff, const double optimal,
    const double antiparallel_overlap);

// Compute the energy between an actin and a myosin segment in 3D.
real am_energy1(const ArrayXreal& center1, const double& length1, const ArrayXreal& dir1,
    const ArrayXreal& center2, const double& length2, const ArrayXreal& dir2,
    const std::vector<double>& box, const std::array<bool,3>& periodic, const double k_am, const double kappa_am,
    const double cutoff, const double optimal);

// Compute the energy between two segments based solely on their angles in 3D.
real am_energy(const ArrayXreal& dir1, const ArrayXreal& dir2, const double kappa_am);

// Compute forces and energy for actin–actin interaction in 3D.
std::vector<double> compute_aa_force_and_energy(Filament& actin,
                                                int& actin1_index, int& actin2_index,
                                                const std::vector<double>& box,
                                                const double k_aa, const double kappa_aa,
                                                const double cutoff, const double optimal,
                                                const double antiparallel_overlap = 0.0);

struct LocalAALinkResult {
    vec force_on_first {0.0, 0.0, 0.0};
    vec torque_on_first {0.0, 0.0, 0.0};
    vec torque_on_second {0.0, 0.0, 0.0};
    vec connector {0.0, 0.0, 0.0};
    double length = 0.0;
    double energy = 0.0;
    double registry_error = 0.0;
    double direction_dot = 1.0;
};

// One implicit ACTN2 molecule attached at fixed material coordinates on two
// actin rods. The effective connector joins the rod centerlines; its rest
// length therefore includes the actin radii and projected ACTN2 span.
LocalAALinkResult compute_local_aa_link(
    const Filament& actin, int first, int second,
    double s_first, double s_second,
    const std::vector<double>& box,
    double radial_stiffness, double perpendicular_stiffness,
    double alignment_stiffness,
    double preferred_centerline,
    double antiparallel_overlap = 0.0,
    double axial_stiffness_ratio = 1.0,
    double axial_saturation_error = 0.05,
    // >0 switches the barbed-end register term from the Gaussian-attenuated
    // form (whose restoring force vanishes far from the target) to a
    // saturating form whose force approaches this ceiling in pN. The
    // small-error stiffness is unchanged, so a registered pair behaves
    // exactly as before; only badly mis-registered pairs feel a difference.
    double axial_force_max = 0.0);


// Compute forces and energy for actin–myosin interaction in 3D.
std::vector<double> compute_am_force_and_energy(Filament& actin, Myosin& myosin,
                                                int& actin_index, int& myosin_index,
                                                const std::vector<double>& box,
                                                const double k_am, const double kappa_am,
                                                const double cutoff, const double optimal);

struct AMTetherResult {
    vec force_on_actin {0.0, 0.0, 0.0};
    vec force_on_myosin {0.0, 0.0, 0.0};
    vec actin_attachment {0.0, 0.0, 0.0};
    vec myosin_attachment {0.0, 0.0, 0.0};
    double transverse_distance = 0.0;
    double centerline_distance = 0.0;
    double surface_gap = 0.0;
    double head_surface_distance = 0.0;
    // Full 3-D distance from the spatial myosin-head anchor to the closest
    // finite actin centerline point, and the corresponding surface distance.
    // Unlike transverse_distance these retain axial mismatch and therefore
    // determine whether this particular material head can reach the actin.
    double head_connector_distance = 0.0;
    double head_reach_surface_distance = 0.0;
    // Projection of the material head anchor onto the infinite actin axis.
    // Association may only use this coordinate when it lies on the finite
    // actin segment; actin_coordinate remains the clamped coordinate used to
    // construct endpoint geometry for already-bound heads and diagnostics.
    double unbounded_actin_coordinate = 0.0;
    bool projection_inside_actin = false;
    double actin_coordinate = 0.0;
};

// Passive proximity-only AM tether. It removes the common filament-axis
// component from the closest-point connector so it does not propel or brake
// longitudinal sliding.
AMTetherResult compute_am_proximity_tether(const Filament& actin,
                                           const Myosin& myosin,
                                           int actin_index,
                                           int myosin_index,
                                           const std::vector<double>& box,
                                           double stiffness,
                                           double optimal,
                                           bool directional = false);

// Proximity tether for one explicit head rooted at a fixed material coordinate
// along the myosin axis. The actin attachment is the closest point on its
// segment; only transverse separation is penalized.
AMTetherResult compute_am_head_tether(const Filament& actin,
                                      const Myosin& myosin,
                                      int actin_index,
                                      int myosin_index,
                                      const std::vector<double>& box,
                                      double myosin_anchor_axial,
                                      double stiffness,
                                      double optimal,
                                      double myosin_anchor_azimuth = 0.0,
                                      double myosin_radial_offset = 0.0,
                                      const vec* material_e1 = nullptr,
                                      const vec* material_e2 = nullptr,
                                      double actin_radius = 0.0,
                                      double radial_slack = 0.0);

// Energy-consistent explicit cross-bridge.  The myosin anchor is fixed in the
// transported thick-filament material frame and actin_coordinate is fixed in
// the thin-filament material frame for the lifetime of an attachment.
// Forces and physical torques are negative derivatives of one objective
// two-spring energy, so an isolated head conserves total force and torque.
struct AMTwoSpringHeadResult {
    vec actin_attachment {0.0, 0.0, 0.0};
    vec myosin_attachment {0.0, 0.0, 0.0};
    vec connector {0.0, 0.0, 0.0};
    vec radial_force_on_actin {0.0, 0.0, 0.0};
    vec axial_force_on_actin {0.0, 0.0, 0.0};
    vec force_on_actin {0.0, 0.0, 0.0};
    vec force_on_myosin {0.0, 0.0, 0.0};
    vec physical_torque_on_actin {0.0, 0.0, 0.0};
    vec physical_torque_on_myosin {0.0, 0.0, 0.0};
    double energy = 0.0;
    double transverse_distance = 0.0;
    double head_reach_surface_distance = 0.0;
    double longitudinal_extension = 0.0;
    double productive_strain = 0.0;
    double active_force = 0.0;
    double radial_energy = 0.0;
    double axial_energy = 0.0;
};

AMTwoSpringHeadResult compute_am_two_spring_head(
    const Filament& actin, const Myosin& myosin,
    int actin_index, int myosin_index, const std::vector<double>& box,
    double actin_coordinate, double myosin_anchor_axial,
    double myosin_anchor_azimuth, double myosin_radial_offset,
    const vec* material_e1, const vec* material_e2,
    double radial_preferred, double radial_stiffness, double radial_slack,
    double axial_stiffness, double stroke_progress, double stall_force,
    bool strong, bool bidirectional_strong = false);

// Choose the closest azimuth reachable from a fixed crown anchor. The result
// is clamped to crown_azimuth +/- max_swing and returned in radians.
double select_am_binding_azimuth(
    const Filament& actin, const Myosin& myosin,
    int actin_index, int myosin_index, const std::vector<double>& box,
    double myosin_anchor_axial, double crown_azimuth, double max_swing,
    const vec& material_e1, const vec& material_e2);

vec compute_actin_myosin_radial_repulsion(
    const Filament& actin, const Myosin& myosin,
    int actin_index, int myosin_index,
    const std::vector<double>& box, double excluded_centerline_distance,
    double stiffness);

struct AMOffCenterTorqueResult {
    vec lever_arm {0.0, 0.0, 0.0};
    vec physical_torque {0.0, 0.0, 0.0};
    // Tangent generalized force used by the unit-vector orientation dynamics.
    vec orientation_drive {0.0, 0.0, 0.0};
};

// Convert a force applied at an off-center cross-bridge attachment into the
// physical moment tau = (r_head-r_center) x F and the corresponding tangent
// drive du/dt ~ tau x axis expected by the existing orientation integrator.
AMOffCenterTorqueResult compute_am_offcenter_torque(
    const vec& center,
    const vec& axis,
    const vec& attachment,
    const vec& force,
    const std::vector<double>& box,
    const std::array<bool,3>& periodic);

// Positive values rotate the moving nematic axis along the shortest tangent
// direction toward the target axis; negative values increase misalignment.
double project_nematic_alignment_drive(
    const vec& moving_axis, const vec& target_axis,
    const vec& orientation_drive);

// Weak aggregate angular compliance of an occupied AM pair. Explicit heads
// still supply the active off-center torque; this flat-bottom potential only
// represents passive resistance of the occupied ensemble to large axial
// misalignment.
struct AMNematicTetherResult {
    vec drive_on_actin {0.0, 0.0, 0.0};
    vec drive_on_myosin {0.0, 0.0, 0.0};
    double angle = 0.0;
    double energy = 0.0;
    bool active = false;
};

AMNematicTetherResult compute_am_nematic_flat_bottom_tether(
    const vec& actin_axis, const vec& myosin_axis,
    double stiffness, double slack_angle);

// Result container for pairwise repulsion calculations.
struct RepulsionResult {
    bool applied = false;
    vec force_on_first {0.0, 0.0, 0.0};
    vec force_on_second {0.0, 0.0, 0.0};
    // Tangent director drives generated by applying the forces at the closest
    // points on the two finite rods.  The orientation integrator consumes
    // tau x axis rather than the bare physical moment tau.
    vec torque_on_first {0.0, 0.0, 0.0};
    vec torque_on_second {0.0, 0.0, 0.0};
    double segment_distance = std::numeric_limits<double>::infinity();
    // Populated by the polarity-aware pointed-end contact.  The signed axial
    // gap is positive while the two pointed-end caps remain separated and
    // negative after a finite-step attempted crossing.  Penetration is the
    // unilateral-contact violation after accounting for transverse offset.
    double axial_gap = std::numeric_limits<double>::infinity();
    double transverse_distance = std::numeric_limits<double>::infinity();
    double contact_penetration = 0.0;
};

// Position-level hard-contact correction used after a Brownian/Langevin
// update.  The soft steric potential remains responsible for resolved
// forces and torques; this projection is the overdamped unilateral-contact
// constraint that prevents a finite time step or a force cap from allowing
// two thick-filament capsules to pass through one another.
struct ContactProjectionResult {
    bool applied = false;
    double distance_before = std::numeric_limits<double>::infinity();
    double correction = 0.0;
};

ContactProjectionResult project_myosin_hard_contact(
    Myosin& myosin, int i, int j, const std::vector<double>& box,
    int fix_myosin, double target_distance,
    double translation_weight = 1.0, double rotation_weight = 0.0);

RepulsionResult compute_actin_myosin_radial_repulsion_result(
    const Filament& actin, const Myosin& myosin,
    int actin_index, int myosin_index,
    const std::vector<double>& box, double excluded_centerline_distance,
    double stiffness);

// Myosin–myosin segment repulsion helper.
RepulsionResult compute_myosin_repulsion(const Filament& actin,
                                         const Myosin& myosin,
                                         int i,
                                         int j,
                                         const std::vector<double>& box,
                                         int fix_myosin,
                                         const utils::MoleculeConnection& actinIndicesPerMyosin,
                                         double stiffness);

// Actin–actin segment repulsion helper.
RepulsionResult compute_actin_repulsion(const Filament& actin,
                                        int i,
                                        int j,
                                        const std::vector<double>& box,
                                        double crosslinker_length,
                                        double stiffness);

// Local, status-independent contact between the pointed ends of two
// antiparallel thin filaments occupying the same lattice lane.  Unlike the
// generic capsule repulsion, this retains the signed end-to-end ordering, so
// a finite integration step cannot turn an end-on collision into an arbitrary
// transverse force after the tips cross.  epsilon is the WCA energy scale
// (pN um), contact_diameter is the decorated cap diameter (um), and force_cap
// is an optional coarse-grained buckling/load-sharing ceiling (pN; <=0 means
// uncapped).
RepulsionResult compute_opposing_pointed_end_contact(
    const Filament& actin, int i, int j,
    const std::vector<double>& box,
    double contact_diameter, double epsilon,
    double force_cap = 0.0,
    double antiparallel_dot_max = -0.8);

// Actin–myosin axial end-stop repulsion helper.
vec compute_actin_myosin_repulsion(const Filament& actin,
                                   const Myosin& myosin,
                                   int act_idx,
                                   int myo_idx,
                                   const std::vector<double>& box,
                                   double radius,
                                   double stiffness);

// Convenience wrappers that accumulate repulsion forces.
bool apply_myosin_repulsion(const Filament& actin,
                            const Myosin& myosin,
                            int i,
                            int j,
                            const std::vector<double>& box,
                            int fix_myosin,
                            const utils::MoleculeConnection& actinIndicesPerMyosin,
                            double stiffness,
                            vec& force_on_first,
                            vec& force_on_second,
                            double& segment_distance);

bool apply_actin_repulsion(const Filament& actin,
                           int i,
                           int j,
                           const std::vector<double>& box,
                           double crosslinker_length,
                           double stiffness,
                           vec& force_on_first,
                           vec& force_on_second);

bool apply_actin_myosin_repulsion(const Filament& actin,
                                  const Myosin& myosin,
                                  int act_idx,
                                  int myo_idx,
                                  const std::vector<double>& box,
                                  double radius,
                                  double stiffness,
                                  vec& force_on_actin,
                                  vec& force_on_myosin);

#endif  // INTERACTION_H
