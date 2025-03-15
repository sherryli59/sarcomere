#include "geometry.h"
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <utility>
#include <tuple>
#include <string>


namespace geometry {

const double EPS = 1e-9;


// Helper function to compute the difference between two vectors using PBC
vec pbc_diff(const vec& a, const vec& b, const std::vector<double>& box) {
    vec diff = a - b;
    diff.pbc_wrap(box); // Adjusts diff.x, diff.y, and diff.z using the minimal image convention
    return diff;
}


// A simple clamp helper
double clamp(double x, double lo, double hi) {
    return std::max(lo, std::min(x, hi));
}

// This function computes the shortest distance between segments AB and CD by solving
// the constrained optimization problem: minimize || (A+t*(B-A)) - (C+s*(D-C)) ||^2 subject to t,s in [0,1].
double segment_segment_distance(const vec& A, const vec& B, 
                                              const vec& C, const vec& D, 
                                              const std::vector<double>& box) {
    // Compute direction vectors (PBC applied)
    vec u = pbc_diff(B, A, box);  // Direction of segment AB
    vec v = pbc_diff(D, C, box);  // Direction of segment CD
    vec w = pbc_diff(A, C, box);  // Vector from C to A

    double a = u.dot(u);         // Squared length of AB
    double b = u.dot(v);
    double c = v.dot(v);         // Squared length of CD
    double d_val = u.dot(w);
    double e_val = v.dot(w);

    double denom = a * c - b * b; // Denominator in the unconstrained solution

    // Unconstrained optimum:
    double t_opt = 0.0, s_opt = 0.0;
    if (std::fabs(denom) > EPS) {
        t_opt = (b * e_val - c * d_val) / denom;
        s_opt = (a * e_val - b * d_val) / denom;
    } else {
        // Segments are nearly parallel: choose t_opt = 0 and optimize s.
        t_opt = 0.0;
        s_opt = (v.dot(w)) / c;
    }

    // We'll check the unconstrained optimum and also boundary candidates.
    double best_dist2 = std::numeric_limits<double>::infinity();
    double best_t = 0.0, best_s = 0.0;

    auto evaluate = [&](double t, double s) {
        vec diff = w + u * t - v * s;
        double d2 = diff.norm_squared();
        if (d2 < best_dist2) {
            best_dist2 = d2;
            best_t = t;
            best_s = s;
        }
    };

    // If unconstrained optimum lies within [0,1]^2, use it.
    if (t_opt >= 0.0 && t_opt <= 1.0 && s_opt >= 0.0 && s_opt <= 1.0) {
        evaluate(t_opt, s_opt);
    }

    // Otherwise, check boundaries.
    // Check t = 0 and t = 1, and optimize s.
    for (int i = 0; i < 2; i++) {
        double t_candidate = (i == 0) ? 0.0 : 1.0;
        // For fixed t, f(s) = || w + t*u - s*v ||^2 is quadratic in s.
        // The optimum (unconstrained) is s = (v dot (w + t*u)) / c.
        double s_candidate = clamp(v.dot(w + u * t_candidate) / c, 0.0, 1.0);
        evaluate(t_candidate, s_candidate);
    }
    // Check s = 0 and s = 1, and optimize t.
    for (int j = 0; j < 2; j++) {
        double s_candidate = (j == 0) ? 0.0 : 1.0;
        // For fixed s, f(t) = || w + t*u - s*v ||^2 is quadratic in t.
        // The optimum is t = (u dot (s*v - w)) / a.
        double t_candidate = clamp(u.dot(s_candidate * v - w) / a, 0.0, 1.0);
        evaluate(t_candidate, s_candidate);
    }

    return std::sqrt(best_dist2);
}



// Compute the shortest distance between two line segments by solving
// the constrained minimization problem.
// Returns a pair: (distance, { "start": P, "end": Q, "normal": (P-Q) } ).
std::pair<double, std::map<std::string, vec>> segment_segment_distance_w_normal(
    const vec& A, const vec& B, 
    const vec& C, const vec& D, 
    const std::vector<double>& box)
{
    // Compute direction vectors (with PBC adjustments)
    vec u = pbc_diff(B, A, box);  // AB direction
    vec v = pbc_diff(D, C, box);  // CD direction
    vec w = pbc_diff(A, C, box);  // From C to A

    // Squared lengths and dot products.
    double a = u.dot(u);         // |u|^2
    double b = u.dot(v);
    double c = v.dot(v);         // |v|^2
    double d_val = u.dot(w);
    double e_val = v.dot(w);

    double denom = a * c - b * b; // Denom for unconstrained optimum

    // Unconstrained optimum parameters:
    double t_opt = 0.0, s_opt = 0.0;
    if (std::fabs(denom) > EPS) {
        t_opt = (b * e_val - c * d_val) / denom;
        s_opt = (a * e_val - b * d_val) / denom;
    } else {
        // Segments nearly parallel: fix t = 0 and optimize s.
        t_opt = 0.0;
        s_opt = clamp(v.dot(w) / c, 0.0, 1.0);
    }
    
    // We'll search for the best (t,s) in [0,1]^2 by evaluating:
    double best_dist2 = std::numeric_limits<double>::infinity();
    double best_t = 0.0, best_s = 0.0;
    
    auto evaluate = [&](double t, double s) {
        vec diff = w + u * t - v * s;  // P(t) - Q(s)
        double d2 = diff.norm_squared();
        if (d2 < best_dist2) {
            best_dist2 = d2;
            best_t = t;
            best_s = s;
        }
    };
    
    // If unconstrained optimum lies in [0,1]^2, evaluate it.
    if (t_opt >= 0.0 && t_opt <= 1.0 && s_opt >= 0.0 && s_opt <= 1.0) {
        evaluate(t_opt, s_opt);
    }
    
    // Evaluate boundaries:
    // 1. t = 0 and t = 1, optimizing over s.
    for (int i = 0; i < 2; i++) {
        double t_candidate = (i == 0) ? 0.0 : 1.0;
        double s_candidate = clamp(v.dot(w + u * t_candidate) / c, 0.0, 1.0);
        evaluate(t_candidate, s_candidate);
    }
    // 2. s = 0 and s = 1, optimizing over t.
    for (int j = 0; j < 2; j++) {
        double s_candidate = (j == 0) ? 0.0 : 1.0;
        double t_candidate = clamp(u.dot(s_candidate * v - w) / a, 0.0, 1.0);
        evaluate(t_candidate, s_candidate);
    }
    
    double distance = std::sqrt(best_dist2);
    
    // Compute the closest points on each segment:
    vec P = A + u * best_t;  // Closest point on AB
    vec Q = C + v * best_s;  // Closest point on CD
    
    // Compute the normal vector from Q to P (apply PBC if needed).
    vec normal = pbc_diff(P, Q, box);
    
    // Package the extra information in a map.
    std::map<std::string, vec> info;
    info["start"] = P;
    info["end"] = Q;
    info["normal"] = normal;
    
    return std::make_pair(distance, info);
}


void apply_pbc(vec& actin_left, vec& actin_right, 
                        vec& myosin_left, vec& myosin_right, std::vector<double> box){

    vec actin_mid = (actin_left + actin_right) / 2;
    vec myosin_mid = (myosin_left + myosin_right) / 2;
    vec displacement = actin_mid - myosin_mid;
    vec shift;
    shift.x =  -box[0]*round(displacement.x/box[0]);
    shift.y =  - box[1]*round(displacement.y/box[1]);
    shift.z =  - box[2]*round(displacement.z/box[2]);
    actin_left = actin_left + shift;
    actin_right = actin_right + shift;
}




// Filter roots within the interval [t_start, t_end]
std::vector<double> findRootsInInterval(const std::vector<double>& roots, double t_start, double t_end) {
    std::vector<double> valid_roots;
    for (double t : roots) {
        if (t >= t_start - EPS && t <= t_end + EPS) {
            valid_roots.push_back(t);
        }
    }
    return valid_roots;
}




// Solve quadratic equation: a*t^2 + b*t + c = 0
std::vector<double> solveQuadratic(double a, double b, double c) {
    std::vector<double> roots;
    double discriminant = b * b - 4 * a * c;

    if (discriminant < -EPS) {
        // No real roots
        return roots;
    } else if (std::abs(discriminant) < EPS) {
        // One real root
        double t = -b / (2 * a);
        roots.push_back(t);
    } else {
        // Two real roots
        double sqrt_disc = std::sqrt(discriminant);
        double t1 = (-b - sqrt_disc) / (2 * a);
        double t2 = (-b + sqrt_disc) / (2 * a);
        roots.push_back(t1);
        roots.push_back(t2);
    }
    return roots;
}

std::vector<double> compute_transitions(vec A, vec B, vec C, vec D) {
    vec AB_dir = B - A;
    vec CD_dir = D - C;
    double CD_len_sq = CD_dir.dot(CD_dir);
    
    // Avoid division by zero if CD is a point
    if (CD_len_sq < EPS) return {0.0, 1.0};

    // Compute s(t) coefficients
    vec AC = A - C;
    double numerator_const = AC.dot(CD_dir);
    double numerator_t_coeff = AB_dir.dot(CD_dir);
    
    // Case 1: s(t) = 0
    double t0 = -numerator_const / numerator_t_coeff;
    
    // Case 2: s(t) = 1
    double t1 = (CD_len_sq - numerator_const) / numerator_t_coeff;
    
    // Case 3: Equidistant to C and D
    double dist_C_sq = (A - C).norm_squared();
    double dist_D_sq = (A - D).norm_squared();
    double t2 = (dist_D_sq - dist_C_sq) / (2 * AB_dir.dot(CD_dir));

    // Collect and clamp boundaries
    std::vector<double> boundaries = {0.0, 1.0};
    auto add_clamped = [&](double t) {
        if (t > 0 && t < 1) boundaries.push_back(std::clamp(t, 0.0, 1.0));
    };
    add_clamped(t0);
    add_clamped(t1);
    add_clamped(t2);

    // Remove duplicates
    std::sort(boundaries.begin(), boundaries.end());
    boundaries.erase(std::unique(boundaries.begin(), boundaries.end()), boundaries.end());
    return boundaries;
}

// Main function to find intervals on AB where distance to CD is less than d
std::tuple<double,vec,vec> subsegment_within_distance(vec A, vec B, vec C, vec D, double d) {
    std::vector<std::pair<vec, vec>> intervals;
    double interval_start = 1.0;
    double interval_end = 0.0;
    vec dirAB = B - A; // Direction vector of AB
    vec dirCD = D - C; // Direction vector of CD
    double lenAB_sq = dirAB.norm_squared();
    double lenCD_sq = dirCD.norm_squared();

    //rule out the case where the midpoints of the two segments are more than d apart
    vec AB_mid = (A + B) / 2;
    vec CD_mid = (C + D) / 2;
    if ((AB_mid - CD_mid).norm() > dirAB.norm() + dirCD.norm() + d) {
        return std::make_tuple(0.0, A, A);
    }
    // Find transition points where closest point on CD changes
    std::vector<double> boundaries = {0.0, 1.0};

    if (lenAB_sq > EPS) {
        // Compute transition points
        boundaries = compute_transitions(A, B, C, D);
    }
    // Process each subsegment
    for (size_t i = 0; i < boundaries.size() - 1; ++i) {
        double t_start = boundaries[i];
        double t_end = boundaries[i + 1];
        if (t_end - t_start < EPS) {
            continue; // Skip negligible intervals
        }

        // Determine the closest point on CD for this subsegment
        // Use midpoint to determine whether closest point is along CD or at an endpoint
        //double t_mid = (t_start + t_end) / 2.0;
        double t_mid = t_start+(t_end-t_start)/5.0;
        vec P_mid = A + dirAB * t_mid;

        // Compute s for the projection of P_mid onto CD
        double s = (P_mid - C).dot(dirCD) / lenCD_sq;
        vec closest_point;
        bool closest_is_projection = true;

        if (s <= 0.0) {
            closest_point = C;
            closest_is_projection = false;
        } else if (s >= 1.0) {
            closest_point = D;
            closest_is_projection = false;
        } else {
            closest_point = C + dirCD * s;
        }
        if (closest_is_projection) {
            // Closest point lies along CD
            // Set up the quadratic equation for D(t)^2 = d^2
            vec E = A - C;
            double A_coef = lenCD_sq * dirAB.dot(dirAB) - pow(dirAB.dot(dirCD), 2);
            double B_coef = 2 * (lenCD_sq * dirAB.dot(E) - dirAB.dot(dirCD) * dirCD.dot(E));
            double C_coef = lenCD_sq * E.dot(E) - pow(dirCD.dot(E), 2) - lenCD_sq * d * d;

            // Solve quadratic equation A_coef * t^2 + B_coef * t + C_coef = 0
            std::vector<double> roots = solveQuadratic(A_coef, B_coef, C_coef);

            // Filter roots within [t_start, t_end]
            std::vector<double> valid_roots;
            for (double t_root : roots) {
                if (t_root >= t_start - EPS && t_root <= t_end + EPS) {
                    valid_roots.push_back(t_root);
                }
            }

            // Collect interval points
            std::vector<double> interval_points = {t_start};
            interval_points.insert(interval_points.end(), valid_roots.begin(), valid_roots.end());
            interval_points.push_back(t_end);

            // Sort and remove duplicates
            std::sort(interval_points.begin(), interval_points.end());
            interval_points.erase(std::unique(interval_points.begin(), interval_points.end(), [](double a, double b) {
                return std::abs(a - b) < EPS;
            }), interval_points.end());

            // Determine intervals where D(t) < d
            for (size_t j = 0; j < interval_points.size() - 1; ++j) {
                double t0 = interval_points[j];
                double t1 = interval_points[j + 1];

                if (t1 - t0 < EPS) {
                    continue;
                }

                // Test midpoint
                double t_test = (t0 + t1) / 2.0;
                vec P_test = A + dirAB * t_test;

                // Compute s for P_test
                double s_test = (P_test - C).dot(dirCD) / lenCD_sq;
                s_test = std::max(0.0, std::min(1.0, s_test));
                vec closest_test = C + dirCD * s_test;
                double dist_sq = (P_test - closest_test).norm_squared();

                if (dist_sq < d * d + EPS) {
                    // Interval where distance < d
                    vec start_point = A + dirAB * t0;
                    vec end_point = A + dirAB * t1;
                    if (t0<interval_start){
                        interval_start = t0;
                    }
                    if (t1>interval_end){
                        interval_end = t1;
                    }
                    intervals.push_back({start_point, end_point});
                }
            }
        } else {
            // Closest point is at C or D
            vec Q = closest_point;
            vec M = A - Q;
            vec N = dirAB;

            double a_quad = N.dot(N);
            double b_quad = 2 * M.dot(N);
            double c_quad = M.dot(M) - d * d;

            // Solve quadratic equation a_quad * t^2 + b_quad * t + c_quad = 0
            std::vector<double> roots = solveQuadratic(a_quad, b_quad, c_quad);

            // Filter roots within [t_start, t_end]
            std::vector<double> valid_roots;
            for (double t_root : roots) {
                if (t_root >= t_start - EPS && t_root <= t_end + EPS) {
                    valid_roots.push_back(t_root);
                }
            }

            // Collect interval points
            std::vector<double> interval_points = {t_start};
            interval_points.insert(interval_points.end(), valid_roots.begin(), valid_roots.end());
            interval_points.push_back(t_end);

            // Sort and remove duplicates
            std::sort(interval_points.begin(), interval_points.end());
            interval_points.erase(std::unique(interval_points.begin(), interval_points.end(), [](double a, double b) {
                return std::abs(a - b) < EPS;
            }), interval_points.end());

            // Determine intervals where D(t) < d
            for (size_t j = 0; j < interval_points.size() - 1; ++j) {
                double t0 = interval_points[j];
                double t1 = interval_points[j + 1];
                if (t1 - t0 < EPS) {
                    continue;
                }

                // Test midpoint
                double t_test = (t0 + t1) / 2.0;
                vec P_test = A + dirAB * t_test;
                double dist_sq = (P_test - Q).norm_squared();

                if (dist_sq < d * d + EPS) {
                    // Interval where distance < d
                    vec start_point = A + dirAB * t0;
                    vec end_point = A + dirAB * t1;
                    intervals.push_back({start_point, end_point});
                    if (t0<interval_start){
                        interval_start = t0;
                    }
                    if (t1>interval_end){
                        interval_end = t1;
                    }
                }
            }
        }
    }
    if (interval_start<=interval_end){
        vec inter = A + dirAB * interval_start;
        return {interval_end-interval_start, A + dirAB * interval_start, A+dirAB*interval_end};
    }
    else{
        return {0.0,A,A};
    }
}

// Function to calculate the points on segment 1 that are distance d away from segment 2
am_interaction analyze_am(vec actin_left, vec actin_right, 
                        vec myosin_left, vec myosin_right, double d, std::vector<double> box) {
    apply_pbc(actin_left, actin_right, myosin_left, myosin_right, box);
    vec actin_mid = 0.5 * (actin_left + actin_right);
    vec myosin_mid = 0.5 * (myosin_left + myosin_right);
    std::tuple<double, vec, vec> result = subsegment_within_distance(actin_left, actin_right, myosin_left, myosin_right, d);
    am_interaction interaction;
    interaction.myosin_binding_ratio = std::get<0>(result);
    interaction.myosin_binding_start = std::get<1>(result);
    interaction.myosin_binding_end = std::get<2>(result);
    // if (interaction.myosin_binding_ratio > 0) {
    // printf("actin left (%f, %f, %f)\n"
    //         "actin right (%f, %f, %f)\n"
    //         "myosin left (%f, %f, %f)\n"
    //         "myosin right (%f, %f, %f)\n"
    //         "myosin binding start (%f, %f, %f)\n"
    //         "myosin binding end (%f, %f, %f)\n",
    //         actin_left.x, actin_left.y, actin_left.z,
    //         actin_right.x, actin_right.y, actin_right.z,
    //         myosin_left.x, myosin_left.y, myosin_left.z,
    //         myosin_right.x, myosin_right.y, myosin_right.z,
    //         interaction.myosin_binding_start.x, interaction.myosin_binding_start.y, interaction.myosin_binding_start.z,
    //         interaction.myosin_binding_end.x, interaction.myosin_binding_end.y, interaction.myosin_binding_end.z);
    // }

    //crosslinkable
    interaction.crosslinkable_start = actin_left;
    interaction.crosslinkable_end = interaction.myosin_binding_start;
    interaction.crosslinkable_ratio = (interaction.crosslinkable_end-actin_left).norm()/actin_right.distance(actin_left,box);
    double dot = (actin_mid-actin_left).dot(myosin_mid-myosin_left);
    vec partial_start, partial_end;
    if (dot>0){
        partial_start = myosin_left;
        partial_end = myosin_left + (myosin_right-myosin_left)/3.0;
    }
    else{
        partial_start = myosin_right - (myosin_right-myosin_left)/3.0;
        partial_end = myosin_right;
    }
    // printf("partial_start: %f %f\n",partial_start.x,partial_start.y);
    // printf("partial_end: %f %f\n",partial_end.x,partial_end.y);
    std::tuple<double, vec, vec> result2 = subsegment_within_distance(actin_left, actin_right, partial_start, partial_end, d);
    interaction.partial_binding_ratio = std::get<0>(result2);
    return interaction;   
}

} // namespace geometry








