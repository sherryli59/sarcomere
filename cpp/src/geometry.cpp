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

// Compute the shortest distance between two line segments in 3D with PBC.
// Segment 1: S1(s) = A + s*(B-A), s in [0,1].
// Segment 2: S2(t) = C + t*(D-C), t in [0,1].
double segment_segment_distance(const vec& A, const vec& B, const vec& C, const vec& D, const std::vector<double>& box) {
    // Compute direction vectors with PBC adjustments
    vec u = pbc_diff(B, A, box); // Direction for segment AB
    vec v = pbc_diff(D, C, box); // Direction for segment CD
    vec w = pbc_diff(A, C, box); // Difference between starting points

    double a_val = u.dot(u);      // Squared length of u
    double b_val = u.dot(v);
    double c_val = v.dot(v);      // Squared length of v
    double d_val = u.dot(w);
    double e_val = v.dot(w);
    double D_val = a_val * c_val - b_val * b_val; // Denomimator

    double sN, sD = D_val; // sN/sD = s, parameter for segment AB
    double tN, tD = D_val; // tN/tD = t, parameter for segment CD

    // If segments are nearly parallel
    if (D_val < EPS) {
        sN = 0.0;
        sD = 1.0;
        tN = e_val;
        tD = c_val;
    } else {
        sN = (b_val * e_val - c_val * d_val);
        tN = (a_val * e_val - b_val * d_val);
        if (sN < 0.0) {
            sN = 0.0;
            tN = e_val;
            tD = c_val;
        } else if (sN > sD) {
            sN = sD;
            tN = e_val + b_val;
            tD = c_val;
        }
    }

    // Clamp t to [0, 1]
    if (tN < 0.0) {
        tN = 0.0;
        if (-d_val < 0.0)
            sN = 0.0;
        else if (-d_val > a_val)
            sN = sD;
        else {
            sN = -d_val;
            sD = a_val;
        }
    } else if (tN > tD) {
        tN = tD;
        if (-d_val + b_val < 0.0)
            sN = 0.0;
        else if (-d_val + b_val > a_val)
            sN = sD;
        else {
            sN = -d_val + b_val;
            sD = a_val;
        }
    }

    double sc = (std::abs(sN) < EPS ? 0.0 : sN / sD);
    double tc = (std::abs(tN) < EPS ? 0.0 : tN / tD);

    vec dP = w + u * sc - v * tc;
    return dP.norm();
}

std::pair<double, std::map<std::string, vec>> segment_segment_distance_w_normal(
    const vec& A, const vec& B, const vec& C, const vec& D, const std::vector<double>& box)
{
    // Compute directional vectors (using PBC)
    vec u = pbc_diff(B, A, box); // segment AB direction
    vec v = pbc_diff(D, C, box); // segment CD direction
    vec w = pbc_diff(A, C, box); // vector from C to A

    // Compute scalar coefficients
    double a_val = u.dot(u);      // |u|^2
    double b_val = u.dot(v);
    double c_val = v.dot(v);      // |v|^2
    double d_val = u.dot(w);
    double e_val = v.dot(w);
    double D_val = a_val * c_val - b_val * b_val;  // denominator

    double s, t;

    // If segments are nearly parallel, use default values
    if (std::abs(D_val) < EPS) {
        s = 0.0;
        t = (c_val > EPS ? e_val / c_val : 0.0);
    } else {
        s = (b_val * e_val - c_val * d_val) / D_val;
        t = (a_val * e_val - b_val * d_val) / D_val;
    }

    // Clamp s and t to the [0, 1] interval
    s = std::max(0.0, std::min(1.0, s));
    t = std::max(0.0, std::min(1.0, t));

    // Compute the closest points on each segment:
    vec P = A + u * s;  // Closest point on segment AB
    vec Q = C + v * t;  // Closest point on segment CD

    // Compute the normal vector from Q to P, applying PBC
    vec normal = pbc_diff(P, Q, box);

    double distance = normal.norm();

    // Create map to store extra information
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
        double t1 = (-b + sqrt_disc) / (2 * a);
        double t2 = (-b - sqrt_disc) / (2 * a);
        roots.push_back(t1);
        roots.push_back(t2);
    }
    return roots;
}

// Main function to find intervals on AB where the distance to CD is less than d, extended to 3D with PBC.
std::tuple<double, vec, vec> subsegment_within_distance(vec A, vec B, vec C, vec D, double d, std::vector<double> box) {
    std::vector<std::pair<vec, vec>> intervals;
    double interval_start = 1.0;
    double interval_end = 0.0;

    // Compute directional vectors using PBC.
    vec dirAB = pbc_diff(B, A, box); // Direction from A to B
    vec dirCD = pbc_diff(D, C, box); // Direction from C to D
    double lenAB_sq = dirAB.norm_squared();
    double lenCD_sq = dirCD.norm_squared();

    // Compute midpoints using minimal image convention.
    vec AB_mid = A + dirAB * 0.5;
    vec CD_mid = C + dirCD * 0.5;
    if (pbc_diff(AB_mid, CD_mid, box).norm() > dirAB.norm() + dirCD.norm() + d) {
        return std::make_tuple(0.0, A, A);
    }

    // Find transition points where the closest point on CD changes along AB.
    std::vector<double> boundaries = {0.0, 1.0};

    if (lenAB_sq > EPS) {
        vec AC = pbc_diff(A, C, box);
        vec AD = pbc_diff(A, D, box);
        double t_c = - (dirAB.dot(AC)) / lenAB_sq;
        double t_d = - (dirAB.dot(AD)) / lenAB_sq;

        if (t_c > EPS && t_c < 1.0 - EPS)
            boundaries.push_back(t_c);
        if (t_d > EPS && t_d < 1.0 - EPS)
            boundaries.push_back(t_d);
    }

    std::sort(boundaries.begin(), boundaries.end());
    boundaries.erase(std::unique(boundaries.begin(), boundaries.end(), [](double a, double b) {
        return std::abs(a - b) < EPS;
    }), boundaries.end());

    // Process each subsegment of AB defined by the boundaries.
    for (size_t i = 0; i < boundaries.size() - 1; ++i) {
        double t_start = boundaries[i];
        double t_end = boundaries[i + 1];
        if (t_end - t_start < EPS)
            continue;

        double t_mid = (t_start + t_end) / 2.0;
        vec P_mid = A + dirAB * t_mid;

        // Compute the projection parameter s for P_mid onto CD (using PBC for difference)
        double s = (pbc_diff(P_mid, C, box)).dot(dirCD) / lenCD_sq;

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
            // Set up quadratic equation for distance from P(t) to CD:
            vec E = pbc_diff(A, C, box);
            double A_coef = lenCD_sq * dirAB.dot(dirAB) - std::pow(dirAB.dot(dirCD), 2);
            double B_coef = 2 * (lenCD_sq * dirAB.dot(E) - dirAB.dot(dirCD) * dirCD.dot(E));
            double C_coef = lenCD_sq * E.dot(E) - std::pow(dirCD.dot(E), 2) - lenCD_sq * d * d;

            std::vector<double> roots = solveQuadratic(A_coef, B_coef, C_coef);
            std::vector<double> valid_roots;
            for (double t_root : roots) {
                if (t_root >= t_start - EPS && t_root <= t_end + EPS)
                    valid_roots.push_back(t_root);
            }

            std::vector<double> interval_points = {t_start};
            interval_points.insert(interval_points.end(), valid_roots.begin(), valid_roots.end());
            interval_points.push_back(t_end);
            std::sort(interval_points.begin(), interval_points.end());
            interval_points.erase(std::unique(interval_points.begin(), interval_points.end(), [](double a, double b) {
                return std::abs(a - b) < EPS;
            }), interval_points.end());

            // Determine subintervals where the distance is less than d.
            for (size_t j = 0; j < interval_points.size() - 1; ++j) {
                double t0 = interval_points[j];
                double t1 = interval_points[j + 1];
                if (t1 - t0 < EPS)
                    continue;
                double t_test = (t0 + t1) / 2.0;
                vec P_test = A + dirAB * t_test;
                // Compute s for P_test onto CD using PBC
                double s_test = (pbc_diff(P_test, C, box)).dot(dirCD) / lenCD_sq;
                s_test = std::max(0.0, std::min(1.0, s_test));
                vec closest_test = C + dirCD * s_test;
                double dist_sq = pbc_diff(P_test, closest_test, box).norm_squared();
                if (dist_sq < d * d + EPS) {
                    intervals.push_back({A + dirAB * t0, A + dirAB * t1});
                    interval_start = std::min(interval_start, t0);
                    interval_end = std::max(interval_end, t1);
                }
            }
        } else {
            // Closest point is at an endpoint (C or D)
            vec Q = closest_point;
            vec M = pbc_diff(A, Q, box);
            vec N = dirAB;
            double a_quad = N.dot(N);
            double b_quad = 2 * M.dot(N);
            double c_quad = M.dot(M) - d * d;

            std::vector<double> roots = solveQuadratic(a_quad, b_quad, c_quad);
            std::vector<double> valid_roots;
            for (double t_root : roots) {
                if (t_root >= t_start - EPS && t_root <= t_end + EPS)
                    valid_roots.push_back(t_root);
            }
            std::vector<double> interval_points = {t_start};
            interval_points.insert(interval_points.end(), valid_roots.begin(), valid_roots.end());
            interval_points.push_back(t_end);
            std::sort(interval_points.begin(), interval_points.end());
            interval_points.erase(std::unique(interval_points.begin(), interval_points.end(), [](double a, double b) {
                return std::abs(a - b) < EPS;
            }), interval_points.end());

            for (size_t j = 0; j < interval_points.size() - 1; ++j) {
                double t0 = interval_points[j];
                double t1 = interval_points[j + 1];
                if (t1 - t0 < EPS)
                    continue;
                double t_test = (t0 + t1) / 2.0;
                vec P_test = A + dirAB * t_test;
                double dist_sq = pbc_diff(P_test, Q, box).norm_squared();
                if (dist_sq < d * d + EPS) {
                    intervals.push_back({A + dirAB * t0, A + dirAB * t1});
                    interval_start = std::min(interval_start, t0);
                    interval_end = std::max(interval_end, t1);
                }
            }
        }
    }

    if (interval_start <= interval_end) {
        return std::make_tuple(interval_end - interval_start,
                               A + dirAB * interval_start,
                               A + dirAB * interval_end);
    } else {
        return std::make_tuple(0.0, A, A);
    }
}




// Function to calculate the points on segment 1 that are distance d away from segment 2
am_interaction analyze_am(vec actin_left, vec actin_right, 
                        vec myosin_left, vec myosin_right, double d, std::vector<double> box) {
    apply_pbc(actin_left, actin_right, myosin_left, myosin_right, box);
    vec actin_mid = 0.5 * (actin_left + actin_right);
    vec myosin_mid = 0.5 * (myosin_left + myosin_right);
    std::tuple<double, vec, vec> result = subsegment_within_distance(actin_left, actin_right, myosin_left, myosin_right, d, box);
    am_interaction interaction;
    interaction.myosin_binding_ratio = std::get<0>(result);
    interaction.myosin_binding_start = std::get<1>(result);
    interaction.myosin_binding_end = std::get<2>(result);
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
    std::tuple<double, vec, vec> result2 = subsegment_within_distance(actin_left, actin_right, partial_start, partial_end, d, box);
    interaction.partial_binding_ratio = std::get<0>(result2);
    return interaction;   
}

} // namespace geometry








