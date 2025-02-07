#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <utility>
#include <tuple>
#include <string>

#include "utils.h"  // This header must define utils::vec
using vec = utils::vec;

namespace geometry {


    // Function declarations

    // Compute the shortest distance between a point x and a segment ab (with PBC).
    double point_segment_distance(vec x, vec a, vec b, std::vector<double> box);

    // Compute the shortest distance (and the corresponding normal vector) between a point x and a segment ab (with PBC).
    std::pair<double, vec> point_segment_distance_w_normal(vec x, vec a, vec b, std::vector<double> box);

    // Compute the orientation of three points with periodic boundary conditions.
    int orientation(const vec& p, const vec& q, const vec& r, const std::vector<double>& box);

    // Compute the orientation of three points (without PBC).
    int orientation(const vec& p, const vec& q, const vec& r);

    // Compute the shortest distance between two segments ab and cd (with PBC).
    double segment_segment_distance(const vec& a, const vec& b, const vec& c, const vec& d, const std::vector<double>& box);

    // Compute the shortest distance between two segments and return a map containing a normal vector and associated information.
    std::pair<double, std::map<std::string, vec>> segment_segment_distance_w_normal(const vec& a, const vec& b, const vec& c, const vec& d, const std::vector<double>& box);

    // Compute the dot product of two 2D vectors.
    double dotProduct(double x1, double y1, double x2, double y2);

    // Structure to store actin–myosin interaction information.
    struct am_interaction {
        double myosin_binding_ratio;
        vec myosin_binding_start;
        vec myosin_binding_end;
        vec crosslinkable_start;
        vec crosslinkable_end;
        double partial_binding_ratio;
        double crosslinkable_ratio;
    };

    // Apply periodic boundary conditions to align four endpoints.
    void apply_pbc(vec& actin_left, vec& actin_right,
                   vec& myosin_left, vec& myosin_right, std::vector<double> box);

    // Return only the roots within [t_start, t_end] (used in quadratic solving).
    std::vector<double> findRootsInInterval(const std::vector<double>& roots, double t_start, double t_end);

    // Solve a quadratic equation: a*t^2 + b*t + c = 0.
    std::vector<double> solveQuadratic(double a, double b, double c);

    // A simple structure to represent an interval.
    struct Interval {
        double t_start;
        double t_end;
    };

    // Find a subsegment on AB where the distance to CD is less than d.
    std::tuple<double, vec, vec> subsegment_within_distance(vec A, vec B, vec C, vec D, double d);

    // Overloaded version: first applies PBC using the provided box.
    std::tuple<double, vec, vec> subsegment_within_distance(vec A, vec B, vec C, vec D, double d, std::vector<double> box);

    // Analyze actin–myosin interaction geometry and return an interaction struct.
    am_interaction analyze_am(vec actin_left, vec actin_right,
                              vec myosin_left, vec myosin_right, double d, std::vector<double> box);
}

#endif // GEOMETRY_H
