#pragma once

#include <cmath>

namespace lattice_spacing {

constexpr double physiological_aa_am_ratio = 0.92;

inline double am_preferred_centerline(double myosin_radius,
                                      double actin_radius,
                                      double surface_gap) {
    return myosin_radius + actin_radius + surface_gap;
}

inline double aa_preferred_centerline(double am_centerline,
                                      double aa_am_ratio) {
    return aa_am_ratio * am_centerline;
}

inline double effective_ratio(double aa_centerline,
                              double am_centerline) {
    return am_centerline > 0.0 ? aa_centerline / am_centerline : 0.0;
}

inline double projected_surface_span(double centerline_distance,
                                     double filament_radius) {
    return centerline_distance - 2.0 * filament_radius;
}

}  // namespace lattice_spacing
