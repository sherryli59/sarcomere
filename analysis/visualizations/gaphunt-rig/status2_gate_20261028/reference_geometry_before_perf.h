#ifndef ACTN2_REFERENCE_GEOMETRY_H
#define ACTN2_REFERENCE_GEOMETRY_H

#include "utils.h"
#include <algorithm>

namespace actn2_reference {

// Independent of soft-tip force support, taper, and axial onset. Units: um.
constexpr double exclusion_radius = 0.130;
// Retention hysteresis: tolerate 10 nm penetration from either boundary.
constexpr double deletion_buffer = 0.010;
constexpr double deletion_radius = exclusion_radius - deletion_buffer;

// Both predicates must hold at the SAME connector point. Clip the segment
// to the finite axial slab before minimizing its squared radial distance.
// End planes and the 130 nm radial boundary count as excluded.
inline bool intersects_axial_slab(
    const utils::vec& first, const utils::vec& connector,
    const utils::vec& unit_axis, double half_length, double radius) {
    const double z0 = first.dot(unit_axis);
    const double dz = connector.dot(unit_axis);
    double lo = 0.0, hi = 1.0;
    if (dz == 0.0) {
        if (std::abs(z0) > half_length) return false;
    } else {
        double entry = (-half_length - z0) / dz;
        double exit = (half_length - z0) / dz;
        if (entry > exit) std::swap(entry, exit);
        lo = std::max(lo, entry);
        hi = std::min(hi, exit);
        if (lo > hi) return false;
    }
    const utils::vec radial_first = first - z0 * unit_axis;
    const utils::vec radial_delta = connector - dz * unit_axis;
    const double radial_speed_sq = radial_delta.norm_squared();
    const double t = radial_speed_sq > 0.0
        ? std::clamp(-radial_first.dot(radial_delta) / radial_speed_sq, lo, hi)
        : lo;
    return (radial_first + t * radial_delta).norm_squared() <= radius * radius;
}

// Connector endpoints use the same shortest periodic image as ACTN2
// mechanics. Enumerate all myosin images whose bounding sphere can meet
// that segment; a long or tilted cylinder need not use the nearest image.
inline bool connector_overlaps_myosin(
    const utils::vec& first, const utils::vec& second,
    const utils::vec& center, utils::vec axis, double length,
    const std::vector<double>& box, const std::array<bool, 3>& periodic,
    double radius = exclusion_radius, double axial_inset = 0.0) {
    length -= 2.0 * axial_inset;
    if (!(axis.norm_squared() > 0.0) || !(length > 0.0) || !(radius > 0.0))
        return false;
    axis.normalize();
    utils::vec connector = second - first;
    connector.pbc_wrap(box, periodic);
    utils::vec midpoint = first + 0.5 * connector - center;
    midpoint.pbc_wrap(box, periodic);
    const double reach = std::hypot(0.5 * length, radius) + 0.5 * connector.norm();
    const double coord[3]{midpoint.x, midpoint.y, midpoint.z};
    double period[3]{0.0, 0.0, 0.0};
    int lower[3]{0, 0, 0}, upper[3]{0, 0, 0};
    for (int k = 0; k < 3; ++k) {
        if (periodic[k] && k < static_cast<int>(box.size()) && box[k] > 0.0) {
            period[k] = box[k];
            lower[k] = static_cast<int>(std::ceil((-reach - coord[k]) / box[k]));
            upper[k] = static_cast<int>(std::floor((reach - coord[k]) / box[k]));
            if (lower[k] > upper[k]) return false;
        } else if (std::abs(coord[k]) > reach) {
            return false;
        }
    }
    for (int x = lower[0]; x <= upper[0]; ++x)
        for (int y = lower[1]; y <= upper[1]; ++y)
            for (int z = lower[2]; z <= upper[2]; ++z) {
                const utils::vec image_midpoint = midpoint +
                    utils::vec{x * period[0], y * period[1], z * period[2]};
                if (image_midpoint.norm_squared() > reach * reach) continue;
                if (intersects_axial_slab(
                        image_midpoint - 0.5 * connector, connector, axis,
                        0.5 * length, radius)) return true;
            }
    return false;
}

}  // namespace actn2_maturation
#endif
