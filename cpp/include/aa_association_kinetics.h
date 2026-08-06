#ifndef AA_ASSOCIATION_KINETICS_H
#define AA_ASSOCIATION_KINETICS_H

#include <algorithm>
#include <cmath>

namespace aa_association_kinetics {

inline double elastic_cost(
    double distance, double preferred_distance, double axis_dot,
    double radial_stiffness, double angular_stiffness) {
    const double extension = distance - preferred_distance;
    const double alignment = std::clamp(axis_dot, -1.0, 1.0);
    return 0.5 * std::max(0.0, radial_stiffness) * extension * extension +
           0.5 * std::max(0.0, angular_stiffness) *
               (1.0 - alignment * alignment);
}

inline double transition_state_on_rate(
    double attempt_rate, double beta, double elastic_energy, double theta) {
    if (!(attempt_rate > 0.0)) return 0.0;
    if (!(beta > 0.0) || !(elastic_energy > 0.0)) return attempt_rate;
    const double split = std::clamp(theta, 0.0, 1.0);
    return attempt_rate *
           std::exp(-std::min(700.0, split * beta * elastic_energy));
}

inline double transition_state_off_rate(
    double reference_rate, double beta, double elastic_energy, double theta) {
    if (!(reference_rate > 0.0)) return 0.0;
    if (!(beta > 0.0) || !(elastic_energy > 0.0)) return reference_rate;
    const double complement = 1.0 - std::clamp(theta, 0.0, 1.0);
    return reference_rate *
           std::exp(std::min(700.0, complement * beta * elastic_energy));
}

// Backward-compatible name for the original all-on-barrier implementation.
inline double gated_rate(double attempt_rate, double beta, double elastic_energy) {
    return transition_state_on_rate(attempt_rate, beta, elastic_energy, 1.0);
}

inline double event_probability(double rate, double dt) {
    if (!(rate > 0.0) || !(dt > 0.0)) return 0.0;
    return -std::expm1(-rate * dt);
}

inline double overlap_extensive_factor(double accessible_overlap,
                                       double reference_length) {
    if (!(accessible_overlap > 0.0) || !(reference_length > 0.0)) return 0.0;
    return accessible_overlap / reference_length;
}

}  // namespace aa_association_kinetics

#endif
