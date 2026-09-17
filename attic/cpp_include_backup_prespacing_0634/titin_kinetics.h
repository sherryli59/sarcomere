#ifndef TITIN_KINETICS_H
#define TITIN_KINETICS_H

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace titin_kinetics {

inline double extension(double distance, double slack_length) {
    return std::max(0.0, distance - slack_length);
}

inline double energy(double distance, double slack_length, double stiffness) {
    const double ext = extension(distance, slack_length);
    return 0.5 * stiffness * ext * ext;
}

inline double force(double distance, double slack_length, double stiffness) {
    return stiffness * extension(distance, slack_length);
}

inline double distance_limited_association_rate(double base_rate,
                                                double distance,
                                                double capture_distance) {
    if (base_rate < 0.0 || distance < 0.0 || !(capture_distance > 0.0)) {
        throw std::invalid_argument("invalid titin association kinetics");
    }
    return distance <= capture_distance ? base_rate : 0.0;
}

inline double dissociation_rate(double base_rate, double beta, double theta,
                                double elastic_energy) {
    if (base_rate < 0.0 || beta < 0.0 || theta < 0.0 || theta > 1.0 ||
        elastic_energy < 0.0) {
        throw std::invalid_argument("invalid titin dissociation kinetics");
    }
    const double exponent =
        std::min(50.0, beta * (1.0 - theta) * elastic_energy);
    return base_rate * std::exp(exponent);
}

inline double transition_probability(double rate, double dt) {
    if (!(rate > 0.0) || !(dt > 0.0)) return 0.0;
    return -std::expm1(-rate * dt);
}

}  // namespace titin_kinetics

#endif
