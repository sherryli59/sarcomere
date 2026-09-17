#ifndef ROD_MOBILITY_H
#define ROD_MOBILITY_H

#include "utils.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace rod_mobility {

// Unit convention: um, s, pN, pN um.  The exact SI Boltzmann constant is
// 1.380649e-23 J/K and 1 pN um = 1e-18 J.
constexpr double kBoltzmannPnUmPerK = 1.380649e-5;

inline double kbt(double temperature_kelvin) {
    if (!(temperature_kelvin > 0.0) || !std::isfinite(temperature_kelvin)) {
        throw std::invalid_argument("temperature must be finite and positive");
    }
    return kBoltzmannPnUmPerK * temperature_kelvin;
}

inline double beta(double temperature_kelvin) {
    return 1.0 / kbt(temperature_kelvin);
}

struct Diffusion {
    bool enabled = false;
    double length = 0.0;
    double diameter = 0.0;
    double viscosity = 0.0;
    double crowding_parallel = 1.0;
    double crowding_perpendicular = 1.0;
    double crowding_tumble = 1.0;
    double crowding_roll = 1.0;
    double parallel = 0.0;
    double perpendicular = 0.0;
    double tumble = 0.0;
    double roll = 0.0;
};

inline Diffusion finite_cylinder(
    double length,
    double diameter,
    double temperature_kelvin,
    double viscosity_pas,
    double crowding_parallel = 1.0,
    double crowding_perpendicular = 1.0,
    double crowding_tumble = 1.0,
    double crowding_roll = 1.0) {
    if (!(length > 0.0) || !(diameter > 0.0) || !(viscosity_pas > 0.0) ||
        !(crowding_parallel > 0.0) || !(crowding_perpendicular > 0.0) ||
        !(crowding_tumble > 0.0) || !(crowding_roll > 0.0)) {
        throw std::invalid_argument("rod dimensions, viscosity, and crowding factors must be positive");
    }
    const double p = length / diameter;
    if (!(p >= 2.0)) {
        throw std::invalid_argument("finite-cylinder mobility requires aspect ratio >= 2");
    }
    const double inverse_p = 1.0 / p;
    const double inverse_p2 = inverse_p * inverse_p;
    const double nu_parallel = -0.207 + 0.980 * inverse_p - 0.133 * inverse_p2;
    const double nu_perpendicular = 0.839 + 0.185 * inverse_p + 0.233 * inverse_p2;
    const double nu_tumble = -0.662 + 0.917 * inverse_p - 0.050 * inverse_p2;
    const double thermal_energy = kbt(temperature_kelvin);

    Diffusion result;
    result.enabled = true;
    result.length = length;
    result.diameter = diameter;
    result.viscosity = viscosity_pas;
    result.crowding_parallel = crowding_parallel;
    result.crowding_perpendicular = crowding_perpendicular;
    result.crowding_tumble = crowding_tumble;
    result.crowding_roll = crowding_roll;
    result.parallel = thermal_energy * (std::log(p) + nu_parallel) /
                      (2.0 * M_PI * viscosity_pas * length * crowding_parallel);
    result.perpendicular = thermal_energy * (std::log(p) + nu_perpendicular) /
                           (4.0 * M_PI * viscosity_pas * length * crowding_perpendicular);
    result.tumble = 3.0 * thermal_energy * (std::log(p) + nu_tumble) /
                    (M_PI * viscosity_pas * length * length * length * crowding_tumble);
    // Axial rotational drag of a slender cylinder: zeta_roll = pi eta L d^2.
    result.roll = thermal_energy /
                  (M_PI * viscosity_pas * length * diameter * diameter * crowding_roll);
    return result;
}

inline utils::vec project_parallel(const utils::vec& value, const utils::vec& axis) {
    return value.dot(axis) * axis;
}

inline utils::vec project_perpendicular(const utils::vec& value, const utils::vec& axis) {
    return value - project_parallel(value, axis);
}

inline utils::vec apply_diffusion_tensor(
    const utils::vec& value, const utils::vec& axis, const Diffusion& diffusion) {
    return diffusion.perpendicular * project_perpendicular(value, axis) +
           diffusion.parallel * project_parallel(value, axis);
}

inline utils::vec apply_sqrt_diffusion_tensor(
    const utils::vec& gaussian, const utils::vec& axis, const Diffusion& diffusion) {
    return std::sqrt(diffusion.perpendicular) * project_perpendicular(gaussian, axis) +
           std::sqrt(diffusion.parallel) * project_parallel(gaussian, axis);
}

// Exponential map on the unit sphere for a tangent displacement.  This keeps
// the rod axis normalized without treating a direction-space gradient as a
// conventional axial torque.
inline utils::vec sphere_exp(const utils::vec& axis, const utils::vec& tangent) {
    const double magnitude = tangent.norm();
    if (magnitude <= 1e-14) {
        utils::vec result = axis + tangent;
        result.normalize();
        return result;
    }
    utils::vec tangent_direction = tangent / magnitude;
    utils::vec result = std::cos(magnitude) * axis +
                        std::sin(magnitude) * tangent_direction;
    result.normalize();
    return result;
}

// Scale a realized director-space increment (deterministic drift plus thermal
// noise) to the configured angular trust radius. Returning the scale makes it
// possible to apply the identical factor to coupled material-frame roll.
inline double angular_increment_scale(
    const utils::vec& increment, double limit, double slack = 1.0) {
    if (!(limit > 0.0) || !std::isfinite(limit) ||
        !(slack > 0.0) || !std::isfinite(slack)) {
        return 1.0;
    }
    const double magnitude = increment.norm();
    const double allowed = slack * limit;
    if (!(magnitude > allowed) || magnitude <= 1e-14) return 1.0;
    return allowed / magnitude;
}

}  // namespace rod_mobility

#endif
