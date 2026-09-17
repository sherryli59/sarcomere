#ifndef CATCH_BOND_KINETICS_H
#define CATCH_BOND_KINETICS_H

#include <algorithm>
#include <limits>

namespace catch_bond_kinetics {

// By default, status 1 is an ordinary short-lived contact and only status 2
// receives load-dependent stabilization. The runtime toggle restores the
// original rule, in which lambda*f_load also stabilizes status 1. The
// directional lifetime bonus remains status-2-only in both modes.
inline double lifetime_for_status(
    int status,
    double base_lifetime,
    double directional_base_lifetime,
    double lifetime_coeff,
    double f_load,
    bool stabilize_status1 = false) {
    if (status <= 0) {
        return 0.0;
    }

    double lifetime = std::max(0.0, base_lifetime);
    if (status == 2) {
        lifetime += std::max(0.0, directional_base_lifetime);
    }
    if (status == 2 || (status == 1 && stabilize_status1)) {
        lifetime += std::max(0.0, lifetime_coeff) * std::max(0.0, f_load);
    }
    return lifetime;
}

inline double break_probability(double dt, double lifetime) {
    if (dt <= 0.0) {
        return 0.0;
    }
    if (lifetime <= 0.0) {
        return 1.0;
    }
    return std::clamp(dt / lifetime, 0.0, 1.0);
}

enum class DetachmentReason {
    none = 0,
    geometric_contour = 1,
    stochastic = 2,
};

struct DetachmentDecision {
    bool geometric = false;
    bool stochastic = false;
    double probability = 0.0;
    // NaN means the geometric short circuit fired and no random number was
    // consumed. This distinction is important for exact trajectory replay.
    double draw = std::numeric_limits<double>::quiet_NaN();

    bool detaches() const { return geometric || stochastic; }

    DetachmentReason reason() const {
        if (geometric) return DetachmentReason::geometric_contour;
        if (stochastic) return DetachmentReason::stochastic;
        return DetachmentReason::none;
    }
};

// Classify ACTN2 detachment without changing the existing RNG semantics:
// connector overextension short-circuits before a stochastic draw is made.
template <typename UniformDraw>
inline DetachmentDecision decide_detachment(
    double connector_length, double rupture_length,
    double dt, double lifetime, UniformDraw&& uniform_draw) {
    DetachmentDecision result;
    result.geometric = rupture_length > 0.0 &&
                       connector_length > rupture_length;
    result.probability = break_probability(dt, lifetime);
    if (result.geometric) return result;

    result.draw = uniform_draw();
    result.stochastic = result.draw < result.probability;
    return result;
}

}  // namespace catch_bond_kinetics

#endif  // CATCH_BOND_KINETICS_H
