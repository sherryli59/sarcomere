#ifndef ACTN2_ARM_KINETICS_H
#define ACTN2_ARM_KINETICS_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

namespace actn2_arm_kinetics {

// A bridge has both actin-binding domains closed.  A tether retains exactly
// one material attachment; the free domain may rebind locally without the
// ACTN2 molecule first entering the soluble reservoir.
enum class State : int {
    weak_bridge = 1,
    mature_bridge = 2,
    tethered_to_i = 3,
    tethered_to_j = 4,
};

inline bool is_bridge(State state) {
    return state == State::weak_bridge || state == State::mature_bridge;
}

inline bool is_tether(State state) {
    return state == State::tethered_to_i || state == State::tethered_to_j;
}

inline bool arm_i_bound(State state) {
    return state != State::tethered_to_j;
}

inline bool arm_j_bound(State state) {
    return state != State::tethered_to_i;
}

inline State bridge_state_for_status(int status) {
    return status == 2 ? State::mature_bridge : State::weak_bridge;
}

inline double event_probability(double rate, double dt) {
    if (!(rate > 0.0) || !(dt > 0.0) || !std::isfinite(rate) ||
        !std::isfinite(dt)) {
        return 0.0;
    }
    return -std::expm1(-rate * dt);
}

// Bell slip pathway for opening one ACTN2-actin interface.  Only nonnegative
// tensile/contact load is admitted by the caller; compression is not silently
// reinterpreted as stabilizing tension. Lengths are in um, force in pN, and
// beta in (pN um)^-1.
inline double opening_rate(
    double zero_force_rate, double load, double bell_distance,
    double beta, double max_exponent = 50.0) {
    if (!(zero_force_rate > 0.0) || !std::isfinite(zero_force_rate)) return 0.0;
    const double favorable_load =
        std::isfinite(load) ? std::max(0.0, load) : 0.0;
    const double exponent = std::clamp(
        beta > 0.0 && std::isfinite(beta) && bell_distance > 0.0 &&
                std::isfinite(bell_distance)
            ? beta * favorable_load * bell_distance
            : 0.0,
        0.0, std::max(0.0, max_exponent));
    return zero_force_rate * std::exp(exponent);
}

struct OpeningDecision {
    bool opens_i = false;
    bool opens_j = false;
    double rate_i = 0.0;
    double rate_j = 0.0;
    double probability = 0.0;
};

// Competing hazards: at most one arm opens in one timestep.  The first draw
// decides whether any opening occurs, and the second selects the arm in
// proportion to its rate.  This remains well behaved if dt is not infinitesimal.
template <class UniformDraw>
inline OpeningDecision decide_bridge_opening(
    double rate_i, double rate_j, double dt, UniformDraw&& draw) {
    OpeningDecision result;
    result.rate_i = std::max(0.0, rate_i);
    result.rate_j = std::max(0.0, rate_j);
    const double total = result.rate_i + result.rate_j;
    result.probability = event_probability(total, dt);
    if (!(total > 0.0) || draw() >= result.probability) return result;
    const double split = result.rate_i / total;
    if (draw() < split) {
        result.opens_i = true;
    } else {
        result.opens_j = true;
    }
    return result;
}

inline double elastic_rebinding_rate(
    double attempt_rate, double bridge_energy, double beta,
    double transition_state_fraction) {
    if (!(attempt_rate > 0.0) || !std::isfinite(attempt_rate)) return 0.0;
    const double energy =
        std::isfinite(bridge_energy) ? std::max(0.0, bridge_energy) : 0.0;
    const double theta = std::clamp(transition_state_fraction, 0.0, 1.0);
    const double exponent = std::clamp(
        -(beta > 0.0 && std::isfinite(beta) ? beta * theta * energy : 0.0),
        -50.0, 0.0);
    return attempt_rate * std::exp(exponent);
}

inline bool contour_fully_occluded(double crosslinkable_ratio,
                                   double epsilon = 1e-12) {
    return !std::isfinite(crosslinkable_ratio) ||
        crosslinkable_ratio <= std::max(0.0, epsilon);
}

}  // namespace actn2_arm_kinetics

#endif
