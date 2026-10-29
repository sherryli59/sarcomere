#ifndef CATCH_BOND_KINETICS_H
#define CATCH_BOND_KINETICS_H

#include <algorithm>
#include <cmath>
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
    am_overlap_occlusion = 3,
    actn2_myosin_collision = 4,
    mature_connector_exclusion = 5,
};

struct DetachmentDecision {
    bool geometric = false;
    bool stochastic = false;
    bool am_overlap_occluded = false;
    bool actn2_myosin_collision = false;
    bool mature_connector_excluded = false;
    double probability = 0.0;
    // NaN means the geometric short circuit fired and no random number was
    // consumed. This distinction is important for exact trajectory replay.
    double draw = std::numeric_limits<double>::quiet_NaN();

    bool detaches() const {
        return geometric || stochastic || am_overlap_occluded ||
               actn2_myosin_collision || mature_connector_excluded;
    }

    DetachmentReason reason() const {
        if (mature_connector_excluded) {
            return DetachmentReason::mature_connector_exclusion;
        }
        if (am_overlap_occluded) {
            return DetachmentReason::am_overlap_occlusion;
        }
        if (actn2_myosin_collision) {
            return DetachmentReason::actn2_myosin_collision;
        }
        if (geometric) return DetachmentReason::geometric_contour;
        if (stochastic) return DetachmentReason::stochastic;
        return DetachmentReason::none;
    }
};

// Reason 3 is reserved for actual penetration of an actin barbed endpoint
// through the myosin capture capsule.  Mere AM reach/contour overlap is not a
// rupture condition.  The geometric short circuit consumes no RNG draw.
inline DetachmentDecision decide_barbed_tip_capsule_penetration(
    double penetration_i, double penetration_j,
    double rupture_penetration) {
    DetachmentDecision result;
    result.am_overlap_occluded =
        rupture_penetration > 0.0 &&
        (penetration_i >= rupture_penetration ||
         penetration_j >= rupture_penetration);
    result.probability = result.am_overlap_occluded ? 1.0 : 0.0;
    return result;
}

// Burial as a HAZARD RATE rather than a certainty.  The probability-1 form
// above is the only ACTN2 death channel that is a certainty rather than a
// rate, and pairing it with a SOFT repulsion is unwinnable: a spring in a
// thermal bath holds a DISTRIBUTION of penetration depths (measured median
// 8.3 nm against a 10 nm threshold), so any fixed line inside that
// distribution is crossed eventually and the first excursion executes a
// bond whose natural lifetime is orders of magnitude longer.  Measured
// burial episodes last a median 20 us with exit/entry = 1.00; stiffening
// the spring moves the median but never removes the tail.
//
// The hazard is Bell-like in the depth past the threshold, with the
// threshold itself as the scale so no extra length parameter is needed:
//     k(p) = rate0 * exp((p - p_thr) / p_thr)   for p > 0
// k = rate0 at p = p_thr and grows e-fold per further p_thr of burial, so
// a transient dip costs a small probability while sustained deep burial
// still kills promptly.  The two barbed tips are independent channels and
// their hazards add.  rate0 <= 0 is the legacy probability-1 veto.
template <typename UniformDraw>
inline DetachmentDecision decide_barbed_tip_capsule_penetration_rate(
    double penetration_i, double penetration_j,
    double rupture_penetration, double rupture_rate, double dt,
    UniformDraw&& uniform_draw) {
    DetachmentDecision result;
    if (!(rupture_penetration > 0.0) || !(rupture_rate > 0.0) ||
        !(dt > 0.0)) {
        return result;
    }
    auto hazard = [&](double penetration) {
        if (!std::isfinite(penetration) || !(penetration > 0.0)) return 0.0;
        return rupture_rate * std::exp(
            (penetration - rupture_penetration) / rupture_penetration);
    };
    const double rate = hazard(penetration_i) + hazard(penetration_j);
    if (!(rate > 0.0)) return result;
    result.probability = std::clamp(rate * dt, 0.0, 1.0);
    result.draw = uniform_draw();
    result.am_overlap_occluded = result.draw < result.probability;
    return result;
}

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
