#ifndef ACTN2_LOAD_KINETICS_H
#define ACTN2_LOAD_KINETICS_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>

namespace actn2_load_kinetics {

// Legacy reproduces the original strong, distinct-myosin proxy. ForceProjection
// classifies an existing ACTN2 crosslink from the signed force that the
// crosslink transmits to both actins.
enum class StatusMode {
    legacy_motor_identity = 0,
    force_projection = 1,
    motor_only = 2,
    transmitted_only = 3,
    opposed_and_transmitted = 4,
    productive_opposed = 5,
    arm_load_productive = 6,
    productive_register = 7,
    motor_transmitted_register = 8,
};

inline StatusMode parse_status_mode(const std::string& value) {
    if (value == "legacy_motor_identity" || value == "legacy") {
        return StatusMode::legacy_motor_identity;
    }
    if (value == "aa_force_projection" || value == "force_projection") {
        return StatusMode::force_projection;
    }
    if (value == "motor_only") return StatusMode::motor_only;
    if (value == "transmitted_only") return StatusMode::transmitted_only;
    if (value == "opposed_and_transmitted" || value == "explicit_crossbridge") {
        return StatusMode::opposed_and_transmitted;
    }
    if (value == "productive_opposed" || value == "productive_motor_pair") {
        return StatusMode::productive_opposed;
    }
    if (value == "arm_load_productive" ||
        value == "bilateral_arm_load_productive") {
        return StatusMode::arm_load_productive;
    }
    if (value == "productive_register" ||
        value == "minimal_productive_register") {
        return StatusMode::productive_register;
    }
    if (value == "motor_transmitted_register" ||
        value == "identity_transmitted_register") {
        return StatusMode::motor_transmitted_register;
    }
    throw std::invalid_argument(
        "aa_status_mode must be legacy_motor_identity, aa_force_projection, "
        "motor_only, transmitted_only, opposed_and_transmitted, or "
        "productive_opposed, arm_load_productive, productive_register, or "
        "motor_transmitted_register");
}

inline const char* status_mode_name(StatusMode mode) {
    switch (mode) {
        case StatusMode::legacy_motor_identity:
            return "legacy_motor_identity";
        case StatusMode::force_projection:
            return "aa_force_projection";
        case StatusMode::motor_only:
            return "motor_only";
        case StatusMode::transmitted_only:
            return "transmitted_only";
        case StatusMode::opposed_and_transmitted:
            return "opposed_and_transmitted";
        case StatusMode::productive_opposed:
            return "productive_opposed";
        case StatusMode::arm_load_productive:
            return "arm_load_productive";
        case StatusMode::productive_register:
            return "productive_register";
        case StatusMode::motor_transmitted_register:
            return "motor_transmitted_register";
    }
    return "unknown";
}

struct ArmLoads {
    // Signed force exerted by actin on its ACTN2 arm, projected from the
    // barbed (+) end toward the pointed (-) end. Positive is the experimentally
    // favorable direction. The force routine supplies force_on_first_actin;
    // force on the second actin is its negative by translational invariance.
    double first = 0.0;
    double second = 0.0;

    double bilateral() const {
        return std::min(first, second);
    }

    double favorable_bilateral() const {
        return std::max(0.0, bilateral());
    }
};

inline double dot(const std::array<double, 3>& a,
                  const std::array<double, 3>& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline double norm(const std::array<double, 3>& value) {
    return std::sqrt(dot(value, value));
}

// Portion of a source-resolved motor force that directly balances the force
// applied to the same actin by one ACTN2 connector.  A passive ACTN2 spring,
// AA prestress, and steric forces cannot create this quantity by themselves.
inline double opposing_support(
    const std::array<double, 3>& motor_force_on_actin,
    const std::array<double, 3>& connector_force_on_actin) {
    const double connector_norm = norm(connector_force_on_actin);
    if (!(connector_norm > 0.0)) return 0.0;
    return std::max(
        0.0, -dot(motor_force_on_actin, connector_force_on_actin) /
                 connector_norm);
}

// Both ACTN2 arms must carry favorable connector load and both selected motor
// sources must balance that same connector.  The weakest of those four
// requirements is the motor-transmitted evidence for this particular pair.
inline double motor_transmitted_load(const ArmLoads& connector_arm_loads,
                                     double first_motor_support,
                                     double second_motor_support) {
    return std::min({
        std::max(0.0, connector_arm_loads.first),
        std::max(0.0, connector_arm_loads.second),
        std::max(0.0, first_motor_support),
        std::max(0.0, second_motor_support)});
}

inline bool distinct_motor_sources(int first_myosin, int second_myosin) {
    return first_myosin >= 0 && second_myosin >= 0 &&
           first_myosin != second_myosin;
}

inline double update_support_ema(double previous, double instantaneous,
                                 double dt, double averaging_window) {
    previous = std::max(0.0, previous);
    instantaneous = std::max(0.0, instantaneous);
    if (!(averaging_window > 0.0)) return instantaneous;
    const double alpha = std::clamp(
        std::max(0.0, dt) / averaging_window, 0.0, 1.0);
    return previous + alpha * (instantaneous - previous);
}

// A mature ACTN2 conformation is stabilized by recently transmitted motor
// support but is not an irreversible record of a past loading event.  The
// unloaded registered complex relaxes at zero_load_rate; support suppresses
// that relaxation, while registry error beyond the retention window
// accelerates it.  The exponential form is the coarse-grained Bell analogue
// for a reversible conformational transition rather than molecular
// detachment.  Exponent clamping keeps extreme diagnostic inputs finite.
inline double mature_unlock_rate(
    double zero_load_rate, double filtered_motor_support,
    double support_scale, double registry_error,
    double registry_retention, double registry_scale) {
    zero_load_rate = std::max(0.0, zero_load_rate);
    if (!(zero_load_rate > 0.0)) return 0.0;
    filtered_motor_support = std::max(0.0, filtered_motor_support);
    support_scale = std::max(support_scale, 1e-12);
    registry_error = std::max(0.0, registry_error);
    registry_retention = std::max(0.0, registry_retention);
    registry_scale = std::max(registry_scale, 1e-12);
    const double support_exponent = -filtered_motor_support / support_scale;
    const double registry_exponent =
        std::max(0.0, registry_error - registry_retention) / registry_scale;
    return zero_load_rate * std::exp(std::clamp(
        support_exponent + registry_exponent, -50.0, 50.0));
}

inline double mature_unlock_probability(double rate, double dt) {
    if (!(rate > 0.0) || !(dt > 0.0)) return 0.0;
    return -std::expm1(-std::min(rate * dt, 700.0));
}

inline bool replace_motor_pair(double best_alternative_support,
                               double current_pair_support,
                               double switch_ratio) {
    switch_ratio = std::max(1.0, switch_ratio);
    return best_alternative_support >
           switch_ratio * std::max(0.0, current_pair_support);
}

// A source-resolved motor force is a finite budget.  When the same source can
// oppose several ACTN2 connector arms, distribute that budget in proportion
// to their raw projected demands instead of crediting the full force to every
// link independently.
inline double allocated_source_support(double raw_support,
                                       double total_raw_support,
                                       double source_force_budget) {
    raw_support = std::max(0.0, raw_support);
    total_raw_support = std::max(0.0, total_raw_support);
    source_force_budget = std::max(0.0, source_force_budget);
    if (!(raw_support > 0.0) || !(total_raw_support > 0.0) ||
        !(source_force_budget > 0.0)) {
        return 0.0;
    }
    return raw_support * std::min(1.0,
        source_force_budget / total_raw_support);
}

struct MotorPairIdentity {
    int myosin_i = -1;
    int myosin_j = -1;
    std::int8_t half_i = 0;
    std::int8_t half_j = 0;

    bool valid() const {
        return distinct_motor_sources(myosin_i, myosin_j);
    }
};

inline bool same_motor_pair(const MotorPairIdentity& first,
                            const MotorPairIdentity& second) {
    return first.myosin_i == second.myosin_i &&
           first.myosin_j == second.myosin_j &&
           first.half_i == second.half_i &&
           first.half_j == second.half_j;
}

struct MotorPairState {
    MotorPairIdentity identity;
    std::uint64_t since_step = 0;
    std::uint64_t last_seen_step = 0;
    double support_ema = 0.0;
    double support_instantaneous = 0.0;
    MotorPairIdentity challenger_identity;
    std::uint64_t challenger_since_step = 0;
    double challenger_support_ema = 0.0;
    bool locked = false;
};

enum class MotorPairEvent : int {
    none = 0,
    selected = 1,
    switched = 2,
    dropped = 3,
};

struct MotorPairUpdate {
    MotorPairEvent event = MotorPairEvent::none;
    double candidate_age = 0.0;
    bool committed = false;
};

// Update one weak ACTN2 molecule's selected motor pair.  Short gaps in strong
// head support are allowed for head turnover, but an identity that has supplied
// no support for dropout_time is forgotten and must earn a fresh commitment
// interval.  A positive alternative may replace the current identity before
// that timeout when it clears the switch hysteresis.
inline MotorPairUpdate update_motor_pair_candidate(
    MotorPairState& state, const MotorPairIdentity& best_identity,
    double best_support, double current_support, std::uint64_t current_step,
    double dt, double support_window, double commit_time,
    double switch_ratio, double exchange_time, double dropout_time) {
    MotorPairUpdate result;
    best_support = std::max(0.0, best_support);
    current_support = std::max(0.0, current_support);
    if (state.locked) {
        // The motor identity is locked, not its mechanical evidence.  Continue
        // filtering the support carried by that exact pair so a reversible
        // mature ACTN2 state can bridge brief head turnover yet eventually
        // unlock after sustained unloading.
        state.support_instantaneous = current_support;
        state.support_ema = update_support_ema(
            state.support_ema, current_support, dt, support_window);
        if (current_support > 0.0) state.last_seen_step = current_step;
        result.candidate_age = current_step >= state.since_step
            ? (current_step - state.since_step) * std::max(0.0, dt) : 0.0;
        result.committed = state.identity.valid();
        return result;
    }

    const bool had_candidate = state.identity.valid();
    const bool best_is_current = had_candidate &&
        same_motor_pair(best_identity, state.identity);

    // Do not let ordinary head-turnover gaps turn the ratio test into
    // "anything beats zero".  During the dropout grace interval the incumbent
    // identity remains the candidate and a transient alternative cannot steal
    // it.  If the gap reaches dropout_time, forget the incumbent first; a new
    // identity may be selected on the following update and must earn a fresh
    // commitment interval.
    if (had_candidate && !(current_support > 0.0)) {
        state.support_instantaneous = 0.0;
        state.challenger_identity = MotorPairIdentity{};
        state.challenger_since_step = 0;
        state.challenger_support_ema = 0.0;
        const double unseen_time = current_step >= state.last_seen_step
            ? (current_step - state.last_seen_step) * std::max(0.0, dt)
            : 0.0;
        if (unseen_time + 1e-15 >= std::max(0.0, dropout_time)) {
            state = MotorPairState{};
            result.event = MotorPairEvent::dropped;
            return result;
        }
        state.support_ema = update_support_ema(
            state.support_ema, 0.0, dt, support_window);
        result.candidate_age = current_step >= state.since_step
            ? (current_step - state.since_step) * std::max(0.0, dt) : 0.0;
        result.committed = result.candidate_age + 1e-15 >=
            std::max(0.0, commit_time);
        return result;
    }

    const bool qualifying_alternative = had_candidate &&
        best_identity.valid() && best_support > 0.0 && !best_is_current &&
        replace_motor_pair(best_support, current_support, switch_ratio);
    bool replace = !had_candidate && best_identity.valid() && best_support > 0.0;
    if (qualifying_alternative) {
        if (!same_motor_pair(state.challenger_identity, best_identity)) {
            state.challenger_identity = best_identity;
            state.challenger_since_step = current_step;
            state.challenger_support_ema = 0.0;
        }
        state.challenger_support_ema = update_support_ema(
            state.challenger_support_ema, best_support, dt, support_window);
        const double challenger_age = current_step >= state.challenger_since_step
            ? (current_step - state.challenger_since_step) *
                  std::max(0.0, dt)
            : 0.0;
        replace = challenger_age + 1e-15 >= std::max(0.0, exchange_time);
    } else {
        state.challenger_identity = MotorPairIdentity{};
        state.challenger_since_step = 0;
        state.challenger_support_ema = 0.0;
    }
    if (replace) {
        result.event = had_candidate ? MotorPairEvent::switched
                                     : MotorPairEvent::selected;
        state.identity = best_identity;
        state.since_step = current_step;
        state.last_seen_step = current_step;
        state.support_ema = 0.0;
        state.challenger_identity = MotorPairIdentity{};
        state.challenger_since_step = 0;
        state.challenger_support_ema = 0.0;
        current_support = best_support;
    }

    if (state.identity.valid()) {
        state.support_instantaneous = current_support;
        if (current_support > 0.0) {
            state.last_seen_step = current_step;
        }
        state.support_ema = update_support_ema(
            state.support_ema, current_support, dt, support_window);
    } else {
        state.support_ema = 0.0;
        state.support_instantaneous = 0.0;
    }

    result.candidate_age = state.identity.valid() &&
            current_step >= state.since_step
        ? (current_step - state.since_step) * std::max(0.0, dt) : 0.0;
    result.committed = state.identity.valid() &&
        result.candidate_age + 1e-15 >= std::max(0.0, commit_time);
    return result;
}

inline double opposed_motor_drive(const std::array<double, 3>& first,
                                  const std::array<double, 3>& second) {
    return 0.5 * norm({first[0] - second[0],
                       first[1] - second[1],
                       first[2] - second[2]});
}

// Motor-only cannot reject co-translation; transmitted-only can be contaminated
// by passive prestress. The combined rule requires both opposed active drive
// and force actually carried by the connector.
inline double catch_load(StatusMode mode,
                         const std::array<double, 3>& first_motor_force,
                         const std::array<double, 3>& second_motor_force,
                         double transmitted_force) {
    transmitted_force = std::max(0.0, transmitted_force);
    const double motor_only_load =
        std::min(norm(first_motor_force), norm(second_motor_force));
    if (mode == StatusMode::motor_only) return motor_only_load;
    if (mode == StatusMode::transmitted_only) return transmitted_force;
    if (mode == StatusMode::opposed_and_transmitted) {
        return std::min(opposed_motor_drive(first_motor_force, second_motor_force),
                        transmitted_force);
    }
    return 0.0;
}

inline int classify_existing_bond_scalar(int previous_status,
                                         double load,
                                         double force_on,
                                         double force_off) {
    force_on = std::max(0.0, force_on);
    force_off = std::clamp(force_off, 0.0, force_on);
    const double threshold = previous_status == 2 ? force_off : force_on;
    return load > threshold ? 2 : 1;
}

inline ArmLoads signed_arm_loads(
    const std::array<double, 3>& force_on_first_actin,
    const std::array<double, 3>& first_barbed_to_pointed,
    const std::array<double, 3>& second_barbed_to_pointed) {
    // ACTN2 experiences the reaction to the force it applies to each actin.
    // F_on_second_actin = -F_on_first_actin.
    return {
        -dot(force_on_first_actin, first_barbed_to_pointed),
         dot(force_on_first_actin, second_barbed_to_pointed),
    };
}

inline ArmLoads limit_by_available_load(const ArmLoads& raw,
                                        double first_available,
                                        double second_available) {
    return {
        std::min(std::max(0.0, raw.first), std::max(0.0, first_available)),
        std::min(std::max(0.0, raw.second), std::max(0.0, second_available)),
    };
}

// Allocate the productive myosin load on one actin among all ACTN2 arms that
// are currently capable of transmitting favorable load.  Raw spring force is
// used only as a load-sharing weight and a transmission ceiling; it is not an
// independent source of catch stabilization.
inline double allocated_motor_load(double raw_favorable_arm_load,
                                   double total_raw_favorable_arm_load,
                                   double available_productive_motor_load) {
    const double raw = std::max(0.0, raw_favorable_arm_load);
    const double total = std::max(0.0, total_raw_favorable_arm_load);
    const double available = std::max(0.0, available_productive_motor_load);
    if (!(raw > 0.0) || !(total > 0.0) || !(available > 0.0)) return 0.0;
    return available * raw / total;
}

inline double motor_supported_arm_load(double raw_favorable_arm_load,
                                       double total_raw_favorable_arm_load,
                                       double available_productive_motor_load) {
    const double raw = std::max(0.0, raw_favorable_arm_load);
    return std::min(
        raw,
        allocated_motor_load(raw, total_raw_favorable_arm_load,
                             available_productive_motor_load));
}

// Convert a legacy pair-level stiffness into a fixed stiffness for each
// explicit ACTN2 molecule. The divisor is a reference coarse-graining choice,
// not the instantaneous occupancy: surviving links must not change stiffness
// when a neighbor binds or detaches.
inline double fixed_single_link_stiffness(double legacy_pair_stiffness,
                                          int reference_links) {
    return std::max(0.0, legacy_pair_stiffness) /
           static_cast<double>(std::max(1, reference_links));
}

inline int classify_existing_bond(int previous_status,
                                  const ArmLoads& loads,
                                  double force_on,
                                  double force_off) {
    force_on = std::max(0.0, force_on);
    force_off = std::clamp(force_off, 0.0, force_on);
    const double threshold = previous_status == 2 ? force_off : force_on;
    return (loads.first > threshold && loads.second > threshold) ? 2 : 1;
}

// A molecular ACTN2 catch state requires load at both directional
// ACTN2--actin interfaces. Callers first allocate the productive motor budget
// among transmitting arms and bound it by each arm's raw mechanical capacity;
// productive_opposed_support remains the pair-level direction/source gate.
inline int classify_arm_loaded_productive(
    int previous_status, const ArmLoads& loads,
    double productive_opposed_support, double force_on, double force_off,
    double support_threshold = 0.0) {
    if (!(productive_opposed_support > std::max(0.0, support_threshold))) {
        return 1;
    }
    return classify_existing_bond(
        previous_status, loads, force_on, force_off);
}

inline bool registry_gate(int previous_status, double direction_dot,
                          double registry_error, double delta_on,
                          double delta_off) {
    if (!(direction_dot < 0.0) || !std::isfinite(registry_error)) return false;
    delta_on = std::max(0.0, delta_on);
    delta_off = std::max(delta_on, delta_off);
    const double threshold = previous_status == 2 ? delta_off : delta_on;
    return registry_error <= threshold;
}

// Productive opposed motor load is required only to create status 2. Once
// created, the state persists until stochastic detachment or loss of the
// wider geometric retention gate.
inline int classify_productive_register(
    int previous_status, double productive_opposed_support, double force_on,
    double direction_dot, double registry_error, double delta_on,
    double delta_off) {
    if (previous_status == 2) {
        return registry_gate(
            2, direction_dot, registry_error, delta_on, delta_off) ? 2 : 1;
    }
    const bool registered = registry_gate(
        1, direction_dot, registry_error, delta_on, delta_off);
    return registered &&
           productive_opposed_support > std::max(0.0, force_on) ? 2 : 1;
}

// Identity-resolved variant. Formation additionally requires that one motor
// pair has persisted for the requested commitment time. Retention remains a
// molecular/geometric state, so transient head turnover cannot unlock it.
inline int classify_motor_transmitted_register(
    int previous_status, double selected_pair_support, bool pair_committed,
    double force_on, double direction_dot, double registry_error,
    double delta_on, double delta_off) {
    if (previous_status == 2) {
        return registry_gate(
            2, direction_dot, registry_error, delta_on, delta_off) ? 2 : 1;
    }
    return pair_committed &&
           registry_gate(
               1, direction_dot, registry_error, delta_on, delta_off) &&
           selected_pair_support > std::max(0.0, force_on) ? 2 : 1;
}

// Diagnostics for overlap-based motor capacity. The overlap calculation itself
// remains in geometry::analyze_am; these helpers keep orientation semantics
// independent and unit-testable.
inline double nematic_axis_alignment(double actin_dot_myosin) {
    return std::clamp(std::abs(actin_dot_myosin), 0.0, 1.0);
}

inline int myosin_half(double contact_axial_coordinate) {
    return contact_axial_coordinate >= 0.0 ? 1 : -1;
}

inline double productive_polarity_score(double actin_dot_myosin,
                                        int half_sign) {
    const int half = half_sign >= 0 ? 1 : -1;
    return std::clamp(-static_cast<double>(half) * actin_dot_myosin, -1.0, 1.0);
}

inline int productive_myosin_half(double actin_dot_myosin) {
    return actin_dot_myosin > 0.0 ? -1 : 1;
}

// Status 2 represents an existing AA link whose two antiparallel actins are
// being pulled in opposition by productive strong heads.  The productive
// force vectors have already been filtered by myosin-half polarity.  Requiring
// force on both arms rejects single-sided occupancy and co-translation.
inline double productive_opposed_load(
    const std::array<double, 3>& first_direction,
    const std::array<double, 3>& second_direction,
    const std::array<double, 3>& first_productive_force,
    const std::array<double, 3>& second_productive_force) {
    const double polarity = dot(first_direction, second_direction);
    if (!(polarity < 0.0)) return 0.0;
    const double first_norm = norm(first_productive_force);
    const double second_norm = norm(second_productive_force);
    if (!(first_norm > 0.0) || !(second_norm > 0.0)) return 0.0;
    const double force_opposition = std::clamp(
        -dot(first_productive_force, second_productive_force) /
            (first_norm * second_norm),
        0.0, 1.0);
    return force_opposition * std::min(first_norm, second_norm);
}

}  // namespace actn2_load_kinetics

#endif
