#include "actn2_load_kinetics.h"

#include <cassert>
#include <cmath>
#include <stdexcept>

namespace {
bool near(double a, double b, double tol = 1e-12) {
    return std::abs(a - b) <= tol;
}
}

int main() {
    using namespace actn2_load_kinetics;

    assert(parse_status_mode("legacy") == StatusMode::legacy_motor_identity);
    assert(parse_status_mode("aa_force_projection") == StatusMode::force_projection);
    assert(parse_status_mode("motor_only") == StatusMode::motor_only);
    assert(parse_status_mode("transmitted_only") == StatusMode::transmitted_only);
    assert(parse_status_mode("explicit_crossbridge") == StatusMode::opposed_and_transmitted);
    assert(parse_status_mode("productive_opposed") == StatusMode::productive_opposed);
    assert(parse_status_mode("arm_load_productive") ==
           StatusMode::arm_load_productive);
    assert(parse_status_mode("productive_register") ==
           StatusMode::productive_register);
    assert(parse_status_mode("motor_transmitted_register") ==
           StatusMode::motor_transmitted_register);
    assert(std::string(status_mode_name(StatusMode::productive_register)) ==
           "productive_register");
    assert(std::string(
               status_mode_name(StatusMode::motor_transmitted_register)) ==
           "motor_transmitted_register");
    bool rejected = false;
    try {
        (void)parse_status_mode("motor_count");
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    assert(rejected);

    // Antiparallel actins pulled away from their shared barbed-end crosslink.
    // The AA force on the first actin points back toward the crosslink.
    const std::array<double, 3> force_on_first{-3.0, 0.0, 0.0};
    const std::array<double, 3> p_first{1.0, 0.0, 0.0};
    const std::array<double, 3> p_second{-1.0, 0.0, 0.0};
    const auto anti = signed_arm_loads(force_on_first, p_first, p_second);
    assert(near(anti.first, 3.0));
    assert(near(anti.second, 3.0));
    assert(near(anti.bilateral(), 3.0));
    assert(classify_existing_bond(1, anti, 2.0, 1.0) == 2);
    const auto no_motor = limit_by_available_load(anti, 0.0, 0.0);
    assert(near(no_motor.bilateral(), 0.0));
    assert(classify_existing_bond(1, no_motor, 0.01, 0.005) == 1);
    const auto motor_limited = limit_by_available_load(anti, 0.4, 0.2);
    assert(near(motor_limited.first, 0.4));
    assert(near(motor_limited.second, 0.2));
    assert(near(motor_limited.bilateral(), 0.2));

    // A 6-pN productive motor budget shared by raw 3- and 1-pN ACTN2 arms
    // allocates 4.5 and 1.5 pN, but actual transmitted load remains bounded
    // by each arm's raw mechanical capacity.  Conversely, when the motor
    // budget is limiting, the allocated loads conserve that budget exactly.
    assert(near(allocated_motor_load(3.0, 4.0, 6.0), 4.5));
    assert(near(allocated_motor_load(1.0, 4.0, 6.0), 1.5));
    assert(near(motor_supported_arm_load(3.0, 4.0, 6.0), 3.0));
    assert(near(motor_supported_arm_load(1.0, 4.0, 6.0), 1.0));
    const double limited_first = motor_supported_arm_load(3.0, 4.0, 2.0);
    const double limited_second = motor_supported_arm_load(1.0, 4.0, 2.0);
    assert(near(limited_first, 1.5));
    assert(near(limited_second, 0.5));
    assert(near(limited_first + limited_second, 2.0));
    assert(near(allocated_motor_load(-1.0, 4.0, 2.0), 0.0));
    assert(near(allocated_motor_load(1.0, 0.0, 2.0), 0.0));
    assert(near(fixed_single_link_stiffness(3000.0, 3), 1000.0));
    assert(near(fixed_single_link_stiffness(3000.0, 1), 3000.0));
    assert(near(fixed_single_link_stiffness(-10.0, 3), 0.0));

    // A common-direction parallel pair cannot favorably load both arms.
    const auto parallel = signed_arm_loads(force_on_first, p_first, p_first);
    assert(near(parallel.first, 3.0));
    assert(near(parallel.second, -3.0));
    assert(near(parallel.favorable_bilateral(), 0.0));
    assert(classify_existing_bond(1, parallel, 2.0, 1.0) == 1);

    // One favorable and one transverse arm is status 1.
    const auto one_arm = signed_arm_loads(
        force_on_first, p_first, std::array<double, 3>{0.0, 1.0, 0.0});
    assert(classify_existing_bond(1, one_arm, 2.0, 1.0) == 1);

    // Hysteresis: an existing status-2 bond remains stabilized between the
    // on/off thresholds, while a status-1 bond cannot newly stabilize there.
    const ArmLoads intermediate{1.5, 1.5};
    assert(classify_existing_bond(1, intermediate, 2.0, 1.0) == 1);
    assert(classify_existing_bond(2, intermediate, 2.0, 1.0) == 2);
    const ArmLoads unloaded{0.5, 2.0};
    assert(classify_existing_bond(2, unloaded, 2.0, 1.0) == 1);

    // Bilateral interface load and productive opposed motor support are both
    // necessary. Motor support gates status but does not replace arm load.
    assert(classify_arm_loaded_productive(1, anti, 0.0, 2.0, 1.0) == 1);
    assert(classify_arm_loaded_productive(1, ArmLoads{0.5, 3.0}, 2.0,
                                          2.0, 1.0) == 1);
    assert(classify_arm_loaded_productive(1, anti, 2.0, 2.0, 1.0) == 2);
    assert(classify_arm_loaded_productive(1, parallel, 2.0, 2.0, 1.0) == 1);
    assert(classify_arm_loaded_productive(2, intermediate, 2.0,
                                          2.0, 1.0) == 2);
    assert(classify_arm_loaded_productive(2, intermediate, 0.0,
                                          2.0, 1.0) == 1);

    assert(near(nematic_axis_alignment(-0.8), 0.8));
    assert(myosin_half(0.0) == 1);
    assert(myosin_half(-0.1) == -1);
    // actin.direction is barbed-to-pointed. Productive attachment therefore
    // requires opposite AM axes on the + myosin half and aligned AM axes on
    // the - half. Test every polarity/half combination explicitly.
    assert(near(productive_polarity_score(-1.0, 1), 1.0));
    assert(near(productive_polarity_score(1.0, 1), -1.0));
    assert(near(productive_polarity_score(1.0, -1), 1.0));
    assert(near(productive_polarity_score(-1.0, -1), -1.0));
    assert(productive_myosin_half(1.0) == -1);
    assert(productive_myosin_half(-1.0) == 1);

    const std::array<double, 3> right{2.0, 0.0, 0.0};
    const std::array<double, 3> left{-2.0, 0.0, 0.0};
    assert(near(opposed_motor_drive(right, right), 0.0));
    assert(near(opposed_motor_drive(right, left), 2.0));
    assert(near(catch_load(StatusMode::motor_only, right, right, 0.0), 2.0));
    assert(near(catch_load(StatusMode::transmitted_only, right, left, 0.7), 0.7));
    assert(near(catch_load(StatusMode::opposed_and_transmitted, right, right, 3.0), 0.0));
    assert(near(catch_load(StatusMode::opposed_and_transmitted, right, left, 0.7), 0.7));

    const std::array<double, 3> zero{0.0, 0.0, 0.0};
    assert(near(productive_opposed_load(right, left, right, left), 2.0));
    assert(near(productive_opposed_load(right, left,
                                        {2.0, 0.0, 0.0},
                                        {-0.5, 0.0, 0.0}), 0.5));
    assert(near(productive_opposed_load(right, right, right, left), 0.0));
    assert(near(productive_opposed_load(right, left, right, zero), 0.0));
    assert(near(productive_opposed_load(right, left, right, right), 0.0));
    assert(classify_existing_bond_scalar(1, 0.7, 1.0, 0.5) == 1);
    assert(classify_existing_bond_scalar(2, 0.7, 1.0, 0.5) == 2);

    // Registry is a formation/retention gate independent of the load gate.
    assert(registry_gate(1, -1.0, 0.10, 0.10, 0.15));
    assert(!registry_gate(1, -1.0, 0.100001, 0.10, 0.15));
    assert(registry_gate(2, -1.0, 0.14, 0.10, 0.15));
    assert(!registry_gate(2, -1.0, 0.16, 0.10, 0.15));
    assert(!registry_gate(1, 1.0, 0.0, 0.10, 0.15));

    // Minimal productive-register rule: motor load is a formation condition,
    // not an instantaneous retention condition. No spring or allocated-arm
    // load enters this classifier.
    assert(classify_productive_register(
               1, 1.01, 1.0, -1.0, 0.05, 0.10, 0.15) == 2);
    assert(classify_productive_register(
               1, 1.00, 1.0, -1.0, 0.05, 0.10, 0.15) == 1);
    assert(classify_productive_register(
               1, 2.00, 1.0, -1.0, 0.11, 0.10, 0.15) == 1);
    assert(classify_productive_register(
               1, 2.00, 1.0, 1.0, 0.00, 0.10, 0.15) == 1);
    // An established catch state survives complete transient motor-load loss.
    assert(classify_productive_register(
               2, 0.00, 1.0, -1.0, 0.14, 0.10, 0.15) == 2);
    assert(classify_productive_register(
               2, 0.00, 1.0, -1.0, 0.151, 0.10, 0.15) == 1);
    assert(classify_productive_register(
               2, 10.0, 1.0, 1.0, 0.00, 0.10, 0.15) == 1);

    // Source-resolved transmission rejects connector prestress without motor
    // opposition, co-directed motor force, and unilateral motor support.
    const std::array<double, 3> connector_i{-3.0, 0.0, 0.0};
    const std::array<double, 3> connector_j{3.0, 0.0, 0.0};
    assert(near(opposing_support({2.0, 0.0, 0.0}, connector_i), 2.0));
    assert(near(opposing_support({-2.0, 0.0, 0.0}, connector_i), 0.0));
    assert(near(opposing_support({-1.5, 0.0, 0.0}, connector_j), 1.5));
    assert(near(motor_transmitted_load(anti, 2.0, 1.5), 1.5));
    assert(near(motor_transmitted_load(anti, 2.0, 0.0), 0.0));
    assert(near(motor_transmitted_load({3.0, -3.0}, 2.0, 2.0), 0.0));
    assert(distinct_motor_sources(2, 3));
    assert(!distinct_motor_sources(2, 2));
    assert(!distinct_motor_sources(-1, 2));

    assert(near(update_support_ema(0.0, 2.0, 0.002, 0.010), 0.4));
    assert(near(update_support_ema(0.4, 2.0, 0.002, 0.010), 0.72));
    assert(near(update_support_ema(0.4, 2.0, 0.002, 0.0), 2.0));
    // Mature-state relaxation is finite when unloaded, suppressed by recent
    // motor support, and accelerated only beyond the registry retention band.
    assert(near(mature_unlock_rate(
        20.0, 0.0, 0.1, 0.02, 0.04, 0.02), 20.0));
    assert(mature_unlock_rate(
        20.0, 0.25, 0.1, 0.02, 0.04, 0.02) < 2.0);
    assert(mature_unlock_rate(
        20.0, 0.0, 0.1, 0.08, 0.04, 0.02) > 100.0);
    assert(near(mature_unlock_rate(
        0.0, 0.0, 0.1, 1.0, 0.04, 0.02), 0.0));
    assert(near(mature_unlock_probability(20.0, 0.001),
                1.0 - std::exp(-0.02)));
    assert(near(mature_unlock_probability(0.0, 0.001), 0.0));
    assert(!replace_motor_pair(1.49, 1.0, 1.5));
    assert(replace_motor_pair(1.51, 1.0, 1.5));

    // One motor source cannot be counted at full strength by every ACTN2
    // molecule it happens to oppose.  Three oversubscribed arms share the
    // finite 4-pN source budget proportionally; under-subscription is left
    // unchanged rather than artificially amplifying the source.
    const double source_total = 2.0 + 3.0 + 5.0;
    const double source_a = allocated_source_support(2.0, source_total, 4.0);
    const double source_b = allocated_source_support(3.0, source_total, 4.0);
    const double source_c = allocated_source_support(5.0, source_total, 4.0);
    assert(near(source_a, 0.8));
    assert(near(source_b, 1.2));
    assert(near(source_c, 2.0));
    assert(near(source_a + source_b + source_c, 4.0));
    assert(near(allocated_source_support(2.0, 3.0, 4.0), 2.0));
    assert(near(allocated_source_support(0.0, 3.0, 4.0), 0.0));

    // Deterministic weak-state annealing: select a pair, reject a subcritical
    // competitor, switch only after it clears hysteresis, then commit after a
    // continuous residence interval.  Locking freezes the identities even if
    // a much stronger alternative later appears.
    MotorPairState pair_state;
    const MotorPairIdentity pair_a{1, 2, 1, -1};
    const MotorPairIdentity pair_b{3, 4, -1, 1};
    auto pair_update = update_motor_pair_candidate(
        pair_state, pair_a, 1.0, 0.0, 0, 0.001, 0.0, 0.010, 1.5, 0.0, 0.010);
    assert(pair_update.event == MotorPairEvent::selected);
    assert(same_motor_pair(pair_state.identity, pair_a));
    assert(!pair_update.committed);
    pair_update = update_motor_pair_candidate(
        pair_state, pair_b, 1.49, 1.0, 5, 0.001, 0.0, 0.010, 1.5, 0.0, 0.010);
    assert(pair_update.event == MotorPairEvent::none);
    assert(same_motor_pair(pair_state.identity, pair_a));
    pair_update = update_motor_pair_candidate(
        pair_state, pair_b, 1.51, 1.0, 6, 0.001, 0.0, 0.010, 1.5, 0.0, 0.010);
    assert(pair_update.event == MotorPairEvent::switched);
    assert(same_motor_pair(pair_state.identity, pair_b));
    assert(!pair_update.committed);
    pair_update = update_motor_pair_candidate(
        pair_state, pair_b, 1.51, 1.51, 16, 0.001, 0.0, 0.010, 1.5, 0.0, 0.010);
    assert(pair_update.committed);
    pair_state.locked = true;
    pair_update = update_motor_pair_candidate(
        pair_state, pair_a, 100.0, 1.0, 17, 0.001, 0.0, 0.010, 1.5, 0.0, 0.010);
    assert(pair_update.event == MotorPairEvent::none);
    assert(same_motor_pair(pair_state.identity, pair_b));
    assert(near(pair_state.support_ema, 1.0));
    pair_update = update_motor_pair_candidate(
        pair_state, MotorPairIdentity{}, 0.0, 0.0, 18, 0.001, 0.010,
        0.010, 1.5, 0.0, 0.010);
    assert(pair_update.event == MotorPairEvent::none);
    assert(same_motor_pair(pair_state.identity, pair_b));
    assert(near(pair_state.support_ema, 0.9));

    // A weak pair may bridge a brief motor-head turnover gap, but the
    // residence clock resets after the configured dropout interval.
    MotorPairState dropout_state;
    pair_update = update_motor_pair_candidate(
        dropout_state, pair_a, 1.0, 0.0, 20, 0.001, 0.0, 0.010, 1.5, 0.0, 0.010);
    assert(pair_update.event == MotorPairEvent::selected);
    pair_update = update_motor_pair_candidate(
        dropout_state, MotorPairIdentity{}, 0.0, 0.0, 29, 0.001, 0.0,
        0.010, 1.5, 0.0, 0.010);
    assert(pair_update.event == MotorPairEvent::none);
    assert(dropout_state.identity.valid());
    // Even a strongly supported challenger cannot exploit an incumbent's
    // brief zero-support turnover gap.
    pair_update = update_motor_pair_candidate(
        dropout_state, pair_b, 100.0, 0.0, 29, 0.001, 0.0,
        0.010, 1.5, 0.0, 0.010);
    assert(pair_update.event == MotorPairEvent::none);
    assert(same_motor_pair(dropout_state.identity, pair_a));
    pair_update = update_motor_pair_candidate(
        dropout_state, MotorPairIdentity{}, 0.0, 0.0, 30, 0.001, 0.0,
        0.010, 1.5, 0.0, 0.010);
    assert(pair_update.event == MotorPairEvent::dropped);
    assert(!dropout_state.identity.valid());
    assert(!pair_update.committed);

    // A mechanically stronger challenger must remain the same best identity
    // continuously for the exchange dwell; momentary dominance is not a
    // partner-exchange event.
    MotorPairState dwell_state;
    pair_update = update_motor_pair_candidate(
        dwell_state, pair_a, 1.0, 0.0, 40, 0.001, 0.0, 0.010,
        1.5, 0.002, 0.010);
    assert(pair_update.event == MotorPairEvent::selected);
    pair_update = update_motor_pair_candidate(
        dwell_state, pair_b, 2.0, 1.0, 41, 0.001, 0.0, 0.010,
        1.5, 0.002, 0.010);
    assert(pair_update.event == MotorPairEvent::none);
    assert(same_motor_pair(dwell_state.identity, pair_a));
    assert(same_motor_pair(dwell_state.challenger_identity, pair_b));
    // Losing best-challenger status resets the dwell clock.
    pair_update = update_motor_pair_candidate(
        dwell_state, pair_a, 1.0, 1.0, 42, 0.001, 0.0, 0.010,
        1.5, 0.002, 0.010);
    assert(!dwell_state.challenger_identity.valid());
    pair_update = update_motor_pair_candidate(
        dwell_state, pair_b, 2.0, 1.0, 43, 0.001, 0.0, 0.010,
        1.5, 0.002, 0.010);
    assert(pair_update.event == MotorPairEvent::none);
    pair_update = update_motor_pair_candidate(
        dwell_state, pair_b, 2.0, 1.0, 44, 0.001, 0.0, 0.010,
        1.5, 0.002, 0.010);
    assert(pair_update.event == MotorPairEvent::none);
    pair_update = update_motor_pair_candidate(
        dwell_state, pair_b, 2.0, 1.0, 45, 0.001, 0.0, 0.010,
        1.5, 0.002, 0.010);
    assert(pair_update.event == MotorPairEvent::switched);
    assert(same_motor_pair(dwell_state.identity, pair_b));
    assert(!dwell_state.challenger_identity.valid());

    assert(classify_motor_transmitted_register(
               1, 2.0, false, 1.0, -1.0, 0.05, 0.10, 0.15) == 1);
    assert(classify_motor_transmitted_register(
               1, 2.0, true, 1.0, -1.0, 0.05, 0.10, 0.15) == 2);
    assert(classify_motor_transmitted_register(
               1, 2.0, true, 1.0, -1.0, 0.11, 0.10, 0.15) == 1);
    // Once mature, loss or replacement of strong heads does not relabel the
    // locked molecular state; only the wider registry gate can demote it.
    assert(classify_motor_transmitted_register(
               2, 0.0, false, 1.0, -1.0, 0.14, 0.10, 0.15) == 2);
    assert(classify_motor_transmitted_register(
               2, 10.0, true, 1.0, -1.0, 0.16, 0.10, 0.15) == 1);

    return 0;
}
