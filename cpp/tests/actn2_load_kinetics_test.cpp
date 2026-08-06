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
    assert(std::string(status_mode_name(StatusMode::productive_register)) ==
           "productive_register");
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

    return 0;
}
