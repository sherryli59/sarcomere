#include "catch_bond_kinetics.h"

#include <cassert>
#include <cmath>

namespace {

bool near(double lhs, double rhs, double tolerance = 1e-12) {
    return std::abs(lhs - rhs) <= tolerance;
}

}  // namespace

int main() {
    using catch_bond_kinetics::break_probability;
    using catch_bond_kinetics::decide_barbed_tip_capsule_penetration;
    using catch_bond_kinetics::decide_detachment;
    using catch_bond_kinetics::DetachmentReason;
    using catch_bond_kinetics::lifetime_for_status;

    // Default: status 1 gets only the short baseline and is load-independent.
    assert(near(lifetime_for_status(1, 1e-5, 2e-3, 0.1, 0.9), 1e-5));
    assert(near(lifetime_for_status(1, 0.0, 0.0, 1.0, 1.0), 0.0));

    // Original mode: the same lambda*f_load term stabilizes status 1 and 2.
    assert(near(lifetime_for_status(1, 1e-5, 2e-3, 0.1, 0.9, true), 0.09001));
    assert(near(lifetime_for_status(1, 0.0, 0.0, 1.0, 1.0, true), 1.0));

    // Status 2 gets the baseline, directional bonus, and catch contribution.
    assert(near(lifetime_for_status(2, 1e-5, 2e-3, 0.1, 0.9), 0.09201));
    assert(near(lifetime_for_status(2, 1e-5, 2e-3, 0.1, 0.9, true), 0.09201));

    // A zero-lifetime status-1 contact is removed on its next KMC decision.
    assert(near(break_probability(1e-5, 0.0), 1.0));
    assert(near(break_probability(1e-5, 1e-5), 1.0));
    assert(near(break_probability(1e-5, 1e-3), 0.01));

    // A contour rupture must not consume a random number. Diagnostics can
    // therefore be added without perturbing the simulated trajectory.
    int draws = 0;
    const auto geometric = decide_detachment(
        0.121, 0.12, 1e-5, 1e-3, [&]() {
            ++draws;
            return 0.0;
        });
    assert(geometric.detaches());
    assert(geometric.geometric);
    assert(!geometric.stochastic);
    assert(geometric.reason() == DetachmentReason::geometric_contour);
    assert(std::isnan(geometric.draw));
    assert(draws == 0);

    const auto stochastic = decide_detachment(
        0.119, 0.12, 1e-5, 1e-3, [&]() {
            ++draws;
            return 0.005;
        });
    assert(stochastic.detaches());
    assert(!stochastic.geometric);
    assert(stochastic.stochastic);
    assert(stochastic.reason() == DetachmentReason::stochastic);
    assert(near(stochastic.probability, 0.01));
    assert(near(stochastic.draw, 0.005));
    assert(draws == 1);

    const auto survives = decide_detachment(
        0.119, 0.12, 1e-5, 1e-3, [&]() {
            ++draws;
            return 0.5;
        });
    assert(!survives.detaches());
    assert(survives.reason() == DetachmentReason::none);
    assert(draws == 2);

    // Reason 3 fires only for true barbed-tip penetration, not mere overlap.
    const auto penetrated =
        decide_barbed_tip_capsule_penetration(0.010, 0.0, 0.010);
    assert(penetrated.detaches());
    assert(penetrated.am_overlap_occluded);
    assert(!penetrated.geometric);
    assert(!penetrated.stochastic);
    assert(penetrated.reason() ==
           DetachmentReason::am_overlap_occlusion);
    assert(near(penetrated.probability, 1.0));
    assert(std::isnan(penetrated.draw));
    catch_bond_kinetics::DetachmentDecision connector;
    connector.mature_connector_excluded = true;
    connector.probability = 1.0;
    assert(connector.detaches());
    assert(connector.reason() == DetachmentReason::mature_connector_exclusion);
    assert(std::isnan(connector.draw));
    assert(!connector.am_overlap_occluded);

    const auto outside =
        decide_barbed_tip_capsule_penetration(0.0, 0.009, 0.010);
    assert(!outside.detaches());
    assert(outside.reason() == DetachmentReason::none);

    // A disabled threshold cannot accidentally recreate the old immediate
    // contour-overlap veto.
    const auto disabled =
        decide_barbed_tip_capsule_penetration(1.0, 1.0, 0.0);
    assert(!disabled.detaches());

    return 0;
}
