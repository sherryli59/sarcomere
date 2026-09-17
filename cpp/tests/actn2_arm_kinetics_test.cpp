#include <cassert>
#include <cmath>

#include "actn2_arm_kinetics.h"

int main() {
    using namespace actn2_arm_kinetics;

    assert(is_bridge(State::weak_bridge));
    assert(is_bridge(State::mature_bridge));
    assert(is_tether(State::tethered_to_i));
    assert(arm_i_bound(State::tethered_to_i));
    assert(!arm_j_bound(State::tethered_to_i));

    assert(std::abs(event_probability(2.0, 0.5) -
                    (1.0 - std::exp(-1.0))) < 1e-14);
    const double unloaded = opening_rate(1.0, 0.0, 0.001, 240.0);
    const double loaded = opening_rate(1.0, 10.0, 0.001, 240.0);
    assert(std::abs(unloaded - 1.0) < 1e-14);
    assert(loaded > unloaded);
    // Negative/compressive load must not alter a tensile Bell pathway.
    assert(std::abs(opening_rate(1.0, -10.0, 0.001, 240.0) - 1.0) < 1e-14);

    double draws[] = {0.0, 0.2};
    int at = 0;
    const auto open_i = decide_bridge_opening(
        3.0, 1.0, 1.0, [&]() { return draws[at++]; });
    assert(open_i.opens_i && !open_i.opens_j);
    assert(open_i.probability > 0.98);

    double draws_j[] = {0.0, 0.9};
    at = 0;
    const auto open_j = decide_bridge_opening(
        3.0, 1.0, 1.0, [&]() { return draws_j[at++]; });
    assert(!open_j.opens_i && open_j.opens_j);

    const double registered = elastic_rebinding_rate(100.0, 0.0, 240.0, 0.5);
    const double strained = elastic_rebinding_rate(100.0, 0.02, 240.0, 0.5);
    assert(std::abs(registered - 100.0) < 1e-14);
    assert(strained < registered);

    // Partial A-M overlap is handled by local steric/contact mechanics.  A
    // deterministic arm opening is reserved for complete loss of exposed
    // ACTN2-binding contour.
    assert(!contour_fully_occluded(1.0));
    assert(!contour_fully_occluded(0.2));
    assert(contour_fully_occluded(0.0));
    return 0;
}
