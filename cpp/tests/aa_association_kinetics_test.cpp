#include "aa_association_kinetics.h"

#include <cassert>
#include <cmath>

int main() {
    using namespace aa_association_kinetics;
    assert(std::abs(elastic_cost(0.05, 0.05, 1.0, 3000.0, 50.0)) < 1e-14);
    const double strained = elastic_cost(0.06, 0.05, 1.0, 3000.0, 50.0);
    assert(std::abs(strained - 0.15) < 1e-12);
    assert(gated_rate(5000.0, 233.0, strained) < 1e-10 * 5000.0);
    assert(std::abs(gated_rate(5000.0, 233.0, 0.0) - 5000.0) < 1e-12);
    const double on0 = 5000.0;
    const double off0 = 100.0;
    for (double theta : {0.0, 0.1, 0.25, 0.5, 1.0}) {
        const double on = transition_state_on_rate(on0, 233.0, strained, theta);
        const double off = transition_state_off_rate(off0, 233.0, strained, theta);
        const double expected_ratio =
            (on0 / off0) * std::exp(-233.0 * strained);
        assert(std::abs(on / off - expected_ratio) <
               1e-12 * std::max(1.0, expected_ratio));
    }
    assert(std::abs(transition_state_on_rate(on0, 233.0, strained, 0.0) - on0) < 1e-12);
    assert(std::abs(transition_state_off_rate(off0, 233.0, strained, 1.0) - off0) < 1e-12);
    assert(std::abs(event_probability(100.0, 1e-3) -
                    (1.0 - std::exp(-0.1))) < 1e-14);
    assert(std::abs(overlap_extensive_factor(0.05, 0.1) - 0.5) < 1e-14);
    assert(std::abs(overlap_extensive_factor(0.2, 0.1) - 2.0) < 1e-14);
    assert(overlap_extensive_factor(0.4, 0.1) ==
           2.0 * overlap_extensive_factor(0.2, 0.1));
    assert(overlap_extensive_factor(0.0, 0.1) == 0.0);
    // A 20 nm capture mismatch should remain accessible with the soft
    // encounter stiffness even though the mature 3000 pN/um bridge would be
    // overwhelmingly suppressed. Random orientation is deliberately not
    // penalized when encounter angular stiffness is zero.
    const double encounter = elastic_cost(0.0706, 0.0506, 0.0, 10.0, 0.0);
    const double mature = elastic_cost(0.0706, 0.0506, 0.0, 3000.0, 50.0);
    assert(transition_state_on_rate(5000.0, 233.0, encounter, 0.1) > 0.9 * 5000.0);
    assert(transition_state_on_rate(5000.0, 233.0, mature, 0.1) < 1e-100 * 5000.0);
    return 0;
}
