#include "titin_kinetics.h"

#include <cassert>
#include <cmath>

namespace {

bool near(double first, double second, double tolerance = 1e-12) {
    return std::abs(first - second) <= tolerance;
}

}  // namespace

int main() {
    constexpr double stiffness = 40.0;
    constexpr double slack = 0.50;

    // Titin is tension-only: compression and slack generate no force or
    // elastic energy.
    assert(near(titin_kinetics::extension(0.40, slack), 0.0));
    assert(near(titin_kinetics::force(0.40, slack, stiffness), 0.0));
    assert(near(titin_kinetics::energy(0.40, slack, stiffness), 0.0));
    assert(near(titin_kinetics::force(0.50, slack, stiffness), 0.0));

    assert(near(titin_kinetics::extension(0.60, slack), 0.10));
    assert(near(titin_kinetics::force(0.60, slack, stiffness), 4.0));
    assert(near(titin_kinetics::energy(0.60, slack, stiffness), 0.20));

    constexpr double beta = 0.25;
    constexpr double theta = 0.5;
    const double inside_capture =
        titin_kinetics::distance_limited_association_rate(5000.0, 0.66, 0.80);
    const double outside_capture =
        titin_kinetics::distance_limited_association_rate(5000.0, 0.81, 0.80);
    const double unstrained_off =
        titin_kinetics::dissociation_rate(0.5, beta, theta, 0.0);
    const double strained_off =
        titin_kinetics::dissociation_rate(0.5, beta, theta, 10.0);
    assert(near(inside_capture, 5000.0));
    assert(near(outside_capture, 0.0));
    assert(near(unstrained_off, 0.5));
    assert(strained_off > unstrained_off);

    assert(near(titin_kinetics::transition_probability(0.0, 1e-3), 0.0));
    assert(titin_kinetics::transition_probability(5000.0, 1e-3) > 0.99);
    return 0;
}
