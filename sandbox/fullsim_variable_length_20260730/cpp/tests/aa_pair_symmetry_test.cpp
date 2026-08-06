#include "components.h"
#include "interaction.h"

#include <cassert>
#include <cmath>
#include <gsl/gsl_rng.h>

namespace {
bool near(double lhs, double rhs, double tolerance = 1e-8) {
    return std::abs(lhs - rhs) <= tolerance;
}
}

int main() {
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    const std::vector<double> box{4.0, 4.0, 4.0};
    Filament actin(2, 1.0, box, rng);
    actin.set_periodic_axes({false, false, false});
    actin.center[0] = {0.03, 0.0, 0.0};
    actin.center[1] = {-0.01, 0.075, 0.008};
    actin.direction[0] = {1.0, 0.0, 0.0};
    actin.direction[1] = {1.0, 0.0, 0.0};
    actin.update_endpoints();
    int first = 0;
    int second = 1;
    const auto forward = compute_aa_force_and_energy(
        actin, first, second, box, 10.0, 100.0, 0.2, 0.065, 0.1);
    const auto reverse = compute_aa_force_and_energy(
        actin, second, first, box, 10.0, 100.0, 0.2, 0.065, 0.1);

    // Swapping an unordered AA pair reverses center force and exchanges the
    // two orientation drives, without changing the physical interaction.
    for (int component = 0; component < 3; ++component) {
        assert(near(forward[component], -reverse[component]));
        assert(near(forward[3 + component], reverse[6 + component]));
        assert(near(forward[6 + component], reverse[3 + component]));
    }

    // At the radial and axial rest geometry, the complete AA gradient vanishes.
    actin.center[0] = {0.0, 0.0, 0.0};
    actin.center[1] = {0.0, 0.065, 0.0};
    actin.update_endpoints();
    const auto equilibrium = compute_aa_force_and_energy(
        actin, first, second, box, 10.0, 100.0, 0.2, 0.065, 0.1);
    for (double value : equilibrium) assert(near(value, 0.0, 1e-7));

    // The working coarse-grained Z-disc geometry has a finite antiparallel
    // barbed-end overlap. Unit rods with W=0.1 therefore have centers 0.9
    // apart, while parallel rods above retain coincident barbed-end planes.
    actin.center[0] = {0.5, 0.0, 0.0};
    actin.center[1] = {-0.4, 0.065, 0.0};
    actin.direction[0] = {1.0, 0.0, 0.0};
    actin.direction[1] = {-1.0, 0.0, 0.0};
    actin.update_endpoints();
    const auto antiparallel_equilibrium = compute_aa_force_and_energy(
        actin, first, second, box, 10.0, 100.0, 0.2, 0.065, 0.1);
    for (double value : antiparallel_equilibrium) {
        assert(near(value, 0.0, 1e-7));
    }

    // Center-aligned antiparallel rods have barbed ends one length apart and
    // must not be mistaken for a stress-free Z-disc crosslink.
    actin.center[0] = {0.0, 0.0, 0.0};
    actin.center[1] = {0.0, 0.065, 0.0};
    actin.update_endpoints();
    const auto wrong_registry = compute_aa_force_and_energy(
        actin, first, second, box, 10.0, 100.0, 0.2, 0.065, 0.1);
    double wrong_norm = 0.0;
    for (double value : wrong_registry) wrong_norm += value * value;
    assert(wrong_norm > 1e-6);
    gsl_rng_free(rng);
    return 0;
}
