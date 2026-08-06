#include "components.h"
#include "interaction.h"

#include <cassert>
#include <cmath>
#include <gsl/gsl_rng.h>

namespace {
bool near_vec(const vec& lhs, const vec& rhs, double tolerance = 1e-10) {
    return (lhs - rhs).norm() <= tolerance;
}
}

int main() {
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    const std::vector<double> box{4.0, 4.0, 4.0};
    Filament actin(2, 1.0, box, rng);
    actin.set_periodic_axes({true, true, true});
    actin.center[0] = {0.0, 0.0, 0.0};
    actin.center[1] = {0.0, 0.0, 0.0};
    actin.direction[0] = {1.0, 0.0, 0.0};
    actin.direction[1] = {-1.0, 0.0, 0.0};
    actin.update_endpoints();

    constexpr double diameter = 0.030;
    constexpr double stiffness = 100.0;
    for (int first_status = 0; first_status <= 2; ++first_status) {
        for (int second_status = 0; second_status <= 2; ++second_status) {
            actin.cb_status[0] = first_status;
            actin.cb_status[1] = second_status;
            const auto result = compute_actin_repulsion(
                actin, 0, 1, box, diameter, stiffness);
            assert(result.applied);
            assert(result.force_on_first.norm() > 0.0);
            assert(near_vec(result.force_on_first, -result.force_on_second));

            const auto swapped = compute_actin_repulsion(
                actin, 1, 0, box, diameter, stiffness);
            assert(swapped.applied);
            assert(near_vec(result.force_on_first, swapped.force_on_second));
            assert(near_vec(result.force_on_second, swapped.force_on_first));
        }
    }

    // A resolved transverse separation should use the geometric normal rather
    // than the coincident-rod regularizer.
    actin.center[1] = {0.0, 0.020, 0.0};
    actin.update_endpoints();
    const auto separated = compute_actin_repulsion(
        actin, 0, 1, box, diameter, stiffness);
    assert(separated.applied);
    assert(separated.force_on_first.y < 0.0);
    assert(std::abs(separated.force_on_first.norm() - 1.0) < 1e-10);
    assert(near_vec(separated.force_on_first, -separated.force_on_second));

    gsl_rng_free(rng);
    return 0;
}
