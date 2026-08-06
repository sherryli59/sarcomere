#include "rod_mobility.h"

#include <cassert>
#include <cmath>
#include <limits>

namespace {

bool close(double actual, double expected, double tolerance = 1e-6) {
    return std::abs(actual - expected) <= tolerance * std::max(1.0, std::abs(expected));
}

}

int main() {
    assert(close(rod_mobility::beta(310.0), 233.6442102, 1e-6));

    const auto thin = rod_mobility::finite_cylinder(
        1.0, 0.010, 310.0, 0.001);
    assert(close(thin.parallel, 3.0024, 2e-4));
    assert(close(thin.perpendicular, 1.8547, 2e-4));
    assert(close(thin.tumble, 16.152, 2e-4));

    utils::vec axis{1.0, 0.0, 0.0};
    utils::vec value{2.0, 3.0, 4.0};
    const auto transformed = rod_mobility::apply_diffusion_tensor(value, axis, thin);
    assert(close(transformed.x, 2.0 * thin.parallel));
    assert(close(transformed.y, 3.0 * thin.perpendicular));
    assert(close(transformed.z, 4.0 * thin.perpendicular));

    const auto rotated = rod_mobility::sphere_exp(axis, {0.0, M_PI / 2.0, 0.0});
    assert(std::abs(rotated.x) < 1e-12);
    assert(close(rotated.y, 1.0));
    assert(std::abs(rotated.z) < 1e-12);
    assert(close(rotated.norm(), 1.0));

    // The realized angular increment is left untouched below the trust radius
    // and is scaled exactly to it above the limit. Disabled/infinite limits
    // preserve the original increment.
    const utils::vec angular_step{0.0, 0.012, 0.016};
    assert(close(angular_step.norm(), 0.020));
    assert(close(rod_mobility::angular_increment_scale(
                     angular_step, 0.010, 1.0), 0.5));
    assert(close(rod_mobility::angular_increment_scale(
                     angular_step, 0.020, 1.0), 1.0));
    assert(close(rod_mobility::angular_increment_scale(
                     angular_step, 0.010, 2.0), 1.0));
    assert(close(rod_mobility::angular_increment_scale(
                     angular_step, 0.0, 1.0), 1.0));
    assert(close(rod_mobility::angular_increment_scale(
                     angular_step, std::numeric_limits<double>::infinity(),
                     1.0), 1.0));
    return 0;
}
