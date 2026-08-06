#include "rod_mobility.h"

#include <cassert>
#include <cmath>
#include <random>

int main() {
    rod_mobility::Diffusion diffusion;
    diffusion.enabled = true;
    diffusion.parallel = 0.30;
    diffusion.perpendicular = 0.12;
    diffusion.tumble = 0.08;
    diffusion.roll = 0.20;
    utils::vec axis{1.0, 2.0, -1.0};
    axis.normalize();
    utils::vec e1 = rod_mobility::project_perpendicular({0.0, 0.0, 1.0}, axis);
    e1.normalize();
    utils::vec e2 = axis.cross(e1);
    e2.normalize();

    constexpr int samples = 250000;
    constexpr double dt = 2e-4;
    std::mt19937_64 rng(1701);
    std::normal_distribution<double> normal(0.0, 1.0);
    double parallel2 = 0.0;
    double perpendicular12 = 0.0;
    double perpendicular22 = 0.0;
    double tumble12 = 0.0;
    double tumble22 = 0.0;
    double cross12 = 0.0;
    for (int sample = 0; sample < samples; ++sample) {
        const utils::vec gaussian{normal(rng), normal(rng), normal(rng)};
        const utils::vec displacement = std::sqrt(2.0 * dt) *
            rod_mobility::apply_sqrt_diffusion_tensor(gaussian, axis, diffusion);
        const double xp = displacement.dot(axis);
        const double x1 = displacement.dot(e1);
        const double x2 = displacement.dot(e2);
        parallel2 += xp * xp;
        perpendicular12 += x1 * x1;
        perpendicular22 += x2 * x2;
        cross12 += x1 * x2;

        const utils::vec rotational =
            std::sqrt(2.0 * diffusion.tumble * dt) *
            rod_mobility::project_perpendicular(gaussian, axis);
        const double r1 = rotational.dot(e1);
        const double r2 = rotational.dot(e2);
        tumble12 += r1 * r1;
        tumble22 += r2 * r2;
    }
    const auto relative_error = [](double measured, double expected) {
        return std::abs(measured / expected - 1.0);
    };
    assert(relative_error(parallel2 / samples,
                          2.0 * diffusion.parallel * dt) < 0.015);
    assert(relative_error(perpendicular12 / samples,
                          2.0 * diffusion.perpendicular * dt) < 0.015);
    assert(relative_error(perpendicular22 / samples,
                          2.0 * diffusion.perpendicular * dt) < 0.015);
    assert(std::abs(cross12 / samples) < 0.01 *
           (2.0 * diffusion.perpendicular * dt));
    assert(relative_error(tumble12 / samples,
                          2.0 * diffusion.tumble * dt) < 0.015);
    assert(relative_error(tumble22 / samples,
                          2.0 * diffusion.tumble * dt) < 0.015);

    // Einstein drift in each channel uses the same diffusion tensor as noise.
    const double beta = rod_mobility::beta(310.0);
    const utils::vec force = 2.3 * axis - 1.7 * e1;
    const utils::vec drift = beta * dt *
        rod_mobility::apply_diffusion_tensor(force, axis, diffusion);
    assert(relative_error(drift.dot(axis),
                          beta * dt * diffusion.parallel * 2.3) < 1e-12);
    assert(relative_error(drift.dot(e1),
                          beta * dt * diffusion.perpendicular * -1.7) < 1e-12);
    assert(std::abs(drift.dot(e2)) < 1e-14);
    return 0;
}
