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
    Myosin myosin(1, 1.5, 0.05, box, rng);
    actin.set_periodic_axes({true, true, true});
    myosin.set_periodic_axes({true, true, true});

    actin.direction[0] = {1.0, 0.0, 0.0};
    actin.direction[1] = {-1.0, 0.0, 0.0};
    myosin.center[0] = {0.0, 0.0, 0.0};
    myosin.direction[0] = {1.0, 0.0, 0.0};
    myosin.update_endpoints();

    auto place_connector = [&](double distance) {
        // The two s=0 material sites form a short connector perpendicular to
        // the myosin axis whose exact segment distance is `distance`.
        actin.center[0] = {0.0, distance, -0.01};
        actin.center[1] = {0.0, distance, 0.01};
        actin.update_endpoints();
    };

    constexpr double hard_core_distance = 0.055;
    constexpr double force_cap = 40.0;

    // Width=0 is the original hard-core WCA contact.
    place_connector(0.050);
    const auto contact = compute_actn2_myosin_contact(
        actin, myosin, 0, 1, 0.0, 0.0, 0, box,
        hard_core_distance, 0.00428, force_cap);
    assert(contact.applied);
    assert(std::abs(contact.penetration - 0.005) < 1e-12);
    assert(contact.shell_penetration == 0.0);
    assert(contact.force > 0.0 && contact.force <= force_cap + 1e-12);
    assert(contact.energy > 0.0);
    assert(std::abs(contact.arm_weight_i + contact.arm_weight_j - 1.0) <
           1e-12);
    assert(near_vec(contact.force_on_actin_i + contact.force_on_actin_j,
                    -contact.force_on_myosin));

    place_connector(0.060);
    const auto separated = compute_actn2_myosin_contact(
        actin, myosin, 0, 1, 0.0, 0.0, 0, box,
        hard_core_distance, 0.00428, force_cap);
    assert(!separated.applied);
    assert(separated.force == 0.0);

    // A 15-nm, 5-kBT shoulder acts before the unchanged 55-nm hard core.
    // shell_energy is expressed directly in pN um (5*kBT at 310 K).
    constexpr double shell_width = 0.015;
    constexpr double shell_energy = 5.0 * 0.00428;
    place_connector(0.065);
    const auto shell = compute_actn2_myosin_contact(
        actin, myosin, 0, 1, 0.0, 0.0, 0, box,
        hard_core_distance, 0.00428, force_cap,
        shell_width, shell_energy);
    const double shell_k = 2.0 * shell_energy /
        (shell_width * shell_width);
    assert(shell.applied);
    assert(shell.penetration == 0.0);
    assert(std::abs(shell.shell_penetration - 0.005) < 1e-12);
    assert(std::abs(shell.force - shell_k * 0.005) < 1e-10);
    assert(std::abs(shell.energy - 0.5 * shell_k * 0.005 * 0.005) <
           1e-12);

    // Entering the corona must not move the geometric rupture surface: core
    // penetration remains zero until the connector crosses 55 nm.
    place_connector(0.055);
    const auto at_core = compute_actn2_myosin_contact(
        actin, myosin, 0, 1, 0.0, 0.0, 0, box,
        hard_core_distance, 0.00428, force_cap,
        shell_width, shell_energy);
    assert(at_core.penetration == 0.0);
    assert(std::abs(at_core.energy - shell_energy) < 1e-12);

    place_connector(0.050);
    const auto inside_core = compute_actn2_myosin_contact(
        actin, myosin, 0, 1, 0.0, 0.0, 0, box,
        hard_core_distance, 0.00428, force_cap,
        shell_width, shell_energy);
    assert(std::abs(inside_core.penetration - 0.005) < 1e-12);
    assert(inside_core.energy > shell_energy);

    place_connector(0.071);
    const auto outside_shell = compute_actn2_myosin_contact(
        actin, myosin, 0, 1, 0.0, 0.0, 0, box,
        hard_core_distance, 0.00428, force_cap,
        shell_width, shell_energy);
    assert(!outside_shell.applied);

    // With the current 50-nm thick-filament radius, a connector at the usual
    // 80-nm A-M centerline lies outside the 55+15=70-nm corona.  A 35-nm
    // shoulder reaches 90 nm and therefore supplies a genuine pre-contact
    // force at that geometry.  This test makes the range distinction explicit;
    // it does not turn connector contact into an axial overlap end-stop.
    place_connector(0.080);
    const auto current_range = compute_actn2_myosin_contact(
        actin, myosin, 0, 1, 0.0, 0.0, 0, box,
        hard_core_distance, 0.00428, force_cap,
        0.015, shell_energy);
    assert(!current_range.applied);
    const auto wider_range = compute_actn2_myosin_contact(
        actin, myosin, 0, 1, 0.0, 0.0, 0, box,
        hard_core_distance, 0.00428, force_cap,
        0.035, shell_energy);
    assert(wider_range.applied);
    assert(wider_range.penetration == 0.0);
    assert(std::abs(wider_range.shell_penetration - 0.010) < 1e-12);

    gsl_rng_free(rng);
    return 0;
}
