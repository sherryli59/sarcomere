#include "components.h"
#include "interaction.h"

#include <cassert>
#include <cmath>
#include <gsl/gsl_rng.h>

namespace {
bool near(double a, double b, double tolerance = 1e-8) {
    return std::abs(a - b) <= tolerance;
}

bool near_vec(const vec& a, const vec& b, double tolerance = 1e-8) {
    return near(a.x, b.x, tolerance) && near(a.y, b.y, tolerance) &&
           near(a.z, b.z, tolerance);
}
}

int main() {
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    const std::vector<double> box{4.0, 4.0, 4.0};
    Filament actin(2, 1.0, box, rng);
    actin.set_periodic_axes({true, true, true});
    actin.center[0] = {0.0, 0.0, 0.0};
    actin.center[1] = {0.0, 0.05, 0.0};
    actin.direction[0] = {1.0, 0.0, 0.0};
    actin.direction[1] = {1.0, 0.0, 0.0};

    constexpr double radial_k = 100.0;
    constexpr double perpendicular_k = 10.0;
    constexpr double preferred = 0.05;
    const auto rest = compute_local_aa_link(
        actin, 0, 1, 0.0, 0.0, box, radial_k,
        perpendicular_k, 0.0, preferred);
    assert(rest.energy < 1e-12);
    assert(rest.force_on_first.norm() < 1e-8);
    assert(rest.torque_on_first.norm() < 1e-8);
    assert(rest.torque_on_second.norm() < 1e-8);

    // A stretched connector pulls the two actins together.
    actin.center[1] = {0.0, 0.06, 0.0};
    const auto stretched = compute_local_aa_link(
        actin, 0, 1, 0.0, 0.0, box, radial_k,
        perpendicular_k, 0.0, preferred);
    assert(stretched.force_on_first.y > 0.0);
    assert(near(stretched.force_on_first.y, 1.0, 1e-6));

    // Swapping link ownership preserves energy and reverses force while
    // exchanging the two torques.
    const auto reversed = compute_local_aa_link(
        actin, 1, 0, 0.0, 0.0, box, radial_k,
        perpendicular_k, 0.0, preferred);
    assert(near(stretched.energy, reversed.energy));
    assert(near_vec(stretched.force_on_first, -reversed.force_on_first));
    assert(near_vec(stretched.torque_on_first, reversed.torque_on_second));
    assert(near_vec(stretched.torque_on_second, reversed.torque_on_first));

    // Axially mismatched attachment sites make the linker oblique.  Both the
    // length and perpendicularity terms then generate restoring mechanics.
    actin.center[1] = {0.0, preferred, 0.0};
    const auto oblique = compute_local_aa_link(
        actin, 0, 1, -0.20, 0.20, box, radial_k,
        perpendicular_k, 0.0, preferred);
    assert(oblique.energy > 0.0);
    assert(oblique.force_on_first.x > 0.0);
    assert(oblique.torque_on_first.norm() > 1e-6);
    assert(oblique.torque_on_second.norm() > 1e-6);

    const auto reverse_oblique = compute_local_aa_link(
        actin, 0, 1, 0.20, -0.20, box, radial_k,
        perpendicular_k, 0.0, preferred);
    assert(reverse_oblique.force_on_first.x < 0.0);

    // Minimum-image wrapping gives identical mechanics across a boundary.
    actin.center[0] = {0.0, 1.98, 0.0};
    actin.center[1] = {0.0, -1.97, 0.0};
    const auto pbc = compute_local_aa_link(
        actin, 0, 1, 0.0, 0.0, box, radial_k,
        perpendicular_k, 0.0, preferred);
    assert(pbc.energy < 1e-12);
    assert(pbc.force_on_first.norm() < 1e-8);

    // Connector perpendicularity alone leaves an in-plane hinge. A mature
    // ACTN2 link's nematic term removes that hinge without selecting polarity.
    actin.center[0] = {0.0, 0.0, 0.0};
    actin.center[1] = {0.0, preferred, 0.0};
    actin.direction[0] = {1.0, 0.0, 0.0};
    constexpr double tilt = 0.2;
    actin.direction[1] = {std::cos(tilt), 0.0, std::sin(tilt)};
    const auto nematic = compute_local_aa_link(
        actin, 0, 1, 0.0, 0.0, box, 0.0, 0.0, 10.0, preferred);
    assert(nematic.energy > 0.0);
    assert(nematic.torque_on_first.norm() > 1e-6);
    assert(nematic.torque_on_second.norm() > 1e-6);

    // Status 1 has no global axial register even when a diagnostic target is
    // supplied: setting axial stiffness to zero exactly recovers the local
    // connector mechanics at arbitrary overlap.
    actin.center[0] = {0.0, 0.0, 0.0};
    actin.center[1] = {-0.55, preferred, 0.0};
    actin.direction[0] = {1.0, 0.0, 0.0};
    actin.direction[1] = {-1.0, 0.0, 0.0};
    const auto status1 = compute_local_aa_link(
        actin, 0, 1, -0.275, -0.275, box, radial_k,
        perpendicular_k, 0.0, preferred, 0.1, 0.0, 0.05);
    const auto status1_no_target = compute_local_aa_link(
        actin, 0, 1, -0.275, -0.275, box, radial_k,
        perpendicular_k, 0.0, preferred, 0.0, 1.0, 0.05);
    assert(near(status1.energy, status1_no_target.energy));
    assert(near_vec(status1.force_on_first, status1_no_target.force_on_first));
    assert(near_vec(status1.torque_on_first, status1_no_target.torque_on_first));
    assert(status1.registry_error > 0.3);

    // A status-2 link is force balanced at the 0.1 um antiparallel register.
    actin.center[1] = {-0.9, preferred, 0.0};
    const auto registered = compute_local_aa_link(
        actin, 0, 1, -0.45, -0.45, box, radial_k,
        0.0, 0.0, preferred, 0.1, 1.0, 0.05);
    assert(registered.registry_error < 1e-12);
    assert(registered.energy < 1e-12);
    assert(registered.force_on_first.norm() < 1e-8);
    assert(registered.torque_on_first.norm() < 1e-8);
    assert(registered.torque_on_second.norm() < 1e-8);

    // Over-overlap and under-overlap restore in opposite directions. Material
    // attachment coordinates are chosen so the local connector itself remains
    // at rest and only the global status-2 register contributes axial force.
    actin.center[1] = {-0.8, preferred, 0.0};
    const auto over_overlap = compute_local_aa_link(
        actin, 0, 1, -0.40, -0.40, box, radial_k,
        0.0, 0.0, preferred, 0.1, 1.0, 0.05);
    assert(near(over_overlap.registry_error, 0.1));
    assert(over_overlap.force_on_first.x > 0.0);
    actin.center[1] = {-1.0, preferred, 0.0};
    const auto under_overlap = compute_local_aa_link(
        actin, 0, 1, -0.50, -0.50, box, radial_k,
        0.0, 0.0, preferred, 0.1, 1.0, 0.05);
    assert(near(under_overlap.registry_error, 0.1));
    assert(under_overlap.force_on_first.x < 0.0);

    // The register energy and force remain bounded even for a very large
    // mismatch. For two scalar errors U_ax approaches k*sigma^2.
    actin.center[1] = {-0.3, preferred, 0.0};
    const auto far = compute_local_aa_link(
        actin, 0, 1, -0.15, -0.15, box, radial_k,
        0.0, 0.0, preferred, 0.1, 1.0, 0.05);
    assert(far.energy <= radial_k * 0.05 * 0.05 + 1e-12);
    assert(far.force_on_first.norm() < radial_k * 0.05);

    // Status-2 mechanics are invariant to pair ownership.
    actin.center[1] = {-0.8, preferred, 0.0};
    const auto forward_status2 = compute_local_aa_link(
        actin, 0, 1, -0.40, -0.40, box, radial_k,
        0.0, 0.0, preferred, 0.1, 1.0, 0.05);
    const auto reverse_status2 = compute_local_aa_link(
        actin, 1, 0, -0.40, -0.40, box, radial_k,
        0.0, 0.0, preferred, 0.1, 1.0, 0.05);
    assert(near(forward_status2.energy, reverse_status2.energy));
    assert(near_vec(forward_status2.force_on_first,
                    -reverse_status2.force_on_first));
    assert(near_vec(forward_status2.torque_on_first,
                    reverse_status2.torque_on_second));
    assert(near_vec(forward_status2.torque_on_second,
                    reverse_status2.torque_on_first));

    // Finite-difference translation gradient of the complete bounded energy.
    constexpr double eps = 1e-6;
    actin.center[0].x = eps;
    const double energy_plus = compute_local_aa_link(
        actin, 0, 1, -0.40, -0.40, box, radial_k,
        0.0, 0.0, preferred, 0.1, 1.0, 0.05).energy;
    actin.center[0].x = -eps;
    const double energy_minus = compute_local_aa_link(
        actin, 0, 1, -0.40, -0.40, box, radial_k,
        0.0, 0.0, preferred, 0.1, 1.0, 0.05).energy;
    actin.center[0].x = 0.0;
    const double fd_force_x = -(energy_plus - energy_minus) / (2.0 * eps);
    assert(near(fd_force_x, forward_status2.force_on_first.x, 1e-6));

    // Finite-difference orientational gradient includes the register torque.
    actin.center[1] = {-0.82, preferred, 0.02};
    actin.direction[1] = {-std::cos(0.1), 0.0, std::sin(0.1)};
    actin.direction[0] = {1.0, 0.0, 0.0};
    const auto tilted_register = compute_local_aa_link(
        actin, 0, 1, -0.40, -0.40, box, radial_k,
        0.0, 0.0, preferred, 0.1, 1.0, 0.05);
    actin.direction[0] = {1.0, 0.0, eps};
    actin.direction[0].normalize();
    const double orient_plus = compute_local_aa_link(
        actin, 0, 1, -0.40, -0.40, box, radial_k,
        0.0, 0.0, preferred, 0.1, 1.0, 0.05).energy;
    actin.direction[0] = {1.0, 0.0, -eps};
    actin.direction[0].normalize();
    const double orient_minus = compute_local_aa_link(
        actin, 0, 1, -0.40, -0.40, box, radial_k,
        0.0, 0.0, preferred, 0.1, 1.0, 0.05).energy;
    const double fd_torque_z = -(orient_plus - orient_minus) / (2.0 * eps);
    assert(near(fd_torque_z, tilted_register.torque_on_first.z, 1e-5));

    gsl_rng_free(rng);
    return 0;
}
