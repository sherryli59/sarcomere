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

    // A polarity-aware pointed-end contact retains axial ordering. Generic
    // capsule sterics cannot do that after two collinear rods cross because
    // their closest-point normal is then degenerate.
    actin.direction[0] = {1.0, 0.0, 0.0};
    actin.direction[1] = {-1.0, 0.0, 0.0};
    actin.center[0] = {-0.495, 0.0, 0.0};
    actin.center[1] = {0.495, 0.0, 0.0};
    actin.update_endpoints();
    const auto end_contact = compute_opposing_pointed_end_contact(
        actin, 0, 1, box, 0.012, 0.00428, 10.0);
    assert(end_contact.applied);
    assert(end_contact.axial_gap < 0.0);
    assert(end_contact.contact_penetration > 0.012);
    assert(end_contact.force_on_first.x < 0.0);
    assert(near_vec(end_contact.force_on_first,
                    -end_contact.force_on_second));
    assert(end_contact.force_on_first.norm() <= 10.0 + 1e-12);
    const auto swapped_end_contact = compute_opposing_pointed_end_contact(
        actin, 1, 0, box, 0.012, 0.00428, 10.0);
    assert(near_vec(end_contact.force_on_first,
                    swapped_end_contact.force_on_second));

    actin.center[0] = {-0.510, 0.0, 0.0};
    actin.center[1] = {0.510, 0.0, 0.0};
    actin.update_endpoints();
    const auto separated_caps = compute_opposing_pointed_end_contact(
        actin, 0, 1, box, 0.012, 0.00428, 10.0);
    assert(!separated_caps.applied);

    actin.center[0] = {-0.495, 0.0, 0.0};
    actin.center[1] = {0.495, 0.020, 0.0};
    actin.update_endpoints();
    const auto different_lanes = compute_opposing_pointed_end_contact(
        actin, 0, 1, box, 0.012, 0.00428, 10.0);
    assert(!different_lanes.applied);

    // Off-center nonparallel capsule intersections must not be represented by
    // center forces alone.  The contact force is reciprocal, uses the normal
    // to both rod axes, and rotates the rod whose closest point is off-center.
    actin.center[0] = {0.0, 0.0, 0.0};
    actin.center[1] = {0.25, 0.0, 0.0};
    actin.direction[0] = {1.0, 0.0, 0.0};
    actin.direction[1] = {0.0, 1.0, 0.0};
    actin.update_endpoints();
    const auto crossing = compute_actin_repulsion(
        actin, 0, 1, box, diameter, stiffness);
    assert(crossing.applied);
    assert(crossing.segment_distance < 1e-7);
    assert(std::abs(crossing.force_on_first.z) > 0.0);
    assert(near_vec(crossing.force_on_first, -crossing.force_on_second));
    assert(crossing.torque_on_first.norm() > 0.0);
    const auto crossing_swapped = compute_actin_repulsion(
        actin, 1, 0, box, diameter, stiffness);
    assert(near_vec(crossing.force_on_first,
                    crossing_swapped.force_on_second));
    assert(near_vec(crossing.torque_on_first,
                    crossing_swapped.torque_on_second));
    constexpr double mobility_step = 1e-4;
    actin.center[0] += mobility_step * crossing.force_on_first;
    actin.center[1] += mobility_step * crossing.force_on_second;
    actin.direction[0] += mobility_step * crossing.torque_on_first;
    actin.direction[1] += mobility_step * crossing.torque_on_second;
    actin.direction[0].normalize();
    actin.direction[1].normalize();
    actin.update_endpoints();
    const auto separated_crossing = compute_actin_repulsion(
        actin, 0, 1, box, diameter, stiffness);
    assert(separated_crossing.segment_distance > crossing.segment_distance);

    Myosin myosin(2, 1.5, 0.025, box, rng);
    myosin.set_periodic_axes({true, true, true});
    myosin.center[0] = {0.0, 0.0, 0.0};
    myosin.center[1] = {0.30, 0.0, 0.0};
    myosin.direction[0] = {1.0, 0.0, 0.0};
    myosin.direction[1] = {0.0, 1.0, 0.0};
    myosin.update_endpoints();
    utils::MoleculeConnection actins_per_myosin(2);
    const auto myosin_crossing = compute_myosin_repulsion(
        actin, myosin, 0, 1, box, 0, actins_per_myosin, stiffness);
    assert(myosin_crossing.applied);
    assert(myosin_crossing.segment_distance < 1e-7);
    assert(near_vec(myosin_crossing.force_on_first,
                    -myosin_crossing.force_on_second));
    assert(myosin_crossing.torque_on_first.norm() > 0.0);
    actin.cb_status[0] = 2;
    actins_per_myosin.addConnection(0, 0);
    const auto occupied_myosin_crossing = compute_myosin_repulsion(
        actin, myosin, 0, 1, box, 0, actins_per_myosin, stiffness);
    assert(near_vec(occupied_myosin_crossing.force_on_first,
                    myosin_crossing.force_on_first));
    assert(near_vec(occupied_myosin_crossing.force_on_first,
                    -occupied_myosin_crossing.force_on_second));

    // The force-level penalty is complemented by an overdamped hard-contact
    // projection, so a finite time step cannot leave two thick rods crossed.
    const vec midpoint_before =
        0.5 * (static_cast<vec>(myosin.center[0]) +
               static_cast<vec>(myosin.center[1]));
    const auto projected = project_myosin_hard_contact(
        myosin, 0, 1, box, 0, 2.0 * myosin.radius);
    assert(projected.applied);
    assert(projected.distance_before < 1e-7);
    const auto after_projection = compute_myosin_repulsion(
        actin, myosin, 0, 1, box, 0, actins_per_myosin, stiffness);
    assert(!after_projection.applied);
    const vec midpoint_after =
        0.5 * (static_cast<vec>(myosin.center[0]) +
               static_cast<vec>(myosin.center[1]));
    assert(near_vec(midpoint_before, midpoint_after));

    // An explicitly fixed first rod receives no position correction; the
    // mobile partner receives the full separating displacement.
    myosin.center[0] = {0.0, 0.0, 0.0};
    myosin.center[1] = {0.30, 0.0, 0.0};
    myosin.direction[0] = {1.0, 0.0, 0.0};
    myosin.direction[1] = {0.0, 1.0, 0.0};
    myosin.update_endpoints();
    const vec fixed_center = myosin.center[0];
    const auto fixed_projection = project_myosin_hard_contact(
        myosin, 0, 1, box, 1, 2.0 * myosin.radius);
    assert(fixed_projection.applied);
    assert(near_vec(fixed_center, myosin.center[0]));
    const auto fixed_after = compute_myosin_repulsion(
        actin, myosin, 0, 1, box, 1, actins_per_myosin, stiffness);
    assert(!fixed_after.applied);

    // With rotational mobility enabled, an off-center hard contact also
    // changes the directors while preserving the pair midpoint.
    myosin.center[0] = {0.0, 0.0, 0.0};
    myosin.center[1] = {0.30, 0.0, 0.0};
    myosin.direction[0] = {1.0, 0.0, 0.0};
    myosin.direction[1] = {0.0, 1.0, 0.0};
    myosin.update_endpoints();
    const vec rotating_midpoint_before =
        0.5 * (static_cast<vec>(myosin.center[0]) +
               static_cast<vec>(myosin.center[1]));
    const vec rotating_direction_before = myosin.direction[0];
    const auto rotating_projection = project_myosin_hard_contact(
        myosin, 0, 1, box, 0, 2.0 * myosin.radius, 0.05, 0.06);
    assert(rotating_projection.applied);
    assert((static_cast<vec>(myosin.direction[0]) -
            rotating_direction_before).norm() > 0.0);
    const vec rotating_midpoint_after =
        0.5 * (static_cast<vec>(myosin.center[0]) +
               static_cast<vec>(myosin.center[1]));
    assert(near_vec(rotating_midpoint_before, rotating_midpoint_after));

    // Mixed actin-myosin contacts carry the same reciprocal force and
    // contact-point orientation drive used by the production dynamics.
    myosin.center[0] = {0.0, 0.0, 0.0};
    myosin.direction[0] = {1.0, 0.0, 0.0};
    myosin.update_endpoints();
    actin.center[0] = {0.30, 0.0, 0.0};
    actin.direction[0] = {0.0, 1.0, 0.0};
    actin.update_endpoints();
    const auto mixed_crossing =
        compute_actin_myosin_radial_repulsion_result(
            actin, myosin, 0, 0, box, 0.04, stiffness);
    assert(mixed_crossing.applied);
    assert(mixed_crossing.segment_distance < 1e-7);
    assert(near_vec(mixed_crossing.force_on_first,
                    -mixed_crossing.force_on_second));
    assert(mixed_crossing.torque_on_second.norm() > 0.0);

    gsl_rng_free(rng);
    return 0;
}
