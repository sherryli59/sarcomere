#include "am_motor_kinetics.h"
#include "components.h"
#include "interaction.h"
#include "geometry.h"

#include <cassert>
#include <cmath>
#include <gsl/gsl_rng.h>

namespace {
bool near(double a, double b, double tol = 1e-10) {
    return std::abs(a - b) <= tol;
}

vec rotate_about(const vec& value, vec axis, double angle) {
    axis.normalize();
    return std::cos(angle) * value +
           std::sin(angle) * axis.cross(value) +
           (1.0 - std::cos(angle)) * axis * axis.dot(value);
}
}

int main() {
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, 17);
    const std::vector<double> box{10.0, 10.0, 10.0};
    Filament actin(1, 1.0, box, rng);
    Myosin myosin(1, 1.0, 0.01, box, rng);
    actin.set_periodic_axes({false, false, false});
    myosin.set_periodic_axes({false, false, false});

    actin.center[0] = {0.0, 0.05, 0.0};
    actin.direction[0] = {1.0, 0.0, 0.0};
    myosin.center[0] = {0.0, 0.0, 0.0};
    myosin.direction[0] = {1.0, 0.0, 0.0};
    actin.update_endpoints();
    myosin.update_endpoints();

    const auto tether = compute_am_proximity_tether(
        actin, myosin, 0, 0, box, 10.0, 0.02);
    assert(near(tether.transverse_distance, 0.05));
    assert(near(tether.force_on_actin.x, 0.0));
    assert(tether.force_on_actin.y < 0.0);
    assert(near(tether.force_on_actin.x + tether.force_on_myosin.x, 0.0));
    assert(near(tether.force_on_actin.y + tether.force_on_myosin.y, 0.0));
    assert(near(tether.force_on_actin.z + tether.force_on_myosin.z, 0.0));

    // Axial translation changes neither transverse distance nor tether force.
    actin.center[0] = {0.2, 0.05, 0.0};
    actin.update_endpoints();
    const auto shifted = compute_am_proximity_tether(
        actin, myosin, 0, 0, box, 10.0, 0.02);
    assert(near(shifted.transverse_distance, tether.transverse_distance));
    assert(near(shifted.force_on_actin.x, 0.0));
    assert(near(shifted.force_on_actin.y, tether.force_on_actin.y));

    // The objective two-spring result exposes energy components used by the
    // continuous detachment hazard. Their sum remains the mechanical energy.
    const auto strained_head = compute_am_two_spring_head(
        actin, myosin, 0, 0, box, 0.0, 0.0, 0.0, myosin.radius,
        nullptr, nullptr, 0.03, 20.0, 0.0, 10.0, 0.0, 4.5, false);
    assert(strained_head.radial_energy >= 0.0);
    assert(strained_head.axial_energy >= 0.0);
    assert(near(strained_head.energy,
                strained_head.radial_energy + strained_head.axial_energy));
    assert(strained_head.radial_energy > 0.0);
    assert(strained_head.radial_force_on_actin.norm() > 0.0);
    assert(near((strained_head.force_on_actin +
                 strained_head.force_on_myosin).norm(), 0.0));
    const double unloaded_off_rate = 200.0;
    const double loaded_off_rate =
        am_motor_kinetics::strain_accelerated_rate(
            unloaded_off_rate, 240.0 * strained_head.radial_energy, 0.1);
    assert(loaded_off_rate > unloaded_off_rate);

    // A cross-bridge force applied away from a thick-filament center produces
    // tau = r x F.  The orientation integrator consumes tau x axis, which is
    // tangent to the unit sphere.  A center-applied force produces no moment.
    const vec filament_axis{1.0, 0.0, 0.0};
    const auto offcenter = compute_am_offcenter_torque(
        {0.0, 0.0, 0.0}, filament_axis,
        {0.25, 0.0, 0.0}, {0.0, 2.0, 0.0}, box,
        {false, false, false});
    assert(near(offcenter.lever_arm.x, 0.25));
    assert(near(offcenter.physical_torque.z, 0.5));
    assert(near(offcenter.orientation_drive.x, 0.0));
    assert(near(offcenter.orientation_drive.y, 0.5));
    assert(near(offcenter.orientation_drive.dot(filament_axis), 0.0));

    const auto centered = compute_am_offcenter_torque(
        {0.0, 0.0, 0.0}, filament_axis,
        {0.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, box,
        {false, false, false});
    assert(near(centered.physical_torque.norm(), 0.0));
    assert(near(centered.orientation_drive.norm(), 0.0));

    // Aggregate occupied-head compliance is nematic and flat-bottomed. It
    // contributes no drive inside the allowed band, but both directors rotate
    // toward their nearest parallel/antiparallel alignment outside it.
    constexpr double pi = 3.14159265358979323846;
    const vec small_tilt{
        std::cos(pi / 180.0), std::sin(pi / 180.0), 0.0};
    const auto inside_angular_slack =
        compute_am_nematic_flat_bottom_tether(
            small_tilt, filament_axis, 10.0, 2.0 * pi / 180.0);
    assert(!inside_angular_slack.active);
    assert(near(inside_angular_slack.energy, 0.0));
    assert(near(inside_angular_slack.drive_on_actin.norm(), 0.0));

    const vec large_tilt{
        std::cos(8.0 * pi / 180.0),
        std::sin(8.0 * pi / 180.0), 0.0};
    const auto outside_angular_slack =
        compute_am_nematic_flat_bottom_tether(
            large_tilt, filament_axis, 10.0, 2.0 * pi / 180.0);
    assert(outside_angular_slack.active);
    assert(outside_angular_slack.energy > 0.0);
    assert(project_nematic_alignment_drive(
               large_tilt, filament_axis,
               outside_angular_slack.drive_on_actin) > 0.0);
    assert(project_nematic_alignment_drive(
               filament_axis, large_tilt,
               outside_angular_slack.drive_on_myosin) > 0.0);
    const auto antiparallel_angular_slack =
        compute_am_nematic_flat_bottom_tether(
            -1.0 * large_tilt, filament_axis,
            10.0, 2.0 * pi / 180.0);
    assert(near(antiparallel_angular_slack.energy,
                outside_angular_slack.energy));

    // Two productively oriented heads on opposite halves of a bipolar thick
    // filament both rotate its axis toward a parallel local actin bundle.
    const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
    const vec tilted_axis{inv_sqrt2, inv_sqrt2, 0.0};
    const vec negative_head = -0.3 * tilted_axis;
    const vec positive_head = 0.3 * tilted_axis;
    const auto negative_half = compute_am_offcenter_torque(
        {0.0, 0.0, 0.0}, tilted_axis, negative_head,
        {-1.0, 0.0, 0.0}, box, {false, false, false});
    const auto positive_half = compute_am_offcenter_torque(
        {0.0, 0.0, 0.0}, tilted_axis, positive_head,
        {1.0, 0.0, 0.0}, box, {false, false, false});
    const vec bundle_drive = negative_half.orientation_drive +
                             positive_half.orientation_drive;
    assert(project_nematic_alignment_drive(
               tilted_axis, filament_axis,
               negative_half.orientation_drive) > 0.0);
    assert(project_nematic_alignment_drive(
               tilted_axis, filament_axis,
               positive_half.orientation_drive) > 0.0);
    assert(project_nematic_alignment_drive(
               tilted_axis, filament_axis, -bundle_drive) < 0.0);
    assert(near(project_nematic_alignment_drive(
               filament_axis, filament_axis, bundle_drive), 0.0));
    vec aligned_axis = tilted_axis + 1e-3 * bundle_drive;
    aligned_axis.normalize();
    assert(std::abs(aligned_axis.x) > std::abs(tilted_axis.x));

    // Directional attachment geometry chooses the head-bearing region on the
    // correct half, so the active-force lever arm cannot enter the bare zone.
    actin.center[0] = {-0.25, 0.05, 0.0};
    actin.direction[0] = {1.0, 0.0, 0.0};
    myosin.center[0] = {0.0, 0.0, 0.0};
    myosin.direction[0] = {1.0, 0.0, 0.0};
    actin.update_endpoints();
    myosin.update_endpoints();
    const auto directional_tether = compute_am_proximity_tether(
        actin, myosin, 0, 0, box, 10.0, 0.02, true);
    const vec directional_arm = directional_tether.myosin_attachment -
                                myosin.center[0];
    assert(directional_arm.x <=
           -0.5 * geometry::am_bare_zone_fraction * myosin.length + 1e-10);

    // Independent heads use distinct material anchors along the thick
    // filament. Their passive tethers remain action-reaction pairs and do not
    // create axial propulsion when the filaments are parallel.
    const auto head_inner = compute_am_head_tether(
        actin, myosin, 0, 0, box, -0.20, 10.0, 0.02);
    const auto head_outer = compute_am_head_tether(
        actin, myosin, 0, 0, box, -0.45, 10.0, 0.02);
    assert(near(head_inner.myosin_attachment.x, -0.20));
    assert(near(head_outer.myosin_attachment.x, -0.45));
    assert(near(head_inner.force_on_actin.x, 0.0));
    assert(near(head_outer.force_on_actin.x, 0.0));
    assert(near((head_inner.force_on_actin +
                 head_inner.force_on_myosin).norm(), 0.0));
    assert(near((head_outer.force_on_actin +
                 head_outer.force_on_myosin).norm(), 0.0));

    // Pair-level radial overlap is not sufficient for attachment: a material
    // head whose axial anchor lies beyond the finite actin segment cannot
    // reach it even when the two rod axes are radially close.  The full 3-D
    // reach metric retains that axial mismatch.
    actin.center[0] = {-0.55, 0.05, 0.0};
    actin.update_endpoints();
    const auto axially_reachable = compute_am_head_tether(
        actin, myosin, 0, 0, box, -0.45, 0.0, 0.055,
        0.0, myosin.radius, nullptr, nullptr, 0.015);
    const auto axially_unreachable = compute_am_head_tether(
        actin, myosin, 0, 0, box, 0.45, 0.0, 0.055,
        0.0, myosin.radius, nullptr, nullptr, 0.015);
    assert(axially_reachable.head_reach_surface_distance < 0.06);
    assert(axially_unreachable.head_reach_surface_distance > 0.25);
    assert(axially_reachable.projection_inside_actin);
    assert(!axially_unreachable.projection_inside_actin);
    assert(std::abs(axially_reachable.unbounded_actin_coordinate) <
           0.5 * actin.length);
    assert(std::abs(axially_unreachable.unbounded_actin_coordinate) >
           0.5 * actin.length);
    assert(near(axially_reachable.head_surface_distance,
                axially_unreachable.head_surface_distance, 1e-10));

    // A thick filament is nematic, not polar. Reversing its stored director,
    // swapping the axial half label, and transforming the transported material
    // frame must describe the same physical head and produce the same tether.
    actin.center[0] = {-0.2, 0.06, 0.01};
    actin.direction[0] = {1.0, 0.0, 0.0};
    myosin.center[0] = {0.0, 0.0, 0.0};
    myosin.direction[0] = {1.0, 0.0, 0.0};
    actin.update_endpoints();
    myosin.update_endpoints();
    const vec material_e1{0.0, 1.0, 0.0};
    const vec material_e2{0.0, 0.0, 1.0};
    const double azimuth = 0.37;
    const auto nematic_original = compute_am_head_tether(
        actin, myosin, 0, 0, box, -0.3, 12.0, 0.02,
        azimuth, myosin.radius, &material_e1, &material_e2, 0.01);
    myosin.direction[0] = {-1.0, 0.0, 0.0};
    myosin.update_endpoints();
    const vec flipped_e1 = material_e1;
    const vec flipped_e2 = -1.0 * material_e2;
    const auto nematic_flipped = compute_am_head_tether(
        actin, myosin, 0, 0, box, 0.3, 12.0, 0.02,
        -azimuth, myosin.radius, &flipped_e1, &flipped_e2, 0.01);
    assert((nematic_original.myosin_attachment -
            nematic_flipped.myosin_attachment).norm() < 1e-10);
    assert((nematic_original.actin_attachment -
            nematic_flipped.actin_attachment).norm() < 1e-10);
    assert((nematic_original.force_on_actin -
            nematic_flipped.force_on_actin).norm() < 1e-10);

    // The new radial geometry distinguishes centerline separation, surface
    // gap, and hard-core exclusion. Its attached-head force vanishes at the
    // requested actin-myosin separation and changes sign across it.
    constexpr double actin_radius = 0.015;
    constexpr double preferred_centerline = 0.055;
    actin.center[0] = {-0.25, preferred_centerline, 0.0};
    actin.update_endpoints();
    const auto radial_equilibrium = compute_am_head_tether(
        actin, myosin, 0, 0, box, -0.20, 10.0,
        preferred_centerline, 0.0, myosin.radius, nullptr, nullptr,
        actin_radius);
    assert(near(radial_equilibrium.centerline_distance,
                preferred_centerline));
    assert(near(radial_equilibrium.surface_gap,
                preferred_centerline - myosin.radius - actin_radius));
    assert(radial_equilibrium.force_on_actin.norm() < 1e-10);

    actin.center[0].y = 0.065;
    actin.update_endpoints();
    const auto radial_stretched = compute_am_head_tether(
        actin, myosin, 0, 0, box, -0.20, 10.0,
        preferred_centerline, 0.0, myosin.radius, nullptr, nullptr,
        actin_radius);
    assert(radial_stretched.force_on_actin.y < 0.0);
    assert(near((radial_stretched.force_on_actin +
                 radial_stretched.force_on_myosin).norm(), 0.0));

    // The proximity tether is force-free throughout its flat-bottom band and
    // only restores separation after the configured slack is exhausted.
    const auto flat_bottom_free = compute_am_head_tether(
        actin, myosin, 0, 0, box, 0.20, 10.0,
        preferred_centerline, 0.0, 0.0, nullptr, nullptr,
        actin_radius, 0.015);
    assert(flat_bottom_free.force_on_actin.norm() < 1e-10);
    const auto flat_bottom_loaded = compute_am_head_tether(
        actin, myosin, 0, 0, box, 0.20, 10.0,
        preferred_centerline, 0.0, 0.0, nullptr, nullptr,
        actin_radius, 0.005);
    assert(flat_bottom_loaded.force_on_actin.norm() > 1e-6);

    // At binding, a crown may swivel toward a reachable actin, but never by
    // more than its prescribed angular range. The chosen azimuth is then
    // stored by the head rather than continually tracking the actin.
    myosin.direction[0] = {1.0, 0.0, 0.0};
    myosin.update_endpoints();
    actin.center[0] = {-0.20, 0.0, 0.06};
    actin.direction[0] = {1.0, 0.0, 0.0};
    actin.update_endpoints();
    const double selected_azimuth = select_am_binding_azimuth(
        actin, myosin, 0, 0, box, -0.20, 0.0, M_PI / 6.0,
        material_e1, material_e2);
    assert(near(selected_azimuth, M_PI / 6.0, 1e-10));

    actin.center[0].y = 0.020;
    actin.center[0].z = 0.0;
    actin.update_endpoints();
    const vec hard_core = compute_actin_myosin_radial_repulsion(
        actin, myosin, 0, 0, box,
        myosin.radius + actin_radius, 1000.0);
    assert(hard_core.y > 0.0);

    // The production two-spring head is derived from one objective energy.
    // Its action-reaction pair and physical moments therefore have zero net
    // force and zero net torque about any laboratory origin.
    actin.center[0] = {0.03, 0.075, -0.012};
    actin.direction[0] = {1.0, 0.18, 0.08};
    actin.direction[0].normalize();
    myosin.center[0] = {-0.02, -0.01, 0.015};
    myosin.direction[0] = {1.0, -0.12, 0.05};
    myosin.direction[0].normalize();
    actin.update_endpoints();
    myosin.update_endpoints();
    vec objective_e1 = static_cast<vec>(myosin.direction[0]).cross(
        {0.0, 0.0, 1.0});
    objective_e1.normalize();
    vec objective_e2 = static_cast<vec>(myosin.direction[0]).cross(
        objective_e1);
    objective_e2.normalize();
    constexpr double material_s = 0.11;
    constexpr double anchor_s = -0.22;
    constexpr double head_phi = 0.41;
    constexpr double head_radius = 0.01;
    constexpr double radial_preferred = 0.028;
    constexpr double radial_k = 35.0;
    constexpr double radial_slack = 0.001;
    constexpr double axial_k = 18.0;
    constexpr double stroke_progress = 0.045;
    constexpr double large_stall = 20.0;
    const auto objective_head = compute_am_two_spring_head(
        actin, myosin, 0, 0, box, material_s, anchor_s, head_phi,
        head_radius, &objective_e1, &objective_e2,
        radial_preferred, radial_k, radial_slack, axial_k,
        stroke_progress, large_stall, true);
    assert((objective_head.force_on_actin +
            objective_head.force_on_myosin).norm() < 1e-11);
    const vec arbitrary_origin{0.17, -0.08, 0.04};
    const vec total_moment =
        (static_cast<vec>(actin.center[0]) - arbitrary_origin).cross(
            objective_head.force_on_actin) +
        objective_head.physical_torque_on_actin +
        (static_cast<vec>(myosin.center[0]) - arbitrary_origin).cross(
            objective_head.force_on_myosin) +
        objective_head.physical_torque_on_myosin;
    assert(total_moment.norm() < 1e-10);

    // Finite differences independently verify all force components and both
    // rotational moments. This catches sign errors that a trajectory can hide.
    const vec saved_actin_center = actin.center[0];
    const vec saved_actin_axis = actin.direction[0];
    const vec saved_myosin_center = myosin.center[0];
    const vec saved_myosin_axis = myosin.direction[0];
    auto energy = [&]() {
        return compute_am_two_spring_head(
            actin, myosin, 0, 0, box, material_s, anchor_s, head_phi,
            head_radius, &objective_e1, &objective_e2,
            radial_preferred, radial_k, radial_slack, axial_k,
            stroke_progress, large_stall, true).energy;
    };
    constexpr double fd_eps = 1e-7;
    const std::array<vec, 3> basis{{
        {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
    for (const vec& direction : basis) {
        actin.center[0] = saved_actin_center + fd_eps * direction;
        actin.update_endpoints();
        const double plus = energy();
        actin.center[0] = saved_actin_center - fd_eps * direction;
        actin.update_endpoints();
        const double minus = energy();
        actin.center[0] = saved_actin_center;
        actin.update_endpoints();
        assert(near((plus - minus) / (2.0 * fd_eps),
                    -objective_head.force_on_actin.dot(direction), 2e-7));

        myosin.center[0] = saved_myosin_center + fd_eps * direction;
        myosin.update_endpoints();
        const double myosin_plus = energy();
        myosin.center[0] = saved_myosin_center - fd_eps * direction;
        myosin.update_endpoints();
        const double myosin_minus = energy();
        myosin.center[0] = saved_myosin_center;
        myosin.update_endpoints();
        assert(near((myosin_plus - myosin_minus) / (2.0 * fd_eps),
                    -objective_head.force_on_myosin.dot(direction), 2e-7));
    }
    const vec rotation_axis = vec{0.3, -0.4, 0.7};
    actin.direction[0] = rotate_about(saved_actin_axis, rotation_axis, fd_eps);
    actin.update_endpoints();
    const double actin_rot_plus = energy();
    actin.direction[0] = rotate_about(saved_actin_axis, rotation_axis, -fd_eps);
    actin.update_endpoints();
    const double actin_rot_minus = energy();
    actin.direction[0] = saved_actin_axis;
    actin.update_endpoints();
    vec unit_rotation_axis = rotation_axis;
    unit_rotation_axis.normalize();
    assert(near((actin_rot_plus - actin_rot_minus) / (2.0 * fd_eps),
                -objective_head.physical_torque_on_actin.dot(
                    unit_rotation_axis), 3e-7));

    myosin.direction[0] = rotate_about(saved_myosin_axis, rotation_axis, fd_eps);
    objective_e1 = rotate_about(objective_e1, rotation_axis, fd_eps);
    objective_e2 = rotate_about(objective_e2, rotation_axis, fd_eps);
    myosin.update_endpoints();
    const double myosin_rot_plus = energy();
    myosin.direction[0] = rotate_about(saved_myosin_axis, rotation_axis, -fd_eps);
    objective_e1 = rotate_about(
        rotate_about(objective_e1, rotation_axis, -fd_eps),
        rotation_axis, -fd_eps);
    objective_e2 = rotate_about(
        rotate_about(objective_e2, rotation_axis, -fd_eps),
        rotation_axis, -fd_eps);
    myosin.update_endpoints();
    const double myosin_rot_minus = energy();
    myosin.direction[0] = saved_myosin_axis;
    objective_e1 = saved_myosin_axis.cross({0.0, 0.0, 1.0});
    objective_e1.normalize();
    objective_e2 = saved_myosin_axis.cross(objective_e1);
    objective_e2.normalize();
    myosin.update_endpoints();
    assert(near((myosin_rot_plus - myosin_rot_minus) / (2.0 * fd_eps),
                -objective_head.physical_torque_on_myosin.dot(
                    unit_rotation_axis), 3e-7));

    // Nematic director reversal is only a coordinate change for a bipolar
    // thick filament; the physical head wrench must remain identical.
    myosin.direction[0] = -1.0 * saved_myosin_axis;
    const vec flipped_objective_e1 = objective_e1;
    const vec flipped_objective_e2 = -1.0 * objective_e2;
    myosin.update_endpoints();
    const auto flipped_objective_head = compute_am_two_spring_head(
        actin, myosin, 0, 0, box, material_s, -anchor_s, -head_phi,
        head_radius, &flipped_objective_e1, &flipped_objective_e2,
        radial_preferred, radial_k, radial_slack, axial_k,
        stroke_progress, large_stall, true);
    assert((objective_head.force_on_actin -
            flipped_objective_head.force_on_actin).norm() < 1e-10);
    assert((objective_head.physical_torque_on_actin -
            flipped_objective_head.physical_torque_on_actin).norm() < 1e-10);
    assert((objective_head.physical_torque_on_myosin -
            flipped_objective_head.physical_torque_on_myosin).norm() < 1e-10);
    myosin.direction[0] = saved_myosin_axis;
    myosin.update_endpoints();

    // Weak heads use the same objective mechanics but a passive zero-rest
    // longitudinal spring. They must obey the same balance identities.
    const auto weak_objective_head = compute_am_two_spring_head(
        actin, myosin, 0, 0, box, material_s, anchor_s, head_phi,
        head_radius, &objective_e1, &objective_e2,
        radial_preferred, radial_k, radial_slack, 0.2 * axial_k,
        0.0, large_stall, false);
    assert(near(weak_objective_head.active_force, 0.0));
    assert((weak_objective_head.force_on_actin +
            weak_objective_head.force_on_myosin).norm() < 1e-11);
    const vec weak_total_moment =
        (static_cast<vec>(actin.center[0]) - arbitrary_origin).cross(
            weak_objective_head.force_on_actin) +
        weak_objective_head.physical_torque_on_actin +
        (static_cast<vec>(myosin.center[0]) - arbitrary_origin).cross(
            weak_objective_head.force_on_myosin) +
        weak_objective_head.physical_torque_on_myosin;
    assert(weak_total_moment.norm() < 1e-10);

    // A rigid rotation of the complete physical state rotates the resulting
    // wrench without changing its energy (objectivity/covariance).
    const vec covariance_axis{0.2, 0.6, -0.3};
    constexpr double covariance_angle = 0.63;
    actin.center[0] = rotate_about(
        saved_actin_center, covariance_axis, covariance_angle);
    actin.direction[0] = rotate_about(
        saved_actin_axis, covariance_axis, covariance_angle);
    myosin.center[0] = rotate_about(
        saved_myosin_center, covariance_axis, covariance_angle);
    myosin.direction[0] = rotate_about(
        saved_myosin_axis, covariance_axis, covariance_angle);
    const vec rotated_e1 = rotate_about(
        objective_e1, covariance_axis, covariance_angle);
    const vec rotated_e2 = rotate_about(
        objective_e2, covariance_axis, covariance_angle);
    actin.update_endpoints();
    myosin.update_endpoints();
    const auto rotated_objective_head = compute_am_two_spring_head(
        actin, myosin, 0, 0, box, material_s, anchor_s, head_phi,
        head_radius, &rotated_e1, &rotated_e2,
        radial_preferred, radial_k, radial_slack, axial_k,
        stroke_progress, large_stall, true);
    assert(near(rotated_objective_head.energy, objective_head.energy, 1e-10));
    assert((rotated_objective_head.force_on_actin - rotate_about(
                objective_head.force_on_actin, covariance_axis,
                covariance_angle)).norm() < 1e-10);
    assert((rotated_objective_head.physical_torque_on_actin - rotate_about(
                objective_head.physical_torque_on_actin, covariance_axis,
                covariance_angle)).norm() < 1e-10);
    assert((rotated_objective_head.physical_torque_on_myosin - rotate_about(
                objective_head.physical_torque_on_myosin, covariance_axis,
                covariance_angle)).norm() < 1e-10);
    actin.center[0] = saved_actin_center;
    actin.direction[0] = saved_actin_axis;
    myosin.center[0] = saved_myosin_center;
    myosin.direction[0] = saved_myosin_axis;
    actin.update_endpoints();
    myosin.update_endpoints();

    // Minimal stationarity experiment. With no balancing force, an overdamped
    // actin shortens the active element and unloads it. A fixed actin retains
    // isometric force. An equal external load is stationary by force balance.
    constexpr double stiffness = 10.0;
    constexpr double stroke = 0.05;
    constexpr double stall = 1.0;
    constexpr double mobility_dt = 0.01;
    double free_x = 0.0;
    for (int step = 0; step < 2000; ++step) {
        const double force = am_motor_kinetics::active_crossbridge_force(
            stiffness, stroke - free_x, stall);
        free_x += mobility_dt * force;
    }
    const double free_force = am_motor_kinetics::active_crossbridge_force(
        stiffness, stroke - free_x, stall);
    const double fixed_force = am_motor_kinetics::active_crossbridge_force(
        stiffness, stroke, stall);
    assert(free_x > 0.049);
    assert(free_force < 1e-6);
    assert(near(fixed_force, 0.5));

    double balanced_x = 0.0;
    const double external_force = fixed_force;
    for (int step = 0; step < 100; ++step) {
        const double motor_force = am_motor_kinetics::active_crossbridge_force(
            stiffness, stroke - balanced_x, stall);
        balanced_x += mobility_dt * (motor_force - external_force);
    }
    assert(near(balanced_x, 0.0));

    gsl_rng_free(rng);
    return 0;
}
