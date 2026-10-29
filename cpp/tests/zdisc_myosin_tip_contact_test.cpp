#include "interaction.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include <gsl/gsl_rng.h>

namespace {
bool close(double a, double b, double tol = 1e-9) {
    return std::abs(a - b) <= tol;
}
}

int main() {
    const std::vector<double> box{8.0, 8.0, 8.0};
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    Filament actin(2, 0.50, box, rng);
    Myosin myosin(1, 1.00, 0.05, box, rng);
    actin.center[0] = {0.54, -0.08, 0.0};
    actin.center[1] = {0.54,  0.08, 0.0};
    actin.direction[0] = {0.0, 1.0, 0.0};
    actin.direction[1] = {0.0, 1.0, 0.0};
    myosin.center[0] = {0.0, 0.0, 0.0};
    myosin.direction[0] = {1.0, 0.0, 0.0};
    actin.update_endpoints();
    myosin.update_endpoints();

    const auto full = compute_zdisc_myosin_tip_contact(
        actin, myosin, 0, 1, 0.08, -0.08, 0, 1, box,
        0.05, 1000.0, 20.0, 0.10, 0.13);
    assert(full.applied);
    assert(close(full.axial_gap, 0.04));
    assert(close(full.transverse_distance, 0.0));
    assert(close(full.radial_gate, 1.0));
    assert(close(full.force, 20.0 * std::tanh(0.5)));
    assert(full.force_on_actin_i.x > 0.0);
    assert(full.force_on_myosin.x < 0.0);
    const vec net = full.force_on_actin_i + full.force_on_actin_j +
        full.force_on_myosin;
    assert(net.norm() < 1e-12);

    // The opposite tip cannot claim a patch on the right half.
    const auto wrong_half = compute_zdisc_myosin_tip_contact(
        actin, myosin, 0, 1, 0.08, -0.08, 0, -1, box,
        0.05, 1000.0, 20.0, 0.10, 0.13);
    assert(!wrong_half.applied);

    // Move the patch outside the axial onset.
    actin.center[0].x = actin.center[1].x = 0.56;
    actin.update_endpoints();
    const auto axial_clear = compute_zdisc_myosin_tip_contact(
        actin, myosin, 0, 1, 0.08, -0.08, 0, 1, box,
        0.05, 1000.0, 20.0, 0.10, 0.13);
    assert(!axial_clear.applied);

    // At the midpoint of the transverse taper, the cubic gate is one half.
    actin.center[0] = {0.54, 0.035, 0.0};
    actin.center[1] = {0.54, 0.195, 0.0};
    actin.update_endpoints();
    const auto tapered = compute_zdisc_myosin_tip_contact(
        actin, myosin, 0, 1, 0.08, -0.08, 0, 1, box,
        0.05, 1000.0, 20.0, 0.10, 0.13);
    assert(tapered.applied);
    assert(close(tapered.transverse_distance, 0.115));
    assert(close(tapered.radial_gate, 0.5));
    assert(close(tapered.force, 0.5 * full.force));

    actin.center[0].y += 0.03;
    actin.center[1].y += 0.03;
    actin.update_endpoints();
    const auto radial_clear = compute_zdisc_myosin_tip_contact(
        actin, myosin, 0, 1, 0.08, -0.08, 0, 1, box,
        0.05, 1000.0, 20.0, 0.10, 0.13);
    assert(!radial_clear.applied);

    actin.set_periodic_axes({false, false, false});
    myosin.set_periodic_axes({false, false, false});
    auto contact_at = [&](vec p, vec q, int sign = 1) {
        actin.center[0] = p;
        actin.center[1] = q;
        actin.update_endpoints();
        return compute_zdisc_myosin_tip_contact(
            actin, myosin, 0, 1, 0, 0, 0, sign, box,
            .05, 1000, 20, .10, .13);
    };
    // Midpoint at 140 nm; the 120 nm endpoint still requests repulsion.
    auto endpoint = contact_at({.48, .120, 0}, {.48, .160, 0});
    assert(endpoint.applied);
    assert(close(endpoint.transverse_distance, .120));
    assert(close(endpoint.axial_gap, -.02));
    assert(endpoint.force > 0 && endpoint.force < 20);
    assert(close(endpoint.force_on_actin_i.x, .5 * endpoint.force));
    assert(close(endpoint.force_on_actin_j.x, .5 * endpoint.force));
    assert(endpoint.force_on_actin_i.y == 0); // remains axial only
    const auto reverse = contact_at({.48, .160, 0}, {.48, .120, 0});
    assert(close(reverse.force, endpoint.force));
    // Neither endpoint is close; the connector interior crosses the axis.
    assert(contact_at({.48, -.2, 0}, {.48, .2, 0}).applied);
    assert(!contact_at({.48, .131, 0}, {.48, .16, 0}).applied);
    assert(!contact_at({.48, .130, 0}, {.48, .16, 0}).applied);
    // Radial and axial eligibility at DIFFERENT points must not combine.
    assert(!contact_at({.7, .05, 0}, {.5, .30, 0}).applied);
    // Minimum radius lies on onset: some interior points are still eligible.
    assert(contact_at({.7, .05, 0}, {.5, .15, 0}).applied);
    assert(!contact_at({.55, .1, 0}, {.55, .12, 0}).applied);
    assert(contact_at({.48, .1, 0}, {.48, .1, 0}).applied);
    // One endpoint can be on the opposite half; the eligible part survives.
    assert(contact_at({-.3, .12, 0}, {.1, .12, 0}).applied);
    assert(contact_at({-.48, .12, 0}, {-.48, .16, 0}, -1).applied);
    // The old midpoint-only sphere would reject this endpoint contact.
    assert(contact_at({.54, .12, 0}, {.9, .4, 0}).applied);
    auto patch = compute_actn2_patch_geometry(actin, 0, 1, 0, 0, box);
    assert(actn2_patch_may_contact_myosin_tip(patch, myosin, 0, box, .05, .13));
    myosin.direction[0] = {0, 1, 0};
    myosin.update_endpoints();
    assert(contact_at({.12, .48, 0}, {.16, .48, 0}).applied);
    myosin.direction[0] = {1, 0, 0};
    myosin.center[0] = {0, 3.9, 0};
    myosin.update_endpoints();
    actin.set_periodic_axes({true, true, true});
    myosin.set_periodic_axes({true, true, true});
    const auto periodic = contact_at({.48, -3.98, 0}, {.48, -3.94, 0});
    assert(periodic.applied);
    assert(close(periodic.transverse_distance, .12));
    const auto translated = contact_at({8.48, 4.02, 0}, {8.48, 4.06, 0});
    assert(close(translated.force, periodic.force));
    {
        // A tilted long myosin can contact a non-nearest periodic image.
        const double s = std::sqrt(.5);
        Myosin tilted(1, 1.5, .05, box, rng);
        tilted.center[0] = {0, 0, 0};
        tilted.direction[0] = {s, s, 0};
        tilted.update_endpoints();
        actin.center[0] = {.48, -.48, 0};
        actin.center[1] = {.49, -.48, 0};
        actin.update_endpoints();
        const auto image = compute_zdisc_myosin_tip_contact(
            actin, tilted, 0, 1, 0, 0, 0, -1, {1, 2, 2},
            .05, 1000, 20, .10, .13);
        assert(image.applied);
    }
    // Independent dense sampling: every sampled eligible point must trigger
    // the continuous segment kernel, even when endpoint lists are empty.
    actin.set_periodic_axes({false, false, false});
    myosin.set_periodic_axes({false, false, false});
    myosin.center[0] = {0, 0, 0};
    myosin.update_endpoints();
    // Measured original overshoot geometry: r27zt0np, frame 475,
    // link 179352 / myosin 244, in a local myosin-axis basis.
    assert(contact_at(
        {-.003812085957410629, .11732520591698475, .03345072485325059},
        {.014175370169763528, .09717262963312846, -.03345072485325024}).applied);
    std::mt19937_64 random(20261028);
    std::uniform_real_distribution<double> axial(-.8, .8), radial(-.3, .3);
    for (int trial = 0; trial < 3000; ++trial) {
        vec p{axial(random), radial(random), radial(random)};
        vec q{axial(random), radial(random), radial(random)};
        const auto candidate = contact_at(p, q);
        bool sampled = false;
        for (int k = 0; k <= 200; ++k) {
            vec point = p + (k / 200.0) * (q - p);
            sampled |= point.x >= 0 && point.x < .55 &&
                std::hypot(point.y, point.z) < .129999;
        }
        if (sampled) assert(candidate.applied);
        if (candidate.applied) {
            assert(candidate.transverse_distance < .13);
            assert(candidate.axial_gap < .05);
            assert(candidate.axial_gap >= -.5 - 1e-12);
            const auto bounds = compute_actn2_patch_geometry(
                actin, 0, 1, 0, 0, box);
            assert(actn2_patch_may_contact_myosin_tip(
                bounds, myosin, 0, box, .05, .13));
        }
    }

    // Independent barbed-tip end-stop. With a 1-um myosin and R=0.13,
    // the on-axis capsule end is x=0.63. The actin's left endpoint is its
    // material barbed tip.
    Filament tip_actin(1, 0.50, box, rng);
    tip_actin.direction[0] = {1.0, 0.0, 0.0};
    tip_actin.center[0] = {0.90, 0.0, 0.0};  // barbed x=0.65, g=+0.02
    tip_actin.update_endpoints();
    myosin.center[0] = {0.0, 0.0, 0.0};
    myosin.direction[0] = {1.0, 0.0, 0.0};
    myosin.update_endpoints();
    const auto approach = compute_barbed_tip_myosin_endstop(
        tip_actin, myosin, 0, 0, box, 0.13, 0.05, 1000.0, 20.0);
    assert(approach.applied);
    assert(close(approach.signed_clearance, 0.02));
    assert(close(approach.penetration, 0.0));
    assert(approach.force_on_actin.x > 0.0);
    assert((approach.force_on_actin + approach.force_on_myosin).norm() < 1e-12);

    tip_actin.center[0].x = 0.87;  // barbed x=0.62, g=-0.01
    tip_actin.update_endpoints();
    const auto penetrated = compute_barbed_tip_myosin_endstop(
        tip_actin, myosin, 0, 0, box, 0.13, 0.05, 1000.0, 20.0);
    assert(penetrated.applied);
    assert(close(penetrated.signed_clearance, -0.01));
    assert(close(penetrated.penetration, 0.01));
    assert(penetrated.force > approach.force);

    // The production steric radius is 0.080 um.  At the same barbed-tip
    // position, the 0.130-um head-reach capsule reports penetration whereas
    // the steric capsule retains 0.040 um of positive axial clearance.
    const auto steric_clear = compute_barbed_tip_myosin_endstop(
        tip_actin, myosin, 0, 0, box, 0.080, 0.05, 1000.0, 20.0);
    assert(steric_clear.applied);
    assert(close(steric_clear.signed_clearance, 0.04));
    assert(close(steric_clear.penetration, 0.0));

    // The same axial coordinate is allowed outside the capsule radius: the
    // rule is axial-dominant but remains a finite capsule, not an infinite wall.
    tip_actin.center[0] = {0.87, 0.14, 0.0};
    tip_actin.update_endpoints();
    const auto outside_radius = compute_barbed_tip_myosin_endstop(
        tip_actin, myosin, 0, 0, box, 0.13, 0.05, 1000.0, 20.0);
    assert(!outside_radius.radially_eligible);
    assert(!outside_radius.applied);

    gsl_rng_free(rng);
    std::cout << "zdisc myosin-tip contact tests passed\n";
    return 0;
}
