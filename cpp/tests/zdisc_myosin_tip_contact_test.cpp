#include "interaction.h"

#include <cassert>
#include <cmath>
#include <iostream>
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

    gsl_rng_free(rng);
    std::cout << "zdisc myosin-tip contact tests passed\n";
    return 0;
}
