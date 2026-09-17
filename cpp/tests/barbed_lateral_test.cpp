// Unit test for the barbed-end lateral repulsion.
//
// The term exists because ACTN2 is a purely attractive crosslinker: a
// crosslinked cluster contracts until barbed ends reach the actin hard core
// (measured 41 nm -> 26 nm as links accumulate), while tiling a Z-disc face at
// this filament count needs ~81 nm.  These checks pin the properties that
// matter: it is inert when unconfigured, it pushes transversely only, it is
// reciprocal, it respects the axial gate and the force cap, and it does not
// invent a direction for coincident ends.
#include <cassert>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include <gsl/gsl_rng.h>

#include "components.h"
#include "interaction.h"

namespace {

const std::vector<double> kBox{6.0, 2.0, 2.0};
gsl_rng* g_rng = nullptr;

Filament make_pair(double separation_y, double axial_offset,
                   bool antiparallel) {
    Filament actin(2, 1.0, kBox, g_rng);
    actin.set_periodic_axes({true, true, true});
    // direction points barbed -> pointed, so barbed = centre - L/2 * direction
    actin.direction[0] = vec{1.0, 0.0, 0.0};
    actin.direction[1] = vec{antiparallel ? -1.0 : 1.0, 0.0, 0.0};
    // place barbed end 0 at the origin
    actin.center[0] = vec{0.5, 0.0, 0.0};
    // barbed end 1 offset transversely by separation_y and axially by offset
    const vec barbed1{axial_offset, separation_y, 0.0};
    actin.center[1] = barbed1 + 0.5 * actin.length * actin.direction[1];
    actin.update_endpoints();
    return actin;
}

bool nearly(double a, double b, double tol = 1e-9) {
    return std::abs(a - b) < tol;
}


// _volume_exclusion now skips this kernel when the barbed ends are further
// apart than sqrt(axial_window^2 + range^2).  That is only sound if the kernel
// itself never applies a force beyond that separation, which is what this
// pins: the axial and transverse gates constrain orthogonal components of the
// same vector, so the bound is necessary, not merely usually true.
void reach_bound_is_necessary() {
    const double range = 0.080;
    const double window = 0.250;
    const double reach_squared = window * window + range * range;
    std::mt19937_64 rng(20260904u);
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    long long applied = 0, violations = 0;

    for (int trial = 0; trial < 200000; ++trial) {
        Filament actin(2, 1.0, kBox, g_rng);
        actin.set_periodic_axes({true, true, true});
        for (int k = 0; k < 2; ++k) {
            vec dir{k == 0 ? 1.0 : (u01(rng) < 0.5 ? 1.0 : -1.0),
                    (u01(rng) - 0.5) * 0.8, (u01(rng) - 0.5) * 0.8};
            dir.normalize();
            actin.direction[k] = dir;
        }
        const vec first_center{u01(rng) * 0.9, (u01(rng) - 0.5) * 0.9,
                               (u01(rng) - 0.5) * 0.9};
        // Draw the partner within reach often enough that the kernel fires.
        const vec offset{(u01(rng) - 0.5) * 0.8, (u01(rng) - 0.5) * 0.3,
                         (u01(rng) - 0.5) * 0.3};
        actin.center[0] = first_center;
        actin.center[1] = first_center + offset;
        actin.update_endpoints();

        const auto result = compute_barbed_end_lateral_repulsion(
            actin, 0, 1, kBox, range, 2000.0, 40.0, window);
        if (!result.applied) continue;
        ++applied;
        vec u0 = actin.direction[0], u1 = actin.direction[1];
        u0.normalize(); u1.normalize();
        vec separation = (static_cast<vec>(actin.center[1]) - 0.5 * actin.length * u1) -
                         (static_cast<vec>(actin.center[0]) - 0.5 * actin.length * u0);
        separation.pbc_wrap(kBox, {true, true, true});
        if (!(separation.norm_squared() < reach_squared)) ++violations;
    }
    printf("  reach_bound: %lld applications, %lld outside the gate (must be 0)\n",
           applied, violations);
    assert(applied > 500);
    assert(violations == 0);
}

}  // namespace

int main() {
    g_rng = gsl_rng_alloc(gsl_rng_mt19937);
    reach_bound_is_necessary();
    const double range = 0.080;
    const double k = 500.0;

    // 1. inert when unconfigured - this is what keeps existing runs identical
    {
        Filament actin = make_pair(0.030, 0.0, true);
        assert(!compute_barbed_end_lateral_repulsion(
                   actin, 0, 1, kBox, 0.0, k, 0.0, 0.25).applied);
        assert(!compute_barbed_end_lateral_repulsion(
                   actin, 0, 1, kBox, range, 0.0, 0.0, 0.25).applied);
    }

    // 2. compressed pair pushes apart, transversely, reciprocally
    {
        Filament actin = make_pair(0.030, 0.0, true);
        const auto r = compute_barbed_end_lateral_repulsion(
            actin, 0, 1, kBox, range, k, 0.0, 0.25);
        assert(r.applied);
        assert(nearly(r.transverse_distance, 0.030));
        const double expected = k * (range - 0.030);
        assert(nearly(r.force_on_first.norm(), expected, 1e-6));
        // pushes first away from second: second is at +y, so force is -y
        assert(r.force_on_first.y < 0.0);
        assert(r.force_on_second.y > 0.0);
        // reciprocal: no net momentum
        const vec net = r.force_on_first + r.force_on_second;
        assert(nearly(net.norm(), 0.0, 1e-12));
        // transverse only: nothing along the shared x axis
        assert(nearly(r.force_on_first.x, 0.0, 1e-12));
    }

    // 3. beyond range is inert, and force grows as spacing shrinks
    {
        Filament far = make_pair(range * 1.01, 0.0, true);
        assert(!compute_barbed_end_lateral_repulsion(
                   far, 0, 1, kBox, range, k, 0.0, 0.25).applied);
        Filament near = make_pair(0.020, 0.0, true);
        Filament mid = make_pair(0.060, 0.0, true);
        const double f_near = compute_barbed_end_lateral_repulsion(
            near, 0, 1, kBox, range, k, 0.0, 0.25).force_on_first.norm();
        const double f_mid = compute_barbed_end_lateral_repulsion(
            mid, 0, 1, kBox, range, k, 0.0, 0.25).force_on_first.norm();
        assert(f_near > f_mid);
    }

    // 4. axial gate: ends in different Z-discs must not push on each other
    {
        Filament inside = make_pair(0.030, 0.20, true);
        Filament outside = make_pair(0.030, 0.30, true);
        assert(compute_barbed_end_lateral_repulsion(
                   inside, 0, 1, kBox, range, k, 0.0, 0.25).applied);
        assert(!compute_barbed_end_lateral_repulsion(
                   outside, 0, 1, kBox, range, k, 0.0, 0.25).applied);
    }

    // 5. force cap saturates
    {
        Filament actin = make_pair(0.005, 0.0, true);
        const double cap = 3.0;
        const auto r = compute_barbed_end_lateral_repulsion(
            actin, 0, 1, kBox, range, k, cap, 0.25);
        assert(r.applied);
        assert(nearly(r.force_on_first.norm(), cap, 1e-9));
    }

    // 6. coincident ends have no defined push direction - stay out of it
    {
        Filament actin = make_pair(0.0, 0.0, true);
        assert(!compute_barbed_end_lateral_repulsion(
                   actin, 0, 1, kBox, range, k, 0.0, 0.25).applied);
    }

    // 7. parallel pairs are included too: the collapse is driven by unlinked
    //    neighbours, which have no polarity requirement
    {
        Filament actin = make_pair(0.030, 0.0, false);
        assert(compute_barbed_end_lateral_repulsion(
                   actin, 0, 1, kBox, range, k, 0.0, 0.25).applied);
    }

    std::printf("barbed_end_lateral_repulsion: all checks passed\n");
    gsl_rng_free(g_rng);
    return 0;
}
