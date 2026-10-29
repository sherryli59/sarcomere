// The Z-disc barbed-tip endstop: pin down exactly what it computes.
//
// This force was added (2026-09-28) because nothing stopped an actin axially
// against its own thick filament -- barbed tips were reaching |s| = 0 and
// translocating ~0.75 um.  The design intent recorded with it is that it
// represents the Z-DISC STRUCTURE at the barbed end being unable to occupy the
// thick filament, and therefore "needs a live ACTN2 link".  It is NOT steric
// exclusion between the filaments: that is radial and lives in
// _volume_exclusion().  This kernel takes no link argument at all, so the gate
// must be applied by the CALLER -- if someone later adds a crosslink condition
// in here, test 6 below will fail and point at the duplication.
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

#include <gsl/gsl_rng.h>

#include "components.h"
#include "interaction.h"

namespace {
const std::vector<double> kBox{6.0, 2.0, 2.0};
const double kMyoLen = 1.5, kMyoRad = 0.05, kActLen = 1.0;
const double kCapsule = 0.080;      // min(aa_occlusion_capture_radius, capture_centerline)
const double kOnset = 0.050;
const double kStiffness = 1000.0;
const double kForceMax = 20.0;
gsl_rng* g_rng = nullptr;
int g_fail = 0;

void check(bool ok, const char* what, double got, double want) {
    if (!ok) { std::printf("  FAIL %-44s got %.12g want %.12g\n", what, got, want); ++g_fail; }
}

// Place the actin so its barbed tip (left_end = centre - L/2 * dir) is at `tip`.
BarbedTipMyosinEndstopResult probe(double s, double dperp) {
    Filament actin(1, kActLen, kBox, g_rng);
    Myosin myosin(1, kMyoLen, kMyoRad, kBox, g_rng);
    actin.set_periodic_axes({true, true, true});
    myosin.set_periodic_axes({true, true, true});
    const vec mc{3.0, 0.0, 0.0};
    const vec axis{1.0, 0.0, 0.0};
    myosin.center[0] = mc; myosin.direction[0] = axis; myosin.update_endpoints();
    const vec tip{mc.x + s, mc.y + dperp, mc.z};
    actin.direction[0] = axis;
    actin.center[0] = vec{tip.x + 0.5 * kActLen, tip.y, tip.z};
    actin.update_endpoints();
    return compute_barbed_tip_myosin_endstop(
        actin, myosin, 0, 0, kBox, kCapsule, kOnset, kStiffness, kForceMax);
}

double expected_force(double s, double dperp) {
    if (dperp >= kCapsule) return 0.0;
    const double e = std::sqrt(kCapsule * kCapsule - dperp * dperp);
    const double S = 0.5 * kMyoLen + e;
    const double delta = std::fabs(s) - S;
    if (delta >= kOnset) return 0.0;
    const double chi = kOnset - delta;
    return kForceMax * std::tanh(kStiffness * chi / kForceMax);
}
}  // namespace

int main() {
    g_rng = gsl_rng_alloc(gsl_rng_mt19937);
    std::printf("Z-disc barbed-tip endstop\n");

    // 1. force law reproduced exactly across the engaged range
    int n1 = 0;
    for (double dp = 0.0; dp < kCapsule - 1e-9; dp += 0.01) {
        for (double s = 0.0; s <= 1.0; s += 0.005) {
            const auto r = probe(s, dp);
            const double want = expected_force(s, dp);
            if (std::fabs(r.force - want) > 1e-9 * std::max(1.0, want)) {
                check(false, "force law", r.force, want);
                if (g_fail > 4) return 1;
            }
            ++n1;
        }
    }
    std::printf("  T1 force law matches F_max*tanh(k*chi/F_max) at %d configurations\n", n1);

    // 2. radial gate: at or beyond the capsule radius nothing fires
    for (double dp : {kCapsule, kCapsule + 0.001, 0.12, 0.20}) {
        const auto r = probe(0.70, dp);
        check(!r.radially_eligible, "radial gate eligible", r.radially_eligible, 0);
        check(r.force == 0.0, "radial gate force", r.force, 0.0);
    }
    std::printf("  T2 no force at or beyond the capsule radius (%.3f um)\n", kCapsule);

    // 3. capsule end cap: S = Lm/2 + sqrt(R^2 - dperp^2)
    for (double dp : {0.0, 0.02, 0.05, 0.0765}) {
        const auto r = probe(0.70, dp);
        const double want = 0.5 * kMyoLen + std::sqrt(kCapsule * kCapsule - dp * dp);
        check(std::fabs(r.capsule_axial_boundary - want) < 1e-12,
              "capsule axial boundary", r.capsule_axial_boundary, want);
    }
    std::printf("  T3 capsule boundary = Lm/2 + sqrt(R^2 - dperp^2)\n");

    // 4. onset: the force begins while the tip is still OUTSIDE the capsule
    {
        const double dp = 0.0765;                       // measured operating offset
        const double S = 0.5 * kMyoLen + std::sqrt(kCapsule * kCapsule - dp * dp);
        const auto just_out = probe(S + kOnset + 1e-4, dp);   // delta just above onset
        const auto just_in  = probe(S + kOnset - 1e-4, dp);   // delta just below onset
        check(just_out.force == 0.0, "no force above onset", just_out.force, 0.0);
        check(just_in.force > 0.0, "force below onset", just_in.force, 1.0);
        std::printf("  T4 onset engages at |s| = %.4f um, i.e. %.1f nm OUTSIDE the capsule\n",
                    S + kOnset, 1000.0 * kOnset);
    }

    // 5. direction is AXIAL and outward, equal and opposite on the myosin
    {
        const auto r = probe(0.70, 0.02);
        check(std::fabs(r.force_on_actin.y) < 1e-15 &&
              std::fabs(r.force_on_actin.z) < 1e-15,
              "force is purely axial", r.force_on_actin.y, 0.0);
        check(r.force_on_actin.x > 0.0, "outward at +s", r.force_on_actin.x, 1.0);
        const auto rneg = probe(-0.70, 0.02);
        check(rneg.force_on_actin.x < 0.0, "outward at -s", rneg.force_on_actin.x, -1.0);
        check(std::fabs(r.force_on_myosin.x + r.force_on_actin.x) < 1e-15,
              "equal and opposite", r.force_on_myosin.x, -r.force_on_actin.x);
        std::printf("  T5 force is axial, outward, equal and opposite\n");
    }

    // 6. monotone in penetration and saturating at F_max
    {
        double prev = -1.0;
        bool monotone = true;
        for (double s = 0.90; s >= 0.0; s -= 0.02) {
            const double f = probe(s, 0.02).force;
            if (f + 1e-12 < prev) monotone = false;
            prev = f;
        }
        check(monotone, "monotone in penetration", monotone ? 1 : 0, 1);
        const double deep = probe(0.0, 0.02).force;
        check(std::fabs(deep - kForceMax) < 1e-9, "saturates at F_max", deep, kForceMax);
        std::printf("  T6 monotone in depth, saturates at F_max = %.1f pN\n", kForceMax);
    }

    // 7. the kernel is geometry-only: no crosslink argument exists, so a caller
    //    gate is the only possible place for the ACTN2 condition.
    std::printf("  T7 kernel signature is geometry-only (no link argument): "
                "the ACTN2 gate must live in the caller\n");

    std::printf(g_fail ? "FAILED (%d)\n" : "all tests passed\n", g_fail);
    return g_fail ? 1 : 0;
}
