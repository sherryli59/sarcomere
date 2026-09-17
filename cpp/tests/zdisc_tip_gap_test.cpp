// geometry::-free unit test: the cheap tip predicate must agree exactly with
// the full contact kernel it short-circuits.
//
// resolve_actn2_tip_hard_contacts evaluates the full kernel ~53,000 times per
// step -- smoothstep, tanh, an exp/log1p/log energy and three orientation
// drives -- and then discards all of it unless axial_gap < -tolerance, which
// is true about twice per step.  zdisc_myosin_tip_gap replaces that decision.
// It is only admissible if, for every configuration:
//   (a) axial_gap and transverse_distance match bit for bit, and
//   (b) `eligible && axial_gap < -tol`  <=>  `applied && axial_gap < -tol`.
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>

#include <gsl/gsl_rng.h>

#include "components.h"
#include "interaction.h"

namespace {

const std::vector<double> kBox{6.0, 2.0, 2.0};
const double kOnset = 0.050;
const double kStiffness = 1000.0;
const double kForceMax = 20.0;
const double kFlat = 0.100;
const double kCutoff = 0.130;
const double kTolerance = 1e-9;
gsl_rng* g_rng = nullptr;

bool same_bits(double a, double b) {
    if (std::isnan(a) && std::isnan(b)) return true;
    std::uint64_t x, y;
    std::memcpy(&x, &a, sizeof x);
    std::memcpy(&y, &b, sizeof y);
    return x == y;
}

}  // namespace

int main() {
    g_rng = gsl_rng_alloc(gsl_rng_mt19937);
    std::mt19937_64 rng(20260901u);
    std::uniform_real_distribution<double> u01(0.0, 1.0);

    long long tested = 0, eligible = 0, would_project = 0, mismatch = 0;
    long long decision_mismatch = 0;

    for (int trial = 0; trial < 300000; ++trial) {
        Filament actin(2, 1.0, kBox, g_rng);
        actin.set_periodic_axes({true, true, true});
        Myosin myosin(1, 1.5, 0.05, kBox, g_rng);
        myosin.set_periodic_axes({true, true, true});

        // Put the pair near a myosin tip often enough to exercise the branch.
        const vec mc{u01(rng) * kBox[0], (u01(rng) - 0.5) * 1.5,
                     (u01(rng) - 0.5) * 1.5};
        vec md{1.0, (u01(rng) - 0.5) * 0.5, (u01(rng) - 0.5) * 0.5};
        md.normalize();
        myosin.center[0] = mc;
        myosin.direction[0] = md;
        myosin.update_endpoints();

        const int end_sign = u01(rng) < 0.5 ? -1 : 1;
        const vec tip = mc + static_cast<double>(end_sign) * 0.75 * md;
        for (int k = 0; k < 2; ++k) {
            vec d{u01(rng) < 0.5 ? 1.0 : -1.0, (u01(rng) - 0.5) * 0.4,
                  (u01(rng) - 0.5) * 0.4};
            d.normalize();
            actin.direction[k] = d;
            // scatter around the tip, sometimes well inside the filament
            const vec jitter{(u01(rng) - 0.5) * 1.8, (u01(rng) - 0.5) * 0.4,
                             (u01(rng) - 0.5) * 0.4};
            actin.center[k] = tip + jitter;
        }
        actin.update_endpoints();
        const double s_i = (u01(rng) - 0.5) * 1.0;
        const double s_j = (u01(rng) - 0.5) * 1.0;

        const auto gap = zdisc_myosin_tip_gap(
            actin, myosin, 0, 1, s_i, s_j, 0, end_sign, kBox, kOnset, kCutoff);
        const auto full = compute_zdisc_myosin_tip_contact(
            actin, myosin, 0, 1, s_i, s_j, 0, end_sign, kBox, kOnset,
            kStiffness, kForceMax, kFlat, kCutoff);

        ++tested;
        if (gap.eligible) ++eligible;

        // (a) geometry agrees bit for bit wherever the full kernel got that far
        if (full.eligible_half) {
            if (!same_bits(gap.axial_gap, full.axial_gap) ||
                !same_bits(gap.transverse_distance, full.transverse_distance)) {
                ++mismatch;
            }
        }
        // (b) the decision the solver actually makes is identical
        const bool cheap = gap.eligible && gap.axial_gap < -kTolerance;
        const bool exact = full.applied && full.axial_gap < -kTolerance;
        if (cheap != exact) ++decision_mismatch;
        if (exact) ++would_project;
    }

    std::printf("zdisc_myosin_tip_gap vs compute_zdisc_myosin_tip_contact\n");
    std::printf("  tested %lld, eligible %lld, would-project %lld\n",
                tested, eligible, would_project);
    std::printf("  geometry bit mismatches   %lld (must be 0)\n", mismatch);
    std::printf("  decision disagreements    %lld (must be 0)\n",
                decision_mismatch);
    assert(mismatch == 0);
    assert(decision_mismatch == 0);
    assert(would_project > 200);    // the projecting branch must be exercised
    assert(eligible > 2000);
    std::printf("zdisc_tip_gap_test: all checks passed\n");
    return 0;
}
