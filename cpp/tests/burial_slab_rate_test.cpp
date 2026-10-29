// The three burial fixes, pinned down.
//
// Context.  Reason-3 "burial" deletes an ACTN2 crosslink with probability 1
// when a barbed tip penetrates the myosin capture capsule by
// aa_occlusion_rupture_penetration.  MEASURED in production: 65.3% of
// 5,963,901 ACTN2 removals, and in the per-step rig the repulsion force was
// acting on 0.0% of the 600 us before every burial deletion -- the tip went
// from penetration 0 to past the kill depth in a single step with the force
// arriving at full 20 pN in that same step.  Three causes, three fixes:
//
//   1. The transverse HARD CUTOFF.  In production capsule_radius (0.080) is
//      exactly the resting A-M centerline, so a tip outside it is invisible
//      and meets no opposition at all; when it crosses in it is already deep.
//   2. The CAP EXTENSION makes the axial boundary a function of the
//      transverse distance, so a purely axial force acquires a curl.
//   3. A probability-1 veto reading a SOFT force's thermal tail is
//      unwinnable: the spring holds a distribution of depths, and the first
//      excursion past the line is fatal however stiff the spring is.
//
// Tests 1 and 4 are the inertness guards: with the new flags off, every
// number must be the legacy one.
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

#include <gsl/gsl_rng.h>

#include "catch_bond_kinetics.h"
#include "components.h"
#include "interaction.h"

namespace {
const std::vector<double> kBox{6.0, 2.0, 2.0};
const double kMyoLen = 1.5, kMyoRad = 0.05, kActLen = 1.0;
const double kCapsule = 0.080;   // = myosin_radius + actin_radius + am_surface_gap
const double kOnset = 0.050;
const double kStiffness = 1000.0;
const double kForceMax = 20.0;
const double kWideCutoff = 0.250;
const double kKillOld = 0.010;
const double kKillNew = 0.050;
gsl_rng* g_rng = nullptr;
int g_fail = 0;

void check(bool ok, const char* what, double got, double want) {
    if (!ok) {
        std::printf("  FAIL %-50s got %.12g want %.12g\n", what, got, want);
        ++g_fail;
    }
}

BarbedTipMyosinEndstopResult probe(double s, double dperp, bool slab,
                                   double cutoff) {
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
        actin, myosin, 0, 0, kBox, kCapsule, kOnset, kStiffness, kForceMax,
        slab, cutoff);
}

// The legacy force law, written out independently of the implementation.
double legacy_force(double s, double dperp) {
    if (dperp >= kCapsule) return 0.0;
    const double e = std::sqrt(kCapsule * kCapsule - dperp * dperp);
    const double delta = std::fabs(s) - (0.5 * kMyoLen + e);
    if (delta >= kOnset) return 0.0;
    return kForceMax * std::tanh(kStiffness * (kOnset - delta) / kForceMax);
}
}  // namespace

int main() {
    g_rng = gsl_rng_alloc(gsl_rng_mt19937);
    std::printf("Burial fixes: slab geometry + hazard rate\n");

    // --- 1. DEFAULT INERT: slab=false reproduces the legacy law exactly -----
    int n1 = 0;
    for (double dp = 0.0; dp < 0.120; dp += 0.005) {
        for (double s = 0.60; s < 0.95; s += 0.005) {
            const auto r = probe(s, dp, false, 0.0);
            const double want = legacy_force(s, dp);
            check(std::fabs(r.force - want) < 1e-12,
                  "legacy force law unchanged", r.force, want);
            ++n1;
        }
    }
    std::printf("1. default (slab=false) reproduces the legacy law: %d points\n", n1);

    // --- 2. THE BLIND ANNULUS, reproduced and then closed -------------------
    // A tip sitting beside the barrel, 300 nm inside the A-band axially.
    const double s_in = 0.45;             // |s| < 0.75, i.e. inside the slab
    const auto out_legacy = probe(s_in, 0.090, false, 0.0);   // just outside the gate
    const auto in_legacy  = probe(s_in, 0.079, false, 0.0);   // just inside it
    check(!out_legacy.radially_eligible, "legacy: tip at 90 nm is invisible",
          out_legacy.radially_eligible ? 1 : 0, 0);
    check(out_legacy.force == 0.0, "legacy: no force at 90 nm", out_legacy.force, 0.0);
    check(in_legacy.penetration > kKillNew,
          "legacy: one step later it is already past the kill line",
          in_legacy.penetration, kKillNew);
    std::printf("2. LEGACY at s=%.3f: dperp 0.090 -> eligible=%d force=%.2f pN | "
                "dperp 0.079 -> penetration=%.1f nm force=%.2f pN\n",
                s_in, out_legacy.radially_eligible ? 1 : 0, out_legacy.force,
                1000.0 * in_legacy.penetration, in_legacy.force);
    const auto out_slab = probe(s_in, 0.090, true, kWideCutoff);
    check(out_slab.radially_eligible, "slab: tip at 90 nm is now seen",
          out_slab.radially_eligible ? 1 : 0, 1);
    check(out_slab.force > 0.0, "slab: force acts at 90 nm", out_slab.force, 1.0);
    std::printf("   SLAB   at s=%.3f: dperp 0.090 -> eligible=%d force=%.2f pN "
                "penetration=%.1f nm\n", s_in, out_slab.radially_eligible ? 1 : 0,
                out_slab.force, 1000.0 * out_slab.penetration);

    // --- 3. SLAB BOUNDARY IS INDEPENDENT OF THE TRANSVERSE DISTANCE --------
    // This is what removes the curl: F_z may not depend on r.
    int n3 = 0;
    for (double s = 0.60; s < 0.80; s += 0.01) {
        const double want = std::max(0.0, 0.5 * kMyoLen - std::fabs(s));
        double f0 = -1.0;
        for (double dp : {0.005, 0.040, 0.079, 0.120, 0.200}) {
            const auto r = probe(s, dp, true, kWideCutoff);
            check(std::fabs(r.penetration - want) < 1e-12,
                  "slab penetration = L/2 - |s|, independent of r",
                  r.penetration, want);
            if (f0 < 0.0) f0 = r.force;
            check(std::fabs(r.force - f0) < 1e-12,
                  "slab force independent of r (zero curl)", r.force, f0);
            ++n3;
        }
    }
    std::printf("3. slab boundary and force independent of r: %d points\n", n3);

    // --- 4. THE CURL, measured as closed-loop work -------------------------
    // F has only a z component, so a rectangle in (z, r) has zero work on its
    // two radial legs and W = integral[s1..s2] (F_z(s,r1) - F_z(s,r2)) ds.
    // That is nonzero exactly when F_z depends on r, i.e. exactly when the
    // cap extension is in play.  Pick s where the tanh is NOT saturated at
    // both radii -- deep inside the capsule both legs sit at force_max and the
    // defect is invisible, which is what an earlier version of this test
    // sampled.  Reported in kT at 310 K (4.28 pN*nm) because the question is
    // whether thermal cycling can pump it.
    auto loop_work = [](bool slab, double cut) {
        const double r1 = 0.005, r2 = 0.079;
        const double s1 = 0.78, s2 = 0.88;
        const int N = 20000;
        const double ds = (s2 - s1) / N;
        double w = 0.0;
        for (int k = 0; k < N; ++k) {
            const double s = s1 + (k + 0.5) * ds;
            w += (probe(s, r1, slab, cut).force -
                  probe(s, r2, slab, cut).force) * ds;
        }
        return w;
    };
    const double kT = 0.00428;  // pN*um at 310 K
    const double w_legacy = loop_work(false, 0.0);
    const double w_slab = loop_work(true, kWideCutoff);
    check(std::fabs(w_legacy) > 10.0 * kT,
          "legacy endstop is not a gradient (closed-loop work != 0)",
          std::fabs(w_legacy) / kT, 10.0);
    check(std::fabs(w_slab) < 1e-12,
          "slab endstop is a gradient (closed-loop work == 0)",
          w_slab, 0.0);
    const auto a = probe(0.83, 0.005, false, 0.0);
    const auto b = probe(0.83, 0.079, false, 0.0);
    std::printf("4. LEGACY at s=0.830: F(r=5nm)=%.3f pN vs F(r=79nm)=%.3f pN "
                "at the SAME axial position, with F_r == 0\n", a.force, b.force);
    std::printf("   closed-loop work around (s 0.78-0.88) x (r 0.005-0.079): "
                "legacy %+.4g pN*um = %+.1f kT ; slab %+.3g pN*um\n",
                w_legacy, w_legacy / kT, w_slab);

    // --- 5. THE FORCE NOW GETS RUNWAY --------------------------------------
    // Walk a tip in along the axis and find the first sample at which the kill
    // condition is met.  The force must already have been acting before it.
    for (double dp : {0.079, 0.090, 0.150}) {
        double prev_force = -1.0, prev_pen = -1.0;
        double s_kill = 0.0, force_before = 0.0, runway = 0.0;
        bool found = false;
        double first_force_s = 0.0;
        bool force_seen = false;
        for (double s = 0.90; s > 0.55; s -= 0.0005) {
            const auto r = probe(s, dp, true, kWideCutoff);
            if (!force_seen && r.force > 0.0) { force_seen = true; first_force_s = s; }
            if (!found && r.penetration >= kKillNew) {
                found = true; s_kill = s; force_before = prev_force;
                runway = first_force_s - s;
                check(prev_pen >= 0.0 && prev_pen < kKillNew,
                      "penetration arrives continuously, not in a jump",
                      prev_pen, kKillNew);
                check(force_before > 0.0,
                      "force was already acting before the kill line",
                      force_before, 1.0);
            }
            prev_force = r.force; prev_pen = r.penetration;
        }
        check(found, "kill line is reached somewhere in the sweep", found ? 1 : 0, 1);
        std::printf("5. slab, dperp=%.3f: force first acts at s=%.4f, kill at "
                    "s=%.4f -> %.1f nm of runway, F just before = %.2f pN\n",
                    dp, first_force_s, s_kill, 1000.0 * runway, force_before);
    }
    // The same sweep in legacy mode outside the gate never engages at all.
    {
        bool any = false;
        for (double s = 0.90; s > 0.55; s -= 0.0005)
            if (probe(s, 0.090, false, 0.0).force > 0.0) any = true;
        check(!any, "legacy: no force anywhere along the approach at 90 nm",
              any ? 1 : 0, 0);
        std::printf("   LEGACY, dperp=0.090: force acts at NO point on the "
                    "entire approach\n");
    }

    // --- 6. THE HAZARD RATE -------------------------------------------------
    using catch_bond_kinetics::decide_barbed_tip_capsule_penetration;
    using catch_bond_kinetics::decide_barbed_tip_capsule_penetration_rate;
    const double dt = 1e-5, rate0 = 1.0;
    auto fixed = [](double v) { return [v]() { return v; }; };

    // rate0 = 0 -> no decision at all, caller keeps the legacy veto.
    {
        const auto d = decide_barbed_tip_capsule_penetration_rate(
            0.350, 0.0, kKillNew, 0.0, dt, fixed(0.0));
        check(!d.am_overlap_occluded, "rate0=0 is inert", d.am_overlap_occluded ? 1 : 0, 0);
        check(d.probability == 0.0, "rate0=0 has zero probability", d.probability, 0.0);
    }
    // Zero penetration consumes no draw and never fires.
    {
        const auto d = decide_barbed_tip_capsule_penetration_rate(
            0.0, 0.0, kKillNew, rate0, dt, fixed(0.0));
        check(!d.am_overlap_occluded, "p=0 never fires", d.am_overlap_occluded ? 1 : 0, 0);
        check(std::isnan(d.draw), "p=0 consumes no draw", d.draw, NAN);
    }
    // Functional form, and the two tips adding as independent channels.
    for (double p : {0.001, 0.010, 0.050, 0.100, 0.350}) {
        const double want = std::min(
            1.0, rate0 * std::exp((p - kKillNew) / kKillNew) * dt);
        const auto d = decide_barbed_tip_capsule_penetration_rate(
            p, 0.0, kKillNew, rate0, dt, fixed(1.0));
        check(std::fabs(d.probability - want) < 1e-15,
              "hazard = rate*exp((p-thr)/thr)", d.probability, want);
        const auto d2 = decide_barbed_tip_capsule_penetration_rate(
            p, p, kKillNew, rate0, dt, fixed(1.0));
        check(std::fabs(d2.probability - std::min(1.0, 2.0 * want)) < 1e-15,
              "two buried tips add", d2.probability, std::min(1.0, 2.0 * want));
    }
    // k = rate0 exactly at the threshold.
    {
        const auto d = decide_barbed_tip_capsule_penetration_rate(
            kKillNew, 0.0, kKillNew, rate0, dt, fixed(1.0));
        check(std::fabs(d.probability - rate0 * dt) < 1e-18,
              "k = rate0 at p = threshold", d.probability, rate0 * dt);
    }
    // The draw decides.
    {
        const double p = 0.350;
        const double prob = rate0 * std::exp((p - kKillNew) / kKillNew) * dt;
        check(decide_barbed_tip_capsule_penetration_rate(
                  p, 0.0, kKillNew, rate0, dt, fixed(prob * 0.5))
                  .am_overlap_occluded,
              "draw below probability fires", 1, 1);
        check(!decide_barbed_tip_capsule_penetration_rate(
                   p, 0.0, kKillNew, rate0, dt, fixed(prob * 2.0))
                   .am_overlap_occluded,
              "draw above probability does not fire", 0, 0);
    }
    // The legacy veto, for contrast: probability 1 at the same depths.
    {
        const auto v = decide_barbed_tip_capsule_penetration(0.0501, 0.0, kKillNew);
        check(v.probability == 1.0, "legacy veto is a certainty", v.probability, 1.0);
    }
    std::printf("6. hazard form, additivity, threshold value and draw: checked\n");

    // --- 7. TRANSIENT vs SUSTAINED -----------------------------------------
    // The measured burial episode is a median 20 us (2 steps at dt=1e-5).
    auto survival = [&](double p, int steps) {
        const double q = std::min(
            1.0, rate0 * std::exp((p - kKillNew) / kKillNew) * dt);
        return std::pow(1.0 - q, steps);
    };
    const double transient = 1.0 - survival(0.055, 2);
    const double sustained_1ms = 1.0 - survival(0.350, 100);
    check(transient < 0.01, "a 20 us dip is survivable", transient, 0.01);
    std::printf("7. rate0=%.3g/s: P(death | 20 us dip at 55 nm) = %.3g ; "
                "P(death | 1 ms sustained at 350 nm) = %.4f\n",
                rate0, transient, sustained_1ms);
    std::printf("   legacy veto: both of those are 1.0\n");

    gsl_rng_free(g_rng);
    if (g_fail) { std::printf("FAILED: %d checks\n", g_fail); return 1; }
    std::printf("ALL CHECKS PASSED\n");
    return 0;
}
