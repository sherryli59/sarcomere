// Unit test for the azimuth-selection gate: only the crown spoke that points
// most toward an actin may bind it.
//
// The gate is a BINDING PERMISSION, deliberately not a geometry change: with
// am_head_azimuth = false the head anchor sits on the myosin axis at radius 0,
// so an azimuth used as an anchor offset multiplies zero and is provably inert
// (verified bitwise).  Gating eligibility instead keeps the axis anchor and its
// single 80 nm radial spring, avoiding the surface-anchor conflict where the
// head measures a ~30 nm rest against the pair tether's 80 nm.
//
// Properties under test:
//   T1 partition      - with spokes 2pi/n apart and half_width pi/n, EXACTLY one
//                       spoke admits any given actin azimuth.
//   T2 argmin         - admits(spoke) iff that spoke minimises |delta azimuth|.
//                       This is what "points most toward the actin" means.
//   T3 no dead gap    - the admitting count is never 0 anywhere on the circle.
//   T4 wide limit     - half_width >= pi admits everything (ungated behaviour).
//   T5 rotation       - rigidly rotating actin+myosin+frame changes nothing.
//   T6 pbc            - an actin across a periodic face decides as its unwrapped
//                       image does.
//   T7 degenerate     - an actin through the myosin axis yields NaN offset and
//                       must ADMIT, never silently forbid every head.
//   T8 axial station  - a tilted actin subtends different azimuths at different
//                       crowns, so different spokes must admit.  Guards against
//                       "simplifying" the gate to use the actin centre.
#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

#include <gsl/gsl_rng.h>

#include "components.h"
#include "interaction.h"

namespace {

const std::vector<double> kBox{6.0, 2.0, 2.0};
const double kActinLength = 1.0;
const double kMyosinLength = 1.5;
const double kMyosinRadius = 0.05;
const double kShell = 0.080;          // am_preferred_centerline
int g_failures = 0;

#define CHECK(cond, ...)                                                   \
    do { if (!(cond)) { ++g_failures;                                      \
        std::printf("  FAIL %s:%d  ", __FILE__, __LINE__);                 \
        std::printf(__VA_ARGS__); std::printf("\n"); } } while (0)

struct Rig {
    gsl_rng* rng;
    Filament actin;
    Myosin myosin;
    vec e1, e2;
    Rig() : rng(gsl_rng_alloc(gsl_rng_mt19937)),
            actin(1, kActinLength, kBox, rng),
            myosin(1, kMyosinLength, kMyosinRadius, kBox, rng) {
        actin.set_periodic_axes({true, true, true});
        myosin.set_periodic_axes({true, true, true});
        set_myosin(vec{0.0, 0.0, 0.0}, vec{1.0, 0.0, 0.0});
        e1 = vec{0.0, 1.0, 0.0};
        e2 = vec{0.0, 0.0, 1.0};
    }
    ~Rig() { gsl_rng_free(rng); }
    void set_myosin(const vec& c, const vec& d) {
        myosin.center[0] = c; myosin.direction[0] = d; myosin.update_endpoints();
    }
    void set_actin(const vec& c, const vec& d) {
        actin.center[0] = c; actin.direction[0] = d; actin.update_endpoints();
    }
    // Place an actin parallel to x at azimuth phi on the shell, centred at x0.
    void place_parallel(double phi, double x0 = 0.0, double radius = kShell) {
        set_actin(vec{x0, radius * std::cos(phi), radius * std::sin(phi)},
                  vec{1.0, 0.0, 0.0});
    }
    double offset(double crown_azimuth, double anchor_axial = 0.0) const {
        return am_head_azimuth_offset(actin, myosin, 0, 0, kBox,
                                      anchor_axial, crown_azimuth, e1, e2);
    }
    bool admits(double crown_azimuth, double half_width,
                double anchor_axial = 0.0) const {
        return am_head_azimuth_admits(actin, myosin, 0, 0, kBox, anchor_axial,
                                      crown_azimuth, half_width, e1, e2);
    }
};

std::vector<double> spokes_of(int n) {
    std::vector<double> s;
    for (int k = 0; k < n; ++k) s.push_back(2.0 * M_PI * k / n);
    return s;
}

// T1 + T2 + T3 -------------------------------------------------------------
void partition_and_argmin(int n_spokes) {
    Rig rig;
    const auto spokes = spokes_of(n_spokes);
    const double half = M_PI / n_spokes;          // the partitioning half-width
    int min_admits = 99, max_admits = 0, argmin_mismatches = 0, steps = 0;
    for (int deg = 0; deg < 360; ++deg) {
        // offset the sweep by half a step to avoid sitting exactly on a
        // boundary, which T3b covers separately
        const double phi = (deg + 0.5) * M_PI / 180.0;
        rig.place_parallel(phi);
        int n_admit = 0; int admitted = -1;
        double best = 1e9; int best_k = -1;
        for (int k = 0; k < n_spokes; ++k) {
            if (rig.admits(spokes[k], half)) { ++n_admit; admitted = k; }
            const double d = std::abs(rig.offset(spokes[k]));
            if (d < best) { best = d; best_k = k; }
        }
        min_admits = std::min(min_admits, n_admit);
        max_admits = std::max(max_admits, n_admit);
        if (n_admit == 1 && admitted != best_k) ++argmin_mismatches;
        ++steps;
    }
    std::printf("  n=%d spokes, half_width=%.1f deg: admits per azimuth in [%d,%d]"
                " over %d samples; argmin mismatches %d\n",
                n_spokes, half * 180.0 / M_PI, min_admits, max_admits,
                steps, argmin_mismatches);
    CHECK(min_admits == 1 && max_admits == 1,
          "T1 partition: expected exactly 1 admitting spoke everywhere, got [%d,%d]",
          min_admits, max_admits);
    CHECK(argmin_mismatches == 0,
          "T2 argmin: %d azimuths admitted a spoke that was not the closest",
          argmin_mismatches);
    CHECK(min_admits >= 1, "T3 dead gap: some azimuth admitted no spoke at all");
}

// T3b: exactly on a boundary, never zero ------------------------------------
void boundary_has_no_gap(int n_spokes) {
    Rig rig;
    const auto spokes = spokes_of(n_spokes);
    const double half = M_PI / n_spokes;
    int zero_count = 0, counts_seen_min = 99, counts_seen_max = 0;
    for (int k = 0; k < n_spokes; ++k) {
        const double phi = spokes[k] + half;        // exact boundary
        rig.place_parallel(phi);
        int n_admit = 0;
        for (int q = 0; q < n_spokes; ++q)
            if (rig.admits(spokes[q], half)) ++n_admit;
        if (n_admit == 0) ++zero_count;
        counts_seen_min = std::min(counts_seen_min, n_admit);
        counts_seen_max = std::max(counts_seen_max, n_admit);
    }
    std::printf("  n=%d boundaries: admit count in [%d,%d], zero-admit cases %d\n",
                n_spokes, counts_seen_min, counts_seen_max, zero_count);
    CHECK(zero_count == 0,
          "T3 boundary: %d boundary azimuths admitted NO spoke", zero_count);
    CHECK(counts_seen_max <= 2,
          "T3 boundary: %d spokes admitted at once, expected at most 2",
          counts_seen_max);
}

// T4: wide half-width admits everything -------------------------------------
void wide_limit_admits_all() {
    Rig rig;
    const auto spokes = spokes_of(3);
    int rejects = 0;
    for (int deg = 0; deg < 360; deg += 3) {
        rig.place_parallel(deg * M_PI / 180.0);
        for (double s : spokes) {
            if (!rig.admits(s, M_PI)) ++rejects;
            if (!rig.admits(s, 0.0)) ++rejects;      // <=0 also admits all
        }
    }
    std::printf("  wide/zero half-width: %d rejections (expected 0)\n", rejects);
    CHECK(rejects == 0, "T4 wide limit: gate rejected %d pairs at half_width>=pi"
                        " or <=0", rejects);
}

// T5: rigid rotation invariance ---------------------------------------------
void rotation_invariant() {
    const double phi = 0.7, crown = 1.1, half = M_PI / 3.0;
    Rig a; a.place_parallel(phi);
    const bool ref_admit = a.admits(crown, half);
    const double ref_off = a.offset(crown);
    // Rotate everything by 50 deg about z: myosin axis, material frame, actin.
    const double t = 50.0 * M_PI / 180.0, c = std::cos(t), s = std::sin(t);
    auto rot = [&](const vec& v) {
        return vec{c * v.x - s * v.y, s * v.x + c * v.y, v.z};
    };
    Rig b;
    b.set_myosin(vec{0, 0, 0}, rot(vec{1, 0, 0}));
    b.e1 = rot(vec{0, 1, 0});
    b.e2 = rot(vec{0, 0, 1});
    const vec ac = vec{0.0, kShell * std::cos(phi), kShell * std::sin(phi)};
    b.set_actin(rot(ac), rot(vec{1, 0, 0}));
    const double got_off = b.offset(crown);
    const bool got_admit = b.admits(crown, half);
    std::printf("  rotation: offset %.12f -> %.12f (delta %.3e), admit %d -> %d\n",
                ref_off, got_off, std::abs(got_off - ref_off),
                int(ref_admit), int(got_admit));
    CHECK(std::abs(got_off - ref_off) < 1e-9,
          "T5 rotation: offset changed by %.3e under a rigid rotation",
          std::abs(got_off - ref_off));
    CHECK(got_admit == ref_admit, "T5 rotation: admit decision flipped");
}

// T6: periodic image equivalence --------------------------------------------
void pbc_equivalent() {
    const double phi = 2.1, crown = 2.0, half = M_PI / 3.0;
    Rig a;
    a.set_myosin(vec{0.1, 0, 0}, vec{1, 0, 0});
    a.place_parallel(phi, 0.1);
    const double ref = a.offset(crown);
    Rig b;   // same relative geometry, shifted so the actin sits across x = Lx
    b.set_myosin(vec{kBox[0] - 0.05, 0, 0}, vec{1, 0, 0});
    b.set_actin(vec{kBox[0] + 0.05 - kBox[0],           // wrapped to -0.0.. side
                    kShell * std::cos(phi), kShell * std::sin(phi)},
                vec{1, 0, 0});
    const double got = b.offset(crown);
    std::printf("  pbc: offset %.12f vs wrapped %.12f (delta %.3e)\n",
                ref, got, std::abs(got - ref));
    CHECK(std::abs(got - ref) < 1e-9,
          "T6 pbc: wrapped actin gave a different azimuth offset (delta %.3e)",
          std::abs(got - ref));
}

// T7: degenerate geometry admits --------------------------------------------
void degenerate_admits() {
    Rig rig;
    rig.set_actin(vec{0.0, 0.0, 0.0}, vec{1.0, 0.0, 0.0});   // on the axis
    const double off = rig.offset(0.0);
    const bool admit = rig.admits(0.0, M_PI / 3.0);
    std::printf("  degenerate (actin on the axis): offset %s, admit %d\n",
                std::isfinite(off) ? "finite" : "NaN", int(admit));
    CHECK(!std::isfinite(off), "T7 degenerate: expected a NaN offset on the axis");
    CHECK(admit, "T7 degenerate: the gate FORBADE binding where the azimuth is "
                 "meaningless; that silently blocks every head");
    // a zero-length actin direction is equally degenerate
    Rig z;
    z.set_actin(vec{0.0, kShell, 0.0}, vec{0.0, 0.0, 0.0});
    CHECK(z.admits(0.0, M_PI / 3.0),
          "T7 degenerate: zero-length actin direction must admit");
}

// T8: the gate must use the binding station, not the actin centre ------------
void axial_station_matters() {
    Rig rig;
    // An actin whose centre sits at azimuth 0 but which tilts in z, so at
    // negative x it is below the axis and at positive x above it.  Its azimuth
    // therefore differs strongly between two crowns on either side.
    const double tilt = 0.35;         // radians away from x, into z
    rig.set_actin(vec{0.0, kShell, 0.0},
                  vec{std::cos(tilt), 0.0, std::sin(tilt)});
    const double off_lo = rig.offset(0.0, -0.45);
    const double off_hi = rig.offset(0.0, +0.45);
    std::printf("  axial station: offset at anchor -0.45 = %+.4f rad, "
                "at +0.45 = %+.4f rad, spread %.4f\n",
                off_lo, off_hi, std::abs(off_hi - off_lo));
    CHECK(std::isfinite(off_lo) && std::isfinite(off_hi),
          "T8 axial: offsets should be finite");
    CHECK(std::abs(off_hi - off_lo) > 0.2,
          "T8 axial: the offset barely changed between crowns (%.4f rad); the "
          "gate is probably using the actin centre instead of the binding "
          "station", std::abs(off_hi - off_lo));
}

}  // namespace

int main() {
    std::printf("am_azimuth_gate_test\n");
    std::printf("T1/T2/T3 partition and argmin:\n");
    partition_and_argmin(3);
    partition_and_argmin(4);
    std::printf("T3b boundaries:\n");
    boundary_has_no_gap(3);
    boundary_has_no_gap(4);
    std::printf("T4 wide limit:\n");       wide_limit_admits_all();
    std::printf("T5 rotation invariance:\n"); rotation_invariant();
    std::printf("T6 periodic image:\n");   pbc_equivalent();
    std::printf("T7 degenerate geometry:\n"); degenerate_admits();
    std::printf("T8 axial station:\n");    axial_station_matters();
    if (g_failures) {
        std::printf("am_azimuth_gate_test: %d FAILURES\n", g_failures);
        return 1;
    }
    std::printf("am_azimuth_gate_test: all checks passed\n");
    return 0;
}
