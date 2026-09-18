// Azimuth-OFF radial mode: does it deliver the three properties it claims?
//
// In the new default mode the head anchor sits on the myosin AXIS (radial
// offset 0) and the radial spring is measured against am_preferred_centerline
// instead of actin_radius + am_surface_gap.  The claims are:
//
//   1. AZIMUTH INDEPENDENCE.  rho equals the centreline separation at every
//      azimuth, so a head bound off-axis is not strained by geometry alone.
//      In the legacy mode rho = sqrt(d^2 + R^2 - 2dR cos phi), which is 30 nm
//      at phi = 0 but 44.4 nm at phi = 30 deg against a 30 nm rest -- a 14 nm
//      stretch from nothing but the anchor's position.
//   2. BOUNDED FORCE.  With am_head_radial_force_max > 0 the radial force
//      saturates instead of growing linearly without limit.  The legacy law
//      reaches 243 pN at 140 nm separation, an order of magnitude past the
//      integrator's own displacement clamp.
//   3. THE COUPLE SURVIVES.  Heads at different axial stations on a TILTED
//      actin must see different local separations, with extensions of opposite
//      sign about the pivot -- that difference is the restoring torque.  A
//      version of this fix that used the pair's GLOBAL closest approach would
//      give every head one number and kill the couple; this asserts it does not.
#include "components.h"
#include "interaction.h"
#include "geometry.h"

#include <cmath>
#include <cstdio>
#include <gsl/gsl_rng.h>
#include <vector>

namespace {

constexpr double kMyosinRadius        = 0.050;
constexpr double kActinRadius         = 0.015;
constexpr double kSurfaceGap          = 0.015;
constexpr double kPreferredCenterline = 0.080;
constexpr double kLegacyRest          = kActinRadius + kSurfaceGap;  // 0.030
constexpr double kHeadK               = 600.0;   // folded per-head, pN/um
constexpr double kSlack               = 0.002;   // am_tether_slack
constexpr double kForceMax            = 12.0;    // am_head_radial_force_max, pN

int failures = 0;

void check(bool ok, const char* what) {
    std::printf(ok ? "  PASS  %s\n" : "  FAIL  %s\n", what);
    if (!ok) ++failures;
}

// One head in the NEW mode: anchor on the axis, rest = centreline.
AMTwoSpringHeadResult head_off(const Filament& actin, const Myosin& myosin,
                               const std::vector<double>& box,
                               double actin_coordinate, double anchor_axial,
                               double phi, double force_max) {
    return compute_am_two_spring_head(
        actin, myosin, 0, 0, box, actin_coordinate, anchor_axial,
        /*anchor_azimuth=*/phi, /*radial_offset=*/0.0,
        nullptr, nullptr,
        /*radial_preferred=*/kPreferredCenterline,
        /*radial_stiffness=*/kHeadK, /*radial_slack=*/kSlack,
        /*axial_stiffness=*/0.0,          // isolate the radial channel
        /*stroke_progress=*/0.0, /*stall_force=*/4.5, /*strong=*/false,
        /*bidirectional_strong=*/false, /*radial_force_max=*/force_max);
}

}  // namespace

int main() {
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, 17);
    const std::vector<double> box{10.0, 10.0, 10.0};
    Filament actin(1, 1.0, box, rng);
    Myosin myosin(1, 1.5, kMyosinRadius, box, rng);
    actin.set_periodic_axes({false, false, false});
    myosin.set_periodic_axes({false, false, false});
    myosin.center[0] = {0.0, 0.0, 0.0};
    myosin.direction[0] = {1.0, 0.0, 0.0};
    actin.center[0] = {0.0, -kPreferredCenterline, 0.0};
    actin.direction[0] = {1.0, 0.0, 0.0};
    actin.update_endpoints();
    myosin.update_endpoints();

    // ---- 1. azimuth independence ----------------------------------------
    std::printf("1. AZIMUTH INDEPENDENCE at the preferred spacing (%.0f nm)\n",
                kPreferredCenterline * 1e3);
    std::printf("   %8s %12s %12s %14s %14s\n", "azimuth",
                "rho OFF nm", "|F| OFF pN", "rho LEGACY nm", "|F| LEGACY pN");
    bool all_rest = true, legacy_varies = false;
    for (double phi_deg : {0.0, 15.0, 30.0, 45.0, 90.0, 180.0}) {
        const double phi = phi_deg * M_PI / 180.0;
        const auto off = head_off(actin, myosin, box, 0.0, 0.0, phi, kForceMax);
        const auto leg = compute_am_two_spring_head(
            actin, myosin, 0, 0, box, 0.0, 0.0, phi, kMyosinRadius,
            nullptr, nullptr, kLegacyRest, 100.0, kSlack, 0.0, 0.0, 4.5, false);
        std::printf("   %7.0f  %12.3f %12.4f %14.3f %14.4f\n", phi_deg,
                    off.transverse_distance * 1e3, off.radial_force_on_actin.norm(),
                    leg.transverse_distance * 1e3, leg.radial_force_on_actin.norm());
        if (std::abs(off.transverse_distance - kPreferredCenterline) > 1e-12)
            all_rest = false;
        if (off.radial_force_on_actin.norm() > 1e-9) all_rest = false;
        if (std::abs(leg.transverse_distance - kLegacyRest) > 1e-3)
            legacy_varies = true;
    }
    check(all_rest, "OFF mode: rho equals the centreline and the force vanishes "
                    "at EVERY azimuth");
    check(legacy_varies, "legacy mode does vary with azimuth (guards the "
                         "comparison, so a null result cannot be a silent no-op)");

    // ---- 2. the force saturates ------------------------------------------
    std::printf("\n2. BOUNDED FORCE, cap = %.0f pN\n", kForceMax);
    std::printf("   %12s %14s %16s\n", "separation", "|F| capped pN", "|F| uncapped pN");
    double max_capped = 0.0, uncapped_far = 0.0;
    for (double d_nm : {82.0, 90.0, 100.0, 140.0, 200.0, 400.0}) {
        actin.center[0] = {0.0, -d_nm * 1e-3, 0.0};
        actin.update_endpoints();
        const double fc = head_off(actin, myosin, box, 0.0, 0.0, 0.0,
                                   kForceMax).radial_force_on_actin.norm();
        const double fu = head_off(actin, myosin, box, 0.0, 0.0, 0.0,
                                   0.0).radial_force_on_actin.norm();
        std::printf("   %10.0f nm %14.3f %16.3f\n", d_nm, fc, fu);
        max_capped = std::max(max_capped, fc);
        if (d_nm == 400.0) uncapped_far = fu;
    }
    check(max_capped <= kForceMax + 1e-6,
          "capped radial force never exceeds am_head_radial_force_max");
    // One head at 600 pN/um reaches 190.8 pN by 400 nm; the 243 pN figure
    // quoted elsewhere is for the whole 7-head ensemble, not a single head.
    check(uncapped_far > 10.0 * kForceMax,
          "the uncapped law really does diverge (so the cap is doing work)");

    // ---- 3. the couple survives on a tilted actin ------------------------
    std::printf("\n3. THE COUPLE: a tilted actin must strain heads at different\n"
                "   axial stations by DIFFERENT amounts, with opposite signs\n"
                "   about the pivot.\n");
    // The tilt MUST lie in the plane containing the actin-myosin offset.  The
    // offset here is along -y, so tilting in x-y brings one end closer and the
    // other farther and produces a couple.  Tilting in x-z instead moves both
    // ends perpendicular to the offset, straining them EQUALLY (+82.19 nm at
    // both +/-0.4125 um) -- a symmetric stretch with zero torque.  That
    // degenerate orientation is easy to pick by accident.
    const double tilt = 20.0 * M_PI / 180.0;
    actin.center[0] = {0.0, -kPreferredCenterline, 0.0};
    actin.direction[0] = {std::cos(tilt), std::sin(tilt), 0.0};
    actin.update_endpoints();
    std::printf("   %14s %14s %14s\n", "station (um)", "rho (nm)", "extension (nm)");
    double ext_neg = 0.0, ext_pos = 0.0;
    for (double a : {-0.4125, -0.1875, 0.1875, 0.4125}) {
        const auto h = head_off(actin, myosin, box, a, a, 0.0, kForceMax);
        const double ext = (h.transverse_distance - kPreferredCenterline) * 1e3;
        std::printf("   %14.4f %14.3f %14.4f\n", a,
                    h.transverse_distance * 1e3, ext);
        if (a < 0.0) ext_neg = ext;
        if (a > 0.0) ext_pos = ext;
    }
    check(ext_pos * ext_neg < 0.0,
          "extensions have OPPOSITE SIGNS about the pivot -- a restoring "
          "couple exists, so folding the radial channel per-head did not "
          "destroy the alignment mechanism");
    check(std::abs(ext_pos - ext_neg) > 1e-6,
          "extensions differ between the two sides");

    gsl_rng_free(rng);
    std::printf("\n%s\n", failures ? "FAILURES" : "all checks passed");
    return failures ? 1 : 0;
}
