// Sweep the actin-myosin separation and print the radial force AS IMPLEMENTED.
//
// Every number comes from the shipped interaction functions -- no reimplementation
// of any force law here.  Sign convention: + = OUTWARD (away from the myosin axis).
//
// Geometry: myosin along +x through the origin, actin parallel to it and offset
// along -y by d.  For parallel filaments rho projects out the axial offset, so a
// head's radial force is the same at every crown station and scaling by the bound
// head count is exact.
#include "components.h"
#include "interaction.h"
#include "geometry.h"

#include <cmath>
#include <cstdio>
#include <gsl/gsl_rng.h>
#include <vector>

namespace {
constexpr double kMyosinRadius = 0.050;
constexpr double kActinRadius  = 0.015;
constexpr double kSurfaceGap   = 0.015;
constexpr double kCenterline   = 0.080;
constexpr double kLegacyRest   = kActinRadius + kSurfaceGap;
constexpr double kExcluded     = kMyosinRadius + kActinRadius;   // 65 nm
constexpr double kStericK      = 1000.0;
constexpr double kLegacyHeadK  = 100.0;
constexpr double kTetherK      = 3000.0;
constexpr double kFoldedHeadK  = 600.0;                          // 100 + 3000/6
constexpr double kSlack        = 0.002;
// Derived from the integrator clamp, matching main.cpp:
//   F_clamp = max_actin_displacement/(mu_perp*dt) = 0.005/(43.34*5e-6) = 23.07 pN
//   per head = F_clamp / am_tether_saturation_heads = 23.07/6 = 3.85 pN
constexpr double kCap          = 3.846;
constexpr int    kHeads        = 7;
constexpr int    kSat          = 6;
constexpr double kLegacyAzim   = 30.0 * M_PI / 180.0;            // clamp limit
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
    actin.direction[0] = {1.0, 0.0, 0.0};
    myosin.update_endpoints();

    const vec rhat{0.0, -1.0, 0.0};   // outward: actin sits at -y
    const double occ = std::min(1.0, static_cast<double>(kHeads) /
                                     static_cast<double>(kSat));

    std::printf("d_nm,head_off,head_on,tether_on,hard_core,net_off,net_on\n");
    for (int i = 0; i <= 700; ++i) {
        const double d = 0.055 + i * (0.145 - 0.055) / 700.0;
        actin.center[0] = {0.0, -d, 0.0};
        actin.update_endpoints();

        // new mode: anchor on the axis, centreline rest, folded k, saturating
        const auto off = compute_am_two_spring_head(
            actin, myosin, 0, 0, box, 0.0, 0.0, 0.0, /*radial_offset=*/0.0,
            nullptr, nullptr, kCenterline, kFoldedHeadK, kSlack,
            0.0, 0.0, 4.5, false, false, /*radial_force_max=*/kCap);
        // legacy mode: surface anchor at the clamp limit, 30 nm rest, uncapped
        const auto on = compute_am_two_spring_head(
            actin, myosin, 0, 0, box, 0.0, 0.0, kLegacyAzim, kMyosinRadius,
            nullptr, nullptr, kLegacyRest, kLegacyHeadK, kSlack,
            0.0, 0.0, 4.5, false);
        // legacy aggregate pair tether
        const auto teth = compute_am_head_tether(
            actin, myosin, 0, 0, box, 0.0, kTetherK * occ, kCenterline,
            0.0, 0.0, nullptr, nullptr, kActinRadius, kSlack, 0.0);
        // excluded volume, live in both modes
        const auto hard = compute_actin_myosin_radial_repulsion_result(
            actin, myosin, 0, 0, box, kExcluded, kStericK);

        const double f_off  = off.radial_force_on_actin.dot(rhat) * kHeads;
        const double f_on   = on.radial_force_on_actin.dot(rhat) * kHeads;
        const double f_teth = teth.force_on_actin.dot(rhat);
        const double f_hard = hard.applied ? hard.force_on_first.dot(rhat) : 0.0;
        std::printf("%.4f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                    d * 1e3, f_off, f_on, f_teth, f_hard,
                    f_off + f_hard, f_on + f_teth + f_hard);
    }
    gsl_rng_free(rng);
    return 0;
}
