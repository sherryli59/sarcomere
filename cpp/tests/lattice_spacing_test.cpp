#include "lattice_spacing.h"

#include <array>
#include <cmath>
#include <stdexcept>

namespace {

bool near(double lhs, double rhs, double tolerance = 1e-12) {
    return std::abs(lhs - rhs) <= tolerance;
}

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

}  // namespace

int main() {
    const double am = lattice_spacing::am_preferred_centerline(
        0.025, 0.015, 0.015);
    require(near(am, 0.055), "incorrect AM centerline");

    // Pin the constant explicitly, so an accidental edit is still caught.
    // 0.9875 = 79/80: the ACTN2 rest length must agree with
    // aa_barbed_lateral_range, a one-sided wall at 1.00 * am_preferred_
    // centerline and 6x stiffer than ACTN2.  The old 0.92 lost that contest
    // and left every link ~5.4 nm strained.
    require(near(lattice_spacing::physiological_aa_am_ratio, 0.9875),
            "physiological AA/AM ratio changed unexpectedly");

    // The round-trip assertions below derive from the constant rather than
    // duplicating its value -- duplicating it is what silently broke this
    // test when the ratio moved.
    const double aa_default = lattice_spacing::aa_preferred_centerline(
        am, lattice_spacing::physiological_aa_am_ratio);
    require(near(aa_default, am * lattice_spacing::physiological_aa_am_ratio),
            "incorrect default AA centerline");
    require(near(lattice_spacing::effective_ratio(aa_default, am),
                 lattice_spacing::physiological_aa_am_ratio),
            "incorrect default spacing ratio");

    for (double rho : std::array<double, 3>{0.8, 0.9, 1.0}) {
        const double aa = lattice_spacing::aa_preferred_centerline(am, rho);
        require(near(lattice_spacing::effective_ratio(aa, am), rho),
                "rho sweep does not round-trip");
    }

    // Cardiac cryo-ET/Z-disc values: 26 nm AM and 24 nm AA centerlines.
    require(near(lattice_spacing::effective_ratio(24.0, 26.0),
                 24.0 / 26.0),
            "cardiac spacing ratio is incorrect");
    require(near(lattice_spacing::projected_surface_span(24.0, 4.2), 15.6),
            "projected ACTN2 surface span is incorrect");
    return 0;
}
