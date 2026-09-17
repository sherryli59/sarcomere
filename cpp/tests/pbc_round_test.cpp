// utils::round_half_away replaces std::round in the periodic minimum image,
// which perf showed costing 6.1% of all cycles as a PLT call.  A replacement is
// only admissible if it is bit-identical, including on the sign of zero and on
// halfway cases, so this compares the two exhaustively over the regimes the
// wrap actually sees plus every special value.
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <random>
#include <vector>

#include "utils.h"

namespace {

bool same_bits(double a, double b) {
    if (std::isnan(a) && std::isnan(b)) return true;
    std::uint64_t x, y;
    std::memcpy(&x, &a, sizeof x);
    std::memcpy(&y, &b, sizeof y);
    return x == y;   // catches -0.0 vs +0.0, which == would not
}

void check(double v, const char* what) {
    const double expected = std::round(v);
    const double actual = utils::round_half_away(v);
    if (!same_bits(expected, actual)) {
        std::printf("MISMATCH (%s): x=%.17g std::round=%.17g ours=%.17g\n",
                    what, v, expected, actual);
        assert(false);
    }
}

void edge_cases() {
    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const std::vector<double> values{
        0.0, -0.0, 0.5, -0.5, 1.5, -1.5, 2.5, -2.5, 0.49999999999999994,
        -0.49999999999999994, 0.5000000000000001, -0.5000000000000001,
        1.0, -1.0, 3.0, -3.0,
        std::nextafter(0.5, 0.0), std::nextafter(0.5, 1.0),
        4503599627370495.5, -4503599627370495.5,          // just under 2^52
        4503599627370496.0, -4503599627370496.0,          // 2^52 exactly
        9007199254740993.0, 1e300, -1e300,
        std::numeric_limits<double>::denorm_min(),
        -std::numeric_limits<double>::denorm_min(),
        std::numeric_limits<double>::min(),
        std::numeric_limits<double>::max(),
        inf, -inf, nan};
    for (double v : values) check(v, "edge");
    std::printf("  %zu edge cases match (including +-0, +-inf, NaN, 2^52)\n",
                values.size());
}

void randomized() {
    std::mt19937_64 rng(20260901u);
    // The wrap computes round(disp / L); with |disp| <~ a few boxes the
    // argument is small, but sample far wider than that.
    std::uniform_real_distribution<double> small(-8.0, 8.0);
    std::uniform_real_distribution<double> wide(-1e17, 1e17);
    for (int k = 0; k < 4000000; ++k) check(small(rng), "small");
    for (int k = 0; k < 1000000; ++k) check(wide(rng), "wide");
    // Dense sweep across every halfway point in the range the wrap sees.
    for (int n = -32; n <= 32; ++n) {
        const double half = n + 0.5;
        check(half, "halfway");
        check(std::nextafter(half, -1e18), "just below halfway");
        check(std::nextafter(half, 1e18), "just above halfway");
    }
    std::printf("  5,000,000 random values and 195 halfway neighbours match\n");
}

// The property that actually matters downstream: the wrap itself is unchanged.
void wrap_axis_matches_reference() {
    std::mt19937_64 rng(20260902u);
    std::uniform_real_distribution<double> pos(-20.0, 20.0);
    const double boxes[] = {6.0, 2.0, 0.5, 13.7};
    for (int k = 0; k < 1000000; ++k) {
        const double x = pos(rng);
        for (double L : boxes) {
            const double reference = x - L * std::round(x / L);
            const double actual = utils::wrap_axis(x, L, true);
            assert(same_bits(reference, actual));
        }
    }
    assert(utils::wrap_axis(3.0, 0.0, true) == 3.0);    // zero box is a no-op
    assert(utils::wrap_axis(3.0, 6.0, false) == 3.0);   // aperiodic is a no-op
    std::printf("  4,000,000 wrap_axis evaluations bit-identical\n");
}

}  // namespace

int main() {
    std::printf("round_half_away vs std::round:\n");
    edge_cases();
    randomized();
    wrap_axis_matches_reference();
    std::printf("pbc_round_test: all checks passed\n");
    return 0;
}
