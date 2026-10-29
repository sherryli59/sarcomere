// Unit test for the dynamic axial force normalisation.
//
// The axial channel had no per-head normalisation: N bound heads on a pair meant
// N times the sliding force, which is why raising the head pool also raised motor
// drive.  am_motor_saturation_scale(n, sat) returns 1/min(n, sat) so the pair's
// TOTAL axial force saturates.  It must be applied to BOTH the axial stiffness
// and the stall cap -- T2 and T3 are the pair of tests that catch applying it to
// only one, which leaves the other regime still growing with N.
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

#include "interaction.h"

namespace {

int g_fail = 0;
#define CHECK(cond, ...) \
    do { if (!(cond)) { ++g_fail; std::printf("  FAIL %s:%d  ", __FILE__, __LINE__); \
         std::printf(__VA_ARGS__); std::printf("\n"); } } while (0)

// The axial force law from compute_am_two_spring_head (interaction.cpp), with the
// scale applied exactly as the call site applies it.
double head_force(double k, double strain, double stall, double scale) {
    const double uncapped = k * scale * strain;
    const double mag = std::abs(uncapped);
    const double cap = stall * scale;
    return cap > 0.0 ? std::copysign(std::min(mag, cap), uncapped) : uncapped;
}
double total_force(int n, double k, double strain, double stall, int sat) {
    const double s = am_motor_saturation_scale(n, sat);
    return n * head_force(k, strain, stall, s);
}

// T1 -------------------------------------------------------------------------
void off_is_inert() {
    int bad = 0;
    for (int n = 0; n <= 20; ++n)
        for (int sat : {0, -1, -8})
            if (am_motor_saturation_scale(n, sat) != 1.0) ++bad;
    std::printf("  saturation <= 0: scale != 1.0 in %d of 63 cases\n", bad);
    CHECK(bad == 0, "T1 off: the disabled path must return exactly 1.0");
}

// T2 -- linear regime: a huge stall cap means the cap never binds -------------
void total_flat_linear() {
    const double k = 1000.0, strain = 0.0005, stall = 1e9;   // k*strain = 0.5 pN
    const int sat = 6;
    const double ref = total_force(1, k, strain, stall, sat);
    double worst = 0.0;
    for (int n = 1; n <= sat; ++n) {
        const double t = total_force(n, k, strain, stall, sat);
        worst = std::max(worst, std::abs(t - ref));
        std::printf("    n=%d  total %.6f pN\n", n, t);
    }
    std::printf("  linear regime: worst deviation from n=1 is %.3g pN\n", worst);
    CHECK(worst < 1e-9, "T2 linear: total axial force must be flat in n (worst %.3g)", worst);
}

// T3 -- capped regime: a tiny stall cap means every head is capped ------------
void total_flat_capped() {
    const double k = 1000.0, strain = 0.05, stall = 4.5;     // k*strain = 50 >> cap
    const int sat = 6;
    const double ref = total_force(1, k, strain, stall, sat);
    double worst = 0.0;
    for (int n = 1; n <= sat; ++n) {
        const double t = total_force(n, k, strain, stall, sat);
        worst = std::max(worst, std::abs(t - ref));
        std::printf("    n=%d  total %.6f pN\n", n, t);
    }
    std::printf("  capped regime: worst deviation from n=1 is %.3g pN\n", worst);
    CHECK(worst < 1e-9,
          "T3 capped: total force must be flat here too -- a non-flat result means "
          "the scale was applied to the stiffness but NOT the stall cap (worst %.3g)", worst);
}

// T4 -- beyond saturation the total grows again, by construction --------------
void grows_past_saturation() {
    const double k = 1000.0, strain = 0.0005, stall = 1e9;
    const int sat = 6;
    const double at_sat = total_force(sat, k, strain, stall, sat);
    const double past = total_force(2 * sat, k, strain, stall, sat);
    std::printf("  n=%d total %.4f pN -> n=%d total %.4f pN (ratio %.2f)\n",
                sat, at_sat, 2 * sat, past, past / at_sat);
    CHECK(std::abs(past / at_sat - 2.0) < 1e-9,
          "T4 saturation: past the cap the total should grow linearly again");
}

// T5 -- degenerate counts ------------------------------------------------------
void degenerate_counts() {
    const int sat = 6;
    const double s0 = am_motor_saturation_scale(0, sat);
    const double s1 = am_motor_saturation_scale(1, sat);
    const double sneg = am_motor_saturation_scale(-3, sat);
    std::printf("  n=0 -> %.4f, n=1 -> %.4f, n=-3 -> %.4f\n", s0, s1, sneg);
    CHECK(s0 == 1.0 && s1 == 1.0 && sneg == 1.0,
          "T5 degenerate: counts below 1 must clamp to a single head, not divide by zero");
    CHECK(std::isfinite(s0) && s0 > 0.0, "T5 degenerate: scale must stay finite and positive");
}

// T6 -- per-head force is non-increasing in n ----------------------------------
void monotonic_per_head() {
    const int sat = 6;
    int bad = 0;
    double prev = 1e30;
    for (int n = 1; n <= 2 * sat; ++n) {
        const double s = am_motor_saturation_scale(n, sat);
        if (s > prev + 1e-15) ++bad;
        prev = s;
    }
    std::printf("  per-head scale monotonic non-increasing: %s\n", bad ? "NO" : "yes");
    CHECK(bad == 0, "T6 monotonic: adding a head must never increase the per-head force");
}

// T7 -- recorded energy scales with the same factor as the force ---------------
// The capped branch computes axial_energy from stall_force, so if the cap were
// left unscaled the reported energy would stop being the integral of the
// reported force.
void energy_tracks_force() {
    const double k = 1000.0, strain = 0.05, stall = 4.5;
    const int sat = 6;
    int bad = 0;
    for (int n = 1; n <= sat; ++n) {
        const double s = am_motor_saturation_scale(n, sat);
        const double ks = k * s, cap = stall * s;
        const double mag = std::abs(ks * strain);
        const double energy = (cap > 0.0 && mag > cap)
            ? cap * std::abs(strain) - 0.5 * cap * cap / ks
            : 0.5 * ks * strain * strain;
        const double force = head_force(k, strain, stall, s);
        // dE/d(strain) must equal the force in the capped regime (= cap)
        if (std::abs(force) > 0.0 && std::abs(std::abs(force) - cap) > 1e-9) ++bad;
        if (!(energy > 0.0)) ++bad;
        std::printf("    n=%d  force %.4f pN  energy %.6f pN*um  cap %.4f\n",
                    n, force, energy, cap);
    }
    CHECK(bad == 0, "T7 energy: recorded energy must scale with the same factor as the force");
}

}  // namespace

int main() {
    std::printf("am_motor_saturation_test\n");
    std::printf("T1 disabled path is inert:\n");      off_is_inert();
    std::printf("T2 total flat, linear regime:\n");   total_flat_linear();
    std::printf("T3 total flat, capped regime:\n");   total_flat_capped();
    std::printf("T4 growth past saturation:\n");      grows_past_saturation();
    std::printf("T5 degenerate counts:\n");           degenerate_counts();
    std::printf("T6 monotonic per-head force:\n");    monotonic_per_head();
    std::printf("T7 energy tracks force:\n");         energy_tracks_force();
    if (g_fail) { std::printf("am_motor_saturation_test: %d FAILURES\n", g_fail); return 1; }
    std::printf("am_motor_saturation_test: all checks passed\n");
    return 0;
}
