// Unit test for geometry::segments_may_interact.
//
// The pre-filter exists because the neighbour list is built on a centre-centre
// sphere sized by the longest filament, which at Lx=6 is larger than the box:
// ~183 myosins per actin are handed to _process_actin_myosin_binding and ~2 of
// them actually interact.  The filter is only admissible if it is a strict
// lower bound on the segment-segment distance, so the property under test is
// one-sided: it may keep a pair that does not interact, but it must never
// reject one that does.
#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "geometry.h"

namespace {

const std::vector<double> kBox{6.0, 2.0, 2.0};
const std::array<bool, 3> kPeriodic{true, true, true};
const double kCapture = 0.05 + 0.015 + 0.05 + 0.015;  // 0.13 um
const double kActinLength = 1.0;
const double kMyosinLength = 1.5;

vec unit_from_angles(double cos_theta, double phi) {
    const double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
    return vec{sin_theta * std::cos(phi), sin_theta * std::sin(phi), cos_theta};
}

// Never reject a pair whose segments are genuinely within the capture distance,
// nor one that analyze_am reports as binding.  Sampled over both isotropic and
// nematic orientations, since the bound's tightness depends on sin(theta).
void randomized_no_false_rejection(bool nematic, const char* label) {
    std::mt19937_64 rng(nematic ? 20260901u : 20260902u);
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    long long tested = 0, interacting = 0, rejected = 0, kept = 0;
    long long false_rejects = 0, binding_rejects = 0;

    for (int trial = 0; trial < 400000; ++trial) {
        // Bias the centres together so a useful fraction of samples interact.
        vec actin_center{u01(rng) * kBox[0], u01(rng) * kBox[1], u01(rng) * kBox[2]};
        vec offset{(u01(rng) - 0.5) * 2.0, (u01(rng) - 0.5) * 0.6,
                   (u01(rng) - 0.5) * 0.6};
        vec myosin_center = actin_center + offset;

        vec actin_dir, myosin_dir;
        if (nematic) {
            // Tight cone about x, matching the S ~ 0.95 of a production frame.
            actin_dir = vec{1.0, (u01(rng) - 0.5) * 0.4, (u01(rng) - 0.5) * 0.4};
            myosin_dir = vec{1.0, (u01(rng) - 0.5) * 0.4, (u01(rng) - 0.5) * 0.4};
        } else {
            actin_dir = unit_from_angles(2.0 * u01(rng) - 1.0,
                                         2.0 * M_PI * u01(rng));
            myosin_dir = unit_from_angles(2.0 * u01(rng) - 1.0,
                                          2.0 * M_PI * u01(rng));
        }
        actin_dir.normalize();
        myosin_dir.normalize();

        const vec a_left = actin_center - 0.5 * kActinLength * actin_dir;
        const vec a_right = actin_center + 0.5 * kActinLength * actin_dir;
        const vec m_left = myosin_center - 0.5 * kMyosinLength * myosin_dir;
        const vec m_right = myosin_center + 0.5 * kMyosinLength * myosin_dir;

        const bool may = geometry::segments_may_interact(
            a_left, a_right, m_left, m_right, kCapture, kBox, kPeriodic);
        const double distance = geometry::segment_segment_distance(
            a_left, a_right, m_left, m_right, kBox, kPeriodic);
        const auto full = geometry::analyze_am(a_left, a_right, m_left, m_right,
                                               kCapture, kBox, kPeriodic, true);

        ++tested;
        may ? ++kept : ++rejected;
        if (distance < kCapture) {
            ++interacting;
            if (!may) ++false_rejects;
        }
        if (!may && (full.myosin_binding_ratio != 0.0 ||
                     full.partial_binding_ratio != 0.0)) {
            ++binding_rejects;
        }
    }

    std::printf("  %-10s tested %lld, interacting %lld, rejected %lld (%.1f%%), "
                "false rejects %lld, nonzero-binding rejects %lld\n",
                label, tested, interacting, rejected,
                100.0 * double(rejected) / double(tested),
                false_rejects, binding_rejects);
    assert(false_rejects == 0);
    assert(binding_rejects == 0);
    assert(interacting > 1000);   // the sampler must actually exercise the case
    assert(rejected > tested / 4); // and the filter must actually reject things
}

// A pair straddling a periodic face must be treated as close, not as a box
// apart: the filter min-images the centre separation exactly as analyze_am and
// segment_segment_distance do.
void periodic_wrap_is_respected() {
    const vec actin_dir{1.0, 0.0, 0.0};
    const vec actin_center{0.0, 0.02, 0.0};
    const vec myosin_center{0.0, kBox[1] - 0.02, 0.0};
    const vec a_left = actin_center - 0.5 * kActinLength * actin_dir;
    const vec a_right = actin_center + 0.5 * kActinLength * actin_dir;
    const vec m_left = myosin_center - 0.5 * kMyosinLength * actin_dir;
    const vec m_right = myosin_center + 0.5 * kMyosinLength * actin_dir;
    assert(geometry::segments_may_interact(a_left, a_right, m_left, m_right,
                                          kCapture, kBox, kPeriodic));
    const double distance = geometry::segment_segment_distance(
        a_left, a_right, m_left, m_right, kBox, kPeriodic);
    assert(distance < kCapture);
    std::printf("  periodic wrap: separation %.4f um across the y face, kept\n",
                distance);
}

// A degenerate (zero-length) rod carries no orientation, so the filter must
// decline to reject rather than divide by zero.
void degenerate_rod_falls_through() {
    const vec p{1.0, 1.0, 1.0};
    assert(geometry::segments_may_interact(p, p, p, p, kCapture, kBox, kPeriodic));
    const vec q{3.0, 1.0, 1.0};
    assert(geometry::segments_may_interact(p, p, q, q, kCapture, kBox, kPeriodic));
    std::printf("  degenerate rods fall through to the exact analysis\n");
}

// Far apart along x with both rods on the axis: purely axial rejection.
void axial_rejection() {
    const vec dir{1.0, 0.0, 0.0};
    const vec a_left{0.0, 1.0, 1.0}, a_right{1.0, 1.0, 1.0};
    const vec m_left{2.0, 1.0, 1.0}, m_right{3.5, 1.0, 1.0};
    assert(!geometry::segments_may_interact(a_left, a_right, m_left, m_right,
                                           kCapture, kBox, kPeriodic));
    // Slide the myosin until the ends nearly touch: must be kept.
    const vec n_left{1.05, 1.0, 1.0}, n_right{2.55, 1.0, 1.0};
    assert(geometry::segments_may_interact(a_left, a_right, n_left, n_right,
                                          kCapture, kBox, kPeriodic));
    (void)dir;
    std::printf("  axial gate rejects a 1 um gap and keeps a 50 nm gap\n");
}


// The same bound now gates compute_actin_repulsion, whose criterion is
// segment_segment_distance_w_normal >= 2*actin_radius.  That is a different
// function from segment_segment_distance, so check the bound against both:
// two equal-length 1 um rods, capture 30 nm.
void actin_actin_capsule() {
    const double kCapsule = 2.0 * 0.015;
    const double kLength = 1.0;
    std::mt19937_64 rng(20260903u);
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    long long interacting = 0, rejected = 0, false_rejects = 0, normal_mismatch = 0;

    for (int trial = 0; trial < 400000; ++trial) {
        vec first_center{u01(rng) * kBox[0], u01(rng) * kBox[1], u01(rng) * kBox[2]};
        vec offset{(u01(rng) - 0.5) * 1.6, (u01(rng) - 0.5) * 0.2,
                   (u01(rng) - 0.5) * 0.2};
        vec second_center = first_center + offset;
        // Nematic cone, and half the pairs antiparallel as in a Z-disc.
        vec first_dir{1.0, (u01(rng) - 0.5) * 0.4, (u01(rng) - 0.5) * 0.4};
        vec second_dir{u01(rng) < 0.5 ? 1.0 : -1.0,
                       (u01(rng) - 0.5) * 0.4, (u01(rng) - 0.5) * 0.4};
        first_dir.normalize();
        second_dir.normalize();

        const vec a_left = first_center - 0.5 * kLength * first_dir;
        const vec a_right = first_center + 0.5 * kLength * first_dir;
        const vec b_left = second_center - 0.5 * kLength * second_dir;
        const vec b_right = second_center + 0.5 * kLength * second_dir;

        const bool may = geometry::segments_may_interact(
            a_left, a_right, b_left, b_right, kCapsule, kBox, kPeriodic);
        const double plain = geometry::segment_segment_distance(
            a_left, a_right, b_left, b_right, kBox, kPeriodic);
        const double with_normal = geometry::segment_segment_distance_w_normal(
            a_left, a_right, b_left, b_right, kBox, kPeriodic).first;

        if (std::abs(plain - with_normal) > 1e-12) ++normal_mismatch;
        if (with_normal < kCapsule) {
            ++interacting;
            if (!may) ++false_rejects;
        }
        if (!may) ++rejected;
    }
    std::printf("  %-10s interacting %lld, rejected %lld (%.1f%%), "
                "false rejects %lld, w_normal disagreements %lld\n",
                "capsule", interacting, rejected,
                100.0 * double(rejected) / 400000.0, false_rejects,
                normal_mismatch);
    assert(false_rejects == 0);
    assert(normal_mismatch == 0);   // the two distance functions must agree
    assert(interacting > 500);
    assert(rejected > 100000);
}

}  // namespace

int main() {
    std::printf("segments_may_interact:\n");
    randomized_no_false_rejection(false, "isotropic");
    randomized_no_false_rejection(true, "nematic");
    actin_actin_capsule();
    periodic_wrap_is_respected();
    degenerate_rod_falls_through();
    axial_rejection();
    std::printf("am_prefilter_test: all checks passed\n");
    return 0;
}
