#include "sarcomere.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <string>
#include <vector>

namespace {

int occupied(const Sarcomere& model, int half_index, int actin_index = -2) {
    int count = 0;
    for (const auto& head : model.am_heads_by_myosin[0][half_index]) {
        const bool bound = head.bound_actin >= 0 &&
            head.state != static_cast<Sarcomere::BondFlag>(
                am_motor_kinetics::State::detached) &&
            head.state != static_cast<Sarcomere::BondFlag>(
                am_motor_kinetics::State::off_srx);
        if (!bound) continue;
        if (actin_index == -2 || head.bound_actin == actin_index) ++count;
    }
    return count;
}

void detach_all(Sarcomere& model) {
    for (auto& half_heads : model.am_heads_by_myosin[0]) {
        for (auto& head : half_heads) {
            head.bound_actin = -1;
            head.state = static_cast<Sarcomere::BondFlag>(
                am_motor_kinetics::State::detached);
            head.stroke_progress = 0.0;
        }
    }
}

}  // namespace

int main() {
    int n_actins = 4;
    int n_myosins = 1;
    std::vector<double> box{4.0, 4.0, 4.0};
    double actin_length = 1.0;
    double myosin_length = 1.5;
    double myosin_radius = 0.025;
    double am_cutoff = 0.020;
    double am_optimal = 0.055;
    double aa_cutoff = 0.08;
    double aa_optimal = 0.05;
    double k_on = 0.0;
    double k_off = 0.0;
    double base_lifetime = 0.0;
    double directional_lifetime = 0.0;
    double lifetime_coeff = 0.0;
    double diff_ratio = 1.0;
    double zero = 0.0;
    std::string filename = "/tmp/am_geometry_regression_test.h5";
    int seed = 991;
    int fix_myosin = 0;
    double dt = 1e-3;
    bool directional = true;

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, seed);
    Sarcomere model(
        n_actins, n_myosins, box, actin_length, myosin_length,
        myosin_radius, am_cutoff, am_optimal, aa_cutoff, aa_optimal,
        k_on, k_off, base_lifetime, directional_lifetime, lifetime_coeff,
        diff_ratio, zero, zero, zero, zero, zero, zero,
        filename, rng, seed, fix_myosin, dt, 0.0, 0.0, 0.0,
        directional, 24, zero, zero, zero, zero);
    model.sterics_enabled = false;
    model.set_am_kinetics(
        true, 1e12, 0.0, 0.0, 0.0, 1.0, 1000.0, 0.005,
        4.5, 5.0, 0.1, 0.12, 0.0, 0.02, 0.2, 1.0, 0.0);
    model.set_am_thin_activation(
        false, false, 0.0, 0.5, 3.0, 0.0, 0.02, 2.0, 0.0);
    model.set_am_thick_regulation(false, 0.0, 0.0, 1.0, 1.0, -1.0);
    model.set_am_head_count(12);
    // This test isolates axial retention and detachment position. Disable the
    // separate lattice-neighbor sharing constraint so crown availability does
    // not change which final axial head is being tracked.
    model.set_am_neighbor_head_cap(0);
    model.set_am_proximity_tether(0.0, 0.005, false, 180.0);
    model.myosin.center[0] = {0.0, 0.0, 0.0};
    model.myosin.direction[0] = {1.0, 0.0, 0.0};

    // Full end-to-end KMC truth table. Actins 0 and 2 occupy the productive
    // halves; 1 and 3 have the same respective polarities on the wrong half.
    const double radial = 0.055;
    model.actin.center[0] = {-0.50, radial, 0.0};
    model.actin.direction[0] = {1.0, 0.0, 0.0};
    model.actin.center[1] = {0.50, radial, 0.0};
    model.actin.direction[1] = {1.0, 0.0, 0.0};
    model.actin.center[2] = {0.50, radial, 0.0};
    model.actin.direction[2] = {-1.0, 0.0, 0.0};
    model.actin.center[3] = {-0.50, radial, 0.0};
    model.actin.direction[3] = {-1.0, 0.0, 0.0};
    model.actin.update_endpoints();
    model.myosin.update_endpoints();
    model.update_system();
    assert(occupied(model, 0, 0) > 0);
    assert(occupied(model, 1, 2) > 0);
    assert(occupied(model, 0, 1) == 0);
    assert(occupied(model, 1, 1) == 0);
    assert(occupied(model, 0, 3) == 0);
    assert(occupied(model, 1, 3) == 0);

    // Existing attachments survive outside the association shell but inside
    // the explicit rupture radius. A detached head cannot newly associate at
    // the same radial position.
    for (int i = 1; i < n_actins; ++i) {
        model.actin.center[i] = {0.0, 1.5, 0.0};
    }
    detach_all(model);
    model.actin.center[0] = {-0.50, radial, 0.0};
    model.actin.direction[0] = {1.0, 0.0, 0.0};
    model.actin.update_endpoints();
    model.update_system();
    assert(occupied(model, 0, 0) > 0);
    model.actin.center[0].y = 0.090;
    int actin_zero = 0;
    model.actin.update_endpoints(actin_zero);
    const double capture_centerline = model.myosin.radius +
        model.actin_radius + model.am_surface_gap + am_cutoff;
    const auto outside_capture = geometry::analyze_am(
        model.actin.left_end[0], model.actin.right_end[0],
        model.myosin.left_end[0], model.myosin.right_end[0],
        capture_centerline, box, model.is_periodic, true);
    assert(outside_capture.myosin_binding_ratio == 0.0);
    assert(outside_capture.partial_binding_ratio == 0.0);
    model.update_system();
    assert(occupied(model, 0, 0) > 0);
    detach_all(model);
    model.update_system();
    assert(occupied(model, 0) == 0);
    assert(occupied(model, 1) == 0);

    // Axial passage through the bare zone and onto the wrong half. Repeat for
    // both polarities, multiple time steps, and two neighbor-cell placements.
    const std::array<double, 3> time_steps{0.02, 0.01, 0.005};
    const std::array<double, 2> shifts{0.0, 1.0};
    const double speed = 1.0;
    const double bare_half =
        0.5 * geometry::am_bare_zone_fraction * myosin_length;
    const double analytic_magnitude = 0.5 * actin_length - bare_half;
    std::array<double, 3> baseline_positive{};
    for (int polarity_index = 0; polarity_index < 2; ++polarity_index) {
        const double polarity = polarity_index == 0 ? 1.0 : -1.0;
        const int correct_half = polarity > 0.0 ? 0 : 1;
        const int wrong_half = 1 - correct_half;
        for (int dt_index = 0; dt_index < 3; ++dt_index) {
            for (int shift_index = 0; shift_index < 2; ++shift_index) {
                const double local_dt = time_steps[dt_index];
                const double shift = shifts[shift_index];
                model.dt = local_dt;
                model.myosin.center[0] = {shift, 0.0, 0.0};
                model.myosin.direction[0] = {1.0, 0.0, 0.0};
                for (int i = 1; i < n_actins; ++i) {
                    model.actin.center[i] = {0.0, 1.5, 0.0};
                }
                model.actin.direction[0] = {polarity, 0.0, 0.0};
                model.actin.center[0] = {
                    shift - 0.50 * polarity, radial, 0.0};
                model.actin.update_endpoints();
                model.myosin.update_endpoints();
                // Resolve the physical 14.3-nm crown spacing. A six-head test
                // pool has only two axial crowns and necessarily loses the
                // final sliver of overlap too early.
                model.set_am_head_count(141);
                model.am_forced_rupture_events.clear();
                model.update_system();
                assert(occupied(model, correct_half, 0) > 0);
                assert(occupied(model, wrong_half) == 0);

                const double increment = polarity * speed * local_dt;
                double detached_relative = 0.0;
                bool detached = false;
                for (int step = 0; step < 300; ++step) {
                    model.actin.center[0].x += increment;
                    model.actin.update_endpoints(actin_zero);
                    model.update_system();
                    if (occupied(model, correct_half, 0) == 0) {
                        detached_relative =
                            model.actin.center[0].x - model.myosin.center[0].x;
                        detached = true;
                        break;
                    }
                }
                assert(detached);
                const double expected = polarity * analytic_magnitude;
                if (!(std::abs(detached_relative - expected) <=
                      speed * local_dt + 2e-6)) {
                    std::cerr << "axial detach mismatch polarity=" << polarity
                              << " dt=" << local_dt
                              << " shift=" << shift
                              << " detached=" << detached_relative
                              << " expected=" << expected << '\n';
                }
                assert(std::abs(detached_relative - expected) <=
                       speed * local_dt + 2e-6);
                assert(occupied(model, wrong_half) == 0);
                for (int extra = 0; extra < 5; ++extra) model.update_system();
                assert(occupied(model, wrong_half) == 0);
                assert(!model.am_forced_rupture_events.empty());
                const size_t row = model.am_forced_rupture_events.size() - 15;
                const int mask = static_cast<int>(
                    model.am_forced_rupture_events[row + 6]);
                assert(mask & Sarcomere::am_no_productive_axial_overlap);
                assert(!(mask & Sarcomere::am_not_in_neighbor_list));
                assert(!(mask & Sarcomere::am_unexpected_candidate_loss));
                assert(model.am_forced_rupture_events_v2.size() >= 19);
                const size_t row_v2 =
                    model.am_forced_rupture_events_v2.size() - 19;
                assert(static_cast<int>(
                           model.am_forced_rupture_events_v2[row_v2 + 6]) == mask);
                assert(model.am_forced_rupture_events_v2[row_v2 + 15] >= 0.0);
                assert(model.am_forced_rupture_events_v2[row_v2 + 16] >= 0.0);
                assert(model.am_forced_rupture_events_v2[row_v2 + 18] ==
                       model.am_rupture_distance);
                if (shift_index == 0 && polarity > 0.0) {
                    baseline_positive[dt_index] = detached_relative;
                } else if (shift_index == 1 && polarity > 0.0) {
                    assert(std::abs(detached_relative -
                                    baseline_positive[dt_index]) < 1e-10);
                }
            }
        }
    }
    // Smaller steps cannot increase the axial-exit discretization error.
    assert(std::abs(baseline_positive[2] - analytic_magnitude) <=
           std::abs(baseline_positive[0] - analytic_magnitude) + 1e-10);

    gsl_rng_free(rng);
    return 0;
}
