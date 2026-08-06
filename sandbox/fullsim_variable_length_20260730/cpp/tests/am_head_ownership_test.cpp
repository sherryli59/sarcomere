#include "sarcomere.h"

#include <cassert>
#include <cmath>
#include <set>
#include <string>

int main() {
    int n_actins = 6;
    int n_myosins = 2;
    std::vector<double> box{4.0, 4.0, 4.0};
    double actin_length = 1.0;
    double myosin_length = 1.5;
    double myosin_radius = 0.025;
    double am_cutoff = 0.1;
    double am_optimal = 0.05;
    double aa_cutoff = 0.1;
    double aa_optimal = 0.05;
    double k_on = 0.0;
    double k_off = 0.0;
    double base_lifetime = 0.0;
    double directional_lifetime = 0.0;
    double lifetime_coeff = 0.0;
    double diff_ratio = 1.0;
    double k_aa = 10.0;
    double kappa_aa = 0.0;
    double k_am = 10.0;
    double kappa_am = 0.0;
    double k_mm = 0.0;
    double v_am = 0.0;
    std::string filename = "/tmp/am_head_ownership_test.h5";
    int seed = 17;
    int fix_myosin = 0;
    double dt = 1e-5;
    double tau_rec = 0.0;
    double titin_k = 0.0;
    double titin_rest = 0.0;
    bool directional = true;
    int max_myosin_bonds = 6;
    double no_cap = 0.0;

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, seed);
    Sarcomere model(
        n_actins, n_myosins, box, actin_length, myosin_length,
        myosin_radius, am_cutoff, am_optimal, aa_cutoff, aa_optimal,
        k_on, k_off, base_lifetime, directional_lifetime, lifetime_coeff,
        diff_ratio, k_aa, kappa_aa, k_am, kappa_am, k_mm, v_am,
        filename, rng, seed, fix_myosin, dt, tau_rec, titin_k, titin_rest,
        directional, max_myosin_bonds, no_cap, no_cap, no_cap, no_cap);

    // Scaffold-free sterics must not depend on an ACTN2 bond flag.
    for (int a = 0; a < n_actins; ++a) {
        model.actin.center[a] = {0.0, 0.5 + 0.3 * a, 0.0};
        model.actin.direction[a] = {1.0, 0.0, 0.0};
        model.actin.cb_status[a] = 0;
    }
    model.actin.center[0] = {0.0, 0.0, 0.0};
    model.actin.center[1] = {0.0, 0.01, 0.0};
    for (int m = 0; m < n_myosins; ++m) {
        model.myosin.center[m] = {1.5, 1.5, 1.5};
    }
    model.actin.update_endpoints();
    model.myosin.update_endpoints();
    model.update_system();
    const vec force0 = model.actin.force[0];
    const vec force1 = model.actin.force[1];
    assert(force0.norm() > 0.0);
    assert((force0 + force1).norm() < 1e-10);

    constexpr int heads_per_half = 7;
    model.set_am_head_count(heads_per_half);
    assert(model.am_heads_per_half == heads_per_half);
    assert(model.am_heads_by_myosin.size() == static_cast<size_t>(n_myosins));

    for (const auto& myosin_halves : model.am_heads_by_myosin) {
        int total = 0;
        for (int half_index = 0; half_index < 2; ++half_index) {
            const auto& heads = myosin_halves[static_cast<size_t>(half_index)];
            assert(heads.size() == heads_per_half);
            total += static_cast<int>(heads.size());
            std::set<double> azimuths;
            for (const auto& head : heads) {
                assert(head.bound_actin == -1);
                assert(head.state == static_cast<Sarcomere::BondFlag>(
                    am_motor_kinetics::State::detached));
                assert((half_index == 0 && head.half == -1) ||
                       (half_index == 1 && head.half == 1));
                assert((half_index == 0 && head.anchor_axial < 0.0) ||
                       (half_index == 1 && head.anchor_axial > 0.0));
                const double abs_axial = std::abs(head.anchor_axial);
                assert(abs_axial >
                       0.5 * geometry::am_bare_zone_fraction * myosin_length);
                assert(abs_axial < myosin_length / 2.0);
                azimuths.insert(head.anchor_azimuth);
            }
            assert(azimuths.size() >= 3);
            // The first physical crown is a balanced 120-degree triplet at
            // one axial material coordinate.
            assert(std::abs(heads[0].anchor_axial -
                            heads[1].anchor_axial) < 1e-12);
            assert(std::abs(heads[1].anchor_axial -
                            heads[2].anchor_axial) < 1e-12);
            assert(std::abs(std::remainder(
                heads[1].anchor_azimuth - heads[0].anchor_azimuth,
                2.0 * M_PI) - 2.0 * M_PI / 3.0) < 1e-12);
            assert(std::abs(std::remainder(
                heads[2].anchor_azimuth - heads[1].anchor_azimuth,
                2.0 * M_PI) - 2.0 * M_PI / 3.0) < 1e-12);
        }
        assert(total == 2 * heads_per_half);
    }

    // Several actins can share a thick filament, but each physical head owns
    // exactly one target and the pool size does not change with target count.
    auto& positive_heads = model.am_heads_by_myosin[0][1];
    for (int h = 0; h < heads_per_half; ++h) {
        positive_heads[static_cast<size_t>(h)].bound_actin = h % n_actins;
        positive_heads[static_cast<size_t>(h)].state =
            static_cast<Sarcomere::BondFlag>(am_motor_kinetics::State::weak);
    }
    std::set<int> targets;
    for (const auto& head : positive_heads) targets.insert(head.bound_actin);
    assert(targets.size() == static_cast<size_t>(n_actins));
    assert(positive_heads.size() == heads_per_half);

    model.set_am_head_count(3);

    // Directional overlap, rather than the actin center, decides which
    // bipolar half can supply heads.  On the positive myosin half only an
    // oppositely oriented actin is productive.
    const vec myosin_left{-0.5, 0.0, 0.0};
    const vec myosin_right{0.5, 0.0, 0.0};
    const std::vector<double> test_box{4.0, 4.0, 4.0};
    const auto wrong_polarity = geometry::analyze_am(
        {0.10, 0.05, 0.0}, {0.60, 0.05, 0.0},
        myosin_left, myosin_right, 0.10, test_box, true);
    const auto correct_polarity = geometry::analyze_am(
        {0.60, 0.05, 0.0}, {0.10, 0.05, 0.0},
        myosin_left, myosin_right, 0.10, test_box, true);
    assert(wrong_polarity.myosin_binding_ratio > 0.0);
    assert(wrong_polarity.partial_binding_ratio == 0.0);
    assert(correct_polarity.partial_binding_ratio > 0.0);

    // This actin's center lies on the negative side of the M-line, but its
    // correctly polarized tip still overlaps the positive head-bearing half.
    // A center-based owner would incorrectly give it to the negative half.
    const auto crossed_center = geometry::analyze_am(
        {0.20, 0.05, 0.0}, {-0.40, 0.05, 0.0},
        myosin_left, myosin_right, 0.10, test_box, true);
    assert(crossed_center.partial_binding_ratio > 0.0);

    // The central 10% bare zone geometrically overlaps the thick backbone but
    // has no motor-bearing overlap.
    const auto bare_zone_only = geometry::analyze_am(
        {-0.04, 0.05, 0.0}, {0.04, 0.05, 0.0},
        myosin_left, myosin_right, 0.10, test_box, true);
    assert(bare_zone_only.myosin_binding_ratio > 0.0);
    assert(bare_zone_only.partial_binding_ratio == 0.0);
    for (const auto& myosin_halves : model.am_heads_by_myosin) {
        assert(myosin_halves[0].size() == 3);
        assert(myosin_halves[1].size() == 3);
        for (const auto& half_heads : myosin_halves) {
            for (const auto& head : half_heads) assert(head.bound_actin == -1);
        }
    }

    // The material crown retains axial roll, remains orthonormal after tilt,
    // and does not drift under repeated rotations.
    model.myosin.direction[0] = {1.0, 0.0, 0.0};
    model.reset_myosin_material_frames();
    const vec initial_e1 = model.myosin_material_e1[0];
    const vec initial_e2 = model.myosin_material_e2[0];
    model.transport_myosin_material_frame(
        0, {1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, 0.5 * M_PI);
    assert(model.myosin_material_e1[0].dot(initial_e2) > 1.0 - 1e-10);
    assert(std::abs(model.myosin_material_e1[0].dot(initial_e1)) < 1e-10);
    vec tilted{1.0, 1.0, 0.25};
    tilted.normalize();
    model.transport_myosin_material_frame(
        0, {1.0, 0.0, 0.0}, tilted, 0.0);
    for (int iteration = 0; iteration < 1000; ++iteration) {
        model.transport_myosin_material_frame(0, tilted, tilted, 1e-3);
    }
    assert(std::abs(model.myosin_material_e1[0].norm() - 1.0) < 1e-10);
    assert(std::abs(model.myosin_material_e2[0].norm() - 1.0) < 1e-10);
    assert(std::abs(model.myosin_material_e1[0].dot(tilted)) < 1e-10);
    assert(std::abs(model.myosin_material_e2[0].dot(tilted)) < 1e-10);
    assert(std::abs(model.myosin_material_e1[0].dot(
        model.myosin_material_e2[0])) < 1e-10);

    // Owner-centric competing hazards distribute a conserved half-filament
    // pool across several simultaneously accessible thin filaments. No AM
    // pair is allowed to instantiate another copy of the pool.
    model.directional = false;
    model.set_am_kinetics(
        true, 1.0e12, 0.0, 0.0, 0.0, 1.0, 10.0, 0.005, 1.0, 5.0,
        0.1, 0.2, 0.0, 0.02, 0.2, 1.0, 0.0);
    model.set_am_head_count(12);
    model.set_am_neighbor_head_cap(2);
    model.myosin.center[0] = {0.0, 0.0, 0.0};
    model.myosin.direction[0] = {1.0, 0.0, 0.0};
    model.reset_myosin_material_frames();
    for (int a = 0; a < n_actins; ++a) {
        const double angle = 2.0 * M_PI * static_cast<double>(a) /
                             static_cast<double>(n_actins);
        model.actin.center[a] =
            {0.3, 0.04 * std::cos(angle), 0.04 * std::sin(angle)};
        model.actin.direction[a] = {1.0, 0.0, 0.0};
    }
    model.actin.update_endpoints();
    model.myosin.update_endpoints();
    model.update_system();

    int bound_count = 0;
    std::set<int> runtime_targets;
    std::vector<int> target_counts(static_cast<size_t>(n_actins), 0);
    for (const auto& head : model.am_heads_by_myosin[0][1]) {
        if (head.bound_actin < 0) continue;
        ++bound_count;
        assert(head.bound_actin < n_actins);
        runtime_targets.insert(head.bound_actin);
        ++target_counts[static_cast<size_t>(head.bound_actin)];
    }
    assert(bound_count == 12);
    assert(runtime_targets.size() == static_cast<size_t>(n_actins));
    for (int count : target_counts) assert(count == 2);
    for (const auto& head : model.am_heads_by_myosin[0][0]) {
        assert(head.bound_actin == -1);
    }

    // Strong heads reserve capacity. Deliberately overfill actin 0 with two
    // additional weak heads; the next update must retain its two strong heads
    // and evict/reassign only the weak excess.
    auto& capped_heads = model.am_heads_by_myosin[0][1];
    std::vector<int> on_zero;
    std::vector<int> elsewhere;
    for (int h = 0; h < static_cast<int>(capped_heads.size()); ++h) {
        if (capped_heads[static_cast<size_t>(h)].bound_actin == 0) {
            on_zero.push_back(h);
        } else {
            elsewhere.push_back(h);
        }
    }
    assert(on_zero.size() == 2 && elsewhere.size() >= 2);
    for (int h : on_zero) {
        capped_heads[static_cast<size_t>(h)].state =
            static_cast<Sarcomere::BondFlag>(am_motor_kinetics::State::strong);
    }
    for (int k = 0; k < 2; ++k) {
        auto& head = capped_heads[static_cast<size_t>(elsewhere[k])];
        head.bound_actin = 0;
        head.state = static_cast<Sarcomere::BondFlag>(
            am_motor_kinetics::State::weak);
    }
    model.update_system();
    int zero_strong = 0;
    int zero_weak = 0;
    for (const auto& head : capped_heads) {
        if (head.bound_actin != 0) continue;
        if (head.state == static_cast<Sarcomere::BondFlag>(
                am_motor_kinetics::State::strong)) ++zero_strong;
        if (head.state == static_cast<Sarcomere::BondFlag>(
                am_motor_kinetics::State::weak)) ++zero_weak;
    }
    assert(zero_strong == 2);
    assert(zero_weak == 0);

    gsl_rng_free(rng);
    return 0;
}
