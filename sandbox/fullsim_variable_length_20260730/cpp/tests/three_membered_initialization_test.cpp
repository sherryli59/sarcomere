#include "axial_cell.h"
#include "sarcomere.h"

#include <array>
#include <cassert>
#include <cmath>
#include <set>
#include <string>
#include <utility>

namespace {

int occupied_heads(const Sarcomere& model, int myosin, int half,
                   int actin) {
    int count = 0;
    for (const auto& head : model.am_heads_by_myosin[myosin][half]) {
        if (head.bound_actin == actin &&
            head.state == static_cast<Sarcomere::BondFlag>(
                am_motor_kinetics::State::weak)) {
            ++count;
        }
    }
    return count;
}

double transverse_distance(const vec& first, const vec& second) {
    const double dy = first.y - second.y;
    const double dz = first.z - second.z;
    return std::sqrt(dy * dy + dz * dz);
}

}  // namespace

int main() {
    constexpr double d_am = 0.055;
    int n_actins = 12;
    int n_myosins = 2;
    std::vector<double> box{5.32, 3.0, 3.0};
    double actin_length = 1.0;
    double myosin_length = 1.5;
    double myosin_radius = 0.025;
    double am_cutoff = 0.05;
    double am_optimal = d_am;
    double aa_cutoff = 0.0706;
    double aa_optimal = 0.0506;
    double k_on = 500.0;
    double k_off = 1.0;
    double base_lifetime = 0.05;
    double directional_lifetime = 0.05;
    double lifetime_coeff = 0.01;
    double diff_ratio = 1.0;
    double k_aa = 3000.0;
    double kappa_aa = 50.0;
    double zero = 0.0;
    std::string filename = "/tmp/three_membered_initialization_test.h5";
    int seed = 4517;
    int fixed_myosins = 0;
    double dt = 1e-6;
    bool directional = true;

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, seed);
    Sarcomere model(
        n_actins, n_myosins, box, actin_length, myosin_length,
        myosin_radius, am_cutoff, am_optimal, aa_cutoff, aa_optimal,
        k_on, k_off, base_lifetime, directional_lifetime,
        lifetime_coeff, diff_ratio, k_aa, kappa_aa,
        zero, zero, zero, zero,
        filename, rng, seed, fixed_myosins, dt, 0.002, zero, zero,
        directional, 2, zero, zero, zero, zero);
    model.set_periodicity({true, false, false});
    model.set_local_aa_model(true, 2, 0.0192, 0.1, 50.0, 0);
    model.set_aa_axial_overlap(0.1);
    model.aa_barbed_zone = 0.10;
    model.set_am_radial_geometry(0.015, 0.015, 1000.0);
    model.set_am_kinetics(
        true, 1000.0, 200.0, 200.0, 80.0, 1.0,
        1000.0, 0.005, 4.5, 5.0, 0.1, 0.1,
        0.05, 0.005, 0.03, 0.0, 0.0005);
    model.set_am_proximity_tether(1000.0, 0.005, true, 180.0);
    model.set_am_head_count(35);
    model.init_aa_overlap = 0.1;
    model.three_membered_lattice();
    model.seed_three_membered_heads(6);

    // Every thick-filament half sees exactly three thin filaments on a ring
    // of radius d_AM, separated by 120 degrees.
    const std::array<std::array<int, 5>, 4> populations{{
        {{0, 0, 9, 10, 11}}, {{0, 1, 0, 1, 2}},
        {{1, 0, 3, 4, 5}}, {{1, 1, 6, 7, 8}},
    }};
    for (const auto& population : populations) {
        const int myosin = population[0];
        const int half = population[1];
        for (int offset = 0; offset < 3; ++offset) {
            const int actin = population[2 + offset];
            assert(std::abs(transverse_distance(
                model.actin.center[actin], model.myosin.center[myosin]) -
                d_am) < 1e-10);
            assert(occupied_heads(model, myosin, half, actin) == 6);
        }
        for (int offset = 0; offset < 3; ++offset) {
            const vec first = model.actin.center[population[2 + offset]];
            const vec second =
                model.actin.center[population[2 + (offset + 1) % 3]];
            assert(std::abs(transverse_distance(first, second) -
                            std::sqrt(3.0) * d_am) < 1e-10);
        }
    }

    // The alternating-polarity populations form one closed six-edge ACTN2
    // ring at each Z-disc. Each actin has two distinct ring neighbours.
    std::set<std::pair<int, int>> pairs;
    std::array<int, 12> degree{};
    int active_links = 0;
    for (const auto& link : model.aa_links) {
        if (!link.active) continue;
        ++active_links;
        const auto pair = std::minmax(link.actin_i, link.actin_j);
        if (pairs.insert(pair).second) {
            ++degree[static_cast<size_t>(pair.first)];
            ++degree[static_cast<size_t>(pair.second)];
        }
    }
    assert(active_links == 24);
    assert(pairs.size() == 12);
    for (int value : degree) assert(value == 2);

    // Reinitialize the same populations near contact but without a perfect
    // ring, prebound ACTN2, or seeded heads. This is the emergent-assembly
    // starting condition used by the explicit-titin diagnostic.
    for (auto& link : model.aa_links) link.active = false;
    for (auto& myosin_heads : model.am_heads_by_myosin) {
        for (auto& half_heads : myosin_heads) {
            for (auto& head : half_heads) {
                head.bound_actin = -1;
                head.state = static_cast<Sarcomere::BondFlag>(
                    am_motor_kinetics::State::detached);
            }
        }
    }
    model.three_membered_lattice(false, 0.01, 0.02, 0.015);
    for (const auto& link : model.aa_links) assert(!link.active);
    assert(std::abs(transverse_distance(
        model.actin.center[0], model.myosin.center[0]) - d_am) > 1e-6);

    // Twelve independent titins (three per thick half) capture twelve
    // individual barbed-end domains. No head occupancy or ACTN2 component is
    // available, so successful binding proves the old Z-node gate is absent.
    axial_cell::Config titin;
    titin.mode = axial_cell::Mode::fixed;
    titin.apply_internal_forces_when_fixed = true;
    titin.explicit_titin = true;
    titin.explicit_titin_per_half = 3;
    titin.explicit_titin_k = 220.0;
    titin.explicit_titin_slack_length = 0.375;
    titin.explicit_titin_capture_distance = 0.75;
    titin.explicit_titin_k_on = 1.0e12;
    // This block isolates topology/ownership rather than physical kinetics.
    titin.explicit_titin_strain_gated_kinetics = false;
    axial_cell::initialize_explicit_titins(model, titin);
    assert(titin.explicit_titins.size() == 12);
    for (int actin = 0; actin < model.actin.n; ++actin) {
        model.actin.force[actin] = {0.0, 0.0, 0.0};
        model.actin.torque[actin] = {0.0, 0.0, 0.0};
    }
    for (int myosin = 0; myosin < model.myosin.n; ++myosin) {
        model.myosin.force[myosin] = {0.0, 0.0, 0.0};
        model.myosin.torque[myosin] = {0.0, 0.0, 0.0};
    }
    const auto capture = axial_cell::prepare_step(model, titin);
    assert(capture.titin_candidate_links == 12);
    assert(capture.titin_bound_links == 12);
    std::set<int> bound_actins;
    for (const auto& molecule : titin.explicit_titins) {
        assert(molecule.bound);
        assert(molecule.bound_actin >= 0);
        bound_actins.insert(molecule.bound_actin);
    }
    assert(bound_actins.size() == 12);
    vec total_force{0.0, 0.0, 0.0};
    bool found_endpoint_torque = false;
    for (int actin = 0; actin < model.actin.n; ++actin) {
        total_force += model.actin.force[actin];
        found_endpoint_torque = found_endpoint_torque ||
            static_cast<vec>(model.actin.torque[actin]).norm() > 1e-12;
    }
    for (int myosin = 0; myosin < model.myosin.n; ++myosin) {
        total_force += model.myosin.force[myosin];
    }
    assert(total_force.norm() < 1e-10);
    assert(found_endpoint_torque);

    gsl_rng_free(rng);
    return 0;
}
