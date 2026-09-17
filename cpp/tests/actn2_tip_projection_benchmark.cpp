#include "langevin.h"
#include "sarcomere.h"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace {

struct Fixture {
    std::unique_ptr<Sarcomere> model;
    std::unique_ptr<Langevin> integrator;
    std::vector<utils::vec> actin_centers;
    std::vector<utils::vec> myosin_centers;
};

Fixture make_fixture(const std::string& output, int pair_count,
                     int links_per_pair, int penetrating_stride,
                     gsl_rng* rng, int seed) {
    int n_actins = 2 * pair_count;
    int n_myosins = pair_count;
    std::vector<double> box{
        std::max(4.0, 2.0 * static_cast<double>(pair_count) + 2.0), 4.0, 4.0};
    double actin_length = 1.0;
    double myosin_length = 1.5;
    double myosin_radius = 0.025;
    double am_cutoff = 0.05;
    double am_optimal = 0.055;
    double aa_cutoff = 0.10;
    double aa_optimal = 0.05;
    double zero = 0.0;
    double k_on = 0.0;
    double k_off = 0.0;
    double base_lifetime = 1e9;
    double directional_base_lifetime = 0.0;
    double lifetime_coeff = 0.0;
    double diff_ratio = 1.0;
    int fixed_myosins = 0;
    double dt = 1e-5;
    bool directional = true;
    int max_myosin_bonds = 24;
    std::string filename = output;
    auto model = std::make_unique<Sarcomere>(
        n_actins, n_myosins, box, actin_length, myosin_length,
        myosin_radius, am_cutoff, am_optimal, aa_cutoff, aa_optimal,
        k_on, k_off, base_lifetime, directional_base_lifetime,
        lifetime_coeff, diff_ratio, zero, zero, zero, zero, zero, zero,
        filename, rng, seed, fixed_myosins, dt, zero, zero, zero,
        directional, max_myosin_bonds, zero, zero, zero, zero);
    model->set_periodicity({false, false, false});
    model->dynamic_am = false;
    model->sterics_enabled = false;
    model->aa_local_multilink = true;
    model->aa_arm_state_model = false;
    model->aa_zdisc_tip_contact = true;
    model->aa_zdisc_tip_onset = 0.050;
    model->aa_zdisc_tip_stiffness = 1000.0;
    model->aa_zdisc_tip_force_max = 20.0;
    model->aa_zdisc_tip_transverse_flat = 0.100;
    model->aa_zdisc_tip_transverse_cutoff = 0.130;
    model->aa_links.clear();
    model->aa_active_link_slots.clear();
    model->aa_links_by_actin.assign(static_cast<size_t>(n_actins), {});
    model->actin_neighbors_by_species.resize(static_cast<size_t>(n_actins));

    std::vector<utils::vec> initial_actin(static_cast<size_t>(n_actins));
    std::vector<utils::vec> initial_myosin(static_cast<size_t>(n_myosins));
    for (int pair = 0; pair < pair_count; ++pair) {
        const int first = 2 * pair;
        const int second = first + 1;
        const double x = 2.0 * static_cast<double>(pair) + 1.0;
        model->actin.center[first] = {x, 0.0, 0.0};
        model->actin.center[second] = {x, aa_optimal, 0.0};
        model->actin.direction[first] = {1.0, 0.0, 0.0};
        model->actin.direction[second] = {-1.0, 0.0, 0.0};
        const bool penetrating = penetrating_stride > 0 &&
            pair % penetrating_stride == 0;
        model->myosin.center[pair] = {
            x, penetrating ? -0.715 : -1.0, 0.0};
        model->myosin.direction[pair] = {0.0, 1.0, 0.0};
        model->actin_neighbors_by_species[static_cast<size_t>(first)].second =
            {pair};
        model->actin_neighbors_by_species[static_cast<size_t>(second)].second =
            {pair};
        initial_actin[static_cast<size_t>(first)] = model->actin.center[first];
        initial_actin[static_cast<size_t>(second)] = model->actin.center[second];
        initial_myosin[static_cast<size_t>(pair)] = model->myosin.center[pair];
        for (int ordinal = 0; ordinal < links_per_pair; ++ordinal) {
            Sarcomere::AACrosslink link;
            link.id = model->aa_next_link_id++;
            link.actin_i = first;
            link.actin_j = second;
            link.s_i = 0.02 * static_cast<double>(
                ordinal - (links_per_pair - 1) / 2.0);
            link.s_j = -link.s_i;
            link.status = 1;
            link.active = true;
            const size_t slot = model->aa_links.size();
            model->aa_links.push_back(link);
            model->aa_active_link_slots.push_back(slot);
            model->aa_links_by_actin[static_cast<size_t>(first)].push_back(slot);
            model->aa_links_by_actin[static_cast<size_t>(second)].push_back(slot);
        }
    }
    model->actin.update_endpoints();
    model->myosin.update_endpoints();

    double beta = 1.0;
    double d_actin = 1.0;
    double d_actin_rot = 0.0;
    double d_myosin = 1.0;
    double d_myosin_rot = 0.0;
    int save_every = 1000000;
    bool resume = false;
    std::remove(output.c_str());
    auto integrator = std::make_unique<Langevin>(
        *model, beta, dt, d_actin, d_actin_rot, d_myosin, d_myosin_rot,
        save_every, resume, true, 0.0, 0.0, 0.0, 0.0);
    return {std::move(model), std::move(integrator),
            std::move(initial_actin), std::move(initial_myosin)};
}

void reset_geometry(Fixture& fixture) {
    for (size_t index = 0; index < fixture.actin_centers.size(); ++index) {
        fixture.model->actin.center[index] = fixture.actin_centers[index];
    }
    for (size_t index = 0; index < fixture.myosin_centers.size(); ++index) {
        fixture.model->myosin.center[index] = fixture.myosin_centers[index];
    }
    fixture.model->actin.update_endpoints();
    fixture.model->myosin.update_endpoints();
}

Fixture make_dynamic_discovery_fixture(
    const std::string& output, gsl_rng* rng, int seed) {
    auto fixture = make_fixture(output, 2, 1, 0, rng, seed);
    auto& model = *fixture.model;
    model.aa_links.clear();
    model.aa_active_link_slots.clear();
    model.aa_links_by_actin.assign(static_cast<size_t>(model.actin.n), {});
    for (auto& neighbors : model.actin_neighbors_by_species) {
        neighbors.second.clear();
    }
    model.actin_neighbors_by_species[0].second = {0, 1};
    model.actin_neighbors_by_species[1].second = {0, 1};
    model.myosin.center[0] = {1.0, -0.10, 0.0};
    model.myosin.center[1] = {1.0, 0.95, 0.0};
    model.myosin.direction[0] = {0.0, 1.0, 0.0};
    model.myosin.direction[1] = {0.0, 1.0, 0.0};
    Sarcomere::AACrosslink link;
    link.id = model.aa_next_link_id++;
    link.actin_i = 0;
    link.actin_j = 1;
    link.s_i = 0.0;
    link.s_j = 0.0;
    link.status = 1;
    link.active = true;
    model.aa_links.push_back(link);
    model.aa_active_link_slots.push_back(0);
    model.aa_links_by_actin[0].push_back(0);
    model.aa_links_by_actin[1].push_back(0);
    model.actin.update_endpoints();
    model.myosin.update_endpoints();
    fixture.actin_centers.resize(model.actin.center.size());
    for (size_t index = 0; index < model.actin.center.size(); ++index) {
        fixture.actin_centers[index] = model.actin.center[index];
    }
    fixture.myosin_centers.resize(model.myosin.center.size());
    for (size_t index = 0; index < model.myosin.center.size(); ++index) {
        fixture.myosin_centers[index] = model.myosin.center[index];
    }
    const auto patch = compute_actn2_patch_geometry(
        model.actin, 0, 1, 0.0, 0.0, model.box);
    assert(!actn2_patch_may_contact_myosin_tip(
        patch, model.myosin, 1, model.box, model.aa_zdisc_tip_onset,
        model.aa_zdisc_tip_transverse_cutoff));
    return fixture;
}

void print_profile(const char* label, const Actn2TipProjectionProfile& p,
                   double outer_seconds) {
    std::cout << label
              << " outer_s=" << outer_seconds
              << " calls=" << p.calls
              << " links=" << p.active_links
              << " pairs=" << p.unique_actin_pairs
              << " constraints=" << p.constraints
              << " broadphase_tests=" << p.broadphase_tests
              << " broadphase_candidates=" << p.broadphase_candidates
              << " solve_evaluations=" << p.solve_evaluations
              << " audit_evaluations=" << p.audit_evaluations
              << " corrections=" << p.corrections
              << " sweeps=" << p.sweeps
              << " exhausted=" << p.exhausted_calls
              << " audited=" << p.audited_calls
              << " certified_without_audit="
              << p.certified_without_audit_calls
              << " build_s=" << p.build_seconds
              << " solve_s=" << p.solve_seconds
              << " audit_s=" << p.audit_seconds
              << " total_s=" << p.total_seconds << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    const int pair_count = argc > 1 ? std::stoi(argv[1]) : 1000;
    const int repeats = argc > 2 ? std::stoi(argv[2]) : 20;
    const int links_per_pair = 3;
    const int penetrating_stride = 20;  // Five percent of pairs need repair.
    assert(pair_count > 0 && repeats > 0);

    gsl_rng* full_rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng* sparse_rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(full_rng, 7101);
    gsl_rng_set(sparse_rng, 7101);
    auto full = make_fixture(
        "/tmp/actn2_tip_projection_full_scan.h5", pair_count,
        links_per_pair, penetrating_stride, full_rng, 7101);
    auto sparse = make_fixture(
        "/tmp/actn2_tip_projection_sparse.h5", pair_count,
        links_per_pair, penetrating_stride, sparse_rng, 7101);
    int fixed_myosins = 0;

    auto run = [&](Fixture& fixture, Actn2TipProjectionTraversal traversal) {
        fixture.integrator->reset_actn2_tip_projection_profile();
        const auto start = std::chrono::steady_clock::now();
        for (int repeat = 0; repeat < repeats; ++repeat) {
            reset_geometry(fixture);
            fixture.integrator->resolve_actn2_tip_hard_contacts(
                fixed_myosins, 64, 1e-9, traversal);
        }
        return std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
    };

    const double full_seconds = run(
        full, Actn2TipProjectionTraversal::full_scan);
    const double sparse_seconds = run(
        sparse, Actn2TipProjectionTraversal::spatial_incidence);
    const auto& full_profile =
        full.integrator->actn2_tip_projection_profile();
    const auto& sparse_profile =
        sparse.integrator->actn2_tip_projection_profile();
    print_profile("full_scan", full_profile, full_seconds);
    print_profile("spatial_incidence", sparse_profile, sparse_seconds);

    assert(full.model->actin.center.size() == sparse.model->actin.center.size());
    assert(full.model->myosin.center.size() == sparse.model->myosin.center.size());
    double max_center_difference = 0.0;
    for (size_t index = 0; index < full.model->actin.center.size(); ++index) {
        max_center_difference = std::max(
            max_center_difference,
            (full.model->actin.center[index] -
             sparse.model->actin.center[index]).norm());
    }
    for (size_t index = 0; index < full.model->myosin.center.size(); ++index) {
        max_center_difference = std::max(
            max_center_difference,
            (full.model->myosin.center[index] -
             sparse.model->myosin.center[index]).norm());
    }
    std::cout << "parity max_center_difference_um=" << max_center_difference
              << " correction_match="
              << (full_profile.corrections == sparse_profile.corrections)
              << " speedup=" << full_seconds / sparse_seconds << '\n';
    assert(max_center_difference <= 1e-15);
    assert(full_profile.corrections == sparse_profile.corrections);
    assert(sparse_profile.solve_evaluations <= full_profile.solve_evaluations);

    gsl_rng* dynamic_full_rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng* dynamic_spatial_rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(dynamic_full_rng, 7102);
    gsl_rng_set(dynamic_spatial_rng, 7102);
    auto dynamic_full = make_dynamic_discovery_fixture(
        "/tmp/actn2_tip_dynamic_full.h5", dynamic_full_rng, 7102);
    auto dynamic_spatial = make_dynamic_discovery_fixture(
        "/tmp/actn2_tip_dynamic_spatial.h5", dynamic_spatial_rng, 7102);
    dynamic_full.integrator->resolve_actn2_tip_hard_contacts(
        fixed_myosins, 64, 1e-9,
        Actn2TipProjectionTraversal::full_scan);
    dynamic_spatial.integrator->resolve_actn2_tip_hard_contacts(
        fixed_myosins, 64, 1e-9,
        Actn2TipProjectionTraversal::spatial_incidence);
    double dynamic_difference = 0.0;
    for (size_t index = 0;
         index < dynamic_full.model->actin.center.size(); ++index) {
        dynamic_difference = std::max(
            dynamic_difference,
            (dynamic_full.model->actin.center[index] -
             dynamic_spatial.model->actin.center[index]).norm());
    }
    for (size_t index = 0;
         index < dynamic_full.model->myosin.center.size(); ++index) {
        dynamic_difference = std::max(
            dynamic_difference,
            (dynamic_full.model->myosin.center[index] -
             dynamic_spatial.model->myosin.center[index]).norm());
    }
    const auto& dynamic_full_profile =
        dynamic_full.integrator->actn2_tip_projection_profile();
    const auto& dynamic_spatial_profile =
        dynamic_spatial.integrator->actn2_tip_projection_profile();
    std::cout << "dynamic_discovery parity_um=" << dynamic_difference
              << " full_corrections=" << dynamic_full_profile.corrections
              << " spatial_corrections="
              << dynamic_spatial_profile.corrections
              << " broadphase_candidates="
              << dynamic_spatial_profile.broadphase_candidates
              << std::endl;
    assert(dynamic_difference <= 1e-15);
    assert(dynamic_full_profile.corrections ==
           dynamic_spatial_profile.corrections);
    assert(dynamic_full_profile.corrections >= 2);
    gsl_rng_free(dynamic_full_rng);
    gsl_rng_free(dynamic_spatial_rng);

    gsl_rng_free(full_rng);
    gsl_rng_free(sparse_rng);
    return 0;
}
