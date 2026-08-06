#include "sarcomere.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <set>
#include <string>

namespace {

double axial_coordinate(const Sarcomere& model, int actin, int myosin) {
    vec displacement = model.actin.center[actin] - model.myosin.center[myosin];
    displacement.pbc_wrap(model.box, model.is_periodic);
    return displacement.dot(model.actin.direction[actin]);
}

Sarcomere::AACrosslink* active_link(Sarcomere& model) {
    for (auto& link : model.aa_links) {
        if (link.active) return &link;
    }
    return nullptr;
}

void seed_relaxed_strong_head(Sarcomere& model, int myosin_index,
                              int half_index, int actin_index) {
    auto& head = model.am_heads_by_myosin[myosin_index][half_index][0];
    head.binding_azimuth = head.anchor_azimuth;
    const auto geometry = compute_am_head_tether(
        model.actin, model.myosin, actin_index, myosin_index, model.box,
        head.anchor_axial, 0.0,
        model.myosin.radius + model.actin_radius + model.am_surface_gap,
        head.binding_azimuth,
        model.am_spatial_heads ? model.myosin.radius : 0.0,
        model.am_spatial_heads ? &model.myosin_material_e1[myosin_index]
                               : nullptr,
        model.am_spatial_heads ? &model.myosin_material_e2[myosin_index]
                               : nullptr,
        model.actin_radius);
    const double fraction = std::clamp(
        geometry.actin_coordinate / model.actin.length + 0.5, 0.0,
        std::nextafter(1.0, 0.0));
    const int site = std::clamp(
        static_cast<int>(std::floor(fraction * model.am_sites_per_actin)),
        0, model.am_sites_per_actin - 1);
    head.target_axial = -0.5 * model.actin.length +
        (static_cast<double>(site) + 0.5) * model.actin.length /
            static_cast<double>(model.am_sites_per_actin);
    head.bound_actin = actin_index;
    head.state = static_cast<Sarcomere::BondFlag>(
        am_motor_kinetics::State::strong);
    head.stroke_progress = 0.005;
}

}  // namespace

int main() {
    int n_actins = 2;
    int n_myosins = 2;
    std::vector<double> box{5.0, 3.0, 3.0};
    double actin_length = 1.0;
    double myosin_length = 1.5;
    double myosin_radius = 0.025;
    double am_cutoff = 0.10;
    double am_optimal = 0.055;
    double aa_cutoff = 0.08;
    double aa_optimal = 0.0506;
    double k_on = 0.0;
    double k_off = 0.0;
    double base_lifetime = 1000.0;
    double directional_lifetime = 0.0;
    double lifetime_coeff = 1.0;
    double diff_ratio = 1.0;
    double k_aa = 3000.0;
    double zero = 0.0;
    std::string filename = "/tmp/actn2_motor_removal_test.h5";
    int seed = 1771;
    int fix_myosin = 0;
    double dt = 1e-4;
    bool directional = true;

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, seed);
    Sarcomere model(
        n_actins, n_myosins, box, actin_length, myosin_length,
        myosin_radius, am_cutoff, am_optimal, aa_cutoff, aa_optimal,
        k_on, k_off, base_lifetime, directional_lifetime, lifetime_coeff,
        diff_ratio, k_aa, zero, zero, zero, zero, zero,
        filename, rng, seed, fix_myosin, dt, 0.0, 0.0, 0.0,
        directional, 8, zero, zero, zero, zero);
    model.sterics_enabled = false;
    model.set_local_aa_model(true, 4, 0.0192, 0.1, 50.0, 0);
    model.set_aa_axial_overlap(0.1);
    model.set_am_radial_geometry(0.015, 0.015, 0.0);
    model.set_am_kinetics(
        true, 0.0, 0.0, 0.0, 0.0, 1.0, 1000.0, 0.005,
        4.5, 5.0, 0.1, 0.12, 0.0, 0.02, 0.2, 1.0, 0.0);
    model.set_am_thin_activation(
        false, false, 0.0, 0.5, 3.0, 0.0, 0.02, 2.0, 0.0);
    model.set_am_thick_regulation(false, 0.0, 0.0, 1.0, 1.0, -1.0);
    model.set_am_proximity_tether(0.0, 0.005, false, 180.0);
    model.set_am_head_count(1);
    model.set_aa_status_kinetics(
        actn2_load_kinetics::StatusMode::arm_load_productive, 0.25, 0.125);
    model.minimal_catch();

    std::set<double> first_arm_sites;
    for (const auto& candidate : model.aa_links) {
        if (candidate.active) first_arm_sites.insert(candidate.s_i);
    }
    assert(first_arm_sites.size() == 4);

    auto* link = active_link(model);
    assert(link != nullptr);
    const std::uint64_t link_id = link->id;
    assert(link->status == 1);

    // Actin 0 (+ polarity) is driven by the negative half of myosin 0;
    // actin 1 (- polarity) is driven by the positive half of myosin 1.
    seed_relaxed_strong_head(model, 0, 0, 0);
    seed_relaxed_strong_head(model, 1, 1, 1);
    auto& second_head = model.am_heads_by_myosin[1][1][0];

    // Emulate the small outward displacement produced by the two strong heads.
    // Status 2 must respond to the force actually carried by both ACTN2 arms,
    // not to motor occupancy alone.
    model.actin.center[0].x += 0.002;
    model.actin.center[1].x -= 0.002;
    model.actin.update_endpoints();

    model.update_system();
    link = active_link(model);
    assert(link != nullptr && link->id == link_id);
    assert(link->status == 2);
    assert(model.actin["am_active_force_x"][0] > 0.0);
    assert(model.actin["am_active_force_x"][1] < 0.0);
    const vec total_force = model.actin.force[0] + model.actin.force[1] +
                            model.myosin.force[0] + model.myosin.force[1];
    assert(total_force.norm() < 1e-8);

    // Removing one motor eliminates bilateral productive load. The ACTN2
    // molecule remains chemically bound and demotes from status 2 to status 1.
    second_head.bound_actin = -1;
    second_head.state = static_cast<Sarcomere::BondFlag>(
        am_motor_kinetics::State::detached);
    model.update_system();
    link = active_link(model);
    assert(link != nullptr && link->id == link_id);
    assert(link->active);
    assert(link->status == 1);

    gsl_rng_free(rng);
    return 0;
}
