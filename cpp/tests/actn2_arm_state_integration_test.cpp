#include "actn2_arm_kinetics.h"
#include "sarcomere.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>

namespace {
std::unique_ptr<Sarcomere> make_model(
    const std::string& output, gsl_rng* rng, int& seed) {
    const std::vector<double> box{4.0, 4.0, 4.0};
    int n_actins = 2;
    int n_myosins = 1;
    double actin_length = 1.0;
    double myosin_length = 1.5;
    double myosin_radius = 0.05;
    double am_cutoff = 0.05;
    double am_optimal = 0.055;
    double aa_cutoff = 0.08;
    double aa_optimal = 0.05;
    double base_lifetime = 1e9;
    double diff_ratio = 1.0;
    int fixed_myosins = 0;
    double dt = 1e-5;
    bool directional = true;
    double zero = 0.0;
    std::string filename = output;
    auto model = std::make_unique<Sarcomere>(
        n_actins, n_myosins, box, actin_length, myosin_length,
        myosin_radius, am_cutoff, am_optimal, aa_cutoff, aa_optimal,
        zero, zero, base_lifetime, zero, zero, diff_ratio, zero, zero,
        zero, zero, zero, zero, filename, rng, seed, fixed_myosins, dt,
        zero, zero, zero, directional, 24, zero, zero, zero, zero);
    model->set_periodicity({false, false, false});
    model->dynamic_am = false;
    model->aa_local_multilink = true;
    model->aa_prebundle_links_per_pair = 1;
    model->aa_neighbor_identity_enabled = false;
    model->aa_arm_state_model = true;
    model->aa_status2_locked_until_detach = false;
    model->aa_arm_open_rate_weak = 0.0;
    model->aa_arm_open_rate_mature = 1e12;
    model->aa_arm_bell_distance = 0.0;
    model->aa_tether_off_rate = 0.0;
    model->aa_tether_rebind_rate = 0.0;
    model->aa_link_rupture_length = 0.0;
    model->aa_myosin_contact_epsilon = 0.0;
    model->sterics_enabled = false;

    model->actin.center[0] = {0.0, 0.0, 0.0};
    model->actin.direction[0] = {1.0, 0.0, 0.0};
    model->actin.center[1] = {0.0, 0.05, 0.0};
    model->actin.direction[1] = {-1.0, 0.0, 0.0};
    model->actin.update_endpoints();
    model->myosin.center[0] = {1.5, 0.0, 0.0};
    model->myosin.direction[0] = {1.0, 0.0, 0.0};
    model->myosin.update_endpoints();

    model->prebind_local_aa_pair(0, 1, 1);
    auto& link = model->aa_links[model->aa_active_link_slots.front()];
    link.status = static_cast<Sarcomere::BondFlag>(2);
    link.arm_state = actn2_arm_kinetics::State::mature_bridge;
    return model;
}

size_t active_count(const Sarcomere& model) {
    size_t count = 0;
    for (const auto& link : model.aa_links) count += link.active ? 1 : 0;
    return count;
}
}

int main() {
    const std::string output = "/tmp/actn2_arm_state_integration_test.h5";
    std::remove(output.c_str());
    int seed = 9917;
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, seed);
    auto model = make_model(output, rng, seed);

    // Opening a mature bridge releases exactly one arm; it does not erase the
    // molecule or turn the attachment coordinates into a diffusing pair.
    model->update_system();
    assert(active_count(*model) == 1);
    const auto& tether = model->aa_links[model->aa_active_link_slots.front()];
    assert(actn2_arm_kinetics::is_tether(tether.arm_state));
    assert(actn2_arm_kinetics::arm_i_bound(tether.arm_state) !=
           actn2_arm_kinetics::arm_j_bound(tether.arm_state));
    const auto tether_state = tether.arm_state;
    const double bound_site = actn2_arm_kinetics::arm_i_bound(tether_state)
        ? tether.s_i : tether.s_j;

    // With both tether turnover rates disabled, the same material arm remains
    // bound and exerts no two-filament spring force on the next step.
    model->update_system();
    assert(active_count(*model) == 1);
    const auto& persistent =
        model->aa_links[model->aa_active_link_slots.front()];
    assert(persistent.arm_state == tether_state);
    assert((actn2_arm_kinetics::arm_i_bound(tether_state)
                ? persistent.s_i : persistent.s_j) == bound_site);

    // Checkpoint/resume must preserve which physical arm remains bound.
    model->new_file();
    model->save_state();
    model->save_resume_snapshot();
    gsl_rng* restored_rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(restored_rng, seed);
    auto restored = make_model(output, restored_rng, seed);
    int frames = 0;
    restored->load_state(frames);
    assert(active_count(*restored) == 1);
    const auto& restored_link =
        restored->aa_links[restored->aa_active_link_slots.front()];
    assert(restored_link.arm_state == tether_state);
    assert((actn2_arm_kinetics::arm_i_bound(tether_state)
                ? restored_link.s_i : restored_link.s_j) == bound_site);

    // Put the surviving arm far from the pair's sole canonical association
    // site and tilt the partner slightly.  The old three-site/full-bridge
    // barrier implementation cannot close this tether: its canonical
    // connector is outside capture, and the whole-filament alignment energy
    // is enormous.  A continuous local contour search with a local arm
    // deformation barrier can still reach the nearby physical contour.
    auto& displaced_tether =
        restored->aa_links[restored->aa_active_link_slots.front()];
    if (actn2_arm_kinetics::arm_i_bound(tether_state)) {
        displaced_tether.s_i = -0.49;
    } else {
        displaced_tether.s_j = -0.49;
    }
    const double displaced_bound_site = -0.49;
    restored->actin.direction[1] = {-1.0, 0.05, 0.0};
    restored->actin.direction[1].normalize();
    restored->actin.update_endpoints();
    restored->kappa_aa = 3000.0;
    restored->aa_tether_rebind_rate = 1e12;
    restored->update_system();
    assert(active_count(*restored) == 1);
    const auto& rebound =
        restored->aa_links[restored->aa_active_link_slots.front()];
    assert(rebound.arm_state == actn2_arm_kinetics::State::weak_bridge);
    assert(rebound.status == 1);
    assert((actn2_arm_kinetics::arm_i_bound(tether_state)
                ? rebound.s_i : rebound.s_j) == displaced_bound_site);
    const double free_site = actn2_arm_kinetics::arm_i_bound(tether_state)
        ? rebound.s_j : rebound.s_i;
    assert(std::abs(free_site) > 0.30);

    // A finite molecular pool also constrains deterministic prebinding; a
    // tether or bridge consumes one molecule until complete detachment.
    restored->aa_prebundle_links_per_pair = 3;
    restored->aa_total_link_pool = 1;
    restored->prebind_local_aa_pair(0, 1, 3);
    assert(active_count(*restored) == 1);

    gsl_rng_free(restored_rng);
    gsl_rng_free(rng);
    std::remove(output.c_str());
    return 0;
}
