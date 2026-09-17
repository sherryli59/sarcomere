#include "am_motor_kinetics.h"
#include "catch_bond_kinetics.h"
#include "langevin.h"
#include "sarcomere.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>

namespace {

std::unique_ptr<Sarcomere> make_model(
    const std::string& output, gsl_rng* rng, int& seed,
    bool prebind = true, double association_rate = 0.0) {
    int n_actins = 2;
    int n_myosins = 1;
    std::vector<double> box{4.0, 4.0, 4.0};
    double actin_length = 1.0;
    double myosin_length = 1.5;
    double myosin_radius = 0.025;
    double am_cutoff = 0.05;
    double am_optimal = 0.055;
    double aa_cutoff = 0.08;
    double aa_optimal = 0.05;
    double k_on = association_rate;
    double k_off = 0.0;
    double base_lifetime = 1e9;
    double directional_base_lifetime = 0.0;
    double lifetime_coeff = 0.0;
    double diff_ratio = 1.0;
    double zero = 0.0;
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
    model->aa_local_multilink = true;
    model->aa_prebundle_links_per_pair = 1;
    model->aa_neighbor_identity_enabled = false;
    model->aa_status2_locked_until_detach = true;
    model->aa_link_rupture_length = 0.0;
    model->sterics_enabled = false;

    // Antiparallel actins at the ACTN2 radial rest separation.
    model->actin.center[0] = {0.0, 0.0, 0.0};
    model->actin.direction[0] = {1.0, 0.0, 0.0};
    model->actin.center[1] = {0.0, aa_optimal, 0.0};
    model->actin.direction[1] = {-1.0, 0.0, 0.0};
    model->actin.update_endpoints();

    model->myosin.direction[0] = {1.0, 0.0, 0.0};
    model->myosin.center[0] = {1.5, 0.0, 0.0};
    model->myosin.update_endpoints();

    if (prebind) {
        model->prebind_local_aa_pair(0, 1, 1);
        assert(model->aa_active_link_slots.size() == 1);
        model->aa_links[model->aa_active_link_slots.front()].status =
            static_cast<Sarcomere::BondFlag>(2);
    }
    return model;
}

size_t active_link_count(const Sarcomere& model) {
    size_t count = 0;
    for (const auto& link : model.aa_links) {
        if (link.active) ++count;
    }
    return count;
}

}  // namespace

int main() {
    const std::string output = "/tmp/actn2_am_occlusion_test.h5";
    std::remove(output.c_str());
    int seed = 8123;
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, seed);
    auto model = make_model(output, rng, seed);

    // With the myosin outside AM capture, the pre-existing mature link remains.
    model->update_system();
    assert(active_link_count(*model) == 1);
    assert(model->actin_crosslink_ratio[0] > 0.0);
    assert(model->actin_crosslink_ratio[1] > 0.0);

    // Reversible maturation changes conformation without removing the ACTN2
    // molecule.  An effectively certain unloaded transition returns status 2
    // to the weak state while retaining the same link and material sites.
    const std::string unlock_output =
        "/tmp/actn2_status2_unlock_test.h5";
    std::remove(unlock_output.c_str());
    auto unlock_model = make_model(unlock_output, rng, seed);
    unlock_model->aa_status2_locked_until_detach = false;
    unlock_model->aa_status2_unlock_rate = 1e12;
    unlock_model->aa_status2_unlock_force_scale = 0.1;
    unlock_model->aa_status2_unlock_registry_scale = 0.02;
    const auto unlock_slot = unlock_model->aa_active_link_slots.front();
    const auto unlock_link_id = unlock_model->aa_links[unlock_slot].id;
    const double old_s_i = unlock_model->aa_links[unlock_slot].s_i;
    const double old_s_j = unlock_model->aa_links[unlock_slot].s_j;
    unlock_model->update_system();
    assert(active_link_count(*unlock_model) == 1);
    assert(unlock_model->aa_links[unlock_slot].id == unlock_link_id);
    assert(unlock_model->aa_links[unlock_slot].status == 1);
    assert(unlock_model->aa_links[unlock_slot].s_i == old_s_i);
    assert(unlock_model->aa_links[unlock_slot].s_j == old_s_j);

    // Move the thick filament over both complete actins. The current AM
    // geometry leaves no barbed-end contour on either actin, so the existing
    // status-2 link must rupture in this update despite status-2 locking and an
    // effectively infinite catch-bond lifetime.
    model->myosin.center[0] = {0.0, 0.0, 0.0};
    model->myosin.update_endpoints();
    model->update_system();
    assert(model->actin_crosslink_ratio[0] <= 1e-12);
    assert(model->actin_crosslink_ratio[1] <= 1e-12);
    assert(active_link_count(*model) == 0);

    // The v3 audit row records a deterministic geometric rupture with the
    // dedicated AM-overlap reason code and no stochastic draw.
    constexpr size_t width = 35;
    assert(model->aa_local_link_events_v3.size() >= width);
    const size_t start = model->aa_local_link_events_v3.size() - width;
    assert(static_cast<int>(model->aa_local_link_events_v3[start + 4]) == 2);
    assert(static_cast<int>(model->aa_local_link_events_v3[start + 30]) ==
           static_cast<int>(
               catch_bond_kinetics::DetachmentReason::am_overlap_occlusion));
    assert(std::isnan(model->aa_local_link_events_v3[start + 32]));
    assert(model->aa_local_link_events_v3[start + 33] == 1.0);
    assert(model->aa_local_link_events_v3[start + 34] == 0.0);

    // Diffusion-mode ACTN2 uses the same explicit connector/thick-filament
    // contact.  A perpendicular thick filament can cross the connector while
    // leaving exposed actin contour on both filaments; this isolates the deep
    // geometric-clash rule from the whole-contour AM-occlusion rule above.
    const std::string collision_output =
        "/tmp/actn2_am_diffusion_collision_test.h5";
    std::remove(collision_output.c_str());
    auto collision_model = make_model(collision_output, rng, seed);
    collision_model->aa_arm_state_model = false;
    collision_model->aa_myosin_contact_shell_width = 0.015;
    collision_model->aa_myosin_contact_shell_energy = 0.0214;
    collision_model->aa_occlusion_rupture_penetration = 0.010;
    collision_model->myosin.center[0] = {0.0, 0.025, 0.0};
    collision_model->myosin.direction[0] = {0.0, 1.0, 0.0};
    collision_model->myosin.update_endpoints();
    collision_model->update_system();
    assert(collision_model->actin_crosslink_ratio[0] > 0.0);
    assert(collision_model->actin_crosslink_ratio[1] > 0.0);
    assert(active_link_count(*collision_model) == 0);
    assert(collision_model->aa_local_link_events_v3.size() >= width);
    const size_t collision_start =
        collision_model->aa_local_link_events_v3.size() - width;
    assert(static_cast<int>(collision_model->aa_local_link_events_v3[
               collision_start + 30]) ==
           static_cast<int>(
               catch_bond_kinetics::DetachmentReason::actn2_myosin_collision));

    // A head-free, spatially remote thick filament does not suppress a valid
    // association, and the accepted-birth audit records no occupied tip zone.
    const std::string accessible_birth_output =
        "/tmp/actn2_tip_accessible_birth_test.h5";
    std::remove(accessible_birth_output.c_str());
    auto accessible_birth_model = make_model(
        accessible_birth_output, rng, seed, false, 1e12);
    accessible_birth_model->aa_zdisc_tip_contact = true;
    accessible_birth_model->update_system();
    assert(active_link_count(*accessible_birth_model) == 1);
    assert(accessible_birth_model->aa_tip_accessibility_events.size() == 13);
    assert(accessible_birth_model->aa_tip_accessibility_events[11] == 0.0);

    // A crossed weak bridge is not deleted while exposed contour remains. It
    // retains an annealing route, and the axial end-stop carries the load.
    const std::string tip_output =
        "/tmp/actn2_zdisc_tip_occlusion_test.h5";
    std::remove(tip_output.c_str());
    auto tip_model = make_model(tip_output, rng, seed);
    tip_model->aa_zdisc_tip_contact = true;
    // Radial connector/thick-filament clash rupture is tested independently
    // above. Disable it here so this block isolates the axial accessibility
    // rule for a nearby myosin tip.
    tip_model->aa_myosin_contact_epsilon = 0.0;
    tip_model->aa_occlusion_rupture_penetration = 0.0;
    tip_model->aa_links[tip_model->aa_active_link_slots.front()].status = 1;
    tip_model->myosin.center[0] = {0.0, -0.715, 0.0};
    tip_model->myosin.direction[0] = {0.0, 1.0, 0.0};
    tip_model->myosin.update_endpoints();
    tip_model->set_am_head_count(1);
    auto& tip_head = tip_model->am_heads_by_myosin[0][1][0];
    tip_head.bound_actin = 0;
    tip_head.state = static_cast<Sarcomere::BondFlag>(
        am_motor_kinetics::State::weak);
    tip_model->update_system();
    assert(tip_model->actin_crosslink_ratio[0] > 0.0);
    assert(tip_model->actin_crosslink_ratio[1] > 0.0);
    assert(active_link_count(*tip_model) == 1);
    const auto& retained_tip_link =
        tip_model->aa_links[tip_model->aa_active_link_slots.front()];
    const auto engaged_tip_exclusion =
        tip_model->evaluate_local_aa_myosin_exclusion(
            retained_tip_link.actin_i, retained_tip_link.actin_j,
            retained_tip_link.s_i, retained_tip_link.s_j);
    assert(engaged_tip_exclusion.axial_tip_contact);
    assert(engaged_tip_exclusion.axial_tip_penetration > 0.0);

    // The same occupied geometry cannot nucleate a new weak bridge. The
    // enormous on-rate makes failure to bind deterministic evidence for the
    // local accessibility gate rather than a missed stochastic event.
    const std::string birth_output =
        "/tmp/actn2_zdisc_tip_birth_gate_test.h5";
    std::remove(birth_output.c_str());
    auto birth_model = make_model(
        birth_output, rng, seed, false, 1e12);
    birth_model->aa_zdisc_tip_contact = true;
    birth_model->aa_myosin_contact_epsilon = 0.0;
    birth_model->aa_occlusion_rupture_penetration = 0.0;
    // The tip is 10 nm short of the patch: inside the 50-nm end-stop zone,
    // but not yet geometrically overlapping the attachment site.
    birth_model->myosin.center[0] = {0.0, -0.735, 0.0};
    birth_model->myosin.direction[0] = {0.0, 1.0, 0.0};
    birth_model->myosin.update_endpoints();
    birth_model->set_am_head_count(1);
    auto& birth_head = birth_model->am_heads_by_myosin[0][1][0];
    birth_head.bound_actin = 0;
    birth_head.state = static_cast<Sarcomere::BondFlag>(
        am_motor_kinetics::State::weak);
    birth_model->update_system();
    assert(birth_model->actin_crosslink_ratio[0] > 0.0);
    assert(birth_model->actin_crosslink_ratio[1] > 0.0);
    assert(active_link_count(*birth_model) == 0);

    // Motor occupancy is irrelevant to steric accessibility: the identical
    // candidate geometry remains excluded after removing the bound head.
    const auto head_free_candidate =
        birth_model->evaluate_local_aa_myosin_exclusion(0, 1, 0.0, 0.0);
    birth_head.bound_actin = -1;
    birth_head.state = static_cast<Sarcomere::BondFlag>(
        am_motor_kinetics::State::detached);
    const auto head_free_exclusion =
        birth_model->evaluate_local_aa_myosin_exclusion(0, 1, 0.0, 0.0);
    assert(head_free_candidate.axial_tip_contact);
    assert(head_free_exclusion.axial_tip_contact);

    // Diffusion may retreat from an existing tip corona, but cannot deepen
    // that contact or exchange it for a different physical tip.
    Sarcomere::AAMyosinExclusionEvaluation no_tip;
    Sarcomere::AAMyosinExclusionEvaluation tip_one_shallow;
    tip_one_shallow.axial_tip_contacts.push_back(
        {0, 1, 0.040, 0.080, 0.010, false});
    Sarcomere::AAMyosinExclusionEvaluation tip_one_deep;
    tip_one_deep.axial_tip_contacts.push_back(
        {0, 1, 0.030, 0.080, 0.020, false});
    Sarcomere::AAMyosinExclusionEvaluation tip_two_shallow;
    tip_two_shallow.axial_tip_contacts.push_back(
        {1, -1, 0.040, 0.080, 0.010, false});
    assert(!Sarcomere::aa_tip_diffusion_move_allowed(
        no_tip, tip_one_shallow));
    assert(Sarcomere::aa_tip_diffusion_move_allowed(
        tip_one_deep, tip_one_shallow));
    assert(!Sarcomere::aa_tip_diffusion_move_allowed(
        tip_one_shallow, tip_one_deep));
    assert(!Sarcomere::aa_tip_diffusion_move_allowed(
        tip_one_shallow, tip_two_shallow));

    // Thick-filament occupancy is head independent. The same crossed geometry
    // is an axial exclusion contact even with no motor head on either actin;
    // the pre-existing molecule remains until a physical detachment event.
    const std::string unrelated_output =
        "/tmp/actn2_zdisc_tip_unrelated_test.h5";
    std::remove(unrelated_output.c_str());
    auto unrelated_model = make_model(unrelated_output, rng, seed);
    unrelated_model->aa_zdisc_tip_contact = true;
    unrelated_model->aa_myosin_contact_epsilon = 0.0;
    unrelated_model->aa_occlusion_rupture_penetration = 0.0;
    unrelated_model->aa_links[
        unrelated_model->aa_active_link_slots.front()].status = 1;
    unrelated_model->myosin.center[0] = {0.0, -0.715, 0.0};
    unrelated_model->myosin.direction[0] = {0.0, 1.0, 0.0};
    unrelated_model->myosin.update_endpoints();
    unrelated_model->update_system();
    assert(active_link_count(*unrelated_model) == 1);
    const auto& unrelated_link =
        unrelated_model->aa_links[
            unrelated_model->aa_active_link_slots.front()];
    const auto unrelated_exclusion =
        unrelated_model->evaluate_local_aa_myosin_exclusion(
            unrelated_link.actin_i, unrelated_link.actin_j,
            unrelated_link.s_i, unrelated_link.s_j);
    assert(unrelated_exclusion.axial_tip_contact);
    assert(unrelated_exclusion.axial_tip_penetration > 0.0);
    assert(!unrelated_exclusion.nearest_tip_motor_engaged);

    // Position-level overdamped contact resolves the invalid linked geometry
    // without deleting the ACTN2 molecule. The committed gap is nonnegative.
    double beta = 1.0;
    double dt = 1e-5;
    double d_actin = 1.0;
    double d_actin_rot = 0.0;
    double d_myosin = 1.0;
    double d_myosin_rot = 0.0;
    int save_every = 100;
    bool resume = false;
    int fixed_myosins = 0;
    Langevin projection_integrator(
        *unrelated_model, beta, dt, d_actin, d_actin_rot,
        d_myosin, d_myosin_rot, save_every, resume, true,
        0.0, 0.0, 0.0, 0.0);
    projection_integrator.resolve_actn2_tip_hard_contacts(fixed_myosins);
    assert(active_link_count(*unrelated_model) == 1);
    const auto projected_exclusion =
        unrelated_model->evaluate_local_aa_myosin_exclusion(
            unrelated_link.actin_i, unrelated_link.actin_j,
            unrelated_link.s_i, unrelated_link.s_j);
    assert(projected_exclusion.axial_tip_contact);
    for (const auto& contact : projected_exclusion.axial_tip_contacts) {
        assert(contact.axial_gap >= -1e-9);
    }

    gsl_rng_free(rng);
    return 0;
}
