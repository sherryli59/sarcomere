#include "actn2_maturation_geometry.h"
#include "sarcomere.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>

using utils::vec;

namespace {
const std::vector<double> box{6.0, 2.0, 2.0};
const std::array<bool, 3> no_pbc{false, false, false};
const std::array<bool, 3> pbc{true, true, true};

bool overlaps(vec p, vec q, vec c = {0, 0, 0}, vec u = {1, 0, 0},
              const std::array<bool, 3>& periodic = no_pbc,
              const std::vector<double>& domain = box) {
    return actn2_maturation::connector_overlaps_myosin(
        p, q, c, u, 1.5, domain, periodic);
}

std::unique_ptr<Sarcomere> make_model(gsl_rng* rng, int& seed) {
    int na = 2, nm = 1, fixed = 0;
    std::vector<double> domain = box;
    double la = 1.0, lm = 1.5, rm = 0.05, am_cut = 0.05;
    double am_opt = 0.08, aa_cut = 0.08, aa_opt = 0.04;
    double zero = 0.0, life = 1e9, ratio = 1.0, dt = 1e-5;
    bool directional = true;
    std::string filename = "/tmp/status2_connector_gate_test.h5";
    auto model = std::make_unique<Sarcomere>(
        na, nm, domain, la, lm, rm, am_cut, am_opt, aa_cut, aa_opt,
        zero, zero, life, zero, zero, ratio,
        zero, zero, zero, zero, zero, zero, filename, rng, seed, fixed,
        dt, zero, zero, zero, directional, 24, zero, zero, zero, zero);
    model->set_periodicity(no_pbc);
    model->aa_local_multilink = true;
    model->aa_prebundle_links_per_pair = 1;
    model->aa_neighbor_identity_enabled = false;
    model->aa_barbed_zone = 0.0;
    // These antiparallel 1 um actins have full axial overlap. Set the
    // registry reference accordingly, so the existing maturation gate passes.
    model->aa_antiparallel_overlap = 1.0;
    model->aa_link_rupture_length = 0.0;
    model->aa_overlap_occlusion_veto = false;
    model->aa_myosin_contact_epsilon = 0.0;
    model->aa_myosin_contact_shell_energy = 0.0;
    model->aa_zdisc_tip_contact = false;
    model->sterics_enabled = false;
    model->aa_status_mode = actn2_load_kinetics::StatusMode::productive_register;
    model->aa_status2_force_on = 0.5;
    model->aa_status2_locked_until_detach = true;
    model->aa_status2_support_window = 1.0;
    model->actin.center[0] = {0, 0.120, 0};
    model->actin.center[1] = {0, 0.160, 0};
    model->actin.direction[0] = {1, 0, 0};
    model->actin.direction[1] = {-1, 0, 0};
    model->actin.update_endpoints();
    model->myosin.center[0] = {0, 0, 0};
    model->myosin.direction[0] = {1, 0, 0};
    model->myosin.update_endpoints();
    model->prebind_local_aa_pair(0, 1, 1);
    assert(model->aa_active_link_slots.size() == 1);
    auto& link = model->aa_links[model->aa_active_link_slots.front()];
    link.status = 1;
    // Represents productive motor support accumulated over the existing
    // maturation window; keeps this test independent of head binding RNG.
    link.support_ema = 2.0;
    return model;
}
}

int main() {
    // r27zt0np_seed20260805.h5, frame 475 / step 237500:
    // link 179352, actins 1316 and 233, myosin 244. The endpoints below
    // are the stored pose in a local basis with myosin axis = x.
    assert(overlaps(
        {-.003812085957410629, .11732520591698475, .03345072485325059},
        {.014175370169763528, .09717262963312846, -.03345072485325024}));
    // The observed loophole: midpoint at 140 nm, connector end at 120 nm.
    assert(overlaps({0, .120, 0}, {0, .160, 0}));
    assert(!overlaps({0, .131, 0}, {0, .160, 0}));
    assert(overlaps({0, .130, 0}, {0, .160, 0}));
    assert(overlaps({0, -.2, 0}, {0, .2, 0})); // interior minimum
    assert(overlaps({0, .1, 0}, {0, .1, 0})); // zero-length connector
    assert(!overlaps({0, .14, 0}, {0, .14, 0}));
    assert(overlaps({-.8, .1, 0}, {.8, .1, 0})); // axial clipping
    assert(overlaps({.70, .1, 0}, {.80, .1, 0}));
    assert(overlaps({.75, .1, 0}, {.75, .12, 0})); // end plane
    assert(!overlaps({.750001, .1, 0}, {.80, .1, 0}));
    assert(!overlaps({-.80, .1, 0}, {-.750001, .1, 0}));
    // Radial proximity outside the slab cannot be combined with axial
    // overlap at a different, radially distant segment point.
    assert(!overlaps({.9, .05, 0}, {.7, .3, 0}));
    assert(overlaps({0, .12, 0}, {0, .16, 0}, {0, 0, 0}, {-1, 0, 0}));
    assert(overlaps({.12, 0, 0}, {.16, 0, 0}, {0, 0, 0}, {0, 1, 0}));
    const double s = std::sqrt(0.5);
    assert(overlaps({-.12*s, .12*s, 0}, {-.16*s, .16*s, 0},
                    {0, 0, 0}, {s, s, 0}));
    // Wrapping a connector across y and translating either filament by a
    // whole box must not change the result.
    assert(overlaps({0, .98, 0}, {0, -.98, 0}, {0, -.9, 0}, {1, 0, 0}, pbc));
    assert(overlaps({6, .98, 0}, {6, 1.02, 0}, {0, -.9, 0}, {1, 0, 0}, pbc));
    assert(!overlaps({0, .98, 0}, {0, 1.02, 0}, {0, -.9, 0}));
    // The nearest cylinder center is not always the intersecting image
    // when a tilted filament is comparable to the box size.
    assert(overlaps({.48, -.48, 0}, {.49, -.48, 0}, {0, 0, 0}, {s, s, 0},
                    pbc, {1.0, 2.0, 2.0}));
    assert(!overlaps({.48, -.48, 0}, {.49, -.48, 0},
                    {0, 0, 0}, {s, s, 0}));

    auto deletion_overlap = [&](vec p, vec q,
                               const std::array<bool, 3>& periodic = no_pbc,
                               vec center = {0, 0, 0}) {
        return actn2_maturation::connector_overlaps_myosin(
            p, q, center, {1, 0, 0}, 1.5, box, periodic,
            actn2_maturation::deletion_radius,
            actn2_maturation::deletion_buffer);
    };
    assert(overlaps({0, .125, 0}, {0, .165, 0}));
    assert(!deletion_overlap({0, .125, 0}, {0, .165, 0}));
    assert(deletion_overlap({0, .119, 0}, {0, .159, 0}));
    assert(deletion_overlap({0, .120, 0}, {0, .160, 0}));
    assert(!deletion_overlap({.745, .1, 0}, {.745, .14, 0}));
    assert(deletion_overlap({.739, .1, 0}, {.739, .14, 0}));
    assert(deletion_overlap({.740, .1, 0}, {.740, .14, 0}));
    assert(!deletion_overlap({.9, .05, 0}, {.7, .3, 0}));
    assert(deletion_overlap({0, .98, 0}, {0, -.98, 0}, pbc, {0, -.9, 0}));

    int seed = 20261028;
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, seed);
    {
        auto model = make_model(rng, seed);
        const size_t slot = model->aa_active_link_slots.front();
        model->update_system();
        assert(model->aa_links[slot].active);
        assert(model->aa_links[slot].status == 1);
        assert(!model->aa_links[slot].candidate_locked);
        assert(model->aa_links[slot].status2_memory_until == 0);
        // A move just outside the physical tip plane now permits maturation,
        // even though it is still inside the soft force's 50 nm onset.
        model->actin.center[0].x = .750001;
        model->actin.center[1].x = .750001;
        model->actin.update_endpoints();
        model->update_system();
        assert(model->aa_links[slot].active);
        assert(model->aa_links[slot].status == 2);
        // Established status 2 may penetrate 5 nm through the end plane.
        model->actin.center[0].x = .745;
        model->actin.center[1].x = .745;
        model->actin.update_endpoints();
        model->update_system();
        assert(model->aa_links[slot].status == 2);
        model->new_file();
        model->save_resume_snapshot();
        {
            H5::H5File output(model->filename, H5F_ACC_RDONLY);
            for (const char* path : {"/meta", "/resume"}) {
                auto group = output.openGroup(path);
                double radius = 0.0;
                group.openAttribute("aa_status2_connector_exclusion_radius").read(
                    H5::PredType::NATIVE_DOUBLE, &radius);
                assert(radius == .130);
                assert(group.attrExists("aa_status2_connector_exclusion_geometry"));
                group.openAttribute("aa_status2_connector_deletion_radius").read(
                    H5::PredType::NATIVE_DOUBLE, &radius);
                assert(std::abs(radius - .120) < 1e-12);
                group.openAttribute("aa_status2_connector_deletion_axial_inset").read(
                    H5::PredType::NATIVE_DOUBLE, &radius);
                assert(radius == .010);
                assert(group.attrExists("aa_status2_connector_deletion_rule"));
                assert(group.attrExists("aa_zdisc_tip_contact_geometry"));
            }
        }
        // Full-connector deletion despite midpoint outside the 120 nm
        // retention cylinder; force and legacy burial veto remain disabled.
        model->actin.center[0] = {0, .119, 0};
        model->actin.center[1] = {0, .159, 0};
        model->actin.update_endpoints();
        model->update_system();
        assert(!model->aa_links[slot].active);
        assert(model->aa_local_link_events_v3.size() >= 35);
        assert(model->aa_local_link_events_v3[
            model->aa_local_link_events_v3.size() - 5] == 5.0);
        assert(model->aa_local_link_events_v3[
            model->aa_local_link_events_v3.size() - 4] == 1.0);
    }
    {
        auto model = make_model(rng, seed);
        const size_t slot = model->aa_active_link_slots.front();
        model->aa_status_mode = actn2_load_kinetics::StatusMode::productive_opposed;
        model->aa_links[slot].support_ema = 0;
        model->aa_links[slot].status2_memory_until = 1000;
        model->am_status2_memory_steps = 1000;
        model->update_system();
        assert(model->aa_links[slot].active);
        assert(model->aa_links[slot].status == 1);
        assert(model->aa_links[slot].status2_memory_until == 0);
    }
    // Inherited mature links obey the rule on their first update, both in
    // the arm-state model and when support drops them back to status 1.
    for (int mode = 0; mode < 3; ++mode) {
        auto model = make_model(rng, seed);
        const size_t slot = model->aa_active_link_slots.front();
        auto& link = model->aa_links[slot];
        link.status = 2;
        if (mode == 1) {
            model->aa_arm_state_model = true;
            link.arm_state = actn2_arm_kinetics::State::mature_bridge;
        } else if (mode == 2) {
            model->aa_status2_locked_until_detach = false;
            link.support_ema = 0.0;
        }
        model->actin.center[0].y = .119;
        model->actin.center[1].y = .159;
        model->actin.update_endpoints();
        model->update_system();
        assert(!model->aa_links[slot].active);
        assert(model->aa_local_link_events_v3[
            model->aa_local_link_events_v3.size() - 5] == 5.0);
    }
    {
        auto model = make_model(rng, seed);
        model->aa_zdisc_tip_contact = true;
        // Deliberately long bridge, with contour rupture disabled by the
        // fixture, isolates coverage beyond the endpoint neighbor spheres.
        model->actin.center[0] = {.5, -2., 0};
        model->actin.center[1] = {.5, 2., 0};
        model->actin.update_endpoints();
        model->update_system();
        // Endpoints remain far from the myosin; their connector crosses it.
        // This must not depend on either actin's A-M neighbor candidates.
        assert(model->actin_neighbors_by_species[0].second.empty());
        assert(model->actin_neighbors_by_species[1].second.empty());
        assert(model->myosin["aa_zdisc_tip_force_plus"][0] > 0);
        assert(model->actin["aa_zdisc_tip_force"][0] > 0);
        assert(model->actin["aa_zdisc_tip_force"][1] > 0);
    }
    gsl_rng_free(rng);
    std::cout << "status2_connector_gate: geometry, maturation, memory and buffered deletion passed\n";
}
