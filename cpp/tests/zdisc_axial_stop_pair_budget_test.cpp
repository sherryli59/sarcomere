#include "sarcomere.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>

// Force-pass selection for the axial-stop connector force: every distinct
// actin pair contributes, redundant links on one pair share the budget of
// the link with the greatest energy.

using utils::vec;

namespace {
const std::vector<double> box{6.0, 2.0, 2.0};
const std::array<bool, 3> no_pbc{false, false, false};

std::unique_ptr<Sarcomere> make_model(gsl_rng* rng, int& seed, int links_per_pair) {
    int na = 4, nm = 1, fixed = 0;
    std::vector<double> domain = box;
    double la = 1.0, lm = 1.5, rm = 0.05, am_cut = 0.05;
    double am_opt = 0.08, aa_cut = 0.08, aa_opt = 0.04;
    double zero = 0.0, life = 1e9, ratio = 1.0, dt = 1e-5;
    bool directional = true;
    std::string filename = "/tmp/zdisc_axial_stop_pair_budget_test.h5";
    auto model = std::make_unique<Sarcomere>(
        na, nm, domain, la, lm, rm, am_cut, am_opt, aa_cut, aa_opt,
        zero, zero, life, zero, zero, ratio,
        zero, zero, zero, zero, zero, zero, filename, rng, seed, fixed,
        dt, zero, zero, zero, directional, 24, zero, zero, zero, zero);
    model->set_periodicity(no_pbc);
    model->aa_local_multilink = true;
    model->aa_prebundle_links_per_pair = links_per_pair;
    model->aa_neighbor_identity_enabled = false;
    model->aa_barbed_zone = 0.0;
    model->aa_antiparallel_overlap = 1.0;
    model->aa_link_rupture_length = 0.0;
    model->aa_overlap_occlusion_veto = false;
    model->aa_myosin_contact_epsilon = 0.0;
    model->aa_myosin_contact_shell_energy = 0.0;
    model->aa_zdisc_tip_contact = true;
    model->sterics_enabled = false;
    // Pair A on the +y side, pair B on the -y side, both bridging radially
    // near the myosin + tip (x = 0.75).  Antiparallel actins; the inner
    // (80 nm) actins' left (barbed) ends lie beyond the tip, outside the
    // myosin capsule.
    const double ys[4]{0.08, 0.12, -0.08, -0.12};
    for (int a = 0; a < 4; ++a) {
        const bool inner = a % 2 == 0;
        model->actin.direction[a] = {inner ? -1.0 : 1.0, 0, 0};
        model->actin.center[a] = {inner ? 0.83 : 0.63, ys[a], 0};
    }
    model->actin.update_endpoints();
    model->myosin.center[0] = {0, 0, 0};
    model->myosin.direction[0] = {1, 0, 0};
    model->myosin.update_endpoints();
    return model;
}

// Put a link's two sites at axial position x on its two actins.
void set_sites(Sarcomere& model, size_t slot, double x) {
    auto& link = model.aa_links[slot];
    const double ci = model.actin.center[link.actin_i].x;
    const double cj = model.actin.center[link.actin_j].x;
    const double di = model.actin.direction[link.actin_i].x;
    const double dj = model.actin.direction[link.actin_j].x;
    link.s_i = (x - ci) / di;
    link.s_j = (x - cj) / dj;
}

struct LinkContact { double energy = 0.0, plus_force = 0.0; };

LinkContact kernel(const Sarcomere& model, size_t slot) {
    const auto& link = model.aa_links[slot];
    LinkContact out;
    for (int end_sign : {-1, 1}) {
        const auto c = compute_zdisc_axial_stop_contact(
            model.actin, model.myosin, link.actin_i, link.actin_j,
            link.s_i, link.s_j, 0, end_sign, model.box,
            model.aa_zdisc_tip_onset, model.aa_zdisc_tip_stop_stiffness,
            model.aa_zdisc_tip_stop_force_max,
            model.aa_zdisc_tip_stop_gate_radius,
            model.aa_zdisc_tip_stop_gate_cutoff,
            model.aa_zdisc_tip_stop_epsilon);
        if (!c.applied) continue;
        out.energy += c.energy;
        if (end_sign > 0) out.plus_force += c.axial_force;
    }
    return out;
}

double plus_force(Sarcomere& model) {
    return model.myosin["aa_zdisc_tip_force_plus"][0];
}
}  // namespace

int main() {
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    int seed = 7;
    {
        // One link on pair A.
        auto model = make_model(rng, seed, 1);
        model->prebind_local_aa_pair(0, 1, 1);
        assert(model->aa_active_link_slots.size() == 1);
        const size_t slot = model->aa_active_link_slots.front();
        set_sites(*model, slot, 0.73);
        model->update_system();
        const auto expected = kernel(*model, slot);
        assert(expected.plus_force > 29.0 && expected.plus_force <= 30.0);
        assert(std::abs(plus_force(*model) - expected.plus_force) < 1e-9);
        // Equal axial shares on the two actins of this radial connector.
        assert(std::abs(model->actin["aa_zdisc_tip_force"][0] -
                        0.5 * expected.plus_force) < 1e-9);
        assert(std::abs(model->actin["aa_zdisc_tip_force"][1] -
                        0.5 * expected.plus_force) < 1e-9);
    }
    for (bool deeper_first : {true, false}) {
        // Two redundant links on pair A: only the greater-energy link acts,
        // whichever comes first in canonical order.
        auto model = make_model(rng, seed, 2);
        model->prebind_local_aa_pair(0, 1, 2);
        assert(model->aa_active_link_slots.size() == 2);
        const size_t a = model->aa_active_link_slots[0];
        const size_t b = model->aa_active_link_slots[1];
        set_sites(*model, a, deeper_first ? 0.73 : 0.79);
        set_sites(*model, b, deeper_first ? 0.79 : 0.73);
        model->update_system();
        const auto first = kernel(*model, a), second = kernel(*model, b);
        assert(first.plus_force > 0.0 && second.plus_force > 0.0);
        const auto& winner = first.energy > second.energy ? first : second;
        const auto& loser = first.energy > second.energy ? second : first;
        assert(winner.plus_force > loser.plus_force + 5.0);
        assert(std::abs(plus_force(*model) - winner.plus_force) < 1e-9);
    }
    {
        // Two distinct pairs: both contribute.
        auto model = make_model(rng, seed, 1);
        model->prebind_local_aa_pair(0, 1, 1);
        model->prebind_local_aa_pair(2, 3, 1);
        assert(model->aa_active_link_slots.size() == 2);
        for (size_t slot : model->aa_active_link_slots) set_sites(*model, slot, 0.73);
        model->update_system();
        double expected = 0.0;
        for (size_t slot : model->aa_active_link_slots) {
            expected += kernel(*model, slot).plus_force;
        }
        assert(expected > 58.0);
        assert(std::abs(plus_force(*model) - expected) < 1e-9);
        // Deepest contacting point is recorded, 20 nm past the physical tip.
        assert(std::isfinite(model->myosin["aa_zdisc_tip_gap_plus"][0]));
        assert(model->myosin["aa_zdisc_tip_gap_plus"][0] < 0.0);
    }
    gsl_rng_free(rng);
    std::cout << "zdisc axial-stop pair budget: single, redundant (both orders) and distinct pairs passed\n";
    return 0;
}
