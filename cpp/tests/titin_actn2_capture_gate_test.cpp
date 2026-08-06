// Deterministic semantics of the ACTN2-support capture gate
// (--titin_require_live_actn2).
//
// The tension-only titin force law is exactly zero over the whole slack
// region, so elastic gating cannot discriminate between candidates there and
// every barbed end inside the 0.85 um capture shell receives the full k_on.
// In a dense box that filled the monovalent ownership cap to 100% within the
// first analyzed millisecond.  The gate makes a *new* association additionally
// require local evidence of a nascent Z-disc complex: at least one live ACTN2
// A-A link on that individual barbed end.
//
// The six properties asserted here are exactly the ones that distinguish this
// rule from the two neighbouring rules that would be wrong -- a status-2
// requirement (which would couple titin to mature motor load) and a
// continuously enforced requirement (which would make the complex brittle
// under ordinary ACTN2 turnover).
//
// The model is driven with association certain and dissociation impossible
// (k_on huge, k_off zero) so that every assertion is a statement about the
// gate rather than about a random draw.

#include "sarcomere.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>

namespace {

// Two antiparallel actins whose barbed ends sit inside the capture shell of a
// single myosin's positive half.  aa_local_multilink is left off so that the
// ACTN2 link state is set explicitly by the test and is never overwritten by
// the stochastic A-A chemistry: the gate's predicate is what is under test,
// not the kinetics that populate it.
std::unique_ptr<Sarcomere> make_model(
    const std::string& output, gsl_rng* rng, int& seed, bool gate,
    double k_off_unsupported = 0.0, bool rebind_requires_actn2 = false) {
    int n_actins = 2;
    int n_myosins = 1;
    std::vector<double> box{4.0, 4.0, 4.0};
    double actin_length = 1.0;
    double myosin_length = 1.5;
    double myosin_radius = 0.025;
    double am_cutoff = 0.10;
    double am_optimal = 0.055;
    double aa_cutoff = 0.08;
    double aa_optimal = 0.0506;
    double zero = 0.0;
    double diff_ratio = 1.0;
    int fixed_myosins = 0;
    double dt = 1e-4;
    bool directional = true;
    int max_myosin_bonds = 6;
    std::string filename = output;
    auto model = std::make_unique<Sarcomere>(
        n_actins, n_myosins, box, actin_length, myosin_length,
        myosin_radius, am_cutoff, am_optimal, aa_cutoff, aa_optimal,
        zero, zero, zero, zero, zero, diff_ratio, zero, zero, zero,
        zero, zero, zero, filename, rng, seed, fixed_myosins, dt, zero,
        40.0, 0.50, directional, max_myosin_bonds, zero, zero, zero,
        zero);
    model->set_periodicity({false, false, false});
    model->dynamic_am = false;
    model->aa_local_multilink = false;
    model->sterics_enabled = false;
    model->titin_mode = Sarcomere::TitinMode::explicit_dynamic;
    // k_off = 0 so a bond, once formed, can only be removed by the code paths
    // under test.  rebind memory 0.02 s with a x100 multiplier matches the
    // promoted production parameters.
    model->configure_explicit_titin(3, 0.70, 1e12, 0.0, false, 0.5,
                                    0.02, 100.0, gate, k_off_unsupported,
                                    rebind_requires_actn2);
    model->aa_links_by_actin.assign(static_cast<size_t>(n_actins), {});

    model->myosin.center[0] = {0.0, 0.0, 0.0};
    model->myosin.direction[0] = {1.0, 0.0, 0.0};
    model->reset_myosin_material_frames();
    model->myosin.update_endpoints();

    // Both barbed ends 0.60 um behind track zero on the positive half, so both
    // are geometrically eligible and only the gate can separate them.  Actin 1
    // is displaced transversely by well under the capture distance.
    const auto positive_track = static_cast<size_t>(3);
    const vec anchor = model->explicit_titin_thick_anchor(
        model->explicit_titins[positive_track]);
    for (int i = 0; i < n_actins; ++i) {
        model->actin.direction[i] = {1.0, 0.0, 0.0};
        model->actin.center[i] = anchor - vec{0.60, 0.0, 0.0} +
                                 vec{0.50, 0.0, 0.0} +
                                 vec{0.0, 0.02 * i, 0.0};
    }
    model->actin.update_endpoints();
    return model;
}

// Install one live ACTN2 A-A link between the two actins at the requested
// status, mirroring what _rebuild_local_aa_indices() would produce.
void add_actn2_link(Sarcomere& model, int status) {
    Sarcomere::AACrosslink link;
    link.id = 1;
    link.actin_i = 0;
    link.actin_j = 1;
    link.s_i = -0.45;
    link.s_j = -0.45;
    link.status = static_cast<Sarcomere::BondFlag>(status);
    link.active = true;
    const size_t slot = model.aa_links.size();
    model.aa_links.push_back(link);
    model.aa_active_link_slots.push_back(slot);
    model.aa_links_by_actin[0].push_back(slot);
    model.aa_links_by_actin[1].push_back(slot);
}

void clear_actn2_links(Sarcomere& model) {
    for (auto& link : model.aa_links) link.active = false;
    model.aa_active_link_slots.clear();
    for (auto& bucket : model.aa_links_by_actin) bucket.clear();
}

int titin_owner_of(const Sarcomere& model, int actin_index) {
    return model.explicit_titin_owner_by_actin[
        static_cast<size_t>(actin_index)];
}

size_t bound_titin_count(const Sarcomere& model) {
    size_t count = 0;
    for (const auto& titin : model.explicit_titins) {
        if (titin.bound) ++count;
    }
    return count;
}

}  // namespace

int main() {
    const std::string output = "/tmp/titin_actn2_capture_gate_test.h5";
    int seed = 4127;

    // 1. Ungated control. Both isolated actins are captured immediately, which
    //    is the saturation mechanism the gate exists to remove. Without this
    //    baseline a later "nothing bound" assertion would not distinguish the
    //    gate from a broken geometry.
    {
        std::remove(output.c_str());
        gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng, seed);
        auto model = make_model(output, rng, seed, false);
        model->update_system();
        assert(bound_titin_count(*model) == 2);
        assert(titin_owner_of(*model, 0) >= 0);
        assert(titin_owner_of(*model, 1) >= 0);
        assert(model->explicit_titin_supported_associations == 0);
        assert(model->explicit_titin_unsupported_associations == 2);
        assert(model->explicit_titin_rebind_associations == 0);
        gsl_rng_free(rng);
    }

    // 2. Gated, no ACTN2. An isolated actin inside the capture shell cannot
    //    bind, and no unsupported association is ever recorded.
    {
        std::remove(output.c_str());
        gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng, seed);
        auto model = make_model(output, rng, seed, true);
        for (int step = 0; step < 5; ++step) model->update_system();
        assert(bound_titin_count(*model) == 0);
        assert(titin_owner_of(*model, 0) < 0);
        assert(titin_owner_of(*model, 1) < 0);
        assert(model->explicit_titin_supported_associations == 0);
        assert(model->explicit_titin_unsupported_associations == 0);
        gsl_rng_free(rng);
    }

    // 3. A live status-1 ACTN2 link is sufficient: the same geometry that was
    //    rejected above now binds, and the association is scored as supported.
    //    Status 2 must NOT be required -- it is gated on mature motor load, so
    //    requiring it would reproduce the original low-capture failure.
    {
        std::remove(output.c_str());
        gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng, seed);
        auto model = make_model(output, rng, seed, true);
        model->update_system();
        assert(bound_titin_count(*model) == 0);

        add_actn2_link(*model, 1);
        model->update_system();
        assert(bound_titin_count(*model) == 2);
        assert(titin_owner_of(*model, 0) >= 0);
        assert(titin_owner_of(*model, 1) >= 0);
        assert(model->explicit_titin_supported_associations == 2);
        assert(model->explicit_titin_unsupported_associations == 0);
        gsl_rng_free(rng);
    }

    // 4. Status 2 also works, confirming the predicate is "any live status"
    //    rather than "status 1 only".
    {
        std::remove(output.c_str());
        gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng, seed);
        auto model = make_model(output, rng, seed, true);
        add_actn2_link(*model, 2);
        model->update_system();
        assert(bound_titin_count(*model) == 2);
        assert(model->explicit_titin_supported_associations == 2);
        gsl_rng_free(rng);
    }

    // 5. Losing the ACTN2 link does not release an already bound titin. The
    //    gate is a formation rule; ordinary ACTN2 turnover under an
    //    established titin must not dissolve the complex.
    {
        std::remove(output.c_str());
        gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng, seed);
        auto model = make_model(output, rng, seed, true);
        add_actn2_link(*model, 1);
        model->update_system();
        assert(bound_titin_count(*model) == 2);
        const int owner0 = titin_owner_of(*model, 0);
        const int owner1 = titin_owner_of(*model, 1);
        assert(owner0 >= 0 && owner1 >= 0);

        clear_actn2_links(*model);
        for (int step = 0; step < 20; ++step) model->update_system();
        assert(bound_titin_count(*model) == 2);
        assert(titin_owner_of(*model, 0) == owner0);
        assert(titin_owner_of(*model, 1) == owner1);
        // No new associations were counted while support was absent.
        assert(model->explicit_titin_supported_associations == 2);
        assert(model->explicit_titin_unsupported_associations == 0);
        gsl_rng_free(rng);
    }

    // 6. Same-partner recapture inside the memory window is exempt from the
    //    gate, and once the window has expired a new association again
    //    requires live ACTN2. The two halves are asserted against the same
    //    forced detachment so the only difference is elapsed time.
    //
    //    dt = 1e-4 s and the memory is 0.02 s, i.e. 200 steps.
    {
        std::remove(output.c_str());
        gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng, seed);
        auto model = make_model(output, rng, seed, true);
        add_actn2_link(*model, 1);
        model->update_system();
        assert(bound_titin_count(*model) == 2);
        const int owner0 = titin_owner_of(*model, 0);

        // Force a molecular detachment of actin 0's titin and remove all ACTN2
        // support, then step once while still well inside the memory window.
        clear_actn2_links(*model);
        auto& titin = model->explicit_titins[static_cast<size_t>(owner0)];
        titin.bound = false;
        titin.bound_actin = -1;
        titin.candidate_actin = -1;
        titin.last_bound_actin = 0;
        titin.last_detach_step =
            static_cast<decltype(titin.last_detach_step)>(model->current_step);
        model->explicit_titin_owner_by_actin[0] = -1;
        const auto supported_before =
            model->explicit_titin_supported_associations;

        model->update_system();
        assert(titin_owner_of(*model, 0) == owner0);
        assert(model->explicit_titins[static_cast<size_t>(owner0)]
                   .bound_actin == 0);
        assert(model->explicit_titin_rebind_associations == 1);
        // A recapture is never scored as supported or unsupported: it bypasses
        // the gate by design and must not read as gate leakage.
        assert(model->explicit_titin_supported_associations ==
               supported_before);
        assert(model->explicit_titin_unsupported_associations == 0);

        // Now detach again and let the memory window expire. With no live
        // ACTN2 the identical geometry must stay unbound.
        auto& again = model->explicit_titins[static_cast<size_t>(owner0)];
        again.bound = false;
        again.bound_actin = -1;
        again.candidate_actin = -1;
        again.last_bound_actin = 0;
        again.last_detach_step =
            static_cast<decltype(again.last_detach_step)>(model->current_step);
        model->explicit_titin_owner_by_actin[0] = -1;
        const size_t detach_step = model->current_step;

        for (int step = 0; step < 260; ++step) {
            model->update_system();
            const double elapsed = model->dt * static_cast<double>(
                model->current_step - detach_step);
            if (elapsed > model->explicit_titin_rebind_memory) break;
        }
        // Everything up to the expiry was a recapture; drive it unbound one
        // last time and confirm the post-expiry search finds nothing.
        auto& expired = model->explicit_titins[static_cast<size_t>(owner0)];
        expired.bound = false;
        expired.bound_actin = -1;
        expired.candidate_actin = -1;
        model->explicit_titin_owner_by_actin[0] = -1;
        const auto rebinds_before =
            model->explicit_titin_rebind_associations;

        for (int step = 0; step < 10; ++step) model->update_system();
        assert(titin_owner_of(*model, 0) < 0);
        assert(model->explicit_titin_rebind_associations == rebinds_before);
        assert(model->explicit_titin_unsupported_associations == 0);

        // Restoring ACTN2 support reopens the ordinary path.
        add_actn2_link(*model, 1);
        model->update_system();
        assert(titin_owner_of(*model, 0) >= 0);
        assert(model->explicit_titin_supported_associations ==
               supported_before + 1);
        gsl_rng_free(rng);
    }

    // 7. Serial-anchor dissociation. Titin reaches actin through ACTN2, so a
    //    bound titin whose barbed end has lost every live ACTN2 link runs at
    //    the unsupported rate instead of its own k_off. This is what gives
    //    occupancy a fixed point: without it, a permanent latch fed by
    //    transient encounters ratchets to full occupancy.
    //
    //    Rebinding memory is disabled here (k_off_unsupported is tested
    //    against a fresh model built with rebind memory 0.02 s, so the run is
    //    kept shorter than the window would matter for a *different* partner,
    //    and the same-partner recapture is asserted separately below).
    {
        std::remove(output.c_str());
        gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng, seed);
        // Certain release when unsupported; k_off stays 0 when supported.
        auto model = make_model(output, rng, seed, true, 1e12);
        add_actn2_link(*model, 1);
        model->update_system();
        assert(bound_titin_count(*model) == 2);
        assert(model->explicit_titin_unsupported_releases == 0);

        // Supported titins are untouched by the unsupported rate however many
        // steps elapse: the anchor is intact, so k_off = 0 governs.
        for (int step = 0; step < 10; ++step) model->update_system();
        assert(bound_titin_count(*model) == 2);
        assert(model->explicit_titin_unsupported_releases == 0);

        // Remove the ACTN2 link. The serial anchor is now absent, so both
        // titins release -- stochastically, at the unsupported rate.
        clear_actn2_links(*model);
        model->update_system();
        assert(model->explicit_titin_unsupported_releases == 2);
        assert(titin_owner_of(*model, 0) < 0);
        assert(titin_owner_of(*model, 1) < 0);
        gsl_rng_free(rng);
    }

    // 8. The serial anchor must not be brittle. Inside the rebinding memory
    //    window a released titin can still recapture the same actin by
    //    identity, so a momentary ACTN2 gap anneals instead of dissolving the
    //    complex permanently.
    //
    //    k_off_unsupported is 1e12 here, i.e. release is certain on every step
    //    while support is absent. Because binding and release are evaluated in
    //    the same sweep, a recaptured titin is released again before the step
    //    ends, so the observable is that the recapture *happened* -- the
    //    counter -- rather than the ownership surviving the step. At a
    //    physical rate (order 100 /s, p ~ 5e-5 per step) that coincidence is
    //    negligible and ownership simply persists.
    {
        std::remove(output.c_str());
        gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng, seed);
        auto model = make_model(output, rng, seed, true, 1e12);
        add_actn2_link(*model, 1);
        model->update_system();
        const int owner0 = titin_owner_of(*model, 0);
        assert(owner0 >= 0);

        clear_actn2_links(*model);
        model->update_system();
        assert(model->explicit_titin_unsupported_releases >= 1);

        // Still inside the 0.02 s window and with ACTN2 absent, the only path
        // that can fire is same-partner recapture. Its counter proves the
        // memory window bypasses the gate under the serial-anchor law too.
        const auto rebinds_before = model->explicit_titin_rebind_associations;
        model->update_system();
        assert(model->explicit_titin_rebind_associations > rebinds_before);

        // Restoring support ends the gap and the actin is anchored again.
        add_actn2_link(*model, 1);
        model->update_system();
        assert(titin_owner_of(*model, 0) >= 0);
        gsl_rng_free(rng);
    }

    // 9. Zero disables the law exactly: a bound titin with no ACTN2 support
    //    keeps its own k_off, which is the pre-existing single-rate behaviour.
    {
        std::remove(output.c_str());
        gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng, seed);
        auto model = make_model(output, rng, seed, true, 0.0);
        add_actn2_link(*model, 1);
        model->update_system();
        assert(bound_titin_count(*model) == 2);
        clear_actn2_links(*model);
        for (int step = 0; step < 50; ++step) model->update_system();
        assert(bound_titin_count(*model) == 2);
        assert(model->explicit_titin_unsupported_releases == 0);
        gsl_rng_free(rng);
    }

    // 10. Recapture under the serial anchor. With the memory exempt from the
    //     ACTN2 requirement, an unsupported release is undone on the very next
    //     step and the serial law has no net effect -- this was measured at
    //     94% of releases in the 640-actin scan. Requiring support for
    //     recapture too keeps the memory's rate advantage without its
    //     permanence, so the actin stays free until support returns.
    {
        std::remove(output.c_str());
        gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng, seed);
        auto model = make_model(output, rng, seed, true, 1e12, true);
        add_actn2_link(*model, 1);
        model->update_system();
        assert(bound_titin_count(*model) == 2);

        clear_actn2_links(*model);
        model->update_system();
        assert(model->explicit_titin_unsupported_releases >= 1);
        const auto rebinds_after_release =
            model->explicit_titin_rebind_associations;

        // No recapture may occur while support is absent, however long we wait
        // inside the 0.02 s memory window.
        for (int step = 0; step < 30; ++step) model->update_system();
        assert(titin_owner_of(*model, 0) < 0);
        assert(titin_owner_of(*model, 1) < 0);
        assert(model->explicit_titin_rebind_associations ==
               rebinds_after_release);

        // Support returns inside the window: recapture fires immediately, and
        // it is scored as a rebind rather than as a fresh association.
        add_actn2_link(*model, 1);
        model->update_system();
        assert(model->explicit_titin_rebind_associations >
               rebinds_after_release);
        gsl_rng_free(rng);
    }

    std::remove(output.c_str());
    return 0;
}
