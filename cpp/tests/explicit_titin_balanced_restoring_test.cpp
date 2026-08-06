#include "sarcomere.h"

#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <string>

namespace {

int bound_titins(const Sarcomere& model) {
    int count = 0;
    for (const auto& titin : model.explicit_titins) {
        if (titin.bound) ++count;
    }
    return count;
}

}  // namespace

int main() {
    std::vector<double> box{5.32, 3.0, 3.0};
    std::string filename = "/tmp/explicit_titin_balanced_restoring_test.h5";
    std::remove(filename.c_str());
    int seed = 20260804;
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, seed);

    int n_actins = 4;
    int n_myosins = 2;
    double actin_length = 1.0;
    double myosin_length = 1.5;
    double myosin_radius = 0.025;
    double am_cutoff = 0.10;
    double am_optimal = 0.055;
    double aa_cutoff = 0.10;
    double aa_optimal = 0.055;
    double aa_k_on = 1000.0;
    double diff_ratio = 1.0;
    int fixed_myosins = 0;
    double dt = 5.0e-7;
    bool directional = true;
    double zero = 0.0;
    Sarcomere model(
        n_actins, n_myosins, box, actin_length, myosin_length, myosin_radius,
        am_cutoff, am_optimal, aa_cutoff, aa_optimal,
        aa_k_on, zero, zero, zero, zero,
        diff_ratio, zero, zero, zero, zero, zero, zero,
        filename, rng, seed, fixed_myosins, dt, zero,
        100.0, 0.665, directional, 24, zero, zero, zero, zero);
    model.set_periodicity({true, false, false});
    model.set_am_radial_geometry(0.015, 0.015, zero);
    model.sterics_enabled = false;
    model.dynamic_am = false;
    model.set_local_aa_model(true, 3, 0.0192, 0.1, zero, 0, 3);
    model.aa_link_rupture_length = 0.12;
    model.aa_barbed_zone = 0.20;
    model.set_aa_axial_overlap(0.16, 0.0, 0.1, 0.15, 0.05);
    model.init_aa_overlap = 0.16;
    model.minimal_catch_variant("balanced");
    // Rebuild the neighbour list after the exact geometry replaces the
    // constructor's random coordinates.
    model.set_periodicity({true, false, false});

    model.titin_mode = Sarcomere::TitinMode::explicit_dynamic;
    model.configure_explicit_titin(
        3, 0.80, 0.0, 0.0, false, 0.5,
        0.0, 1.0, true, 0.0, false);

    // Candidate recognition and ACTN2-gated capture have their own exhaustive
    // regression. Here we freeze the four expected molecular owners so this
    // test isolates the mechanical Jacobian of the bilateral topology.
    const std::array<int, 4> owner_myosin{{0, 1, 0, 1}};
    const std::array<int, 4> owner_half{{0, 1, 1, 0}};
    for (int a = 0; a < model.actin.n; ++a) {
        size_t best = model.explicit_titins.size();
        double best_distance = 1e100;
        for (size_t i = 0; i < model.explicit_titins.size(); ++i) {
            const auto& titin = model.explicit_titins[i];
            if (titin.myosin != owner_myosin[a] ||
                titin.half != owner_half[a]) {
                continue;
            }
            vec displacement = model.explicit_titin_thick_anchor(titin) -
                               static_cast<vec>(model.actin.left_end[a]);
            displacement.pbc_wrap(model.box, model.is_periodic);
            if (displacement.norm() < best_distance) {
                best = i;
                best_distance = displacement.norm();
            }
        }
        assert(best < model.explicit_titins.size());
        auto& titin = model.explicit_titins[best];
        titin.bound = true;
        titin.bound_actin = a;
        titin.last_bound_actin = a;
        model.explicit_titin_owner_by_actin[static_cast<size_t>(a)] =
            static_cast<int>(best);
    }
    model.update_system();

    // One individual barbed-end domain is captured on every actin. Extra
    // molecular tracks remain available; no actin is multiply owned.
    const int initial_bound = bound_titins(model);
    assert(initial_bound == 4);
    for (int a = 0; a < model.actin.n; ++a) {
        assert(model.explicit_titin_owner_by_actin[static_cast<size_t>(a)] >= 0);
    }

    const std::array<vec, 2> base_centers{{
        model.myosin.center[0], model.myosin.center[1]}};
    const double dx = 0.025;
    auto force_after_shift = [&](double shift) {
        model.myosin.center[0] = base_centers[0] + vec{shift, 0.0, 0.0};
        model.myosin.center[1] = base_centers[1];
        model.myosin.update_endpoints();
        model.update_system();
        vec total{0.0, 0.0, 0.0};
        for (int a = 0; a < model.actin.n; ++a) total += model.actin.force[a];
        for (int m = 0; m < model.myosin.n; ++m) total += model.myosin.force[m];
        assert(total.norm() < 1e-9);
        return model.myosin.force[0].x;
    };

    const double force_plus = force_after_shift(dx);
    const double force_minus = force_after_shift(-dx);
    const double slope = (force_plus - force_minus) / (2.0 * dx);
    std::printf("explicit titin local slope: F(+dx)=%.6f pN "
                "F(-dx)=%.6f pN dF/dx=%.6f pN/um\n",
                force_plus, force_minus, slope);

    // Moving a thick filament toward either Z-side must increase tension on
    // the opposite molecule and pull it back. This is the local mode that the
    // ACTN2-only cell leaves unconstrained.
    assert(force_plus < 0.0);
    assert(force_minus > 0.0);
    assert(slope < -50.0);

    // The deliberately buried geometry is not represented by an invented
    // compression spring.  After a large shift, the near-side molecule must
    // be exactly slack while the already-bound far-side molecule carries the
    // restoring tension.  The ownership is preserved across the displacement.
    const double buried_force = force_after_shift(0.55);
    int slack_links = 0;
    int tension_links = 0;
    for (const auto& titin : model.explicit_titins) {
        if (!titin.bound || titin.myosin != 0) continue;
        if (titin.force <= 1e-12) ++slack_links;
        else ++tension_links;
    }
    assert(buried_force < 0.0);
    assert(slack_links == 1);
    assert(tension_links == 1);
    assert(bound_titins(model) == initial_bound);

    gsl_rng_free(rng);
    std::remove(filename.c_str());
    return 0;
}
