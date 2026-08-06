#include "sarcomere.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>

namespace {

std::unique_ptr<Sarcomere> make_model(
    const std::string& output, gsl_rng* rng, int& seed) {
    int n_actins = 1;
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
    model->configure_explicit_titin(3, 0.70, 1e12, 0.0, false, 0.5);
    model->myosin.center[0] = {0.0, 0.0, 0.0};
    model->myosin.direction[0] = {1.0, 0.0, 0.0};
    model->reset_myosin_material_frames();
    model->myosin.update_endpoints();

    // Put the individual actin barbed endpoint 0.60 um behind track zero on
    // the positive half, using the transported thick-filament material frame.
    const auto positive_track = static_cast<size_t>(3);
    const vec anchor = model->explicit_titin_thick_anchor(
        model->explicit_titins[positive_track]);
    model->actin.direction[0] = {1.0, 0.0, 0.0};
    model->actin.center[0] = anchor - vec{0.60, 0.0, 0.0} +
                             vec{0.50, 0.0, 0.0};
    model->actin.update_endpoints();
    return model;
}

size_t bound_titin(const Sarcomere& model) {
    for (size_t index = 0; index < model.explicit_titins.size(); ++index) {
        if (model.explicit_titins[index].bound) return index;
    }
    return model.explicit_titins.size();
}

}  // namespace

int main() {
    const std::string output = "/tmp/explicit_titin_test.h5";
    std::remove(output.c_str());
    int seed = 9731;
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, seed);
    auto model = make_model(output, rng, seed);

    assert(model->explicit_titins.size() == 6);
    model->update_system();
    const size_t owner = bound_titin(*model);
    assert(owner < model->explicit_titins.size());
    const auto& stretched = model->explicit_titins[owner];
    assert(stretched.myosin == 0);
    assert(stretched.half == 1);
    assert(stretched.bound_actin == 0);
    assert(std::abs(stretched.distance - 0.60) < 1e-10);
    assert(std::abs(stretched.extension - 0.10) < 1e-10);
    assert(std::abs(stretched.force - 4.0) < 1e-10);
    assert(std::abs(model->actin.force[0].norm() - 4.0) < 1e-10);
    assert((static_cast<vec>(model->actin.force[0]) +
            static_cast<vec>(model->myosin.force[0])).norm() < 1e-10);
    assert(model->explicit_titin_owner_by_actin[0] ==
           static_cast<int>(owner));

    // Compression retains molecular ownership but contributes no mechanical
    // force, so titin cannot push the actin away from the thick filament.
    const vec anchor = model->explicit_titin_thick_anchor(stretched);
    model->actin.center[0] = anchor - vec{0.40, 0.0, 0.0} +
                             vec{0.50, 0.0, 0.0};
    model->actin.update_endpoints();
    model->update_system();
    assert(model->explicit_titins[owner].bound);
    assert(std::abs(model->explicit_titins[owner].distance - 0.40) < 1e-10);
    assert(std::abs(model->explicit_titins[owner].force) < 1e-12);
    assert(model->actin.force[0].norm() < 1e-12);
    assert(model->myosin.force[0].norm() < 1e-12);

    // Both checkpoint layers must preserve the molecular owner exactly.
    model->new_file();
    model->save_state();
    model->save_resume_snapshot();
    gsl_rng* restored_rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(restored_rng, seed);
    auto restored = make_model(output, restored_rng, seed);
    int frames = 0;
    restored->load_state(frames);
    const size_t restored_owner = bound_titin(*restored);
    assert(restored_owner == owner);
    assert(restored->explicit_titins[restored_owner].bound_actin == 0);
    assert(restored->explicit_titin_owner_by_actin[0] ==
           static_cast<int>(restored_owner));
    assert(std::abs(restored->explicit_titins[restored_owner].distance - 0.40) <
           1e-10);

    gsl_rng_free(restored_rng);
    gsl_rng_free(rng);
    std::remove(output.c_str());
    return 0;
}
