// Myomesin M-band spring.
//
// This exists because the spring was previously computed and then discarded:
// _apply_myomesin_spring filled local force vectors that were never written
// back, so --k_mm silently did nothing and the M-band mechanism had never run
// in any trajectory. Every assertion below is about the force actually
// ARRIVING on the myosins, which is exactly what the old code failed to do and
// what no existing test checked.
//
// Physically the spring draws the centres of two transversely adjacent thick
// filaments into a common axial plane. That is an M-band; combined with ACTN2
// at the barbed ends it supplies the second registration plane a sarcomere
// needs.

#include "sarcomere.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>

namespace {

std::unique_ptr<Sarcomere> make_model(
    const std::string& output, gsl_rng* rng, int& seed, double k_mm) {
    int n_actins = 1;
    int n_myosins = 2;
    std::vector<double> box{8.0, 8.0, 8.0};
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
    double dt = 1e-6;
    bool directional = true;
    int max_myosin_bonds = 6;
    std::string filename = output;   // constructor takes a non-const reference
    auto model = std::make_unique<Sarcomere>(
        n_actins, n_myosins, box, actin_length, myosin_length,
        myosin_radius, am_cutoff, am_optimal, aa_cutoff, aa_optimal,
        zero, zero, zero, zero, zero, diff_ratio, zero, zero, zero,
        zero, k_mm, zero, filename, rng, seed, fixed_myosins, dt, zero,
        zero, zero, directional, max_myosin_bonds, zero, zero, zero,
        zero);
    model->set_periodicity({false, false, false});
    model->dynamic_am = false;
    model->aa_local_multilink = false;
    model->sterics_enabled = true;
    // Park the actin far away so only the myosin pair is under test.
    model->actin.center[0] = {0.0, 3.5, 3.5};
    model->actin.direction[0] = {1.0, 0.0, 0.0};
    model->actin.update_endpoints();
    return model;
}

// Two parallel thick filaments, transversely adjacent, offset axially by
// `axial`. The caller gates on segment distance <= 2*am_cutoff = 0.20 um.
void place(Sarcomere& model, double axial, double transverse) {
    model.myosin.direction[0] = {1.0, 0.0, 0.0};
    model.myosin.direction[1] = {1.0, 0.0, 0.0};
    model.myosin.center[0] = {0.0, 0.0, 0.0};
    model.myosin.center[1] = {axial, transverse, 0.0};
    model.reset_myosin_material_frames();
    model.myosin.update_endpoints();
}

}  // namespace

int main() {
    const std::string output = "/tmp/myomesin_mband_test.h5";
    int seed = 5501;
    constexpr double K = 100.0;
    constexpr double TRANSVERSE = 0.08;   // inside 2*am_cutoff

    // 1. The force arrives, is attractive along the axis, and is reciprocal.
    //    A regression on the discarded-force bug: with the old code every one
    //    of these would read exactly zero.
    {
        std::remove(output.c_str());
        gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng, seed);
        auto model = make_model(output, rng, seed, K);
        place(*model, 0.40, TRANSVERSE);
        model->update_system();
        const vec f0 = model->myosin.force[0];
        const vec f1 = model->myosin.force[1];
        assert(f0.norm() > 1e-9);
        assert(f1.norm() > 1e-9);
        // myosin 1 sits at +x of myosin 0, so myosin 0 is pulled toward +x.
        assert(f0.x > 0.0);
        assert(f1.x < 0.0);
        // Reciprocal: the pair exerts no net force on the system.
        assert(std::abs(f0.x + f1.x) < 1e-9);
        gsl_rng_free(rng);
    }

    // 2. k_mm = 0 leaves the myosins untouched, so the term is opt-in and the
    //    previous behaviour is exactly recoverable.
    {
        std::remove(output.c_str());
        gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng, seed);
        auto model = make_model(output, rng, seed, 0.0);
        place(*model, 0.40, TRANSVERSE);
        model->update_system();
        assert(model->myosin.force[0].norm() < 1e-12);
        assert(model->myosin.force[1].norm() < 1e-12);
        gsl_rng_free(rng);
    }

    // 3. Magnitude grows with axial offset -- the law is a zero-rest-length
    //    spring, so coplanar filaments feel nothing and offset ones are pulled
    //    back together.
    {
        std::remove(output.c_str());
        gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng, seed);
        auto model = make_model(output, rng, seed, K);
        place(*model, 0.20, TRANSVERSE);
        model->update_system();
        const double near = model->myosin.force[0].x;
        place(*model, 0.60, TRANSVERSE);
        model->update_system();
        const double far = model->myosin.force[0].x;
        assert(far > near);
        assert(near > 0.0);
        gsl_rng_free(rng);
    }

    // 4. Beyond the axial cutoff (one myosin length) the pair does not couple.
    //    This is what keeps the spring from collapsing every thick filament
    //    onto a single plane: separate M-bands stay independent.
    {
        std::remove(output.c_str());
        gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng, seed);
        auto model = make_model(output, rng, seed, K);
        place(*model, 1.90, TRANSVERSE);   // > myosin_length = 1.5
        model->update_system();
        assert(model->myosin.force[0].norm() < 1e-12);
        gsl_rng_free(rng);
    }

    // 5. Transversely distant filaments do not couple either, so the spring is
    //    a local neighbour rule rather than a global centroid device.
    {
        std::remove(output.c_str());
        gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng, seed);
        auto model = make_model(output, rng, seed, K);
        place(*model, 0.40, 1.20);   // far beyond 2*am_cutoff
        model->update_system();
        assert(model->myosin.force[0].norm() < 1e-12);
        gsl_rng_free(rng);
    }

    std::remove(output.c_str());
    return 0;
}
