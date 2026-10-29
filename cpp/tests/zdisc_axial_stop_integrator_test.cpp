#include "langevin.h"
#include "sarcomere.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <omp.h>

// The axial-stop connector force through the real update_system() +
// Langevin::sample_step() loop, in the controlled rig of
// analysis/standalone_actn2_axial_contact.py: connectors on four radial
// actin pairs at the X = 0.75 um plane, one 1.5 um myosin starting 25 nm
// past it, a constant 40 pN axial drive.  Noise is off, actins are held
// fixed and myosin rotation is off (the standalone rig has no rotations).
// Then: a scene evaluated with 1 and 4 OpenMP threads gives identical forces.

using utils::vec;

namespace {
const std::vector<double> box{6.0, 2.0, 2.0};
constexpr double kPi = 3.14159265358979323846;
constexpr double kPlane = 0.75, kDrive = 40.0;

std::unique_ptr<Sarcomere> make_model(gsl_rng* rng, int& seed, int n_actins,
                                      int n_myosins, const std::string& name) {
    int na = n_actins, nm = n_myosins, fixed = 0;
    std::vector<double> domain = box;
    double la = 1.0, lm = 1.5, rm = 0.05, am_cut = 0.05;
    double am_opt = 0.08, aa_cut = 0.08, aa_opt = 0.04;
    double zero = 0.0, life = 1e9, ratio = 1.0, dt = 1e-6;
    bool directional = true;
    std::string filename = "/tmp/zdisc_axial_stop_integrator_" + name + ".h5";
    auto model = std::make_unique<Sarcomere>(
        na, nm, domain, la, lm, rm, am_cut, am_opt, aa_cut, aa_opt,
        zero, zero, life, zero, zero, ratio,
        zero, zero, zero, zero, zero, zero, filename, rng, seed, fixed,
        dt, zero, zero, zero, directional, 24, zero, zero, zero, zero);
    model->set_periodicity({false, false, false});
    model->aa_local_multilink = true;
    model->aa_prebundle_links_per_pair = 2;
    model->aa_neighbor_identity_enabled = false;
    model->aa_barbed_zone = 0.0;
    model->aa_antiparallel_overlap = 1.0;
    model->aa_link_rupture_length = 0.0;
    model->aa_overlap_occlusion_veto = false;
    model->aa_myosin_contact_epsilon = 0.0;
    model->aa_myosin_contact_shell_energy = 0.0;
    model->aa_zdisc_tip_contact = true;
    model->aa_zdisc_tip_hard_projection = false;
    model->sterics_enabled = false;
    return model;
}

// One antiparallel pair at azimuth theta near a tip whose outward axis is
// tip_sign * x: the inner (80 nm) actin's left (barbed) end lies beyond the
// tip, outside the myosin capsule; the outer (120 nm) actin is antiparallel.
void place_pair(Sarcomere& model, int first, double theta,
                double x_center = 0.73, double tip_sign = 1.0) {
    const double c = std::cos(theta), s = std::sin(theta);
    model.actin.direction[first] = {-tip_sign, 0, 0};
    model.actin.center[first] = {x_center + 0.1 * tip_sign, 0.080 * c, 0.080 * s};
    model.actin.direction[first + 1] = {tip_sign, 0, 0};
    model.actin.center[first + 1] = {x_center - 0.1 * tip_sign, 0.120 * c, 0.120 * s};
}

void set_sites(Sarcomere& model, size_t slot, double x) {
    auto& link = model.aa_links[slot];
    link.s_i = (x - model.actin.center[link.actin_i].x) /
               model.actin.direction[link.actin_i].x;
    link.s_j = (x - model.actin.center[link.actin_j].x) /
               model.actin.direction[link.actin_j].x;
}

// Axial force of one radial connector whose points all sit at axial
// coordinate z from the myosin center (gate 1 inside 130 nm).
double connector_force(const Sarcomere& m, double z) {
    const double z_eps = std::hypot(z, m.aa_zdisc_tip_stop_epsilon);
    const double depth = 0.5 * 1.5 + m.aa_zdisc_tip_onset - z_eps;
    if (depth <= 0.0) return 0.0;
    return std::min(m.aa_zdisc_tip_stop_stiffness * depth,
                    m.aa_zdisc_tip_stop_force_max) * z / z_eps;
}

struct RigResult {
    vec center{0, 0, 0}, direction{1, 0, 0};
    double max_lateral = 0.0;            // |center y,z| over the run
    double max_contact_off_axis = 0.0;   // contact force normal to own axis
    double final_contact_axial = 0.0;
    int active_links = 0;
};

// Run the rig with `pairs` radial connectors (azimuths k*90 deg, k < pairs).
RigResult run_rig(gsl_rng* rng, int& seed, int pairs, bool rotate,
                  int steps, const std::string& name) {
    omp_set_num_threads(1);
    auto model = make_model(rng, seed, 2 * pairs, 1, name);
    model->aa_prebundle_links_per_pair = 1;
    for (int k = 0; k < pairs; ++k) place_pair(*model, 2 * k, k * kPi / 2);
    model->actin.update_endpoints();
    model->myosin.center[0] = {0.025, 0, 0};  // tip 25 nm past the plane
    model->myosin.direction[0] = {1, 0, 0};
    model->myosin.update_endpoints();
    for (int k = 0; k < pairs; ++k) model->prebind_local_aa_pair(2 * k, 2 * k + 1, 1);
    assert(static_cast<int>(model->aa_active_link_slots.size()) == pairs);
    for (size_t slot : model->aa_active_link_slots) set_sites(*model, slot, kPlane);

    // Myosin mobility 19 um/(pN s) as in the standalone rig; actins fixed.
    double beta = 1.0, dt = 1e-6, D_at = 0.0, D_ar = 0.0;
    double D_mt = 19.0, D_mr = rotate ? 2.0 : 0.0;
    int save_every = 1 << 30, fix_myosin = 0;
    bool resume = false;
    Langevin langevin(*model, beta, dt, D_at, D_ar, D_mt, D_mr, save_every,
                      resume, true, 1.0, 1.0, 1.0, 1.0);
    langevin.set_thermal_noise_scale(0.0);
    RigResult out;
    for (int step = 0; step < steps; ++step) {
        model->update_system();
        // No other force acts on the myosin in this fixture: this is the
        // contact force alone, before the drive is added.
        vec axis = model->myosin.direction[0];
        axis.normalize();
        const vec contact = model->myosin.force[0];
        out.max_contact_off_axis = std::max(out.max_contact_off_axis,
            (contact - contact.dot(axis) * axis).norm());
        out.final_contact_axial = contact.dot(axis);
        model->myosin.force[0] += vec{kDrive, 0, 0};
        langevin.sample_step(dt, rng, fix_myosin);
        const vec c = model->myosin.center[0];
        out.max_lateral = std::max(out.max_lateral, std::hypot(c.y, c.z));
    }
    out.center = model->myosin.center[0];
    out.direction = model->myosin.direction[0];
    out.active_links = static_cast<int>(model->aa_active_link_slots.size());
    // Sites are fixed material coordinates through the run.
    for (size_t slot : model->aa_active_link_slots) {
        const auto& link = model->aa_links[slot];
        const double xi = model->actin.center[link.actin_i].x +
            link.s_i * model->actin.direction[link.actin_i].x;
        assert(std::abs(xi - kPlane) < 1e-12);
    }
    return out;
}

// Myosin center where n radial connectors at the plane balance the drive.
double equilibrium_center(const Sarcomere& m, int n) {
    double lo = 0.70, hi = 0.79999;  // z = plane - center; force decreases in z
    for (int it = 0; it < 200; ++it) {
        const double mid = 0.5 * (lo + hi);
        (n * connector_force(m, mid) > kDrive ? lo : hi) = mid;
    }
    return kPlane - 0.5 * (lo + hi);
}
}  // namespace

int main() {
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    int seed = 11;
    // Parameter defaults for the analytic reference.
    auto reference = make_model(rng, seed, 2, 1, "reference");

    // 1. Four connectors hold the myosin short of the plane, where the
    //    contact exactly balances the drive; nothing moves it sideways.
    const auto four = run_rig(rng, seed, 4, false, 4000, "four");
    const double expected4 = equilibrium_center(*reference, 4);
    assert(four.active_links == 4);
    assert(std::abs(four.center.x - expected4) < 1e-7);
    assert(std::abs(four.final_contact_axial + kDrive) < 1e-6);
    assert(four.max_lateral == 0.0);
    assert(four.max_contact_off_axis == 0.0);

    // 2. One connector lost (three left, asymmetric): still holds, deeper,
    //    and the myosin moves exactly zero sideways.  Two left also hold.
    const auto three = run_rig(rng, seed, 3, false, 4000, "three");
    const double expected3 = equilibrium_center(*reference, 3);
    assert(three.active_links == 3);
    assert(std::abs(three.center.x - expected3) < 1e-7);
    assert(three.max_lateral == 0.0);
    assert(three.max_contact_off_axis == 0.0);
    const auto two = run_rig(rng, seed, 2, false, 4000, "two");
    assert(std::abs(two.center.x - equilibrium_center(*reference, 2)) < 1e-7);
    assert(two.max_lateral == 0.0);
    // One connector (30 pN capacity < 40 pN drive) gives way straight
    // along the axis.
    const auto one = run_rig(rng, seed, 1, false, 4000, "one");
    assert(one.center.x > 0.1);
    assert(one.max_lateral == 0.0);
    // With myosin rotation on, the asymmetric contact tilts the myosin, but
    // the contact force stays along the myosin's own axis.
    const auto tilting = run_rig(rng, seed, 3, true, 4000, "tilting");
    vec axis = tilting.direction;
    axis.normalize();
    const double tilt_deg = std::acos(std::min(1.0, std::abs(axis.x))) * 180.0 / kPi;
    assert(tilting.max_contact_off_axis < 1e-9);

    // 3. Thread-count determinism: identical scenes evaluated with 1 and 4
    //    threads give identical contact forces, torques and tip records.
    //    Two myosins, ten actin pairs with two redundant links each.
    struct Snapshot {
        std::vector<vec> actin_force, actin_torque, myosin_force, myosin_torque;
        std::vector<double> tip_force, tip_gap;
    };
    auto evaluate = [&](int threads) {
        omp_set_num_threads(threads);
        int local_seed = 23;
        auto model = make_model(rng, local_seed, 20, 2,
                                "threads" + std::to_string(threads));
        for (int k = 0; k < 10; ++k) {
            // Five pairs around each myosin tip: +x tip of myosin 0 at the
            // origin, -x tip of myosin 1 centered at x = -1.6.
            const double theta = (k % 5) * 2 * kPi / 5 + 0.3 * (k / 5);
            place_pair(*model, 2 * k, theta, k < 5 ? 0.73 : -2.33,
                       k < 5 ? 1.0 : -1.0);
        }
        model->actin.update_endpoints();
        model->myosin.center[0] = {0.0, 0.0, 0.0};
        model->myosin.direction[0] = {1, 0, 0};
        model->myosin.center[1] = {-1.6, 0.0, 0.0};
        model->myosin.direction[1] = {1, 0, 0};
        model->myosin.update_endpoints();
        for (int k = 0; k < 10; ++k) model->prebind_local_aa_pair(2 * k, 2 * k + 1, 2);
        assert(model->aa_active_link_slots.size() == 20);
        int index = 0;
        for (size_t slot : model->aa_active_link_slots) {
            const auto& link = model->aa_links[slot];
            const bool near_plus = link.actin_i < 10;
            // Two redundant links per pair at different depths.
            const double depth = (index++ % 2 == 0) ? 0.02 : 0.035;
            set_sites(*model, slot, near_plus ? 0.75 - depth : -2.35 + depth);
        }
        model->update_system();
        Snapshot snap;
        for (int a = 0; a < model->actin.n; ++a) {
            snap.actin_force.push_back(model->actin.force[a]);
            snap.actin_torque.push_back(model->actin.torque[a]);
        }
        for (int m = 0; m < model->myosin.n; ++m) {
            snap.myosin_force.push_back(model->myosin.force[m]);
            snap.myosin_torque.push_back(model->myosin.torque[m]);
            for (const char* name : {"aa_zdisc_tip_force_plus", "aa_zdisc_tip_force_minus"})
                snap.tip_force.push_back(model->myosin[name][m]);
            for (const char* name : {"aa_zdisc_tip_gap_plus", "aa_zdisc_tip_gap_minus"})
                snap.tip_gap.push_back(model->myosin[name][m]);
        }
        return snap;
    };
    const auto serial = evaluate(1);
    const auto parallel = evaluate(4);
    double worst_thread = 0.0;
    for (size_t a = 0; a < serial.actin_force.size(); ++a) {
        worst_thread = std::max({worst_thread,
            (serial.actin_force[a] - parallel.actin_force[a]).norm(),
            (serial.actin_torque[a] - parallel.actin_torque[a]).norm()});
    }
    for (size_t m = 0; m < serial.myosin_force.size(); ++m) {
        worst_thread = std::max({worst_thread,
            (serial.myosin_force[m] - parallel.myosin_force[m]).norm(),
            (serial.myosin_torque[m] - parallel.myosin_torque[m]).norm()});
    }
    assert(worst_thread < 1e-12);
    for (size_t k = 0; k < serial.tip_force.size(); ++k) {
        assert(serial.tip_force[k] == parallel.tip_force[k]);
        const bool both_nan = std::isnan(serial.tip_gap[k]) && std::isnan(parallel.tip_gap[k]);
        assert(both_nan || serial.tip_gap[k] == parallel.tip_gap[k]);
    }
    // Both tips are in contact, each with five pairs' worth of force.
    assert(serial.tip_force[0] > 100.0 && serial.tip_force[3] > 100.0);

    gsl_rng_free(rng);
    std::cout << "zdisc axial-stop integrator: 4 connectors stop at "
              << 1000 * (four.center.x) << " nm (analytic " << 1000 * expected4
              << "), 3 at " << 1000 * three.center.x << " nm (analytic "
              << 1000 * expected3 << "), 1 slides to " << 1000 * one.center.x
              << " nm, sideways 0; tilting run " << tilt_deg
              << " deg with on-axis contact; 1 vs 4 threads max diff "
              << worst_thread << "\n";
    return 0;
}
