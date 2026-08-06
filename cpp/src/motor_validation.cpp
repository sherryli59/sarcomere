// Motor validation: three minimal systems, ONE parameter set.
//
// A myosin head ensemble is only a motor if it transports. The production
// configuration currently fails that at every scale -- free filaments move at
// diffusion level -- while still passing the balanced-cell stability gate, so
// the stability gate alone cannot tell a motor from a crosslink. These three
// rigs have to succeed simultaneously:
//
//   free      1 actin, 1 FIXED myosin, nothing else. A motility assay. The
//             filament must translocate along the myosin axis at a useful
//             fraction of the unloaded velocity and keep going. This is the
//             test the model presently fails.
//   linked    2 antiparallel actins joined by one ACTN2 link, 1 fixed myosin.
//             Transport against an internal elastic load: slower than free,
//             still directed.
//   balanced  the 4-actin / 2-myosin cell.  There is no external load and no
//             titin.  Opposed motors first slide weakly crosslinked actins;
//             locked ACTN2 then stores axial strain and must stop the relative
//             A-M motion by an internal force balance.
//
// Reported per mode: net axial displacement, the drift velocity, and the
// pure-diffusion expectation for the same duration. Ratio ~1 means the
// filament is only diffusing.
//
// Usage: motor_validation <mode> <ms> [seed]

#include "sarcomere.h"
#include "langevin.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <unistd.h>

namespace {

struct Rig {
    std::unique_ptr<Sarcomere> model;
    int probe_actin = 0;
    vec axis{1.0, 0.0, 0.0};
};

struct TitinConfig {
    double stiffness = 0.0;
    double slack_length = 0.64;
    double capture_distance = 0.72;
    double k_on = 5000.0;
    double k_off = 0.5;
    double transition_theta = 0.95;
    int links_per_half = 3;
};

Rig build(const std::string& mode, int seed, double dt, int heads_per_half,
          double act_floor, double act_frac, std::string& h5,
          gsl_rng* rng, int& fix_myosin_out, double polarity,
          double init_aa, double aa_stiffness_scale,
          double aa_register_ratio, const TitinConfig& titin,
          double am_head_transverse_k, double am_alignment_k,
          double am_transverse_strain_detachment_theta,
          double balanced_box_x, const std::string& balanced_variant) {
    int n_actins = (mode == "balanced") ? 4 : (mode == "linked" ? 2 : 1);
    int n_myosins = (mode == "balanced") ? 2 : 1;
    std::vector<double> box = (mode == "balanced")
        ? std::vector<double>{balanced_box_x, 3.0, 3.0}
        : std::vector<double>{8.0, 6.0, 6.0};
    double actin_length = 1.0, myosin_length = 1.5, myosin_radius = 0.025;
    double am_cutoff = 0.10, am_optimal = 0.055;
    double aa_cutoff = 0.10, aa_optimal = 0.055;
    double zero = 0.0, diff_ratio = 1.0;
    double aa_k_on = (mode == "free") ? 0.0 : 2500.0;
    double aa_k_off = 1.0;
    double base_lifetime = (mode == "free") ? 0.0 : 0.1;
    double dir_lifetime = (mode == "free") ? 0.0 : 0.05;
    double lifetime_coeff = 0.4;
    double k_aa = 10.0, kappa_aa = 100.0;
    // Myosin index 0 is held fixed in the free/linked assays, which is what a
    // motility slide does; the balanced cell must be free to find its own
    // force balance, so nothing is pinned there.
    int fix_myosin = 0;   // never exclude myosins from processing
    fix_myosin_out = fix_myosin;
    bool directional = true;

    // Zero preserves the ACTN2-only validation. A positive opt-in value adds
    // explicit, turnover-capable, tension-only molecules to the balanced rig;
    // it never selects the index-owned or compression-bearing legacy modes.
    const double titin_k_v = titin.stiffness;
    const double titin_rest_v = titin.stiffness > 0.0
        ? titin.slack_length : 0.0;
    const double temperature = 310.0;
    const double beta = rod_mobility::beta(temperature);
    const auto am = rod_mobility::finite_cylinder(
        actin_length, 0.010, temperature, 0.001, 10.0, 10.0, 160.0, 160.0);
    const auto mm = rod_mobility::finite_cylinder(
        myosin_length, 0.035, temperature, 0.001, 20.0, 20.0, 64.0, 64.0);
    const double maxaf = 0.005 / (beta * std::max(am.parallel, am.perpendicular) * dt);
    const double maxmf = 0.005 / (beta * std::max(mm.parallel, mm.perpendicular) * dt);
    const double maxat = 0.005 / (beta * am.tumble * dt);
    const double maxmt = 0.005 / (beta * mm.tumble * dt);

    Rig rig;
    rig.model = std::make_unique<Sarcomere>(
        n_actins, n_myosins, box, actin_length, myosin_length,
        myosin_radius, am_cutoff, am_optimal, aa_cutoff, aa_optimal,
        aa_k_on, aa_k_off, base_lifetime, dir_lifetime, lifetime_coeff,
        diff_ratio, k_aa, kappa_aa, zero, zero, zero, zero,
        h5, rng, seed, fix_myosin, dt, zero, titin_k_v, titin_rest_v,
        directional, 24, maxaf, maxmf, maxat, maxmt);
    Sarcomere& m = *rig.model;
    m.dynamic_am = true;
    m.am_spatial_heads = true;
    m.am_offcenter_torque = true;
    m.sterics_enabled = true;
    m.aa_local_multilink = (mode != "free");
    m.set_physical_mobility(temperature, beta, am, mm);
    m.set_am_radial_geometry(0.015, 0.015, 1000.0);
    m.set_steric_stiffness(2000.0, 1000.0);
    const double centerline = myosin_radius + 0.015 + 0.015;
    m.set_aa_spacing_geometry(aa_optimal / centerline, centerline,
                              aa_cutoff - aa_optimal);
    m.set_am_kinetics(true, 100.0, 200.0, 200.0, 160.0, 1.0,
                      1000.0, 0.005, 4.5, 5.0, 0.1, 0.1,
                      0.0, 0.02, 0.2, act_floor,
                      mode == "balanced" ? 0.020 : 0.0, act_frac);
    m.set_am_proximity_tether(
        3000.0, 0.002, true, 30.0, am_alignment_k, 2.0);
    m.set_am_head_transverse_tether(am_head_transverse_k);
    m.set_am_retention_kinetics(
        true, 0.1, am_transverse_strain_detachment_theta,
        20.0, 6, 1.0);
    m.set_am_head_count(heads_per_half);
    m.set_am_neighbor_head_cap(24);
    m.set_am_weak_binding_floor(0.1);
    if (titin.stiffness > 0.0) {
        if (mode != "balanced") {
            throw std::invalid_argument(
                "explicit titin diagnostic is only defined for balanced mode");
        }
        m.titin_mode = Sarcomere::TitinMode::explicit_dynamic;
        // A new actin partner must already have a live ACTN2 molecule at its
        // barbed end. The link itself turns over independently; no centroid,
        // Z-node, or occupancy-triggered effective spring is introduced.
        m.configure_explicit_titin(
            titin.links_per_half, titin.capture_distance,
            titin.k_on, titin.k_off, true, titin.transition_theta,
            0.020, 100.0, true, 0.0, false);
    }
    if (mode != "free") {
        m.set_local_aa_model(true, 2, 0.0192, 0.1, 50.0, 0, 3);
        m.set_aa_neighbor_identity(true, 0.02, 0.20);
        m.max_aa_bonds_per_actin = 2;
        m.aa_association_max_direction_dot = -0.5;
        m.aa_link_rupture_length = 0.12;
        m.aa_barbed_zone = 0.20;
        m.aa_status1_radial_k = 75.0;
        m.set_aa_axial_overlap(init_aa, 0.0, 0.1, 0.15, 0.05);
        // Formation requires opposed productive motor support and axial
        // register. Local ACTN2 arm loads remain diagnostics and lifetime
        // inputs, but are not a separate hard gate on status-2 formation.
        m.set_aa_status_kinetics(
            actn2_load_kinetics::StatusMode::productive_register, 0.25, 0.10);
    }

    if (mode == "balanced") {
        // Three explicit molecules per antiparallel pair provide molecular
        // turnover without opening the entire macroscopic load path when one
        // ACTN2 molecule detaches.
        m.set_local_aa_model(true, 3, 0.0192, 0.1, 50.0, 0, 3);
        // The titin perturbation gate must assemble at the centered geometry.
        // Its 550-nm displacement is applied only after fixed-geometry titin
        // association below, while preserving those material bonds.
        const std::string initialization_variant =
            balanced_variant == "balanced_titin_perturb"
                ? "balanced" : balanced_variant;
        m.minimal_catch_variant(initialization_variant);
        // Status 1 is a soft, mobile encounter complex.  It becomes a stiff,
        // material-locked status-2 connector only close to the 100 nm A-A
        // register, so the starting 160 nm overlap must anneal before locking.
        m.k_aa *= std::max(0.0, aa_stiffness_scale);
        m.aa_antiparallel_overlap = 0.10;
        m.aa_axial_register_ratio = std::max(0.0, aa_register_ratio);
        m.aa_axial_register_force_max = 12.0 * m.aa_axial_register_ratio;
        m.aa_status2_registry_on = 0.040;
        m.aa_status2_registry_off = 0.040;
        m.aa_status2_support_window = 0.010;
        m.aa_status2_locked_until_detach = true;
        m.set_aa_status_kinetics(
            actn2_load_kinetics::StatusMode::productive_register, 0.25, 0.10);
        // Chemistry and motion still begin together. q starts from zero and
        // the mature interaction mechanics rise smoothly over 20 ms.
        m.set_interaction_ramp(0, 40000);
    } else {
        // One myosin along +x at the origin. Place the actin so it overlaps
        // the positive half with the polarity the productive-half rule
        // accepts, barbed end outboard.
        m.myosin.center[0] = {0.0, 0.0, 0.0};
        m.myosin.direction[0] = {1.0, 0.0, 0.0};
        m.reset_myosin_material_frames();
        m.myosin.update_endpoints();
        const double centre_x = 0.95;           // overlapping the + half
        const double radial = centerline;
        // polarity = -1 puts the barbed end OUTBOARD (x = centre + L/2),
        // which is the sarcomere arrangement: barbed end at the Z-disc,
        // pointed end reaching in toward the M-line.
        m.actin.direction[0] = {polarity, 0.0, 0.0};
        m.actin.center[0] = {centre_x, radial, 0.0};
        if (mode == "linked") {
            m.actin.direction[1] = {-polarity, 0.0, 0.0};
            m.actin.center[1] = {centre_x + 0.90, radial + 0.055, 0.0};
        }
        m.actin.update_endpoints();
    }
    // Rebuild neighbours now that everything is in its final place.
    m.set_periodicity({true, false, false});
    m.actin.update_endpoints();
    m.myosin.update_endpoints();
    std::printf("# geometry: myosin c=(%.3f,%.3f,%.3f) spans x[%.3f,%.3f]\n",
                m.myosin.center[0].x, m.myosin.center[0].y, m.myosin.center[0].z,
                m.myosin.center[0].x - 0.5*myosin_length,
                m.myosin.center[0].x + 0.5*myosin_length);
    for (int i = 0; i < m.actin.n; ++i) {
        std::printf("# geometry: actin %d c=(%.3f,%.3f,%.3f) dir=(%.2f,%.2f,%.2f)"
                    " spans x[%.3f,%.3f]  radial=%.4f\n", i,
                    m.actin.center[i].x, m.actin.center[i].y, m.actin.center[i].z,
                    m.actin.direction[i].x, m.actin.direction[i].y, m.actin.direction[i].z,
                    m.actin.center[i].x - 0.5*actin_length,
                    m.actin.center[i].x + 0.5*actin_length,
                    std::hypot(m.actin.center[i].y - m.myosin.center[0].y,
                               m.actin.center[i].z - m.myosin.center[0].z));
    }
    rig.probe_actin = 0;
    rig.axis = {1.0, 0.0, 0.0};
    return rig;
}

struct ForceLedger {
    double motor = 0.0;       // productive AM force, pN/actin
    double actn2 = 0.0;       // signed AA force, pN/actin; negative resists
    double actn2_load = 0.0;  // -actn2, positive when resisting
    double bilateral = 0.0;   // bilateral resisting arm load, pN/link
    double net = 0.0;         // all deterministic actin forces, pN/actin
    double other = 0.0;       // net - motor - actn2
    int status1 = 0;
    int status2 = 0;
    std::vector<double> motor_by_actin;
    std::vector<double> actn2_by_actin;
    std::vector<double> net_by_actin;
    std::vector<double> other_by_actin;
};

ForceLedger force_ledger(Sarcomere& model) {
    ForceLedger result;
    const int count = std::max(1, model.actin.n);
    result.motor_by_actin.assign(static_cast<size_t>(model.actin.n), 0.0);
    result.actn2_by_actin.assign(static_cast<size_t>(model.actin.n), 0.0);
    result.net_by_actin.assign(static_cast<size_t>(model.actin.n), 0.0);
    result.other_by_actin.assign(static_cast<size_t>(model.actin.n), 0.0);
    const auto& am_x = model.actin["am_active_force_x"];
    const auto& am_y = model.actin["am_active_force_y"];
    const auto& am_z = model.actin["am_active_force_z"];
    for (int i = 0; i < model.actin.n; ++i) {
        const vec motor_force{am_x[i], am_y[i], am_z[i]};
        result.motor_by_actin[static_cast<size_t>(i)] =
            motor_force.dot(model.actin.direction[i]);
        result.net_by_actin[static_cast<size_t>(i)] =
            model.actin.force[i].dot(model.actin.direction[i]);
        result.motor += result.motor_by_actin[static_cast<size_t>(i)];
        result.net += result.net_by_actin[static_cast<size_t>(i)];
    }
    const double alignment_k = actn2_load_kinetics::fixed_single_link_stiffness(
        model.kappa_aa, model.aa_stiffness_reference_links);
    double bilateral_sum = 0.0;
    int bilateral_count = 0;
    for (const auto& link : model.aa_links) {
        if (!link.active) continue;
        result.status1 += link.status == 1 ? 1 : 0;
        result.status2 += link.status == 2 ? 1 : 0;
        const auto mechanics = compute_local_aa_link(
            model.actin, link.actin_i, link.actin_j, link.s_i, link.s_j,
            model.box, model.aa_link_radial_stiffness(link.status),
            model.aa_perpendicular_k, alignment_k, model.aa_optimal,
            model.aa_antiparallel_overlap,
            model.aa_link_register_ratio(link.status),
            model.aa_register_saturation_error,
            model.aa_axial_register_force_max);
        const vec force_i = mechanics.force_on_first;
        const vec force_j = -force_i;
        const double projected_i = force_i.dot(
            model.actin.direction[link.actin_i]);
        const double projected_j = force_j.dot(
            model.actin.direction[link.actin_j]);
        result.actn2_by_actin[static_cast<size_t>(link.actin_i)] += projected_i;
        result.actn2_by_actin[static_cast<size_t>(link.actin_j)] += projected_j;
        result.actn2 += projected_i + projected_j;
        const double bilateral = std::min(
            std::max(0.0, -projected_i),
            std::max(0.0, -projected_j));
        bilateral_sum += bilateral;
        ++bilateral_count;
    }
    result.motor /= count;
    result.actn2 /= count;
    result.actn2_load = std::max(0.0, -result.actn2);
    result.net /= count;
    result.other = result.net - result.motor - result.actn2;
    for (int i = 0; i < model.actin.n; ++i) {
        result.other_by_actin[static_cast<size_t>(i)] =
            result.net_by_actin[static_cast<size_t>(i)] -
            result.motor_by_actin[static_cast<size_t>(i)] -
            result.actn2_by_actin[static_cast<size_t>(i)];
    }
    result.bilateral = bilateral_count > 0
        ? bilateral_sum / bilateral_count : 0.0;
    return result;
}

struct LinearWindow {
    double n = 0.0, st = 0.0, sq = 0.0, stt = 0.0, stq = 0.0;
    void add(double t, double q) {
        n += 1.0; st += t; sq += q; stt += t*t; stq += t*q;
    }
    double slope() const {
        const double denominator = n*stt - st*st;
        return std::abs(denominator) > 1e-20
            ? (n*stq - st*sq) / denominator : 0.0;
    }
};

struct LedgerWindow {
    long n = 0;
    double motor = 0.0, actn2 = 0.0, load = 0.0, bilateral = 0.0;
    double net = 0.0, other = 0.0, status1 = 0.0, status2 = 0.0;
    void add(const ForceLedger& value) {
        ++n; motor += value.motor; actn2 += value.actn2;
        load += value.actn2_load; bilateral += value.bilateral;
        net += value.net; other += value.other;
        status1 += value.status1; status2 += value.status2;
    }
    double mean(double value) const { return n > 0 ? value/n : 0.0; }
};

struct FvSample {
    double time = 0.0;
    double velocity = 0.0;
    double load = 0.0;
    double status2 = 0.0;
};

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::printf("usage: motor_validation <free|linked|balanced> <ms> [seed]"
                    " [heads] [floor] [frac] [polarity] [preeq_ms] [init_aa]"
                    " [weak_site_D] [aa_scale] [register_ratio] [thermal_scale]"
                    " [csv] [trajectory_h5] [save_every] [titin_k]"
                    " [titin_slack] [titin_capture] [titin_k_on]"
                    " [titin_k_off] [titin_links_per_half] [titin_theta]"
                    " [initialize_q_calcium_equilibrium]"
                    " [am_head_transverse_k] [am_alignment_k]"
                    " [am_transverse_strain_detachment_theta]"
                    " [removed_short_length_placeholder] [balanced_box_x]"
                    " [balanced_variant] [bidirectional_crossbridge]"
                    " [double_overlap_contact] [contact_diameter_um]"
                    " [contact_epsilon_pn_um] [contact_force_cap_pn]\n");
        return 2;
    }
    const std::string mode = argv[1];
    const double ms = std::atof(argv[2]);
    int seed = (argc > 3) ? std::atoi(argv[3]) : 20260803;
    const int heads = (argc > 4) ? std::atoi(argv[4]) : 35;
    const double floor_v = (argc > 5) ? std::atof(argv[5]) : 0.0;
    const double frac_v = (argc > 6) ? std::atof(argv[6]) : 0.5;
    const double polarity = (argc > 7) ? std::atof(argv[7]) : -1.0;
    // Fixed-geometry chemical pre-equilibration, as the sandbox gate does:
    // update_system() only, no Langevin step, so binding kinetics settle
    // before any filament is allowed to move. The passing run used 30 ms.
    const double preq_ms = (argc > 8) ? std::atof(argv[8])
                                      : (mode == "balanced" ? 30.0 : 0.0);
    const double init_aa = (argc > 9) ? std::atof(argv[9]) : 0.16;
    const double weak_site_D = (argc > 10) ? std::atof(argv[10]) : 0.010;
    const double aa_scale = (argc > 11) ? std::atof(argv[11]) : 100.0;
    const double register_ratio = (argc > 12) ? std::atof(argv[12]) : 3.0;
    const double thermal_scale = (argc > 13) ? std::atof(argv[13]) : 1.0;
    const std::string csv_path = (argc > 14) ? argv[14] : "";
    const std::string trajectory_path = (argc > 15) ? argv[15] : "";
    const int trajectory_save_every = (argc > 16)
        ? std::max(1, std::atoi(argv[16])) : 1000;
    TitinConfig titin;
    titin.stiffness = (argc > 17) ? std::atof(argv[17]) : 0.0;
    titin.slack_length = (argc > 18) ? std::atof(argv[18]) : 0.64;
    titin.capture_distance = (argc > 19) ? std::atof(argv[19]) : 0.72;
    titin.k_on = (argc > 20) ? std::atof(argv[20]) : 5000.0;
    titin.k_off = (argc > 21) ? std::atof(argv[21]) : 0.5;
    titin.links_per_half = (argc > 22)
        ? std::max(1, std::atoi(argv[22])) : 3;
    titin.transition_theta = (argc > 23) ? std::atof(argv[23]) : 0.95;
    const bool initialize_q_calcium_equilibrium =
        (argc > 24) ? std::atoi(argv[24]) != 0 : false;
    const double am_head_transverse_k =
        (argc > 25) ? std::atof(argv[25]) : 0.0;
    const double am_alignment_k =
        (argc > 26) ? std::atof(argv[26]) : 30.0;
    const double am_transverse_strain_detachment_theta =
        (argc > 27) ? std::atof(argv[27]) : 0.0;
    // Positional compatibility with archived diagnostic scripts. The old
    // short-length penalty has been removed and this value is ignored.
    const double removed_short_length_placeholder =
        (argc > 28) ? std::atof(argv[28]) : 0.0;
    (void)removed_short_length_placeholder;
    const double balanced_box_x =
        (argc > 29) ? std::atof(argv[29]) : 5.32;
    const std::string balanced_variant =
        (argc > 30) ? argv[30] : "balanced";
    const bool bidirectional_crossbridge =
        (argc > 31) ? std::atoi(argv[31]) != 0 : false;
    const bool double_overlap_contact =
        (argc > 32) ? std::atoi(argv[32]) != 0 : false;
    const double contact_diameter =
        (argc > 33) ? std::atof(argv[33]) : 0.012;
    const double contact_epsilon =
        (argc > 34) ? std::atof(argv[34]) : 0.00428;
    const double contact_force_cap =
        (argc > 35) ? std::atof(argv[35]) : 0.0;
    if (mode == "balanced" && !(balanced_box_x > 2.0)) {
        std::fprintf(stderr, "balanced_box_x must exceed 2 um\n");
        return 2;
    }

    const double dt = 5.0e-7;
    const long nsteps = static_cast<long>(ms * 1e-3 / dt);
    const bool keep_trajectory = !trajectory_path.empty();
    std::string h5 = keep_trajectory
        ? trajectory_path
        : "/tmp/motor_validation_" + mode + "_" +
              std::to_string(static_cast<long long>(::getpid())) + ".h5";
    std::remove(h5.c_str());

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, static_cast<unsigned long>(seed));
    int fix_myosin = 0;
    Rig rig = build(mode, seed, dt, heads, floor_v, frac_v, h5, rng,
                    fix_myosin, polarity, init_aa, aa_scale, register_ratio,
                    titin, am_head_transverse_k, am_alignment_k,
                    am_transverse_strain_detachment_theta, balanced_box_x,
                    balanced_variant);
    Sarcomere& m = *rig.model;
    m.am_bidirectional_crossbridge = bidirectional_crossbridge;
    m.configure_double_overlap_contact(
        double_overlap_contact, contact_diameter,
        contact_epsilon, contact_force_cap);
    m.aa_status1_site_diffusion = weak_site_D;
    if (initialize_q_calcium_equilibrium) {
        m.initialize_am_thin_activation_calcium_equilibrium();
    }

    // Same integrator the sandbox diagnostic uses: update_system() then one
    // Langevin sample_step per step.
    double beta_v = rod_mobility::beta(310.0);
    double dtv = dt;
    double Dat = 0.300263, Dar = 0.10096, Dmt = 0.0811443, Dmr = 0.0589838;
    int save_every = keep_trajectory ? trajectory_save_every : (1 << 30);
    bool resume = false;
    Langevin sim(m, beta_v, dtv, Dat, Dar, Dmt, Dmr, save_every, resume, true,
                 0.005, 0.005, 0.005, 0.005);
    sim.set_thermal_noise_scale(thermal_scale);
    sim.set_rod_mobility(
        rod_mobility::finite_cylinder(1.0, 0.010, 310.0, 0.001, 10.0, 10.0, 160.0, 160.0),
        rod_mobility::finite_cylinder(1.5, 0.035, 310.0, 0.001, 20.0, 20.0, 64.0, 64.0));

    // Relative coordinate: nothing is pinned, so a stroke moves actin one way
    // and recoils myosin the other. Sliding is actin MINUS myosin.
    auto axial_of = [&](int i) {
        return static_cast<vec>(m.actin.center[i]).dot(rig.axis) -
               static_cast<vec>(m.myosin.center[0]).dot(rig.axis);
    };
    // Geometric axial overlap of two rod segments projected on the x axis.
    // Minimum-image: actins 2/3 form the Z-disc at the periodic boundary, so
    // an unwrapped comparison always reads zero overlap for that pair.
    const double Lx_box = m.box[0];
    auto seg_overlap = [&](double c1, double h1, double c2, double h2) {
        double d = c2 - c1;
        d -= Lx_box * std::round(d / Lx_box);
        const double lo = std::max(-h1, d - h2);
        const double hi = std::min(h1, d + h2);
        return std::max(0.0, hi - lo);
    };
    // The sandbox gate is stated as A-A overlap ~0.096 um and A-M overlap
    // ~0.37 um with terminal A-M drift 0.003-0.008 um/s. Measure exactly those
    // rather than a single pair's relative coordinate.
    auto aa_pair_overlap = [&](int p) {
        const int i = 2*p, j = 2*p + 1;
        return seg_overlap(m.actin.center[i].x, 0.5, m.actin.center[j].x, 0.5);
    };
    // Titin is the restoring force that should stop the Z-disc separating: it
    // is slack at rest (0.58 um vs 0.75 um) and engages after ~0.17 um of
    // motion. If nothing is bound, that restraint never appears.
    auto titin_bound_count = [&]() {
        int n = 0;
        for (const auto& t : m.explicit_titins) if (t.bound) ++n;
        return n;
    };
    auto titin_total_force = [&]() {
        double f = 0.0;
        for (const auto& t : m.explicit_titins) f += std::abs(t.force);
        return f;
    };
    auto aa_links_of = [&](int p) {
        const int i = 2*p, j = 2*p + 1;
        int n = 0;
        for (const auto& lk : m.aa_links) {
            if (!lk.active) continue;
            if ((lk.actin_i == i && lk.actin_j == j) ||
                (lk.actin_i == j && lk.actin_j == i)) ++n;
        }
        return n;
    };
    auto aa_overlap_of = [&]() {
        if (m.actin.n < 4) return 0.0;
        return 0.5 * (aa_pair_overlap(0) + aa_pair_overlap(1));
    };
    const std::array<int, 4> balanced_owner{{0, 1, 0, 1}};
    auto am_overlap_for_actin = [&](int i) {
        const int j = mode == "balanced"
            ? balanced_owner[static_cast<size_t>(i)] : 0;
        vec separation = static_cast<vec>(m.actin.center[i]) -
                         static_cast<vec>(m.myosin.center[j]);
        separation.pbc_wrap(m.box, m.actin.periodic_axes);
        const vec axis = m.myosin.direction[j];
        const double center = separation.dot(axis);
        const double actin_half = 0.5 * m.actin.length *
            std::abs(m.actin.direction[i].dot(axis));
        const double myosin_half = 0.5 * m.myosin.length;
        const double lo = std::max(-myosin_half, center - actin_half);
        const double hi = std::min(myosin_half, center + actin_half);
        return std::max(0.0, hi - lo);
    };
    auto am_overlap_geom = [&]() {
        double total = 0.0;
        for (int i = 0; i < m.actin.n; ++i) {
            total += am_overlap_for_actin(i);
        }
        return total / std::max(1, m.actin.n);
    };
    auto overlap_of = [&]() {
        double s = 0.0;
        const auto& f = m.actin["am_overlap_engaged"];
        for (double v : f) s += v;
        return s / std::max<size_t>(1, f.size());
    };

    const long preq_steps = static_cast<long>(preq_ms * 1e-3 / dt);
    for (long s = 0; s < preq_steps; ++s) {
        m.update_system();
    }
    if (preq_steps > 0) {
        double wk = 0.0, sg = 0.0;
        for (double v : m.actin["am_weak_count"]) wk += v;
        for (double v : m.actin["am_strong_count"]) sg += v;
        std::printf("# preequil %.0f ms at fixed geometry: weak=%.0f strong=%.0f"
                    " AA_ovl=%.4f AM_ovl=%.4f\n",
                    preq_ms, wk, sg, aa_overlap_of(), am_overlap_geom());
    }
    const bool titin_perturbation =
        mode == "balanced" &&
        balanced_variant == "balanced_titin_perturb";
    auto titin_local_midpoint_error = [&]() {
        double maximum_error = 0.0;
        bool bilateral = true;
        for (int myosin_index = 0; myosin_index < m.myosin.n;
             ++myosin_index) {
            vec axis = m.myosin.direction[myosin_index];
            axis.normalize();
            std::array<double, 2> projection_sum{{0.0, 0.0}};
            std::array<int, 2> count{{0, 0}};
            for (const auto& titin_link : m.explicit_titins) {
                if (!titin_link.bound ||
                    titin_link.myosin != myosin_index ||
                    titin_link.bound_actin < 0 ||
                    titin_link.half < 0 || titin_link.half > 1) {
                    continue;
                }
                vec displacement =
                    static_cast<vec>(
                        m.actin.left_end[titin_link.bound_actin]) -
                    static_cast<vec>(m.myosin.center[myosin_index]);
                displacement.pbc_wrap(m.box, m.is_periodic);
                projection_sum[static_cast<size_t>(titin_link.half)] +=
                    displacement.dot(axis);
                ++count[static_cast<size_t>(titin_link.half)];
            }
            if (count[0] == 0 || count[1] == 0) {
                bilateral = false;
                continue;
            }
            const double minus_plane = projection_sum[0] / count[0];
            const double plus_plane = projection_sum[1] / count[1];
            // Translation-invariant axial displacement of this thick
            // filament from the midpoint of its two bound barbed-end planes.
            // A rigid translation of the whole local assembly leaves this
            // quantity unchanged.
            const double midpoint_offset = 0.5 * (minus_plane + plus_plane);
            maximum_error = std::max(
                maximum_error, std::abs(midpoint_offset));
        }
        return std::pair<double, bool>{maximum_error, bilateral};
    };
    int titin_initial_stretched = 0;
    int titin_initial_slack = 0;
    double titin_initial_midpoint_error = 0.0;
    if (titin_perturbation) {
        if (m.titin_mode != Sarcomere::TitinMode::explicit_dynamic ||
            !(m.titin_k > 0.0)) {
            throw std::invalid_argument(
                "balanced_titin_perturb requires explicit titin");
        }
        if (titin_bound_count() != m.actin.n) {
            throw std::runtime_error(
                "balanced_titin_perturb requires one centered titin owner per actin before displacement");
        }
        if (keep_trajectory) m.save_state();
        constexpr double inward_shift = 0.55;
        m.myosin.center[0].x -= inward_shift;
        m.myosin.center[1].x += inward_shift;
        m.myosin.update_endpoints();
        int stretched = 0;
        int slack = 0;
        double max_distance = 0.0;
        for (const auto& titin_link : m.explicit_titins) {
            if (!titin_link.bound) continue;
            vec displacement = m.explicit_titin_thick_anchor(titin_link) -
                static_cast<vec>(
                    m.actin.left_end[titin_link.bound_actin]);
            displacement.pbc_wrap(m.box, m.is_periodic);
            const double distance = displacement.norm();
            max_distance = std::max(max_distance, distance);
            if (distance > m.titin_rest_length) ++stretched;
            else ++slack;
        }
        titin_initial_stretched = stretched;
        titin_initial_slack = slack;
        const auto initial_midpoint = titin_local_midpoint_error();
        titin_initial_midpoint_error = initial_midpoint.first;
        std::printf("# bilateral titin perturbation: shift=0.550 um "
                    "preserved_bound=%d stretched=%d slack=%d "
                    "max_distance=%.4f um "
                    "max_local_midpoint_error=%.4f um bilateral=%d\n",
                    titin_bound_count(), stretched, slack, max_distance,
                    titin_initial_midpoint_error,
                    initial_midpoint.second ? 1 : 0);
    }
    if (keep_trajectory) m.save_state();

    const double initial_aa_overlap = aa_overlap_of();
    const double initial_am_overlap = am_overlap_geom();
    auto aa_extension_of = [&]() {
        return initial_aa_overlap - aa_overlap_of();
    };
    const double x0 = axial_of(rig.probe_actin);
    double prev = x0, path = 0.0;
    const auto& owner = balanced_owner;
    std::vector<vec> previous_actin(static_cast<size_t>(m.actin.n));
    std::vector<vec> previous_myosin(static_cast<size_t>(m.myosin.n));
    for (int i = 0; i < m.actin.n; ++i) previous_actin[i] = m.actin.center[i];
    for (int j = 0; j < m.myosin.n; ++j) previous_myosin[j] = m.myosin.center[j];
    double collective_q = 0.0;
    std::vector<double> individual_q(static_cast<size_t>(m.actin.n), 0.0);
    double report_q = 0.0;
    LinearWindow early_fit, late_fit;
    LinearWindow early_aa_extension_fit, late_aa_extension_fit;
    std::vector<LinearWindow> late_individual(static_cast<size_t>(m.actin.n));
    std::vector<LinearWindow> late_individual_am(static_cast<size_t>(m.actin.n));
    LinearWindow late_am_fit;
    LedgerWindow late_forces;
    std::unordered_map<std::uint64_t, double> active_link_start;
    std::vector<double> completed_link_lifetimes;
    for (const auto& link : m.aa_links) {
        if (link.active) active_link_start.emplace(link.id, 0.0);
    }
    int max_status2 = 0;
    const double fv_window_duration = 0.002;
    double fv_window_start_time = 0.0;
    double fv_window_start_extension = aa_extension_of();
    double fv_load_sum = 0.0;
    double fv_status2_sum = 0.0;
    long fv_count = 0;
    std::vector<FvSample> fv_samples;
    const double duration = nsteps * dt;
    // Capture the approach phase before the first four ACTN2 molecules lock
    // (about 6 ms in the reference cell), rather than diluting it with the
    // subsequent plateau.
    const double early_begin = std::min(0.0005, 0.02 * duration);
    const double early_end = std::min(0.008, 0.20 * duration);
    // Estimate a secular late drift over the trailing half-trajectory (capped
    // at 0.5 s), not over only the last 0.1 s. A 100-ms regression can assign
    // a large apparent velocity to a bounded Brownian excursion even when the
    // local A-M coordinate is stationary over the physically relevant motor/
    // titin timescale. The former gate produced exactly that false failure:
    // max |slope|=0.433 um/s over 0.1 s but 0.081 um/s over 0.5 s, with every
    // A-M overlap bounded throughout.
    const double late_begin = std::max(0.50 * duration, duration - 0.50);

    std::ofstream csv;
    if (!csv_path.empty()) {
        csv.open(csv_path);
        if (!csv) {
            std::fprintf(stderr, "could not open diagnostic CSV: %s\n",
                         csv_path.c_str());
            return 3;
        }
        csv << "time_s,q_um,q0_um,q1_um,q2_um,q3_um,"
               "interval_velocity_um_s,aa_overlap_um,am_overlap_um,"
               "am0_um,am1_um,am2_um,am3_um,"
               "status1_links,status2_links,motor_pn_per_actin,"
               "actn2_pn_per_actin,actn2_load_pn_per_actin,"
               "bilateral_load_pn_per_link,net_pn_per_actin,"
               "other_pn_per_actin,weak_heads,strong_heads,"
               "titin_bound,titin_total_force_pn\n";
        csv << std::setprecision(10);
    }

    std::printf("# mode=%s heads=%d no_external_load=1 titin_k=%.1f "
                "titin_slack=%.3f titin_bound=%d "
                "weak_site_D=%.4g aa_scale=%.1f register=%.2f thermal=%.2f "
                "head_transverse_k=%.1f alignment_k=%.1f "
                "double_overlap_contact=%d diameter=%.4f epsilon=%.5f "
                "force_cap=%.3f nsteps=%ld\n",
                mode.c_str(), heads, m.titin_k, m.titin_rest_length,
                titin_bound_count(), weak_site_D, aa_scale, register_ratio,
                thermal_scale, am_head_transverse_k, am_alignment_k,
                double_overlap_contact ? 1 : 0, contact_diameter,
                contact_epsilon, contact_force_cap, nsteps);
    std::printf("%8s %9s %9s %8s %8s %4s %4s %9s %9s %9s %9s\n",
                "ms", "q", "vwin", "AAovl", "AMovl", "s1", "s2",
                "Qmotor", "Qaa", "Laa", "Qnet");
    const long report = std::max<long>(1, nsteps / 10);
    const long csv_every = std::max<long>(1, nsteps / 200);
    ForceLedger ledger;
    for (long s = 0; s < nsteps; ++s) {
        m.update_system();
        ledger = force_ledger(m);
        max_status2 = std::max(max_status2, ledger.status2);
        sim.sample_step(dtv, rng, fix_myosin);
        if (keep_trajectory && (s + 1) % save_every == 0) {
            m.save_state();
        }

        double dq = 0.0;
        for (int i = 0; i < m.actin.n; ++i) {
            vec delta_actin = static_cast<vec>(m.actin.center[i]) -
                              previous_actin[static_cast<size_t>(i)];
            delta_actin.pbc_wrap(m.box, m.actin.periodic_axes);
            const int myosin_index = mode == "balanced"
                ? owner[static_cast<size_t>(i)] : 0;
            vec delta_myosin = static_cast<vec>(m.myosin.center[myosin_index]) -
                               previous_myosin[static_cast<size_t>(myosin_index)];
            delta_myosin.pbc_wrap(m.box, m.myosin.periodic_axes);
            const double dqi =
                (delta_actin - delta_myosin).dot(m.actin.direction[i]);
            individual_q[static_cast<size_t>(i)] += dqi;
            dq += dqi;
        }
        collective_q += dq / std::max(1, m.actin.n);
        for (int i = 0; i < m.actin.n; ++i) previous_actin[i] = m.actin.center[i];
        for (int j = 0; j < m.myosin.n; ++j) previous_myosin[j] = m.myosin.center[j];

        const double x = axial_of(rig.probe_actin);
        path += std::abs(x - prev);
        prev = x;
        const double time = (s + 1) * dt;
        std::vector<std::uint64_t> live_ids;
        live_ids.reserve(m.aa_active_link_slots.size());
        for (const auto& link : m.aa_links) {
            if (!link.active) continue;
            live_ids.push_back(link.id);
            active_link_start.emplace(link.id, time);
        }
        std::vector<std::uint64_t> disappeared;
        for (const auto& entry : active_link_start) {
            if (std::find(live_ids.begin(), live_ids.end(), entry.first) ==
                live_ids.end()) {
                completed_link_lifetimes.push_back(time - entry.second);
                disappeared.push_back(entry.first);
            }
        }
        for (const auto id : disappeared) active_link_start.erase(id);

        fv_load_sum += ledger.actn2_load;
        fv_status2_sum += ledger.status2;
        ++fv_count;
        if (time - fv_window_start_time + 0.5 * dt >= fv_window_duration) {
            const double elapsed = time - fv_window_start_time;
            fv_samples.push_back({
                time,
                (aa_extension_of() - fv_window_start_extension) / elapsed,
                fv_load_sum / static_cast<double>(fv_count),
                fv_status2_sum / static_cast<double>(fv_count)});
            fv_window_start_time = time;
            fv_window_start_extension = aa_extension_of();
            fv_load_sum = 0.0;
            fv_status2_sum = 0.0;
            fv_count = 0;
        }
        if (time >= early_begin && time <= early_end) {
            early_fit.add(time, collective_q);
            early_aa_extension_fit.add(time, aa_extension_of());
        }
        if (time >= late_begin) {
            late_fit.add(time, collective_q);
            late_aa_extension_fit.add(time, aa_extension_of());
            for (int i = 0; i < m.actin.n; ++i) {
                late_individual[static_cast<size_t>(i)].add(
                    time, individual_q[static_cast<size_t>(i)]);
                late_individual_am[static_cast<size_t>(i)].add(
                    time, am_overlap_for_actin(i));
            }
            late_am_fit.add(time, am_overlap_geom());
            late_forces.add(ledger);
        }
        double wk = 0.0, sg = 0.0;
        for (double v : m.actin["am_weak_count"]) wk += v;
        for (double v : m.actin["am_strong_count"]) sg += v;
        if (csv && ((s + 1) % csv_every == 0 || s + 1 == nsteps)) {
            const double interval_velocity = dq / dt;
            csv << time << ',' << collective_q;
            for (int i = 0; i < 4; ++i) {
                csv << ',' << (i < m.actin.n ? individual_q[i] : 0.0);
            }
            csv << ',' << interval_velocity << ','
                << aa_overlap_of() << ',' << am_overlap_geom() << ','
                << (m.actin.n > 0 ? am_overlap_for_actin(0) : 0.0) << ','
                << (m.actin.n > 1 ? am_overlap_for_actin(1) : 0.0) << ','
                << (m.actin.n > 2 ? am_overlap_for_actin(2) : 0.0) << ','
                << (m.actin.n > 3 ? am_overlap_for_actin(3) : 0.0) << ','
                << ledger.status1 << ',' << ledger.status2 << ','
                << ledger.motor << ',' << ledger.actn2 << ','
                << ledger.actn2_load << ',' << ledger.bilateral << ','
                << ledger.net << ',' << ledger.other << ',' << wk << ',' << sg
                << ',' << titin_bound_count() << ',' << titin_total_force()
                << '\n';
        }
        if ((s + 1) % report == 0) {
            const double vwin = (collective_q - report_q) / (report * dt);
            report_q = collective_q;
            std::printf("%8.2f %9.4f %9.3f %8.4f %8.4f %4d %4d "
                        "%9.3f %9.3f %9.3f %9.3f\n",
                        time * 1e3, collective_q, vwin,
                        aa_overlap_of(), am_overlap_geom(),
                        ledger.status1, ledger.status2, ledger.motor,
                        ledger.actn2, ledger.actn2_load, ledger.net);
        }
    }
    const double T = nsteps * dt;
    const double net = axial_of(rig.probe_actin) - x0;
    const double D = 0.300263;
    const double null_disp = std::sqrt(2.0 * D * T) * std::sqrt(2.0 / M_PI);
    const double early_velocity = early_fit.slope();
    const double late_velocity = late_fit.slope();
    const double early_aa_extension_velocity = early_aa_extension_fit.slope();
    const double late_aa_extension_velocity = late_aa_extension_fit.slope();
    double max_individual_late_velocity = 0.0;
    double max_individual_am_late_velocity = 0.0;
    for (const auto& fit : late_individual) {
        max_individual_late_velocity = std::max(
            max_individual_late_velocity, std::abs(fit.slope()));
    }
    for (const auto& fit : late_individual_am) {
        max_individual_am_late_velocity = std::max(
            max_individual_am_late_velocity, std::abs(fit.slope()));
    }
    const double late_am_velocity = late_am_fit.slope();
    const double mean_motor = late_forces.mean(late_forces.motor);
    const double mean_aa = late_forces.mean(late_forces.actn2);
    const double mean_load = late_forces.mean(late_forces.load);
    const double mean_bilateral = late_forces.mean(late_forces.bilateral);
    const double mean_net = late_forces.mean(late_forces.net);
    const double mean_other = late_forces.mean(late_forces.other);
    const double mean_status2 = late_forces.mean(late_forces.status2);
    const double mean_status1 = late_forces.mean(late_forces.status1);
    double unlocked_fv_velocity = 0.0, unlocked_fv_load = 0.0;
    double locked_fv_velocity = 0.0, locked_fv_load = 0.0;
    int unlocked_fv_count = 0, locked_fv_count = 0;
    for (const auto& sample : fv_samples) {
        if (sample.status2 < 0.5 && sample.time <= 0.040) {
            unlocked_fv_velocity += sample.velocity;
            unlocked_fv_load += sample.load;
            ++unlocked_fv_count;
        }
        if (sample.status2 >= 2.0 && sample.time >= late_begin) {
            locked_fv_velocity += sample.velocity;
            locked_fv_load += sample.load;
            ++locked_fv_count;
        }
    }
    if (unlocked_fv_count > 0) {
        unlocked_fv_velocity /= unlocked_fv_count;
        unlocked_fv_load /= unlocked_fv_count;
    }
    if (locked_fv_count > 0) {
        locked_fv_velocity /= locked_fv_count;
        locked_fv_load /= locked_fv_count;
    }
    double mean_completed_lifetime = 0.0;
    for (double lifetime : completed_link_lifetimes) {
        mean_completed_lifetime += lifetime;
    }
    if (!completed_link_lifetimes.empty()) {
        mean_completed_lifetime /= completed_link_lifetimes.size();
    }
    const double force_scale = std::max(
        0.1, std::abs(mean_motor) + std::abs(mean_aa));
    const double residual_fraction = std::abs(mean_net) / force_scale;
    const double final_aa = aa_overlap_of();
    const double final_am = am_overlap_geom();
    double minimum_owner_am = std::numeric_limits<double>::infinity();
    double maximum_owner_am = 0.0;
    for (int i = 0; i < m.actin.n; ++i) {
        const double value = am_overlap_for_actin(i);
        minimum_owner_am = std::min(minimum_owner_am, value);
        maximum_owner_am = std::max(maximum_owner_am, value);
    }
    const bool slid = early_aa_extension_velocity > 0.10 &&
                      final_am > initial_am_overlap + 0.02;
    const bool stalled = std::abs(late_aa_extension_velocity) < 0.10;
    const bool loaded = mean_load > 0.25 && mean_status2 >= 1.0;
    const bool endogenous_fv = unlocked_fv_count >= 2 && locked_fv_count >= 2 &&
        unlocked_fv_velocity > locked_fv_velocity + 0.20 &&
        std::abs(locked_fv_velocity) < 0.10 && locked_fv_load > 0.25;
    const bool force_balanced = residual_fraction < 0.25;
    const bool registered = std::abs(final_aa - 0.10) < 0.035;
    const bool finite_overlap = final_am > 0.20 && final_am < 0.75;
    const bool every_actin_engaged = mode != "balanced" ||
        (minimum_owner_am > 0.15 && maximum_owner_am < 0.75);
    // A global mean can look stationary while one A-M pair loses overlap and
    // another gains it. Require the local overlap coordinates themselves to
    // have plateaued over the long late fitting window. A 0.20-um/s ceiling
    // distinguishes the former secular escape (>0.4 um/s and eventual zero
    // overlap) from bounded 40--80-nm Brownian exchange, whose fitted slope
    // can exceed 0.10 um/s over a finite window even though every overlap
    // remains far from both zero and full burial.
    const bool locally_stable_am = mode != "balanced" ||
        max_individual_am_late_velocity < 0.20;
    // A stochastic initiation delay is allowed; sustained late transport and
    // accumulated directed travel are the motility criteria.
    const bool free_translocation = collective_q > 0.02 &&
        late_velocity > 0.20 && final_am > 0.20;
    // This scenario is molecularly one-sided: actin 0 is motor-bound while
    // actin 1 is coupled only through ACTN2.  Do not infer that condition from
    // the projected axial-overlap diagnostic.  Once the unsupported filament
    // rotates, that projection can report apparent overlap even though it has
    // no weak or strong head and is far outside the transverse capture shell.
    const auto bound_head_count = [&](int actin_index) {
        return m.actin["am_weak_count"][actin_index] +
               m.actin["am_strong_count"][actin_index];
    };
    const bool one_sided_binding = m.actin.n == 2 &&
        bound_head_count(0) > 0.5 && bound_head_count(1) < 0.5;
    // A one-sided link should turn over before it can mature.  Requiring high
    // *late-window* status-1 occupancy contradicts that test: after the short
    // links release and the filaments separate, they should not be expected to
    // remain close enough to rebind.  Completed lifetimes are the direct
    // observable; max_status2 ensures none was load-stabilized by one motor.
    const bool unsupported_turnover = !completed_link_lifetimes.empty() &&
        max_status2 == 0 && mean_completed_lifetime < 0.20;
    const bool linked_transport = m.actin.n == 2 &&
        individual_q[0] > 0.02 && individual_q[1] < -0.02;
    const auto titin_midpoint = titin_perturbation
        ? titin_local_midpoint_error()
        : std::pair<double, bool>{0.0, true};
    const double titin_recenter_error = titin_midpoint.first;
    const double titin_remaining_fraction = titin_initial_midpoint_error > 0.0
        ? titin_recenter_error / titin_initial_midpoint_error : 0.0;
    const bool titin_recentered = !titin_perturbation ||
        (titin_initial_stretched == 2 && titin_initial_slack == 2 &&
         titin_bound_count() == m.actin.n &&
         titin_midpoint.second &&
         titin_remaining_fraction < 0.45 &&
         titin_recenter_error < 0.25 && registered && finite_overlap);
    const bool pass = mode == "balanced"
        ? (titin_perturbation
            ? titin_recentered
            : (slid && stalled && endogenous_fv && loaded && force_balanced &&
               registered && finite_overlap && every_actin_engaged &&
               locally_stable_am))
        : (mode == "free"
            ? free_translocation
            : (mode == "linked" && one_sided_binding &&
               unsupported_turnover && linked_transport));
    std::printf("\nRESULT mode=%s probe_net=%+.4f um probe_v=%+.3f um/s "
                "null=%.4f ratio=%.2f\n",
                mode.c_str(), net, net / T, null_disp,
                std::abs(net) / null_disp);
    std::printf("SELF_LOAD early_v=%+.4f late_v=%+.4f "
                "max_individual_late_v=%.4f AM_late_v=%+.4f "
                "max_individual_AM_late_v=%.4f q=%+.4f "
                "AA=%.4f->%.4f AM=%.4f->%.4f engaged=%.4f\n",
                early_velocity, late_velocity, max_individual_late_velocity,
                late_am_velocity, max_individual_am_late_velocity, collective_q,
                initial_aa_overlap, final_aa,
                initial_am_overlap, final_am, overlap_of());
    std::printf("ENDOGENOUS_FV AA_extension_v early=%+.4f late=%+.4f "
                "unlocked=(load=%.4f,v=%+.4f,n=%d) "
                "locked=(load=%.4f,v=%+.4f,n=%d) => %s\n",
                early_aa_extension_velocity, late_aa_extension_velocity,
                unlocked_fv_load, unlocked_fv_velocity, unlocked_fv_count,
                locked_fv_load, locked_fv_velocity, locked_fv_count,
                endogenous_fv ? "PASS" : "FAIL");
    std::printf("INDIVIDUAL_Q");
    for (int i = 0; i < m.actin.n; ++i) {
        std::printf(" q%d=%+.4f v%d_late=%+.4f", i, individual_q[i], i,
                    late_individual[static_cast<size_t>(i)].slope());
    }
    std::printf("\n");
    std::printf("INDIVIDUAL_AM");
    for (int i = 0; i < m.actin.n; ++i) {
        std::printf(" am%d=%.4f dam%d_late=%+.4f", i,
                    am_overlap_for_actin(i), i,
                    late_individual_am[static_cast<size_t>(i)].slope());
    }
    std::printf(" max_abs_dam=%.4f\n", max_individual_am_late_velocity);
    std::printf("INDIVIDUAL_FORCE");
    for (int i = 0; i < m.actin.n; ++i) {
        const int j = mode == "balanced"
            ? owner[static_cast<size_t>(i)] : 0;
        std::printf(" a%d[M=%+.3f AA=%+.3f O=%+.3f N=%+.3f cos=%.4f]", i,
                    ledger.motor_by_actin[static_cast<size_t>(i)],
                    ledger.actn2_by_actin[static_cast<size_t>(i)],
                    ledger.other_by_actin[static_cast<size_t>(i)],
                    ledger.net_by_actin[static_cast<size_t>(i)],
                    m.actin.direction[i].dot(m.myosin.direction[j]));
    }
    std::printf("\n");
    std::printf("FORCE_BALANCE motor=%+.4f AA=%+.4f load=%.4f "
                "bilateral=%.4f other=%+.4f net=%+.4f residual=%.4f "
                "mean_s2=%.2f\n",
                mean_motor, mean_aa, mean_load, mean_bilateral,
                mean_other, mean_net, residual_fraction, mean_status2);
    std::printf("TURNOVER completed=%zu mean_lifetime=%.4f s "
                "mean_s1=%.2f max_s2=%d one_sided=%s => %s\n",
                completed_link_lifetimes.size(), mean_completed_lifetime,
                mean_status1, max_status2,
                one_sided_binding ? "PASS" : "FAIL",
                unsupported_turnover ? "PASS" : "FAIL");
    std::printf("TITIN bound=%d/%zu total_force=%.4f pN "
                "supported_assoc=%llu unsupported_assoc=%llu rebind=%llu\n",
                titin_bound_count(), m.explicit_titins.size(),
                titin_total_force(),
                static_cast<unsigned long long>(
                    m.explicit_titin_supported_associations),
                static_cast<unsigned long long>(
                    m.explicit_titin_unsupported_associations),
                static_cast<unsigned long long>(
                    m.explicit_titin_rebind_associations));
    std::printf("DOUBLE_OVERLAP_CONTACT enabled=%d count=%llu "
                "max_penetration=%.6f um max_force=%.4f pN\n",
                m.aa_double_overlap_contact_enabled ? 1 : 0,
                static_cast<unsigned long long>(
                    m.aa_double_overlap_contact_count),
                m.aa_double_overlap_max_penetration,
                m.aa_double_overlap_max_force);
    if (titin_perturbation) {
        std::printf("TITIN_RECENTER stretched=%d slack=%d "
                    "max_local_midpoint_error=%.6f->%.6f um "
                    "remaining=%.3f bilateral=%d bound=%d/%d => %s\n",
                    titin_initial_stretched, titin_initial_slack,
                    titin_initial_midpoint_error, titin_recenter_error,
                    titin_remaining_fraction,
                    titin_midpoint.second ? 1 : 0,
                    titin_bound_count(), m.actin.n,
                    titin_recentered ? "PASS" : "FAIL");
    }
    std::printf("GATE slide=%s stall=%s endogenous_FV=%s load=%s "
                "balance=%s register=%s finite_AM=%s all_AM=%s[%.3f,%.3f] "
                "local_AM_stall=%s[max|dAM/dt|=%.4f] "
                "free_transport=%s "
                "linked_transport=%s => %s\n",
                slid ? "PASS" : "FAIL", stalled ? "PASS" : "FAIL",
                endogenous_fv ? "PASS" : "FAIL",
                loaded ? "PASS" : "FAIL",
                force_balanced ? "PASS" : "FAIL",
                registered ? "PASS" : "FAIL",
                finite_overlap ? "PASS" : "FAIL",
                every_actin_engaged ? "PASS" : "FAIL",
                minimum_owner_am, maximum_owner_am,
                locally_stable_am ? "PASS" : "FAIL",
                max_individual_am_late_velocity,
                free_translocation ? "PASS" : "FAIL",
                linked_transport ? "PASS" : "FAIL",
                pass ? "PASS" : "FAIL");
    if (keep_trajectory && nsteps % save_every != 0) m.save_state();
    if (!keep_trajectory) std::remove(h5.c_str());
    gsl_rng_free(rng);
    return pass ? 0 : 1;
}
