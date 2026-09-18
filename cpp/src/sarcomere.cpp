#include "sarcomere.h"
#include "aa_refractory.h"
#include "aa_association_kinetics.h"
#include "catch_bond_kinetics.h"
#include "counter_rng.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <thread>
#include <tuple>
#include <unordered_set>
#include <gsl/gsl_randist.h>

#include "sarcomere_internal_gates.h"


// Parameterized Constructor
Sarcomere::Sarcomere(int& n_actins, int& n_myosins, vector box0, double& actin_length, double& myosin_length,
        double& myosin_radius, double& am_cutoff, double& am_optimal, double& aa_cutoff, double& aa_optimal,
        double& k_on, double& k_off,
        double& base_lifetime, double& directional_base_lifetime, double& lifetime_coeff,
        double& diff_coeff_ratio, double& k_aa, double& kappa_aa, double& k_am, double& kappa_am, double& k_mm, double& v_am,
        std::string& filename, gsl_rng* rng, int& seed, int& fix_myosin, double& dt, double tau_rec,
        double titin_k, double titin_rest_length, bool& directional, int max_myosin_bonds,
        double max_actin_force_param, double max_myosin_force_param,
        double max_actin_torque_param, double max_myosin_torque_param)
            : actin(n_actins, actin_length, box0, rng),
              myosin(n_myosins, myosin_length, myosin_radius, box0, rng),
              myosinIndicesPerActin(n_actins),
              actinIndicesPerMyosin(n_myosins),
              neighbor_list(0.0, box0, 0.0, is_periodic),
                actin_actin_bonds(n_actins, std::vector<BondFlag>(n_actins, 0)),
                actin_actin_status(n_actins, std::vector<BondFlag>(n_actins, 0)),
                actin_actin_status_prev(n_actins, std::vector<BondFlag>(n_actins, 0)),
                actin_actin_lifetime(n_actins, std::vector<int>(n_actins, 0)),
                actin_recovery_until(n_actins, std::vector<RecoveryStep>(n_actins, 0)),
                aa_status2_memory_until(n_actins, std::vector<RecoveryStep>(n_actins, 0)),
                am_bonds(n_actins, std::vector<BondFlag>(n_myosins, 0)),
                am_thin_activation(n_actins, 0.0),
                am_half_tension(n_myosins, {0.0, 0.0}),
                myosin_material_e1(n_myosins),
                myosin_material_e2(n_myosins),
                am_state(n_actins, std::vector<BondFlag>(n_myosins, 0)),
                am_seen(n_actins, std::vector<BondFlag>(n_myosins, 0)),
                am_target_axial(n_actins, std::vector<double>(n_myosins, 0.0)),
                am_heads_by_myosin(n_myosins),
                actin_forces_temp(omp_get_max_threads(), std::vector<vec>(n_actins, {0, 0, 0})),
                myosin_forces_temp(omp_get_max_threads(), std::vector<vec>(n_myosins, {0, 0, 0})),
                myosin_velocities_temp(omp_get_max_threads(), std::vector<vec>(n_myosins, {0, 0, 0})),
                actin_torques_temp(omp_get_max_threads(), std::vector<vec>(n_actins, {0, 0, 0})),
                myosin_torques_temp(omp_get_max_threads(), std::vector<vec>(n_myosins, {0, 0, 0})),
                actin_cb_status_temp(omp_get_max_threads(), std::vector<int>(n_actins, 0)),
                myosin_f_load_temp(omp_get_max_threads(), std::vector<std::array<double, 2>>(n_myosins, {0.0, 0.0})),
                cb_breakage_events_temp(omp_get_max_threads()),
                aa_completed_lifetimes_temp(omp_get_max_threads()),
                am_transition_events_temp(omp_get_max_threads()),
                aa_state_events_temp(omp_get_max_threads()),
                actin_chem_entropy_delta_temp(omp_get_max_threads(), std::vector<double>(n_actins, 0.0)),
                actin_chem_input_delta_temp(omp_get_max_threads(), std::vector<double>(n_actins, 0.0)),
                actin_chem_binding_delta_temp(omp_get_max_threads(), std::vector<double>(n_actins, 0.0)),
                actin_chem_unbinding_delta_temp(omp_get_max_threads(), std::vector<double>(n_actins, 0.0)),
                myosin_f_load(n_myosins, {0.0, 0.0}),
                actinIndicesPerMyosin_temp(omp_get_max_threads(), utils::MoleculeConnection(n_myosins)),
                // KMC randomness is counter-based and object-owned. Retain
                // one legacy GSL stream solely for backward-compatible
                // checkpoint fields, independent of the OpenMP team size.
                rng_engines(1, nullptr),
                active_aa_pairs_temp(omp_get_max_threads()),
                active_am_pairs_temp(omp_get_max_threads()),
                recovery_pairs_temp(omp_get_max_threads()),
                actin_f_load_computed(n_actins, false),
                actin_f_load_mutex(n_actins),
                max_actin_force(max_actin_force_param),
                max_myosin_force(max_myosin_force_param),
                max_actin_torque(max_actin_torque_param),
                max_myosin_torque(max_myosin_torque_param)

            {
            if (!(max_actin_force > 0.0)) {
                max_actin_force = std::numeric_limits<double>::infinity();
            }
            if (!(max_myosin_force > 0.0)) {
                max_myosin_force = std::numeric_limits<double>::infinity();
            }
            if (!(max_actin_torque > 0.0)) {
                max_actin_torque = std::numeric_limits<double>::infinity();
            }
            if (!(max_myosin_torque > 0.0)) {
                max_myosin_torque = std::numeric_limits<double>::infinity();
            }
            box.resize(3);
            box[0] = box0[0];
            box[1] = box0[1];
            box[2] = box0[2];
            this->k_aa = k_aa;
            this->kappa_aa = kappa_aa;
            this->k_on = k_on;
            this->k_off = k_off;
            this->k_am = k_am;
            this->kappa_am = kappa_am;
            this->v_am = v_am;
            this->k_mm = k_mm;
            this->skin_distance = skin_distance;
            this->filename = filename;
            this->rng = rng;
            this->initial_seed = seed;
            this->fix_myosin = fix_myosin;
            this->dt = dt;
            this->am_cutoff = am_cutoff;
            this->am_optimal = am_optimal;
            this->aa_cutoff = aa_cutoff;
            this->aa_optimal = aa_optimal;
            this->base_lifetime = base_lifetime;
            this->directional_base_lifetime = directional_base_lifetime;
            this->lifetime_coeff = lifetime_coeff;
            this->diff_coeff_ratio = diff_coeff_ratio;
            this->titin_k = titin_k;
            this->titin_rest_length = titin_rest_length;
            this->directional = directional;
            this->max_myosin_bonds = max_myosin_bonds;
            this->max_strong_actin_bonds = 2;
            if (tau_rec > 0) {
                bond_recovery_steps = static_cast<size_t>(std::ceil(tau_rec / dt));
            } else {
                bond_recovery_steps = 0;
            }
            double neighbor_cutoff = std::max(am_cutoff, aa_cutoff);
            if (k_mm > 0.0) {
                neighbor_cutoff = std::max(
                    neighbor_cutoff,
                    std::max(myomesin_axial_capture,
                             myomesin_radial_rest + myomesin_radial_capture));
            }
            cutoff_radius = std::max(actin_length, myosin_length) + neighbor_cutoff;
            double skin_distance = 0.15 * cutoff_radius;
            neighbor_list = NeighborList(cutoff_radius + skin_distance, box, skin_distance / 2, is_periodic);
            neighbor_list.initialize(actin.center_x, actin.center_y, actin.center_z,
                                    myosin.center_x, myosin.center_y, myosin.center_z);
            actin_actin_bonds_prev = actin_actin_bonds;
            actin_actin_lifetime_prev = actin_actin_lifetime;
            // Initialize attach-step tracker for completed lifetime measurements
            aa_attach_step.assign(n_actins, std::vector<int>(n_actins, -1));
            am_bonds_prev = am_bonds;
            actin_basic_tension.resize(n_actins);
            actin_crosslink_ratio.resize(n_actins);
            actin_crosslink_start.resize(n_actins);
            actin_crosslink_end.resize(n_actins);
            actin_n_bonds.resize(n_actins);
            actin_strong_cb_count.resize(n_actins);
            previous_actin_cb_status.assign(n_actins, 0);
            aa_reserved_degree.assign(n_actins, 0);
            aa_links_by_actin.resize(static_cast<size_t>(n_actins));
            n_myosins_per_actin.resize(n_actins);
            am_interaction.resize(n_actins);
            for (int i = 0; i < n_actins; i++) {
                am_interaction[i].resize(n_myosins);
            }
            actin.register_feature("f_load_cb");
            actin_f_load_cb = &actin["f_load_cb"];
            actin.register_feature("cumulative_chem_entropy");
            actin.register_feature("cumulative_chem_input");
            actin.register_feature("cumulative_chem_binding_events");
            actin.register_feature("cumulative_chem_unbinding_events");
            actin.register_feature("myosin_binding_ratio");
            actin.register_feature("crosslink_ratio");
            actin.register_feature("partial_binding_ratio");
            actin.register_feature("am_weak_count");
            actin.register_feature("am_strong_count");
            actin.register_feature("am_motor_force");
            actin.register_feature("am_active_force_x");
            actin.register_feature("am_active_force_y");
            actin.register_feature("am_active_force_z");
            actin.register_feature("am_tether_force");
            // Separate from am_tether_force: until 2026-09-17 BOTH the per-head
            // transverse spring (am_head_transverse_k) and the aggregate pair
            // tether (am_tether_k) summed into am_tether_force, so the recorded
            // value could not be attributed to either and the two radial
            // channels could not be told apart from a run.
            actin.register_feature("am_head_transverse_force");
            actin.register_feature("am_activation");
            actin.register_feature("am_calcium");
            actin.register_feature("am_q");
            actin.register_feature("turnover_age");
            actin.register_feature("cumulative_turnover_events");
            // Overlap is stored as a length. The two weighted orientation sums
            // can be divided by am_overlap_available during analysis.
            actin.register_feature("am_overlap_available");
            actin.register_feature("am_overlap_engaged");
            actin.register_feature("am_overlap_axis_weighted");
            actin.register_feature("am_overlap_productive_weighted");
            actin.register_feature("aa_arm_load_max");
            actin.register_feature("aa_bilateral_load_max");
            actin.register_feature("aa_raw_arm_load_max");
            actin.register_feature("aa_local_link_count");
            actin.register_feature("aa_local_status1_count");
            actin.register_feature("aa_local_status2_count");
            actin.register_feature("aa_local_tension_max");
            actin.register_feature("aa_zdisc_tip_force");
            actin.register_feature("aa_zdisc_tip_contact_count");
            actin.register_feature("cumulative_aa_zdisc_tip_projection_events");
            actin.register_feature("cumulative_aa_zdisc_tip_projection_correction");
            actin.register_feature("titin_bound");
            actin.register_feature("titin_force");
            actin.register_feature("cumulative_force_cap_events");
            actin.register_feature("cumulative_torque_cap_events");
            actin.register_feature("cumulative_rotation_cap_events");
            myosin.register_feature("am_off_count_minus");
            myosin.register_feature("am_off_count_plus");
            myosin.register_feature("am_on_count_minus");
            myosin.register_feature("am_on_count_plus");
            myosin.register_feature("am_tension_minus");
            myosin.register_feature("am_tension_plus");
            myosin.register_feature("titin_bound_minus");
            myosin.register_feature("titin_bound_plus");
            myosin.register_feature("titin_force_minus");
            myosin.register_feature("titin_force_plus");
            for (const char* suffix : {"minus", "plus"}) {
                myosin.register_feature(std::string("am_align_drive_") + suffix);
                myosin.register_feature(std::string("am_align_abs_drive_") + suffix);
                myosin.register_feature(std::string("am_align_active_drive_") + suffix);
                myosin.register_feature(std::string("am_align_tether_drive_") + suffix);
                myosin.register_feature(std::string("am_align_strong_drive_") + suffix);
                myosin.register_feature(std::string("am_align_head_count_") + suffix);
                myosin.register_feature(std::string("am_align_positive_count_") + suffix);
                myosin.register_feature(std::string("am_align_negative_count_") + suffix);
            }
            myosin.register_feature("material_e1_x");
            myosin.register_feature("material_e1_y");
            myosin.register_feature("material_e1_z");
            myosin.register_feature("material_e2_x");
            myosin.register_feature("material_e2_y");
            myosin.register_feature("material_e2_z");
            myosin.register_feature("cumulative_force_cap_events");
            myosin.register_feature("cumulative_torque_cap_events");
            myosin.register_feature("cumulative_rotation_cap_events");
            myosin.register_feature("cumulative_hard_contact_events");
            myosin.register_feature("cumulative_hard_contact_correction");
            myosin.register_feature("cumulative_aa_zdisc_tip_projection_events");
            myosin.register_feature("cumulative_aa_zdisc_tip_projection_correction");
            for (const char* suffix : {"minus", "plus"}) {
                myosin.register_feature(
                    std::string("aa_zdisc_tip_force_") + suffix);
                myosin.register_feature(
                    std::string("aa_zdisc_tip_gap_") + suffix);
                myosin.register_feature(
                    std::string("aa_zdisc_tip_rho_") + suffix);
                myosin.register_feature(
                    std::string("aa_zdisc_tip_gate_") + suffix);
            }
            for (int m = 0; m < n_myosins; ++m) {
                vec axis = myosin.direction[m];
                axis.normalize();
                const vec reference = std::abs(axis.z) < 0.9
                    ? vec{0.0, 0.0, 1.0} : vec{0.0, 1.0, 0.0};
                myosin_material_e1[m] = axis.cross(reference);
                if (myosin_material_e1[m].norm() <= 1e-12) {
                    myosin_material_e1[m] = {1.0, 0.0, 0.0};
                }
                myosin_material_e1[m].normalize();
                myosin_material_e2[m] = axis.cross(myosin_material_e1[m]);
                myosin_material_e2[m].normalize();
                transport_myosin_material_frame(m, axis, axis, 0.0);
            }
            if (actin_f_load_cb) {
                actin_f_load_cb->assign(n_actins, 0.0);
            }
            // Initialize the canonical legacy checkpoint stream.
            for (int t = 0; t < static_cast<int>(rng_engines.size()); ++t) {
                rng_engines[t] = gsl_rng_alloc(gsl_rng_mt19937); 
                gsl_rng_set(rng_engines[t], seed + t);           
            }
            set_periodicity(is_periodic);
            myosin_bond_matrix.assign(n_myosins, std::vector<BondFlag>(n_myosins, 0));
            has_myosin_bond_pairs = false;
            actin_last_turnover_step.assign(static_cast<size_t>(n_actins), 0);
        }

// Destructor
Sarcomere::~Sarcomere() {}

void Sarcomere::reset_release_diagnostics() {
    am_transition_events.clear();
    am_forced_rupture_events.clear();
    am_forced_rupture_events_v2.clear();
    aa_local_link_events.clear();
    aa_local_link_events_v2.clear();
    aa_local_link_events_v3.clear();
    aa_tip_accessibility_events.clear();
    aa_motor_pair_events.clear();
    aa_arm_events.clear();
    explicit_titin_events.clear();
    myomesin_events.clear();
    for (auto& values : am_transition_events_temp) values.clear();

    const auto reset_cap_counters = [](Filament& filament) {
        for (const char* name : {
                 "cumulative_force_cap_events",
                 "cumulative_torque_cap_events",
                 "cumulative_rotation_cap_events",
                 "cumulative_aa_zdisc_tip_projection_events",
                 "cumulative_aa_zdisc_tip_projection_correction"}) {
            auto found = filament.custom_features.find(name);
            if (found != filament.custom_features.end()) {
                std::fill(found->second.begin(), found->second.end(), 0.0);
            }
        }
    };
    reset_cap_counters(actin);
    reset_cap_counters(myosin);
}

void Sarcomere::configure_explicit_titin(
    int links_per_half, double capture_distance, double association_rate,
    double dissociation_rate, bool strain_gated, double transition_theta,
    double rebind_memory, double rebind_rate_multiplier,
    bool require_live_actn2, double k_off_unsupported,
    bool rebind_requires_actn2) {
    if (links_per_half <= 0 || !(capture_distance > 0.0) ||
        association_rate < 0.0 || dissociation_rate < 0.0 ||
        transition_theta < 0.0 || transition_theta > 1.0 ||
        rebind_memory < 0.0 || !std::isfinite(rebind_memory) ||
        rebind_rate_multiplier < 0.0 ||
        !std::isfinite(rebind_rate_multiplier) ||
        !std::isfinite(k_off_unsupported)) {
        throw std::invalid_argument("invalid explicit titin configuration");
    }
    explicit_titin_per_half = links_per_half;
    explicit_titin_capture_distance = capture_distance;
    explicit_titin_k_on = association_rate;
    explicit_titin_k_off = dissociation_rate;
    explicit_titin_strain_gated_kinetics = strain_gated;
    explicit_titin_transition_theta = transition_theta;
    explicit_titin_rebind_memory = rebind_memory;
    explicit_titin_rebind_rate_multiplier = rebind_rate_multiplier;
    explicit_titin_require_live_actn2 = require_live_actn2;
    explicit_titin_k_off_unsupported = k_off_unsupported;
    explicit_titin_rebind_requires_actn2 = rebind_requires_actn2;
    initialize_explicit_titins();
}

void Sarcomere::initialize_explicit_titins() {
    explicit_titins.clear();
    explicit_titin_owner_by_actin.assign(static_cast<size_t>(actin.n), -1);
    if (titin_mode != TitinMode::explicit_dynamic) return;
    if (explicit_titin_per_half <= 0 || myosin.n <= 0) {
        throw std::invalid_argument(
            "explicit titin requires positive links per myosin half");
    }
    explicit_titins.reserve(static_cast<size_t>(
        2 * myosin.n * explicit_titin_per_half));
    const int count = explicit_titin_per_half;
    for (int m = 0; m < myosin.n; ++m) {
        for (int half = 0; half < 2; ++half) {
            const double phase = count > 1 && half == 0
                ? M_PI / static_cast<double>(count) : 0.0;
            for (int track = 0; track < count; ++track) {
                ExplicitTitin titin;
                titin.myosin = m;
                titin.half = half;
                titin.track = track;
                titin.anchor_azimuth = phase +
                    2.0 * M_PI * static_cast<double>(track) /
                        static_cast<double>(count);
                explicit_titins.push_back(titin);
            }
        }
    }
}

void Sarcomere::configure_myomesin_bridges(
    int sites_per_filament, double association_rate,
    double dissociation_rate, double axial_capture, double radial_rest,
    double radial_capture, double radial_stiffness, double alignment_cos,
    double bell_distance, double force_max) {
    if (sites_per_filament < 1 || association_rate < 0.0 ||
        dissociation_rate < 0.0 || !(axial_capture > 0.0) ||
        !(radial_rest > 0.0) || radial_capture < 0.0 ||
        radial_stiffness < 0.0 || alignment_cos < 0.0 ||
        alignment_cos > 1.0 || bell_distance < 0.0 || force_max < 0.0 ||
        !std::isfinite(association_rate) ||
        !std::isfinite(dissociation_rate) ||
        !std::isfinite(axial_capture) || !std::isfinite(radial_rest) ||
        !std::isfinite(radial_capture) ||
        !std::isfinite(radial_stiffness) ||
        !std::isfinite(alignment_cos) || !std::isfinite(bell_distance) ||
        !std::isfinite(force_max)) {
        throw std::invalid_argument("invalid myomesin bridge configuration");
    }
    myomesin_sites_per_filament = sites_per_filament;
    myomesin_k_on = association_rate;
    myomesin_k_off = dissociation_rate;
    myomesin_axial_capture = axial_capture;
    myomesin_radial_rest = radial_rest;
    myomesin_radial_capture = radial_capture;
    myomesin_radial_k = radial_stiffness;
    myomesin_alignment_cos = alignment_cos;
    myomesin_bell_distance = bell_distance;
    myomesin_force_max = force_max;
    myomesin_bridges.clear();
    myomesin_events.clear();
    myomesin_next_id = 1;
}

int Sarcomere::myomesin_degree(int myosin_index) const {
    if (myosin_index < 0 || myosin_index >= myosin.n) return 0;
    int degree = 0;
    for (const auto& bridge : myomesin_bridges) {
        if (bridge.myosin_i == myosin_index ||
            bridge.myosin_j == myosin_index) {
            ++degree;
        }
    }
    return degree;
}

bool Sarcomere::prebind_myomesin_bridge(int first, int second) {
    if (k_mm <= 0.0 || first < 0 || second < 0 || first >= myosin.n ||
        second >= myosin.n || first == second) {
        return false;
    }
    if (first > second) std::swap(first, second);
    for (const auto& bridge : myomesin_bridges) {
        if (bridge.myosin_i == first && bridge.myosin_j == second) {
            return false;
        }
    }
    if (myomesin_degree(first) >= myomesin_sites_per_filament ||
        myomesin_degree(second) >= myomesin_sites_per_filament) {
        return false;
    }
    vec axis, radial;
    double axial = 0.0;
    double radial_distance = 0.0;
    if (!_myomesin_pair_geometry(first, second, axis, axial, radial,
                                 radial_distance)) {
        return false;
    }
    MyomesinBridge bridge;
    bridge.id = myomesin_next_id++;
    bridge.myosin_i = first;
    bridge.myosin_j = second;
    bridge.attach_step = static_cast<RecoveryStep>(current_step);
    bridge.axial_offset = axial;
    bridge.radial_distance = radial_distance;
    myomesin_bridges.push_back(bridge);
    return true;
}

vec Sarcomere::explicit_titin_thick_anchor(
    const ExplicitTitin& titin) const {
    if (titin.myosin < 0 || titin.myosin >= myosin.n ||
        (titin.half != 0 && titin.half != 1)) {
        throw std::invalid_argument("invalid explicit titin anchor");
    }
    vec axis = myosin.direction[titin.myosin];
    axis.normalize();
    vec e1 = myosin_material_e1[static_cast<size_t>(titin.myosin)];
    e1 -= axis * e1.dot(axis);
    if (e1.norm() <= 1e-12) {
        const vec reference = std::abs(axis.z) < 0.9
            ? vec{0.0, 0.0, 1.0} : vec{0.0, 1.0, 0.0};
        e1 = axis.cross(reference);
    }
    e1.normalize();
    vec e2 = axis.cross(e1);
    e2.normalize();
    const double half_sign = titin.half == 0 ? -1.0 : 1.0;
    const vec radial = std::cos(titin.anchor_azimuth) * e1 +
                       std::sin(titin.anchor_azimuth) * e2;
    return static_cast<vec>(myosin.center[titin.myosin]) +
           half_sign * 0.5 * myosin.length * axis +
           myosin.radius * radial;
}

const char* Sarcomere::profile_section_name(int section) {
    switch (section) {
        case kSecSetToZero:         return "set_to_zero";
        case kSecAmBinding:         return "am_binding";
        case kSecThinActivation:    return "thin_activation_and_heads";
        case kSecBondLimit:         return "myosin_bond_limit";
        case kSecConcatConnections: return "concat_connections";
        case kSecAaChemistry:       return "actn2_chemistry";
        case kSecChemReduce:        return "chem_reduce";
        case kSecBondMatrix:        return "myosin_bond_matrix";
        case kSecAmForce:           return "am_force_velocity";
        case kSecAaStateTitin:      return "aa_state_titin_myomesin";
        case kSecVolumeExclusion:   return "volume_exclusion";
        case kSecReduce:            return "force_reduction";
        case kSecTail:              return "titin_bundling_forcecap";
        default:                    return "unknown";
    }
}

void Sarcomere::update_system() {
    using ProfileClock = std::chrono::steady_clock;
    const auto system_start = ProfileClock::now();
    // Advance global step counter each time the system is updated
    current_step++;
    if (titin_k > 0.0 && titin_rest_length > 0.0 &&
        titin_mode == TitinMode::balanced_owning &&
        (actin.n != 4 || myosin.n != 2)) {
        throw std::invalid_argument(
            "balanced_owning titin requires exactly four actins and two myosins");
    }
    am_activation = am_motor_kinetics::activation(
        static_cast<double>(current_step) * dt, am_cycle_period,
        am_activation_rise, am_activation_decay, am_activation_floor,
        am_activation_active_fraction);
    _process_actin_turnover();
    const double bundle_strength = _current_bundle_strength();
    auto& cumulative_chem_entropy = actin["cumulative_chem_entropy"];
    auto& cumulative_chem_input = actin["cumulative_chem_input"];
    auto& cumulative_chem_binding_events = actin["cumulative_chem_binding_events"];
    auto& cumulative_chem_unbinding_events = actin["cumulative_chem_unbinding_events"];
    auto& actin_force_cap_events = actin["cumulative_force_cap_events"];
    auto& actin_torque_cap_events = actin["cumulative_torque_cap_events"];
    auto& myosin_force_cap_events = myosin["cumulative_force_cap_events"];
    auto& myosin_torque_cap_events = myosin["cumulative_torque_cap_events"];
    const auto neighbor_start = ProfileClock::now();
    _update_neighbors();
    system_neighbor_seconds += std::chrono::duration<double>(
        ProfileClock::now() - neighbor_start).count();
    const auto parallel_start = ProfileClock::now();
    if (parallel_section_seconds.size() <
        static_cast<size_t>(omp_get_max_threads())) {
        parallel_section_seconds.resize(omp_get_max_threads(), {});
    if (aa_phase_seconds.size() <
            static_cast<size_t>(omp_get_max_threads())) {
        aa_phase_seconds.resize(omp_get_max_threads(), {});
    }
    }
    #pragma omp parallel
    {
        // Per-thread section stamps.  No barriers are added: `mark` simply
        // advances at each section boundary, so a thread's entry covers its
        // own work plus its wait at that section's implicit barrier.
        auto& section_seconds =
            parallel_section_seconds[omp_get_thread_num()];
        auto section_mark = ProfileClock::now();
        const auto stamp_section = [&](int section) {
            const auto now = ProfileClock::now();
            section_seconds[section] +=
                std::chrono::duration<double>(now - section_mark).count();
            section_mark = now;
        };
        _set_to_zero();
        stamp_section(kSecSetToZero);
        //Step 2: Compute actin-myosin binding
        #pragma omp for schedule(runtime)
        for (int i = 0; i < actin.n; i++) {
            _process_actin_myosin_binding(i);
        }
        stamp_section(kSecAmBinding);
        if (dynamic_am) {
          // _update_thin_filament_activation() now contains its own orphaned
          // worksharing constructs (a `single` for the resize guard and an
          // `omp for` over actins), so every thread in the enclosing team
          // participates instead of parking at a barrier.
          _update_thin_filament_activation();
          _process_dynamic_am_heads();
        }
        stamp_section(kSecThinActivation);
        #pragma omp single
        {
          _enforce_myosin_bond_limit();
        }
        stamp_section(kSecBondLimit);

        // Step 3: Concatenate actinIndicesPerMyosin connections
        #pragma omp for
        for (int i = 0; i < myosin.n; ++i) {
            for (int t = 0; t < omp_get_num_threads(); ++t) {
                const auto& indices =
                    actinIndicesPerMyosin_temp[t].getConnections(i);
                for (int j = 0; j < indices.size(); j++) {
                    actinIndicesPerMyosin.addConnection(i, indices[j]);
                }
            }
        }
        stamp_section(kSecConcatConnections);

        if (aa_local_multilink) {
        // _update_local_aa_links() now contains its own orphaned worksharing
        // constructs (`single` blocks around the shared mutations, `omp for`
        // over actins and links, and an `ordered` region for link
        // creation/destruction), so every thread in the enclosing team
        // participates instead of parking at a barrier.
        _update_local_aa_links();
        } else if (base_lifetime > 0 || directional_base_lifetime > 0 || lifetime_coeff > 0) {
        // Step 4: Compute catch bonds
        #pragma omp for schedule(runtime)
        for (int i = 0; i < actin.n; i++) {
            _process_catch_bonds(i);
        }
        }
        stamp_section(kSecAaChemistry);

        // Reduce chemistry and catch-bond status together. Both use the same
        // actin index partition and retain their original thread summation
        // order, so the fused pass is bitwise equivalent.
        #pragma omp for
        for (int i = 0; i < actin.n; ++i) {
            for (int t = 0; t < omp_get_num_threads(); ++t) {
                cumulative_chem_entropy[i] +=
                    actin_chem_entropy_delta_temp[t][i];
                cumulative_chem_input[i] +=
                    actin_chem_input_delta_temp[t][i];
                cumulative_chem_binding_events[i] +=
                    actin_chem_binding_delta_temp[t][i];
                cumulative_chem_unbinding_events[i] +=
                    actin_chem_unbinding_delta_temp[t][i];
                actin.cb_status[i] = std::max(actin.cb_status[i], actin_cb_status_temp[t][i]);
            }
        }
        stamp_section(kSecChemReduce);

        #pragma omp single
        {
            _update_myosin_bond_matrix();
        }
        stamp_section(kSecBondMatrix);

        // Step 6: Compute actin-myosin forces
        #pragma omp for schedule(runtime)
        for (int i = 0; i < actin.n; i++) {
            _calc_am_force_velocity(i);
        }
        stamp_section(kSecAmForce);

        #pragma omp single
        {
            _record_aa_state_changes();
        }

        if (titin_k > 0 && titin_rest_length > 0 &&
            titin_mode != TitinMode::explicit_dynamic) {
            #pragma omp for schedule(runtime)
            for (int i = 0; i < actin.n; ++i) {
                _apply_titin_forces(i);
            }
    }

        // Preserve the historical k_mm=0 execution graph.  An unconditional
        // `single` has an implicit team barrier even when the callee returns
        // immediately; that changed which worker reached later stochastic
        // work first in otherwise identical multi-threaded controls.
        if (k_mm > 0.0) {
            #pragma omp single
            {
                _update_myomesin_bridges();
            }
        }
        stamp_section(kSecAaStateTitin);
        _volume_exclusion();
        stamp_section(kSecVolumeExclusion);
    // double k_theta = 1.0;
    // _apply_cb_alignment_bias(k_theta);

        // Step 8: Reduce actin forces and angular forces
        reduce_arrays(
            actin_forces_temp, actin.force,
            actin_torques_temp, actin.torque);

        // Step 9: Reduce myosin forces, velocities, and angular forces
        reduce_arrays(
            myosin_forces_temp, myosin.force,
            myosin_velocities_temp, myosin.velocity,
            myosin_torques_temp, myosin.torque);
        stamp_section(kSecReduce);

        if (titin_mode == TitinMode::explicit_dynamic && titin_k > 0.0) {
            #pragma omp single
            {
                _update_explicit_titin();
            }
        }


        if (bundle_strength > 0.0) {
            #pragma omp single
            {
                _apply_transverse_bundling(bundle_strength);
            }
        }

        if (probe_confine_k > 0.0) {
            #pragma omp single
            {
                _apply_probe_confinement();
            }
        }

        // Cap forces based on specified maxima
        #pragma omp for
        for (int i = 0; i < std::max(actin.n, myosin.n); ++i) {
            if (i < actin.n) {
                double cap = (actin.cb_status[i] < 2)
                    ? max_actin_force : max_myosin_force;
                if (std::isfinite(cap) && cap > 0.0) {
                    double mag = actin.force[i].norm();
                    if (mag > cap && mag > 0.0) {
                        actin.force[i] *= (cap / mag);
                        actin_force_cap_events[i] += 1.0;
                    }
                }
                const double torque_cap = max_actin_torque;
                if (std::isfinite(torque_cap) && torque_cap > 0.0) {
                    double mag = actin.torque[i].norm();
                    if (mag > torque_cap && mag > 0.0) {
                        actin.torque[i] *= (torque_cap / mag);
                        actin_torque_cap_events[i] += 1.0;
                    }
                }
            }
            if (i < myosin.n) {
                const double force_cap = max_myosin_force;
                if (std::isfinite(force_cap) && force_cap > 0.0) {
                    double mag = myosin.force[i].norm();
                    if (mag > force_cap && mag > 0.0) {
                        myosin.force[i] *= (force_cap / mag);
                        myosin_force_cap_events[i] += 1.0;
                    }
                }
                const double torque_cap = max_myosin_torque;
                if (std::isfinite(torque_cap) && torque_cap > 0.0) {
                    double mag = myosin.torque[i].norm();
                    if (mag > torque_cap && mag > 0.0) {
                        myosin.torque[i] *= (torque_cap / mag);
                        myosin_torque_cap_events[i] += 1.0;
                    }
                }
            }
        }

        #pragma omp for
        for (int myosin_idx = 0; myosin_idx < myosin.n; ++myosin_idx) {
            std::array<double, 2> max_load = myosin_f_load_temp[0][myosin_idx];
            for (int thread_idx = 1; thread_idx < omp_get_max_threads(); ++thread_idx) {
                const auto& thread_load = myosin_f_load_temp[thread_idx][myosin_idx];
                max_load[0] = std::max(max_load[0], thread_load[0]);
                max_load[1] = std::max(max_load[1], thread_load[1]);
            }
            myosin_f_load[myosin_idx] = max_load;
            double dot = myosin.velocity[myosin_idx].dot(
                myosin.direction[myosin_idx]);
            size_t dir_idx;
            if (std::abs(dot) < EPS) {
                dir_idx = (myosin_f_load[myosin_idx][0] >=
                           myosin_f_load[myosin_idx][1]) ? 0 : 1;
            } else {
                dir_idx = (dot >= 0.0) ? 0 : 1;
            }
            double load_factor = 1.0 - myosin_f_load[myosin_idx][dir_idx];
            if (load_factor < 0.0) {
                load_factor = 0.0;
            }
            myosin.velocity[myosin_idx].normalize();
            myosin.velocity[myosin_idx] *= v_am * load_factor;
        }
    }
    system_parallel_seconds += std::chrono::duration<double>(
        ProfileClock::now() - parallel_start).count();

    const auto finalize_start = ProfileClock::now();
    _finalize_active_pair_lists();

    for (auto& local_events : cb_breakage_events_temp) {
        if (!local_events.empty()) {
            cb_breakage_events.insert(cb_breakage_events.end(),
                                      local_events.begin(), local_events.end());
            local_events.clear();
        }
    }
    for (auto& local_lifetimes : aa_completed_lifetimes_temp) {
        if (!local_lifetimes.empty()) {
            aa_completed_lifetimes.insert(aa_completed_lifetimes.end(),
                                          local_lifetimes.begin(), local_lifetimes.end());
            local_lifetimes.clear();
        }
    }
    for (auto& local_events : am_transition_events_temp) {
        if (!local_events.empty()) {
            am_transition_events.insert(am_transition_events.end(),
                                        local_events.begin(), local_events.end());
            local_events.clear();
        }
    }
    for (auto& local_events : aa_state_events_temp) {
        if (!local_events.empty()) {
            aa_state_events.insert(aa_state_events.end(),
                                   local_events.begin(), local_events.end());
            local_events.clear();
        }
    }
    const auto system_end = ProfileClock::now();
    system_finalize_seconds += std::chrono::duration<double>(
        system_end - finalize_start).count();
    system_total_seconds += std::chrono::duration<double>(
        system_end - system_start).count();
    ++system_profile_steps;
}


void Sarcomere::update_system_sterics_only() {
    // Advance step counter as well for sterics-only updates
    current_step++;
    const double bundle_strength = _current_bundle_strength();
    auto& actin_torque_cap_events = actin["cumulative_torque_cap_events"];
    auto& myosin_torque_cap_events = myosin["cumulative_torque_cap_events"];
    _update_neighbors();
    #pragma omp parallel
    {   
        _set_to_zero();  
        // Relax all steric overlaps before any AM/AA chemistry or mature
        // springs are enabled. The old pre-phase separated only myosins.
        _volume_exclusion();
        // Step 7: Reduce actin forces and angular forces
        reduce_arrays(
            actin_forces_temp, actin.force,
            actin_torques_temp, actin.torque);
        // Step 8: Reduce myosin forces, velocities, and angular forces
        reduce_arrays(
            myosin_forces_temp, myosin.force,
            myosin_torques_temp, myosin.torque);

        if (bundle_strength > 0.0) {
            #pragma omp single
            {
                _apply_transverse_bundling(bundle_strength);
            }
        }

        if (probe_confine_k > 0.0) {
            #pragma omp single
            {
                _apply_probe_confinement();
            }
        }

        #pragma omp for
        for (int i = 0; i < actin.n; ++i) {
            double cap = max_actin_torque;
            if (std::isfinite(cap) && cap > 0.0) {
                double mag = actin.torque[i].norm();
                if (mag > cap && mag > 0.0) {
                    actin.torque[i] *= (cap / mag);
                    actin_torque_cap_events[i] += 1.0;
                }
            }
        }
        #pragma omp for
        for (int i = 0; i < myosin.n; ++i) {
            double cap = max_myosin_torque;
            if (std::isfinite(cap) && cap > 0.0) {
                double mag = myosin.torque[i].norm();
                if (mag > cap && mag > 0.0) {
                    myosin.torque[i] *= (cap / mag);
                    myosin_torque_cap_events[i] += 1.0;
                }
            }
        }
    }
    _finalize_active_pair_lists();
}


int Sarcomere::load_initial_geometry(const std::string& source_filename,
                                     int frame_index) {
    if (source_filename.empty()) {
        throw std::invalid_argument(
            "initial geometry source filename must not be empty");
    }
    if (frame_index < -1) {
        throw std::invalid_argument(
            "initial geometry frame must be -1 or a non-negative index");
    }

    try {
        H5::H5File source(source_filename, H5F_ACC_RDONLY);

        struct DatasetShape {
            H5::DataSet dataset;
            hsize_t frames = 0;
            hsize_t filaments = 0;
        };
        auto inspect_dataset = [&](const char* path,
                                   int expected_filaments) -> DatasetShape {
            if (!source.nameExists(path)) {
                throw std::runtime_error(
                    std::string("missing geometry dataset ") + path);
            }
            H5::DataSet dataset = source.openDataSet(path);
            H5::DataSpace space = dataset.getSpace();
            if (space.getSimpleExtentNdims() != 3) {
                throw std::runtime_error(
                    std::string("geometry dataset must have rank 3: ") + path);
            }
            hsize_t dims[3] = {0, 0, 0};
            space.getSimpleExtentDims(dims, nullptr);
            if (dims[0] == 0 ||
                dims[1] != static_cast<hsize_t>(expected_filaments) ||
                dims[2] != 3) {
                throw std::runtime_error(
                    std::string("incompatible geometry dataset shape for ") + path);
            }
            return {dataset, dims[0], dims[1]};
        };

        DatasetShape actin_center = inspect_dataset("/actin/center", actin.n);
        DatasetShape actin_direction =
            inspect_dataset("/actin/direction", actin.n);
        DatasetShape myosin_center =
            inspect_dataset("/myosin/center", myosin.n);
        DatasetShape myosin_direction =
            inspect_dataset("/myosin/direction", myosin.n);

        hsize_t usable_frames = std::min(
            std::min(actin_center.frames, actin_direction.frames),
            std::min(myosin_center.frames, myosin_direction.frames));
        if (source.nameExists("/meta/committed_frames")) {
            long long committed = 0;
            H5::DataSet committed_dataset =
                source.openDataSet("/meta/committed_frames");
            committed_dataset.read(&committed, H5::PredType::NATIVE_LLONG);
            if (committed <= 0) {
                throw std::runtime_error(
                    "source trajectory has no committed geometry frames");
            }
            usable_frames = std::min(
                usable_frames, static_cast<hsize_t>(committed));
        }
        if (usable_frames == 0) {
            throw std::runtime_error(
                "source trajectory has no complete geometry frames");
        }

        const hsize_t resolved_frame = frame_index < 0
            ? usable_frames - 1
            : static_cast<hsize_t>(frame_index);
        if (resolved_frame >= usable_frames) {
            throw std::out_of_range(
                "requested initial geometry frame is outside the committed source trajectory");
        }

        // A box mismatch would change periodic images during the transplant.
        // Older trajectories without box metadata remain importable.
        if (source.nameExists("/meta")) {
            H5::Group meta = source.openGroup("/meta");
            const char* names[3] = {"box_x", "box_y", "box_z"};
            for (int axis = 0; axis < 3; ++axis) {
                if (!meta.attrExists(names[axis])) {
                    continue;
                }
                double source_box = 0.0;
                H5::Attribute attr = meta.openAttribute(names[axis]);
                attr.read(H5::PredType::NATIVE_DOUBLE, &source_box);
                const double scale = std::max(1.0, std::abs(box[axis]));
                if (!std::isfinite(source_box) ||
                    std::abs(source_box - box[axis]) > 1e-9 * scale) {
                    throw std::runtime_error(
                        std::string("source/current box mismatch on ") +
                        names[axis]);
                }
            }
        }

        auto read_frame = [&](DatasetShape& shaped) {
            std::vector<double> values(
                static_cast<size_t>(shaped.filaments) * 3);
            H5::DataSpace file_space = shaped.dataset.getSpace();
            hsize_t start[3] = {resolved_frame, 0, 0};
            hsize_t count[3] = {1, shaped.filaments, 3};
            file_space.selectHyperslab(H5S_SELECT_SET, count, start);
            H5::DataSpace memory_space(3, count);
            shaped.dataset.read(values.data(), H5::PredType::NATIVE_DOUBLE,
                                memory_space, file_space);
            return values;
        };

        const auto actin_center_values = read_frame(actin_center);
        const auto actin_direction_values = read_frame(actin_direction);
        const auto myosin_center_values = read_frame(myosin_center);
        const auto myosin_direction_values = read_frame(myosin_direction);

        auto import_species = [&](Filament& filament,
                                  const std::vector<double>& centers,
                                  const std::vector<double>& directions,
                                  const char* species_name) {
            for (int i = 0; i < filament.n; ++i) {
                const size_t offset = static_cast<size_t>(i) * 3;
                vec center{centers[offset], centers[offset + 1],
                           centers[offset + 2]};
                vec direction{directions[offset], directions[offset + 1],
                              directions[offset + 2]};
                const double norm = direction.norm();
                if (!std::isfinite(center.x) || !std::isfinite(center.y) ||
                    !std::isfinite(center.z) || !std::isfinite(norm) ||
                    norm < 1e-12) {
                    throw std::runtime_error(
                        std::string("non-finite or degenerate ") + species_name +
                        " geometry at filament " + std::to_string(i));
                }
                filament.center[static_cast<size_t>(i)] = center;
                filament.direction[static_cast<size_t>(i)] = direction / norm;
            }
            filament.update_endpoints();
            std::fill(filament.force_x.begin(), filament.force_x.end(), 0.0);
            std::fill(filament.force_y.begin(), filament.force_y.end(), 0.0);
            std::fill(filament.force_z.begin(), filament.force_z.end(), 0.0);
            std::fill(filament.torque_x.begin(), filament.torque_x.end(), 0.0);
            std::fill(filament.torque_y.begin(), filament.torque_y.end(), 0.0);
            std::fill(filament.torque_z.begin(), filament.torque_z.end(), 0.0);
            std::fill(filament.velocity_x.begin(), filament.velocity_x.end(), 0.0);
            std::fill(filament.velocity_y.begin(), filament.velocity_y.end(), 0.0);
            std::fill(filament.velocity_z.begin(), filament.velocity_z.end(), 0.0);
        };

        import_species(actin, actin_center_values, actin_direction_values,
                       "actin");
        import_species(myosin, myosin_center_values, myosin_direction_values,
                       "myosin");

        neighbor_list.set_species_positions(
            actin.center_x, actin.center_y, actin.center_z,
            myosin.center_x, myosin.center_y, myosin.center_z);
        neighbor_list.rebuild_neighbor_list();

        initial_geometry_source = source_filename;
        initial_geometry_source_frame_requested = frame_index;
        initial_geometry_source_frame_loaded =
            static_cast<int>(resolved_frame);
        return initial_geometry_source_frame_loaded;
    } catch (const H5::Exception& error) {
        throw std::runtime_error(
            "failed to read initial geometry from " + source_filename +
            ": " + error.getDetailMsg());
    }
}


void Sarcomere::set_periodicity(const std::array<bool,3>& periodic_axes) {
    is_periodic = periodic_axes;
    actin.set_periodic_axes(is_periodic);
    myosin.set_periodic_axes(is_periodic);
    neighbor_list.set_periodic_axes(is_periodic);
    neighbor_list.rebuild_neighbor_list();
}


void Sarcomere::_update_neighbors() {
    neighbor_list.set_species_positions(actin.center_x, actin.center_y, actin.center_z,
                                       myosin.center_x, myosin.center_y, myosin.center_z);
    // Clear previous neighbors data and prepare to store new neighbors
    actin_neighbors_by_species.clear();
    actin_neighbors_by_species.resize(actin.center.size());
    // Sized here, outside the parallel region, for the same reason as the list
    // above.  Guarding a resize with `if (size != n) { #pragma omp single ... }`
    // inside the region is a race: one thread can complete the resize before
    // another evaluates the condition, so that thread never encounters the
    // `single` and its implicit barrier, and it reads .size() while the vector
    // is being reallocated.
    actin_unit_direction.resize(actin.center.size());
    if (neighbor_list.needs_rebuild()) {
        neighbor_list.rebuild_neighbor_list();
        actin_close_valid = false;
    }
    if (!actin_close_valid ||
        actin_close_neighbors.size() != actin.center.size()) {
        _refresh_actin_close_neighbors();
    }
}

// Segment-distance-filtered actin-actin candidates for the sterics loop.
//
// The shared broadphase indexes filament CENTRES, so its cutoff is
// max(filament length) + interaction range = 1.6 um: each actin carries ~857
// candidates of which ~4 are ever in range.  Filtering has to be CACHED to be
// worth anything - a rejected pair already costs only one segment-distance
// evaluation inside compute_actin_repulsion, so re-filtering every step would
// save nothing.  The list is therefore rebuilt only when the cell list is, and
// carries the full skin so it stays conservative until then.
//
// Called from _update_neighbors(), i.e. OUTSIDE the omp parallel region, and
// reads the cell list directly rather than actin_neighbors_by_species (which
// is not populated until _set_to_zero runs inside that region).
void Sarcomere::_refresh_actin_close_neighbors() {
    const double skin = 2.0 * neighbor_list.rebuild_threshold();
    double reach = 2.0 * actin_radius;                       // capsule sterics
    if (aa_double_overlap_contact_enabled) {
        reach = std::max(reach, aa_double_overlap_contact_diameter);
    }
    if (aa_barbed_lateral_range > 0.0 && aa_barbed_lateral_k > 0.0) {
        // furthest reach is the corner of the lateral term's envelope
        reach = std::max(reach, std::hypot(aa_barbed_lateral_axial_window,
                                           aa_barbed_lateral_range));
    }
    // compute_opposing_pointed_end_contact places NO lower bound on the
    // signed axial gap, so it can fire for an antiparallel pair whose
    // pointed ends are transversely aligned at any axial separation.  No
    // segment-distance bound is therefore rigorous for it, and a filament
    // half-length is added as an empirical margin: a cutoff sweep showed
    // bitwise-identical trajectories from 0.60 um upward and divergence at
    // 0.27 um, with no measurable cost difference between them.
    actin_close_cutoff = reach + skin + 0.5 * actin.length;
    if (const char* override_cutoff = std::getenv("SARC_CLOSE_CUTOFF")) {
        actin_close_cutoff = std::atof(override_cutoff);
    }

    const size_t n = actin.center.size();
    actin_close_neighbors.assign(n, {});
    #pragma omp parallel for schedule(dynamic, 32)
    for (size_t i = 0; i < n; ++i) {
        auto& keep = actin_close_neighbors[i];
        for (const auto& [neighbor_index, neighbor_type] :
             neighbor_list.get_neighbors(static_cast<int>(i))) {
            if (neighbor_type != ParticleType::Actin) continue;
            const int j = neighbor_index;
            // Filtering preserves relative order and a dropped pair would have
            // contributed exactly zero, so the sterics loop is unchanged.
            bool keep_pair =
                geometry::segment_segment_distance(
                    actin.left_end[i], actin.right_end[i],
                    actin.left_end[j], actin.right_end[j],
                    box, is_periodic) <= actin_close_cutoff;
            // The pointed-end contact is NOT bounded by segment distance: it
            // fires whenever two antiparallel pointed ends are transversely
            // aligned, with no lower bound on the (signed) axial gap, so a
            // deeply crossed pair qualifies at any axial separation.  Test its
            // own geometry rather than truncating it.
            if (!keep_pair && aa_double_overlap_contact_enabled) {
                vec ui = actin.direction[i];
                vec uj = actin.direction[j];
                if (ui.norm() > 0.0 && uj.norm() > 0.0) {
                    ui.normalize();
                    uj.normalize();
                    if (ui.dot(uj) <= -0.8 + 1e-9) {
                        vec e = ui - uj;
                        if (e.norm() > 1e-12) {
                            e.normalize();
                            vec sep = (static_cast<vec>(actin.center[j]) +
                                       0.5 * actin.length * uj) -
                                      (static_cast<vec>(actin.center[i]) +
                                       0.5 * actin.length * ui);
                            sep.pbc_wrap(box, actin.periodic_axes);
                            const vec transverse = sep - sep.dot(e) * e;
                            keep_pair = transverse.norm() <=
                                aa_double_overlap_contact_diameter + skin;
                        }
                    }
                }
            }
            if (keep_pair) keep.push_back(j);
        }
    }
    actin_close_valid = true;
}

bool Sarcomere::_actin_renewal_position_is_clear(
    int actin_index, const vec& center, const vec& direction) const {
    const vec first = center - 0.5 * actin.length * direction;
    const vec second = center + 0.5 * actin.length * direction;
    const double aa_clearance =
        actin_turnover_clearance_ratio * 2.0 * actin_radius;
    for (int other = 0; other < actin.n; ++other) {
        if (other == actin_index) continue;
        const double distance = geometry::segment_segment_distance(
            first, second, actin.left_end[other], actin.right_end[other],
            box, is_periodic);
        if (distance < aa_clearance - 1e-12) return false;
    }
    const double am_clearance =
        actin_turnover_clearance_ratio * (actin_radius + myosin.radius);
    for (int m = 0; m < myosin.n; ++m) {
        const double distance = geometry::segment_segment_distance(
            first, second, myosin.left_end[m], myosin.right_end[m],
            box, is_periodic);
        if (distance < am_clearance - 1e-12) return false;
    }
    return true;
}

bool Sarcomere::_renew_actin(
    int actin_index, int& placement_attempts, int& released_am,
    int& released_actn2, int& released_titin, vec& new_center,
    vec& new_direction) {
    placement_attempts = 0;
    released_am = 0;
    released_actn2 = 0;
    released_titin = 0;
    const vec old_direction = actin.direction[actin_index];
    const bool old_positive_x = old_direction.x >= 0.0;
    counter_rng::CounterRng placement_rng(
        static_cast<std::uint64_t>(initial_seed),
        static_cast<std::uint64_t>(current_step),
        counter_rng::Stream::actin_turnover,
        static_cast<std::uint64_t>(actin_index), 1);
    bool placed = false;
    for (int attempt = 1; attempt <= actin_turnover_max_placement_attempts;
         ++attempt) {
        placement_attempts = attempt;
        new_center = {
            (placement_rng.uniform() - 0.5) * box[0],
            (placement_rng.uniform() - 0.5) * box[1],
            (placement_rng.uniform() - 0.5) * box[2]};
        const double axial = 2.0 * placement_rng.uniform() - 1.0;
        const double azimuth = 2.0 * M_PI * placement_rng.uniform();
        const double transverse = std::sqrt(std::max(0.0, 1.0 - axial * axial));
        new_direction = {
            axial, transverse * std::cos(azimuth),
            transverse * std::sin(azimuth)};
        if ((new_direction.x >= 0.0) != old_positive_x) {
            new_direction = -1.0 * new_direction;
        }
        if (_actin_renewal_position_is_clear(
                actin_index, new_center, new_direction)) {
            placed = true;
            break;
        }
    }
    if (!placed) return false;

    // Release every explicit motor head owned by a thick filament.
    for (int m = 0; m < myosin.n; ++m) {
        int pair_count = 0;
        for (const auto& half : am_heads_by_myosin[static_cast<size_t>(m)]) {
            for (const auto& head : half) {
                if (head.bound_actin == actin_index) ++pair_count;
            }
        }
        if (pair_count > 0) {
            released_am += pair_count;
            _detach_am_heads_for_pair(
                actin_index, m, AMTransitionReason::actin_turnover);
        }
        am_bonds[actin_index][m] = 0;
        am_bonds_prev[actin_index][m] = 0;
        am_state[actin_index][m] = static_cast<BondFlag>(
            am_motor_kinetics::State::detached);
        am_seen[actin_index][m] = 0;
        am_target_axial[actin_index][m] = 0.0;
        actinIndicesPerMyosin.deleteConnection(m, actin_index);
        for (auto& temporary : actinIndicesPerMyosin_temp) {
            temporary.deleteConnection(m, actin_index);
        }
    }
    myosinIndicesPerActin.deleteAllConnections(actin_index);

    // Release every local ACTN2 molecule, including single-arm tethers.
    for (size_t slot = 0; slot < aa_links.size(); ++slot) {
        auto& link = aa_links[slot];
        if (!link.active ||
            (link.actin_i != actin_index && link.actin_j != actin_index)) {
            continue;
        }
        ++released_actn2;
        link.active = false;
        aa_free_link_slots.push_back(slot);
    }
    for (int other = 0; other < actin.n; ++other) {
        actin_actin_bonds[actin_index][other] = 0;
        actin_actin_bonds[other][actin_index] = 0;
        actin_actin_bonds_prev[actin_index][other] = 0;
        actin_actin_bonds_prev[other][actin_index] = 0;
        actin_actin_status[actin_index][other] = 0;
        actin_actin_status[other][actin_index] = 0;
        actin_actin_status_prev[actin_index][other] = 0;
        actin_actin_status_prev[other][actin_index] = 0;
        actin_actin_lifetime[actin_index][other] = 0;
        actin_actin_lifetime[other][actin_index] = 0;
        actin_actin_lifetime_prev[actin_index][other] = 0;
        actin_actin_lifetime_prev[other][actin_index] = 0;
        actin_recovery_until[actin_index][other] = 0;
        actin_recovery_until[other][actin_index] = 0;
        aa_status2_memory_until[actin_index][other] = 0;
        aa_status2_memory_until[other][actin_index] = 0;
        aa_attach_step[actin_index][other] = -1;
        aa_attach_step[other][actin_index] = -1;
        if (aa_neighbor_identity_enabled) {
            aa_neighbor_slots.erase(actin_index, other);
        }
    }
    auto remove_actin_pairs = [actin_index](std::vector<BondPair>& pairs) {
        pairs.erase(std::remove_if(
            pairs.begin(), pairs.end(), [actin_index](const BondPair& pair) {
                return pair.first == actin_index || pair.second == actin_index;
            }), pairs.end());
    };
    remove_actin_pairs(active_aa_pairs);
    remove_actin_pairs(previous_aa_pairs);
    remove_actin_pairs(active_recovery_pairs);
    remove_actin_pairs(active_am_pairs);
    remove_actin_pairs(previous_am_pairs);
    _rebuild_local_aa_indices();

    // Explicit titin is a molecular attachment and must not follow the new
    // filament geometry.  Candidate-only references are also invalidated.
    for (auto& titin : explicit_titins) {
        if (titin.bound && titin.bound_actin == actin_index) {
            ++released_titin;
            titin.bound = false;
            titin.last_bound_actin = actin_index;
            titin.bound_actin = -1;
            titin.last_detach_step = static_cast<RecoveryStep>(std::min(
                current_step,
                static_cast<size_t>(std::numeric_limits<RecoveryStep>::max())));
            titin.distance = titin.extension = titin.force = 0.0;
        }
        if (titin.candidate_actin == actin_index) titin.candidate_actin = -1;
    }
    if (explicit_titin_owner_by_actin.size() == static_cast<size_t>(actin.n)) {
        explicit_titin_owner_by_actin[static_cast<size_t>(actin_index)] = -1;
    }

    actin.center[actin_index] = new_center;
    actin.direction[actin_index] = new_direction;
    actin.update_endpoints(actin_index);
    actin.force[actin_index] = {0.0, 0.0, 0.0};
    actin.torque[actin_index] = {0.0, 0.0, 0.0};
    actin.velocity[actin_index] = {0.0, 0.0, 0.0};
    actin.f_load[actin_index] = 0.0;
    actin.cb_status[actin_index] = 0;
    previous_actin_cb_status[actin_index] = 0;
    const double reset_q = am_thin_activation_enabled
        ? am_motor_kinetics::thin_activation_calcium_equilibrium(
              am_activation, am_q_calcium_open_rate, am_q_calcium_half,
              am_q_calcium_hill, am_q_close_rate)
        : am_activation;
    am_thin_activation[static_cast<size_t>(actin_index)] = reset_q;
    actin["am_q"][static_cast<size_t>(actin_index)] = reset_q;
    actin["am_activation"][static_cast<size_t>(actin_index)] = am_activation;
    actin["am_calcium"][static_cast<size_t>(actin_index)] = am_activation;
    actin_last_turnover_step[static_cast<size_t>(actin_index)] =
        static_cast<RecoveryStep>(std::min(
            current_step,
            static_cast<size_t>(std::numeric_limits<RecoveryStep>::max())));
    return true;
}

void Sarcomere::_process_actin_turnover() {
    if (actin_last_turnover_step.size() != static_cast<size_t>(actin.n)) {
        actin_last_turnover_step.assign(static_cast<size_t>(actin.n), 0);
    }
    auto& age_feature = actin["turnover_age"];
    for (int i = 0; i < actin.n; ++i) {
        age_feature[static_cast<size_t>(i)] = dt * static_cast<double>(
            current_step - std::min<size_t>(
                current_step,
                actin_last_turnover_step[static_cast<size_t>(i)]));
    }
    if (!actin_turnover_enabled || !(actin_turnover_rate > 0.0) ||
        static_cast<double>(current_step) * dt < actin_turnover_start_time) {
        return;
    }
    const double probability = -std::expm1(-actin_turnover_rate * dt);
    for (int i = 0; i < actin.n; ++i) {
        counter_rng::CounterRng event_rng(
            static_cast<std::uint64_t>(initial_seed),
            static_cast<std::uint64_t>(current_step),
            counter_rng::Stream::actin_turnover,
            static_cast<std::uint64_t>(i), 0);
        if (event_rng.uniform() >= probability) continue;
        const vec old_center = actin.center[i];
        const vec old_direction = actin.direction[i];
        const size_t last_step = std::min<size_t>(
            current_step, actin_last_turnover_step[static_cast<size_t>(i)]);
        int attempts = 0;
        int released_am = 0;
        int released_actn2 = 0;
        int released_titin = 0;
        vec new_center = old_center;
        vec new_direction = old_direction;
        const bool success = _renew_actin(
            i, attempts, released_am, released_actn2, released_titin,
            new_center, new_direction);
        if (success) {
            actin["cumulative_turnover_events"][static_cast<size_t>(i)] += 1.0;
            age_feature[static_cast<size_t>(i)] = 0.0;
        } else {
            new_center = old_center;
            new_direction = old_direction;
        }
        actin_turnover_events.insert(actin_turnover_events.end(), {
            static_cast<double>(i), static_cast<double>(current_step),
            static_cast<double>(current_step - last_step),
            success ? 1.0 : 0.0, static_cast<double>(attempts),
            static_cast<double>(released_am),
            static_cast<double>(released_actn2),
            static_cast<double>(released_titin),
            old_center.x, old_center.y, old_center.z,
            new_center.x, new_center.y, new_center.z,
            old_direction.x, old_direction.y, old_direction.z,
            new_direction.x, new_direction.y, new_direction.z});
    }
}


void Sarcomere::_set_to_zero() {
    #pragma omp single
    {
        for (const auto& [i, j] : previous_aa_pairs) {
            actin_actin_bonds_prev[i][j] = actin_actin_bonds_prev[j][i] = 0;
            actin_actin_status_prev[i][j] = actin_actin_status_prev[j][i] = 0;
            actin_actin_lifetime_prev[i][j] = actin_actin_lifetime_prev[j][i] = 0;
        }
        previous_aa_pairs.clear();
        previous_aa_pairs.reserve(active_aa_pairs.size());
        for (const auto& [i, j] : active_aa_pairs) {
            actin_actin_bonds_prev[i][j] = actin_actin_bonds_prev[j][i] =
                actin_actin_bonds[i][j];
            actin_actin_status_prev[i][j] = actin_actin_status_prev[j][i] =
                actin_actin_status[i][j];
            actin_actin_lifetime_prev[i][j] = actin_actin_lifetime_prev[j][i] =
                actin_actin_lifetime[i][j];
            actin_actin_bonds[i][j] = actin_actin_bonds[j][i] = 0;
            actin_actin_status[i][j] = actin_actin_status[j][i] = 0;
            actin_actin_lifetime[i][j] = actin_actin_lifetime[j][i] = 0;
            previous_aa_pairs.emplace_back(i, j);
        }
        active_aa_pairs.clear();
        std::fill(aa_reserved_degree.begin(), aa_reserved_degree.end(), 0);
        for (const auto& [i, j] : previous_aa_pairs) {
            ++aa_reserved_degree[i];
            ++aa_reserved_degree[j];
        }

        for (const auto& [i, j] : previous_am_pairs) {
            am_bonds_prev[i][j] = 0;
        }
        previous_am_pairs.clear();
        previous_am_pairs.reserve(active_am_pairs.size());
        for (const auto& [i, j] : active_am_pairs) {
            am_bonds_prev[i][j] = am_bonds[i][j];
            am_bonds[i][j] = 0;
            if (dynamic_am) am_seen[i][j] = 0;
            previous_am_pairs.emplace_back(i, j);
        }
        active_am_pairs.clear();

        for (auto& pairs : active_aa_pairs_temp) pairs.clear();
        for (auto& pairs : active_am_pairs_temp) pairs.clear();
        for (auto& pairs : recovery_pairs_temp) pairs.clear();
    }

    #pragma omp for
    for (int t = 0; t < omp_get_max_threads(); ++t) {
        cb_breakage_events_temp[t].clear();
        aa_completed_lifetimes_temp[t].clear();
        am_transition_events_temp[t].clear();
        aa_state_events_temp[t].clear();
        for (int i = 0; i < actin.n; ++i) {
            actin_forces_temp[t][i] = {0, 0, 0};
            actin_torques_temp[t][i] = {0, 0, 0};
            actin_cb_status_temp[t][i] = 0;
            actin_chem_entropy_delta_temp[t][i] = 0.0;
            actin_chem_input_delta_temp[t][i] = 0.0;
            actin_chem_binding_delta_temp[t][i] = 0.0;
            actin_chem_unbinding_delta_temp[t][i] = 0.0;
        }
        for (int i = 0; i < myosin.n; ++i) {
            myosin_forces_temp[t][i] = {0, 0, 0};
            myosin_velocities_temp[t][i] = {0, 0, 0};
            myosin_torques_temp[t][i] = {0, 0, 0};
            myosin_f_load_temp[t][i] = {0.0, 0.0};
            actinIndicesPerMyosin_temp[t].deleteAllConnections(i);
        }
    }

    auto& myosin_binding_ratio_feature = actin["myosin_binding_ratio"];
    auto& crosslink_ratio_feature = actin["crosslink_ratio"];
    auto& partial_binding_ratio_feature = actin["partial_binding_ratio"];
    auto& am_weak_count_feature = actin["am_weak_count"];
    auto& am_strong_count_feature = actin["am_strong_count"];
    auto& am_motor_force_feature = actin["am_motor_force"];
    auto& am_active_force_x_feature = actin["am_active_force_x"];
    auto& am_active_force_y_feature = actin["am_active_force_y"];
    auto& am_active_force_z_feature = actin["am_active_force_z"];
    auto& am_tether_force_feature = actin["am_tether_force"];
    auto& am_head_transverse_force_feature = actin["am_head_transverse_force"];
    auto& am_activation_feature = actin["am_activation"];
    auto& am_calcium_feature = actin["am_calcium"];
    auto& am_q_feature = actin["am_q"];
    auto& am_overlap_available_feature = actin["am_overlap_available"];
    auto& am_overlap_engaged_feature = actin["am_overlap_engaged"];
    auto& am_overlap_axis_weighted_feature =
        actin["am_overlap_axis_weighted"];
    auto& am_overlap_productive_weighted_feature =
        actin["am_overlap_productive_weighted"];
    auto& aa_arm_load_max_feature = actin["aa_arm_load_max"];
    auto& aa_bilateral_load_max_feature = actin["aa_bilateral_load_max"];
    auto& aa_raw_arm_load_max_feature = actin["aa_raw_arm_load_max"];
    auto& aa_local_link_count_feature = actin["aa_local_link_count"];
    auto& aa_local_status1_count_feature = actin["aa_local_status1_count"];
    auto& aa_local_status2_count_feature = actin["aa_local_status2_count"];
    auto& aa_local_tension_max_feature = actin["aa_local_tension_max"];
    auto& aa_zdisc_tip_force_feature = actin["aa_zdisc_tip_force"];
    auto& aa_zdisc_tip_contact_count_feature =
        actin["aa_zdisc_tip_contact_count"];
    #pragma omp for
    for (int i = 0; i < actin.n; i++) {
        actin.update_endpoints(i);
        previous_actin_cb_status[i] = actin.cb_status[i];
        actin.force[i] = {0, 0, 0};
        actin.torque[i] = {0, 0, 0};
        actin.velocity[i] = {0, 0, 0};
        actin.f_load[i] = 0;
        if (actin_f_load_cb) {
            (*actin_f_load_cb)[i] = 0;
        }
        actin.cb_status[i] = 0;
        actin_basic_tension[i] = 0;
        actin_n_bonds[i] = 0;
        actin_strong_cb_count[i] = 0;
        n_myosins_per_actin[i] = 0;
        actin_crosslink_ratio[i] = 1;
        actin_crosslink_start[i] = actin.left_end[i];
        actin_crosslink_end[i] = actin.right_end[i];
        myosin_binding_ratio_feature[i] = 0;
        crosslink_ratio_feature[i] = 1;
        partial_binding_ratio_feature[i] = 0;
        am_weak_count_feature[i] = 0;
        am_strong_count_feature[i] = 0;
        am_motor_force_feature[i] = 0;
        am_active_force_x_feature[i] = 0;
        am_active_force_y_feature[i] = 0;
        am_active_force_z_feature[i] = 0;
        am_tether_force_feature[i] = 0;
        am_head_transverse_force_feature[i] = 0;
        am_activation_feature[i] = am_activation;
        am_calcium_feature[i] = am_activation;
        am_q_feature[i] = am_thin_activation[static_cast<size_t>(i)];
        am_overlap_available_feature[i] = 0;
        am_overlap_engaged_feature[i] = 0;
        am_overlap_axis_weighted_feature[i] = 0;
        am_overlap_productive_weighted_feature[i] = 0;
        aa_arm_load_max_feature[i] = 0;
        aa_bilateral_load_max_feature[i] = 0;
        aa_raw_arm_load_max_feature[i] = 0;
        aa_local_link_count_feature[i] = 0;
        aa_local_status1_count_feature[i] = 0;
        aa_local_status2_count_feature[i] = 0;
        aa_local_tension_max_feature[i] = 0;
        aa_zdisc_tip_force_feature[i] = 0;
        aa_zdisc_tip_contact_count_feature[i] = 0;
        actin_f_load_computed[i] = false;
        myosinIndicesPerActin.deleteAllConnections(i);
    }

    auto& align_drive_minus = myosin["am_align_drive_minus"];
    auto& align_drive_plus = myosin["am_align_drive_plus"];
    auto& align_abs_drive_minus = myosin["am_align_abs_drive_minus"];
    auto& align_abs_drive_plus = myosin["am_align_abs_drive_plus"];
    auto& align_active_drive_minus = myosin["am_align_active_drive_minus"];
    auto& align_active_drive_plus = myosin["am_align_active_drive_plus"];
    auto& align_tether_drive_minus = myosin["am_align_tether_drive_minus"];
    auto& align_tether_drive_plus = myosin["am_align_tether_drive_plus"];
    auto& align_strong_drive_minus = myosin["am_align_strong_drive_minus"];
    auto& align_strong_drive_plus = myosin["am_align_strong_drive_plus"];
    auto& align_head_count_minus = myosin["am_align_head_count_minus"];
    auto& align_head_count_plus = myosin["am_align_head_count_plus"];
    auto& align_positive_count_minus =
        myosin["am_align_positive_count_minus"];
    auto& align_positive_count_plus = myosin["am_align_positive_count_plus"];
    auto& align_negative_count_minus =
        myosin["am_align_negative_count_minus"];
    auto& align_negative_count_plus = myosin["am_align_negative_count_plus"];
    auto& zdisc_force_minus = myosin["aa_zdisc_tip_force_minus"];
    auto& zdisc_force_plus = myosin["aa_zdisc_tip_force_plus"];
    auto& zdisc_gap_minus = myosin["aa_zdisc_tip_gap_minus"];
    auto& zdisc_gap_plus = myosin["aa_zdisc_tip_gap_plus"];
    auto& zdisc_rho_minus = myosin["aa_zdisc_tip_rho_minus"];
    auto& zdisc_rho_plus = myosin["aa_zdisc_tip_rho_plus"];
    auto& zdisc_gate_minus = myosin["aa_zdisc_tip_gate_minus"];
    auto& zdisc_gate_plus = myosin["aa_zdisc_tip_gate_plus"];
    #pragma omp for
    for (int i = 0; i < myosin.n; i++) {
        myosin.update_endpoints(i);
        myosin.force[i] = {0, 0, 0};
        myosin.torque[i] = {0, 0, 0};
        myosin.velocity[i] = {0, 0, 0};
        actinIndicesPerMyosin.deleteAllConnections(i);
        myosin_f_load[i] = {0.0, 0.0};
        align_drive_minus[i] = 0.0;
        align_drive_plus[i] = 0.0;
        align_abs_drive_minus[i] = 0.0;
        align_abs_drive_plus[i] = 0.0;
        align_active_drive_minus[i] = 0.0;
        align_active_drive_plus[i] = 0.0;
        align_tether_drive_minus[i] = 0.0;
        align_tether_drive_plus[i] = 0.0;
        align_strong_drive_minus[i] = 0.0;
        align_strong_drive_plus[i] = 0.0;
        align_head_count_minus[i] = 0.0;
        align_head_count_plus[i] = 0.0;
        align_positive_count_minus[i] = 0.0;
        align_positive_count_plus[i] = 0.0;
        align_negative_count_minus[i] = 0.0;
        align_negative_count_plus[i] = 0.0;
        zdisc_force_minus[i] = 0.0;
        zdisc_force_plus[i] = 0.0;
        zdisc_gap_minus[i] = std::numeric_limits<double>::quiet_NaN();
        zdisc_gap_plus[i] = std::numeric_limits<double>::quiet_NaN();
        zdisc_rho_minus[i] = std::numeric_limits<double>::quiet_NaN();
        zdisc_rho_plus[i] = std::numeric_limits<double>::quiet_NaN();
        zdisc_gate_minus[i] = 0.0;
        zdisc_gate_plus[i] = 0.0;
    }
    #pragma omp for
    for (size_t i = 0; i < actin.center.size(); i++) {
        // Same expression the steric kernels use, so the gates below see the
        // identical bits.  Piggy-backed on this loop rather than a new `omp
        // for`, which would add a team barrier and perturb the execution graph.
        vec unit = actin.direction[i];
        if (unit.norm() > 0.0) unit.normalize();
        actin_unit_direction[i] = unit;
        // Refill in place so stable neighborhoods reuse their vector capacity
        // instead of allocating two temporary vectors per actin per timestep.
        auto& by_species = actin_neighbors_by_species[i];
        by_species.first.clear();
        by_species.second.clear();
        for (const auto& [neighbor_index, neighbor_type] :
             neighbor_list.get_neighbors(static_cast<int>(i))) {
            if (neighbor_type == ParticleType::Actin) {
                by_species.first.push_back(neighbor_index);
            } else {
                by_species.second.push_back(neighbor_index - actin.n);
            }
        }
    }
}

void Sarcomere::_finalize_active_pair_lists() {
    auto merge_pairs = [](auto& destination, auto& per_thread, bool normalize,
                          auto&& is_active) {
        size_t candidate_count = 0;
        for (const auto& pairs : per_thread) candidate_count += pairs.size();
        destination.clear();
        destination.reserve(candidate_count);
        for (const auto& pairs : per_thread) {
            for (auto pair : pairs) {
                if (normalize && pair.first > pair.second) {
                    std::swap(pair.first, pair.second);
                }
                if (is_active(pair.first, pair.second)) {
                    destination.push_back(pair);
                }
            }
        }
        std::sort(destination.begin(), destination.end());
        destination.erase(std::unique(destination.begin(), destination.end()), destination.end());
    };

    merge_pairs(active_aa_pairs, active_aa_pairs_temp, true,
        [&](int i, int j) { return actin_actin_bonds[i][j] == 1; });
    merge_pairs(active_am_pairs, active_am_pairs_temp, false,
        [&](int i, int j) { return am_bonds[i][j] == 1; });
    std::vector<BondPair> recovery_candidates;
    recovery_candidates.reserve(active_recovery_pairs.size());
    for (const auto& pair : active_recovery_pairs) {
        if (actin_recovery_until[pair.first][pair.second] > current_step) {
            recovery_candidates.push_back(pair);
        }
    }
    for (const auto& pairs : recovery_pairs_temp) {
        for (auto pair : pairs) {
            if (pair.first > pair.second) std::swap(pair.first, pair.second);
            recovery_candidates.push_back(pair);
        }
    }
    std::sort(recovery_candidates.begin(), recovery_candidates.end());
    recovery_candidates.erase(
        std::unique(recovery_candidates.begin(), recovery_candidates.end()),
        recovery_candidates.end());
    active_recovery_pairs.clear();
    active_recovery_pairs.reserve(recovery_candidates.size());
    for (auto pair : recovery_candidates) {
        if (pair.first > pair.second) std::swap(pair.first, pair.second);
        if (actin_recovery_until[pair.first][pair.second] > current_step) {
            active_recovery_pairs.push_back(pair);
        }
    }
}

void Sarcomere::_rebuild_active_pair_lists_from_dense() {
    active_aa_pairs.clear();
    previous_aa_pairs.clear();
    active_am_pairs.clear();
    previous_am_pairs.clear();
    active_recovery_pairs.clear();
    for (int i = 0; i < actin.n; ++i) {
        for (int j = i + 1; j < actin.n; ++j) {
            if (actin_actin_bonds[i][j] == 1) {
                active_aa_pairs.emplace_back(i, j);
            }
            if (actin_actin_bonds_prev[i][j] == 1) {
                previous_aa_pairs.emplace_back(i, j);
            }
            if (actin_recovery_until[i][j] > current_step) {
                active_recovery_pairs.emplace_back(i, j);
            }
        }
        for (int j = 0; j < myosin.n; ++j) {
            if (am_bonds[i][j] == 1) {
                active_am_pairs.emplace_back(i, j);
            }
            if (am_bonds_prev[i][j] == 1) {
                previous_am_pairs.emplace_back(i, j);
            }
        }
    }
}


void Sarcomere::_process_actin_myosin_binding(int& i) {
    int thread_id = omp_get_thread_num();
    auto& local_actinIndicesPerMyosin = actinIndicesPerMyosin_temp[thread_id];
    const auto& myosin_neighbors = actin_neighbors_by_species[i].second;
    // Hoist the string-keyed feature lookups out of the neighbor loop: each
    // `actin["name"]` hashes a std::string and walks a bucket chain.
    // operator[] takes a `const std::string&`, so a string literal built a
    // fresh temporary here on every actin of every step; all but the shortest
    // of these names exceed the small-string buffer and therefore hit the
    // allocator.  Holding the keys in function-scope statics removes that
    // churn without changing the lookup itself.
    static const std::string key_overlap_available{"am_overlap_available"};
    static const std::string key_overlap_axis{"am_overlap_axis_weighted"};
    static const std::string key_overlap_productive{
        "am_overlap_productive_weighted"};
    static const std::string key_overlap_engaged{"am_overlap_engaged"};
    static const std::string key_weak_count{"am_weak_count"};
    static const std::string key_strong_count{"am_strong_count"};
    static const std::string key_partial_binding{"partial_binding_ratio"};
    static const std::string key_myosin_binding{"myosin_binding_ratio"};
    static const std::string key_crosslink_ratio{"crosslink_ratio"};
    auto& am_overlap_available_feature = actin[key_overlap_available];
    auto& am_overlap_axis_weighted_feature = actin[key_overlap_axis];
    auto& am_overlap_productive_weighted_feature =
        actin[key_overlap_productive];
    auto& am_overlap_engaged_feature = actin[key_overlap_engaged];
    auto& am_weak_count_feature = actin[key_weak_count];
    auto& am_strong_count_feature = actin[key_strong_count];
    auto& partial_binding_ratio_feature = actin[key_partial_binding];
    auto& myosin_binding_ratio_feature = actin[key_myosin_binding];
    auto& crosslink_ratio_feature = actin[key_crosslink_ratio];
    double f_load_contrib = 0;
    double f_load_contrib_cb = 0;
    for (int index = 0; index < myosin_neighbors.size(); index++) {
        int j = myosin_neighbors[index];
        const double capture_centerline = myosin.radius + actin_radius +
                                          am_surface_gap + am_cutoff;
        if (opt_gates_enabled &&
            !geometry::segments_may_interact(
                                  actin.left_end[i], actin.right_end[i],
                                  myosin.left_end[j], myosin.right_end[j],
                                  capture_centerline, box, is_periodic)) {
            // Out of range: analyze_am would have returned zeros, every
            // overlap-weighted accumulation below would have added zero, and
            // geometry_allowed would have been false.  Store the same zeros so
            // the downstream KMC sees an unchanged am_interaction.
            am_interaction[i][j] = {0.0, 0.0};
            continue;
        }
        const auto interaction_geometry = geometry::analyze_am(
            actin.left_end[i], actin.right_end[i], myosin.left_end[j], myosin.right_end[j],
            capture_centerline, box, is_periodic, directional);
        am_interaction[i][j] = {
            interaction_geometry.myosin_binding_ratio,
            interaction_geometry.partial_binding_ratio};
        const double partial_ratio = am_interaction[i][j].partial_binding_ratio;
        const double binding_ratio = am_interaction[i][j].myosin_binding_ratio;
        const double overlap_length = std::max(0.0, binding_ratio) * actin.length;
        const double actin_myosin_dot = std::clamp(
            actin.direction[i].dot(myosin.direction[j]), -1.0, 1.0);
        vec center_displacement = actin.center[i] - myosin.center[j];
        center_displacement.pbc_wrap(box, is_periodic);
        const double contact_axial = center_displacement.dot(myosin.direction[j]);
        const int half_sign = actn2_load_kinetics::myosin_half(contact_axial);
        const double axis_alignment =
            actn2_load_kinetics::nematic_axis_alignment(actin_myosin_dot);
        const double productive_alignment =
            actn2_load_kinetics::productive_polarity_score(
                actin_myosin_dot, half_sign);
        am_overlap_available_feature[i] += overlap_length;
        am_overlap_axis_weighted_feature[i] +=
            overlap_length * axis_alignment;
        am_overlap_productive_weighted_feature[i] +=
            overlap_length * productive_alignment;
        const bool geometry_allowed = binding_ratio > EPS &&
            (partial_ratio > EPS || !directional);

        if (dynamic_am) {
            // KMC is performed once per conserved myosin head after every
            // actin-myosin geometry has been evaluated. Processing it here
            // would allow the first actin in iteration order to monopolize a
            // shared physical head.
            continue;
        } else if (!geometry_allowed) {
            continue;
        }

        if (geometry_allowed) {
                    if (actin_crosslink_ratio[i] > interaction_geometry.crosslinkable_ratio) {
                        actin_crosslink_ratio[i] = interaction_geometry.crosslinkable_ratio;
                        actin_crosslink_start[i] = interaction_geometry.crosslinkable_start;
                        actin_crosslink_end[i] = interaction_geometry.crosslinkable_end;
                        if (actin_crosslink_ratio[i] < EPS) {
                            const auto& actin_neighbors = actin_neighbors_by_species[i].first;
                            bool prev_catch_bonded = false;
                            for (int idx = 0; idx < actin_neighbors.size(); idx++){
                                if (actin_actin_status_prev[i][actin_neighbors[idx]] > 1){
                                    prev_catch_bonded = true;
                                    break;
                                }
                            }
                            if (prev_catch_bonded) {
                                double angle = std::acos(std::abs(actin.direction[i].dot(myosin.direction[j]))) * 180.0 / M_PI;
                                vec repulsive_force = compute_actin_myosin_repulsion(
                                    actin,
                                    myosin,
                                    i,
                                    j,
                                    box,
                                    am_cutoff*1.1,
                                    actin_steric_k);
                                double force_magnitude = repulsive_force.norm();
                            }
                        }
                    }
        }

        local_actinIndicesPerMyosin.addConnection(j, i);
        myosinIndicesPerActin.addConnection(i, j);
        am_overlap_engaged_feature[i] += overlap_length;
        n_myosins_per_actin[i]++;
        const bool strong = !dynamic_am ||
            am_state[i][j] == static_cast<BondFlag>(am_motor_kinetics::State::strong);
        if (dynamic_am) {
            int weak_heads = 0;
            int strong_heads = 0;
            for (const auto& half_heads : am_heads_by_myosin[static_cast<size_t>(j)]) {
                for (const auto& head : half_heads) {
                    if (head.bound_actin != i) continue;
                    if (head.state == static_cast<BondFlag>(am_motor_kinetics::State::weak)) {
                        ++weak_heads;
                    } else if (head.state == static_cast<BondFlag>(am_motor_kinetics::State::strong)) {
                        ++strong_heads;
                    }
                }
            }
            am_weak_count_feature[i] += static_cast<double>(weak_heads);
            am_strong_count_feature[i] += static_cast<double>(strong_heads);
        }
        double abs_cos = std::abs(actin.direction[i].dot(myosin.direction[j]));
        if (strong && abs_cos > actin_basic_tension[i]) {
            actin_basic_tension[i] = abs_cos;
        }
        if (partial_binding_ratio_feature[i] < partial_ratio) {
            partial_binding_ratio_feature[i] = partial_ratio;
        }
        if (myosin_binding_ratio_feature[i] < binding_ratio) {
            myosin_binding_ratio_feature[i] = binding_ratio;
        }

        am_bonds[i][j] = 1;
        active_am_pairs_temp[thread_id].emplace_back(i, j);
        if (strong) {
            double partial = directional ? partial_ratio : binding_ratio;
            double contrib = std::min(
                partial / geometry::am_head_bearing_fraction, 1.0);
            double contrib_cb = contrib * abs_cos;
            if (contrib > f_load_contrib) f_load_contrib = contrib;
            if (contrib_cb > f_load_contrib_cb) f_load_contrib_cb = contrib_cb;
        }
    }
    double denom = (1 - std::exp(-2));
    //double load_velocity = (1 - std::exp(-2 * f_load_contrib)) / denom;
    double load_velocity = f_load_contrib;
    double load_cb = (1 - std::exp(-2 * f_load_contrib_cb)) / denom;
    if (load_velocity > 1.0) load_velocity = 1.0;
    if (load_cb > 1.0) load_cb = 1.0;
    actin.f_load[i] = load_velocity;
    if (actin_f_load_cb) {
        (*actin_f_load_cb)[i] = load_cb;
    }
    crosslink_ratio_feature[i] = actin_crosslink_ratio[i];

}

void Sarcomere::_process_dynamic_am_heads() {
    using Candidate = AMHeadCandidate;

    struct PairAccountingScratch {
        std::vector<unsigned int> generation;
        std::vector<int> capacity;
        std::vector<int> occupied;
        std::vector<int> strong;
        std::vector<int> candidate_index;
        std::vector<int> touched;
        unsigned int current_generation = 0;

        void begin(size_t count) {
            if (generation.size() < count) {
                generation.resize(count, 0);
                capacity.resize(count);
                occupied.resize(count);
                strong.resize(count);
                candidate_index.resize(count);
            }
            ++current_generation;
            if (current_generation == 0) {
                std::fill(generation.begin(), generation.end(), 0);
                current_generation = 1;
            }
            touched.clear();
        }

        void touch(int actin_index) {
            const size_t index = static_cast<size_t>(actin_index);
            if (generation[index] == current_generation) return;
            generation[index] = current_generation;
            capacity[index] = -1;
            occupied[index] = 0;
            strong[index] = 0;
            candidate_index[index] = -1;
            touched.push_back(actin_index);
        }
    };
    static thread_local PairAccountingScratch pair_scratch;

    const int thread_id = omp_in_parallel() ? omp_get_thread_num() : 0;
    auto& local_actin_indices = actinIndicesPerMyosin_temp[thread_id];
    auto& candidates = am_head_candidates_;
    // Hoist string-keyed feature lookups out of the per-myosin / per-actin
    // loops below.  The half-suffixed myosin features previously built a fresh
    // std::string on every myosin half every step.
    auto& am_off_count_minus_feature = myosin["am_off_count_minus"];
    auto& am_off_count_plus_feature = myosin["am_off_count_plus"];
    auto& am_on_count_minus_feature = myosin["am_on_count_minus"];
    auto& am_on_count_plus_feature = myosin["am_on_count_plus"];
    auto& am_tension_minus_feature = myosin["am_tension_minus"];
    auto& am_tension_plus_feature = myosin["am_tension_plus"];
    auto& am_overlap_engaged_feature = actin["am_overlap_engaged"];
    auto& am_weak_count_feature = actin["am_weak_count"];
    auto& am_strong_count_feature = actin["am_strong_count"];
    auto& partial_binding_ratio_feature = actin["partial_binding_ratio"];
    auto& myosin_binding_ratio_feature = actin["myosin_binding_ratio"];
    auto& crosslink_ratio_feature = actin["crosslink_ratio"];

    #pragma omp single
    {
    candidates.resize(static_cast<size_t>(myosin.n));
    for (auto& myosin_candidates : candidates) {
        // Preserve capacity because local filament neighborhoods usually
        // change slowly from one timestep to the next.
        myosin_candidates.clear();
    }
    }

    // Candidate geometry is independent for each myosin.  Assigning complete
    // myosin candidate lists to threads avoids a serial all-actin scan and
    // requires no locks or cross-thread merge.  Sorting by actin index below
    // preserves the exact ordering of the former actin-major serial loop.
    #pragma omp for schedule(dynamic, 1)
    for (int j = 0; j < myosin.n; ++j) {
        auto& myosin_candidates = candidates[static_cast<size_t>(j)];
        for (const auto& [neighbor_index, neighbor_type] :
             neighbor_list.get_neighbors(j + actin.n)) {
            if (neighbor_type != ParticleType::Actin) continue;
            const int i = neighbor_index;
            const double binding_ratio = am_interaction[i][j].myosin_binding_ratio;
            const double partial_ratio = am_interaction[i][j].partial_binding_ratio;
            const double distance = geometry::segment_segment_distance(
                actin.left_end[i], actin.right_end[i],
                myosin.left_end[j], myosin.right_end[j], box, is_periodic);
            const double association_centerline = myosin.radius + actin_radius +
                am_surface_gap + am_cutoff;
            if (distance > association_centerline) continue;
            const double actin_dot_myosin = std::clamp(
                actin.direction[i].dot(myosin.direction[j]), -1.0, 1.0);
            vec center_displacement = actin.center[i] - myosin.center[j];
            center_displacement.pbc_wrap(box, is_periodic);
            const int productive_half_sign = directional
                ? actn2_load_kinetics::productive_myosin_half(
                      actin_dot_myosin)
                : actn2_load_kinetics::myosin_half(
                      center_displacement.dot(myosin.direction[j]));

            // Compute axial overlap with both physical head-bearing halves.
            // The polarity-compatible half can power stroke. If the actin has
            // crossed the M-line far enough to reach the opposite head-bearing
            // half, those heads may form weak nonproductive encounters that
            // consume real myosin heads and actin sites but never become
            // strong motors. The central bare zone itself remains force-free.
            vec axis_m = myosin.direction[j];
            axis_m.normalize();
            vec left_disp = actin.left_end[i] - myosin.center[j];
            vec right_disp = actin.right_end[i] - myosin.center[j];
            left_disp.pbc_wrap(box, is_periodic);
            right_disp.pbc_wrap(box, is_periodic);
            const double actin_axial_min = std::min(
                left_disp.dot(axis_m), right_disp.dot(axis_m));
            const double actin_axial_max = std::max(
                left_disp.dot(axis_m), right_disp.dot(axis_m));
            const double bare_half =
                0.5 * geometry::am_bare_zone_fraction * myosin.length;
            for (int half_index = 0; half_index < 2; ++half_index) {
                const int half_sign = half_index == 0 ? -1 : 1;
                const double head_axial_min = half_sign < 0
                    ? -0.5 * myosin.length : bare_half;
                const double head_axial_max = half_sign < 0
                    ? -bare_half : 0.5 * myosin.length;
                const double axial_overlap = std::max(
                    0.0, std::min(actin_axial_max, head_axial_max) -
                         std::max(actin_axial_min, head_axial_min));

                bool existing_attachment = false;
                for (const auto& head :
                     am_heads_by_myosin[static_cast<size_t>(j)]
                                        [static_cast<size_t>(half_index)]) {
                    if (head.bound_actin == i &&
                        head.state != static_cast<BondFlag>(
                            am_motor_kinetics::State::detached) &&
                        head.state != static_cast<BondFlag>(
                            am_motor_kinetics::State::off_srx)) {
                        existing_attachment = true;
                        break;
                    }
                }
                const bool association_geometry =
                    binding_ratio > EPS && axial_overlap > EPS;
                const bool retention_geometry =
                    existing_attachment && axial_overlap > EPS;
                if (!association_geometry && !retention_geometry) continue;
                // Reaching here guarantees axial_overlap > EPS for this
                // half, which is exactly what _am_actin_reaches_half() tests:
                // the actin has come past the bare zone into this half's head
                // band.  So under pull-through the opposite half is productive.
                const bool crossover_productive =
                    am_crossover_pull_through &&
                    half_sign != productive_half_sign;
                myosin_candidates.push_back({
                    i,
                    half_index,
                    !directional || half_sign == productive_half_sign ||
                        crossover_productive,
                    axial_overlap,
                    binding_ratio,
                    partial_ratio});
            }
        }
        std::sort(
            myosin_candidates.begin(), myosin_candidates.end(),
            [](const AMHeadCandidate& left, const AMHeadCandidate& right) {
                return left.actin_index < right.actin_index;
            });
    }

    #pragma omp single
    {
    // Clean legacy/resumed duplicates before chemistry, then expose one
    // shared site table to all myosins for this timestep's associations.
    _enforce_global_am_site_ownership();
    _rebuild_am_site_claims();
    }

    #pragma omp for schedule(dynamic, 1)
    for (int j = 0; j < myosin.n; ++j) {
        auto& myosin_candidates = candidates[static_cast<size_t>(j)];
        for (int half_index = 0; half_index < 2; ++half_index) {
            auto& heads = am_heads_by_myosin[static_cast<size_t>(j)]
                                              [static_cast<size_t>(half_index)];
            // A half-filament owns one conserved head pool, but a single
            // neighboring actin may occupy only its overlap-scaled share.
            // Strong heads reserve slots first. If sliding lowers capacity,
            // discard only the most strained excess weak heads; strong heads
            // retain their normal load-dependent KMC lifetime.
            pair_scratch.begin(static_cast<size_t>(actin.n));
            const double head_bearing_length =
                geometry::am_head_bearing_fraction * myosin.length;
            for (int candidate_number = 0;
                 candidate_number < static_cast<int>(myosin_candidates.size());
                 ++candidate_number) {
                const auto& candidate =
                    myosin_candidates[static_cast<size_t>(candidate_number)];
                if (candidate.half_index != half_index) continue;
                pair_scratch.touch(candidate.actin_index);
                const int capacity =
                    am_motor_kinetics::overlap_limited_neighbor_capacity(
                        am_heads_per_half, am_max_heads_per_actin_half,
                        candidate.axial_overlap, head_bearing_length);
                const size_t actin_index =
                    static_cast<size_t>(candidate.actin_index);
                pair_scratch.capacity[actin_index] =
                    std::max(pair_scratch.capacity[actin_index], capacity);
                pair_scratch.candidate_index[actin_index] = candidate_number;
            }
            for (const auto& head : heads) {
                if (head.bound_actin < 0) continue;
                const auto state =
                    static_cast<am_motor_kinetics::State>(head.state);
                if (state != am_motor_kinetics::State::weak &&
                    state != am_motor_kinetics::State::strong) continue;
                pair_scratch.touch(head.bound_actin);
                ++pair_scratch.occupied[static_cast<size_t>(
                    head.bound_actin)];
                if (state == am_motor_kinetics::State::strong) {
                    ++pair_scratch.strong[static_cast<size_t>(
                        head.bound_actin)];
                }
            }
            for (int i : pair_scratch.touched) {
                const size_t actin_index = static_cast<size_t>(i);
                const int capacity = pair_scratch.capacity[actin_index];
                if (capacity < 0) continue;
                const int weak_allowance =
                    am_motor_kinetics::weak_slots_after_strong(
                        capacity, pair_scratch.strong[actin_index]);
                int excess_weak = std::max(
                    0, pair_scratch.occupied[actin_index] -
                       pair_scratch.strong[actin_index] - weak_allowance);
                if (excess_weak == 0) continue;
                std::vector<std::pair<double, int>> weak_heads;
                for (int h = 0; h < static_cast<int>(heads.size()); ++h) {
                    const auto& head = heads[static_cast<size_t>(h)];
                    if (head.bound_actin != i ||
                        head.state != static_cast<BondFlag>(
                            am_motor_kinetics::State::weak)) continue;
                    weak_heads.emplace_back(
                        _am_head_strain_energy(
                            i, j, _am_head_id(half_index, h),
                            am_motor_kinetics::State::weak), h);
                }
                std::sort(weak_heads.begin(), weak_heads.end(),
                          [](const auto& first, const auto& second) {
                              if (first.first != second.first) {
                                  return first.first > second.first;
                              }
                              return first.second > second.second;
                          });
                for (const auto& [energy, h] : weak_heads) {
                    (void)energy;
                    if (excess_weak-- <= 0) break;
                    auto& head = heads[static_cast<size_t>(h)];
                    _record_am_transition(
                        i, j, _am_head_id(half_index, h),
                        am_motor_kinetics::State::weak,
                        am_motor_kinetics::State::detached,
                        AMTransitionReason::capacity_limit);
                    _release_am_site_claim(i, head.target_axial);
                    head.state = static_cast<BondFlag>(
                        am_motor_kinetics::State::detached);
                    head.bound_actin = -1;
                    --pair_scratch.occupied[actin_index];
                }
            }
            double half_tension = 0.0;
            if (am_thick_force_clamp >= 0.0) {
                half_tension = am_thick_force_clamp;
            } else {
                for (int h = 0; h < static_cast<int>(heads.size()); ++h) {
                    const auto& head = heads[static_cast<size_t>(h)];
                    if (head.bound_actin < 0 ||
                        head.state != static_cast<BondFlag>(
                            am_motor_kinetics::State::strong)) continue;
                    half_tension += _am_head_active_force(
                        head.bound_actin, j, _am_head_id(half_index, h));
                }
            }
            am_half_tension[static_cast<size_t>(j)]
                           [static_cast<size_t>(half_index)] = half_tension;
            std::vector<double> attachment_rates;
            std::vector<const Candidate*> eligible_candidates;
            attachment_rates.reserve(myosin_candidates.size());
            eligible_candidates.reserve(myosin_candidates.size());
            for (int h = 0; h < static_cast<int>(heads.size()); ++h) {
                auto& head = heads[static_cast<size_t>(h)];
                const int head_id = _am_head_id(half_index, h);
                counter_rng::CounterRng head_rng(
                    static_cast<std::uint64_t>(initial_seed),
                    static_cast<std::uint64_t>(current_step),
                    counter_rng::Stream::am_head,
                    static_cast<std::uint64_t>(j),
                    (static_cast<std::uint64_t>(half_index) << 32) |
                        static_cast<std::uint32_t>(h));
                auto state = static_cast<am_motor_kinetics::State>(head.state);
                bool just_recruited = false;

                if (!am_thick_feedback_enabled &&
                    state == am_motor_kinetics::State::off_srx) {
                    head.state = static_cast<BondFlag>(
                        am_motor_kinetics::State::detached);
                    state = am_motor_kinetics::State::detached;
                } else if (am_thick_feedback_enabled &&
                           state == am_motor_kinetics::State::off_srx) {
                    const double recruit_rate =
                        am_motor_kinetics::thick_recruitment_rate(
                            am_thick_recruit_rate, half_tension,
                            am_thick_tension_scale);
                    if (head_rng.uniform() <
                        am_motor_kinetics::transition_probability(
                            recruit_rate, dt)) {
                        _record_am_thick_transition(
                            j, half_index, head_id, state,
                            am_motor_kinetics::State::detached, half_tension);
                        head.state = static_cast<BondFlag>(
                            am_motor_kinetics::State::detached);
                        state = am_motor_kinetics::State::detached;
                        just_recruited = true;
                    } else {
                        continue;
                    }
                }

                const Candidate* bound_candidate = nullptr;
                if (head.bound_actin >= 0) {
                    pair_scratch.touch(head.bound_actin);
                    const int candidate_number =
                        pair_scratch.candidate_index[static_cast<size_t>(
                            head.bound_actin)];
                    if (candidate_number >= 0) {
                        bound_candidate = &myosin_candidates[
                            static_cast<size_t>(candidate_number)];
                    }
                    // Candidate lists govern association only. A material
                    // cross-bridge is evaluated directly until it detaches or
                    // truly loses productive overlap; neighbor-cell changes
                    // cannot erase its chemical state.
                    if (!_am_head_retention_reachable(
                            head.bound_actin, j, head) ||
                        state == am_motor_kinetics::State::detached) {
                        const int old_actin = head.bound_actin;
                        const bool counted =
                            state == am_motor_kinetics::State::weak ||
                            state == am_motor_kinetics::State::strong;
                        if (state != am_motor_kinetics::State::detached) {
                            _record_am_forced_rupture(
                                head.bound_actin, j, head_id, state,
                                AMTransitionReason::geometric_rupture,
                                bound_candidate == nullptr);
                            _record_am_transition(
                                head.bound_actin, j, head_id, state,
                                am_motor_kinetics::State::detached,
                                AMTransitionReason::geometric_rupture);
                        }
                        _release_am_site_claim(
                            old_actin, head.target_axial);
                        head.state = static_cast<BondFlag>(
                            am_motor_kinetics::State::detached);
                        head.bound_actin = -1;
                        state = am_motor_kinetics::State::detached;
                        if (counted && old_actin >= 0) {
                            --pair_scratch.occupied[static_cast<size_t>(
                                old_actin)];
                        }
                    }
                }

                if (head.bound_actin >= 0 &&
                    state == am_motor_kinetics::State::weak) {
                    const bool productive_attachment =
                        _am_head_is_productive(head.bound_actin, j, head);
                    const double k_ps = productive_attachment
                        ? am_motor_kinetics::load_slowed_rate(
                              am_k_powerstroke, 0.0, am_force_scale) *
                              std::clamp(am_thin_activation[static_cast<size_t>(
                                  head.bound_actin)], 0.0, 1.0)
                        : 0.0;
                    const double k_weak_det = _am_head_detachment_rate(
                        head.bound_actin, j, head_id, state, am_k_weak_off);
                    const double total_rate = k_weak_det + k_ps;
                    if (head_rng.uniform() <
                        am_motor_kinetics::transition_probability(total_rate, dt)) {
                        const double choice = head_rng.uniform();
                        if (total_rate <= 0.0 ||
                            choice < k_weak_det / total_rate) {
                            const int old_actin = head.bound_actin;
                            _record_am_transition(
                                old_actin, j, head_id, state,
                                am_motor_kinetics::State::detached,
                                AMTransitionReason::weak_detachment);
                            _release_am_site_claim(
                                old_actin, head.target_axial);
                            head.state = static_cast<BondFlag>(
                                am_motor_kinetics::State::detached);
                            head.bound_actin = -1;
                            state = am_motor_kinetics::State::detached;
                            --pair_scratch.occupied[static_cast<size_t>(
                                old_actin)];
                        } else {
                            // Enter the force-producing state at zero stroke
                            // progress. The internal coordinate advances over
                            // subsequent steps according to force velocity,
                            // rather than receiving the full stroke instantly.
                            head.stroke_progress = 0.0;
                            _record_am_transition(
                                head.bound_actin, j, head_id, state,
                                am_motor_kinetics::State::strong,
                                AMTransitionReason::powerstroke);
                            head.state = static_cast<BondFlag>(
                                am_motor_kinetics::State::strong);
                            state = am_motor_kinetics::State::strong;
                        }
                    }
                } else if (head.bound_actin >= 0 &&
                           state == am_motor_kinetics::State::strong) {
                    const double force_before_step = _am_head_active_force(
                        head.bound_actin, j, head_id);
                    const double remaining_stroke = std::max(
                        0.0, am_power_stroke - head.stroke_progress);
                    // Legacy sliding throttle, applied where the legacy model
                    // applied it: to the sliding velocity itself.  The prescribed
                    // model used v = (1 - f_load(overlap)) * v0 with
                    // f_load = (1-exp(-2r))/(1-exp(-2)), so sliding halted once the
                    // engaged AM overlap reached its target.  A bound head keeps its
                    // elastic hold; it simply stops stroking further.
                    double overlap_throttle = 1.0;
                    if (am_overlap_target_ratio > 0.0) {
                        const double engaged =
                            am_interaction[head.bound_actin][j].partial_binding_ratio;
                        const double r = std::clamp(
                            engaged / am_overlap_target_ratio, 0.0, 1.0);
                        const double f_load =
                            (1.0 - std::exp(-2.0 * r)) / (1.0 - std::exp(-2.0));
                        overlap_throttle = std::clamp(1.0 - f_load, 0.0, 1.0);
                    }
                    if (!am_stroke_detailed_balance) {
                        const double stroke_velocity = overlap_throttle *
                            am_motor_kinetics::working_stroke_velocity(
                                am_unloaded_velocity, force_before_step,
                                am_stall_force);
                        const double stroke_increment = std::min(
                            remaining_stroke, stroke_velocity * dt);
                        head.stroke_progress += stroke_increment;
                    } else {
                        // Discrete working stroke with detailed balance.  The
                        // rest coordinate is two-valued (0 or am_power_stroke);
                        // forward and reverse rates depend on the strain the
                        // transition would create, so the force-velocity
                        // relation emerges rather than being imposed.
                        const double beta_value = thermal_beta > 0.0
                            ? thermal_beta
                            : rod_mobility::beta(temperature_kelvin);
                        const double d = am_power_stroke;
                        const double z = head.stroke_progress -
                            _am_productive_strain(head.bound_actin, j, head_id);
                        const double dG_mech =
                            am_motor_kinetics::stroke_mechanical_free_energy(
                                am_motor_k, d, z);
                        const double dG_chem = am_stroke_free_energy >= 0.0
                            ? am_stroke_free_energy
                            : d * (am_stall_force - 0.5 * am_motor_k * d);
                        const auto rates =
                            am_motor_kinetics::stroke_transition_rates(
                                am_stroke_base_rate * overlap_throttle,
                                beta_value, dG_mech, dG_chem, am_stroke_theta,
                                static_cast<double>(
                                    am_strain_detachment_max_exponent));
                        const bool pre_stroke = head.stroke_progress < 0.5 * d;
                        const double rate = pre_stroke ? rates.forward
                                                       : rates.reverse;
                        if (rate > 0.0 &&
                            head_rng.uniform() < 1.0 - std::exp(-rate * dt)) {
                            head.stroke_progress = pre_stroke ? d : 0.0;
                        }
                    }
                    const double chemical_k_det =
                        am_motor_kinetics::load_slowed_rate_with_floor(
                        am_k_strong_off,
                        _am_head_active_force(head.bound_actin, j, head_id),
                        am_force_scale, am_strong_off_floor_fraction);
                    const double k_det = _am_head_detachment_rate(
                        head.bound_actin, j, head_id, state, chemical_k_det);
                    if (head_rng.uniform() <
                        am_motor_kinetics::transition_probability(k_det, dt)) {
                        const int old_actin = head.bound_actin;
                        _record_am_transition(
                            old_actin, j, head_id, state,
                            am_motor_kinetics::State::detached,
                            AMTransitionReason::strong_detachment);
                        _release_am_site_claim(
                            old_actin, head.target_axial);
                        head.state = static_cast<BondFlag>(
                            am_motor_kinetics::State::detached);
                        head.bound_actin = -1;
                        state = am_motor_kinetics::State::detached;
                        --pair_scratch.occupied[static_cast<size_t>(
                            old_actin)];
                    }
                }

                if (head.bound_actin >= 0 ||
                    state != am_motor_kinetics::State::detached) continue;

                if (am_thick_feedback_enabled && !just_recruited) {
                    const double fold_rate =
                        am_motor_kinetics::thick_folding_rate(
                            am_thick_fold_rate, half_tension,
                            am_thick_tension_scale);
                    if (head_rng.uniform() <
                        am_motor_kinetics::transition_probability(
                            fold_rate, dt)) {
                        _record_am_thick_transition(
                            j, half_index, head_id, state,
                            am_motor_kinetics::State::off_srx, half_tension);
                        head.state = static_cast<BondFlag>(
                            am_motor_kinetics::State::off_srx);
                        continue;
                    }
                }

                // Exact rejection thinning: spatial_gate is bounded by one,
                // so the activation-only sum is an upper rate. Most 1-us
                // steps contain no proposal and can skip all expensive head
                // geometry while preserving the continuous-time KMC law.
                double upper_attach_rate = 0.0;
                for (const auto& candidate : myosin_candidates) {
                    if (candidate.half_index != half_index) continue;
                    const size_t actin_index =
                        static_cast<size_t>(candidate.actin_index);
                    const int capacity = pair_scratch.capacity[actin_index];
                    if (capacity >= 0 &&
                        pair_scratch.occupied[actin_index] >= capacity) {
                        continue;
                    }
                    upper_attach_rate += am_k_on *
                        am_motor_kinetics::weak_attachment_gate(
                            am_thin_activation[actin_index],
                            am_weak_binding_floor);
                }
                if (!(upper_attach_rate > 0.0)) continue;
                double proposal_remaining = dt;
                const double first_wait = -std::log(
                    head_rng.uniform_pos()) /
                    upper_attach_rate;
                if (first_wait >= proposal_remaining) continue;
                proposal_remaining -= first_wait;

                attachment_rates.clear();
                eligible_candidates.clear();
                double total_attach_rate = 0.0;
                for (const auto& candidate : myosin_candidates) {
                    if (candidate.half_index != half_index) continue;
                    const size_t actin_index =
                        static_cast<size_t>(candidate.actin_index);
                    const int capacity = pair_scratch.capacity[actin_index];
                    if (capacity >= 0 &&
                        pair_scratch.occupied[actin_index] >= capacity) continue;
                    const double candidate_azimuth =
                        am_binding_azimuth_selection && am_spatial_heads
                        ? select_am_binding_azimuth(
                              actin, myosin, candidate.actin_index, j, box,
                              head.anchor_axial, head.anchor_azimuth,
                              am_binding_azimuth_half_width,
                              myosin_material_e1[static_cast<size_t>(j)],
                              myosin_material_e2[static_cast<size_t>(j)])
                        : head.anchor_azimuth;
                    const double spatial_gate = _am_head_association_gate(
                        candidate.actin_index, j, head, candidate_azimuth);
                    const double rate = am_k_on *
                                        am_motor_kinetics::weak_attachment_gate(
                                            am_thin_activation[static_cast<size_t>(
                                                candidate.actin_index)],
                                            am_weak_binding_floor) *
                                        spatial_gate;
                    if (!(rate > 0.0)) continue;
                    eligible_candidates.push_back(&candidate);
                    attachment_rates.push_back(rate);
                    total_attach_rate += rate;
                }
                if (!(total_attach_rate > 0.0)) continue;
                const double acceptance = std::min(
                    1.0, total_attach_rate / upper_attach_rate);
                bool accepted = false;
                while (true) {
                    if (head_rng.uniform() < acceptance) {
                        accepted = true;
                        break;
                    }
                    const double wait = -std::log(
                        head_rng.uniform_pos()) /
                        upper_attach_rate;
                    if (wait >= proposal_remaining) break;
                    proposal_remaining -= wait;
                }
                if (!accepted) continue;

                double draw = head_rng.uniform() * total_attach_rate;
                size_t chosen = attachment_rates.size() - 1;
                for (size_t k = 0; k < attachment_rates.size(); ++k) {
                    draw -= attachment_rates[k];
                    if (draw <= 0.0) {
                        chosen = k;
                        break;
                    }
                }
                const int target_actin =
                    eligible_candidates[chosen]->actin_index;
                head.binding_azimuth =
                    am_binding_azimuth_selection && am_spatial_heads
                    ? select_am_binding_azimuth(
                          actin, myosin, target_actin, j, box,
                          head.anchor_axial, head.anchor_azimuth,
                          am_binding_azimuth_half_width,
                          myosin_material_e1[static_cast<size_t>(j)],
                          myosin_material_e2[static_cast<size_t>(j)])
                    : head.anchor_azimuth;
                const double target_site = _claim_free_am_site_coordinate(
                    target_actin, j, head, head.binding_azimuth);
                if (!std::isfinite(target_site)) continue;
                head.bound_actin = target_actin;
                head.stroke_progress = 0.0;
                // This is a material coordinate on actin, fixed until the
                // head detaches.  Sliding now stretches the axial head spring
                // instead of silently moving the attachment point.
                head.target_axial = target_site;
                _record_am_transition(
                    target_actin, j, head_id,
                    am_motor_kinetics::State::detached,
                    am_motor_kinetics::State::weak,
                    AMTransitionReason::attachment);
                head.state = static_cast<BondFlag>(am_motor_kinetics::State::weak);
                ++pair_scratch.occupied[static_cast<size_t>(target_actin)];
            }

            int off_count = 0;
            int on_count = 0;
            for (const auto& head : heads) {
                if (head.state == static_cast<BondFlag>(
                        am_motor_kinetics::State::off_srx)) {
                    ++off_count;
                } else {
                    ++on_count;
                }
            }
            (half_index == 0 ? am_off_count_minus_feature
                             : am_off_count_plus_feature)[j] = off_count;
            (half_index == 0 ? am_on_count_minus_feature
                             : am_on_count_plus_feature)[j] = on_count;
            (half_index == 0 ? am_tension_minus_feature
                             : am_tension_plus_feature)[j] = half_tension;
        }

        // Reused across myosins so the per-myosin summary does not allocate;
        // the loop is a worksharing construct, so the buffer is thread_local
        // like pair_scratch above.
        static thread_local std::vector<int> bound_actins;
        bound_actins.clear();
        for (const auto& half_heads :
             am_heads_by_myosin[static_cast<size_t>(j)]) {
            for (const auto& head : half_heads) {
                if (head.bound_actin >= 0 &&
                    head.state != static_cast<BondFlag>(
                        am_motor_kinetics::State::detached)) {
                    bound_actins.push_back(head.bound_actin);
                }
            }
        }
        std::sort(bound_actins.begin(), bound_actins.end());
        bound_actins.erase(
            std::unique(bound_actins.begin(), bound_actins.end()),
            bound_actins.end());

        for (int i : bound_actins) {
            _sync_am_pair_summary(i, j);
            if (am_state[i][j] == static_cast<BondFlag>(
                    am_motor_kinetics::State::detached)) continue;
            am_bonds[i][j] = 1;
            am_seen[i][j] = 1;
            active_am_pairs_temp[thread_id].emplace_back(i, j);
            local_actin_indices.addConnection(j, i);
            // Only the per-actin accumulators below need serializing: several
            // myosins may reduce into the same actin concurrently.  The head
            // scan and analyze_am read nothing the lock protects, so they are
            // evaluated before it is taken -- holding a single global lock
            // across a 2x141-head scan plus a segment-geometry solve made this
            // region the dominant Amdahl term past ~8 threads.
            const double binding_ratio = am_interaction[i][j].myosin_binding_ratio;
            const double partial_ratio = am_interaction[i][j].partial_binding_ratio;
            int weak_heads = 0;
            int strong_heads = 0;
            for (const auto& half_heads :
                 am_heads_by_myosin[static_cast<size_t>(j)]) {
                for (const auto& head : half_heads) {
                    if (head.bound_actin != i) continue;
                    if (head.state == static_cast<BondFlag>(
                            am_motor_kinetics::State::weak)) ++weak_heads;
                    if (head.state == static_cast<BondFlag>(
                            am_motor_kinetics::State::strong)) ++strong_heads;
                }
            }
            const double abs_cos = std::abs(
                actin.direction[i].dot(myosin.direction[j]));
            const auto interaction_geometry = geometry::analyze_am(
                actin.left_end[i], actin.right_end[i],
                myosin.left_end[j], myosin.right_end[j],
                // geometry::analyze_am expects a centerline capture distance.
                // Passing only am_cutoff here omitted both filament radii and
                // the preferred surface gap.  In the explicit-head path that
                // could collapse crosslinkable_ratio to zero even while heads
                // were visibly bound, permanently suppressing replacement
                // ACTN2 links after a prebound link detached.
                myosin.radius + actin_radius + am_surface_gap + am_cutoff,
                box, is_periodic, directional);

            #pragma omp critical(dynamic_am_actin_summary)
            {
            myosinIndicesPerActin.addConnection(i, j);

            am_overlap_engaged_feature[i] +=
                std::max(0.0, binding_ratio) * actin.length;
            ++n_myosins_per_actin[i];

            am_weak_count_feature[i] += static_cast<double>(weak_heads);
            am_strong_count_feature[i] += static_cast<double>(strong_heads);

            if (strong_heads > 0 && abs_cos > actin_basic_tension[i]) {
                actin_basic_tension[i] = abs_cos;
            }
            partial_binding_ratio_feature[i] = std::max(
                partial_binding_ratio_feature[i], partial_ratio);
            myosin_binding_ratio_feature[i] = std::max(
                myosin_binding_ratio_feature[i], binding_ratio);

            if (actin_crosslink_ratio[i] >
                interaction_geometry.crosslinkable_ratio) {
                actin_crosslink_ratio[i] =
                    interaction_geometry.crosslinkable_ratio;
                actin_crosslink_start[i] =
                    interaction_geometry.crosslinkable_start;
                actin_crosslink_end[i] =
                    interaction_geometry.crosslinkable_end;
            }

            if (strong_heads > 0) {
                const double partial = directional ? partial_ratio : binding_ratio;
                const double contribution = std::min(
                    1.0, partial / geometry::am_head_bearing_fraction);
                actin.f_load[i] = std::max(actin.f_load[i], contribution);
                if (actin_f_load_cb) {
                    const double contribution_cb = contribution * abs_cos;
                    const double load_cb = (1.0 - std::exp(-2.0 * contribution_cb)) /
                                           (1.0 - std::exp(-2.0));
                    (*actin_f_load_cb)[i] = std::max(
                        (*actin_f_load_cb)[i], std::min(1.0, load_cb));
                }
            }
            }
        }
    }

    #pragma omp single
    {
    for (int i = 0; i < actin.n; ++i) {
        crosslink_ratio_feature[i] = actin_crosslink_ratio[i];
    }
    }
}

void Sarcomere::_update_thin_filament_activation() {
    // Called from inside the `omp parallel` region of update_system().  Each
    // iteration below writes only am_thin_activation[i] and the per-i entries
    // of two actin features, and reads everything else (am_interaction,
    // am_heads_by_myosin, neighbor lists, filament geometry) read-only, so the
    // loop is an orphaned worksharing construct rather than a `single` block.
    // The resize guard is the one shared write and stays serialized; the
    // implicit barrier of `single` keeps it ordered before the loop.
    #pragma omp single
    {
    if (am_thin_activation.size() != static_cast<size_t>(actin.n)) {
        am_thin_activation.assign(static_cast<size_t>(actin.n),
                                  am_activation_floor);
    }
    }
    auto& am_calcium_feature = actin["am_calcium"];
    auto& am_q_feature = actin["am_q"];
    #pragma omp for schedule(static)
    for (int i = 0; i < actin.n; ++i) {
        int strong_heads = 0;
        int accessible_heads = 0;
        for (int j : actin_neighbors_by_species[static_cast<size_t>(i)].second) {
            const double binding_ratio = am_interaction[i][j].myosin_binding_ratio;
            const double partial_ratio = am_interaction[i][j].partial_binding_ratio;
            const bool association_geometry =
                binding_ratio > EPS && (partial_ratio > EPS || !directional);
            const int half_sign = directional
                ? actn2_load_kinetics::productive_myosin_half(
                      actin.direction[i].dot(myosin.direction[j]))
                : actn2_load_kinetics::myosin_half([&]() {
                      vec displacement = actin.center[i] - myosin.center[j];
                      displacement.pbc_wrap(box, is_periodic);
                      return displacement.dot(myosin.direction[j]);
                  }());
            // Under pull-through the far half also carries productive heads
            // once the actin has reached past the bare zone into it, so the
            // activation feedback must count both halves or it will disagree
            // with the binding path that actually created those heads.
            const int candidate_halves[2] = {half_sign, -half_sign};
            const int half_count =
                (directional && am_crossover_pull_through &&
                 _am_actin_reaches_half(i, j, -half_sign)) ? 2 : 1;
            for (int hs_slot = 0; hs_slot < half_count; ++hs_slot) {
            const int half_index = candidate_halves[hs_slot] < 0 ? 0 : 1;
            const auto& heads =
                am_heads_by_myosin[static_cast<size_t>(j)]
                                   [static_cast<size_t>(half_index)];
            int retained_heads = 0;
            int retained_strong_heads = 0;
            for (const auto& head : heads) {
                if (head.bound_actin != i ||
                    !_am_head_retention_reachable(i, j, head)) continue;
                ++retained_heads;
                if (head.state == static_cast<BondFlag>(
                        am_motor_kinetics::State::strong)) {
                    ++retained_strong_heads;
                }
            }
            int available_slots = 0;
            if (association_geometry) {
                const double overlap = actin.length * std::clamp(
                    directional ? partial_ratio : binding_ratio, 0.0, 1.0);
                available_slots =
                    am_motor_kinetics::overlap_limited_neighbor_capacity(
                        am_heads_per_half, am_max_heads_per_actin_half,
                        overlap,
                        geometry::am_head_bearing_fraction * myosin.length);
            }
            // Use the same productive-overlap capacity enforced by the KMC;
            // retained heads remain accessible throughout their wider shell.
            accessible_heads += std::max(retained_heads, available_slots);
            strong_heads += retained_strong_heads;
            }
        }
        accessible_heads = std::min(
            accessible_heads, std::max(1, am_max_heads_per_actin));
        strong_heads = std::min(strong_heads, accessible_heads);
        const double strong_fraction = accessible_heads > 0
            ? std::clamp(static_cast<double>(strong_heads) /
                         static_cast<double>(accessible_heads), 0.0, 1.0)
            : 0.0;
        if (am_thin_activation_enabled) {
            am_thin_activation[static_cast<size_t>(i)] =
                am_motor_kinetics::thin_activation_step(
                    am_thin_activation[static_cast<size_t>(i)],
                    am_activation, strong_fraction, dt,
                    am_q_calcium_open_rate, am_q_calcium_half,
                    am_q_calcium_hill,
                    am_thin_feedback_enabled
                        ? am_q_strong_feedback_rate : 0.0,
                    am_q_strong_half, am_q_strong_hill,
                    am_q_close_rate);
        } else {
            am_thin_activation[static_cast<size_t>(i)] = am_activation;
        }
        am_calcium_feature[i] = am_activation;
        am_q_feature[i] = am_thin_activation[static_cast<size_t>(i)];
    }
}

double Sarcomere::_am_axial_coordinate(int actin_index, int myosin_index) const {
    vec displacement = actin.center[actin_index] - myosin.center[myosin_index];
    displacement.pbc_wrap(box, is_periodic);
    return displacement.dot(actin.direction[actin_index]);
}

void Sarcomere::_initialize_am_head_pools() {
    am_heads_by_myosin.assign(
        static_cast<size_t>(myosin.n),
        std::array<std::vector<AMHeadState>, 2>{});
    for (int m = 0; m < myosin.n; ++m) {
        for (int half_index = 0; half_index < 2; ++half_index) {
            auto& heads = am_heads_by_myosin[static_cast<size_t>(m)]
                                               [static_cast<size_t>(half_index)];
            heads.assign(static_cast<size_t>(am_heads_per_half), AMHeadState{});
            const int half_sign = half_index == 0 ? -1 : 1;
            const int crown_count = std::max(1, (am_heads_per_half + 2) / 3);
            for (int h = 0; h < am_heads_per_half; ++h) {
                auto& head = heads[static_cast<size_t>(h)];
                const bool initially_on = am_motor_kinetics::initially_on_head(
                    h, am_heads_per_half, am_thick_initial_on_fraction);
                head.state = static_cast<BondFlag>(
                    !am_thick_feedback_enabled || initially_on
                        ? am_motor_kinetics::State::detached
                        : am_motor_kinetics::State::off_srx);
                const int crown = h / 3;
                const int spoke = h % 3;
                const double fraction =
                    (static_cast<double>(crown) + 0.5) /
                    static_cast<double>(crown_count);
                head.half = static_cast<int8_t>(half_sign);
                const double bare_half =
                    0.5 * geometry::am_bare_zone_fraction * myosin.length;
                head.anchor_axial = static_cast<double>(half_sign) *
                    (bare_half + fraction *
                        geometry::am_head_bearing_fraction * myosin.length);
                head.anchor_azimuth =
                    2.0 * M_PI * static_cast<double>(spoke) / 3.0 +
                    static_cast<double>(crown) * (40.0 * M_PI / 180.0);
                head.binding_azimuth = head.anchor_azimuth;
                head.stroke_progress = 0.0;
            }
        }
    }
}

int Sarcomere::_am_site_index(double actin_coordinate) const {
    const int count = std::max(1, am_sites_per_actin);
    const double length = std::max(actin.length, 1e-12);
    const double fraction = std::clamp(
        actin_coordinate / length + 0.5, 0.0,
        std::nextafter(1.0, 0.0));
    return std::clamp(static_cast<int>(std::floor(fraction * count)),
                      0, count - 1);
}

double Sarcomere::_am_site_coordinate(int site_index) const {
    const int count = std::max(1, am_sites_per_actin);
    const int site = std::clamp(site_index, 0, count - 1);
    return -0.5 * actin.length +
        (static_cast<double>(site) + 0.5) * actin.length /
            static_cast<double>(count);
}

double Sarcomere::_nearest_free_am_site_coordinate(
    int actin_index, int myosin_index, const AMHeadState& head,
    double binding_azimuth) const {
    if (actin_index < 0 || actin_index >= actin.n ||
        myosin_index < 0 || myosin_index >= myosin.n) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const auto continuous = compute_am_head_tether(
        actin, myosin, actin_index, myosin_index, box,
        head.anchor_axial, 0.0,
        myosin.radius + actin_radius + am_surface_gap,
        binding_azimuth, am_spatial_heads ? myosin.radius : 0.0,
        am_spatial_heads
            ? &myosin_material_e1[static_cast<size_t>(myosin_index)] : nullptr,
        am_spatial_heads
            ? &myosin_material_e2[static_cast<size_t>(myosin_index)] : nullptr,
        actin_radius);
    if (!continuous.projection_inside_actin) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const int count = std::max(1, am_sites_per_actin);
    double best_coordinate = std::numeric_limits<double>::quiet_NaN();
    double best_axial_error = std::numeric_limits<double>::infinity();
    for (int site = 0; site < count; ++site) {
        bool occupied = false;
        for (const auto& myosin_heads : am_heads_by_myosin) {
            for (const auto& half_heads : myosin_heads) {
                for (const auto& other : half_heads) {
                    if (&other == &head || other.bound_actin != actin_index) {
                        continue;
                    }
                    const auto state =
                        static_cast<am_motor_kinetics::State>(other.state);
                    if (state != am_motor_kinetics::State::weak &&
                        state != am_motor_kinetics::State::strong) {
                        continue;
                    }
                    if (_am_site_index(other.target_axial) == site) {
                        occupied = true;
                        break;
                    }
                }
                if (occupied) break;
            }
            if (occupied) break;
        }
        if (occupied) continue;

        const double coordinate = _am_site_coordinate(site);
        const auto geometry = compute_am_two_spring_head(
            actin, myosin, actin_index, myosin_index, box,
            coordinate, head.anchor_axial, binding_azimuth,
            am_spatial_heads ? myosin.radius : 0.0,
            am_spatial_heads
                ? &myosin_material_e1[static_cast<size_t>(myosin_index)] : nullptr,
            am_spatial_heads
                ? &myosin_material_e2[static_cast<size_t>(myosin_index)] : nullptr,
            actin_radius + am_surface_gap, 0.0, 0.0, 0.0,
            0.0, am_stall_force, false);
        if (geometry.head_reach_surface_distance > am_cutoff) continue;
        const double axial_error =
            std::abs(coordinate - continuous.actin_coordinate);
        if (axial_error < best_axial_error) {
            best_axial_error = axial_error;
            best_coordinate = coordinate;
        }
    }
    return best_coordinate;
}

void Sarcomere::_rebuild_am_site_claims() {
    const int site_count = std::max(1, am_sites_per_actin);
    am_site_claims_.assign(
        static_cast<size_t>(actin.n) * static_cast<size_t>(site_count), -1);
    am_site_claim_count_.assign(static_cast<size_t>(actin.n), 0);
    for (int j = 0; j < myosin.n; ++j) {
        for (int half = 0; half < 2; ++half) {
            const auto& heads = am_heads_by_myosin[static_cast<size_t>(j)]
                                                   [static_cast<size_t>(half)];
            for (int h = 0; h < static_cast<int>(heads.size()); ++h) {
                const auto& head = heads[static_cast<size_t>(h)];
                const auto state =
                    static_cast<am_motor_kinetics::State>(head.state);
                if (head.bound_actin < 0 ||
                    (state != am_motor_kinetics::State::weak &&
                     state != am_motor_kinetics::State::strong)) {
                    continue;
                }
                const int site = _am_site_index(head.target_axial);
                const size_t index =
                    static_cast<size_t>(head.bound_actin) *
                        static_cast<size_t>(site_count) +
                    static_cast<size_t>(site);
                if (am_site_claims_[index] < 0) {
                    am_site_claims_[index] =
                        ((j * 2 + half) * am_heads_per_half) + h;
                    ++am_site_claim_count_[static_cast<size_t>(
                        head.bound_actin)];
                }
            }
        }
    }
}

double Sarcomere::_claim_free_am_site_coordinate(
    int actin_index, int myosin_index, const AMHeadState& head,
    double binding_azimuth) {
    double claimed_coordinate = std::numeric_limits<double>::quiet_NaN();
    #pragma omp critical(am_site_claims)
    {
        const int site_count = std::max(1, am_sites_per_actin);
        if (actin_index >= 0 && actin_index < actin.n &&
            myosin_index >= 0 && myosin_index < myosin.n &&
            am_site_claims_.size() ==
                static_cast<size_t>(actin.n) *
                    static_cast<size_t>(site_count) &&
            am_site_claim_count_.size() == static_cast<size_t>(actin.n) &&
            am_site_claim_count_[static_cast<size_t>(actin_index)] <
                std::max(1, am_max_heads_per_actin)) {
            const auto continuous = compute_am_head_tether(
                actin, myosin, actin_index, myosin_index, box,
                head.anchor_axial, 0.0,
                myosin.radius + actin_radius + am_surface_gap,
                binding_azimuth,
                am_spatial_heads ? myosin.radius : 0.0,
                am_spatial_heads
                    ? &myosin_material_e1[static_cast<size_t>(myosin_index)]
                    : nullptr,
                am_spatial_heads
                    ? &myosin_material_e2[static_cast<size_t>(myosin_index)]
                    : nullptr,
                actin_radius);
            if (continuous.projection_inside_actin) {
                int best_site = -1;
                double best_axial_error =
                    std::numeric_limits<double>::infinity();
                const size_t offset = static_cast<size_t>(actin_index) *
                                      static_cast<size_t>(site_count);
                for (int site = 0; site < site_count; ++site) {
                    if (am_site_claims_[offset + static_cast<size_t>(site)] >= 0) {
                        continue;
                    }
                    const double coordinate = _am_site_coordinate(site);
                    const auto geometry = compute_am_two_spring_head(
                        actin, myosin, actin_index, myosin_index, box,
                        coordinate, head.anchor_axial, binding_azimuth,
                        am_spatial_heads ? myosin.radius : 0.0,
                        am_spatial_heads
                            ? &myosin_material_e1[static_cast<size_t>(
                                  myosin_index)]
                            : nullptr,
                        am_spatial_heads
                            ? &myosin_material_e2[static_cast<size_t>(
                                  myosin_index)]
                            : nullptr,
                        actin_radius + am_surface_gap, 0.0, 0.0, 0.0,
                        0.0, am_stall_force, false);
                    if (geometry.head_reach_surface_distance > am_cutoff) {
                        continue;
                    }
                    const double axial_error = std::abs(
                        coordinate - continuous.actin_coordinate);
                    if (axial_error < best_axial_error) {
                        best_axial_error = axial_error;
                        best_site = site;
                    }
                }
                if (best_site >= 0) {
                    am_site_claims_[offset + static_cast<size_t>(best_site)] = 1;
                    ++am_site_claim_count_[static_cast<size_t>(actin_index)];
                    claimed_coordinate = _am_site_coordinate(best_site);
                }
            }
        }
    }
    return claimed_coordinate;
}

void Sarcomere::_release_am_site_claim(
    int actin_index, double actin_coordinate) {
    #pragma omp critical(am_site_claims)
    {
        const int site_count = std::max(1, am_sites_per_actin);
        if (actin_index >= 0 && actin_index < actin.n &&
            am_site_claims_.size() ==
                static_cast<size_t>(actin.n) *
                    static_cast<size_t>(site_count) &&
            am_site_claim_count_.size() == static_cast<size_t>(actin.n)) {
            const size_t index = static_cast<size_t>(actin_index) *
                    static_cast<size_t>(site_count) +
                static_cast<size_t>(_am_site_index(actin_coordinate));
            if (am_site_claims_[index] >= 0) {
                am_site_claims_[index] = -1;
                am_site_claim_count_[static_cast<size_t>(actin_index)] =
                    std::max(0, am_site_claim_count_[static_cast<size_t>(
                                    actin_index)] - 1);
            }
        }
    }
}

bool Sarcomere::_am_actin_reaches_half(
    int actin_index, int myosin_index, int half) const {
    if (actin_index < 0 || actin_index >= actin.n ||
        myosin_index < 0 || myosin_index >= myosin.n) return false;
    vec axis_m = myosin.direction[myosin_index];
    if (axis_m.norm() <= 1e-12) return false;
    axis_m.normalize();
    vec left_disp = actin.left_end[actin_index] - myosin.center[myosin_index];
    vec right_disp = actin.right_end[actin_index] - myosin.center[myosin_index];
    left_disp.pbc_wrap(box, is_periodic);
    right_disp.pbc_wrap(box, is_periodic);
    const double a_lo = std::min(left_disp.dot(axis_m), right_disp.dot(axis_m));
    const double a_hi = std::max(left_disp.dot(axis_m), right_disp.dot(axis_m));
    const double bare_half =
        0.5 * geometry::am_bare_zone_fraction * myosin.length;
    const double h_lo = half < 0 ? -0.5 * myosin.length : bare_half;
    const double h_hi = half < 0 ? -bare_half : 0.5 * myosin.length;
    return std::min(a_hi, h_hi) - std::max(a_lo, h_lo) > EPS;
}

bool Sarcomere::_am_half_is_productive(
    int actin_index, int myosin_index, int half) const {
    if (!directional) return true;
    if (actin_index < 0 || actin_index >= actin.n ||
        myosin_index < 0 || myosin_index >= myosin.n) return false;
    const int productive_half =
        actn2_load_kinetics::productive_myosin_half(
            actin.direction[actin_index].dot(
                myosin.direction[myosin_index]));
    if (half == productive_half) return true;
    // Opposite half: productive only once the actin has actually reached past
    // the bare zone into this half. Pull direction comes from actin polarity,
    // so these heads continue the motion rather than opposing it.
    return am_crossover_pull_through &&
           _am_actin_reaches_half(actin_index, myosin_index, half);
}

bool Sarcomere::_am_head_is_productive(
    int actin_index, int myosin_index, const AMHeadState& head) const {
    return _am_half_is_productive(actin_index, myosin_index,
                                  static_cast<int>(head.half));
}

void Sarcomere::_enforce_global_am_site_ownership() {
    const int site_count = std::max(1, am_sites_per_actin);
    const size_t total_sites = static_cast<size_t>(actin.n) *
                               static_cast<size_t>(site_count);
    std::vector<int> winner(total_sites, -1);
    std::vector<int> winner_state(total_sites, -1);
    std::vector<double> winner_energy(
        total_sites, std::numeric_limits<double>::infinity());
    std::vector<unsigned char> retained_site(total_sites, 0);

    auto owner_code = [&](int myosin_index, int half_index, int head_index) {
        return ((myosin_index * 2 + half_index) * am_heads_per_half) +
               head_index;
    };

    // Build one deterministic owner for every material actin site. Strong
    // chemistry wins over weak chemistry; ties keep the lower-strain head,
    // then the stable myosin/head identifier. Initializers and simultaneous
    // KMC proposals therefore cannot duplicate a site across myosins.
    for (int j = 0; j < myosin.n; ++j) {
        for (int half = 0; half < 2; ++half) {
            auto& heads = am_heads_by_myosin[static_cast<size_t>(j)]
                                              [static_cast<size_t>(half)];
            for (int h = 0; h < static_cast<int>(heads.size()); ++h) {
                auto& head = heads[static_cast<size_t>(h)];
                const auto state =
                    static_cast<am_motor_kinetics::State>(head.state);
                if (head.bound_actin < 0 ||
                    (state != am_motor_kinetics::State::weak &&
                     state != am_motor_kinetics::State::strong)) continue;
                const int site = _am_site_index(head.target_axial);
                head.target_axial = _am_site_coordinate(site);
                const size_t global_site =
                    static_cast<size_t>(head.bound_actin) *
                        static_cast<size_t>(site_count) +
                    static_cast<size_t>(site);
                const int state_priority =
                    state == am_motor_kinetics::State::strong ? 1 : 0;
                const int code = owner_code(j, half, h);
                const double energy = _am_head_strain_energy(
                    head.bound_actin, j, _am_head_id(half, h), state);
                const bool better = winner[global_site] < 0 ||
                    state_priority > winner_state[global_site] ||
                    (state_priority == winner_state[global_site] &&
                     (energy < winner_energy[global_site] - 1e-15 ||
                      (std::abs(energy - winner_energy[global_site]) <= 1e-15 &&
                       code < winner[global_site])));
                if (better) {
                    winner[global_site] = code;
                    winner_state[global_site] = state_priority;
                    winner_energy[global_site] = energy;
                }
            }
        }
    }

    // Site exclusion alone prevents two heads using the same monomer-scale
    // location. A separate actin-level cap prevents every neighboring myosin
    // from recreating another 24-head allowance on different sites. The same
    // chemistry/strain ordering makes this global arbitration deterministic.
    const int global_cap = std::max(1, am_max_heads_per_actin);
    for (int a = 0; a < actin.n; ++a) {
        std::vector<size_t> occupied_sites;
        occupied_sites.reserve(static_cast<size_t>(site_count));
        const size_t offset = static_cast<size_t>(a) *
                              static_cast<size_t>(site_count);
        for (int site = 0; site < site_count; ++site) {
            const size_t index = offset + static_cast<size_t>(site);
            if (winner[index] >= 0) occupied_sites.push_back(index);
        }
        std::sort(
            occupied_sites.begin(), occupied_sites.end(),
            [&](size_t left, size_t right) {
                if (winner_state[left] != winner_state[right]) {
                    return winner_state[left] > winner_state[right];
                }
                if (winner_energy[left] != winner_energy[right]) {
                    return winner_energy[left] < winner_energy[right];
                }
                return winner[left] < winner[right];
            });
        const int retained = std::min(
            global_cap, static_cast<int>(occupied_sites.size()));
        for (int index = 0; index < retained; ++index) {
            retained_site[occupied_sites[static_cast<size_t>(index)]] = 1;
        }
    }

    for (int j = 0; j < myosin.n; ++j) {
        for (int half = 0; half < 2; ++half) {
            auto& heads = am_heads_by_myosin[static_cast<size_t>(j)]
                                              [static_cast<size_t>(half)];
            for (int h = 0; h < static_cast<int>(heads.size()); ++h) {
                auto& head = heads[static_cast<size_t>(h)];
                const auto state =
                    static_cast<am_motor_kinetics::State>(head.state);
                if (head.bound_actin < 0 ||
                    (state != am_motor_kinetics::State::weak &&
                     state != am_motor_kinetics::State::strong)) continue;
                const int site = _am_site_index(head.target_axial);
                const size_t global_site =
                    static_cast<size_t>(head.bound_actin) *
                        static_cast<size_t>(site_count) +
                    static_cast<size_t>(site);
                if (winner[global_site] == owner_code(j, half, h) &&
                    retained_site[global_site] != 0) continue;
                const int old_actin = head.bound_actin;
                _record_am_transition(
                    old_actin, j, _am_head_id(half, h), state,
                    am_motor_kinetics::State::detached,
                    AMTransitionReason::site_occupied);
                head.state = static_cast<BondFlag>(
                    am_motor_kinetics::State::detached);
                head.bound_actin = -1;
                head.stroke_progress = 0.0;
            }
        }
    }
}

double Sarcomere::_am_head_association_gate(
    int actin_index, int myosin_index, const AMHeadState& head,
    double binding_azimuth) const {
    if (actin_index < 0 || actin_index >= actin.n ||
        myosin_index < 0 || myosin_index >= myosin.n) return 0.0;
    const auto geometry = compute_am_head_tether(
        actin, myosin, actin_index, myosin_index, box,
        head.anchor_axial, 0.0,
        myosin.radius + actin_radius + am_surface_gap,
        binding_azimuth, am_spatial_heads ? myosin.radius : 0.0,
        am_spatial_heads
            ? &myosin_material_e1[static_cast<size_t>(myosin_index)] : nullptr,
        am_spatial_heads
            ? &myosin_material_e2[static_cast<size_t>(myosin_index)] : nullptr,
        actin_radius);
    // A head is an association candidate only when its material anchor
    // projects onto the finite actin segment.  Clamping an outside projection
    // to an endpoint is useful for measuring an already-bound connector, but
    // must never turn a head beyond the filament end into a new binding event.
    if (!geometry.projection_inside_actin) return 0.0;
    const double site_coordinate = _am_site_coordinate(
        _am_site_index(geometry.actin_coordinate));
    const auto site_geometry = compute_am_two_spring_head(
        actin, myosin, actin_index, myosin_index, box,
        site_coordinate, head.anchor_axial, binding_azimuth,
        am_spatial_heads ? myosin.radius : 0.0,
        am_spatial_heads
            ? &myosin_material_e1[static_cast<size_t>(myosin_index)] : nullptr,
        am_spatial_heads
            ? &myosin_material_e2[static_cast<size_t>(myosin_index)] : nullptr,
        actin_radius + am_surface_gap, 0.0, 0.0, 0.0,
        0.0, am_stall_force, false);
    const double capture_reach = am_cutoff;
    if (!(site_geometry.head_reach_surface_distance <= capture_reach)) {
        return 0.0;
    }
    const double spatial_scale = std::max(am_optimal, 1e-12);
    return std::exp(-0.5 * std::pow(
        site_geometry.head_reach_surface_distance /
            spatial_scale,
        2.0));
}

bool Sarcomere::_am_head_retention_reachable(
    int actin_index, int myosin_index, const AMHeadState& head) const {
    if (actin_index < 0 || actin_index >= actin.n ||
        myosin_index < 0 || myosin_index >= myosin.n) return false;
    vec axis_m = myosin.direction[myosin_index];
    axis_m.normalize();
    vec left_disp = actin.left_end[actin_index] - myosin.center[myosin_index];
    vec right_disp = actin.right_end[actin_index] - myosin.center[myosin_index];
    left_disp.pbc_wrap(box, is_periodic);
    right_disp.pbc_wrap(box, is_periodic);
    const double actin_axial_min = std::min(
        left_disp.dot(axis_m), right_disp.dot(axis_m));
    const double actin_axial_max = std::max(
        left_disp.dot(axis_m), right_disp.dot(axis_m));
    const double bare_half =
        0.5 * geometry::am_bare_zone_fraction * myosin.length;
    const double head_axial_min = head.half < 0
        ? -0.5 * myosin.length : bare_half;
    const double head_axial_max = head.half < 0
        ? -bare_half : 0.5 * myosin.length;
    if (directional && !(std::min(actin_axial_max, head_axial_max) -
                         std::max(actin_axial_min, head_axial_min) > EPS)) {
        return false;
    }

    const auto geometry = compute_am_two_spring_head(
        actin, myosin, actin_index, myosin_index, box,
        head.target_axial, head.anchor_axial, head.binding_azimuth,
        am_spatial_heads ? myosin.radius : 0.0,
        am_spatial_heads
            ? &myosin_material_e1[static_cast<size_t>(myosin_index)] : nullptr,
        am_spatial_heads
            ? &myosin_material_e2[static_cast<size_t>(myosin_index)] : nullptr,
        actin_radius + am_surface_gap, 0.0, 0.0, 0.0,
        head.stroke_progress, am_stall_force, false);
    // Measure extension relative to the relaxed two-spring connector, not
    // absolute distance from the thick-filament anchor to the actin surface.
    // The former is zero at the preferred AM spacing.
    const double excess_extension = std::max(
        0.0, geometry.head_reach_surface_distance);
    return excess_extension <= am_rupture_distance;
}

int Sarcomere::_am_bound_head_count(
    int actin_index, int myosin_index) const {
    if (myosin_index < 0 || myosin_index >= myosin.n) return 0;
    int count = 0;
    for (const auto& half_heads :
         am_heads_by_myosin[static_cast<size_t>(myosin_index)]) {
        for (const auto& head : half_heads) {
            if (head.bound_actin == actin_index &&
                head.state != static_cast<BondFlag>(
                    am_motor_kinetics::State::detached) &&
                head.state != static_cast<BondFlag>(
                    am_motor_kinetics::State::off_srx)) {
                ++count;
            }
        }
    }
    return count;
}

double Sarcomere::_am_per_head_tether_stiffness(
    int bound_head_count) const {
    return am_motor_kinetics::occupancy_scaled_per_head_stiffness(
        am_tether_k, bound_head_count, am_tether_saturation_heads);
}

double Sarcomere::_am_head_strain_energy(
    int actin_index, int myosin_index, int head_index,
    am_motor_kinetics::State state) const {
    const auto* head = _am_head(myosin_index, head_index);
    if (head == nullptr || head->bound_actin != actin_index ||
        (state != am_motor_kinetics::State::weak &&
         state != am_motor_kinetics::State::strong)) return 0.0;
    // Pair-level centerline containment remains separate. This optional term
    // is the physical transverse compliance of the individual, material-point
    // head connector and therefore contributes to its stochastic off-rate.
    const double radial_k = am_head_transverse_k;
    const bool strong = state == am_motor_kinetics::State::strong;
    const double axial_k = am_motor_k *
        (strong ? 1.0 : am_weak_stiffness_ratio);
    const auto mechanics = compute_am_two_spring_head(
        actin, myosin, actin_index, myosin_index, box,
        head->target_axial, head->anchor_axial, head->binding_azimuth,
        am_spatial_heads ? myosin.radius : 0.0,
        am_spatial_heads
            ? &myosin_material_e1[static_cast<size_t>(myosin_index)] : nullptr,
        am_spatial_heads
            ? &myosin_material_e2[static_cast<size_t>(myosin_index)] : nullptr,
        actin_radius + am_surface_gap, radial_k, am_tether_slack,
        axial_k, head->stroke_progress, am_stall_force, strong);
    if (!strong) return mechanics.energy;

    // Productive strong-head strain is already catch-stabilized through force.
    // Only transverse distortion and negative (dragged-ahead) axial strain
    // accelerate strong-state loss.
    const double negative_strain = std::max(
        0.0, mechanics.longitudinal_extension - head->stroke_progress);
    const double negative_energy = 0.5 * am_motor_k *
        am_weak_stiffness_ratio * negative_strain * negative_strain;
    return mechanics.radial_energy + negative_energy;
}

std::array<double, 2> Sarcomere::_am_head_detachment_energies(
    int actin_index, int myosin_index, int head_index,
    am_motor_kinetics::State state) const {
    const auto* head = _am_head(myosin_index, head_index);
    if (head == nullptr || head->bound_actin != actin_index ||
        (state != am_motor_kinetics::State::weak &&
         state != am_motor_kinetics::State::strong)) {
        return {0.0, 0.0};
    }
    const bool strong = state == am_motor_kinetics::State::strong;
    const double axial_k = am_motor_k *
        (strong ? 1.0 : am_weak_stiffness_ratio);
    const auto mechanics = compute_am_two_spring_head(
        actin, myosin, actin_index, myosin_index, box,
        head->target_axial, head->anchor_axial, head->binding_azimuth,
        am_spatial_heads ? myosin.radius : 0.0,
        am_spatial_heads
            ? &myosin_material_e1[static_cast<size_t>(myosin_index)] : nullptr,
        am_spatial_heads
            ? &myosin_material_e2[static_cast<size_t>(myosin_index)] : nullptr,
        actin_radius + am_surface_gap, am_head_transverse_k,
        am_tether_slack, axial_k, head->stroke_progress,
        am_stall_force, strong);
    // Weak binding is non-stereospecific in this coarse-grained model. Its
    // chemical off-rate is not accelerated by axial sliding before locking.
    double negative_axial_energy = 0.0;
    if (strong) {
        const double negative_strain = std::max(
            0.0, mechanics.longitudinal_extension - head->stroke_progress);
        negative_axial_energy = 0.5 * am_motor_k *
            am_weak_stiffness_ratio * negative_strain * negative_strain;
    }
    return {negative_axial_energy, mechanics.radial_energy};
}

double Sarcomere::_am_head_detachment_rate(
    int actin_index, int myosin_index, int head_index,
    am_motor_kinetics::State state, double chemical_rate) const {
    if (!am_strain_dependent_detachment) return chemical_rate;
    const double beta_value = thermal_beta > 0.0
        ? thermal_beta : rod_mobility::beta(temperature_kelvin);
    const auto energies = _am_head_detachment_energies(
        actin_index, myosin_index, head_index, state);
    return am_motor_kinetics::anisotropic_strain_accelerated_rate(
        chemical_rate, beta_value,
        energies[0], am_strain_detachment_theta,
        energies[1], am_transverse_strain_detachment_theta,
        am_strain_detachment_max_exponent);
}

int Sarcomere::_am_head_id(int half_index, int head_index) const {
    return half_index * am_heads_per_half + head_index;
}

Sarcomere::AMHeadState* Sarcomere::_am_head(int myosin_index, int head_id) {
    if (myosin_index < 0 || myosin_index >= myosin.n || head_id < 0 ||
        am_heads_per_half <= 0) return nullptr;
    const int half_index = head_id / am_heads_per_half;
    const int head_index = head_id % am_heads_per_half;
    if (half_index < 0 || half_index >= 2) return nullptr;
    auto& heads = am_heads_by_myosin[static_cast<size_t>(myosin_index)]
                                       [static_cast<size_t>(half_index)];
    if (head_index >= static_cast<int>(heads.size())) return nullptr;
    return &heads[static_cast<size_t>(head_index)];
}

const Sarcomere::AMHeadState* Sarcomere::_am_head(
    int myosin_index, int head_id) const {
    if (myosin_index < 0 || myosin_index >= myosin.n || head_id < 0 ||
        am_heads_per_half <= 0) return nullptr;
    const int half_index = head_id / am_heads_per_half;
    const int head_index = head_id % am_heads_per_half;
    if (half_index < 0 || half_index >= 2) return nullptr;
    const auto& heads = am_heads_by_myosin[static_cast<size_t>(myosin_index)]
                                             [static_cast<size_t>(half_index)];
    if (head_index >= static_cast<int>(heads.size())) return nullptr;
    return &heads[static_cast<size_t>(head_index)];
}

double Sarcomere::_am_productive_strain(
    int actin_index, int myosin_index, int head_index) const {
    double target = am_target_axial[actin_index][myosin_index];
    if (head_index >= 0) {
        const auto* head = _am_head(myosin_index, head_index);
        if (head == nullptr || head->bound_actin != actin_index) return 0.0;
        const auto geometry = compute_am_two_spring_head(
            actin, myosin, actin_index, myosin_index, box,
            head->target_axial, head->anchor_axial, head->binding_azimuth,
            am_spatial_heads ? myosin.radius : 0.0,
            am_spatial_heads
                ? &myosin_material_e1[static_cast<size_t>(myosin_index)] : nullptr,
            am_spatial_heads
                ? &myosin_material_e2[static_cast<size_t>(myosin_index)] : nullptr,
            actin_radius + am_surface_gap, 0.0, 0.0, 0.0,
            head->stroke_progress, am_stall_force, true);
        return geometry.productive_strain;
    }
    return target -
           _am_axial_coordinate(actin_index, myosin_index);
}

double Sarcomere::_am_head_active_force(
    int actin_index, int myosin_index, int head_index) const {
    const auto* head = _am_head(myosin_index, head_index);
    if (head == nullptr || head->bound_actin != actin_index ||
        head->state != static_cast<BondFlag>(am_motor_kinetics::State::strong)) return 0.0;
    return am_motor_kinetics::active_crossbridge_force(
        am_motor_k, _am_productive_strain(actin_index, myosin_index, head_index),
        am_stall_force, am_bidirectional_crossbridge);
}

double Sarcomere::_am_active_force(int actin_index, int myosin_index) const {
    if (!dynamic_am) return 0.0;
    double force = 0.0;
    for (int half_index = 0; half_index < 2; ++half_index) {
        const auto& heads = am_heads_by_myosin[static_cast<size_t>(myosin_index)]
                                               [static_cast<size_t>(half_index)];
        for (int h = 0; h < static_cast<int>(heads.size()); ++h) {
            force += _am_head_active_force(
                actin_index, myosin_index, _am_head_id(half_index, h));
        }
    }
    return force;
}

double Sarcomere::_am_resisting_force(int actin_index, int myosin_index) const {
    if (!dynamic_am) {
        return 0.0;
    }
    const auto state = static_cast<am_motor_kinetics::State>(am_state[actin_index][myosin_index]);
    if (state == am_motor_kinetics::State::detached) return 0.0;
    if (state != am_motor_kinetics::State::strong) return 0.0;
    return _am_active_force(actin_index, myosin_index);
}

void Sarcomere::_record_am_transition(
    int actin_index, int myosin_index, int head_index,
    am_motor_kinetics::State old_state,
    am_motor_kinetics::State new_state,
    AMTransitionReason reason) {
    if (!dynamic_am || old_state == new_state) return;
    const auto diagnostic_state = new_state != am_motor_kinetics::State::detached
        ? new_state : old_state;
    const double strain = _am_productive_strain(
        actin_index, myosin_index, head_index);
    // The passive proximity tether is deliberately excluded. Weak heads carry
    // no active axial force, and strong heads are tension-only and stall-capped.
    const double active_force = diagnostic_state == am_motor_kinetics::State::strong
        ? am_motor_kinetics::active_crossbridge_force(
              am_motor_k, strain, am_stall_force)
        : 0.0;
    const int thread_id = omp_in_parallel() ? omp_get_thread_num() : 0;
    auto& events = am_transition_events_temp[thread_id];
    events.insert(events.end(), {
        static_cast<double>(actin_index),
        static_cast<double>(myosin_index),
        static_cast<double>(current_step),
        static_cast<double>(static_cast<int>(old_state)),
        static_cast<double>(static_cast<int>(new_state)),
        static_cast<double>(static_cast<int>(reason)),
        am_activation,
        strain,
        active_force,
        static_cast<double>(previous_actin_cb_status[actin_index]),
        static_cast<double>(head_index)});
}

Sarcomere::AMRetentionDiagnostic Sarcomere::_am_retention_diagnostic(
    int actin_index, int myosin_index, int head_index,
    bool candidate_missing) const {
    AMRetentionDiagnostic diagnostic;
    diagnostic.rupture_distance = am_rupture_distance;
    if (actin_index < 0 || actin_index >= actin.n ||
        myosin_index < 0 || myosin_index >= myosin.n) {
        diagnostic.failure_mask = am_unexpected_candidate_loss;
        return diagnostic;
    }

    const auto& myosin_neighbors =
        actin_neighbors_by_species[static_cast<size_t>(actin_index)].second;
    const bool neighbor_present =
        std::find(myosin_neighbors.begin(), myosin_neighbors.end(),
                  myosin_index) != myosin_neighbors.end();
    if (!neighbor_present) {
        diagnostic.failure_mask |= am_not_in_neighbor_list;
    }

    diagnostic.segment_distance = geometry::segment_segment_distance(
        actin.left_end[actin_index], actin.right_end[actin_index],
        myosin.left_end[myosin_index], myosin.right_end[myosin_index],
        box, is_periodic);
    if (diagnostic.segment_distance > am_rupture_distance) {
        diagnostic.failure_mask |= am_beyond_rupture_distance;
    }

    diagnostic.actin_myosin_dot = std::clamp(
        actin.direction[actin_index].dot(myosin.direction[myosin_index]),
        -1.0, 1.0);
    diagnostic.expected_half = directional
        ? actn2_load_kinetics::productive_myosin_half(
              diagnostic.actin_myosin_dot)
        : actn2_load_kinetics::myosin_half([&]() {
              vec displacement = actin.center[actin_index] -
                                 myosin.center[myosin_index];
              displacement.pbc_wrap(box, is_periodic);
              return displacement.dot(myosin.direction[myosin_index]);
          }());
    const auto* head = _am_head(myosin_index, head_index);
    diagnostic.head_half = head == nullptr ? 0 : static_cast<int>(head->half);
    // Under pull-through a far-half head is legitimate once the actin has
    // reached into it, so it must not be reported as a wrong-half failure.
    const bool crossover_ok = am_crossover_pull_through && head != nullptr &&
        _am_actin_reaches_half(actin_index, myosin_index,
                               static_cast<int>(head->half));
    if (diagnostic.head_half != diagnostic.expected_half && !crossover_ok) {
        diagnostic.failure_mask |= am_wrong_productive_half;
    }
    diagnostic.capture_reach = am_surface_gap + am_cutoff;
    diagnostic.retention_reach = am_rupture_distance;
    if (head != nullptr) {
        const auto head_geometry = compute_am_two_spring_head(
            actin, myosin, actin_index, myosin_index, box,
            head->target_axial, head->anchor_axial, head->binding_azimuth,
            am_spatial_heads ? myosin.radius : 0.0,
            am_spatial_heads
                ? &myosin_material_e1[static_cast<size_t>(myosin_index)] : nullptr,
            am_spatial_heads
                ? &myosin_material_e2[static_cast<size_t>(myosin_index)] : nullptr,
            actin_radius + am_surface_gap, 0.0, 0.0, 0.0,
            head->stroke_progress, am_stall_force, false);
        diagnostic.head_connector_distance =
            head_geometry.connector.norm();
        diagnostic.head_surface_distance =
            diagnostic.head_connector_distance - actin_radius;
        diagnostic.head_excess_extension = std::max(
            0.0, head_geometry.head_reach_surface_distance);
        diagnostic.strain_energy = _am_head_strain_energy(
            actin_index, myosin_index, head_index,
            static_cast<am_motor_kinetics::State>(head->state));
        if (diagnostic.head_excess_extension > diagnostic.retention_reach) {
            diagnostic.failure_mask |= am_head_beyond_retention_reach;
        }
    }

    vec axis_m = myosin.direction[myosin_index];
    axis_m.normalize();
    vec left_disp = actin.left_end[actin_index] - myosin.center[myosin_index];
    vec right_disp = actin.right_end[actin_index] - myosin.center[myosin_index];
    left_disp.pbc_wrap(box, is_periodic);
    right_disp.pbc_wrap(box, is_periodic);
    const double actin_axial_min = std::min(
        left_disp.dot(axis_m), right_disp.dot(axis_m));
    const double actin_axial_max = std::max(
        left_disp.dot(axis_m), right_disp.dot(axis_m));
    const double bare_half =
        0.5 * geometry::am_bare_zone_fraction * myosin.length;
    const double head_axial_min = diagnostic.expected_half < 0
        ? -0.5 * myosin.length : bare_half;
    const double head_axial_max = diagnostic.expected_half < 0
        ? -bare_half : 0.5 * myosin.length;
    diagnostic.productive_axial_overlap = std::max(
        0.0, std::min(actin_axial_max, head_axial_max) -
             std::max(actin_axial_min, head_axial_min));
    if (directional && !(diagnostic.productive_axial_overlap > EPS)) {
        diagnostic.failure_mask |= am_no_productive_axial_overlap;
    }

    diagnostic.binding_ratio =
        am_interaction[actin_index][myosin_index].myosin_binding_ratio;
    diagnostic.partial_ratio =
        am_interaction[actin_index][myosin_index].partial_binding_ratio;
    const bool association_geometry = diagnostic.binding_ratio > EPS &&
        (diagnostic.partial_ratio > EPS || !directional);
    if (!association_geometry) {
        diagnostic.failure_mask |= am_association_gate_closed;
    }

    const bool retention_geometry =
        diagnostic.head_half == diagnostic.expected_half &&
        (!directional || diagnostic.productive_axial_overlap > EPS) &&
        diagnostic.head_excess_extension <= diagnostic.retention_reach;
    if (candidate_missing && retention_geometry) {
        diagnostic.failure_mask |= am_unexpected_candidate_loss;
    }
    return diagnostic;
}

void Sarcomere::_record_am_forced_rupture(
    int actin_index, int myosin_index, int head_index,
    am_motor_kinetics::State old_state, AMTransitionReason reason,
    bool candidate_missing) {
    const auto diagnostic = _am_retention_diagnostic(
        actin_index, myosin_index, head_index, candidate_missing);
    // Both buffers are shared members appended from inside the `omp for` in
    // _process_dynamic_am_heads; concurrent vector reallocation corrupts the
    // heap.  Serialise the appends (rare event, no measurable cost).
    #pragma omp critical(am_rupture_event_log)
    {
    am_forced_rupture_events.insert(am_forced_rupture_events.end(), {
        static_cast<double>(actin_index),
        static_cast<double>(myosin_index),
        static_cast<double>(current_step),
        static_cast<double>(head_index),
        static_cast<double>(static_cast<int>(old_state)),
        static_cast<double>(static_cast<int>(reason)),
        static_cast<double>(diagnostic.failure_mask),
        diagnostic.segment_distance,
        diagnostic.rupture_distance,
        diagnostic.actin_myosin_dot,
        static_cast<double>(diagnostic.expected_half),
        static_cast<double>(diagnostic.head_half),
        diagnostic.productive_axial_overlap,
        diagnostic.binding_ratio,
        diagnostic.partial_ratio});
    am_forced_rupture_events_v2.insert(am_forced_rupture_events_v2.end(), {
        static_cast<double>(actin_index), static_cast<double>(myosin_index),
        static_cast<double>(current_step), static_cast<double>(head_index),
        static_cast<double>(static_cast<int>(old_state)),
        static_cast<double>(static_cast<int>(reason)),
        static_cast<double>(diagnostic.failure_mask), diagnostic.segment_distance,
        diagnostic.rupture_distance, diagnostic.actin_myosin_dot,
        static_cast<double>(diagnostic.expected_half),
        static_cast<double>(diagnostic.head_half),
        diagnostic.productive_axial_overlap, diagnostic.binding_ratio,
        diagnostic.partial_ratio, diagnostic.head_connector_distance,
        diagnostic.head_surface_distance, diagnostic.capture_reach,
        diagnostic.retention_reach});
    }
}

void Sarcomere::_record_am_thick_transition(
    int myosin_index, int half_index, int head_id,
    am_motor_kinetics::State old_state,
    am_motor_kinetics::State new_state,
    double half_tension) {
    if (!dynamic_am || old_state == new_state) return;
    // Shared member appended from inside the parallel region -- see
    // am_rupture_event_log above.
    #pragma omp critical(am_thick_transition_log)
    am_thick_transition_events.insert(am_thick_transition_events.end(), {
        static_cast<double>(myosin_index),
        static_cast<double>(half_index == 0 ? -1 : 1),
        static_cast<double>(head_id),
        static_cast<double>(current_step),
        static_cast<double>(static_cast<int>(old_state)),
        static_cast<double>(static_cast<int>(new_state)),
        half_tension});
}

void Sarcomere::_sync_am_pair_summary(int actin_index, int myosin_index) {
    int summary_state = 0;
    double target_sum = 0.0;
    int target_count = 0;
    for (const auto& half_heads :
         am_heads_by_myosin[static_cast<size_t>(myosin_index)]) {
        for (const auto& head : half_heads) {
            if (head.bound_actin != actin_index) continue;
            summary_state = std::max(summary_state, static_cast<int>(head.state));
            if (head.state != static_cast<BondFlag>(am_motor_kinetics::State::detached)) {
                target_sum += head.target_axial;
                ++target_count;
            }
        }
    }
    am_state[actin_index][myosin_index] = static_cast<BondFlag>(summary_state);
    am_target_axial[actin_index][myosin_index] = target_count > 0
        ? target_sum / static_cast<double>(target_count) : 0.0;
}

void Sarcomere::_detach_am_heads_for_pair(
    int actin_index, int myosin_index, AMTransitionReason reason) {
    for (int half_index = 0; half_index < 2; ++half_index) {
        auto& heads = am_heads_by_myosin[static_cast<size_t>(myosin_index)]
                                         [static_cast<size_t>(half_index)];
        for (int h = 0; h < static_cast<int>(heads.size()); ++h) {
            auto& head = heads[static_cast<size_t>(h)];
            if (head.bound_actin != actin_index) continue;
            const auto state = static_cast<am_motor_kinetics::State>(head.state);
            if (state != am_motor_kinetics::State::detached) {
                _record_am_forced_rupture(
                    actin_index, myosin_index, _am_head_id(half_index, h),
                    state, reason, reason == AMTransitionReason::neighbor_loss);
                _record_am_transition(
                    actin_index, myosin_index, _am_head_id(half_index, h),
                    state, am_motor_kinetics::State::detached, reason);
            }
            head.state = static_cast<BondFlag>(am_motor_kinetics::State::detached);
            head.bound_actin = -1;
        }
    }
    _sync_am_pair_summary(actin_index, myosin_index);
}

void Sarcomere::_drop_unseen_am_bonds() {
    for (const auto& [i, j] : previous_am_pairs) {
        if (am_state[i][j] != static_cast<BondFlag>(am_motor_kinetics::State::detached) &&
            am_seen[i][j] == 0) {
            _detach_am_heads_for_pair(i, j, AMTransitionReason::neighbor_loss);
            am_bonds[i][j] = 0;
            myosinIndicesPerActin.deleteConnection(i, j);
            for (auto& temp_conn : actinIndicesPerMyosin_temp) {
                temp_conn.deleteConnection(j, i);
            }
        }
    }
}



void Sarcomere::_enforce_myosin_bond_limit() {
    for (int j = 0; j < myosin.n; ++j) {
        struct BoundInfo {
            int actin_idx;
            bool cb;
            bool prev;
            int motor_state;
            double dot;
            bool active;
        };

        std::vector<BoundInfo> bound;
        bound.reserve(actin.n);
        for (int i = 0; i < actin.n; ++i) {
            if (am_bonds[i][j] == 1) {
                double dot = actin.direction[i].dot(myosin.direction[j]);
                bool cb = actin.cb_status[i] >= 2;
                bool prev = am_bonds_prev[i][j] == 1;
                int motor_state = dynamic_am ? static_cast<int>(am_state[i][j]) : 2;
                bound.push_back({i, cb, prev, motor_state, dot, true});
            }
        }

        if (bound.empty()) {
            continue;
        }

        auto remove_connection = [&](BoundInfo& info) {
            if (!info.active) {
                return;
            }
            am_bonds[info.actin_idx][j] = 0;
            if (dynamic_am) {
                int weak_heads = 0;
                int strong_heads = 0;
                for (const auto& half_heads :
                     am_heads_by_myosin[static_cast<size_t>(j)]) {
                    for (const auto& head : half_heads) {
                        if (head.bound_actin != info.actin_idx) continue;
                        if (head.state == static_cast<BondFlag>(am_motor_kinetics::State::weak)) {
                            ++weak_heads;
                        } else if (head.state == static_cast<BondFlag>(am_motor_kinetics::State::strong)) {
                            ++strong_heads;
                        }
                    }
                }
                actin["am_weak_count"][info.actin_idx] = std::max(
                    0.0, actin["am_weak_count"][info.actin_idx] - weak_heads);
                actin["am_strong_count"][info.actin_idx] = std::max(
                    0.0, actin["am_strong_count"][info.actin_idx] - strong_heads);
                _detach_am_heads_for_pair(
                    info.actin_idx, j, AMTransitionReason::capacity_limit);
            }
            myosinIndicesPerActin.deleteConnection(info.actin_idx, j);
            for (auto& temp_conn : actinIndicesPerMyosin_temp) {
                temp_conn.deleteConnection(j, info.actin_idx);
            }
            info.active = false;
        };

        auto cmp_priority = [](const BoundInfo* lhs, const BoundInfo* rhs) {
            auto score = [](const BoundInfo* entry) {
                int value = 0;
                value += 4 * entry->motor_state;
                if (entry->cb) value += 2;
                if (entry->prev) value += 1;
                double align = std::abs(entry->dot);
                return std::make_tuple(value, align, -entry->actin_idx);
            };
            auto l = score(lhs);
            auto r = score(rhs);
            if (std::get<0>(l) != std::get<0>(r)) {
                return std::get<0>(l) > std::get<0>(r);
            }
            if (std::get<1>(l) != std::get<1>(r)) {
                return std::get<1>(l) > std::get<1>(r);
            }
            return std::get<2>(l) > std::get<2>(r);
        };

        int active_count = static_cast<int>(bound.size());

        if (active_count > max_myosin_bonds) {
            for (auto& info : bound) {
                if (active_count <= max_myosin_bonds) {
                    break;
                }
                if (!info.active || info.cb) {
                    continue;
                }
                remove_connection(info);
                --active_count;
            }
        }

        int per_group_limit = 0;
        if (max_myosin_bonds >= 2) {
            per_group_limit = max_myosin_bonds / 2;
        } else {
            per_group_limit = max_myosin_bonds;
        }

        // The legacy pair-owned model needed a polarity-group cap to avoid
        // duplicating an unlimited number of bonds. Conserved myosin-half
        // head pools already impose that capacity explicitly and may
        // distribute physical heads among all accessible neighboring actins.
        if (!dynamic_am && directional && per_group_limit > 0) {
            std::vector<BoundInfo*> positive;
            std::vector<BoundInfo*> negative;
            positive.reserve(bound.size());
            negative.reserve(bound.size());
            for (auto& info : bound) {
                if (!info.active) continue;
                if (info.dot >= 0.0) {
                    positive.push_back(&info);
                } else {
                    negative.push_back(&info);
                }
            }

            auto trim_group = [&](std::vector<BoundInfo*>& group) {
                if (static_cast<int>(group.size()) <= per_group_limit) {
                    return;
                }
                std::sort(group.begin(), group.end(), cmp_priority);
                for (size_t idx = per_group_limit; idx < group.size(); ++idx) {
                    if (!group[idx]->active) {
                        continue;
                    }
                    remove_connection(*group[idx]);
                    --active_count;
                }
            };

            trim_group(positive);
            trim_group(negative);
        }

        if (active_count > max_myosin_bonds) {
            std::vector<BoundInfo*> survivors;
            survivors.reserve(bound.size());
            for (auto& info : bound) {
                if (info.active) {
                    survivors.push_back(&info);
                }
            }
            std::sort(survivors.begin(), survivors.end(), cmp_priority);
            for (size_t idx = max_myosin_bonds; idx < survivors.size(); ++idx) {
                remove_connection(*survivors[idx]);
            }
        }
    }
}

void Sarcomere::_calc_am_force_velocity(int& i) {
    int thread_id = omp_get_thread_num();
    // Thread-local temporary lists for forces and velocities
    auto& local_actin_forces = actin_forces_temp[thread_id];
    auto& local_actin_torques = actin_torques_temp[thread_id];
    auto& local_myosin_forces = myosin_forces_temp[thread_id];
    auto& local_myosin_f_load = myosin_f_load_temp[thread_id];
    auto& local_myosin_torques = myosin_torques_temp[thread_id];
    auto& local_myosin_velocities = myosin_velocities_temp[thread_id];
    // getConnections() already returns a const reference and nothing below
    // mutates the connection list, so binding it avoids a per-actin heap copy.
    const std::vector<int>& myosin_indices =
        myosinIndicesPerActin.getConnections(i);
    const std::vector<int>& myosin_neighbors = actin_neighbors_by_species[i].second;
    vec velocity = v_am * actin.direction[i];
    const double interaction_ramp = _current_interaction_ramp();
    // Hoist string-keyed feature lookups out of the per-neighbor / per-head
    // loops below.  The keys themselves are constants, so they are built once
    // rather than allocated afresh for every actin of every step.
    static const std::string key_motor_force{"am_motor_force"};
    static const std::string key_active_force_x{"am_active_force_x"};
    static const std::string key_active_force_y{"am_active_force_y"};
    static const std::string key_active_force_z{"am_active_force_z"};
    static const std::string key_tether_force{"am_tether_force"};
    static const std::string key_head_transverse_force{"am_head_transverse_force"};
    auto& am_motor_force_feature = actin[key_motor_force];
    auto& am_active_force_x_feature = actin[key_active_force_x];
    auto& am_active_force_y_feature = actin[key_active_force_y];
    auto& am_active_force_z_feature = actin[key_active_force_z];
    auto& am_tether_force_feature = actin[key_tether_force];
    auto& am_head_transverse_force_feature = actin[key_head_transverse_force];
    // The half-suffixed myosin alignment diagnostics were looked up once per
    // bound head, each time building a fresh std::string.  Resolve both halves
    // once per actin instead.
    std::array<std::vector<double>*, 2> align_drive_feature{nullptr, nullptr};
    std::array<std::vector<double>*, 2> align_abs_feature{nullptr, nullptr};
    std::array<std::vector<double>*, 2> align_active_feature{nullptr, nullptr};
    std::array<std::vector<double>*, 2> align_tether_feature{nullptr, nullptr};
    std::array<std::vector<double>*, 2> align_strong_feature{nullptr, nullptr};
    std::array<std::vector<double>*, 2> align_count_feature{nullptr, nullptr};
    std::array<std::vector<double>*, 2> align_positive_feature{nullptr, nullptr};
    std::array<std::vector<double>*, 2> align_negative_feature{nullptr, nullptr};
    if (am_offcenter_torque) {
        // Both halves of all eight diagnostics are fixed names.  Concatenating
        // them per actin allocated and freed sixteen strings on every actin of
        // every step; the lookups below are now pure hash-table probes.
        static const std::array<std::array<std::string, 2>, 8> align_keys{{
            {{"am_align_drive_minus", "am_align_drive_plus"}},
            {{"am_align_abs_drive_minus", "am_align_abs_drive_plus"}},
            {{"am_align_active_drive_minus", "am_align_active_drive_plus"}},
            {{"am_align_tether_drive_minus", "am_align_tether_drive_plus"}},
            {{"am_align_strong_drive_minus", "am_align_strong_drive_plus"}},
            {{"am_align_head_count_minus", "am_align_head_count_plus"}},
            {{"am_align_positive_count_minus", "am_align_positive_count_plus"}},
            {{"am_align_negative_count_minus", "am_align_negative_count_plus"}},
        }};
        for (int half = 0; half < 2; ++half) {
            const auto h = static_cast<size_t>(half);
            align_drive_feature[half] = &myosin[align_keys[0][h]];
            align_abs_feature[half] = &myosin[align_keys[1][h]];
            align_active_feature[half] = &myosin[align_keys[2][h]];
            align_tether_feature[half] = &myosin[align_keys[3][h]];
            align_strong_feature[half] = &myosin[align_keys[4][h]];
            align_count_feature[half] = &myosin[align_keys[5][h]];
            align_positive_feature[half] = &myosin[align_keys[6][h]];
            align_negative_feature[half] = &myosin[align_keys[7][h]];
        }
    }
    for (int index = 0; index < myosin_indices.size(); index++) {
        int j = myosin_indices[index];
        if (am_bonds[i][j] != 1) {
            continue;
        }
        //scale by angle between actin and myosin
        double abs_cos_angle = std::abs(actin.direction[i].dot(myosin.direction[j]));
        double active_ratio = directional ? am_interaction[i][j].partial_binding_ratio
                                          : am_interaction[i][j].myosin_binding_ratio;
        double normalized_partial_ratio = std::min(
            active_ratio / geometry::am_head_bearing_fraction, 1.0);
        //apply exponential scaling to partial binding ratio
        normalized_partial_ratio = (1 - std::exp(-2 * normalized_partial_ratio)) / (1 - std::exp(-2));
        const bool strong = !dynamic_am ||
            am_state[i][j] == static_cast<BondFlag>(am_motor_kinetics::State::strong);
        const double state_stiffness = strong ? 1.0 : am_weak_stiffness_ratio;
        vec am_force_on_actin {0.0, 0.0, 0.0};
        vector force_vec;
        if (dynamic_am) {
            int bound_head_count = 0;
            for (const auto& half_heads :
                 am_heads_by_myosin[static_cast<size_t>(j)]) {
                for (const auto& head : half_heads) {
                    if (head.bound_actin == i &&
                        head.state != static_cast<BondFlag>(
                            am_motor_kinetics::State::detached) &&
                        head.state != static_cast<BondFlag>(
                            am_motor_kinetics::State::off_srx)) {
                        ++bound_head_count;
                    }
                }
            }
            // Spatial heads carry axial cross-bridge mechanics, optional
            // passive transverse compliance, and reciprocal off-center torque.
            // The aggregate pair tether below remains a separate coarse-grained
            // centerline containment interaction in both comparison arms.
            const double per_head_tether_k =
                interaction_ramp * am_head_transverse_k;
            for (int half_index = 0; half_index < 2; ++half_index) {
                const auto& heads =
                    am_heads_by_myosin[static_cast<size_t>(j)]
                                       [static_cast<size_t>(half_index)];
                for (int h = 0; h < static_cast<int>(heads.size()); ++h) {
                    const auto& head = heads[static_cast<size_t>(h)];
                    if (head.bound_actin != i) continue;
                    const auto state = static_cast<am_motor_kinetics::State>(head.state);
                    if (state == am_motor_kinetics::State::detached) continue;
                    const bool head_is_strong =
                        state == am_motor_kinetics::State::strong;
                    const double axial_stiffness = interaction_ramp *
                        am_motor_k * (head_is_strong
                            ? 1.0 : am_weak_stiffness_ratio);
                    const auto head_mechanics = compute_am_two_spring_head(
                        actin, myosin, i, j, box,
                        head.target_axial, head.anchor_axial,
                        head.binding_azimuth,
                        am_spatial_heads ? myosin.radius : 0.0,
                        am_spatial_heads
                            ? &myosin_material_e1[static_cast<size_t>(j)] : nullptr,
                        am_spatial_heads
                            ? &myosin_material_e2[static_cast<size_t>(j)] : nullptr,
                        actin_radius + am_surface_gap,
                        per_head_tether_k, am_tether_slack,
                        axial_stiffness, head.stroke_progress,
                        am_stall_force, head_is_strong,
                        am_bidirectional_crossbridge);
                    // Legacy sliding throttle: the prescribed-velocity model slowed
                    // filament sliding as the engaged AM overlap saturated
                    //   f_load = (1-exp(-2r))/(1-exp(-2)),  v = (1-f_load)*v0
                    // so sliding stopped once the optimal overlap was reached.  The
                    // explicit-crossbridge rewrite dropped it, which lets thin
                    // filaments slide through to the M-line.  Restore it here as a
                    // scale on the force-producing (active) term.
                    double overlap_throttle = 1.0;
                    if (am_overlap_target_ratio > 0.0) {
                        const double r = std::clamp(
                            active_ratio / am_overlap_target_ratio, 0.0, 1.0);
                        const double f_load =
                            (1.0 - std::exp(-2.0 * r)) / (1.0 - std::exp(-2.0));
                        overlap_throttle = std::clamp(1.0 - f_load, 0.0, 1.0);
                    }
                    const double active_scalar = head_mechanics.active_force;
                    const vec active_force = head_is_strong
                        ? head_mechanics.axial_force_on_actin
                        : vec{0.0, 0.0, 0.0};
                    const vec head_force_on_actin =
                        head_mechanics.force_on_actin;
                    am_force_on_actin += head_force_on_actin;

                    if (am_offcenter_torque) {
                        vec actin_axis = actin.direction[i];
                        vec myosin_axis = myosin.direction[j];
                        actin_axis.normalize();
                        myosin_axis.normalize();
                        const vec actin_orientation_drive =
                            head_mechanics.physical_torque_on_actin.cross(
                                actin_axis);
                        const vec myosin_orientation_drive =
                            head_mechanics.physical_torque_on_myosin.cross(
                                myosin_axis);
                        const double signed_drive =
                            project_nematic_alignment_drive(
                                myosin.direction[j], actin.direction[i],
                                myosin_orientation_drive);
                        const double active_drive = head_is_strong
                            ? signed_drive : 0.0;
                        const double tether_drive = head_is_strong
                            ? 0.0 : signed_drive;
                        auto& drive_feature = *align_drive_feature[half_index];
                        auto& abs_feature = *align_abs_feature[half_index];
                        auto& active_feature = *align_active_feature[half_index];
                        auto& tether_feature = *align_tether_feature[half_index];
                        auto& strong_feature = *align_strong_feature[half_index];
                        auto& count_feature = *align_count_feature[half_index];
                        auto& positive_feature = *align_positive_feature[half_index];
                        auto& negative_feature = *align_negative_feature[half_index];
                        #pragma omp atomic update
                        drive_feature[j] += signed_drive;
                        #pragma omp atomic update
                        abs_feature[j] += std::abs(signed_drive);
                        #pragma omp atomic update
                        active_feature[j] += active_drive;
                        #pragma omp atomic update
                        tether_feature[j] += tether_drive;
                        if (state == am_motor_kinetics::State::strong) {
                            #pragma omp atomic update
                            strong_feature[j] += signed_drive;
                        }
                        #pragma omp atomic update
                        count_feature[j] += 1.0;
                        if (signed_drive > 0.0) {
                            #pragma omp atomic update
                            positive_feature[j] += 1.0;
                        } else if (signed_drive < 0.0) {
                            #pragma omp atomic update
                            negative_feature[j] += 1.0;
                        }
                        local_actin_torques[i] += actin_orientation_drive;
                        const vec axial_drive =
                            head_mechanics.physical_torque_on_myosin.dot(
                                myosin_axis) *
                            myosin_axis;
                        // The director integrator consumes tau x u for tumble;
                        // the transported material frame additionally consumes
                        // the physical axial torque for roll.
                        local_myosin_torques[j] +=
                            myosin_orientation_drive + axial_drive;
                    }
                    am_motor_force_feature[i] += active_scalar;
                    am_active_force_x_feature[i] += active_force.x;
                    am_active_force_y_feature[i] += active_force.y;
                    am_active_force_z_feature[i] += active_force.z;
                    am_head_transverse_force_feature[i] +=
                        head_mechanics.radial_force_on_actin.norm();
                }
            }
            // Explicit heads supply active physical moments. A real occupied
            // cross-bridge ensemble additionally resists large axial
            // misalignment even though individual heads swivel. Represent
            // that unresolved passive compliance with an occupancy-scaled,
            // force-free angular band.
            const double angular_occupancy = std::min(
                1.0, static_cast<double>(bound_head_count) /
                         static_cast<double>(
                             std::max(1, am_tether_saturation_heads)));
            // A constitutive share of the radial containment, independent of
            // whether any head is currently attached.  The occupancy gate above
            // is keyed to bound heads and therefore vanishes exactly when they
            // let go, which is when an unconstrained pair drifts apart.  This
            // term is a pure centre force with no torque (see the accumulation
            // below), so it constrains separation without touching orientation.
            //
            // Reach-gated: am_tether_force_max defaults to 0, meaning an
            // unbounded quadratic, and the neighbour list is effectively
            // all-pairs, so an ungated floor would tether every pair in the box.
            // Guarded so the default (floor 0) costs nothing and is bitwise inert.
            double radial_occupancy = angular_occupancy;
            if (am_tether_occupancy_floor > 0.0) {
                const double pair_separation = geometry::segment_segment_distance(
                    actin.left_end[i], actin.right_end[i],
                    myosin.left_end[j], myosin.right_end[j],
                    box, actin.periodic_axes);
                if (pair_separation <= am_preferred_centerline +
                                       std::max(0.0, am_tether_occupancy_floor_reach)) {
                    radial_occupancy = std::max(
                        radial_occupancy,
                        std::min(1.0, am_tether_occupancy_floor));
                }
            }
            const auto radial_tether = compute_am_head_tether(
                actin, myosin, i, j, box, 0.0,
                interaction_ramp * am_tether_k * radial_occupancy,
                am_preferred_centerline, 0.0, 0.0,
                nullptr, nullptr, actin_radius, am_tether_slack,
                interaction_ramp * am_tether_force_max);
            am_force_on_actin += radial_tether.force_on_actin;
            am_tether_force_feature[i] +=
                radial_tether.force_on_actin.norm();
            const auto angular_tether =
                compute_am_nematic_flat_bottom_tether(
                    actin.direction[i], myosin.direction[j],
                    interaction_ramp * am_alignment_k * angular_occupancy,
                    am_alignment_slack);
            local_actin_torques[i] += angular_tether.drive_on_actin;
            local_myosin_torques[j] += angular_tether.drive_on_myosin;
        } else {
            force_vec = compute_am_force_and_energy(
                actin, myosin, i, j, box,
                interaction_ramp * k_am * normalized_partial_ratio * state_stiffness,
                interaction_ramp * kappa_am * state_stiffness,
                am_preferred_centerline + am_cutoff,
                am_preferred_centerline);
            am_force_on_actin = {force_vec[0], force_vec[1], force_vec[2]};
            local_actin_torques[i].x += force_vec[3];
            local_actin_torques[i].y += force_vec[4];
            local_actin_torques[i].z += force_vec[5];
            local_myosin_torques[j].x += force_vec[6];
            local_myosin_torques[j].y += force_vec[7];
            local_myosin_torques[j].z += force_vec[8];
        }
        local_actin_forces[i] += am_force_on_actin;
        local_myosin_forces[j] -= am_force_on_actin;

        if (!dynamic_am && actin.cb_status[i] == 2) {
            double f_load_am = (1 - std::exp(-2 * normalized_partial_ratio)) / (1 - std::exp(-2));
            vec delta_velocity = -(1 - f_load_am) * velocity;
            double dot = -velocity.dot(myosin.direction[j]);
            size_t dir_idx;
            dir_idx = (dot >= 0.0) ? 0 : 1;
            // printf("Actin %d (CB) and Myosin %d: partial ratio %f,  f_load_am %f, dir_idx %d, actin center (%f, %f, %f), actin direction (%f, %f, %f)\n",
            //        i, j, am_interaction[i][j].partial_binding_ratio, f_load_am, dir_idx,
            //          actin.center[i].x, actin.center[i].y, actin.center[i].z,
            //             actin.direction[i].x, actin.direction[i].y, actin.direction[i].z);
            if (f_load_am > local_myosin_f_load[j][dir_idx]) {
                local_myosin_f_load[j][dir_idx] = f_load_am;
            }
            local_myosin_velocities[j] += delta_velocity;
        }
    }

    if (dynamic_am) {
        // Dynamic motors move the particles through their explicit axial spring;
        // no additional prescribed sliding velocity is applied.
        actin.velocity[i] = {0, 0, 0};
    } else if (actin.cb_status[i] == 2) {
        if (myosin_indices.size() == 0) {
        }
        // Keep the prescribed-velocity compatibility path under the same
        // explicit overlap-throttle control as the dynamic-head path.  A
        // non-positive target disables the empirical overlap arrest; a
        // positive target retains the historical status-2 slowdown.
        const double overlap_throttle = am_overlap_target_ratio > 0.0
            ? std::clamp(1.0 - actin.f_load[i], 0.0, 1.0)
            : 1.0;
        actin.velocity[i] = overlap_throttle * velocity;
        // printf("Actin %d (CB): f_load %f, velocity set to (%f, %f, %f)\n",
        //        i, actin.f_load[i],
        //        actin.velocity[i].x,
        //        actin.velocity[i].y,
        //        actin.velocity[i].z);
        if (has_myosin_bond_pairs) {
            for (int j : myosin_neighbors) {
                if (am_bonds[i][j] == 1) {
                    continue;
                }
                bool bonded_pair = false;
                for (int bound_idx : myosin_indices) {
                    if (am_interaction[i][bound_idx].partial_binding_ratio <= EPS) {
                        continue;
                    }
                    if (_myosin_pair_bonded(bound_idx, j)) {
                        bonded_pair = true;
                        break;
                    }
                }
                if (!bonded_pair) {
                    continue;
                }
                apply_actin_myosin_repulsion(
                    actin,
                    myosin,
                    i,
                    j,
                    box,
                    am_cutoff,
                    myosin_steric_k,
                    local_actin_forces[i],
                    local_myosin_forces[j]);
            }
        }
    } else {
        if (myosin_indices.size() == 0) {
            actin.velocity[i] = {0, 0, 0};
        } else {
            actin.velocity[i] = velocity;
        }
    }
}


int Sarcomere::_best_explicit_titin_candidate(
    size_t titin_index, const vec& thick_anchor) const {
    if (titin_index >= explicit_titins.size()) return -1;
    const auto& titin = explicit_titins[titin_index];
    // Same-material-partner recapture. This branch is deliberately exempt from
    // explicit_titin_require_live_actn2: the identity and the short memory
    // window are themselves the specificity, so momentary ACTN2 turnover under
    // an established titin lets the complex anneal instead of dissolving.
    if (explicit_titin_rebind_memory > 0.0 &&
        titin.last_bound_actin >= 0 &&
        titin.last_bound_actin < actin.n &&
        explicit_titin_owner_by_actin[
            static_cast<size_t>(titin.last_bound_actin)] < 0 &&
        (!explicit_titin_rebind_requires_actn2 ||
         actin_has_live_actn2(titin.last_bound_actin))) {
        const double elapsed = dt * static_cast<double>(
            current_step >= titin.last_detach_step
                ? current_step - titin.last_detach_step : 0);
        if (elapsed <= explicit_titin_rebind_memory) {
            vec displacement = thick_anchor -
                static_cast<vec>(actin.left_end[titin.last_bound_actin]);
            displacement.pbc_wrap(box, is_periodic);
            if (displacement.norm() <= explicit_titin_capture_distance) {
                return titin.last_bound_actin;
            }
        }
    }
    vec myosin_axis = myosin.direction[titin.myosin];
    myosin_axis.normalize();
    int best_actin = -1;
    double best_distance = explicit_titin_capture_distance;
    const auto& neighbors =
        neighbor_list.get_neighbors(titin.myosin + actin.n);
    for (const auto& [neighbor_index, neighbor_type] : neighbors) {
        if (neighbor_type != ParticleType::Actin) continue;
        const int a = neighbor_index;
        if (a < 0 || a >= actin.n) continue;
        if (explicit_titin_owner_by_actin[static_cast<size_t>(a)] >= 0 &&
            explicit_titin_owner_by_actin[static_cast<size_t>(a)] !=
                static_cast<int>(titin_index)) {
            continue;
        }
        // Ordinary new-partner search. Unlike the same-partner recapture above,
        // it must find local ACTN2 support on this individual barbed end.
        if (explicit_titin_require_live_actn2 && !actin_has_live_actn2(a)) {
            continue;
        }
        vec from_myosin = static_cast<vec>(actin.left_end[a]) -
                          static_cast<vec>(myosin.center[titin.myosin]);
        from_myosin.pbc_wrap(box, is_periodic);
        const double half_projection = from_myosin.dot(myosin_axis);
        if ((titin.half == 0 && half_projection >= 0.0) ||
            (titin.half == 1 && half_projection <= 0.0)) {
            continue;
        }
        vec displacement = thick_anchor -
                           static_cast<vec>(actin.left_end[a]);
        displacement.pbc_wrap(box, is_periodic);
        vec actin_axis = actin.direction[a];
        actin_axis.normalize();
        if (actin_axis.dot(displacement) <= 0.0) continue;
        const double distance = displacement.norm();
        if (distance > best_distance) continue;
        if (best_actin < 0 || distance < best_distance ||
            (distance == best_distance && a < best_actin)) {
            best_actin = a;
            best_distance = distance;
        }
    }
    return best_actin;
}

void Sarcomere::_update_explicit_titin() {
    if (titin_mode != TitinMode::explicit_dynamic) return;
    if (!(titin_k > 0.0) || !(titin_rest_length > 0.0) ||
        !(explicit_titin_capture_distance > 0.0) ||
        explicit_titin_k_on < 0.0 || explicit_titin_k_off < 0.0 ||
        explicit_titin_transition_theta < 0.0 ||
        explicit_titin_transition_theta > 1.0 ||
        explicit_titin_rebind_memory < 0.0 ||
        !std::isfinite(explicit_titin_rebind_memory) ||
        explicit_titin_rebind_rate_multiplier < 0.0 ||
        !std::isfinite(explicit_titin_rebind_rate_multiplier) ||
        !std::isfinite(explicit_titin_k_off_unsupported) ||
        (explicit_titin_strain_gated_kinetics && !(thermal_beta > 0.0))) {
        throw std::invalid_argument("invalid explicit titin configuration");
    }
    if (explicit_titins.empty()) initialize_explicit_titins();
    if (explicit_titin_owner_by_actin.size() !=
        static_cast<size_t>(actin.n)) {
        throw std::runtime_error("invalid explicit titin ownership state");
    }

    auto& actin_bound = actin["titin_bound"];
    auto& actin_force_feature = actin["titin_force"];
    auto& bound_minus = myosin["titin_bound_minus"];
    auto& bound_plus = myosin["titin_bound_plus"];
    auto& force_minus = myosin["titin_force_minus"];
    auto& force_plus = myosin["titin_force_plus"];
    std::fill(actin_bound.begin(), actin_bound.end(), 0.0);
    std::fill(actin_force_feature.begin(), actin_force_feature.end(), 0.0);
    std::fill(bound_minus.begin(), bound_minus.end(), 0.0);
    std::fill(bound_plus.begin(), bound_plus.end(), 0.0);
    std::fill(force_minus.begin(), force_minus.end(), 0.0);
    std::fill(force_plus.begin(), force_plus.end(), 0.0);

    auto record_event = [&](int event_code, size_t titin_index, int a,
                            double distance, double elastic_energy,
                            double force_value) {
        const auto& titin = explicit_titins[titin_index];
        const std::array<double, 10> row{{
            static_cast<double>(current_step),
            static_cast<double>(event_code),
            static_cast<double>(titin_index),
            static_cast<double>(titin.myosin),
            static_cast<double>(titin.half),
            static_cast<double>(titin.track),
            static_cast<double>(a), distance, elastic_energy, force_value}};
        explicit_titin_events.insert(
            explicit_titin_events.end(), row.begin(), row.end());
    };

    for (size_t index = 0; index < explicit_titins.size(); ++index) {
        auto& titin = explicit_titins[index];
        const vec thick_anchor = explicit_titin_thick_anchor(titin);
        if (!titin.bound) {
            titin.distance = 0.0;
            titin.extension = 0.0;
            titin.force = 0.0;
            titin.candidate_distance = 0.0;
            titin.candidate_energy = 0.0;
            titin.candidate_on_rate = 0.0;
            titin.candidate_actin =
                _best_explicit_titin_candidate(index, thick_anchor);
            if (titin.candidate_actin < 0) continue;
            vec candidate_displacement = thick_anchor -
                static_cast<vec>(actin.left_end[titin.candidate_actin]);
            candidate_displacement.pbc_wrap(box, is_periodic);
            titin.candidate_distance = candidate_displacement.norm();
            titin.candidate_energy = titin_kinetics::energy(
                titin.candidate_distance, titin_rest_length, titin_k);
            double candidate_base_rate = explicit_titin_k_on;
            const bool same_partner_recapture =
                explicit_titin_rebind_memory > 0.0 &&
                titin.candidate_actin == titin.last_bound_actin &&
                current_step >= titin.last_detach_step &&
                dt * static_cast<double>(
                    current_step - titin.last_detach_step) <=
                    explicit_titin_rebind_memory;
            if (same_partner_recapture) {
                candidate_base_rate *=
                    explicit_titin_rebind_rate_multiplier;
            }
            // Capture is a local recognition event.  Do not charge the elastic
            // energy of the fully formed tension-bearing spring before the
            // molecular association exists.  Once captured, the existing
            // tension-only force law and (optionally) strain-dependent off-rate
            // govern the bound state.
            titin.candidate_on_rate =
                titin_kinetics::distance_limited_association_rate(
                    candidate_base_rate, titin.candidate_distance,
                    explicit_titin_capture_distance);
            counter_rng::CounterRng binding_rng(
                static_cast<std::uint64_t>(initial_seed),
                static_cast<std::uint64_t>(current_step),
                counter_rng::Stream::explicit_titin_binding,
                static_cast<std::uint64_t>(index),
                static_cast<std::uint64_t>(titin.candidate_actin));
            if (binding_rng.uniform() >= titin_kinetics::transition_probability(
                    titin.candidate_on_rate, dt)) {
                continue;
            }
            const int a = titin.candidate_actin;
            if (explicit_titin_owner_by_actin[static_cast<size_t>(a)] >= 0) {
                continue;
            }
            // Provenance of this association, recorded before any mechanics.
            // Recaptures are counted on their own because they bypass the gate
            // by design and must not be read as gate leakage.
            if (same_partner_recapture) {
                ++explicit_titin_rebind_associations;
            } else if (actin_has_live_actn2(a)) {
                ++explicit_titin_supported_associations;
            } else {
                ++explicit_titin_unsupported_associations;
            }
            titin.bound = true;
            titin.bound_actin = a;
            titin.last_bound_actin = a;
            titin.candidate_actin = -1;
            titin.attach_step = static_cast<RecoveryStep>(current_step);
            titin.attachment_distance = titin.candidate_distance;
            titin.attachment_energy = titin.candidate_energy;
            titin.attachment_force = titin_kinetics::force(
                titin.attachment_distance, titin_rest_length, titin_k);
            explicit_titin_owner_by_actin[static_cast<size_t>(a)] =
                static_cast<int>(index);
            record_event(1, index, a, titin.attachment_distance,
                         titin.attachment_energy, titin.attachment_force);
        }

        const int a = titin.bound_actin;
        if (a < 0 || a >= actin.n ||
            explicit_titin_owner_by_actin[static_cast<size_t>(a)] !=
                static_cast<int>(index)) {
            throw std::runtime_error("invalid explicit titin binding state");
        }
        const vec actin_anchor = actin.left_end[a];
        vec displacement = thick_anchor - actin_anchor;
        displacement.pbc_wrap(box, is_periodic);
        titin.distance = displacement.norm();
        titin.extension = titin_kinetics::extension(
            titin.distance, titin_rest_length);
        titin.force = titin_kinetics::force(
            titin.distance, titin_rest_length, titin_k);
        const double elastic_energy = titin_kinetics::energy(
            titin.distance, titin_rest_length, titin_k);

        // Serial anchor: actin -- ACTN2 -- titin. With no live ACTN2 on this
        // barbed end the composite link is only as good as its missing middle
        // element, so dissociation runs at the unsupported rate. This is still
        // a stochastic, strain-gated release rather than a forced detachment,
        // and the same-partner rebinding memory can carry the molecule across
        // a momentary gap.
        const bool anchor_supported = actin_has_live_actn2(a);
        const double base_off_rate =
            (explicit_titin_k_off_unsupported > 0.0 && !anchor_supported)
                ? explicit_titin_k_off_unsupported
                : explicit_titin_k_off;
        double off_rate = base_off_rate;
        if (explicit_titin_strain_gated_kinetics) {
            off_rate = titin_kinetics::dissociation_rate(
                base_off_rate, thermal_beta,
                explicit_titin_transition_theta, elastic_energy);
        }
        counter_rng::CounterRng release_rng(
            static_cast<std::uint64_t>(initial_seed),
            static_cast<std::uint64_t>(current_step),
            counter_rng::Stream::explicit_titin_release,
            static_cast<std::uint64_t>(index),
            static_cast<std::uint64_t>(a));
        if (release_rng.uniform() < titin_kinetics::transition_probability(
                off_rate, dt)) {
            // Event code 3 marks a release that happened while the serial
            // anchor was absent, so anchor-loss releases stay separable from
            // titin's own chemistry in the event stream.
            record_event(anchor_supported ? 2 : 3, index, a, titin.distance,
                         elastic_energy, titin.force);
            if (!anchor_supported) ++explicit_titin_unsupported_releases;
            explicit_titin_owner_by_actin[static_cast<size_t>(a)] = -1;
            titin.bound = false;
            titin.last_bound_actin = a;
            titin.bound_actin = -1;
            titin.candidate_actin = -1;
            titin.last_detach_step = static_cast<RecoveryStep>(current_step);
            titin.distance = 0.0;
            titin.extension = 0.0;
            titin.force = 0.0;
            continue;
        }

        actin_bound[static_cast<size_t>(a)] = 1.0;
        actin_force_feature[static_cast<size_t>(a)] = titin.force;
        auto& half_bound = titin.half == 0 ? bound_minus : bound_plus;
        auto& half_force = titin.half == 0 ? force_minus : force_plus;
        half_bound[static_cast<size_t>(titin.myosin)] += 1.0;
        half_force[static_cast<size_t>(titin.myosin)] += titin.force;
        if (!(titin.extension > 0.0) || titin.distance <= 1e-12) continue;

        const vec force_on_actin =
            (titin.force / titin.distance) * displacement;
        const vec force_on_myosin = -force_on_actin;
        actin.force[a] += force_on_actin;
        myosin.force[titin.myosin] += force_on_myosin;
        const auto actin_torque = compute_am_offcenter_torque(
            actin.center[a], actin.direction[a], actin_anchor,
            force_on_actin, box, is_periodic);
        actin.torque[a] += actin_torque.orientation_drive;
        const auto myosin_torque = compute_am_offcenter_torque(
            myosin.center[titin.myosin], myosin.direction[titin.myosin],
            thick_anchor, force_on_myosin, box, is_periodic);
        vec myosin_axis = myosin.direction[titin.myosin];
        myosin_axis.normalize();
        const vec axial_roll =
            myosin_torque.physical_torque.dot(myosin_axis) * myosin_axis;
        myosin.torque[titin.myosin] +=
            myosin_torque.orientation_drive + axial_roll;
    }
}

void Sarcomere::_apply_titin_forces(int& i) {
    int thread_id = omp_get_thread_num();
    auto& local_actin_forces = actin_forces_temp[thread_id];
    auto& local_myosin_forces = myosin_forces_temp[thread_id];

    if (titin_mode == TitinMode::balanced_owning) {
        // The minimal balanced cell has two periodic Z-discs. A0/A2 belong
        // to M0 and A1/A3 belong to M1. This spring is always present and
        // reciprocal; it does not disappear when either AM heads or ACTN2
        // molecules turn over.
        constexpr std::array<int, 4> owning_myosin{{0, 1, 0, 1}};
        const int m = owning_myosin[static_cast<size_t>(i)];
        vec diff = myosin.center[m] - actin.left_end[i];
        diff.pbc_wrap(box, is_periodic);
        const double dist = diff.norm();
        if (dist > EPS) {
            const vec force =
                (titin_k * (dist - titin_rest_length) / dist) * diff;
            local_actin_forces[i] += force;
            local_myosin_forces[m] -= force;
        }
        return;
    }

    for (int j = i + 1; j < actin.n; ++j) {
        if (actin_actin_status[i][j] != 2) {
            continue;
        }

        auto apply_spring = [&](int act_a, int act_b) {
            const auto& myosin_indices = myosinIndicesPerActin.getConnections(act_b);
            vec anchor = actin.left_end[act_a];
            for (int m : myosin_indices) {
                if (am_bonds[act_b][m] != 1) {
                    continue;
                }
                vec myo_center = myosin.center[m];
                vec diff = myo_center - anchor;
                diff.pbc_wrap(box, is_periodic);
                double dist = diff.norm();
                if (dist < EPS) {
                    continue;
                }
                double stretch = dist - titin_rest_length;
                vec force = (titin_k * stretch / dist) * diff;
                if (force.norm_squared() < EPS) {
                    continue;
                }
                local_actin_forces[act_a] += force;
                local_myosin_forces[m] -= force;
                double dot_prod = actin.direction[act_a].dot(diff)/dist;
                if (stretch < 0) {
                    dot_prod = -dot_prod;
                }
                double am_dot_prod = actin.direction[act_b].dot(myosin.direction[m]);
            }
        };

        apply_spring(i, j);
        apply_spring(j, i);
    }
}


bool Sarcomere::_myomesin_pair_geometry(
    int i, int j, vec& axis, double& axial, vec& radial,
    double& radial_distance) const {
    if (i < 0 || j < 0 || i >= myosin.n || j >= myosin.n || i == j) {
        return false;
    }
    vec di = myosin.direction[i];
    vec dj = myosin.direction[j];
    if (di.norm() <= 1e-12 || dj.norm() <= 1e-12) return false;
    di.normalize();
    dj.normalize();
    if (di.dot(dj) < 0.0) dj = -dj;
    axis = di + dj;
    if (axis.norm() <= 1e-12) return false;
    axis.normalize();
    // Vector from i to j (pbc_diff_masked returns first minus second).
    const vec diff = utils::pbc_diff_masked(
        myosin.center[j], myosin.center[i], box, is_periodic);
    axial = diff.dot(axis);
    radial = diff - axial * axis;
    radial_distance = radial.norm();
    return std::isfinite(axial) && std::isfinite(radial_distance);
}

void Sarcomere::_update_myomesin_bridges() {
    if (k_mm <= 0.0) {
        myomesin_bridges.clear();
        return;
    }

    auto record_event = [&](int event_code, const MyomesinBridge& bridge) {
        const std::array<double, 9> row{{
            static_cast<double>(current_step),
            static_cast<double>(event_code),
            static_cast<double>(bridge.id),
            static_cast<double>(bridge.myosin_i),
            static_cast<double>(bridge.myosin_j),
            bridge.axial_offset, bridge.radial_distance, bridge.force,
            dt * static_cast<double>(current_step - bridge.attach_step)}};
        myomesin_events.insert(myomesin_events.end(), row.begin(), row.end());
    };

    // Bound molecules turn over independently.  The optional Bell distance
    // accelerates release under load but never forces an actomyosin head off.
    std::vector<MyomesinBridge> retained;
    retained.reserve(myomesin_bridges.size());
    for (auto bridge : myomesin_bridges) {
        vec axis, radial;
        double axial = 0.0;
        double radial_distance = 0.0;
        if (!_myomesin_pair_geometry(bridge.myosin_i, bridge.myosin_j,
                                     axis, axial, radial,
                                     radial_distance)) {
            record_event(3, bridge);
            continue;
        }
        const double radial_k = myomesin_radial_k > 0.0
            ? myomesin_radial_k : k_mm;
        const double axial_force = k_mm * axial;
        const double radial_force = radial_k *
            (radial_distance - myomesin_radial_rest);
        bridge.axial_offset = axial;
        bridge.radial_distance = radial_distance;
        bridge.force = std::hypot(axial_force, radial_force);
        double off_rate = myomesin_k_off;
        if (myomesin_bell_distance > 0.0 && thermal_beta > 0.0) {
            const double exponent = std::min(
                20.0, thermal_beta * myomesin_bell_distance * bridge.force);
            off_rate *= std::exp(exponent);
        }
        counter_rng::CounterRng release_rng(
            static_cast<std::uint64_t>(initial_seed),
            static_cast<std::uint64_t>(current_step),
            counter_rng::Stream::myomesin_release, bridge.id,
            counter_rng::canonical_pair_key(
                bridge.myosin_i, bridge.myosin_j));
        if (release_rng.uniform() < 1.0 - std::exp(-off_rate * dt)) {
            record_event(2, bridge);
            continue;
        }
        retained.push_back(bridge);
    }
    myomesin_bridges.swap(retained);

    if (!(myomesin_k_on > 0.0)) return;
    std::vector<int> degree(static_cast<size_t>(myosin.n), 0);
    std::unordered_set<std::uint64_t> occupied;
    for (const auto& bridge : myomesin_bridges) {
        ++degree[static_cast<size_t>(bridge.myosin_i)];
        ++degree[static_cast<size_t>(bridge.myosin_j)];
        occupied.insert(counter_rng::canonical_pair_key(
            bridge.myosin_i, bridge.myosin_j));
    }

    struct Candidate {
        int i = -1;
        int j = -1;
        double axial = 0.0;
        double radial = 0.0;
        double priority = 0.0;
    };
    std::vector<Candidate> accepted;
    const double probability = 1.0 - std::exp(-myomesin_k_on * dt);
    for (int i = 0; i < myosin.n; ++i) {
        for (const auto& [neighbor_index, neighbor_type] :
             neighbor_list.get_neighbors(i + actin.n)) {
            if (neighbor_type != ParticleType::Myosin) continue;
            const int j = neighbor_index - actin.n;
            if (i >= j) continue;
            const std::uint64_t pair_key =
                counter_rng::canonical_pair_key(i, j);
            if (occupied.count(pair_key) != 0) continue;
            vec axis, radial;
            double axial = 0.0;
            double radial_distance = 0.0;
            vec di = myosin.direction[i];
            vec dj = myosin.direction[j];
            if (di.norm() <= 1e-12 || dj.norm() <= 1e-12) continue;
            di.normalize();
            dj.normalize();
            if (!_myomesin_pair_geometry(i, j, axis, axial, radial,
                                         radial_distance) ||
                std::abs(di.dot(dj)) + 1e-12 < myomesin_alignment_cos ||
                std::abs(axial) > myomesin_axial_capture ||
                std::abs(radial_distance - myomesin_radial_rest) >
                    myomesin_radial_capture) {
                continue;
            }
            counter_rng::CounterRng binding_rng(
                static_cast<std::uint64_t>(initial_seed),
                static_cast<std::uint64_t>(current_step),
                counter_rng::Stream::myomesin_binding, pair_key);
            if (binding_rng.uniform() >= probability) continue;
            accepted.push_back(
                Candidate{i, j, axial, radial_distance,
                          binding_rng.uniform()});
        }
    }
    // Successful encounters compete in a random, counter-based order so a
    // saturated filament does not systematically favor the lowest pair ID.
    std::sort(accepted.begin(), accepted.end(),
              [](const Candidate& a, const Candidate& b) {
                  if (a.priority != b.priority) return a.priority < b.priority;
                  return std::tie(a.i, a.j) < std::tie(b.i, b.j);
              });
    for (const auto& candidate : accepted) {
        if (degree[static_cast<size_t>(candidate.i)] >=
                myomesin_sites_per_filament ||
            degree[static_cast<size_t>(candidate.j)] >=
                myomesin_sites_per_filament) {
            continue;
        }
        MyomesinBridge bridge;
        bridge.id = myomesin_next_id++;
        bridge.myosin_i = candidate.i;
        bridge.myosin_j = candidate.j;
        bridge.attach_step = static_cast<RecoveryStep>(current_step);
        bridge.axial_offset = candidate.axial;
        bridge.radial_distance = candidate.radial;
        myomesin_bridges.push_back(bridge);
        ++degree[static_cast<size_t>(candidate.i)];
        ++degree[static_cast<size_t>(candidate.j)];
        record_event(1, myomesin_bridges.back());
    }
}

void Sarcomere::_apply_myomesin_bridge(
    MyomesinBridge& bridge, std::vector<vec>& local_myosin_forces) {
    if (k_mm <= 0.0) return;
    vec axis, radial;
    double axial = 0.0;
    double radial_distance = 0.0;
    if (!_myomesin_pair_geometry(bridge.myosin_i, bridge.myosin_j,
                                 axis, axial, radial, radial_distance)) {
        return;
    }
    const double ramp = _current_interaction_ramp();
    const double radial_k = myomesin_radial_k > 0.0
        ? myomesin_radial_k : k_mm;
    vec force_on_i = ramp * k_mm * axial * axis;
    if (radial_distance > 1e-12) {
        force_on_i += ramp * radial_k *
            (radial_distance - myomesin_radial_rest) /
            radial_distance * radial;
    }
    double magnitude = force_on_i.norm();
    if (myomesin_force_max > 0.0 && magnitude > myomesin_force_max) {
        force_on_i = (myomesin_force_max / magnitude) * force_on_i;
        magnitude = myomesin_force_max;
    }
    bridge.axial_offset = axial;
    bridge.radial_distance = radial_distance;
    bridge.force = magnitude;
    local_myosin_forces[bridge.myosin_i] += force_on_i;
    local_myosin_forces[bridge.myosin_j] -= force_on_i;
}



// The lateral kernel fires only when the barbed-end separation has axial
// component within aa_barbed_lateral_axial_window AND transverse component
// below aa_barbed_lateral_range.  Those are orthogonal parts of the same
// vector, so |separation|^2 < window^2 + range^2 is a strict necessary
// condition -- and it needs no normalization, because the unit directions were
// cached in _set_to_zero using the kernel's own expression.
bool Sarcomere::_barbed_ends_within_reach(int i, int j) const {
    const vec barbed_i = static_cast<vec>(actin.center[i]) -
                         0.5 * actin.length * actin_unit_direction[
                             static_cast<size_t>(i)];
    const vec barbed_j = static_cast<vec>(actin.center[j]) -
                         0.5 * actin.length * actin_unit_direction[
                             static_cast<size_t>(j)];
    vec separation = barbed_j - barbed_i;
    separation.pbc_wrap(box, is_periodic);
    const double reach_squared =
        aa_barbed_lateral_axial_window * aa_barbed_lateral_axial_window +
        aa_barbed_lateral_range * aa_barbed_lateral_range;
    return separation.norm_squared() < reach_squared;
}

void Sarcomere::_volume_exclusion(){
    #pragma omp for schedule(static)
    for (size_t bridge_index = 0;
         bridge_index < myomesin_bridges.size(); ++bridge_index) {
        auto& local_myosin_forces =
            myosin_forces_temp[omp_get_thread_num()];
        _apply_myomesin_bridge(
            myomesin_bridges[bridge_index], local_myosin_forces);
    }
    // Myomesin is a molecular interaction, not a steric force.  Turning off
    // volume exclusion must not silently turn off already occupied bridges.
    if (!sterics_enabled) return;
    const double EPS_FORCE = 1e-9;
    double myosin_cutoff = 2.0 * myosin.radius;
    #pragma omp for schedule(runtime)
    for (int i = 0; i<myosin.n; i++){
        auto& local_myosin_forces = myosin_forces_temp[omp_get_thread_num()];
        auto& local_myosin_torques = myosin_torques_temp[omp_get_thread_num()];
        for (const auto& [neighbor_index, neighbor_type] :
             neighbor_list.get_neighbors(i + actin.n)) {
            if (neighbor_type != ParticleType::Myosin) continue;
            const int j = neighbor_index - actin.n;
            if (i<j){
                const auto repulsion = compute_myosin_repulsion(
                    actin,
                    myosin,
                    i,
                    j,
                    box,
                    fix_myosin,
                    actinIndicesPerMyosin,
                    myosin_steric_k);
                if (repulsion.applied) {
                    local_myosin_forces[i] += repulsion.force_on_first;
                    local_myosin_forces[j] += repulsion.force_on_second;
                    local_myosin_torques[i] += repulsion.torque_on_first;
                    local_myosin_torques[j] += repulsion.torque_on_second;
                }
            }
        }
    }

    #pragma omp for schedule(runtime)
    for (int i = 0; i < actin.n; i++){
        // Segment-distance-filtered subset: every interaction in this loop is
        // short range, so the ~857-entry centre-based list is ~4x more than
        // needed.  Dropped pairs would have contributed exactly zero.
        const auto& actin_indices =
            actin_close_neighbors[static_cast<size_t>(i)];
        auto& local_actin_forces = actin_forces_temp[omp_get_thread_num()];
        auto& local_actin_torques = actin_torques_temp[omp_get_thread_num()];
        for (int j : actin_indices) {
            if (i<j){
                RepulsionResult repulsion{};
                if (aa_double_overlap_contact_enabled) {
                    // Deliberately ungated.  This kernel is transverse-gated
                    // with no bound on the axial gap, so no distance test can
                    // bracket it; its own early-out is already the cheap path.
                    repulsion = compute_opposing_pointed_end_contact(
                        actin, i, j, box,
                        aa_double_overlap_contact_diameter,
                        aa_double_overlap_contact_epsilon,
                        aa_double_overlap_contact_force_cap);
                }
                // Retain ordinary capsule sterics for every non-end-on pair.
                // Do not add it to an active pointed-end contact: after an
                // attempted collinear crossing its zero-distance regularizer
                // chooses an arbitrary transverse direction and would mask
                // the physically resolved axial reaction above.
                if (!repulsion.applied &&
                    (!opt_gates_enabled ||
                     geometry::segments_may_interact(
                        actin.left_end[i], actin.right_end[i],
                        actin.left_end[j], actin.right_end[j],
                        2.0 * actin_radius, box, is_periodic))) {
                    // The kernel returns no force once the segment distance
                    // reaches 2*actin_radius, and the gate is a strict lower
                    // bound on that same distance, so a rejection here is
                    // exactly the zero result the call would have produced.
                    // This is what removes segment_segment_distance_w_normal
                    // from ~99% of pairs; perf had it at 15.3% of all cycles.
                    repulsion = compute_actin_repulsion(
                        actin, i, j, box,
                        2.0 * actin_radius,
                        actin_steric_k);
                }
                // Barbed-end lateral pressure.  Additive to the capsule
                // sterics: longer range, transverse only, and it is what
                // keeps a crosslinked cluster from collapsing to the hard
                // core.  Inert unless range and stiffness are both set.
                if (aa_barbed_lateral_range > 0.0 &&
                    aa_barbed_lateral_k > 0.0 &&
                    (!opt_gates_enabled ||
                     _barbed_ends_within_reach(i, j))) {
                    const auto lateral = compute_barbed_end_lateral_repulsion(
                        actin, i, j, box,
                        aa_barbed_lateral_range,
                        aa_barbed_lateral_k,
                        aa_barbed_lateral_force_max,
                        aa_barbed_lateral_axial_window);
                    if (lateral.applied) {
                        local_actin_forces[i] += lateral.force_on_first;
                        local_actin_forces[j] += lateral.force_on_second;
                        local_actin_torques[i] += lateral.torque_on_first;
                        local_actin_torques[j] += lateral.torque_on_second;
                    }
                }
                if (repulsion.applied) {
                    local_actin_forces[i] += repulsion.force_on_first;
                    local_actin_forces[j] += repulsion.force_on_second;
                    local_actin_torques[i] += repulsion.torque_on_first;
                    local_actin_torques[j] += repulsion.torque_on_second;
                    if (repulsion.contact_penetration > 0.0) {
                        #pragma omp atomic update
                        aa_double_overlap_contact_count += 1;
                        #pragma omp critical(double_overlap_contact_maxima)
                        {
                            aa_double_overlap_max_penetration = std::max(
                                aa_double_overlap_max_penetration,
                                repulsion.contact_penetration);
                            aa_double_overlap_max_force = std::max(
                                aa_double_overlap_max_force,
                                repulsion.force_on_first.norm());
                        }
                    }
                }
            }
        }
    }

    #pragma omp for schedule(runtime)
    for (int i = 0; i < actin.n; ++i) {
        auto& local_actin_forces = actin_forces_temp[omp_get_thread_num()];
        auto& local_myosin_forces = myosin_forces_temp[omp_get_thread_num()];
        auto& local_actin_torques = actin_torques_temp[omp_get_thread_num()];
        auto& local_myosin_torques = myosin_torques_temp[omp_get_thread_num()];
        for (int j :
             actin_neighbors_by_species[static_cast<size_t>(i)].second) {
            // Same unfiltered ~183-entry list as the binding loop.
            if (opt_gates_enabled &&
                !geometry::segments_may_interact(
                    actin.left_end[i], actin.right_end[i],
                    myosin.left_end[j], myosin.right_end[j],
                    myosin.radius + actin_radius, box, is_periodic)) {
                continue;
            }
            const auto repulsion =
                compute_actin_myosin_radial_repulsion_result(
                actin, myosin, i, j, box,
                myosin.radius + actin_radius, am_steric_k);
            if (repulsion.applied) {
                local_actin_forces[i] += repulsion.force_on_first;
                local_myosin_forces[j] += repulsion.force_on_second;
                local_actin_torques[i] += repulsion.torque_on_first;
                local_myosin_torques[j] += repulsion.torque_on_second;
            }
        }
    }
}

void Sarcomere::_myosin_exclusion(){
    const double EPS_FORCE = 1e-9;
    double myosin_cutoff = 2.0 * myosin.radius;
    #pragma omp for schedule(runtime)
    for (int i = 0; i<myosin.n; i++){
        auto& local_myosin_forces = myosin_forces_temp[omp_get_thread_num()];
        for (const auto& [neighbor_index, neighbor_type] :
             neighbor_list.get_neighbors(i + actin.n)) {
            if (neighbor_type != ParticleType::Myosin) continue;
            const int j = neighbor_index - actin.n;
            if (i<j){
                double seg_distance = std::numeric_limits<double>::infinity();
                apply_myosin_repulsion(
                    actin,
                    myosin,
                    i,
                    j,
                    box,
                    fix_myosin,
                    actinIndicesPerMyosin,
                    myosin_steric_k,
                    local_myosin_forces[i],
                    local_myosin_forces[j],
                    seg_distance,
                    myosin_exclusion_distance);
            }
        }
    }
}

actn2_load_kinetics::ArmLoads Sarcomere::_aa_signed_arm_loads(int i, int j) {
    // Catch-bond pairs are evaluated by _process_catch_bonds with i > j.  The
    // AA energy is pair-symmetric; retain canonical ordering so diagnostics
    // and classification follow the same numerical path.
    int first = std::max(i, j);
    int second = std::min(i, j);
    const vector force_vec = compute_aa_force_and_energy(
        actin, first, second, box,
        _current_interaction_ramp() * k_aa,
        _current_interaction_ramp() * kappa_aa,
        aa_cutoff, aa_optimal, aa_antiparallel_overlap);
    if (force_vec.size() < 3) {
        return {};
    }
    const std::array<double, 3> force_on_first{
        force_vec[0], force_vec[1], force_vec[2]};
    const std::array<double, 3> pointed_first{
        actin.direction[first].x, actin.direction[first].y, actin.direction[first].z};
    const std::array<double, 3> pointed_second{
        actin.direction[second].x, actin.direction[second].y, actin.direction[second].z};
    const auto canonical_loads = actn2_load_kinetics::signed_arm_loads(
        force_on_first, pointed_first, pointed_second);
    if (i == first) {
        return canonical_loads;
    }
    return {canonical_loads.second, canonical_loads.first};
}

double Sarcomere::_actin_motor_load(int actin_index) const {
    if (!dynamic_am) {
        return actin_f_load_cb
            ? std::max(0.0, (*actin_f_load_cb)[actin_index])
            : std::max(0.0, actin.f_load[actin_index]);
    }
    double load = 0.0;
    for (int myosin_index : myosinIndicesPerActin.getConnections(actin_index)) {
        if (am_bonds[actin_index][myosin_index] != 1 ||
            am_state[actin_index][myosin_index] !=
                static_cast<BondFlag>(am_motor_kinetics::State::strong)) {
            continue;
        }
        load += _am_resisting_force(actin_index, myosin_index);
    }
    return load;
}

std::array<double, 3> Sarcomere::_actin_motor_force_vector(int actin_index) const {
    std::array<double, 3> result{0.0, 0.0, 0.0};
    if (!dynamic_am) return result;
    for (int myosin_index : myosinIndicesPerActin.getConnections(actin_index)) {
        const double force = _am_active_force(actin_index, myosin_index);
        result[0] += force * actin.direction[actin_index].x;
        result[1] += force * actin.direction[actin_index].y;
        result[2] += force * actin.direction[actin_index].z;
    }
    return result;
}

std::array<double, 3> Sarcomere::_actin_productive_motor_force_vector(
    int actin_index) const {
    std::array<double, 3> result{0.0, 0.0, 0.0};
    if (!dynamic_am) return result;
    for (int myosin_index :
         myosinIndicesPerActin.getConnections(actin_index)) {
        const double actin_dot_myosin = std::clamp(
            actin.direction[actin_index].dot(
                myosin.direction[myosin_index]),
            -1.0, 1.0);
        for (int half_index = 0; half_index < 2; ++half_index) {
            const auto& heads =
                am_heads_by_myosin[static_cast<size_t>(myosin_index)]
                                   [static_cast<size_t>(half_index)];
            for (int h = 0; h < static_cast<int>(heads.size()); ++h) {
                const auto& head = heads[static_cast<size_t>(h)];
                if (head.bound_actin != actin_index ||
                    head.state != static_cast<BondFlag>(
                        am_motor_kinetics::State::strong)) {
                    continue;
                }
                const double productive_alignment =
                    actn2_load_kinetics::productive_polarity_score(
                        actin_dot_myosin, static_cast<int>(head.half));
                if (!(productive_alignment > 0.0)) continue;
                const double force = productive_alignment *
                    _am_head_active_force(
                        actin_index, myosin_index,
                        _am_head_id(half_index, h));
                result[0] += force * actin.direction[actin_index].x;
                result[1] += force * actin.direction[actin_index].y;
                result[2] += force * actin.direction[actin_index].z;
            }
        }
    }
    return result;
}

std::vector<Sarcomere::ProductiveMotorContribution>
Sarcomere::_actin_productive_motor_contributions(int actin_index) const {
    std::vector<ProductiveMotorContribution> result;
    if (!dynamic_am) return result;
    for (int myosin_index :
         myosinIndicesPerActin.getConnections(actin_index)) {
        const double actin_dot_myosin = std::clamp(
            actin.direction[actin_index].dot(
                myosin.direction[myosin_index]),
            -1.0, 1.0);
        for (int half_index = 0; half_index < 2; ++half_index) {
            ProductiveMotorContribution contribution;
            contribution.myosin = myosin_index;
            const auto& heads =
                am_heads_by_myosin[static_cast<size_t>(myosin_index)]
                                   [static_cast<size_t>(half_index)];
            for (int h = 0; h < static_cast<int>(heads.size()); ++h) {
                const auto& head = heads[static_cast<size_t>(h)];
                if (head.bound_actin != actin_index ||
                    head.state != static_cast<BondFlag>(
                        am_motor_kinetics::State::strong)) {
                    continue;
                }
                const double productive_alignment =
                    actn2_load_kinetics::productive_polarity_score(
                        actin_dot_myosin, static_cast<int>(head.half));
                if (!(productive_alignment > 0.0)) continue;
                contribution.half = head.half;
                const double force = productive_alignment *
                    _am_head_active_force(
                        actin_index, myosin_index,
                        _am_head_id(half_index, h));
                contribution.force[0] +=
                    force * actin.direction[actin_index].x;
                contribution.force[1] +=
                    force * actin.direction[actin_index].y;
                contribution.force[2] +=
                    force * actin.direction[actin_index].z;
            }
            if (actn2_load_kinetics::norm(contribution.force) > 0.0) {
                result.push_back(contribution);
            }
        }
    }
    return result;
}

double Sarcomere::_aa_productive_opposed_load(int i, int j) const {
    const std::array<double, 3> direction_i{
        actin.direction[i].x, actin.direction[i].y, actin.direction[i].z};
    const std::array<double, 3> direction_j{
        actin.direction[j].x, actin.direction[j].y, actin.direction[j].z};
    return actn2_load_kinetics::productive_opposed_load(
        direction_i, direction_j,
        _actin_productive_motor_force_vector(i),
        _actin_productive_motor_force_vector(j));
}

double Sarcomere::_aa_transmitted_load(int i, int j) {
    int first = i;
    int second = j;
    const vector force_vec = compute_aa_force_and_energy(
        actin, first, second, box,
        _current_interaction_ramp() * k_aa,
        _current_interaction_ramp() * kappa_aa,
        aa_cutoff, aa_optimal, aa_antiparallel_overlap);
    if (force_vec.size() < 3) return 0.0;
    return std::sqrt(force_vec[0] * force_vec[0] +
                     force_vec[1] * force_vec[1] +
                     force_vec[2] * force_vec[2]);
}

actn2_load_kinetics::ArmLoads Sarcomere::_aa_effective_arm_loads(int i, int j) {
    if (aa_status_mode ==
            actn2_load_kinetics::StatusMode::productive_opposed ||
        aa_status_mode ==
            actn2_load_kinetics::StatusMode::productive_register) {
        const double load = _aa_productive_opposed_load(i, j);
        return {load, load};
    }
    if (aa_status_mode == actn2_load_kinetics::StatusMode::motor_only ||
        aa_status_mode == actn2_load_kinetics::StatusMode::transmitted_only ||
        aa_status_mode == actn2_load_kinetics::StatusMode::opposed_and_transmitted) {
        const double load = actn2_load_kinetics::catch_load(
            aa_status_mode,
            _actin_motor_force_vector(i),
            _actin_motor_force_vector(j),
            _aa_transmitted_load(i, j));
        return {load, load};
    }
    const auto raw = _aa_signed_arm_loads(i, j);
    if (aa_status_mode ==
        actn2_load_kinetics::StatusMode::arm_load_productive) {
        return raw;
    }
    // The current AA potential contains passive geometric prestress. Limit the
    // directionally compatible connector reaction by the favorable motor load
    // available on each actin so passive AA relaxation alone cannot create a
    // force-stabilized state. This is an actual-strain load for dynamic AM and
    // the documented overlap proxy for legacy AM.
    return actn2_load_kinetics::limit_by_available_load(
        raw, _actin_motor_load(i), _actin_motor_load(j));
}

int Sarcomere::determine_cb_status(int& i, int& j){
    auto clear_status2_memory = [&]() {
        aa_status2_memory_until[i][j] = 0;
        aa_status2_memory_until[j][i] = 0;
    };
    auto refresh_status2_memory = [&]() {
        if (!dynamic_am || am_status2_memory_steps == 0) return;
        const size_t max_checkpoint_step =
            static_cast<size_t>(std::numeric_limits<int>::max());
        if (current_step > max_checkpoint_step ||
            am_status2_memory_steps > max_checkpoint_step - current_step) {
            throw std::overflow_error("AM status-2 memory step exceeds checkpoint range");
        }
        const RecoveryStep until =
            static_cast<RecoveryStep>(current_step + am_status2_memory_steps);
        aa_status2_memory_until[i][j] = until;
        aa_status2_memory_until[j][i] = until;
    };
    auto retains_status2 = [&]() {
        return dynamic_am && am_status2_memory_steps > 0 &&
            actin_actin_bonds_prev[i][j] == 1 &&
            actin_actin_status_prev[i][j] == 2 &&
            static_cast<size_t>(aa_status2_memory_until[i][j]) >= current_step;
    };
    double crosslink_i = std::clamp(actin["crosslink_ratio"][i], 0.0, 1.0);
    double crosslink_j = std::clamp(actin["crosslink_ratio"][j], 0.0, 1.0);
    // if (crosslink_i <= EPS || crosslink_j <= EPS) {
    //     return 0;
    // }

    double distance;
    if (directional) {
        vec crosslink_point_i = actin.left_end[i];
        crosslink_point_i += actin.direction[i] * (actin.length * crosslink_i);
        vec crosslink_point_j = actin.left_end[j];
        crosslink_point_j += actin.direction[j] * (actin.length * crosslink_j);
        // Compute geometric metrics using the first binding-zone points as endpoints
        distance = geometry::segment_segment_distance(
            actin.left_end[i], crosslink_point_i, actin.left_end[j], crosslink_point_j, box, is_periodic);
    } else {
        // Use the stored crosslinkable segment endpoints (from the chosen myosin)
        vec start_i = actin_crosslink_start[i];
        vec end_i = actin_crosslink_end[i];
        vec start_j = actin_crosslink_start[j];
        vec end_j = actin_crosslink_end[j];
        distance = geometry::segment_segment_distance(start_i, end_i, start_j, end_j, box, is_periodic);
    }
    double cos_angle = actin.direction[i].dot(actin.direction[j]);

    bool was_strong = (actin_actin_status_prev[i][j] == 2);
    int thread_id = omp_get_thread_num();
    auto& local_breakage_events = cb_breakage_events_temp[thread_id];
    auto& local_completed_lifetimes = aa_completed_lifetimes_temp[thread_id];

    auto record_break = [&](void){
        auto& myosin_indices_i = myosinIndicesPerActin.getConnections(i);
        auto& myosin_indices_j = myosinIndicesPerActin.getConnections(j);
        double tension_i = actin_basic_tension[i];
        double tension_j = actin_basic_tension[j];
        local_breakage_events.insert(local_breakage_events.end(),
                                     {static_cast<double>(i), static_cast<double>(j),
                                      static_cast<double>(current_step), distance, cos_angle,
                                      tension_i, tension_j,
                                      actin_crosslink_ratio[i], actin_crosslink_ratio[j],
                                      static_cast<double>(myosin_indices_i.size()),
                                      static_cast<double>(myosin_indices_j.size())});
        for (int k = 0; k < max_myosin_bonds; ++k) {
            local_breakage_events.push_back(
                k < static_cast<int>(myosin_indices_i.size()) ?
                    static_cast<double>(myosin_indices_i[k]) : -1.0);
        }
        for (int k = 0; k < max_myosin_bonds; ++k) {
            local_breakage_events.push_back(
                k < static_cast<int>(myosin_indices_j.size()) ?
                    static_cast<double>(myosin_indices_j[k]) : -1.0);
        }
    };
    bool crosslink = false;
    // if ((crosslink_i > EPS) && (crosslink_j > EPS) || !directional) {
    if (aa_refractory::mechanically_active(
            actin_actin_bonds_prev[i][j] == 1, distance < aa_cutoff)) {
        crosslink = true;
        // Per-pair, per-frame debug print disabled: at high actin density this emits tens of
        // millions of lines (tens of GB of stdout per run), starving HDF5 I/O. Re-enable only
        // for targeted debugging.
        // printf("Actins %d and %d are crosslinked, distance: %f, cos_angle: %f, crosslink ratio: %f, %f\n",
        //        i, j, distance, cos_angle,crosslink_i, crosslink_j);
    }
    //}
    if (!crosslink){
        clear_status2_memory();
        if (was_strong){
            int myosin_i = myosinIndicesPerActin.getConnections(i)[0];
            int myosin_j = myosinIndicesPerActin.getConnections(j)[0];
            // Per-pair, per-frame debug print disabled (see matching note above): tens of GB of
            // stdout per run at high density. Re-enable only for targeted debugging.
            // printf("Actins %d and %d no longer crosslinked, distance: %f, cos_angle: %f, crosslink ratio: %f, %f myosin indices %d, %d\n",
            //      i, j, distance, cos_angle,crosslink_i, crosslink_j, myosin_i, myosin_j);
            (void)myosin_i; (void)myosin_j;
            record_break();
            // Record completed lifetime on geometric loss once
            if (aa_attach_step.size() && aa_attach_step[i][j] >= 0 &&
                actin_actin_status_prev[i][j] == 2) {
                local_completed_lifetimes.push_back((current_step - aa_attach_step[i][j]) * dt);
            }
            if (aa_attach_step.size()) {
                aa_attach_step[i][j] = -1;
                aa_attach_step[j][i] = -1;
            }
        }
        return 0; // return -1 for non-catch bond
    }

    if (aa_status_mode != actn2_load_kinetics::StatusMode::legacy_motor_identity) {
        // A newly associated ACTN2 crosslink begins in the ordinary bonded
        // state. It can become load-stabilized only after an existing bond has
        // actually transmitted favorable force on both arms.
        if (actin_actin_bonds_prev[i][j] != 1) {
            clear_status2_memory();
            return 1;
        }
        const auto loads = _aa_effective_arm_loads(i, j);
        const int status =
            aa_status_mode == actn2_load_kinetics::StatusMode::arm_load_productive
            ? actn2_load_kinetics::classify_arm_loaded_productive(
                  static_cast<int>(actin_actin_status_prev[i][j]), loads,
                  _aa_productive_opposed_load(i, j),
                  aa_status2_force_on, aa_status2_force_off)
            : aa_status_mode == actn2_load_kinetics::StatusMode::force_projection
            ? actn2_load_kinetics::classify_existing_bond(
                  static_cast<int>(actin_actin_status_prev[i][j]), loads,
                  aa_status2_force_on, aa_status2_force_off)
            : actn2_load_kinetics::classify_existing_bond_scalar(
                  static_cast<int>(actin_actin_status_prev[i][j]),
                  loads.favorable_bilateral(),
                  aa_status2_force_on, aa_status2_force_off);
        if (status == 2) {
            refresh_status2_memory();
            return 2;
        }
        // This is a persistence-only memory: it cannot establish status 2,
        // but it can bridge a short loss of mechanical support after an
        // existing status-2 crosslink was loaded.
        if (retains_status2()) {
            return 2;
        }
        clear_status2_memory();
        if (was_strong) {
            record_break();
        }
        return status;
    }

    bool catch_bond =(actin_basic_tension[i]>EPS && actin_basic_tension[j]>EPS);
    if (directional){
        catch_bond = (catch_bond && cos_angle<0);
    }
    if (!catch_bond){
        if (retains_status2()) {
            return 2;
        }
        if (was_strong){
            record_break();
        }
        return 1;
    }
    auto& myosin_indices_i = myosinIndicesPerActin.getConnections(i);
    auto& myosin_indices_j = myosinIndicesPerActin.getConnections(j);
    if (myosin_indices_i.empty() || myosin_indices_j.empty()){
        if (retains_status2()) {
            return 2;
        }
        if (was_strong){
            record_break();
        }
        return 1;
    }
    for (int mi : myosin_indices_i){
        if (am_bonds[i][mi] != 1) continue;
        if (dynamic_am && am_state[i][mi] !=
            static_cast<BondFlag>(am_motor_kinetics::State::strong)) continue;
        for (int mj : myosin_indices_j){
            if (am_bonds[j][mj] != 1) continue;
            if (dynamic_am && am_state[j][mj] !=
                static_cast<BondFlag>(am_motor_kinetics::State::strong)) continue;
            if (mi != mj){
                refresh_status2_memory();
                return 2;
            }
        }
    }
    if (retains_status2()) {
        return 2;
    }
    if (was_strong){
        record_break();
    }
    return 1;
}

void Sarcomere::_record_chem_transition(int& i, int& j, double affinity, bool binding) {
    int thread_id = omp_get_thread_num();
    double split_entropy = 0.5 * affinity;
    double split_input = 0.5 * std::abs(affinity);
    actin_chem_entropy_delta_temp[thread_id][i] += split_entropy;
    actin_chem_entropy_delta_temp[thread_id][j] += split_entropy;
    actin_chem_input_delta_temp[thread_id][i] += split_input;
    actin_chem_input_delta_temp[thread_id][j] += split_input;
    if (binding) {
        actin_chem_binding_delta_temp[thread_id][i] += 0.5;
        actin_chem_binding_delta_temp[thread_id][j] += 0.5;
    } else {
        actin_chem_unbinding_delta_temp[thread_id][i] += 0.5;
        actin_chem_unbinding_delta_temp[thread_id][j] += 0.5;
    }
}

bool Sarcomere::_cb_decide(int& i, int& j, int status){
    if (status == 0){
        return false;
    }
    const bool was_bonded_prev = actin_actin_bonds_prev[i][j] == 1;
    if (!aa_refractory::can_attempt_association(
            was_bonded_prev,
            current_step,
            static_cast<size_t>(actin_recovery_until[i][j]))) {
        return false;
    }
    int thread_id = omp_get_thread_num();
    auto& local_actin_cb_status = actin_cb_status_temp[thread_id];
    auto& local_completed_lifetimes = aa_completed_lifetimes_temp[thread_id];

    counter_rng::CounterRng pair_rng(
        static_cast<std::uint64_t>(initial_seed),
        static_cast<std::uint64_t>(current_step),
        counter_rng::Stream::legacy_aa_pair,
        counter_rng::canonical_pair_key(i, j));
    double rand = pair_rng.uniform();
    auto current_elastic_energy = [&]() {
        const double distance = geometry::segment_segment_distance(
            actin.left_end[i], actin.right_end[i],
            actin.left_end[j], actin.right_end[j], box, is_periodic);
        return aa_association_kinetics::elastic_cost(
            distance, aa_optimal,
            actin.direction[i].dot(actin.direction[j]),
            aa_encounter_k, aa_encounter_kappa);
    };
    double f_load = 0.0;
    if (aa_status_mode != actn2_load_kinetics::StatusMode::legacy_motor_identity) {
        f_load = _aa_effective_arm_loads(i, j).favorable_bilateral();
        if (aa_status_mode ==
                actn2_load_kinetics::StatusMode::arm_load_productive &&
            !(_aa_productive_opposed_load(i, j) > 0.0)) {
            f_load = 0.0;
        }
    } else {
        const double abs_cos_angle = std::abs(actin.direction[i].dot(actin.direction[j]));
        const double f_load_i = actin_f_load_cb ? (*actin_f_load_cb)[i] : actin.f_load[i];
        const double f_load_j = actin_f_load_cb ? (*actin_f_load_cb)[j] : actin.f_load[j];
        f_load = abs_cos_angle * std::min(f_load_i, f_load_j);
    }
    if (actin_actin_bonds_prev[i][j] == 1) {
        const double lifetime = catch_bond_kinetics::lifetime_for_status(
            status,
            base_lifetime,
            directional_base_lifetime,
            lifetime_coeff,
            f_load,
            stabilize_status1);
        double k_off_adjusted =
            catch_bond_kinetics::break_probability(dt, lifetime);
        if (aa_passive_detailed_balance) {
            const double reference_off_rate = lifetime > 0.0
                ? 1.0 / lifetime
                : std::numeric_limits<double>::infinity();
            const double off_rate = aa_association_kinetics::transition_state_off_rate(
                reference_off_rate, thermal_beta, current_elastic_energy(),
                aa_transition_state_theta);
            k_off_adjusted = aa_association_kinetics::event_probability(off_rate, dt);
        }
        if (rand < k_off_adjusted){
            aa_status2_memory_until[i][j] = 0;
            aa_status2_memory_until[j][i] = 0;
            if (actin_actin_status_prev[i][j] == 2) {
                double affinity = std::log(std::max(k_on, EPS) * std::max(lifetime, EPS));
                _record_chem_transition(i, j, -affinity, false);
            }
            if (f_load>0) {
                }
            // If this pair had an attach timestamp, record completed lifetime
            if (aa_attach_step.size() && aa_attach_step[i][j] >= 0 &&
                actin_actin_status_prev[i][j] == 2) {
                local_completed_lifetimes.push_back((current_step - aa_attach_step[i][j]) * dt);
            }
            if (aa_attach_step.size()) {
                aa_attach_step[i][j] = -1;
                aa_attach_step[j][i] = -1;
            }
            if (bond_recovery_steps > 0) {
                const size_t max_checkpoint_step =
                    static_cast<size_t>(std::numeric_limits<int>::max());
                if (current_step > max_checkpoint_step ||
                    bond_recovery_steps > max_checkpoint_step - current_step) {
                    throw std::overflow_error("actin recovery step exceeds checkpoint integer range");
                }
                const RecoveryStep release = static_cast<RecoveryStep>(current_step + bond_recovery_steps);
                actin_recovery_until[i][j] = release;
                actin_recovery_until[j][i] = release;
                recovery_pairs_temp[thread_id].emplace_back(i, j);
            }
            return false;
        }
    }

    else{
        double association_rate = k_on;
        if (aa_energy_gated_association || aa_passive_detailed_balance) {
            association_rate = aa_association_kinetics::transition_state_on_rate(
                k_on, thermal_beta, current_elastic_energy(),
                aa_transition_state_theta);
        }
        if (rand >= aa_association_kinetics::event_probability(
                        association_rate, dt)){
            //if (status==2){
            // printf("Actins %d and %d fail to form catch bond with status %d, rand %f, k_on*dt %f, f_load %f, abs_cos_angle %f\n",
            //     i, j, status, rand, k_on*dt, f_load, abs_cos_angle);
            //}
            return false;
        }
        if (max_aa_bonds_per_actin > 0) {
            bool reserved = false;
            #pragma omp critical(aa_valency_reservation)
            {
                if (aa_reserved_degree[i] < max_aa_bonds_per_actin &&
                    aa_reserved_degree[j] < max_aa_bonds_per_actin) {
                    ++aa_reserved_degree[i];
                    ++aa_reserved_degree[j];
                    reserved = true;
                }
            }
            if (!reserved) return false;
        }
        if (status == 2) {
            const double lifetime = catch_bond_kinetics::lifetime_for_status(
                status,
                base_lifetime,
                directional_base_lifetime,
                lifetime_coeff,
                f_load,
                stabilize_status1);
            double affinity = std::log(
                std::max(association_rate, EPS) * std::max(lifetime, EPS));
            _record_chem_transition(i, j, affinity, true);
        }
    //     if (status==2){
    //     printf("Actins %d and %d form/maintain catch bond with status %d, rand %f, f_load %f, abs_cos_angle %f\n",
    //         i, j, status, rand, f_load, abs_cos_angle);}
     }
    // Bond is accepted or persists: record status locally
    local_actin_cb_status[i] = std::max(local_actin_cb_status[i], status);
    local_actin_cb_status[j] = std::max(local_actin_cb_status[j], status);

    return true;
}

void Sarcomere::_set_cb(int& i, int& j, int status){
    int thread_id = omp_get_thread_num();
    auto& local_actin_forces = actin_forces_temp[thread_id];
    auto& local_actin_torques = actin_torques_temp[thread_id];
    auto& local_actin_cb_status = actin_cb_status_temp[thread_id];
    auto& local_myosin_forces = myosin_forces_temp[thread_id];
    auto& local_myosin_torques = myosin_torques_temp[thread_id];
 
    const double ramp = _current_interaction_ramp();
    vector force_vec = compute_aa_force_and_energy(
        actin, i, j, box, ramp * k_aa, ramp * kappa_aa,
        aa_cutoff, aa_optimal, aa_antiparallel_overlap);
    const std::array<double, 3> force_on_i{
        force_vec[0], force_vec[1], force_vec[2]};
    const std::array<double, 3> pointed_i{
        actin.direction[i].x, actin.direction[i].y, actin.direction[i].z};
    const std::array<double, 3> pointed_j{
        actin.direction[j].x, actin.direction[j].y, actin.direction[j].z};
    const auto arm_loads = actn2_load_kinetics::signed_arm_loads(
        force_on_i, pointed_i, pointed_j);
    const auto effective_arm_loads =
        actn2_load_kinetics::limit_by_available_load(
            arm_loads, _actin_motor_load(i), _actin_motor_load(j));
    const double bilateral_load = effective_arm_loads.favorable_bilateral();
    // Multiple OpenMP workers can own different pairs that share an actin.
    #pragma omp critical(aa_load_feature_update)
    {
        actin["aa_raw_arm_load_max"][i] = std::max(
            actin["aa_raw_arm_load_max"][i], std::max(0.0, arm_loads.first));
        actin["aa_raw_arm_load_max"][j] = std::max(
            actin["aa_raw_arm_load_max"][j], std::max(0.0, arm_loads.second));
        actin["aa_arm_load_max"][i] = std::max(
            actin["aa_arm_load_max"][i], effective_arm_loads.first);
        actin["aa_arm_load_max"][j] = std::max(
            actin["aa_arm_load_max"][j], effective_arm_loads.second);
        actin["aa_bilateral_load_max"][i] = std::max(
            actin["aa_bilateral_load_max"][i], bilateral_load);
        actin["aa_bilateral_load_max"][j] = std::max(
            actin["aa_bilateral_load_max"][j], bilateral_load);
    }
    local_actin_forces[i].x += force_vec[0];
    local_actin_forces[i].y += force_vec[1];
    local_actin_forces[i].z += force_vec[2];
    local_actin_forces[j].x -= force_vec[0];
    local_actin_forces[j].y -= force_vec[1];
    local_actin_forces[j].z -= force_vec[2];
    local_actin_torques[i].x += force_vec[3];
    local_actin_torques[i].y += force_vec[4];
    local_actin_torques[i].z += force_vec[5];
    local_actin_torques[j].x += force_vec[6];
    local_actin_torques[j].y += force_vec[7];
    local_actin_torques[j].z += force_vec[8];
    actin_n_bonds[i] += 1;
    actin_n_bonds[j] += 1;
    if (status == 2){
        actin_strong_cb_count[i] += 1;
        actin_strong_cb_count[j] += 1;
    }
    actin_actin_bonds[i][j] = 1;
    actin_actin_bonds[j][i] = 1;
    actin_actin_status[i][j] = status;
    actin_actin_status[j][i] = status;
    // Update lifetime: increment if bond persisted, reset if new
    actin_actin_lifetime[i][j] = actin_actin_lifetime_prev[i][j] + 1;
    actin_actin_lifetime[j][i] = actin_actin_lifetime[i][j];
    active_aa_pairs_temp[thread_id].emplace_back(i, j);
    // If the bond is newly formed this step, stamp attach time
    if (actin_actin_bonds_prev[i][j] == 0) {
        if (aa_attach_step.size()) {
            aa_attach_step[i][j] = static_cast<int>(current_step);
            aa_attach_step[j][i] = static_cast<int>(current_step);
        }
        if (bond_recovery_steps > 0) {
            if (current_step > static_cast<size_t>(std::numeric_limits<int>::max())) {
                throw std::overflow_error("actin recovery step exceeds checkpoint integer range");
            }
            const RecoveryStep recovery_step = static_cast<RecoveryStep>(current_step);
            actin_recovery_until[i][j] = recovery_step;
            actin_recovery_until[j][i] = recovery_step;
        }
    }
}

void Sarcomere::_set_cb(int& i, std::vector<int> indices, std::vector<int> status){
    int thread_id = omp_get_thread_num();
    bool add_connection;
    for (size_t index = 0; index < indices.size(); index++){
        int j = indices[index];
        _set_cb(i,j,status[index]);
    }
}



std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
    Sarcomere::_extract_bonded_pairs(
    const std::vector<std::vector<BondFlag>>& actin_actin_bonds,
    const std::vector<std::vector<BondFlag>>& actin_actin_status,
    const utils::MoleculeConnection& myosinIndicesPerActin)
{
    std::vector<double> flatActinBonds;
    flatActinBonds.reserve(active_aa_pairs.size() * 2);
    for (const auto& [i, j] : active_aa_pairs) {
        if (actin_actin_bonds[i][j] == 1 &&
            actin_actin_status[i][j] == 2 &&
            actin.cb_status[i] == 2 &&
            actin.cb_status[j] == 2) {
            flatActinBonds.push_back(static_cast<double>(i));
            flatActinBonds.push_back(static_cast<double>(j));
        }
    }

    // Next, process the flattened actin bonds to extract unique myosin–myosin pairs.
    // For each actin–actin bond (i, j), pair every myosin attached to actin i with every
    // myosin attached to actin j.
    std::set<std::pair<int,int>> uniqueMyosinBonds;
    for (size_t k = 0; k < flatActinBonds.size(); k += 2) {
        int i = static_cast<int>(flatActinBonds[k]);
        int j = static_cast<int>(flatActinBonds[k+1]);

        auto connections_i = myosinIndicesPerActin.getConnections(i);
        auto connections_j = myosinIndicesPerActin.getConnections(j);

        // Skip if either actin has no attached myosin.
        if (connections_i.empty() || connections_j.empty())
            continue;

        for (int myosin_i_raw : connections_i) {
            for (int myosin_j_raw : connections_j) {
                if (myosin_i_raw == myosin_j_raw) {
                    continue;
                }
                int myosin_a = myosin_i_raw;
                int myosin_b = myosin_j_raw;
                if (myosin_a > myosin_b) {
                    std::swap(myosin_a, myosin_b);
                }
                uniqueMyosinBonds.insert(std::make_pair(myosin_a, myosin_b));
            }
        }
    }

    // Flatten the set of unique myosin bonds into a vector<double>.
    std::vector<double> flattenedMyosinBonds;
    for (const auto& bond : uniqueMyosinBonds) {
        flattenedMyosinBonds.push_back(static_cast<double>(bond.first));
        flattenedMyosinBonds.push_back(static_cast<double>(bond.second));
    }

    std::vector<double> flatActinMyosinBonds;
    flatActinMyosinBonds.reserve(active_am_pairs.size() * 2);
    for (const auto& [a, m] : active_am_pairs) {
        if (am_bonds[a][m] == 1) {
            flatActinMyosinBonds.push_back(static_cast<double>(a));
            flatActinMyosinBonds.push_back(static_cast<double>(m));
        }
    }

    // Return the triplet: actin bonds, myosin bonds, and actin–myosin bonds.
    return {flatActinBonds, flattenedMyosinBonds, flatActinMyosinBonds};
}


vec Sarcomere::_alignment_torque(const vec& u, double k_bias)
{
    // u must be unit length
    return { k_bias * (1.0 - u.x),
             -k_bias * u.y,
             -k_bias * u.z};
}

void Sarcomere::_apply_cb_alignment_bias(double& k_theta_bias)
{
    #pragma omp for
        for (int i = 0; i < myosin.n; ++i) {
            vec u = myosin.direction[i];
            u.normalize();
            auto idxs = actinIndicesPerMyosin.getConnections(i);

            double acc_cb = 0.1; // accumulated cross-bridge status
            for (int a : idxs) {
                double cb = static_cast<double>(actin.cb_status[a]);
                if (cb > 0.1) {           // active cross-bridge
                    acc_cb += cb;
                    if (acc_cb > 1.0) { acc_cb = 1.0; break; }
                }
            }
            vec tau = _alignment_torque(u, k_theta_bias * acc_cb);
            myosin_torques_temp[omp_get_thread_num()][i] += tau;
        }
    #pragma omp for
        for (int i = 0; i < actin.n; ++i) {
            if (actin.cb_status[i] < 2) continue; // skip if no active cross-bridge
            vec u = actin.direction[i];
            u.normalize();
            vec tau = _alignment_torque(u, k_theta_bias * actin.cb_status[i]);
            actin_torques_temp[omp_get_thread_num()][i] += tau;
        }
}

double Sarcomere::_current_bundle_strength() const {
    if (k_bundle_max <= 0.0) {
        return 0.0;
    }
    if (bundle_ramp_steps <= 0) {
        return k_bundle_max;
    }
    size_t step_for_ramp = (current_step > 0) ? (current_step - 1) : 0;
    double ramp = static_cast<double>(step_for_ramp) / static_cast<double>(bundle_ramp_steps);
    ramp = std::clamp(ramp, 0.0, 1.0);
    return k_bundle_max * ramp * ramp;
}

void Sarcomere::set_interaction_ramp(size_t start_step, size_t ramp_steps) {
    interaction_ramp_start_step = start_step;
    interaction_ramp_steps = ramp_steps;
}

void Sarcomere::set_steric_stiffness(double actin_stiffness,
                                     double myosin_stiffness,
                                     double myosin_exclusion) {
    actin_steric_k = std::max(0.0, actin_stiffness);
    myosin_steric_k = std::max(0.0, myosin_stiffness);
    myosin_exclusion_distance = std::max(0.0, myosin_exclusion);
}

void Sarcomere::configure_double_overlap_contact(
    bool enabled,
    double contact_diameter,
    double epsilon,
    double force_cap) {
    aa_double_overlap_contact_enabled = enabled;
    aa_double_overlap_contact_diameter =
        std::max(0.0, contact_diameter);
    aa_double_overlap_contact_epsilon = std::max(0.0, epsilon);
    aa_double_overlap_contact_force_cap = std::max(0.0, force_cap);
}

double Sarcomere::_current_interaction_ramp() const {
    if (interaction_ramp_steps == 0) return 1.0;
    if (current_step <= interaction_ramp_start_step) return 0.0;
    const double x = std::clamp(
        static_cast<double>(current_step - interaction_ramp_start_step) /
            static_cast<double>(interaction_ramp_steps),
        0.0, 1.0);
    // Smoothstep avoids a force derivative discontinuity at both ends.
    return x * x * (3.0 - 2.0 * x);
}

void Sarcomere::_apply_probe_confinement() {
    if (!(probe_confine_k > 0.0) || actin.n <= 0) return;
    // Lazy capture: the first force evaluation precedes any motion, so this is
    // the initial configuration.  NOTE this makes the trap anchor to the
    // RESUMED position if a run is ever resumed; the probe does not resume.
    if (probe_confine_anchor.size() != static_cast<size_t>(actin.n)) {
        probe_confine_anchor.assign(static_cast<size_t>(actin.n),
                                    vec{0.0, 0.0, 0.0});
        for (int i = 0; i < actin.n; ++i) {
            probe_confine_anchor[static_cast<size_t>(i)] =
                static_cast<vec>(actin.center[i]);
        }
    }
    for (int i = 0; i < actin.n; ++i) {
        vec offset = static_cast<vec>(actin.center[i]) -
                     probe_confine_anchor[static_cast<size_t>(i)];
        offset.pbc_wrap(box, actin.periodic_axes);
        // Force at the centre: no lever arm, hence no torque contribution.
        actin.force[i].x -= probe_confine_k * offset.x;
        actin.force[i].y -= probe_confine_k * offset.y;
        actin.force[i].z -= probe_confine_k * offset.z;
    }
}

void Sarcomere::_apply_transverse_bundling(double k_bundle) {
    if (k_bundle <= 0.0) {
        return;
    }
    double sum_y = 0.0;
    double sum_z = 0.0;
    size_t participant_count = 0;

    for (int i = fix_myosin; i < myosin.n; ++i) {
        sum_y += myosin.center[i].y;
        sum_z += myosin.center[i].z;
        ++participant_count;
    }
    for (int i = 0; i < actin.n; ++i) {
        sum_y += actin.center[i].y;
        sum_z += actin.center[i].z;
        ++participant_count;
    }

    if (participant_count == 0) {
        return;
    }

    double com_y = sum_y / static_cast<double>(participant_count);
    double com_z = sum_z / static_cast<double>(participant_count);

    for (int i = fix_myosin; i < myosin.n; ++i) {
        double dy = myosin.center[i].y - com_y;
        double dz = myosin.center[i].z - com_z;
        myosin.force[i].y -= k_bundle * dy;
        myosin.force[i].z -= k_bundle * dz;
    }
    for (int i = 0; i < actin.n; ++i) {
        double dy = actin.center[i].y - com_y;
        double dz = actin.center[i].z - com_z;
        actin.force[i].y -= k_bundle * dy;
        actin.force[i].z -= k_bundle * dz;
    }
}

void Sarcomere::_update_myosin_bond_matrix() {
    if (myosin.n <= 0) {
        myosin_bond_matrix.clear();
        has_myosin_bond_pairs = false;
        return;
    }
    if (myosin_bond_matrix.size() != static_cast<size_t>(myosin.n)) {
        myosin_bond_matrix.assign(myosin.n, std::vector<BondFlag>(myosin.n, 0));
    } else {
        for (auto& row : myosin_bond_matrix) {
            std::fill(row.begin(), row.end(), 0);
        }
    }
    has_myosin_bond_pairs = false;
    std::vector<std::vector<int>> actin_attached(actin.n);
    for (int act_idx = 0; act_idx < actin.n; ++act_idx) {
        const auto& attachments = myosinIndicesPerActin.getConnections(act_idx);
        actin_attached[act_idx].reserve(attachments.size());
        for (int m_idx : attachments) {
            if (m_idx < 0 || m_idx >= myosin.n) {
                continue;
            }
            if (am_interaction[act_idx][m_idx].partial_binding_ratio > EPS) {
                actin_attached[act_idx].push_back(m_idx);
            }
        }
    }
    for (int a = 0; a < actin.n; ++a) {
        if (actin.cb_status[a] != 2) {
            continue;
        }
        const auto& attached_a = actin_attached[a];
        if (attached_a.empty()) {
            continue;
        }
        for (int b = a + 1; b < actin.n; ++b) {
            if (actin.cb_status[b] != 2) {
                continue;
            }
            if (actin_actin_status[a][b] != 2) {
                continue;
            }
            const auto& attached_b = actin_attached[b];
            if (attached_b.empty()) {
                continue;
            }
            for (int mi : attached_a) {
                for (int mj : attached_b) {
                    if (mi == mj) {
                        continue;
                    }
                    if (myosin_bond_matrix[mi][mj] == 0) {
                        has_myosin_bond_pairs = true;
                    }
                    myosin_bond_matrix[mi][mj] = 1;
                    myosin_bond_matrix[mj][mi] = 1;
                }
            }
        }
    }
}

bool Sarcomere::_myosin_pair_bonded(int mi, int mj) const {
    if (myosin_bond_matrix.empty()) {
        return false;
    }
    if (mi < 0 || mj < 0) {
        return false;
    }
    if (mi >= static_cast<int>(myosin_bond_matrix.size())) {
        return false;
    }
    if (mj >= static_cast<int>(myosin_bond_matrix[mi].size())) {
        return false;
    }
    return myosin_bond_matrix[mi][mj] != 0;
}


// Sarcomere::new_file() now lives in src/sarcomere_io.cpp (moved 2026-09-17).



// save_state, save_resume_snapshot and load_state now live in
// src/sarcomere_io.cpp (moved 2026-09-17), together with new_file and the
// HDF5 packing helpers they use.


// The ACTN2 subsystem (every _local_aa_*/_aa_*/actn2 function, the
// AAPairAxis and AAAssociationProposal helpers, _process_catch_bonds and
// _update_local_aa_links) moved to src/sarcomere_actn2.cpp on 2026-09-17.

