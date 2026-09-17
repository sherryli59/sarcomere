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

namespace {

double boltzmann_exclusion_factor(double beta, double energy) {
    if (!(energy > 0.0) || !std::isfinite(energy)) {
        return std::isfinite(energy) ? 1.0 : 0.0;
    }
    const double exponent = std::clamp(beta * energy, 0.0, 50.0);
    return std::exp(-exponent);
}

void delete_link_if_exists(H5::H5File& file, const std::string& full_path) {
    if (file.nameExists(full_path)) {
        H5Ldelete(file.getId(), full_path.c_str(), H5P_DEFAULT);
    }
}

void move_link(H5::H5File& file, const std::string& src_path, const std::string& dst_path) {
    H5Lmove(file.getId(), src_path.c_str(), file.getId(), dst_path.c_str(), H5P_DEFAULT, H5P_DEFAULT);
}

void promote_temp_group(H5::H5File& file,
                        const std::string& final_path,
                        const std::string& temp_path,
                        const std::string& backup_path) {
    delete_link_if_exists(file, backup_path);
    if (file.nameExists(final_path)) {
        move_link(file, final_path, backup_path);
    }
    move_link(file, temp_path, final_path);
    delete_link_if_exists(file, backup_path);
    file.flush(H5F_SCOPE_GLOBAL);
}

void maybe_test_sleep(const char* env_name, const char* label) {
    const char* raw = std::getenv(env_name);
    if (raw == nullptr || raw[0] == '\0') {
        return;
    }
    const int sleep_seconds = std::max(0, std::atoi(raw));
    if (sleep_seconds <= 0) {
        return;
    }
    fflush(stdout);
    std::this_thread::sleep_for(std::chrono::seconds(sleep_seconds));
}

void truncate_dataset_first_dim(H5::Group& group, const std::string& dataset_name, hsize_t keep_frames) {
    if (!group.nameExists(dataset_name)) {
        return;
    }
    H5::DataSet dataset = group.openDataSet(dataset_name);
    H5::DataSpace filespace = dataset.getSpace();
    const int rank = filespace.getSimpleExtentNdims();
    if (rank <= 0) {
        return;
    }
    std::vector<hsize_t> dims(rank, 0);
    filespace.getSimpleExtentDims(dims.data(), nullptr);
    if (dims.empty() || dims[0] <= keep_frames) {
        return;
    }
    dims[0] = keep_frames;
    dataset.extend(dims.data());
}

struct PackedAABondState {
    std::vector<int> pairs;
    std::vector<int> status;
    std::vector<int> lifetime;
};

struct PackedRecoveryState {
    std::vector<int> pairs;
    std::vector<int> until;
};

template <typename BondMatrix, typename StatusMatrix>
PackedAABondState pack_aa_bond_state(
    const std::vector<std::pair<int, int>>& active_pairs,
    const BondMatrix& bonds,
    const StatusMatrix& status,
    const std::vector<std::vector<int>>& lifetime) {
    PackedAABondState packed;
    packed.pairs.reserve(active_pairs.size() * 2);
    packed.status.reserve(active_pairs.size());
    packed.lifetime.reserve(active_pairs.size());
    for (const auto& [i, j] : active_pairs) {
        if (bonds[i][j] != 1) {
            continue;
        }
        packed.pairs.push_back(i);
        packed.pairs.push_back(j);
        packed.status.push_back(status[i][j]);
        packed.lifetime.push_back(lifetime[i][j]);
    }
    return packed;
}

template <typename BondMatrix>
std::vector<int> pack_am_bond_pairs(
    const std::vector<std::pair<int, int>>& active_pairs,
    const BondMatrix& bonds) {
    std::vector<int> packed;
    packed.reserve(active_pairs.size() * 2);
    for (const auto& [i, j] : active_pairs) {
        if (bonds[i][j] != 1) {
            continue;
        }
        packed.push_back(i);
        packed.push_back(j);
    }
    return packed;
}

template <typename RecoveryMatrix>
PackedRecoveryState pack_recovery_state(
    const std::vector<std::pair<int, int>>& active_pairs,
    const RecoveryMatrix& recovery_until,
    size_t current_step) {
    PackedRecoveryState packed;
    packed.pairs.reserve(active_pairs.size() * 2);
    packed.until.reserve(active_pairs.size());
    for (const auto& [i, j] : active_pairs) {
        const size_t until = static_cast<size_t>(recovery_until[i][j]);
        if (until <= current_step) {
            continue;
        }
        packed.pairs.push_back(i);
        packed.pairs.push_back(j);
        packed.until.push_back(
            static_cast<int>(std::min(until, static_cast<size_t>(std::numeric_limits<int>::max()))));
    }
    return packed;
}

void pad_int_vector(std::vector<int>& values, size_t target_size, int pad_value) {
    // Fit to exactly target_size: pad with pad_value if shorter, TRUNCATE if longer.
    // Truncation lets the fixed-width compact /state bond snapshot degrade gracefully when the
    // bond count exceeds capacity (dense, low-myosin networks) instead of overflowing the dataset.
    // Trajectory data lives in separate datasets and is unaffected.
    values.resize(target_size, pad_value);
}

std::vector<double> flatten_vec_values(const std::vector<vec>& values) {
    std::vector<double> flat;
    flat.reserve(3 * values.size());
    for (const auto& value : values) {
        flat.push_back(value.x);
        flat.push_back(value.y);
        flat.push_back(value.z);
    }
    return flat;
}

}


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

namespace {
// Escape hatch for A/B timing only.  Every gate it controls is exact -- each
// rejects solely cases where the kernel it guards provably contributes zero --
// so switching it off changes runtime and nothing else.  It exists because
// node-to-node speed varies ~10%, which makes a cumulative speedup stitched
// across jobs meaningless; with this, one job measures the whole set.
const bool opt_gates_enabled = [] {
    for (const char* name : {"SARC_OPT_GATES", "SARC_AM_PREFILTER"}) {
        const char* v = std::getenv(name);
        if (v && std::string(v) == "0") return false;
    }
    return true;
}();
}  // namespace

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

std::uint64_t Sarcomere::_aa_pair_key(int i, int j) const {
    const auto first = static_cast<std::uint32_t>(std::min(i, j));
    const auto second = static_cast<std::uint32_t>(std::max(i, j));
    return (static_cast<std::uint64_t>(first) << 32) |
           static_cast<std::uint64_t>(second);
}

void Sarcomere::_rebuild_local_aa_indices() {
    aa_active_link_slots.clear();
    aa_links_by_pair.clear();
    // assign() would destroy and reallocate every per-actin bucket twice per
    // step.  Clearing in place leaves the same empty buckets but keeps their
    // capacity, which is what the incremental add_link() path relies on too.
    if (aa_links_by_actin.size() != static_cast<size_t>(actin.n)) {
        aa_links_by_actin.assign(static_cast<size_t>(actin.n), {});
    } else {
        for (auto& bucket : aa_links_by_actin) bucket.clear();
    }
    for (size_t slot = 0; slot < aa_links.size(); ++slot) {
        const auto& link = aa_links[slot];
        if (!link.active) continue;
        aa_active_link_slots.push_back(slot);
        const bool arm_bridge = !aa_arm_state_model ||
            actn2_arm_kinetics::is_bridge(link.arm_state);
        if (arm_bridge) {
            aa_links_by_pair[_aa_pair_key(link.actin_i, link.actin_j)].push_back(slot);
            aa_links_by_actin[static_cast<size_t>(link.actin_i)].push_back(slot);
            aa_links_by_actin[static_cast<size_t>(link.actin_j)].push_back(slot);
        } else if (link.arm_state ==
                   actn2_arm_kinetics::State::tethered_to_i) {
            aa_links_by_actin[static_cast<size_t>(link.actin_i)].push_back(slot);
        } else if (link.arm_state ==
                   actn2_arm_kinetics::State::tethered_to_j) {
            aa_links_by_actin[static_cast<size_t>(link.actin_j)].push_back(slot);
        }
    }
}

void Sarcomere::set_aa_neighbor_identity(
    bool enabled, double status1_memory, double status2_memory) {
    aa_neighbor_identity_enabled = enabled;
    aa_neighbor_slot_memory_steps = status1_memory > 0.0
        ? static_cast<size_t>(std::ceil(status1_memory / dt)) : 0;
    aa_status2_slot_memory_steps = status2_memory > 0.0
        ? static_cast<size_t>(std::ceil(status2_memory / dt)) : 0;
    if (!enabled) aa_neighbor_slots.clear();
}

bool Sarcomere::reserve_aa_neighbor_pair(
    int first, int second, bool pinned) {
    if (!aa_neighbor_identity_enabled) return true;
    return aa_neighbor_slots.reserve(
        first, second, static_cast<aa_neighbor_identity::Step>(current_step),
        max_aa_bonds_per_actin, pinned);
}

int Sarcomere::aa_neighbor_degree(int actin_index) const {
    return aa_neighbor_identity_enabled
        ? aa_neighbor_slots.degree(actin_index) : 0;
}

// The axial frame of a filament pair does not depend on which ACTN2 site of
// that pair is being probed, so it is factored out of the ordinal loop in
// _update_local_aa_links().  Every expression below is copied verbatim from
// the original _local_aa_attachment_candidate() so the split is bitwise
// neutral.
Sarcomere::AAPairAxis Sarcomere::_local_aa_pair_axis(int i, int j) const {
    AAPairAxis geometry;
    vec ui = actin.direction[i];
    vec uj = actin.direction[j];
    ui.normalize();
    uj.normalize();
    if (ui.dot(uj) < 0.0) uj = -uj;
    vec axis = ui + uj;
    if (axis.norm() <= 1e-10) axis = ui;
    axis.normalize();
    const double pi = actin.center[i].dot(axis);
    vec center_displacement = actin.center[j] - actin.center[i];
    center_displacement.pbc_wrap(box, is_periodic);
    const double pj = pi + center_displacement.dot(axis);
    const double di = actin.direction[i].dot(axis);
    const double dj = actin.direction[j].dot(axis);
    const double hi = 0.5 * actin.length * std::abs(di);
    const double hj = 0.5 * actin.length * std::abs(dj);
    double lower = std::max(pi - hi, pj - hj);
    double upper = std::min(pi + hi, pj + hj);
    // ACTN2 is confined to the Z-disc: a zone of width aa_barbed_zone measured
    // from each filament's BARBED end, which is filament coordinate s = -L/2.
    // A crosslink can only form where BOTH filaments' zones overlap, so the
    // mutual overlap is intersected with each.  The surviving length is the
    // accessible binding length and scales the association rate through
    // overlap_extensive_factor; side-by-side pairs, whose zones never meet,
    // are left with zero accessible length and cannot bind at all.
    if (aa_barbed_zone > 0.0) {
        const double half = 0.5 * actin.length;
        const double ai = pi - half * di;
        const double bi = pi + (aa_barbed_zone - half) * di;
        const double aj = pj - half * dj;
        const double bj = pj + (aa_barbed_zone - half) * dj;
        lower = std::max(lower, std::max(std::min(ai, bi), std::min(aj, bj)));
        upper = std::min(upper, std::min(std::max(ai, bi), std::max(aj, bj)));
    }
    geometry.axis = axis;
    geometry.pi = pi;
    geometry.pj = pj;
    geometry.lower = lower;
    geometry.overlap = std::max(0.0, upper - lower);
    geometry.denominator_i = di;
    geometry.denominator_j = dj;
    return geometry;
}

std::pair<double, double> Sarcomere::_local_aa_attachment_site(
    const AAPairAxis& geometry, int ordinal, double& overlap) const {
    overlap = geometry.overlap;
    if (overlap <= 0.0) return {0.0, 0.0};
    // Place every coarse-grained ACTN2 molecule at its own material site.
    // The former two-slot expression repeated 0.35/0.65 for ordinal >= 2,
    // so increasing the molecular redundancy created duplicate springs rather
    // than additional binding sites.  Even spacing gives N distinct sites and
    // reduces to the same central pair (1/3, 2/3) when N == 2.
    const int site_count = std::max(1, aa_prebundle_links_per_pair);
    const int site_index = std::clamp(ordinal, 0, site_count - 1);
    const double fraction = static_cast<double>(site_index + 1) /
                            static_cast<double>(site_count + 1);
    const double target = geometry.lower + fraction * overlap;
    const double denominator_i = geometry.denominator_i;
    const double denominator_j = geometry.denominator_j;
    if (std::abs(denominator_i) <= 1e-8 || std::abs(denominator_j) <= 1e-8) {
        overlap = 0.0;
        return {0.0, 0.0};
    }
    const double s_i = std::clamp(
        (target - geometry.pi) / denominator_i,
        -0.5 * actin.length, 0.5 * actin.length);
    const double s_j = std::clamp(
        (target - geometry.pj) / denominator_j,
        -0.5 * actin.length, 0.5 * actin.length);
    return {s_i, s_j};
}

std::pair<double, double> Sarcomere::_local_aa_attachment_candidate(
    int i, int j, int ordinal, double& overlap) const {
    return _local_aa_attachment_site(
        _local_aa_pair_axis(i, j), ordinal, overlap);
}

bool Sarcomere::_local_aa_site_outside_am_overlap(
    int actin_index, double axial_coordinate) const {
    if (actin_index < 0 || actin_index >= actin.n) return false;
    // Legacy directional geometry defined the crosslinkable segment from the
    // barbed end (left_end, s=-L/2) up to the beginning of AM overlap.  With no
    // engaged myosin the ratio is one and the full actin remains available,
    // preserving prebundling before thick-filament incorporation.
    const double ratio = std::clamp(
        actin_crosslink_ratio[actin_index], 0.0, 1.0);
    if (ratio <= EPS) return false;
    const double allowed_end =
        -0.5 * actin.length + ratio * actin.length;
    return axial_coordinate <= allowed_end + 1e-12;
}

// Local evidence that this individual barbed end is part of a nascent Z-disc
// complex: at least one live ACTN2 molecule of any status holds it.  Status 1
// counts, because status 2 is gated on mature motor load and would recreate the
// original low-capture failure.  The per-actin index is rebuilt at the end of
// _update_local_aa_links(), which runs earlier in the same step than
// _update_explicit_titin(), so the buckets are current here.  A size mismatch
// means the local ACTN2 model is not active at all; that is reported as "no
// support" so the gate cannot silently pass on a stale index.
bool Sarcomere::actin_has_live_actn2(int actin_index) const {
    if (actin_index < 0 || actin_index >= actin.n) return false;
    if (aa_links_by_actin.size() != static_cast<size_t>(actin.n)) return false;
    for (size_t slot : aa_links_by_actin[static_cast<size_t>(actin_index)]) {
        if (slot >= aa_links.size() || !aa_links[slot].active) continue;
        const auto& link = aa_links[slot];
        if (!aa_arm_state_model ||
            (link.actin_i == actin_index &&
             actn2_arm_kinetics::arm_i_bound(link.arm_state)) ||
            (link.actin_j == actin_index &&
             actn2_arm_kinetics::arm_j_bound(link.arm_state))) {
            return true;
        }
    }
    return false;
}

int Sarcomere::_local_aa_distinct_neighbors(int actin_index) const {
    std::vector<int> neighbors;
    for (size_t slot : aa_links_by_actin[static_cast<size_t>(actin_index)]) {
        const auto& link = aa_links[slot];
        if (aa_arm_state_model &&
            !actn2_arm_kinetics::is_bridge(link.arm_state)) continue;
        neighbors.push_back(link.actin_i == actin_index
            ? link.actin_j : link.actin_i);
    }
    std::sort(neighbors.begin(), neighbors.end());
    neighbors.erase(std::unique(neighbors.begin(), neighbors.end()),
                    neighbors.end());
    return static_cast<int>(neighbors.size());
}

bool Sarcomere::_local_aa_site_available(
    int actin_index, int partner, double coordinate) const {
    for (size_t slot : aa_links_by_actin[static_cast<size_t>(actin_index)]) {
        const auto& link = aa_links[slot];
        if (!link.active) continue;
        // Tethers are indexed only on their still-bound actin, so their
        // surviving material site remains occupied while the free arm is
        // available to close.  Treating the entire tether as absent would let
        // a second molecule occupy the same pair-facing site during annealing.
        const int other = link.actin_i == actin_index
            ? link.actin_j : link.actin_i;
        // The reduced model has no ACTN2 azimuth.  Links to distinct
        // neighboring filaments may therefore share an axial level and
        // occupy different azimuthal faces; contour exclusion applies
        // only to duplicate links within the same filament pair.
        if (other != partner) continue;
        const double occupied = link.actin_i == actin_index ? link.s_i : link.s_j;
        if (std::abs(occupied - coordinate) < aa_site_exclusion) return false;
    }
    return true;
}

bool Sarcomere::_local_aa_neighbor_pair_available(
    int first, int second) const {
    return !aa_neighbor_identity_enabled ||
        aa_neighbor_slots.can_reserve(
            first, second, max_aa_bonds_per_actin);
}

std::vector<std::pair<double, double>>
Sarcomere::_local_aa_tether_capture_intervals(
    int bound_actin, int candidate_actin, double bound_coordinate) const {
    std::vector<std::pair<double, double>> intervals;
    if (bound_actin < 0 || bound_actin >= actin.n ||
        candidate_actin < 0 || candidate_actin >= actin.n ||
        bound_actin == candidate_actin) {
        return intervals;
    }

    // The free domain searches the candidate's continuous exposed contour,
    // not the canonical sites used to seed a new two-filament encounter.  Its
    // remaining arm stays at its material coordinate, making this a genuinely
    // local intramolecular capture problem.
    const double contour_lower = -0.5 * actin.length;
    double contour_upper = 0.5 * actin.length;
    if (aa_barbed_zone > 0.0) {
        contour_upper = std::min(
            contour_upper, contour_lower + aa_barbed_zone);
    }
    const double crosslink_ratio = std::clamp(
        actin_crosslink_ratio[candidate_actin], 0.0, 1.0);
    contour_upper = std::min(
        contour_upper, contour_lower + crosslink_ratio * actin.length);
    if (!(contour_upper > contour_lower)) return intervals;

    vec candidate_direction = actin.direction[candidate_actin];
    if (!(candidate_direction.norm() > 1e-12)) return intervals;
    candidate_direction.normalize();
    const vec bound_center = actin.center[bound_actin];
    const vec candidate_center = actin.center[candidate_actin];
    const vec bound_point = bound_center +
        bound_coordinate * actin.direction[bound_actin];
    vec center_from_bound = candidate_center - bound_point;
    center_from_bound.pbc_wrap(box, is_periodic);
    const double closest_coordinate =
        -center_from_bound.dot(candidate_direction);
    const vec transverse = center_from_bound +
        closest_coordinate * candidate_direction;
    const double radial_squared = transverse.dot(transverse);
    const double capture_squared =
        aa_tether_capture_distance * aa_tether_capture_distance;
    if (radial_squared >= capture_squared) return intervals;
    const double axial_reach = std::sqrt(std::max(
        0.0, capture_squared - radial_squared));
    const double reachable_lower = std::max(
        contour_lower, closest_coordinate - axial_reach);
    const double reachable_upper = std::min(
        contour_upper, closest_coordinate + axial_reach);
    if (!(reachable_upper > reachable_lower)) return intervals;
    intervals.push_back({reachable_lower, reachable_upper});

    // Remove pair-facing contour already occupied by another ACTN2 molecule.
    // Distinct neighboring filaments may still bind at the same axial level,
    // matching _local_aa_site_available().
    if (!(aa_site_exclusion > 0.0)) return intervals;
    for (size_t slot :
         aa_links_by_actin[static_cast<size_t>(candidate_actin)]) {
        const auto& link = aa_links[slot];
        if (!link.active) continue;
        const int other = link.actin_i == candidate_actin
            ? link.actin_j : link.actin_i;
        if (other != bound_actin) continue;
        const double occupied = link.actin_i == candidate_actin
            ? link.s_i : link.s_j;
        const double excluded_lower = occupied - aa_site_exclusion;
        const double excluded_upper = occupied + aa_site_exclusion;
        std::vector<std::pair<double, double>> remaining;
        remaining.reserve(intervals.size() + 1);
        for (const auto& interval : intervals) {
            if (excluded_upper <= interval.first ||
                excluded_lower >= interval.second) {
                remaining.push_back(interval);
                continue;
            }
            if (excluded_lower > interval.first) {
                remaining.push_back({
                    interval.first,
                    std::min(interval.second, excluded_lower)});
            }
            if (excluded_upper < interval.second) {
                remaining.push_back({
                    std::max(interval.first, excluded_upper),
                    interval.second});
            }
        }
        intervals.swap(remaining);
        if (intervals.empty()) break;
    }
    return intervals;
}

Sarcomere::AAMyosinExclusionEvaluation
Sarcomere::_local_aa_myosin_exclusion(
    int actin_i, int actin_j, double s_i, double s_j) const {
    AAMyosinExclusionEvaluation evaluation;
    const bool radial_contact_enabled =
        aa_actn2_contact_radius > 0.0 &&
        aa_myosin_contact_shell_width > 0.0 &&
        aa_myosin_contact_shell_energy > 0.0;
    if (!radial_contact_enabled && !aa_zdisc_tip_contact) {
        return evaluation;
    }
    const auto& first_neighbors =
        actin_neighbors_by_species[static_cast<size_t>(actin_i)].second;
    const auto& second_neighbors =
        actin_neighbors_by_species[static_cast<size_t>(actin_j)].second;
    const double hard_core_distance =
        myosin.radius + aa_actn2_contact_radius;
    const double radial_outer_distance = hard_core_distance +
        (radial_contact_enabled ? aa_myosin_contact_shell_width : 0.0);
    const auto patch = compute_actn2_patch_geometry(
        actin, actin_i, actin_j, s_i, s_j, box);
    auto evaluate_myosin = [&](int myosin_index) {
        const bool may_contact_body = radial_contact_enabled &&
            actn2_patch_may_contact_myosin_body(
                patch, myosin, myosin_index, box, radial_outer_distance);
        const bool may_contact_tip = aa_zdisc_tip_contact &&
            actn2_patch_may_contact_myosin_tip(
                patch, myosin, myosin_index, box, aa_zdisc_tip_onset,
                aa_zdisc_tip_transverse_cutoff);
        if (!may_contact_body && !may_contact_tip) return;
        if (may_contact_body) {
            const auto contact = compute_actn2_myosin_contact(
                actin, myosin, actin_i, actin_j, s_i, s_j, myosin_index,
                box, hard_core_distance, aa_myosin_contact_epsilon,
                aa_myosin_contact_force_cap, aa_myosin_contact_shell_width,
                aa_myosin_contact_shell_energy);
            evaluation.energy += contact.energy;
            evaluation.hard_core_overlap = evaluation.hard_core_overlap ||
                contact.penetration > 0.0;
        }
        // Steric accessibility is a property of the thick-filament body, not
        // of motor occupancy.  Motor engagement is retained only as an audit
        // field; it never decides whether ACTN2 may occupy this patch.
        bool motor_engages_pair = false;
        if (may_contact_tip && myosin_index >= 0 &&
            myosin_index < static_cast<int>(am_heads_by_myosin.size())) {
            for (const auto& half_heads :
                 am_heads_by_myosin[static_cast<size_t>(myosin_index)]) {
                for (const auto& head : half_heads) {
                    if (head.state == static_cast<BondFlag>(
                            am_motor_kinetics::State::detached)) continue;
                    if (head.bound_actin == actin_i ||
                        head.bound_actin == actin_j) {
                        motor_engages_pair = true;
                        break;
                    }
                }
                if (motor_engages_pair) break;
            }
        }
        if (may_contact_tip) {
            for (int end_sign : {-1, 1}) {
                const auto tip = compute_zdisc_myosin_tip_contact(
                    actin, myosin, actin_i, actin_j, s_i, s_j,
                    myosin_index, end_sign, box, aa_zdisc_tip_onset,
                    aa_zdisc_tip_stiffness, aa_zdisc_tip_force_max,
                    aa_zdisc_tip_transverse_flat,
                    aa_zdisc_tip_transverse_cutoff);
                if (tip.eligible_half &&
                    tip.transverse_distance < aa_zdisc_tip_transverse_cutoff &&
                    tip.axial_gap < evaluation.nearest_tip_gap) {
                    evaluation.nearest_tip_myosin = myosin_index;
                    evaluation.nearest_tip_sign = end_sign;
                    evaluation.nearest_tip_gap = tip.axial_gap;
                    evaluation.nearest_tip_rho = tip.transverse_distance;
                    evaluation.nearest_tip_motor_engaged =
                        motor_engages_pair;
                }
                if (!tip.applied) continue;
                evaluation.axial_tip_contact = true;
                evaluation.axial_tip_penetration = std::max(
                    evaluation.axial_tip_penetration, tip.compression);
                evaluation.axial_tip_contacts.push_back({
                    myosin_index, end_sign, tip.axial_gap,
                    tip.transverse_distance, tip.compression,
                    motor_engages_pair});
            }
        }
    };
    for (int myosin_index : first_neighbors) {
        evaluate_myosin(myosin_index);
    }
    for (int myosin_index : second_neighbors) {
        if (std::find(first_neighbors.begin(), first_neighbors.end(),
                      myosin_index) == first_neighbors.end()) {
            evaluate_myosin(myosin_index);
        }
    }
    return evaluation;
}

Sarcomere::AAMyosinExclusionEvaluation
Sarcomere::evaluate_local_aa_myosin_exclusion(
    int actin_i, int actin_j, double s_i, double s_j) const {
    return _local_aa_myosin_exclusion(actin_i, actin_j, s_i, s_j);
}

bool Sarcomere::aa_tip_diffusion_move_allowed(
    const AAMyosinExclusionEvaluation& old_evaluation,
    const AAMyosinExclusionEvaluation& proposed_evaluation,
    double tolerance) {
    tolerance = std::max(0.0, tolerance);
    for (const auto& proposed_contact :
         proposed_evaluation.axial_tip_contacts) {
        const auto old_contact = std::find_if(
            old_evaluation.axial_tip_contacts.begin(),
            old_evaluation.axial_tip_contacts.end(),
            [&](const AATipExclusionContact& contact) {
                return contact.myosin == proposed_contact.myosin &&
                       contact.end_sign == proposed_contact.end_sign;
            });
        if (old_contact == old_evaluation.axial_tip_contacts.end() ||
            proposed_contact.penetration >
                old_contact->penetration + tolerance) {
            return false;
        }
    }
    return true;
}

// Pure evaluation of one association proposal against the current link state.
// It performs no mutation, so it may run concurrently for many pairs; the
// caller commits the surviving proposals serially in canonical order and
// re-evaluates any pair whose actins were touched by an earlier commit.
// Necessary conditions for _evaluate_aa_association to return a firing
// proposal.  Both are geometry-only, so they cannot change as links commit
// during the step -- which is what makes it safe to drop a pair from the
// candidate list rather than re-test it in the serial loop.
//
//   1. it opens with
//        if (dir_i . dir_j >= aa_association_max_direction_dot) return {};
//   2. it ends the geometric section with
//        if (candidate.length >= aa_cutoff) return {};
//      where length is |connector| between two attachment points
//      centre + s*direction with |s| <= L/2, i.e. points ON the two segments.
//      So segment_segment_distance <= length, and a pair whose segments are
//      further apart than aa_cutoff can never satisfy it.
bool Sarcomere::_aa_pair_may_associate(int i, int j) const {
    if (actin.direction[i].dot(actin.direction[j]) >=
        aa_association_max_direction_dot) {
        return false;
    }
    return geometry::segments_may_interact(
        actin.left_end[i], actin.right_end[i],
        actin.left_end[j], actin.right_end[j],
        aa_cutoff, box, is_periodic);
}

Sarcomere::AAAssociationProposal Sarcomere::_evaluate_aa_association(
    int i, int j) const {
    AAAssociationProposal proposal;
    if (actin.direction[i].dot(actin.direction[j]) >=
        aa_association_max_direction_dot) {
        return proposal;
    }
    const auto key = _aa_pair_key(i, j);
    const auto pair_it = aa_links_by_pair.find(key);
    const int occupied = pair_it == aa_links_by_pair.end()
        ? 0
        : static_cast<int>(pair_it->second.size());
    if (occupied >= aa_prebundle_links_per_pair) return proposal;
    if (aa_total_link_pool > 0 &&
        static_cast<int>(aa_active_link_slots.size()) >= aa_total_link_pool) {
        return proposal;
    }
    const bool has_neighbor_slot =
        aa_neighbor_identity_enabled && aa_neighbor_slots.contains(i, j);
    if (max_aa_bonds_per_actin > 0 && occupied == 0 && !has_neighbor_slot) {
        const int degree_i = aa_neighbor_identity_enabled
            ? aa_neighbor_slots.degree(i) : _local_aa_distinct_neighbors(i);
        const int degree_j = aa_neighbor_identity_enabled
            ? aa_neighbor_slots.degree(j) : _local_aa_distinct_neighbors(j);
        if (degree_i >= max_aa_bonds_per_actin ||
            degree_j >= max_aa_bonds_per_actin) return proposal;
    }
    double overlap = 0.0;
    std::pair<double, double> attachment{0.0, 0.0};
    bool found_site = false;
    // Occupancy alone does not identify which of the two axial slots
    // remains after a detachment.  Probe canonical slots explicitly
    // so either missing site can repopulate.  The pair axis is shared by
    // every ordinal, so it is built once.
    const AAPairAxis geometry = _local_aa_pair_axis(i, j);
    for (int ordinal = 0; ordinal < aa_prebundle_links_per_pair; ++ordinal) {
        double candidate_overlap = 0.0;
        const auto candidate_attachment =
            _local_aa_attachment_site(geometry, ordinal, candidate_overlap);
        if (candidate_overlap <= 0.0) continue;
        if (!_local_aa_site_available(i, j, candidate_attachment.first) ||
            !_local_aa_site_available(j, i, candidate_attachment.second)) {
            continue;
        }
        if (!_local_aa_site_outside_am_overlap(
                i, candidate_attachment.first) ||
            !_local_aa_site_outside_am_overlap(
                j, candidate_attachment.second)) {
            continue;
        }
        overlap = candidate_overlap;
        attachment = candidate_attachment;
        found_site = true;
        break;
    }
    if (!found_site) return proposal;
    const auto candidate = compute_local_aa_link(
        actin, i, j, attachment.first, attachment.second, box,
        aa_encounter_k, aa_encounter_kappa, 0.0, aa_optimal);
    if (candidate.length >= aa_cutoff) return proposal;
    const auto exclusion = _local_aa_myosin_exclusion(
        i, j, attachment.first, attachment.second);
    proposal.nearest_tip_myosin = exclusion.nearest_tip_myosin;
    proposal.nearest_tip_sign = exclusion.nearest_tip_sign;
    proposal.nearest_tip_gap = exclusion.nearest_tip_gap;
    proposal.nearest_tip_rho = exclusion.nearest_tip_rho;
    proposal.nearest_tip_motor_engaged =
        exclusion.nearest_tip_motor_engaged;
    // A molecular bridge cannot be born through the thick-filament hard core.
    // Configurations in the soft corona remain possible but pay their actual
    // conservative Boltzmann weight below.
    if (exclusion.hard_core_overlap || exclusion.axial_tip_contact) {
        return proposal;
    }
    const double overlap_factor =
        aa_association_kinetics::overlap_extensive_factor(
            overlap, aa_overlap_capture);
    double association_rate =
        (aa_prebundle_links_per_pair - occupied) * k_on * overlap_factor;
    if (aa_total_link_pool > 0) {
        association_rate *= std::max(0.0,
            1.0 - static_cast<double>(aa_active_link_slots.size()) /
                      static_cast<double>(aa_total_link_pool));
    }
    if (aa_energy_gated_association || aa_passive_detailed_balance) {
        association_rate = aa_association_kinetics::transition_state_on_rate(
            association_rate, thermal_beta, candidate.energy,
            aa_transition_state_theta);
    }
    association_rate *= boltzmann_exclusion_factor(
        thermal_beta, exclusion.energy);
    counter_rng::CounterRng association_rng(
        static_cast<std::uint64_t>(initial_seed),
        static_cast<std::uint64_t>(current_step),
        counter_rng::Stream::aa_pair_association,
        counter_rng::canonical_pair_key(i, j));
    proposal.s_i = attachment.first;
    proposal.s_j = attachment.second;
    proposal.fire = association_rng.uniform() <
        aa_association_kinetics::event_probability(association_rate, dt);
    return proposal;
}

const char* Sarcomere::aa_phase_name(int phase) {
    switch (phase) {
        case kAaRebuildIndices: return "rebuild_indices";
        case kAaArmKinetics: return "arm_kinetics";
        case kAaSiteDiffusion: return "site_diffusion";
        case kAaMotorForces: return "motor_forces";
        case kAaArmLoads: return "arm_loads";
        case kAaLinkLoop: return "link_loop_ordered";
        case kAaReindex: return "reindex_expire";
        case kAaSortNeighbors: return "sort_neighbors";
        case kAaPrefixSum: return "pair_prefix_sum";
        case kAaScatterPairs: return "scatter_pairs";
        case kAaEvalProposals: return "evaluate_proposals";
        case kAaCommitLinks: return "commit_links";
        case kAaTipContact: return "tip_contact";
        case kAaTipAndScatter: return "mechanics_scatter";
        default: return "unknown";
    }
}

void Sarcomere::_update_local_aa_links() {
    auto& phase_seconds = aa_phase_seconds[omp_get_thread_num()];
    auto phase_mark = std::chrono::steady_clock::now();
    int phase_index = 0;
    const auto stamp_aa_phase = [&](int next) {
        const auto now = std::chrono::steady_clock::now();
        phase_seconds[phase_index] +=
            std::chrono::duration<double>(now - phase_mark).count();
        phase_mark = now;
        phase_index = next;
    };
    // Called from inside the `omp parallel` region of update_system() by every
    // thread of the team.  The routine is a sequence of orphaned worksharing
    // constructs separated by `single` blocks that own the shared mutations
    // (index rebuilds, link creation/destruction, event emission).  All
    // cross-phase state therefore lives in aa_scratch_* members rather than
    // function locals, which would otherwise be private per thread.
    const double ramp = _current_interaction_ramp();
    const double single_link_alignment_k = ramp *
        actn2_load_kinetics::fixed_single_link_stiffness(
            kappa_aa, aa_stiffness_reference_links);
    const bool motor_transmitted_mode =
        aa_status_mode ==
        actn2_load_kinetics::StatusMode::motor_transmitted_register;

    stamp_aa_phase(kAaRebuildIndices);
    #pragma omp single
    {
    if (motor_transmitted_mode) aa_motor_update_started = omp_get_wtime();
    _rebuild_local_aa_indices();
    // Existing single-arm tethers turn over before bridge mechanics.  A free
    // arm closes by a local KMC event; the bound material coordinate never
    // slides or teleports.  Bridges that open later in this step first become
    // eligible to rebind on the next step, giving every molecular transition a
    // well-defined order.
    if (aa_arm_state_model) {
        const auto tether_slots = aa_active_link_slots;
        for (size_t slot : tether_slots) {
            auto& link = aa_links[slot];
            if (!link.active ||
                !actn2_arm_kinetics::is_tether(link.arm_state)) continue;
            const auto old_state = link.arm_state;
            const bool bound_to_i = old_state ==
                actn2_arm_kinetics::State::tethered_to_i;
            const int bound_actin = bound_to_i ? link.actin_i : link.actin_j;
            const int old_partner = bound_to_i ? link.actin_j : link.actin_i;
            const double bound_s = bound_to_i ? link.s_i : link.s_j;
            counter_rng::CounterRng detach_rng(
                static_cast<std::uint64_t>(initial_seed),
                static_cast<std::uint64_t>(current_step),
                counter_rng::Stream::aa_tether_detachment, link.id);
            const double off_probability =
                actn2_arm_kinetics::event_probability(aa_tether_off_rate, dt);
            const bool stochastic_detachment =
                detach_rng.uniform() < off_probability;
            const bool no_crosslinkable_contour =
                actn2_arm_kinetics::contour_fully_occluded(
                    actin_crosslink_ratio[bound_actin], EPS);
            if (stochastic_detachment || no_crosslinkable_contour) {
                aa_arm_events.insert(aa_arm_events.end(), {
                    static_cast<double>(link.id),
                    static_cast<double>(current_step), 3.0,
                    static_cast<double>(static_cast<int>(old_state)), 0.0,
                    static_cast<double>(bound_actin),
                    static_cast<double>(old_partner), -1.0, bound_s,
                    std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0, 0.0,
                    no_crosslinkable_contour
                        ? std::numeric_limits<double>::infinity()
                        : aa_tether_off_rate});
                link.active = false;
                aa_free_link_slots.push_back(slot);
                continue;
            }

            struct RebindCandidate {
                int actin = -1;
                double lower_s = 0.0;
                double upper_s = 0.0;
                double rate = 0.0;
                double energy = 0.0;
            };
            std::vector<RebindCandidate> candidates;
            double total_rate = 0.0;
            for (int candidate_actin :
                 actin_neighbors_by_species[static_cast<size_t>(bound_actin)].first) {
                if (candidate_actin == bound_actin) continue;
                if (actin.direction[bound_actin].dot(
                        actin.direction[candidate_actin]) >=
                    aa_association_max_direction_dot) continue;
                if (!_local_aa_neighbor_pair_available(
                        bound_actin, candidate_actin)) continue;
                const auto pair_it = aa_links_by_pair.find(
                    _aa_pair_key(bound_actin, candidate_actin));
                if (pair_it != aa_links_by_pair.end() &&
                    static_cast<int>(pair_it->second.size()) >=
                        aa_prebundle_links_per_pair) continue;
                const auto capture_intervals =
                    _local_aa_tether_capture_intervals(
                        bound_actin, candidate_actin, bound_s);
                for (const auto& interval : capture_intervals) {
                    const double interval_width =
                        interval.second - interval.first;
                    if (!(interval_width > 0.0)) continue;
                    // Midpoint quadrature resolves the continuous contour at
                    // no coarser than 10 nm (and at half the site-exclusion
                    // scale when that is smaller).  Each bin carries its
                    // physical contour measure, so refining the quadrature
                    // cannot multiply the molecular attempt rate.
                    const double bin_target = std::max(
                        1e-6, std::min(
                            0.010, aa_site_exclusion > 0.0
                                ? 0.5 * aa_site_exclusion : 0.010));
                    const int bin_count = std::max(
                        1, static_cast<int>(std::ceil(
                            interval_width / bin_target)));
                    const double bin_width = interval_width / bin_count;
                    for (int bin = 0; bin < bin_count; ++bin) {
                    const double lower_s = interval.first + bin * bin_width;
                    const double upper_s = lower_s + bin_width;
                    const double free_s = 0.5 * (lower_s + upper_s);
                    const auto candidate_mechanics = compute_local_aa_link(
                        actin,
                        bound_to_i ? bound_actin : candidate_actin,
                        bound_to_i ? candidate_actin : bound_actin,
                        bound_to_i ? bound_s : free_s,
                        bound_to_i ? free_s : bound_s,
                        box, aa_link_radial_stiffness(static_cast<BondFlag>(1)),
                        aa_link_perpendicular_stiffness(static_cast<BondFlag>(1)),
                        0.0, aa_optimal, 0.0, 0.0,
                        aa_register_saturation_error,
                        aa_axial_register_force_max, aa_radial_force_max);
                    if (candidate_mechanics.length >
                        aa_tether_capture_distance) continue;
                    const int candidate_i = bound_to_i
                        ? bound_actin : candidate_actin;
                    const int candidate_j = bound_to_i
                        ? candidate_actin : bound_actin;
                    const double candidate_s_i = bound_to_i
                        ? bound_s : free_s;
                    const double candidate_s_j = bound_to_i
                        ? free_s : bound_s;
                    const auto exclusion = _local_aa_myosin_exclusion(
                        candidate_i, candidate_j,
                        candidate_s_i, candidate_s_j);
                    if (exclusion.hard_core_overlap ||
                        exclusion.axial_tip_contact) continue;
                    const double local_closing_rate =
                        actn2_arm_kinetics::elastic_rebinding_rate(
                            aa_tether_rebind_rate,
                            candidate_mechanics.energy, thermal_beta,
                            aa_arm_rebind_theta);
                    const double rate = local_closing_rate *
                        boltzmann_exclusion_factor(
                            thermal_beta, exclusion.energy) *
                        bin_width / std::max(1e-12, aa_overlap_capture);
                    if (!(rate > 0.0)) continue;
                    candidates.push_back({candidate_actin, lower_s, upper_s, rate,
                                          candidate_mechanics.energy});
                    total_rate += rate;
                    }
                }
            }
            if (!(total_rate > 0.0)) continue;
            counter_rng::CounterRng rebind_rng(
                static_cast<std::uint64_t>(initial_seed),
                static_cast<std::uint64_t>(current_step),
                counter_rng::Stream::aa_arm_rebinding, link.id);
            if (rebind_rng.uniform() >=
                actn2_arm_kinetics::event_probability(total_rate, dt)) continue;
            double target = rebind_rng.uniform() * total_rate;
            const RebindCandidate* selected = &candidates.back();
            for (const auto& candidate : candidates) {
                target -= candidate.rate;
                if (target <= 0.0) {
                    selected = &candidate;
                    break;
                }
            }
            const double selected_free_s = selected->lower_s +
                (selected->upper_s - selected->lower_s) * rebind_rng.uniform();
            const int selected_i = bound_to_i
                ? bound_actin : selected->actin;
            const int selected_j = bound_to_i
                ? selected->actin : bound_actin;
            const double selected_s_i = bound_to_i
                ? bound_s : selected_free_s;
            const double selected_s_j = bound_to_i
                ? selected_free_s : bound_s;
            const auto selected_exclusion = _local_aa_myosin_exclusion(
                selected_i, selected_j, selected_s_i, selected_s_j);
            if (selected_exclusion.hard_core_overlap ||
                selected_exclusion.axial_tip_contact) continue;
            if (aa_zdisc_tip_contact &&
                selected_exclusion.nearest_tip_myosin >= 0 &&
                selected_exclusion.nearest_tip_gap <
                    aa_zdisc_tip_onset - 1e-12) {
                throw std::runtime_error(
                    "ACTN2 tether rebinding attempted to commit inside a "
                    "myosin-tip exclusion zone");
            }
            if (aa_neighbor_identity_enabled &&
                !aa_neighbor_slots.reserve(
                    bound_actin, selected->actin,
                    static_cast<aa_neighbor_identity::Step>(current_step),
                    max_aa_bonds_per_actin)) {
                continue;
            }
            link.last_partner = old_partner;
            if (bound_to_i) {
                link.actin_j = selected->actin;
                link.s_j = selected_free_s;
            } else {
                link.actin_i = selected->actin;
                link.s_i = selected_free_s;
            }
            link.status = 1;
            link.arm_state = actn2_arm_kinetics::State::weak_bridge;
            link.arm_state_since = static_cast<RecoveryStep>(std::min(
                current_step,
                static_cast<size_t>(std::numeric_limits<RecoveryStep>::max())));
            link.attach_step = link.arm_state_since;
            // Make this closure visible to later tether events in the same
            // serial KMC sweep.  Without the incremental index update, two
            // free arms could select the same pair-facing contour before the
            // end-of-phase rebuild.
            aa_links_by_pair[_aa_pair_key(
                bound_actin, selected->actin)].push_back(slot);
            aa_links_by_actin[static_cast<size_t>(selected->actin)]
                .push_back(slot);
            aa_tip_accessibility_events.insert(
                aa_tip_accessibility_events.end(), {
                    static_cast<double>(link.id),
                    static_cast<double>(current_step), 2.0,
                    static_cast<double>(link.actin_i),
                    static_cast<double>(link.actin_j), link.s_i, link.s_j,
                    static_cast<double>(
                        selected_exclusion.nearest_tip_myosin),
                    static_cast<double>(
                        selected_exclusion.nearest_tip_sign),
                    selected_exclusion.nearest_tip_gap,
                    selected_exclusion.nearest_tip_rho, 0.0,
                    selected_exclusion.nearest_tip_motor_engaged ? 1.0
                                                                  : 0.0});
            aa_arm_events.insert(aa_arm_events.end(), {
                static_cast<double>(link.id),
                static_cast<double>(current_step), 2.0,
                static_cast<double>(static_cast<int>(old_state)),
                static_cast<double>(static_cast<int>(link.arm_state)),
                static_cast<double>(bound_actin),
                static_cast<double>(old_partner),
                static_cast<double>(selected->actin), bound_s,
                selected_free_s, 0.0, 0.0, 0.0, selected->rate});
        }
        _rebuild_local_aa_indices();
    }
    aa_scratch_pair_tension.clear();
    aa_scratch_mechanics.assign(aa_links.size(), LocalAALinkResult{});
    aa_scratch_arm_loads.assign(aa_links.size(),
                                actn2_load_kinetics::ArmLoads{});
    aa_scratch_contact_arm_loads.assign(
        aa_links.size(), std::array<double, 2>{0.0, 0.0});
    aa_scratch_contact_penetration.assign(aa_links.size(), 0.0);
    aa_scratch_contact_results.assign(aa_links.size(), {});
    aa_scratch_contact_evaluated.assign(aa_links.size(), 0);
    aa_scratch_total_arm_load.assign(static_cast<size_t>(actin.n), 0.0);
    aa_scratch_productive_load.assign(static_cast<size_t>(actin.n), 0.0);
    aa_scratch_motor_force.assign(static_cast<size_t>(actin.n),
                                  std::array<double, 3>{0.0, 0.0, 0.0});
    aa_scratch_productive_force.assign(static_cast<size_t>(actin.n),
                                       std::array<double, 3>{0.0, 0.0, 0.0});
    if (aa_scratch_productive_by_myosin.size() !=
        static_cast<size_t>(actin.n)) {
        aa_scratch_productive_by_myosin.resize(static_cast<size_t>(actin.n));
    }
    if (aa_scratch_productive_source_demand.size() !=
        static_cast<size_t>(actin.n)) {
        aa_scratch_productive_source_demand.resize(
            static_cast<size_t>(actin.n));
    }
    }

    auto& mechanics = aa_scratch_mechanics;
    auto& raw_arm_loads = aa_scratch_arm_loads;
    auto& contact_arm_loads = aa_scratch_contact_arm_loads;
    auto& contact_penetration = aa_scratch_contact_penetration;
    auto& contact_results = aa_scratch_contact_results;
    auto& contact_evaluated = aa_scratch_contact_evaluated;
    auto& pair_tension_sum = aa_scratch_pair_tension;
    auto& total_raw_favorable_arm_load = aa_scratch_total_arm_load;
    auto& productive_motor_load = aa_scratch_productive_load;
    const auto& active_slots = aa_active_link_slots;
    const long long n_active = static_cast<long long>(active_slots.size());

    // A weak ACTN2 encounter complex is not stereospecifically locked: each
    // binding domain diffuses in one dimension along its actin's barbed-end
    // association zone. The connector force supplies the Einstein drift.
    // Counter-based noise keeps each physical molecule's trajectory invariant
    // to OpenMP scheduling. Status-2 coordinates are deliberately untouched.
    if (!aa_arm_state_model && aa_status1_site_diffusion > 0.0 &&
        thermal_beta > 0.0 && dt > 0.0) {
        const double lower = -0.5 * actin.length;
        const double zone = aa_barbed_zone > 0.0
            ? std::min(actin.length, aa_barbed_zone) : actin.length;
        const double upper = lower + zone;
        const double sigma = std::sqrt(2.0 * aa_status1_site_diffusion * dt);
        auto reflect = [lower, upper](double value) {
            if (!(upper > lower)) return lower;
            const double width = upper - lower;
            double folded = std::fmod(value - lower, 2.0 * width);
            if (folded < 0.0) folded += 2.0 * width;
            return lower + (folded <= width ? folded
                                             : 2.0 * width - folded);
        };
        stamp_aa_phase(kAaArmKinetics);
        #pragma omp for schedule(static)
        for (long long index = 0; index < n_active; ++index) {
            const size_t slot = active_slots[static_cast<size_t>(index)];
            auto& link = aa_links[slot];
            if (!link.active || link.status != 1) continue;
            const auto weak_mechanics = compute_local_aa_link(
                actin, link.actin_i, link.actin_j, link.s_i, link.s_j, box,
                ramp * aa_link_radial_stiffness(link.status),
                ramp * aa_link_perpendicular_stiffness(link.status), single_link_alignment_k,
                aa_optimal, aa_antiparallel_overlap,
                aa_link_register_ratio(link.status),
                aa_register_saturation_error, aa_axial_register_force_max,
            ramp * aa_radial_force_max);
            counter_rng::CounterRng site_rng(
                static_cast<std::uint64_t>(initial_seed),
                static_cast<std::uint64_t>(current_step),
                counter_rng::Stream::aa_status1_site_diffusion,
                static_cast<std::uint64_t>(link.id));
            const double radius = std::sqrt(
                -2.0 * std::log(site_rng.uniform_pos()));
            const double angle = 6.2831853071795864769 * site_rng.uniform();
            const double gaussian_i = radius * std::cos(angle);
            const double gaussian_j = radius * std::sin(angle);
            const vec force_i = weak_mechanics.force_on_first;
            const vec force_j = -force_i;
            const double drift_i = aa_status1_site_diffusion * thermal_beta *
                force_i.dot(actin.direction[link.actin_i]) * dt;
            const double drift_j = aa_status1_site_diffusion * thermal_beta *
                force_j.dot(actin.direction[link.actin_j]) * dt;
            const double proposed_i = reflect(
                link.s_i + drift_i + sigma * gaussian_i);
            const double proposed_j = reflect(
                link.s_j + drift_j + sigma * gaussian_j);
            const auto old_exclusion = _local_aa_myosin_exclusion(
                link.actin_i, link.actin_j, link.s_i, link.s_j);
            const auto proposed_exclusion = _local_aa_myosin_exclusion(
                link.actin_i, link.actin_j, proposed_i, proposed_j);
            // A no-flux Brownian boundary is enforced independently for each
            // physical myosin tip. A move cannot enter a new tip zone, and a
            // site already inside a zone may only reduce that same tip's
            // penetration. Keying by (myosin,end_sign) prevents a move from
            // leaving one excluded zone by silently entering another.
            const bool exits_or_reduces_axial_contact =
                aa_tip_diffusion_move_allowed(
                    old_exclusion, proposed_exclusion);
            if (!proposed_exclusion.hard_core_overlap &&
                exits_or_reduces_axial_contact) {
                link.s_i = proposed_i;
                link.s_j = proposed_j;
            }
        }
    }

    // Phase 0.  The motor-force vectors of an actin are read-only for the whole
    // routine, so caching them per actin replaces the former two evaluations
    // per incident link (each of which rescanned every head of every bound
    // myosin) with one evaluation per actin.
    double motor_support_build_thread_seconds = 0.0;
    stamp_aa_phase(kAaSiteDiffusion);
    #pragma omp for schedule(static)
    for (int i = 0; i < actin.n; ++i) {
        aa_scratch_motor_force[static_cast<size_t>(i)] =
            _actin_motor_force_vector(i);
        if (motor_transmitted_mode) {
            const double started = omp_get_wtime();
            auto& contributions =
                aa_scratch_productive_by_myosin[static_cast<size_t>(i)];
            contributions = _actin_productive_motor_contributions(i);
            aa_scratch_productive_source_demand[
                static_cast<size_t>(i)].assign(contributions.size(), 0.0);
            std::array<double, 3> aggregate{0.0, 0.0, 0.0};
            for (const auto& contribution : contributions) {
                aggregate[0] += contribution.force[0];
                aggregate[1] += contribution.force[1];
                aggregate[2] += contribution.force[2];
            }
            aa_scratch_productive_force[static_cast<size_t>(i)] = aggregate;
            motor_support_build_thread_seconds += omp_get_wtime() - started;
        } else {
            aa_scratch_productive_by_myosin[static_cast<size_t>(i)].clear();
            aa_scratch_productive_source_demand[
                static_cast<size_t>(i)].clear();
            aa_scratch_productive_force[static_cast<size_t>(i)] =
                _actin_productive_motor_force_vector(i);
        }
        productive_motor_load[static_cast<size_t>(i)] =
            actn2_load_kinetics::norm(
                aa_scratch_productive_force[static_cast<size_t>(i)]);
    }
    if (motor_transmitted_mode) {
        #pragma omp atomic update
        aa_motor_support_build_seconds +=
            motor_support_build_thread_seconds;
    }

    // Phase 1a.  Unloaded mechanics and signed arm loads are independent per
    // link, so they are computed by the whole team.
    stamp_aa_phase(kAaMotorForces);
    #pragma omp for schedule(static)
    for (long long index = 0; index < n_active; ++index) {
        const size_t slot = active_slots[static_cast<size_t>(index)];
        const auto& link = aa_links[slot];
        if (aa_arm_state_model &&
            !actn2_arm_kinetics::is_bridge(link.arm_state)) {
            mechanics[slot] = LocalAALinkResult{};
            raw_arm_loads[slot] = actn2_load_kinetics::ArmLoads{};
            continue;
        }
        const double link_radial_k =
            ramp * aa_link_radial_stiffness(link.status);
        mechanics[slot] = compute_local_aa_link(
            actin, link.actin_i, link.actin_j, link.s_i, link.s_j, box,
            link_radial_k, ramp * aa_link_perpendicular_stiffness(link.status),
            single_link_alignment_k, aa_optimal, aa_antiparallel_overlap,
            aa_link_register_ratio(link.status),
            aa_register_saturation_error, aa_axial_register_force_max,
            ramp * aa_radial_force_max);
        const std::array<double, 3> link_force{
            mechanics[slot].force_on_first.x,
            mechanics[slot].force_on_first.y,
            mechanics[slot].force_on_first.z};
        const std::array<double, 3> direction_i{
            actin.direction[link.actin_i].x,
            actin.direction[link.actin_i].y,
            actin.direction[link.actin_i].z};
        const std::array<double, 3> direction_j{
            actin.direction[link.actin_j].x,
            actin.direction[link.actin_j].y,
            actin.direction[link.actin_j].z};
        raw_arm_loads[slot] = actn2_load_kinetics::signed_arm_loads(
            link_force, direction_i, direction_j);
        if ((aa_myosin_contact_epsilon > 0.0 ||
             (aa_myosin_contact_shell_width > 0.0 &&
              aa_myosin_contact_shell_energy > 0.0)) &&
            aa_actn2_contact_radius > 0.0) {
            contact_evaluated[slot] = 1;
            std::vector<int> contact_myosins =
                actin_neighbors_by_species[static_cast<size_t>(link.actin_i)].second;
            const auto& second_neighbors =
                actin_neighbors_by_species[static_cast<size_t>(link.actin_j)].second;
            contact_myosins.insert(contact_myosins.end(),
                                    second_neighbors.begin(), second_neighbors.end());
            std::sort(contact_myosins.begin(), contact_myosins.end());
            contact_myosins.erase(
                std::unique(contact_myosins.begin(), contact_myosins.end()),
                contact_myosins.end());
            const auto patch = compute_actn2_patch_geometry(
                actin, link.actin_i, link.actin_j, link.s_i, link.s_j, box);
            const double outer_contact_distance =
                myosin.radius + aa_actn2_contact_radius +
                aa_myosin_contact_shell_width;
            for (int myosin_index : contact_myosins) {
                if (!actn2_patch_may_contact_myosin_body(
                        patch, myosin, myosin_index, box,
                        outer_contact_distance)) continue;
                const auto contact = compute_actn2_myosin_contact(
                    actin, myosin, link.actin_i, link.actin_j,
                    link.s_i, link.s_j, myosin_index, box,
                    myosin.radius + aa_actn2_contact_radius,
                    aa_myosin_contact_epsilon,
                    aa_myosin_contact_force_cap,
                    aa_myosin_contact_shell_width,
                    aa_myosin_contact_shell_energy);
                if (!contact.applied) continue;
                contact_results[slot].push_back({myosin_index, contact});
                contact_arm_loads[slot][0] +=
                    contact.arm_weight_i * contact.force;
                contact_arm_loads[slot][1] +=
                    contact.arm_weight_j * contact.force;
                contact_penetration[slot] = std::max(
                    contact_penetration[slot], contact.penetration);
            }
        }
    }

    // Phase 1b.  The two accumulators are floating-point reductions over
    // overlapping actin/pair keys.  They stay serial and in slot order so the
    // summation order is exactly the original one; the work is a handful of
    // adds per link and does not scale with head count.
    stamp_aa_phase(kAaArmLoads);
    #pragma omp single
    {
    for (size_t slot : active_slots) {
        const auto& link = aa_links[slot];
        if (aa_arm_state_model &&
            !actn2_arm_kinetics::is_bridge(link.arm_state)) continue;
        pair_tension_sum[_aa_pair_key(link.actin_i, link.actin_j)] +=
            std::max(0.0, ramp * aa_link_radial_stiffness(link.status) *
                (mechanics[slot].length - aa_optimal));
        total_raw_favorable_arm_load[static_cast<size_t>(link.actin_i)] +=
            std::max(0.0, raw_arm_loads[slot].first);
        total_raw_favorable_arm_load[static_cast<size_t>(link.actin_j)] +=
            std::max(0.0, raw_arm_loads[slot].second);
        if (motor_transmitted_mode) {
            const std::array<double, 3> connector_i{
                mechanics[slot].force_on_first.x,
                mechanics[slot].force_on_first.y,
                mechanics[slot].force_on_first.z};
            const std::array<double, 3> connector_j{
                -connector_i[0], -connector_i[1], -connector_i[2]};
            const auto& sources_i = aa_scratch_productive_by_myosin[
                static_cast<size_t>(link.actin_i)];
            const auto& sources_j = aa_scratch_productive_by_myosin[
                static_cast<size_t>(link.actin_j)];
            auto& demand_i = aa_scratch_productive_source_demand[
                static_cast<size_t>(link.actin_i)];
            auto& demand_j = aa_scratch_productive_source_demand[
                static_cast<size_t>(link.actin_j)];
            for (size_t source = 0; source < sources_i.size(); ++source) {
                demand_i[source] += actn2_load_kinetics::opposing_support(
                    sources_i[source].force, connector_i);
            }
            for (size_t source = 0; source < sources_j.size(); ++source) {
                demand_j[source] += actn2_load_kinetics::opposing_support(
                    sources_j[source].force, connector_j);
            }
        }
    }
    }

    // Phase 2.  Detachment and catch-state update are evaluated per physical
    // link. Each actin's productive motor-force budget is shared across every
    // favorable ACTN2 arm incident on that actin, including links to distinct
    // partners.  Everything up to the `ordered` region is independent per link
    // and runs concurrently; link destruction, free-slot recycling and event
    // emission are order-sensitive and run inside `ordered`, which replays them
    // in exactly the original slot sequence.
    double motor_candidate_thread_seconds = 0.0;
    std::uint64_t motor_candidate_thread_pairs = 0;
    stamp_aa_phase(kAaLinkLoop);
    // Serial by measurement, not by necessity: the loop holds only ~0.1 s of
    // work per 100 steps at production density, and `ordered` lock-step cost
    // 0.27 s on top of that at 8-32 threads.  Too cheap to parallelise.
    #pragma omp single
    {
    for (long long index = 0; index < n_active; ++index) {
        const size_t slot = active_slots[static_cast<size_t>(index)];
        auto& link = aa_links[slot];
        if (aa_arm_state_model &&
            !actn2_arm_kinetics::is_bridge(link.arm_state)) {
            continue;
        }
        const auto& mech = mechanics[slot];
        // This scalar connector tension is retained only for the legacy
        // shared-load model.  The production arm_load_productive classifier
        // below uses the signed force delivered to each actin arm together
        // with independently measured opposed productive motor support.
        const double link_tension = std::max(
            0.0, ramp * aa_link_radial_stiffness(link.status) *
                (mech.length - aa_optimal));
        const auto& motor_i =
            aa_scratch_motor_force[static_cast<size_t>(link.actin_i)];
        const auto& motor_j =
            aa_scratch_motor_force[static_cast<size_t>(link.actin_j)];
        const double dx = motor_i[0] - motor_j[0];
        const double dy = motor_i[1] - motor_j[1];
        const double dz = motor_i[2] - motor_j[2];
        const double opposed = 0.5 * std::sqrt(dx * dx + dy * dy + dz * dz);
        // Read-only lookup: phase 1b already inserted every active pair key,
        // and `find` (unlike `operator[]`) cannot mutate the shared map from
        // several threads at once.
        const auto tension_it =
            pair_tension_sum.find(_aa_pair_key(link.actin_i, link.actin_j));
        const double denominator = tension_it == pair_tension_sum.end()
            ? 0.0 : tension_it->second;
        const double weight = denominator > 0.0
            ? link_tension / denominator : 0.0;
        const auto arm_loads = raw_arm_loads[slot];
        double motor_transmitted_support = 0.0;
        bool motor_pair_committed = false;
        const actn2_load_kinetics::MotorPairIdentity old_candidate_identity{
            link.candidate_myosin_i, link.candidate_myosin_j,
            link.candidate_half_i, link.candidate_half_j};
        const bool candidate_was_locked = link.candidate_locked;
        actn2_load_kinetics::MotorPairEvent motor_pair_event =
            actn2_load_kinetics::MotorPairEvent::none;
        double motor_pair_best_support = 0.0;
        double motor_pair_current_support = 0.0;
        double motor_pair_candidate_age = 0.0;
        if (motor_transmitted_mode) {
            const double started = omp_get_wtime();
            const std::array<double, 3> connector_i{
                mech.force_on_first.x, mech.force_on_first.y,
                mech.force_on_first.z};
            const std::array<double, 3> connector_j{
                -connector_i[0], -connector_i[1], -connector_i[2]};
            const auto& sources_i = aa_scratch_productive_by_myosin[
                static_cast<size_t>(link.actin_i)];
            const auto& sources_j = aa_scratch_productive_by_myosin[
                static_cast<size_t>(link.actin_j)];
            const auto& demand_i = aa_scratch_productive_source_demand[
                static_cast<size_t>(link.actin_i)];
            const auto& demand_j = aa_scratch_productive_source_demand[
                static_cast<size_t>(link.actin_j)];
            double best_support = 0.0;
            double current_support = 0.0;
            double best_raw_i = 0.0;
            double best_raw_j = 0.0;
            double best_allocated_i = 0.0;
            double best_allocated_j = 0.0;
            double current_raw_i = 0.0;
            double current_raw_j = 0.0;
            double current_allocated_i = 0.0;
            double current_allocated_j = 0.0;
            int best_myosin_i = -1;
            int best_myosin_j = -1;
            int8_t best_half_i = 0;
            int8_t best_half_j = 0;
            auto identity_less = [](int mi, int8_t hi, int mj, int8_t hj,
                                    int other_mi, int8_t other_hi,
                                    int other_mj, int8_t other_hj) {
                return std::tie(mi, hi, mj, hj) <
                       std::tie(other_mi, other_hi, other_mj, other_hj);
            };
            for (size_t source_i_index = 0;
                 source_i_index < sources_i.size(); ++source_i_index) {
                const auto& source_i = sources_i[source_i_index];
                const double raw_support_i =
                    actn2_load_kinetics::opposing_support(
                        source_i.force, connector_i);
                const double support_i =
                    actn2_load_kinetics::allocated_source_support(
                        raw_support_i, demand_i[source_i_index],
                        actn2_load_kinetics::norm(source_i.force));
                for (size_t source_j_index = 0;
                     source_j_index < sources_j.size(); ++source_j_index) {
                    const auto& source_j = sources_j[source_j_index];
                    // The same thick filament pulling both thin filaments is
                    // not an M-band registration signal.
                    if (!actn2_load_kinetics::distinct_motor_sources(
                            source_i.myosin, source_j.myosin)) {
                        continue;
                    }
                    ++motor_candidate_thread_pairs;
                    const double raw_support_j =
                        actn2_load_kinetics::opposing_support(
                            source_j.force, connector_j);
                    const double support_j =
                        actn2_load_kinetics::allocated_source_support(
                            raw_support_j, demand_j[source_j_index],
                            actn2_load_kinetics::norm(source_j.force));
                    const double support =
                        actn2_load_kinetics::motor_transmitted_load(
                            arm_loads, support_i, support_j);
                    const bool is_current =
                        source_i.myosin == link.candidate_myosin_i &&
                        source_j.myosin == link.candidate_myosin_j &&
                        source_i.half == link.candidate_half_i &&
                        source_j.half == link.candidate_half_j;
                    if (is_current) {
                        current_support = support;
                        current_raw_i = raw_support_i;
                        current_raw_j = raw_support_j;
                        current_allocated_i = support_i;
                        current_allocated_j = support_j;
                    }
                    if (support > best_support ||
                        (support == best_support && support > 0.0 &&
                         (best_myosin_i < 0 || identity_less(
                              source_i.myosin, source_i.half,
                              source_j.myosin, source_j.half,
                              best_myosin_i, best_half_i,
                              best_myosin_j, best_half_j)))) {
                        best_support = support;
                        best_myosin_i = source_i.myosin;
                        best_myosin_j = source_j.myosin;
                        best_half_i = source_i.half;
                        best_half_j = source_j.half;
                        best_raw_i = raw_support_i;
                        best_raw_j = raw_support_j;
                        best_allocated_i = support_i;
                        best_allocated_j = support_j;
                    }
                }
            }

            actn2_load_kinetics::MotorPairState candidate_state;
            candidate_state.identity = old_candidate_identity;
            candidate_state.since_step = link.candidate_since;
            candidate_state.last_seen_step = link.candidate_last_seen;
            candidate_state.support_ema = link.candidate_support_ema;
            candidate_state.support_instantaneous =
                link.candidate_support_instantaneous;
            candidate_state.challenger_identity = {
                link.challenger_myosin_i, link.challenger_myosin_j,
                link.challenger_half_i, link.challenger_half_j};
            candidate_state.challenger_since_step = link.challenger_since;
            candidate_state.challenger_support_ema =
                link.challenger_support_ema;
            candidate_state.locked = link.candidate_locked;
            const actn2_load_kinetics::MotorPairIdentity best_identity{
                best_myosin_i, best_myosin_j, best_half_i, best_half_j};
            const auto candidate_update =
                actn2_load_kinetics::update_motor_pair_candidate(
                    candidate_state, best_identity, best_support,
                    current_support, static_cast<std::uint64_t>(current_step),
                    dt, aa_motor_pair_support_window,
                    aa_motor_pair_commit_time, aa_motor_pair_switch_ratio,
                    aa_motor_pair_exchange_time, aa_motor_pair_dropout_time);
            link.candidate_myosin_i = candidate_state.identity.myosin_i;
            link.candidate_myosin_j = candidate_state.identity.myosin_j;
            link.candidate_half_i = candidate_state.identity.half_i;
            link.candidate_half_j = candidate_state.identity.half_j;
            link.candidate_since = static_cast<RecoveryStep>(std::min(
                candidate_state.since_step,
                static_cast<std::uint64_t>(
                    std::numeric_limits<RecoveryStep>::max())));
            link.candidate_last_seen = static_cast<RecoveryStep>(std::min(
                candidate_state.last_seen_step,
                static_cast<std::uint64_t>(
                    std::numeric_limits<RecoveryStep>::max())));
            link.candidate_support_ema = candidate_state.support_ema;
            link.candidate_support_instantaneous =
                candidate_state.support_instantaneous;
            link.challenger_myosin_i =
                candidate_state.challenger_identity.myosin_i;
            link.challenger_myosin_j =
                candidate_state.challenger_identity.myosin_j;
            link.challenger_half_i =
                candidate_state.challenger_identity.half_i;
            link.challenger_half_j =
                candidate_state.challenger_identity.half_j;
            link.challenger_since = static_cast<RecoveryStep>(std::min(
                candidate_state.challenger_since_step,
                static_cast<std::uint64_t>(
                    std::numeric_limits<RecoveryStep>::max())));
            link.challenger_support_ema =
                candidate_state.challenger_support_ema;
            motor_pair_event = candidate_update.event;
            motor_pair_candidate_age = candidate_update.candidate_age;
            motor_pair_committed = candidate_update.committed;
            motor_pair_best_support = best_support;
            motor_pair_current_support = current_support;
            if (motor_pair_event ==
                    actn2_load_kinetics::MotorPairEvent::selected ||
                motor_pair_event ==
                    actn2_load_kinetics::MotorPairEvent::switched) {
                link.candidate_raw_support_i = best_raw_i;
                link.candidate_raw_support_j = best_raw_j;
                link.candidate_allocated_support_i = best_allocated_i;
                link.candidate_allocated_support_j = best_allocated_j;
            } else if (candidate_state.identity.valid() &&
                       !candidate_state.locked) {
                link.candidate_raw_support_i = current_raw_i;
                link.candidate_raw_support_j = current_raw_j;
                link.candidate_allocated_support_i = current_allocated_i;
                link.candidate_allocated_support_j = current_allocated_j;
            } else if (!candidate_state.identity.valid()) {
                link.candidate_raw_support_i = 0.0;
                link.candidate_raw_support_j = 0.0;
                link.candidate_allocated_support_i = 0.0;
                link.candidate_allocated_support_j = 0.0;
            }
            motor_transmitted_support = link.candidate_support_ema;
            motor_candidate_thread_seconds += omp_get_wtime() - started;
        }
        double allocated_motor_i = 0.0;
        double allocated_motor_j = 0.0;
        auto kinetic_arm_loads = arm_loads;
        if (aa_status_mode ==
            actn2_load_kinetics::StatusMode::arm_load_productive) {
            allocated_motor_i = actn2_load_kinetics::allocated_motor_load(
                arm_loads.first,
                total_raw_favorable_arm_load[static_cast<size_t>(link.actin_i)],
                productive_motor_load[static_cast<size_t>(link.actin_i)]);
            allocated_motor_j = actn2_load_kinetics::allocated_motor_load(
                arm_loads.second,
                total_raw_favorable_arm_load[static_cast<size_t>(link.actin_j)],
                productive_motor_load[static_cast<size_t>(link.actin_j)]);
            kinetic_arm_loads = {
                std::min(std::max(0.0, arm_loads.first), allocated_motor_i),
                std::min(std::max(0.0, arm_loads.second), allocated_motor_j)};
        }
        // Identical to _aa_productive_opposed_load(i, j) with the two
        // productive force vectors taken from the phase-0 cache.
        const double instantaneous_support =
            actn2_load_kinetics::productive_opposed_load(
                std::array<double, 3>{actin.direction[link.actin_i].x,
                                      actin.direction[link.actin_i].y,
                                      actin.direction[link.actin_i].z},
                std::array<double, 3>{actin.direction[link.actin_j].x,
                                      actin.direction[link.actin_j].y,
                                      actin.direction[link.actin_j].z},
                aa_scratch_productive_force[
                    static_cast<size_t>(link.actin_i)],
                aa_scratch_productive_force[
                    static_cast<size_t>(link.actin_j)]);
        // Status-2 formation previously required registry AND bilateral opposed
        // motor load in the SAME timestep.  With ~6 ms strong-head lifetimes that
        // coincidence is rare, so links stalled at status 1 indefinitely.  Average
        // the productive support over a motor timescale instead: a registry-
        // qualified pair may mature on load accumulated across a short window.
        // The catch-bond load still uses the instantaneous value.  The EMA may
        // additionally suppress the opt-in mature-to-weak relaxation below,
        // allowing it to bridge brief head-turnover gaps without storing load
        // indefinitely.
        if (aa_status2_support_window > 0.0) {
            const double alpha = std::min(
                1.0, dt / std::max(aa_status2_support_window, 1e-12));
            link.support_ema += alpha * (instantaneous_support - link.support_ema);
        } else {
            link.support_ema = instantaneous_support;
        }
        const double productive_support =
            (link.status == 2) ? instantaneous_support
                               : std::max(instantaneous_support, link.support_ema);
        double effective_load = std::min(link_tension, weight * opposed);
        if (aa_status_mode ==
                actn2_load_kinetics::StatusMode::productive_opposed ||
            aa_status_mode ==
                actn2_load_kinetics::StatusMode::productive_register) {
            effective_load = productive_support;
        } else if (aa_status_mode ==
                   actn2_load_kinetics::StatusMode::arm_load_productive) {
            effective_load = kinetic_arm_loads.favorable_bilateral();
        }
        const BondFlag old_status = link.status;
        const bool registry_ok = actn2_load_kinetics::registry_gate(
            static_cast<int>(old_status), mech.direction_dot,
            mech.registry_error, aa_status2_registry_on,
            aa_status2_registry_off);
        const bool registry_retention_ok = actn2_load_kinetics::registry_gate(
            2, mech.direction_dot, mech.registry_error,
            aa_status2_registry_on, aa_status2_registry_off);
        bool instantaneous_status2_gate = false;
        if (aa_status_mode ==
            actn2_load_kinetics::StatusMode::productive_register) {
            const int classified =
                actn2_load_kinetics::classify_productive_register(
                    static_cast<int>(old_status), productive_support,
                    aa_status2_force_on, mech.direction_dot,
                    mech.registry_error, aa_status2_registry_on,
                    aa_status2_registry_off);
            instantaneous_status2_gate = old_status != 2 && classified == 2;
            link.status = static_cast<BondFlag>(classified);
            if (link.status != 2) effective_load = 0.0;
        } else if (motor_transmitted_mode) {
            const int classified =
                actn2_load_kinetics::classify_motor_transmitted_register(
                    static_cast<int>(old_status),
                    motor_transmitted_support, motor_pair_committed,
                    aa_status2_force_on, mech.direction_dot,
                    mech.registry_error, aa_status2_registry_on,
                    aa_status2_registry_off);
            instantaneous_status2_gate = old_status != 2 && classified == 2;
            link.status = static_cast<BondFlag>(classified);
            if (link.status == 2) {
                link.candidate_locked = true;
                link.challenger_myosin_i = -1;
                link.challenger_myosin_j = -1;
                link.challenger_half_i = 0;
                link.challenger_half_j = 0;
                link.challenger_since = 0;
                link.challenger_support_ema = 0.0;
            }
            effective_load = link.status == 2
                ? arm_loads.favorable_bilateral() : 0.0;
        } else if (aa_status_mode ==
                   actn2_load_kinetics::StatusMode::arm_load_productive) {
            const int classified =
                actn2_load_kinetics::classify_arm_loaded_productive(
                    static_cast<int>(link.status), kinetic_arm_loads,
                    productive_support, aa_status2_force_on,
                    aa_status2_force_off);
            instantaneous_status2_gate = classified == 2 && registry_ok;
            link.status = static_cast<BondFlag>(
                classified == 2 && registry_ok ? 2 : 1);
            if (link.status != 2) effective_load = 0.0;
        } else {
            const double threshold = link.status == 2
                ? aa_status2_force_off : aa_status2_force_on;
            link.status = effective_load > threshold && registry_ok ? 2 : 1;
        }
        bool status2_unlocked = false;
        if (old_status == 2 && aa_status2_unlock_rate > 0.0) {
            // Status 2 is a reversible ACTN2 conformation.  It survives brief
            // head turnover through the existing motor-timescale filter, but
            // sustained loss of source-resolved support allows a stochastic
            // return to the weak, annealable state without detaching ACTN2.
            const double filtered_support = motor_transmitted_mode
                ? link.candidate_support_ema : link.support_ema;
            const double unlock_rate = actn2_load_kinetics::mature_unlock_rate(
                aa_status2_unlock_rate, filtered_support,
                aa_status2_unlock_force_scale, mech.registry_error,
                aa_status2_registry_off, aa_status2_unlock_registry_scale);
            const double unlock_probability =
                actn2_load_kinetics::mature_unlock_probability(unlock_rate, dt);
            counter_rng::CounterRng unlock_rng(
                static_cast<std::uint64_t>(initial_seed),
                static_cast<std::uint64_t>(current_step),
                counter_rng::Stream::aa_status2_unlock, link.id);
            status2_unlocked = unlock_rng.uniform() < unlock_probability;
            if (status2_unlocked) {
                link.status = 1;
                link.status2_memory_until = 0;
                link.candidate_locked = false;
                effective_load = 0.0;
            } else {
                link.status = 2;
                if (aa_status_mode ==
                        actn2_load_kinetics::StatusMode::productive_opposed ||
                    aa_status_mode ==
                        actn2_load_kinetics::StatusMode::productive_register) {
                    effective_load = instantaneous_support;
                } else if (motor_transmitted_mode) {
                    effective_load = arm_loads.favorable_bilateral();
                } else if (aa_status_mode ==
                           actn2_load_kinetics::StatusMode::arm_load_productive) {
                    effective_load = kinetic_arm_loads.favorable_bilateral();
                }
            }
        }
        if (aa_status2_unlock_rate <= 0.0 &&
            (aa_status2_locked_until_detach || aa_arm_state_model) &&
            old_status == 2) {
            // Stereospecific locking is a molecular state, not an
            // instantaneous force threshold.  Head turnover may remove the
            // load that created the lock, but it cannot make the binding
            // domains diffuse again without an actual ACTN2 unbinding event.
            // The status-2 lifetime below still responds to the current load,
            // and contour rupture remains an immediate escape path.
            link.status = 2;
        }
        if (aa_arm_state_model) {
            if (link.status == 2) {
                if (old_status != 2) {
                    link.arm_state_since = static_cast<RecoveryStep>(std::min(
                        current_step,
                        static_cast<size_t>(
                            std::numeric_limits<RecoveryStep>::max())));
                }
                link.arm_state = actn2_arm_kinetics::State::mature_bridge;
            } else {
                if (old_status == 2) {
                    link.arm_state_since = static_cast<RecoveryStep>(std::min(
                        current_step,
                        static_cast<size_t>(
                            std::numeric_limits<RecoveryStep>::max())));
                }
                link.arm_state = actn2_arm_kinetics::State::weak_bridge;
            }
        }
        const bool register_retention_mode =
            aa_status_mode ==
                actn2_load_kinetics::StatusMode::productive_register ||
            motor_transmitted_mode;
        if (!register_retention_mode &&
            link.status == 2 && am_status2_memory_steps > 0) {
            link.status2_memory_until = static_cast<RecoveryStep>(std::min(
                static_cast<size_t>(std::numeric_limits<RecoveryStep>::max()),
                current_step + am_status2_memory_steps));
        } else if (!register_retention_mode &&
                   link.status == 1 && registry_retention_ok &&
                   static_cast<size_t>(link.status2_memory_until) >= current_step) {
            if (aa_status_mode !=
                    actn2_load_kinetics::StatusMode::arm_load_productive ||
                actn2_load_kinetics::classify_existing_bond(
                    2, kinetic_arm_loads, aa_status2_force_on,
                    aa_status2_force_off) == 2) {
                // In the arm-resolved mode memory may bridge a brief loss of
                // strong-head support only while the physical link continues
                // to carry bilateral load. It cannot stabilize an unloaded
                // link or create status 2 from status 1.
                link.status = 2;
                effective_load = kinetic_arm_loads.favorable_bilateral();
            }
        }
        mechanics[slot] = compute_local_aa_link(
            actin, link.actin_i, link.actin_j, link.s_i, link.s_j, box,
            ramp * aa_link_radial_stiffness(link.status),
            ramp * aa_link_perpendicular_stiffness(link.status),
            single_link_alignment_k, aa_optimal, aa_antiparallel_overlap,
            aa_link_register_ratio(link.status),
            aa_register_saturation_error, aa_axial_register_force_max,
            ramp * aa_radial_force_max);
        const double lifetime = catch_bond_kinetics::lifetime_for_status(
            link.status, base_lifetime, directional_base_lifetime,
            lifetime_coeff, effective_load, stabilize_status1);
        counter_rng::CounterRng detachment_rng(
            static_cast<std::uint64_t>(initial_seed),
            static_cast<std::uint64_t>(current_step),
            counter_rng::Stream::aa_link_detachment,
            static_cast<std::uint64_t>(link.id));
        const double registry_threshold = old_status == 2
            ? aa_status2_registry_off : aa_status2_registry_on;
        const double arm_load_threshold = old_status == 2
            ? aa_status2_force_off : aa_status2_force_on;
        const bool arm_load_ok =
            kinetic_arm_loads.first > arm_load_threshold &&
            kinetic_arm_loads.second > arm_load_threshold;
        const bool productive_support_ok = productive_support > 0.0;
        // Association already rejects ACTN2 sites inside AM overlap. Apply the
        // corresponding retention rule as well: if either actin has no
        // crosslinkable barbed-end contour left, the existing molecular link is
        // geometrically occluded and ruptures immediately, including status 2.
        // Short-circuiting here deliberately avoids consuming a stochastic draw.
        catch_bond_kinetics::DetachmentDecision detachment;
        if (!aa_arm_state_model) {
            detachment = catch_bond_kinetics::decide_am_overlap_occlusion(
                actin_crosslink_ratio[link.actin_i],
                actin_crosslink_ratio[link.actin_j], EPS);
        }
        actn2_arm_kinetics::OpeningDecision arm_opening;
        const auto old_arm_state = link.arm_state;
        const double tensile_arm_i = link_tension + contact_arm_loads[slot][0];
        const double tensile_arm_j = link_tension + contact_arm_loads[slot][1];
        link.arm_force_i = tensile_arm_i;
        link.arm_force_j = tensile_arm_j;
        link.contact_force = contact_arm_loads[slot][0] +
                             contact_arm_loads[slot][1];
        link.contact_penetration = contact_penetration[slot];
        if (!detachment.detaches() &&
            aa_occlusion_rupture_penetration > 0.0 &&
            contact_penetration[slot] >= aa_occlusion_rupture_penetration) {
            detachment.actn2_myosin_collision = true;
            detachment.probability = 1.0;
        }
        if (aa_arm_state_model && !detachment.detaches() &&
            aa_link_rupture_length > 0.0 &&
            mechanics[slot].length > aa_link_rupture_length) {
            detachment.geometric = true;
            detachment.probability = 1.0;
        }
        const bool attachment_occluded_i = aa_arm_state_model &&
            actn2_arm_kinetics::contour_fully_occluded(
                actin_crosslink_ratio[link.actin_i], EPS);
        const bool attachment_occluded_j = aa_arm_state_model &&
            actn2_arm_kinetics::contour_fully_occluded(
                actin_crosslink_ratio[link.actin_j], EPS);
        if (aa_arm_state_model && !detachment.detaches() &&
            (attachment_occluded_i || attachment_occluded_j)) {
            // A scalar partial-overlap boundary is not itself a molecular
            // collision.  Partial overlap is handled by the explicit local
            // ACTN2-myosin contact above, whose force enters the Bell arm
            // hazard.  Deterministic opening remains only when an actin has
            // no crosslinkable contour at all, as required by the geometric
            // model.  If both contours vanish, deterministic link-id ordering
            // chooses i and the remaining tether is removed next step.
            arm_opening.opens_i = attachment_occluded_i;
            arm_opening.opens_j = !attachment_occluded_i && attachment_occluded_j;
            arm_opening.rate_i = attachment_occluded_i
                ? std::numeric_limits<double>::infinity() : 0.0;
            arm_opening.rate_j = attachment_occluded_j
                ? std::numeric_limits<double>::infinity() : 0.0;
            arm_opening.probability = 1.0;
        } else if (aa_arm_state_model && !detachment.detaches()) {
            const double zero_force_rate = link.status == 2
                ? aa_arm_open_rate_mature : aa_arm_open_rate_weak;
            const double rate_i = actn2_arm_kinetics::opening_rate(
                zero_force_rate, tensile_arm_i, aa_arm_bell_distance,
                thermal_beta);
            const double rate_j = actn2_arm_kinetics::opening_rate(
                zero_force_rate, tensile_arm_j, aa_arm_bell_distance,
                thermal_beta);
            counter_rng::CounterRng opening_rng(
                static_cast<std::uint64_t>(initial_seed),
                static_cast<std::uint64_t>(current_step),
                counter_rng::Stream::aa_arm_opening, link.id);
            arm_opening = actn2_arm_kinetics::decide_bridge_opening(
                rate_i, rate_j, dt, [&]() { return opening_rng.uniform(); });
        } else if (!detachment.detaches()) {
            detachment = catch_bond_kinetics::decide_detachment(
                mechanics[slot].length, aa_link_rupture_length, dt, lifetime,
                [&]() { return detachment_rng.uniform(); });
        }
        const bool contour_rupture = detachment.geometric;
        const bool overlap_occlusion = detachment.am_overlap_occluded;
        const bool collision_rupture = detachment.actn2_myosin_collision;
        const bool geometric_rupture =
            contour_rupture || overlap_occlusion || collision_rupture;
        const bool arm_opens = arm_opening.opens_i || arm_opening.opens_j;
        const bool detaches = detachment.detaches();
        // Everything from here mutates shared, order-sensitive state: the
        // event streams are appended in slot order, aa_neighbor_slots is a
        // single hash map, and the free-slot stack decides which slot the next
        // association reuses.  `ordered` replays this block in the original
        // sequential order regardless of how the iterations were scheduled.
        {
        auto record_motor_pair_event = [&actn2_load_kinetics_unused =
                                             motor_transmitted_mode,
                                         this, &link,
                                         motor_pair_current_support,
                                         motor_pair_best_support,
                                         motor_pair_candidate_age](
            int event_code,
            const actn2_load_kinetics::MotorPairIdentity& old_identity,
            const actn2_load_kinetics::MotorPairIdentity& new_identity) {
            if (!actn2_load_kinetics_unused) return;
            aa_motor_pair_events.insert(aa_motor_pair_events.end(), {
                static_cast<double>(link.id >> 32),
                static_cast<double>(link.id & 0xffffffffULL),
                static_cast<double>(current_step),
                static_cast<double>(event_code),
                static_cast<double>(old_identity.myosin_i),
                static_cast<double>(old_identity.myosin_j),
                static_cast<double>(old_identity.half_i),
                static_cast<double>(old_identity.half_j),
                static_cast<double>(new_identity.myosin_i),
                static_cast<double>(new_identity.myosin_j),
                static_cast<double>(new_identity.half_i),
                static_cast<double>(new_identity.half_j),
                motor_pair_current_support, motor_pair_best_support,
                link.candidate_support_ema, motor_pair_candidate_age,
                link.candidate_raw_support_i,
                link.candidate_raw_support_j,
                link.candidate_allocated_support_i,
                link.candidate_allocated_support_j});
        };
        const actn2_load_kinetics::MotorPairIdentity new_candidate_identity{
            link.candidate_myosin_i, link.candidate_myosin_j,
            link.candidate_half_i, link.candidate_half_j};
        if (motor_pair_event !=
            actn2_load_kinetics::MotorPairEvent::none) {
            record_motor_pair_event(
                static_cast<int>(motor_pair_event), old_candidate_identity,
                new_candidate_identity);
        }
        if (!candidate_was_locked && link.candidate_locked) {
            record_motor_pair_event(4, new_candidate_identity,
                                    new_candidate_identity);
        }
        if (status2_unlocked && candidate_was_locked &&
            new_candidate_identity.valid()) {
            record_motor_pair_event(6, new_candidate_identity,
                                    new_candidate_identity);
        }
        if (link.status != old_status) {
            aa_local_link_events.insert(aa_local_link_events.end(), {
                static_cast<double>(link.id),
                static_cast<double>(link.actin_i),
                static_cast<double>(link.actin_j),
                static_cast<double>(current_step), 3.0,
                static_cast<double>(old_status),
                static_cast<double>(link.status), effective_load, lifetime,
                link.s_i, link.s_j});
            aa_local_link_events_v2.insert(aa_local_link_events_v2.end(), {
                static_cast<double>(link.id), static_cast<double>(link.actin_i),
                static_cast<double>(link.actin_j), static_cast<double>(current_step),
                3.0, static_cast<double>(old_status),
                static_cast<double>(link.status), effective_load, lifetime,
                link.s_i, link.s_j, kinetic_arm_loads.first,
                kinetic_arm_loads.second,
                kinetic_arm_loads.favorable_bilateral(), productive_support,
                instantaneous_status2_gate ? 1.0 : 0.0, link_tension,
                arm_loads.first, arm_loads.second,
                allocated_motor_i, allocated_motor_j});
            aa_local_link_events_v3.insert(aa_local_link_events_v3.end(), {
                static_cast<double>(link.id), static_cast<double>(link.actin_i),
                static_cast<double>(link.actin_j), static_cast<double>(current_step),
                3.0, static_cast<double>(old_status),
                static_cast<double>(link.status), effective_load, lifetime,
                link.s_i, link.s_j, kinetic_arm_loads.first,
                kinetic_arm_loads.second,
                kinetic_arm_loads.favorable_bilateral(), productive_support,
                instantaneous_status2_gate ? 1.0 : 0.0, link_tension,
                arm_loads.first, arm_loads.second,
                allocated_motor_i, allocated_motor_j,
                mechanics[slot].length, aa_link_rupture_length,
                mechanics[slot].direction_dot, mechanics[slot].registry_error,
                registry_threshold, registry_ok ? 1.0 : 0.0,
                arm_load_ok ? 1.0 : 0.0,
                productive_support_ok ? 1.0 : 0.0,
                instantaneous_status2_gate ? 1.0 : 0.0,
                0.0, detachment.probability, detachment.draw, 0.0, 0.0});
        }
        if (aa_neighbor_identity_enabled && !arm_opens) {
            aa_neighbor_slots.observe(
                link.actin_i, link.actin_j, static_cast<int>(link.status),
                static_cast<aa_neighbor_identity::Step>(current_step),
                static_cast<aa_neighbor_identity::Step>(
                    aa_status2_slot_memory_steps));
        }
        if (arm_opens && !detaches) {
            const int old_i = link.actin_i;
            const int old_j = link.actin_j;
            const int bound_actin = arm_opening.opens_i ? old_j : old_i;
            const int opened_partner = arm_opening.opens_i ? old_i : old_j;
            const double bound_s = arm_opening.opens_i ? link.s_j : link.s_i;
            const double opened_s = arm_opening.opens_i ? link.s_i : link.s_j;
            link.last_partner = opened_partner;
            link.status = 1;
            link.arm_state = arm_opening.opens_i
                ? actn2_arm_kinetics::State::tethered_to_j
                : actn2_arm_kinetics::State::tethered_to_i;
            link.arm_state_since = static_cast<RecoveryStep>(std::min(
                current_step,
                static_cast<size_t>(std::numeric_limits<RecoveryStep>::max())));
            if (new_candidate_identity.valid()) {
                record_motor_pair_event(
                    5, new_candidate_identity,
                    actn2_load_kinetics::MotorPairIdentity{});
            }
            link.candidate_locked = false;
            link.candidate_myosin_i = -1;
            link.candidate_myosin_j = -1;
            link.challenger_myosin_i = -1;
            link.challenger_myosin_j = -1;
            link.candidate_support_ema = 0.0;
            link.challenger_support_ema = 0.0;
            if (aa_neighbor_identity_enabled) {
                aa_neighbor_slots.erase(old_i, old_j);
            }
            aa_arm_events.insert(aa_arm_events.end(), {
                static_cast<double>(link.id),
                static_cast<double>(current_step), 1.0,
                static_cast<double>(static_cast<int>(old_arm_state)),
                static_cast<double>(static_cast<int>(link.arm_state)),
                static_cast<double>(bound_actin),
                static_cast<double>(opened_partner), -1.0, bound_s, opened_s,
                arm_opening.opens_i ? tensile_arm_i : tensile_arm_j,
                link.contact_force, link.contact_penetration,
                arm_opening.opens_i ? arm_opening.rate_i : arm_opening.rate_j});
        }
        if (detaches) {
            // Code 5 means that an ACTN2 detachment actually released a held
            // motor-pair identity.  Do not emit a release row for links that
            // never selected a candidate; those are ordinary ACTN2
            // detachments and otherwise inflate the partner-release count.
            if (new_candidate_identity.valid()) {
                record_motor_pair_event(
                    5, new_candidate_identity,
                    actn2_load_kinetics::MotorPairIdentity{});
            }
            if (geometric_rupture && aa_neighbor_identity_enabled) {
                aa_neighbor_slots.erase(link.actin_i, link.actin_j);
            }
            if (link.status == 2) {
                // Slot 0 always exists and, unlike omp_get_thread_num() inside
                // an `ordered` region, does not depend on which thread happened
                // to reach this iteration.  At one thread this is the index the
                // routine already used.
                aa_completed_lifetimes_temp[0].push_back(
                    (current_step - static_cast<size_t>(link.attach_step)) * dt);
            }
            aa_local_link_events.insert(aa_local_link_events.end(), {
                static_cast<double>(link.id),
                static_cast<double>(link.actin_i),
                static_cast<double>(link.actin_j),
                static_cast<double>(current_step), 2.0,
                static_cast<double>(link.status), 0.0, effective_load,
                lifetime, link.s_i, link.s_j});
            aa_local_link_events_v2.insert(aa_local_link_events_v2.end(), {
                static_cast<double>(link.id), static_cast<double>(link.actin_i),
                static_cast<double>(link.actin_j), static_cast<double>(current_step),
                2.0, static_cast<double>(link.status), 0.0, effective_load,
                lifetime, link.s_i, link.s_j, kinetic_arm_loads.first,
                kinetic_arm_loads.second,
                kinetic_arm_loads.favorable_bilateral(), productive_support,
                instantaneous_status2_gate ? 1.0 : 0.0, link_tension,
                arm_loads.first, arm_loads.second,
                allocated_motor_i, allocated_motor_j});
            aa_local_link_events_v3.insert(aa_local_link_events_v3.end(), {
                static_cast<double>(link.id), static_cast<double>(link.actin_i),
                static_cast<double>(link.actin_j), static_cast<double>(current_step),
                2.0, static_cast<double>(link.status), 0.0, effective_load,
                lifetime, link.s_i, link.s_j, kinetic_arm_loads.first,
                kinetic_arm_loads.second,
                kinetic_arm_loads.favorable_bilateral(), productive_support,
                instantaneous_status2_gate ? 1.0 : 0.0, link_tension,
                arm_loads.first, arm_loads.second,
                allocated_motor_i, allocated_motor_j,
                mechanics[slot].length, aa_link_rupture_length,
                mechanics[slot].direction_dot, mechanics[slot].registry_error,
                registry_threshold, registry_ok ? 1.0 : 0.0,
                arm_load_ok ? 1.0 : 0.0,
                productive_support_ok ? 1.0 : 0.0,
                instantaneous_status2_gate ? 1.0 : 0.0,
                static_cast<double>(static_cast<int>(detachment.reason())),
                detachment.probability, detachment.draw,
                geometric_rupture ? 1.0 : 0.0,
                detachment.stochastic ? 1.0 : 0.0});
            link.active = false;
            aa_free_link_slots.push_back(slot);
        }
        }  // formerly omp ordered
    }
    }  // omp single
    if (motor_transmitted_mode) {
        #pragma omp atomic update
        aa_motor_candidate_seconds += motor_candidate_thread_seconds;
        #pragma omp atomic update
        aa_motor_candidate_pairs_evaluated += motor_candidate_thread_pairs;
    }

    stamp_aa_phase(kAaReindex);
    #pragma omp single
    {
    _rebuild_local_aa_indices();
    if (aa_neighbor_identity_enabled) {
        aa_neighbor_slots.expire_unseen(
            static_cast<aa_neighbor_identity::Step>(current_step),
            static_cast<aa_neighbor_identity::Step>(
                aa_neighbor_slot_memory_steps));
    }
    if (aa_scratch_sorted_neighbors.size() != static_cast<size_t>(actin.n)) {
        aa_scratch_sorted_neighbors.assign(static_cast<size_t>(actin.n), {});
    }
    aa_scratch_pair_offset.assign(static_cast<size_t>(actin.n) + 1, 0);
    aa_scratch_dirty_actin.assign(static_cast<size_t>(actin.n), 0);
    }

    // Canonical proposal order: outer actin ascending, partner ascending,
    // restricted to i > j.  Materializing that list lets the (pure) evaluation
    // of every pair run across the team before the serial commit.  Building it
    // is itself O(N log k) and was the largest remaining serial block -- worse,
    // one thread rebuilding a list the whole team had just dirtied got *slower*
    // as the team grew.  It is now a per-actin sort in parallel, an O(N) prefix
    // sum, and a parallel scatter into exact offsets, which reproduces the same
    // ordering the serial loop produced.
    stamp_aa_phase(kAaSortNeighbors);
    #pragma omp for schedule(static)
    for (int i = 0; i < actin.n; ++i) {
        auto& sorted = aa_scratch_sorted_neighbors[static_cast<size_t>(i)];
        // Filter before sorting.  The surviving list is a subsequence of the
        // old one in the same ascending order, so the pair list, its ordering,
        // and therefore the greedy commit order below are all unchanged; only
        // pairs that provably cannot fire are absent.
        sorted.clear();
        // Source from the segment-distance-filtered list (~539 entries) rather
        // than the raw centre-based one (~890).  Its cutoff is 0.77 um and a
        // viable ACTN2 pair needs segment distance below aa_cutoff = 0.12 um,
        // so nothing that could associate is missing from it.
        if (opt_gates_enabled) {
            for (int j : actin_close_neighbors[static_cast<size_t>(i)]) {
                if (_aa_pair_may_associate(i, j)) sorted.push_back(j);
            }
        } else {
            sorted.assign(actin_neighbors_by_species[i].first.begin(),
                          actin_neighbors_by_species[i].first.end());
        }
        std::sort(sorted.begin(), sorted.end());
        size_t count = 0;
        for (int j : sorted) if (i > j) ++count;
        aa_scratch_pair_offset[static_cast<size_t>(i) + 1] = count;
    }

    stamp_aa_phase(kAaPrefixSum);
    #pragma omp single
    {
    for (int i = 0; i < actin.n; ++i) {
        aa_scratch_pair_offset[static_cast<size_t>(i) + 1] +=
            aa_scratch_pair_offset[static_cast<size_t>(i)];
    }
    const size_t total_pairs = aa_scratch_pair_offset[
        static_cast<size_t>(actin.n)];
    // resize(), not assign(): every entry is overwritten by the evaluation
    // loop below, so re-initializing the whole buffer each step is wasted.
    if (aa_scratch_pairs.size() < total_pairs) {
        aa_scratch_pairs.resize(total_pairs);
        aa_scratch_proposals.resize(total_pairs);
    }
    aa_scratch_pair_count = total_pairs;
    }

    stamp_aa_phase(kAaScatterPairs);
    #pragma omp for schedule(static)
    for (int i = 0; i < actin.n; ++i) {
        size_t at = aa_scratch_pair_offset[static_cast<size_t>(i)];
        for (int j : aa_scratch_sorted_neighbors[static_cast<size_t>(i)]) {
            if (i <= j) continue;
            aa_scratch_pairs[at++] = {i, j};
        }
    }

    // Phase 3a.  Evaluate every candidate against the post-detachment state.
    // _evaluate_aa_association() is const, so the whole team can run it.
    const long long n_pairs = static_cast<long long>(aa_scratch_pair_count);
    stamp_aa_phase(kAaEvalProposals);
    #pragma omp for schedule(static)
    for (long long k = 0; k < n_pairs; ++k) {
        aa_scratch_proposals[static_cast<size_t>(k)] = _evaluate_aa_association(
            aa_scratch_pairs[static_cast<size_t>(k)].first,
            aa_scratch_pairs[static_cast<size_t>(k)].second);
    }

    stamp_aa_phase(kAaCommitLinks);
    #pragma omp single
    {
    auto add_link = [&](int i, int j, const AAAssociationProposal& proposal) {
        const double s_i = proposal.s_i;
        const double s_j = proposal.s_j;
        if (aa_zdisc_tip_contact && proposal.nearest_tip_myosin >= 0 &&
            proposal.nearest_tip_gap < aa_zdisc_tip_onset - 1e-12) {
            throw std::runtime_error(
                "ACTN2 association attempted to commit inside a myosin-tip "
                "exclusion zone");
        }
        AACrosslink link;
        link.id = aa_next_link_id++;
        link.actin_i = i;
        link.actin_j = j;
        link.s_i = s_i;
        link.s_j = s_j;
        link.status = 1;
        link.arm_state = actn2_arm_kinetics::State::weak_bridge;
        link.attach_step = static_cast<RecoveryStep>(std::min(
            current_step,
            static_cast<size_t>(std::numeric_limits<RecoveryStep>::max())));
        link.arm_state_since = link.attach_step;
        link.active = true;
        aa_tip_accessibility_events.insert(
            aa_tip_accessibility_events.end(), {
                static_cast<double>(link.id),
                static_cast<double>(current_step), 1.0,
                static_cast<double>(i), static_cast<double>(j), s_i, s_j,
                static_cast<double>(proposal.nearest_tip_myosin),
                static_cast<double>(proposal.nearest_tip_sign),
                proposal.nearest_tip_gap, proposal.nearest_tip_rho, 0.0,
                proposal.nearest_tip_motor_engaged ? 1.0 : 0.0});
        aa_local_link_events.insert(aa_local_link_events.end(), {
            static_cast<double>(link.id), static_cast<double>(i),
            static_cast<double>(j), static_cast<double>(current_step),
            1.0, 0.0, 1.0, 0.0, base_lifetime, s_i, s_j});
        aa_local_link_events_v2.insert(aa_local_link_events_v2.end(), {
            static_cast<double>(link.id), static_cast<double>(i),
            static_cast<double>(j), static_cast<double>(current_step),
            1.0, 0.0, 1.0, 0.0, base_lifetime, s_i, s_j,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0});
        size_t slot = 0;
        if (!aa_free_link_slots.empty()) {
            slot = aa_free_link_slots.back();
            aa_free_link_slots.pop_back();
            aa_links[slot] = link;
        } else {
            slot = aa_links.size();
            aa_links.push_back(link);
            mechanics.emplace_back();
            contact_results.emplace_back();
            contact_evaluated.push_back(0);
        }
        contact_results[slot].clear();
        contact_evaluated[slot] = 0;
        mechanics[slot] = compute_local_aa_link(
            actin, i, j, s_i, s_j, box,
            ramp * aa_link_radial_stiffness(link.status),
            ramp * aa_link_perpendicular_stiffness(link.status),
            single_link_alignment_k, aa_optimal, aa_antiparallel_overlap,
            aa_link_register_ratio(link.status),
            aa_register_saturation_error, aa_axial_register_force_max,
            ramp * aa_radial_force_max);
        const auto& initial_mechanics = mechanics[slot];
        const bool initial_registry_ok = actn2_load_kinetics::registry_gate(
            1, initial_mechanics.direction_dot, initial_mechanics.registry_error,
            aa_status2_registry_on, aa_status2_registry_off);
        const double diagnostic_nan =
            std::numeric_limits<double>::quiet_NaN();
        aa_local_link_events_v3.insert(aa_local_link_events_v3.end(), {
            static_cast<double>(link.id), static_cast<double>(i),
            static_cast<double>(j), static_cast<double>(current_step),
            1.0, 0.0, 1.0, 0.0, base_lifetime, s_i, s_j,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            initial_mechanics.length, aa_link_rupture_length,
            initial_mechanics.direction_dot, initial_mechanics.registry_error,
            aa_status2_registry_on, initial_registry_ok ? 1.0 : 0.0,
            0.0, 0.0, 0.0,
            0.0, diagnostic_nan, diagnostic_nan, 0.0, 0.0});
        // Maintain sparse indices incrementally. Rebuilding every index after
        // one local KMC association made a successful event O(Nlinks).
        aa_active_link_slots.push_back(slot);
        aa_links_by_pair[_aa_pair_key(i, j)].push_back(slot);
        aa_links_by_actin[static_cast<size_t>(i)].push_back(slot);
        aa_links_by_actin[static_cast<size_t>(j)].push_back(slot);
    };

    // Association proposals are committed in canonical actin order, exactly as
    // the original single serial loop did.  A proposal only depends on state
    // that add_link()/reserve() touch at the two actins of the committed pair
    // (pair occupancy, occupied axial sites, neighbor degree), so a cached
    // verdict is valid unless one of its actins has already been committed
    // this step; those few pairs are re-evaluated here.  The global link pool,
    // when enabled, makes every rate depend on the running link count, so that
    // configuration re-evaluates unconditionally.
    const bool pool_limited = aa_total_link_pool > 0;
    for (size_t k = 0; k < aa_scratch_pair_count; ++k) {
        const int i = aa_scratch_pairs[k].first;
        const int j = aa_scratch_pairs[k].second;
        const bool stale = pool_limited ||
            aa_scratch_dirty_actin[static_cast<size_t>(i)] ||
            aa_scratch_dirty_actin[static_cast<size_t>(j)];
        const AAAssociationProposal proposal = stale
            ? _evaluate_aa_association(i, j)
            : aa_scratch_proposals[k];
        if (!proposal.fire) continue;
        if (aa_neighbor_identity_enabled &&
            !aa_neighbor_slots.reserve(
                i, j,
                static_cast<aa_neighbor_identity::Step>(current_step),
                max_aa_bonds_per_actin)) {
            continue;
        }
        add_link(i, j, proposal);
        aa_scratch_dirty_actin[static_cast<size_t>(i)] = 1;
        aa_scratch_dirty_actin[static_cast<size_t>(j)] = 1;
    }
    }  // omp single

    // Apply mechanics and create the old pair summaries as derived state.
    // This phase does no filament mechanics; it is a scatter of already
    // computed quantities into shared per-actin accumulators and the dense
    // pair matrices, so it stays serialized rather than paying for privatized
    // copies of those matrices.
    stamp_aa_phase(kAaTipContact);
    // ---- parallel tip-contact scan -------------------------------------
    // Read-only: patch geometry, the conservative sphere test and the contact
    // kernel.  Survivors go to per-link buckets; the order-sensitive reduction
    // runs serially below over the same link sequence.
    if (aa_zdisc_tip_contact) {
        const size_t n_slots = aa_active_link_slots.size();
        #pragma omp single
        {
        if (aa_tip_scan_buckets.size() < n_slots) {
            aa_tip_scan_buckets.resize(n_slots);
        }
        const size_t threads = static_cast<size_t>(omp_get_max_threads());
        if (aa_tip_stamp_pool.size() < threads) {
            aa_tip_stamp_pool.resize(threads);
        }
        for (auto& stamps : aa_tip_stamp_pool) {
            if (stamps.size() != static_cast<size_t>(myosin.n)) {
                stamps.assign(static_cast<size_t>(myosin.n), -1);
            }
        }
        }
        #pragma omp for schedule(dynamic, 8)
        for (long long index = 0; index < static_cast<long long>(n_slots);
             ++index) {
            auto& bucket = aa_tip_scan_buckets[static_cast<size_t>(index)];
            bucket.clear();
            const size_t slot = aa_active_link_slots[static_cast<size_t>(index)];
            const auto& link = aa_links[slot];
            if (aa_arm_state_model &&
                !actn2_arm_kinetics::is_bridge(link.arm_state)) continue;
            auto& stamps = aa_tip_stamp_pool[
                static_cast<size_t>(omp_get_thread_num())];
            // The link index is globally unique and each thread owns its own
            // stamp array, so it doubles as the visit token with no counter.
            const long long token = index;
            const auto patch = compute_actn2_patch_geometry(
                actin, link.actin_i, link.actin_j, link.s_i, link.s_j, box);
            for (int pass = 0; pass < 2; ++pass) {
                const auto& candidates = actin_neighbors_by_species[
                    static_cast<size_t>(pass == 0 ? link.actin_i
                                                  : link.actin_j)].second;
                for (int m : candidates) {
                    if (m < 0 || m >= myosin.n) continue;
                    if (stamps[static_cast<size_t>(m)] == token) continue;
                    stamps[static_cast<size_t>(m)] = token;
                    if (!actn2_patch_may_contact_myosin_tip(
                            patch, myosin, m, box, aa_zdisc_tip_onset,
                            aa_zdisc_tip_transverse_cutoff)) continue;
                    for (int end_sign : {-1, 1}) {
                        auto contact = compute_zdisc_myosin_tip_contact(
                            actin, myosin, link.actin_i, link.actin_j,
                            link.s_i, link.s_j, m, end_sign, box,
                            aa_zdisc_tip_onset, aa_zdisc_tip_stiffness,
                            aa_zdisc_tip_force_max,
                            aa_zdisc_tip_transverse_flat,
                            aa_zdisc_tip_transverse_cutoff);
                        if (!contact.applied) continue;
                        bucket.push_back({slot, m, end_sign, contact});
                    }
                }
            }
        }
    }
    #pragma omp single
    {
    const int thread_id = omp_get_thread_num();
    auto& local_forces = actin_forces_temp[thread_id];
    auto& local_torques = actin_torques_temp[thread_id];
    auto& local_myosin_forces = myosin_forces_temp[thread_id];
    auto& local_myosin_torques = myosin_torques_temp[thread_id];
    auto& local_status = actin_cb_status_temp[thread_id];
    auto& aa_local_link_count_feature = actin["aa_local_link_count"];
    auto& aa_local_status1_count_feature = actin["aa_local_status1_count"];
    auto& aa_local_status2_count_feature = actin["aa_local_status2_count"];
    auto& aa_local_tension_max_feature = actin["aa_local_tension_max"];
    auto& aa_zdisc_tip_force_feature = actin["aa_zdisc_tip_force"];
    auto& aa_zdisc_tip_contact_count_feature =
        actin["aa_zdisc_tip_contact_count"];
    if (aa_zdisc_tip_contact) {
        struct SelectedTipContact {
            size_t slot = std::numeric_limits<size_t>::max();
            int myosin_index = -1;
            ZDiscMyosinTipContactResult contact;
        };
        std::vector<SelectedTipContact> selected(
            static_cast<size_t>(2 * myosin.n));
        // Replay the scan results in canonical link order.  Identical
        // comparison sequence to the old serial loop, so the winner is the
        // same even where the epsilon band makes the rule non-transitive.
        for (size_t index = 0; index < aa_active_link_slots.size(); ++index) {
            for (const auto& hit : aa_tip_scan_buckets[index]) {
                const auto& link = aa_links[hit.slot];
                auto& best = selected[static_cast<size_t>(
                    2 * hit.myosin_index + (hit.end_sign > 0 ? 1 : 0))];
                if (best.myosin_index < 0 ||
                    hit.contact.force > best.contact.force + 1e-12 ||
                    (std::abs(hit.contact.force - best.contact.force) <= 1e-12 &&
                     link.id < aa_links[best.slot].id)) {
                    best = {hit.slot, hit.myosin_index, hit.contact};
                }
            }
        }
        // One selected patch per tip provides one finite force budget even
        // when several ACTN2 molecules occupy the same local Z-disc patch.
        for (size_t tip_index = 0; tip_index < selected.size(); ++tip_index) {
            const auto& best = selected[tip_index];
            if (best.myosin_index < 0) continue;
            const auto& link = aa_links[best.slot];
            const auto& contact = best.contact;
            local_forces[link.actin_i] += contact.force_on_actin_i;
            local_forces[link.actin_j] += contact.force_on_actin_j;
            local_myosin_forces[best.myosin_index] += contact.force_on_myosin;
            local_torques[link.actin_i] += contact.torque_on_actin_i;
            local_torques[link.actin_j] += contact.torque_on_actin_j;
            local_myosin_torques[best.myosin_index] += contact.torque_on_myosin;
            aa_zdisc_tip_force_feature[link.actin_i] += 0.5 * contact.force;
            aa_zdisc_tip_force_feature[link.actin_j] += 0.5 * contact.force;
            aa_zdisc_tip_contact_count_feature[link.actin_i] += 1.0;
            aa_zdisc_tip_contact_count_feature[link.actin_j] += 1.0;
            const bool plus = (tip_index % 2) == 1;
            myosin[plus ? "aa_zdisc_tip_force_plus" :
                           "aa_zdisc_tip_force_minus"][best.myosin_index] =
                contact.force;
            myosin[plus ? "aa_zdisc_tip_gap_plus" :
                           "aa_zdisc_tip_gap_minus"][best.myosin_index] =
                contact.axial_gap;
            myosin[plus ? "aa_zdisc_tip_rho_plus" :
                           "aa_zdisc_tip_rho_minus"][best.myosin_index] =
                contact.transverse_distance;
            myosin[plus ? "aa_zdisc_tip_gate_plus" :
                           "aa_zdisc_tip_gate_minus"][best.myosin_index] =
                contact.radial_gate;
        }
    }
    stamp_aa_phase(kAaTipAndScatter);
    for (size_t slot : aa_active_link_slots) {
        const auto& link = aa_links[slot];
        if (aa_arm_state_model &&
            !actn2_arm_kinetics::is_bridge(link.arm_state)) continue;
        const auto& mech = mechanics[slot];
        local_forces[link.actin_i] += mech.force_on_first;
        local_forces[link.actin_j] -= mech.force_on_first;
        local_torques[link.actin_i] += mech.torque_on_first;
        local_torques[link.actin_j] += mech.torque_on_second;
        if ((aa_myosin_contact_epsilon > 0.0 ||
             (aa_myosin_contact_shell_width > 0.0 &&
              aa_myosin_contact_shell_energy > 0.0)) &&
            aa_actn2_contact_radius > 0.0) {
            auto apply_contact = [&](
                int myosin_index,
                const ACTN2MyosinContactResult& contact) {
                local_forces[link.actin_i] += contact.force_on_actin_i;
                local_forces[link.actin_j] += contact.force_on_actin_j;
                local_myosin_forces[myosin_index] += contact.force_on_myosin;
                local_torques[link.actin_i] += contact.torque_on_actin_i;
                local_torques[link.actin_j] += contact.torque_on_actin_j;
                local_myosin_torques[myosin_index] += contact.torque_on_myosin;
            };
            if (slot < contact_evaluated.size() &&
                contact_evaluated[slot]) {
                for (const auto& entry : contact_results[slot]) {
                    apply_contact(entry.first, entry.second);
                }
            } else {
                // A bridge associated after phase 1 has no cached contact yet.
                // This is the only path that evaluates contact during force
                // application; pre-existing bridges reuse the kinetic pass.
                std::vector<int> contact_myosins =
                    actin_neighbors_by_species[
                        static_cast<size_t>(link.actin_i)].second;
                const auto& second_neighbors =
                    actin_neighbors_by_species[
                        static_cast<size_t>(link.actin_j)].second;
                contact_myosins.insert(
                    contact_myosins.end(), second_neighbors.begin(),
                    second_neighbors.end());
                std::sort(contact_myosins.begin(), contact_myosins.end());
                contact_myosins.erase(
                    std::unique(contact_myosins.begin(), contact_myosins.end()),
                    contact_myosins.end());
                const auto patch = compute_actn2_patch_geometry(
                    actin, link.actin_i, link.actin_j,
                    link.s_i, link.s_j, box);
                const double outer_contact_distance =
                    myosin.radius + aa_actn2_contact_radius +
                    aa_myosin_contact_shell_width;
                for (int myosin_index : contact_myosins) {
                    if (!actn2_patch_may_contact_myosin_body(
                            patch, myosin, myosin_index, box,
                            outer_contact_distance)) continue;
                    const auto contact = compute_actn2_myosin_contact(
                        actin, myosin, link.actin_i, link.actin_j,
                        link.s_i, link.s_j, myosin_index, box,
                        myosin.radius + aa_actn2_contact_radius,
                        aa_myosin_contact_epsilon,
                        aa_myosin_contact_force_cap,
                        aa_myosin_contact_shell_width,
                        aa_myosin_contact_shell_energy);
                    if (contact.applied) apply_contact(myosin_index, contact);
                }
            }
        }
        local_status[link.actin_i] = std::max(local_status[link.actin_i],
                                              static_cast<int>(link.status));
        local_status[link.actin_j] = std::max(local_status[link.actin_j],
                                              static_cast<int>(link.status));
        ++actin_n_bonds[link.actin_i];
        ++actin_n_bonds[link.actin_j];
        aa_local_link_count_feature[link.actin_i] += 1.0;
        aa_local_link_count_feature[link.actin_j] += 1.0;
        auto& status_feature = link.status == 2
            ? aa_local_status2_count_feature : aa_local_status1_count_feature;
        status_feature[link.actin_i] += 1.0;
        status_feature[link.actin_j] += 1.0;
        const double link_tension = std::max(
            0.0, ramp * aa_link_radial_stiffness(link.status) *
                (mech.length - aa_optimal));
        aa_local_tension_max_feature[link.actin_i] = std::max(
            aa_local_tension_max_feature[link.actin_i], link_tension);
        aa_local_tension_max_feature[link.actin_j] = std::max(
            aa_local_tension_max_feature[link.actin_j], link_tension);
        if (link.status == 2) {
            ++actin_strong_cb_count[link.actin_i];
            ++actin_strong_cb_count[link.actin_j];
        }
        const int first = std::min(link.actin_i, link.actin_j);
        const int second = std::max(link.actin_i, link.actin_j);
        actin_actin_bonds[first][second] = actin_actin_bonds[second][first] = 1;
        actin_actin_status[first][second] = actin_actin_status[second][first] =
            std::max(actin_actin_status[first][second], link.status);
        actin_actin_lifetime[first][second] = actin_actin_lifetime[second][first] =
            std::max(actin_actin_lifetime[first][second],
                static_cast<int>(current_step - link.attach_step + 1));
    active_aa_pairs_temp[thread_id].emplace_back(first, second);
    }
    if (motor_transmitted_mode) {
        aa_motor_local_update_seconds += omp_get_wtime() -
            aa_motor_update_started;
        ++aa_motor_profile_steps;
    }
    }  // omp single
    stamp_aa_phase(0);   // close the final phase
}

void Sarcomere::_process_catch_bonds(int& i) {
    std::vector <int> actin_neighbors = actin_neighbors_by_species[i].first;
    // A capture neighbor list is sufficient for new associations, but an
    // occupied spring must continue to be evaluated after it stretches beyond
    // capture range.  KMC, not neighbor-list loss, owns detachment.
    for (const auto& pair : previous_aa_pairs) {
        if (pair.first == i) actin_neighbors.push_back(pair.second);
        if (pair.second == i) actin_neighbors.push_back(pair.first);
    }
    std::sort(actin_neighbors.begin(), actin_neighbors.end());
    actin_neighbors.erase(
        std::unique(actin_neighbors.begin(), actin_neighbors.end()),
        actin_neighbors.end());
    std::vector<int> statuses;
    std::vector<int> cb_indices;
    for (int index = 0; index < actin_neighbors.size(); index++){
            int j = actin_neighbors[index];
            if (i>j){
                int status = determine_cb_status(i,j);
                if (status > 0){
                    bool acc = _cb_decide(i,j,status);
                    if (acc){
                        cb_indices.push_back(j);
                        statuses.push_back(status);
                    }
                }
            }
    }
    _set_cb(i, cb_indices,statuses);
}

void Sarcomere::_record_aa_state_changes() {
    // The candidate list is rebuilt from scratch every step but its size barely
    // moves, so a persistent buffer keeps the capacity instead of reallocating
    // (and copying previous_aa_pairs into a fresh allocation) each time.
    std::vector<BondPair>& candidates = aa_state_change_candidates;
    candidates.assign(previous_aa_pairs.begin(), previous_aa_pairs.end());
    for (const auto& local_pairs : active_aa_pairs_temp) {
        for (auto pair : local_pairs) {
            if (pair.first > pair.second) std::swap(pair.first, pair.second);
            candidates.push_back(pair);
        }
    }
    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());

    const int thread_id = omp_in_parallel() ? omp_get_thread_num() : 0;
    auto& events = aa_state_events_temp[thread_id];
    auto& am_motor_force_feature = actin["am_motor_force"];
    auto& am_strong_count_feature = actin["am_strong_count"];
    for (const auto& [i, j] : candidates) {
        const int old_bond = static_cast<int>(actin_actin_bonds_prev[i][j]);
        const int new_bond = static_cast<int>(actin_actin_bonds[i][j]);
        const int old_status = old_bond
            ? static_cast<int>(actin_actin_status_prev[i][j]) : 0;
        const int new_status = new_bond
            ? static_cast<int>(actin_actin_status[i][j]) : 0;
        if (old_bond == new_bond && old_status == new_status) continue;

        const double load_i = actin_f_load_cb ? (*actin_f_load_cb)[i] : actin.f_load[i];
        const double load_j = actin_f_load_cb ? (*actin_f_load_cb)[j] : actin.f_load[j];
        events.insert(events.end(), {
            static_cast<double>(i),
            static_cast<double>(j),
            static_cast<double>(current_step),
            static_cast<double>(old_bond),
            static_cast<double>(new_bond),
            static_cast<double>(old_status),
            static_cast<double>(new_status),
            am_motor_force_feature[i],
            am_motor_force_feature[j],
            am_strong_count_feature[i],
            am_strong_count_feature[j],
            load_i,
            load_j});

        const auto force_i = _actin_motor_force_vector(i);
        const auto force_j = _actin_motor_force_vector(j);
        const double motor_only = actn2_load_kinetics::catch_load(
            actn2_load_kinetics::StatusMode::motor_only,
            force_i, force_j, 0.0);
        const double opposed_drive =
            actn2_load_kinetics::opposed_motor_drive(force_i, force_j);
        const double transmitted = _aa_transmitted_load(i, j);
        const double combined_candidate = std::min(opposed_drive, transmitted);
        const double selected =
            aa_status_mode == actn2_load_kinetics::StatusMode::force_projection
            ? _aa_effective_arm_loads(i, j).favorable_bilateral()
            : actn2_load_kinetics::catch_load(
                  aa_status_mode, force_i, force_j, transmitted);
        aa_load_events.insert(aa_load_events.end(), {
            static_cast<double>(i),
            static_cast<double>(j),
            static_cast<double>(current_step),
            static_cast<double>(old_status),
            static_cast<double>(new_status),
            motor_only,
            opposed_drive,
            transmitted,
            combined_candidate,
            selected,
            actin.direction[i].dot(actin.direction[j])});
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
    auto& am_motor_force_feature = actin[key_motor_force];
    auto& am_active_force_x_feature = actin[key_active_force_x];
    auto& am_active_force_y_feature = actin[key_active_force_y];
    auto& am_active_force_z_feature = actin[key_active_force_z];
    auto& am_tether_force_feature = actin[key_tether_force];
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
                    am_tether_force_feature[i] +=
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


void Sarcomere::new_file(){
    create_file(filename, actin, myosin, max_myosin_bonds);
    H5::H5File file(filename, H5F_ACC_RDWR);
    H5::Group meta = file.openGroup("/meta");
    H5::DataSpace scalar(H5S_SCALAR);
    const int mode_value = static_cast<int>(aa_status_mode);
    H5::Attribute mode_attr = meta.createAttribute(
        "aa_status_mode", H5::PredType::NATIVE_INT, scalar);
    mode_attr.write(H5::PredType::NATIVE_INT, &mode_value);
    H5::Attribute on_attr = meta.createAttribute(
        "aa_status2_force_on", H5::PredType::NATIVE_DOUBLE, scalar);
    on_attr.write(H5::PredType::NATIVE_DOUBLE, &aa_status2_force_on);
    H5::Attribute off_attr = meta.createAttribute(
        "aa_status2_force_off", H5::PredType::NATIVE_DOUBLE, scalar);
    off_attr.write(H5::PredType::NATIVE_DOUBLE, &aa_status2_force_off);
    const int aa_load_schema = 2;
    H5::Attribute load_schema_attr = meta.createAttribute(
        "aa_load_transition_schema", H5::PredType::NATIVE_INT, scalar);
    load_schema_attr.write(H5::PredType::NATIVE_INT, &aa_load_schema);
    const int am_model_version = dynamic_am ? 22 : 0;
    H5::Attribute am_model_attr = meta.createAttribute(
        "am_model_version", H5::PredType::NATIVE_INT, scalar);
    am_model_attr.write(H5::PredType::NATIVE_INT, &am_model_version);
    H5::Attribute stall_attr = meta.createAttribute(
        "am_stall_force", H5::PredType::NATIVE_DOUBLE, scalar);
    stall_attr.write(H5::PredType::NATIVE_DOUBLE, &am_stall_force);
    H5::Attribute motor_k_attr = meta.createAttribute(
        "am_motor_k", H5::PredType::NATIVE_DOUBLE, scalar);
    motor_k_attr.write(H5::PredType::NATIVE_DOUBLE, &am_motor_k);
    H5::Attribute stroke_attr = meta.createAttribute(
        "am_power_stroke", H5::PredType::NATIVE_DOUBLE, scalar);
    stroke_attr.write(H5::PredType::NATIVE_DOUBLE, &am_power_stroke);
    H5::Attribute heads_attr = meta.createAttribute(
        "am_heads_per_half", H5::PredType::NATIVE_INT, scalar);
    heads_attr.write(H5::PredType::NATIVE_INT, &am_heads_per_half);
    auto write_int_attr = [&](const char* name, int value) {
        H5::Attribute attr = meta.createAttribute(
            name, H5::PredType::NATIVE_INT, scalar);
        attr.write(H5::PredType::NATIVE_INT, &value);
    };
    auto write_double_attr = [&](const char* name, double value) {
        H5::Attribute attr = meta.createAttribute(
            name, H5::PredType::NATIVE_DOUBLE, scalar);
        attr.write(H5::PredType::NATIVE_DOUBLE, &value);
    };
    auto write_string_attr = [&](const char* name, const std::string& value) {
        H5::StrType type(H5::PredType::C_S1, value.size() + 1);
        type.setStrpad(H5T_STR_NULLTERM);
        H5::Attribute attr = meta.createAttribute(name, type, scalar);
        attr.write(type, value);
    };
    write_string_attr("cli_invocation", cli_invocation);
    write_string_attr("cli_defaults", cli_defaults);
    write_double_attr("box_x", box.size() > 0 ? box[0] : 0.0);
    write_double_attr("box_y", box.size() > 1 ? box[1] : 0.0);
    write_double_attr("box_z", box.size() > 2 ? box[2] : 0.0);
    write_int_attr("periodic_x", is_periodic[0] ? 1 : 0);
    write_int_attr("periodic_y", is_periodic[1] ? 1 : 0);
    write_int_attr("periodic_z", is_periodic[2] ? 1 : 0);
    write_int_attr("actin_turnover_enabled",
                   actin_turnover_enabled ? 1 : 0);
    write_double_attr("actin_turnover_rate", actin_turnover_rate);
    write_double_attr("actin_turnover_start_time",
                      actin_turnover_start_time);
    write_int_attr("actin_turnover_max_placement_attempts",
                   actin_turnover_max_placement_attempts);
    write_double_attr("actin_turnover_clearance_ratio",
                      actin_turnover_clearance_ratio);
    write_string_attr(
        "actin_turnover_semantics",
        "fixed_population_whole_filament_reseed_uniform_hazard_"
        "steric_rejection_preserve_x_polarity");
    write_string_attr(
        "actin_turnover_events_columns",
        "actin,step,age_steps,success,placement_attempts,released_am_heads,"
        "released_actn2,released_titin,old_center_x,old_center_y,"
        "old_center_z,new_center_x,new_center_y,new_center_z,old_direction_x,"
        "old_direction_y,old_direction_z,new_direction_x,new_direction_y,"
        "new_direction_z");
    if (!initial_geometry_source.empty()) {
        write_string_attr("initial_geometry_source", initial_geometry_source);
        write_int_attr("initial_geometry_source_frame_requested",
                       initial_geometry_source_frame_requested);
        write_int_attr("initial_geometry_source_frame_loaded",
                       initial_geometry_source_frame_loaded);
        write_string_attr("initial_geometry_import_scope",
                          "filament_centers_and_directions_only");
    }
    write_double_attr("am_tether_k", am_tether_k);
    write_double_attr("am_tether_occupancy_floor", am_tether_occupancy_floor);
    write_double_attr("am_tether_occupancy_floor_reach",
                      am_tether_occupancy_floor_reach);
    write_double_attr("probe_tilt_deg", probe_tilt_deg);
    write_double_attr("probe_confine_k", probe_confine_k);
    write_double_attr("am_tether_slack", am_tether_slack);
    write_double_attr("am_head_transverse_k", am_head_transverse_k);
    write_double_attr("am_alignment_k", am_alignment_k);
    write_double_attr("am_alignment_slack_rad", am_alignment_slack);
    write_int_attr("am_tether_saturation_heads", am_tether_saturation_heads);
    write_int_attr("am_max_heads_per_actin_half",
                   am_max_heads_per_actin_half);
    write_int_attr("am_sites_per_actin", am_sites_per_actin);
    write_int_attr("am_max_heads_per_actin", am_max_heads_per_actin);
    write_string_attr(
        "am_site_ownership_law",
        "global_material_actin_sites_shared_across_all_myosins_strong_then_low_strain_priority");
    write_string_attr(
        "am_double_overlap_law",
        "wrong_half_heads_may_bind_weakly_but_cannot_powerstroke");
    write_double_attr("am_rupture_distance", am_rupture_distance);
    write_int_attr("am_strain_dependent_detachment",
                   am_strain_dependent_detachment ? 1 : 0);
    write_double_attr("am_strain_detachment_theta",
                      am_strain_detachment_theta);
    write_double_attr("am_transverse_strain_detachment_theta",
                      am_transverse_strain_detachment_theta);
    write_double_attr("am_strain_detachment_max_exponent",
                      am_strain_detachment_max_exponent);
    write_double_attr("am_strong_off_floor_fraction",
                      am_strong_off_floor_fraction);
    write_double_attr("am_cycle_period", am_cycle_period);
    write_double_attr("am_activation_rise", am_activation_rise);
    write_double_attr("am_activation_decay", am_activation_decay);
    write_double_attr("am_activation_floor", am_activation_floor);
    write_double_attr("am_activation_active_fraction",
                      am_activation_active_fraction);
    write_double_attr("am_k_on", am_k_on);
    write_double_attr("am_weak_binding_floor", am_weak_binding_floor);
    write_double_attr("am_k_weak_off", am_k_weak_off);
    write_double_attr("am_k_powerstroke", am_k_powerstroke);
    write_double_attr("am_k_strong_off", am_k_strong_off);
    write_double_attr("am_binding_azimuth_half_width_rad",
                      am_binding_azimuth_half_width);
    write_int_attr("am_binding_azimuth_selection",
                   am_binding_azimuth_selection ? 1 : 0);
    write_double_attr("am_unloaded_velocity", am_unloaded_velocity);
    write_string_attr("am_head_target_semantics",
                      "fixed_actin_material_coordinate_um");
    write_string_attr("am_head_mechanics",
                      "objective_two_spring_anisotropic_strain_kmc_v4");
    write_string_attr(
        "am_head_transverse_law",
        "per_bound_material_head_flat_bottom_transverse_spring_with_independent_escape_coupling");
    write_string_attr(
        "am_head_emergency_rupture_law",
        "connector_excess_beyond_relaxed_reach_only_after_stochastic_kmc");
    write_string_attr("am_head_crown_layout",
                      "axial_crowns_with_three_120deg_sites_and_40deg_twist");
    write_string_attr(
        "am_head_association_axial_gate",
        "material_head_projection_must_lie_inside_finite_actin_segment");
    write_double_attr("am_bare_zone_fraction",
                      geometry::am_bare_zone_fraction);
    write_string_attr("unit_length", "um");
    write_string_attr("unit_time", "s");
    write_string_attr("unit_force", "pN");
    write_string_attr("unit_energy", "pN um");
    write_string_attr("unit_stiffness", "pN/um");
    write_string_attr("unit_torque", "pN um");
    write_int_attr("am_thin_activation_enabled",
                   am_thin_activation_enabled ? 1 : 0);
    write_int_attr("am_spatial_heads", am_spatial_heads ? 1 : 0);
    write_int_attr("sterics_enabled", sterics_enabled ? 1 : 0);
    write_double_attr("myosin_radius", myosin.radius);
    write_double_attr("actin_radius", actin_radius);
    write_double_attr("am_surface_gap", am_surface_gap);
    write_double_attr("am_steric_k", am_steric_k);
    write_double_attr("actin_steric_k", actin_steric_k);
    write_double_attr("myosin_steric_k", myosin_steric_k);
    write_double_attr("myosin_exclusion_distance", myosin_exclusion_distance);
    write_int_attr("myomesin_schema_version", k_mm > 0.0 ? 1 : 0);
    write_double_attr("k_mm", k_mm);
    write_int_attr("myomesin_sites_per_filament",
                   myomesin_sites_per_filament);
    write_double_attr("myomesin_k_on", myomesin_k_on);
    write_double_attr("myomesin_k_off", myomesin_k_off);
    write_double_attr("myomesin_axial_capture", myomesin_axial_capture);
    write_double_attr("myomesin_radial_rest", myomesin_radial_rest);
    write_double_attr("myomesin_radial_capture", myomesin_radial_capture);
    write_double_attr("myomesin_radial_k", myomesin_radial_k);
    write_double_attr("myomesin_alignment_cos", myomesin_alignment_cos);
    write_double_attr("myomesin_bell_distance", myomesin_bell_distance);
    write_double_attr("myomesin_force_max", myomesin_force_max);
    write_string_attr(
        "myomesin_model",
        k_mm > 0.0
            ? "sparse_local_dynamic_finite_valency_bare_zone_modules"
            : "disabled");
    write_string_attr(
        "myomesin_force_law",
        "coarse_module_pairwise_axial_zero_register_plus_transverse_nonzero_rest_spacing");
    write_string_attr(
        "myomesin_bridge_state_columns",
        "id_hi,id_lo,myosin_i,myosin_j,attach_step,axial_offset,"
        "radial_distance,force");
    write_string_attr(
        "myomesin_events_columns",
        "step,event_code,id,myosin_i,myosin_j,axial_offset,"
        "radial_distance,force,age");
    write_string_attr("myomesin_event_codes",
                      "1=association,2=stochastic_release,"
                      "3=invalid_geometry_release");
    write_int_attr("aa_double_overlap_contact_enabled",
                   aa_double_overlap_contact_enabled ? 1 : 0);
    write_double_attr("aa_double_overlap_contact_diameter",
                      aa_double_overlap_contact_diameter);
    write_double_attr("aa_double_overlap_contact_epsilon",
                      aa_double_overlap_contact_epsilon);
    write_double_attr("aa_double_overlap_contact_force_cap",
                      aa_double_overlap_contact_force_cap);
    write_string_attr(
        "aa_double_overlap_contact_law",
        "local_antiparallel_pointed_end_unilateral_WCA_equal_opposite_force_and_torque");
    write_double_attr("am_preferred_centerline", am_preferred_centerline);
    write_double_attr("aa_preferred_centerline", aa_optimal);
    write_double_attr("aa_spacing_ratio", aa_spacing_ratio);
    write_double_attr("aa_capture_margin", aa_capture_margin);
    write_double_attr("aa_antiparallel_overlap", aa_antiparallel_overlap);
    write_double_attr("aa_axial_register_ratio", aa_axial_register_ratio);
    write_int_attr("aa_register_all_statuses",
                   aa_register_all_statuses ? 1 : 0);
    write_double_attr("aa_axial_register_force_max",
                      aa_axial_register_force_max);
    write_double_attr("aa_barbed_zone", aa_barbed_zone);
    write_double_attr("aa_status2_registry_on", aa_status2_registry_on);
    write_double_attr("aa_status2_registry_off", aa_status2_registry_off);
    write_double_attr(
        "aa_register_saturation_error", aa_register_saturation_error);
    write_int_attr("aa_local_multilink", aa_local_multilink ? 1 : 0);
    write_int_attr("aa_prebundle_links_per_pair", aa_prebundle_links_per_pair);
    write_double_attr("aa_site_exclusion", aa_site_exclusion);
    write_double_attr("aa_overlap_capture", aa_overlap_capture);
    write_double_attr("aa_perpendicular_k", aa_perpendicular_k);
    write_double_attr("aa_status1_radial_k", aa_status1_radial_k);
    write_double_attr("aa_status1_site_diffusion", aa_status1_site_diffusion);
    write_int_attr("aa_status2_locked_until_detach",
                   aa_status2_locked_until_detach ? 1 : 0);
    write_double_attr("aa_status2_unlock_rate", aa_status2_unlock_rate);
    write_double_attr("aa_status2_unlock_force_scale",
                      aa_status2_unlock_force_scale);
    write_double_attr("aa_status2_unlock_registry_scale",
                      aa_status2_unlock_registry_scale);
    write_int_attr("aa_arm_state_model", aa_arm_state_model ? 1 : 0);
    write_double_attr("aa_arm_open_rate_weak", aa_arm_open_rate_weak);
    write_double_attr("aa_arm_open_rate_mature", aa_arm_open_rate_mature);
    write_double_attr("aa_arm_bell_distance", aa_arm_bell_distance);
    write_double_attr("aa_tether_off_rate", aa_tether_off_rate);
    write_double_attr("aa_tether_rebind_rate", aa_tether_rebind_rate);
    write_double_attr("aa_tether_capture_distance",
                      aa_tether_capture_distance);
    write_double_attr("aa_arm_rebind_theta", aa_arm_rebind_theta);
    write_double_attr("aa_actn2_contact_radius", aa_actn2_contact_radius);
    write_double_attr("aa_myosin_contact_epsilon",
                      aa_myosin_contact_epsilon);
    write_double_attr("aa_myosin_contact_force_cap",
                      aa_myosin_contact_force_cap);
    write_double_attr("aa_myosin_contact_shell_width",
                      aa_myosin_contact_shell_width);
    write_double_attr("aa_myosin_contact_shell_energy",
                      aa_myosin_contact_shell_energy);
    write_int_attr("aa_zdisc_tip_contact", aa_zdisc_tip_contact ? 1 : 0);
    write_double_attr("aa_zdisc_tip_onset", aa_zdisc_tip_onset);
    write_double_attr("aa_zdisc_tip_stiffness", aa_zdisc_tip_stiffness);
    write_double_attr("aa_zdisc_tip_force_max", aa_zdisc_tip_force_max);
    write_double_attr("aa_barbed_lateral_range", aa_barbed_lateral_range);
    write_double_attr("aa_barbed_lateral_k", aa_barbed_lateral_k);
    write_double_attr("aa_barbed_lateral_force_max",
                      aa_barbed_lateral_force_max);
    write_double_attr("aa_barbed_lateral_axial_window",
                      aa_barbed_lateral_axial_window);
    write_double_attr("aa_zdisc_tip_transverse_flat",
                      aa_zdisc_tip_transverse_flat);
    write_double_attr("aa_zdisc_tip_transverse_cutoff",
                      aa_zdisc_tip_transverse_cutoff);
    write_string_attr("aa_tip_accessibility_scope",
                      "all_nearby_myosin_tips_head_independent");
    write_string_attr("aa_zdisc_tip_nonpenetration",
                      "mobility_weighted_projected_gauss_seidel_g_ge_0");
    write_string_attr("aa_zdisc_tip_projection_traversal",
                      "spatial_incidence_exact_conservative_broadphase");
    write_double_attr("aa_occlusion_rupture_penetration",
                      aa_occlusion_rupture_penetration);
    write_double_attr(
        "aa_status1_radial_k_effective",
        aa_link_radial_stiffness(static_cast<BondFlag>(1)));
    write_double_attr("aa_association_max_direction_dot",
                      aa_association_max_direction_dot);
    write_double_attr("aa_link_rupture_length", aa_link_rupture_length);
    write_int_attr("aa_stiffness_reference_links",
                   aa_stiffness_reference_links);
    write_double_attr(
        "aa_single_link_k",
        actn2_load_kinetics::fixed_single_link_stiffness(
            k_aa, aa_stiffness_reference_links));
    write_double_attr(
        "aa_single_link_kappa",
        actn2_load_kinetics::fixed_single_link_stiffness(
            kappa_aa, aa_stiffness_reference_links));
    write_int_attr("aa_total_link_pool", aa_total_link_pool);
    write_int_attr("aa_neighbor_identity",
                   aa_neighbor_identity_enabled ? 1 : 0);
    write_double_attr("aa_neighbor_slot_memory",
                      aa_neighbor_slot_memory_steps * dt);
    write_double_attr("aa_status2_slot_memory",
                      aa_status2_slot_memory_steps * dt);
    write_int_attr("aa_local_link_schema",
                   aa_local_multilink ? (aa_arm_state_model ? 5 : 4) : 0);
    write_double_attr("titin_k", titin_k);
    write_double_attr("titin_rest_length", titin_rest_length);
    write_int_attr("titin_mode", static_cast<int>(titin_mode));
    // Schema 3 adds the ACTN2-support capture gate and the per-frame
    // association-provenance counters. The bump is deliberate: it makes the
    // resume guard reject every schema-2 trajectory, including the cancelled
    // saturated v4 pilot, rather than silently continuing it under new rules.
    write_int_attr("titin_schema_version",
                   titin_mode == TitinMode::explicit_dynamic ? 4 : 0);
    write_double_attr("titin_k_off_unsupported",
                      explicit_titin_k_off_unsupported);
    write_int_attr("titin_rebind_requires_actn2",
                   explicit_titin_rebind_requires_actn2 ? 1 : 0);
    write_string_attr(
        "titin_anchor_topology",
        titin_mode != TitinMode::explicit_dynamic
            ? "not_applicable"
            : (explicit_titin_k_off_unsupported > 0.0
                ? "serial_actin_actn2_titin_composite_lifetime"
                : "parallel_single_k_off"));
    write_int_attr("titin_require_live_actn2",
                   explicit_titin_require_live_actn2 ? 1 : 0);
    write_string_attr(
        "titin_capture_gate",
        titin_mode != TitinMode::explicit_dynamic
            ? "not_applicable"
            : (explicit_titin_require_live_actn2
                ? "formation_requires_live_actn2_on_barbed_end_any_status"
                : "unrestricted_within_capture_shell"));
    write_string_attr(
        "titin_association_counts_columns",
        "actn2_supported,actn2_unsupported,same_partner_rebind,"
        "unsupported_releases");
    write_string_attr(
        "titin_event_codes",
        "1=association,2=release_with_live_actn2,3=release_without_live_actn2");
    write_int_attr("titin_links_per_half", explicit_titin_per_half);
    write_double_attr("titin_capture_distance",
                      explicit_titin_capture_distance);
    write_double_attr("titin_k_on", explicit_titin_k_on);
    write_double_attr("titin_k_off", explicit_titin_k_off);
    write_int_attr("titin_strain_gated_kinetics",
                   explicit_titin_strain_gated_kinetics ? 1 : 0);
    write_string_attr("titin_association_law",
                      "distance_limited_constant_rate_then_tension_only_spring");
    write_double_attr("titin_transition_theta",
                      explicit_titin_transition_theta);
    write_double_attr("titin_rebind_memory",
                      explicit_titin_rebind_memory);
    write_double_attr("titin_rebind_rate_multiplier",
                      explicit_titin_rebind_rate_multiplier);
    write_string_attr("titin_force_law",
                      titin_mode == TitinMode::explicit_dynamic
                          ? "tension_only_flat_bottom"
                          : "legacy_bidirectional_hookean");
    write_string_attr("titin_capture_semantics",
                      titin_mode == TitinMode::explicit_dynamic
                          ? "material_tip_domain_to_individual_actin_barbed_end"
                          : "legacy_effective_centering");
    write_string_attr(
        "titin_topology",
        titin_mode == TitinMode::explicit_dynamic
            ? "explicit_thick_tip_domains_to_dynamic_actin_barbed_ends"
            : (titin_mode == TitinMode::balanced_owning
                ? "balanced_barbed_end_to_owning_myosin_center"
                : "legacy_status2_crosslinked_partner"));
    write_int_attr("release_step_origin", static_cast<int>(
        std::min(release_step_origin,
                 static_cast<size_t>(std::numeric_limits<int>::max()))));
    // The release origin is exactly the number of fixed-geometry chemistry
    // steps requested by main. Store the protocol name explicitly so a
    // production checkpoint can prove that chemistry and motion were coupled
    // from step zero without relying on that implementation detail.
    write_int_attr("fixed_lattice_preequil_steps", static_cast<int>(
        std::min(release_step_origin,
                 static_cast<size_t>(std::numeric_limits<int>::max()))));
    write_string_attr(
        "aa_local_link_mechanics_columns",
        "connector_x,connector_y,connector_z,length,extension,arm_i,arm_j,"
        "bilateral_load,productive_support,status2_gate,lifetime,single_link_tension,"
        "raw_spring_arm_i,raw_spring_arm_j,allocated_motor_i,allocated_motor_j");
    write_string_attr(
        "aa_local_link_mechanics_v2_columns",
        "connector_x,connector_y,connector_z,length,extension,arm_i,arm_j,"
        "bilateral_load,productive_support,status2_gate,lifetime,single_link_tension,"
        "raw_spring_arm_i,raw_spring_arm_j,allocated_motor_i,allocated_motor_j,"
        "direction_dot,registry_error,registry_threshold,registry_ok,arm_load_ok,"
        "productive_support_ok,actual_status2_gate");
    write_string_attr(
        "aa_local_link_mechanics_v2_frame_offsets_columns",
        "step,identity_start,mechanics_v2_start,count");
    write_string_attr(
        "aa_local_link_arm_state_columns",
        "state,state_since_step,last_partner,arm_i_bound,arm_j_bound,"
        "arm_force_i,arm_force_j,actn2_myosin_contact_force,"
        "actn2_myosin_contact_penetration");
    write_string_attr(
        "aa_arm_state_codes",
        "1=weak_two_arm_bridge,2=mature_two_arm_bridge,"
        "3=tethered_to_actin_i,4=tethered_to_actin_j");
    write_string_attr(
        "aa_arm_events_columns",
        "link_id,step,event_code,old_state,new_state,bound_actin,old_partner,"
        "new_or_opened_actin,bound_site,free_or_opened_site,arm_force,"
        "contact_force,contact_penetration,event_rate");
    write_string_attr(
        "aa_arm_event_codes",
        "1=one_arm_opened,2=free_arm_rebound,3=tether_detached");
    write_string_attr(
        "aa_local_link_motor_support_columns",
        "candidate_myosin_i,candidate_myosin_j,candidate_half_i,"
        "candidate_half_j,support_ema,support_instantaneous,candidate_age,"
        "identity_locked,candidate_since_step,candidate_last_seen_step,"
        "raw_source_support_i,raw_source_support_j,"
        "allocated_source_support_i,allocated_source_support_j,"
        "challenger_myosin_i,challenger_myosin_j,challenger_half_i,"
        "challenger_half_j,challenger_support_ema,challenger_age");
    write_string_attr(
        "aa_motor_pair_events_columns",
        "link_id_hi,link_id_lo,step,event_code,old_myosin_i,old_myosin_j,"
        "old_half_i,old_half_j,new_myosin_i,new_myosin_j,new_half_i,"
        "new_half_j,current_support,best_support,support_ema,candidate_age,"
        "raw_source_support_i,raw_source_support_j,"
        "allocated_source_support_i,allocated_source_support_j");
    write_string_attr(
        "aa_motor_pair_event_codes",
        "1=selected,2=switched,3=dropped_after_support_timeout,"
        "4=locked_at_maturation,5=released_by_actn2_detachment,"
        "6=unlocked_to_weak_state");
    write_string_attr(
        "aa_motor_profile_semantics",
        "local_update_is_wall_time;source_build_and_candidate_scan_are_sums_"
        "of_worker_thread_time;printed_at_process_exit_for_motor_transmitted_register");
    write_string_attr(
        "aa_local_link_events_v2_columns",
        "link_id,actin_i,actin_j,step,event_code,old_status,new_status,"
        "effective_load,lifetime,s_i,s_j,arm_i,arm_j,bilateral_load,"
        "productive_support,status2_gate,single_link_tension,raw_spring_arm_i,"
        "raw_spring_arm_j,allocated_motor_i,allocated_motor_j");
    write_string_attr(
        "aa_local_link_events_v3_columns",
        "link_id,actin_i,actin_j,step,event_code,old_status,new_status,"
        "effective_load,lifetime,s_i,s_j,arm_i,arm_j,bilateral_load,"
        "productive_support,status2_gate,single_link_tension,raw_spring_arm_i,"
        "raw_spring_arm_j,allocated_motor_i,allocated_motor_j,connector_length,"
        "rupture_length,direction_dot,registry_error,registry_threshold,registry_ok,"
        "arm_load_ok,productive_support_ok,actual_status2_gate,detachment_reason,"
        "break_probability,stochastic_draw,geometric_detachment,"
        "stochastic_detachment");
    write_string_attr(
        "aa_local_link_detachment_reason_codes",
        "0=not_a_detachment,1=geometric_contour,2=stochastic,"
        "3=am_overlap_occlusion,4=actn2_myosin_collision");
    write_string_attr(
        "aa_tip_accessibility_events_columns",
        "link_id,step,event_code,actin_i,actin_j,s_i,s_j,nearest_myosin,"
        "tip_sign,axial_gap,transverse_distance,tip_zone_occupied,"
        "motor_engaged");
    write_string_attr(
        "aa_tip_accessibility_event_codes",
        "1=accepted_association,2=accepted_tether_rebinding");
    write_string_attr(
        "aa_catch_load_semantics",
        "per_actin_productive_motor_load_allocated_across_transmitting_arms_"
        "and_bounded_by_raw_arm_force");
    write_string_attr(
        "am_forced_ruptures_v2_columns",
        "actin,myosin,step,head_id,old_state,reason,failure_mask,"
        "segment_distance,rupture_distance,actin_myosin_dot,expected_half,"
        "head_half,productive_axial_overlap,binding_ratio,partial_ratio,"
        "head_connector_distance,head_surface_distance,capture_reach,"
        "retention_reach");
    write_int_attr("physical_rod_mobility", physical_rod_mobility ? 1 : 0);
    write_int_attr("isometric_actin_axial", isometric_actin_axial ? 1 : 0);
    write_double_attr("temperature_kelvin", temperature_kelvin);
    write_double_attr("thermal_beta", thermal_beta);
    write_double_attr("legacy_actin_D_trans", legacy_actin_D_trans);
    write_double_attr("legacy_actin_D_rot", legacy_actin_D_rot);
    write_double_attr("legacy_myosin_D_trans", legacy_myosin_D_trans);
    write_double_attr("legacy_myosin_D_rot", legacy_myosin_D_rot);
    write_double_attr("thermal_noise_scale", thermal_noise_scale_metadata);
    write_double_attr("actin_length", actin.length);
    write_double_attr("myosin_length", myosin.length);
    const double expected_beta = rod_mobility::beta(temperature_kelvin);
    write_int_attr("thermal_beta_overridden",
                   std::abs(thermal_beta - expected_beta) >
                       1e-10 * expected_beta ? 1 : 0);
    write_double_attr("aa_refractory_time", bond_recovery_steps * dt);
    write_double_attr("dt", dt);
    write_double_attr("k_am", k_am);
    write_double_attr("k_aa", k_aa);
    write_double_attr("kappa_am", kappa_am);
    write_double_attr("kappa_aa", kappa_aa);
    write_double_attr("aa_base_lifetime", base_lifetime);
    write_double_attr("aa_directional_base_lifetime", directional_base_lifetime);
    write_double_attr("aa_lifetime_coefficient", lifetime_coeff);
    write_double_attr("aa_motor_pair_support_window",
                      aa_motor_pair_support_window);
    write_double_attr("aa_motor_pair_commit_time",
                      aa_motor_pair_commit_time);
    write_double_attr("aa_motor_pair_switch_ratio",
                      aa_motor_pair_switch_ratio);
    write_double_attr("aa_motor_pair_exchange_time",
                      aa_motor_pair_exchange_time);
    write_double_attr("aa_motor_pair_dropout_time",
                      aa_motor_pair_dropout_time);
    write_string_attr(
        "aa_motor_transmitted_register_semantics",
        "distinct_myosin_pair_opposes_actual_local_actn2_connector_on_both_arms;"
        "weak_state_exchange_then_identity_lock_at_maturation");
    write_int_attr("max_aa_bonds_per_actin", max_aa_bonds_per_actin);
    write_int_attr("aa_energy_gated_association",
                   aa_energy_gated_association ? 1 : 0);
    write_double_attr("aa_transition_state_theta", aa_transition_state_theta);
    write_int_attr("aa_passive_detailed_balance",
                   aa_passive_detailed_balance ? 1 : 0);
    write_double_attr("aa_encounter_k", aa_encounter_k);
    write_double_attr("aa_encounter_kappa", aa_encounter_kappa);
    write_int_attr("interaction_ramp_start_step",
                   static_cast<int>(interaction_ramp_start_step));
    write_int_attr("interaction_ramp_steps",
                   static_cast<int>(interaction_ramp_steps));
    auto write_mobility = [&](const char* prefix,
                              const rod_mobility::Diffusion& mobility) {
        const std::string base(prefix);
        write_double_attr((base + "_hydrodynamic_length").c_str(), mobility.length);
        write_double_attr((base + "_hydrodynamic_diameter").c_str(), mobility.diameter);
        write_double_attr((base + "_viscosity_pas").c_str(), mobility.viscosity);
        write_double_attr((base + "_crowding_parallel").c_str(), mobility.crowding_parallel);
        write_double_attr((base + "_crowding_perpendicular").c_str(), mobility.crowding_perpendicular);
        write_double_attr((base + "_crowding_tumble").c_str(), mobility.crowding_tumble);
        write_double_attr((base + "_crowding_roll").c_str(), mobility.crowding_roll);
        write_double_attr((base + "_D_parallel").c_str(), mobility.parallel);
        write_double_attr((base + "_D_perpendicular").c_str(), mobility.perpendicular);
        write_double_attr((base + "_D_tumble").c_str(), mobility.tumble);
        write_double_attr((base + "_D_roll").c_str(), mobility.roll);
    };
    write_mobility("actin", actin_rod_mobility);
    write_mobility("myosin", myosin_rod_mobility);
    write_int_attr("am_thin_feedback_enabled",
                   am_thin_feedback_enabled ? 1 : 0);
    write_double_attr("am_q_calcium_open_rate", am_q_calcium_open_rate);
    write_double_attr("am_q_calcium_half", am_q_calcium_half);
    write_double_attr("am_q_calcium_hill", am_q_calcium_hill);
    write_double_attr("am_q_strong_feedback_rate", am_q_strong_feedback_rate);
    write_double_attr("am_q_strong_half", am_q_strong_half);
    write_double_attr("am_q_strong_hill", am_q_strong_hill);
    write_double_attr("am_q_close_rate", am_q_close_rate);
    write_int_attr("am_q_initialized_at_calcium_equilibrium",
                   am_q_initialized_at_calcium_equilibrium ? 1 : 0);
    write_double_attr("am_initial_q", am_initial_q);
    write_int_attr("am_thick_feedback_enabled",
                   am_thick_feedback_enabled ? 1 : 0);
    write_double_attr("am_thick_recruit_rate", am_thick_recruit_rate);
    write_double_attr("am_thick_fold_rate", am_thick_fold_rate);
    write_double_attr("am_thick_tension_scale", am_thick_tension_scale);
    write_double_attr("am_thick_initial_on_fraction",
                      am_thick_initial_on_fraction);
    write_double_attr("am_thick_force_clamp", am_thick_force_clamp);
    file.flush(H5F_SCOPE_GLOBAL);
}

void Sarcomere::save_state(){
    auto bondData =
        _extract_bonded_pairs(actin_actin_bonds, actin_actin_status, myosinIndicesPerActin);
    std::vector<double> flatActinBonds = std::get<0>(bondData);
    std::vector<double> flatMyosinBonds = std::get<1>(bondData);
    std::vector<double> flatActinMyosinBonds = std::get<2>(bondData);
    append_to_file(filename, actin, myosin, flatActinBonds,
                   flatMyosinBonds, flatActinMyosinBonds, max_myosin_bonds);
    
    // Save compact bond state for resume using pair lists instead of dense matrices.
    {
        H5::H5File file(filename, H5F_ACC_RDWR);
        H5::Group group_state;
        try {
            group_state = file.openGroup("/state");
        } catch (H5::Exception&) {
            group_state = file.createGroup("/state");
        }
        auto ensure_state_dataset_int = [&](const std::string& dataset_name, hsize_t width) {
            const std::string full_path = "/state/" + dataset_name;
            if (!file.nameExists(full_path)) {
                std::vector<hsize_t> initDims = {0, width};
                std::vector<hsize_t> maxDims = {H5S_UNLIMITED, width};
                std::vector<hsize_t> chunkDims = {10, width};
                create_empty_dataset_int(file, "/state", dataset_name, initDims, maxDims, chunkDims);
            }
        };
        auto ensure_state_dataset_double = [&](const std::string& dataset_name, hsize_t width) {
            const std::string full_path = "/state/" + dataset_name;
            if (!file.nameExists(full_path)) {
                std::vector<hsize_t> initDims = {0, width};
                std::vector<hsize_t> maxDims = {H5S_UNLIMITED, width};
                std::vector<hsize_t> chunkDims = {10, width};
                create_empty_dataset(file, "/state", dataset_name, initDims, maxDims, chunkDims);
            }
        };

        const hsize_t actin_width = static_cast<hsize_t>(actin.n);
        const hsize_t myosin_width = static_cast<hsize_t>(myosin.n);
        const hsize_t max_aa_pairs = static_cast<hsize_t>(actin.n) *
                                     static_cast<hsize_t>(std::max(10, 2 * max_myosin_bonds));
        const hsize_t max_am_pairs = static_cast<hsize_t>(myosin.n) *
                                     static_cast<hsize_t>(max_myosin_bonds);
        const hsize_t max_recovery_pairs = max_aa_pairs;

        auto ensure_state_pair_dataset = [&](const std::string& dataset_name, hsize_t rows, hsize_t cols) {
            const std::string full_path = "/state/" + dataset_name;
            if (!file.nameExists(full_path)) {
                std::vector<hsize_t> initDims = {0, rows, cols};
                std::vector<hsize_t> maxDims = {H5S_UNLIMITED, rows, cols};
                std::vector<hsize_t> chunkDims = {10, rows, cols};
                create_empty_dataset_int(file, "/state", dataset_name, initDims, maxDims, chunkDims);
            }
        };
        auto ensure_state_pair_dataset_double = [&](const std::string& dataset_name, hsize_t rows, hsize_t cols) {
            const std::string full_path = "/state/" + dataset_name;
            if (!file.nameExists(full_path)) {
                std::vector<hsize_t> initDims = {0, rows, cols};
                std::vector<hsize_t> maxDims = {H5S_UNLIMITED, rows, cols};
                std::vector<hsize_t> chunkDims = {10, rows, cols};
                create_empty_dataset(file, "/state", dataset_name, initDims, maxDims, chunkDims);
            }
        };

        ensure_state_dataset_int("current_step", 1);
        ensure_state_dataset_int("state_commit_count", 1);
        ensure_state_dataset_int("aa_pairs_prev_count", 1);
        ensure_state_dataset_int("aa_pairs_current_count", 1);
        ensure_state_dataset_int("am_pairs_prev_count", 1);
        ensure_state_dataset_int("am_pairs_current_count", 1);
        ensure_state_pair_dataset("aa_pairs_prev", max_aa_pairs, 2);
        ensure_state_pair_dataset("aa_status_prev", max_aa_pairs, 1);
        ensure_state_pair_dataset("aa_lifetime_prev", max_aa_pairs, 1);
        ensure_state_pair_dataset("aa_pairs_current", max_aa_pairs, 2);
        ensure_state_pair_dataset("aa_status_current", max_aa_pairs, 1);
        ensure_state_pair_dataset("aa_lifetime_current", max_aa_pairs, 1);
        ensure_state_pair_dataset("aa_attach_step_current", max_aa_pairs, 1);
        ensure_state_pair_dataset("aa_status2_memory_until_current", max_aa_pairs, 1);
        ensure_state_pair_dataset("am_pairs_prev", max_am_pairs, 2);
        ensure_state_pair_dataset("am_pairs_current", max_am_pairs, 2);
        ensure_state_dataset_int("am_state_current", max_am_pairs);
        ensure_state_dataset_double("am_target_axial_current", max_am_pairs);
        const hsize_t heads_per_myosin =
            static_cast<hsize_t>(2 * am_heads_per_half);
        ensure_state_pair_dataset("am_head_state_current", myosin_width,
                                  heads_per_myosin);
        ensure_state_pair_dataset("am_head_bound_actin_current", myosin_width,
                                  heads_per_myosin);
        ensure_state_pair_dataset_double("am_head_target_axial_current", myosin_width,
                                         heads_per_myosin);
        ensure_state_pair_dataset_double("am_head_binding_azimuth_current",
                                         myosin_width, heads_per_myosin);
        ensure_state_pair_dataset_double("am_head_stroke_progress_current",
                                         myosin_width, heads_per_myosin);
        ensure_state_pair_dataset_double("am_head_strain_energy_current",
                                         myosin_width, heads_per_myosin);
        ensure_state_pair_dataset_double("am_head_detachment_rate_current",
                                         myosin_width, heads_per_myosin);
        ensure_state_pair_dataset_double("myosin_material_e1_current",
                                         myosin_width, 3);
        ensure_state_pair_dataset_double("myosin_material_e2_current",
                                         myosin_width, 3);
        ensure_state_dataset_int("aa_recovery_count", 1);
        ensure_state_pair_dataset("aa_recovery_pairs", max_recovery_pairs, 2);
        ensure_state_pair_dataset("aa_recovery_until_values", max_recovery_pairs, 1);
        ensure_state_dataset_double("neighbor_last_actin_x", actin_width);
        ensure_state_dataset_double("neighbor_last_actin_y", actin_width);
        ensure_state_dataset_double("neighbor_last_actin_z", actin_width);
        ensure_state_dataset_double("neighbor_last_myosin_x", myosin_width);
        ensure_state_dataset_double("neighbor_last_myosin_y", myosin_width);
        ensure_state_dataset_double("neighbor_last_myosin_z", myosin_width);
        ensure_state_dataset_double("am_q_current", actin_width);
        ensure_state_dataset_int("actin_last_turnover_step_current",
                                 actin_width);
        const hsize_t titin_width =
            static_cast<hsize_t>(explicit_titins.size());
        if (titin_mode == TitinMode::explicit_dynamic) {
            ensure_state_dataset_int("titin_bound_actin_current", titin_width);
            ensure_state_dataset_int("titin_last_bound_actin_current", titin_width);
            ensure_state_dataset_int("titin_attach_step_current", titin_width);
            ensure_state_dataset_int("titin_last_detach_step_current", titin_width);
            ensure_state_dataset_double("titin_distance_current", titin_width);
            ensure_state_dataset_double("titin_extension_current", titin_width);
            ensure_state_dataset_double("titin_force_current", titin_width);
            ensure_state_dataset_double(
                "titin_attachment_distance_current", titin_width);
            ensure_state_dataset_double(
                "titin_attachment_energy_current", titin_width);
            ensure_state_dataset_double(
                "titin_attachment_force_current", titin_width);
            ensure_state_dataset_int("titin_association_counts", 4);
        }

        // Save current_step (scalar value)
        std::vector<int> step_vec = {static_cast<int>(current_step)};
        append_to_dataset_int(group_state, "current_step", step_vec, {1, 1});

        auto aa_prev = pack_aa_bond_state(
            previous_aa_pairs, actin_actin_bonds_prev,
            actin_actin_status_prev, actin_actin_lifetime_prev);
        auto aa_current = pack_aa_bond_state(
            active_aa_pairs, actin_actin_bonds,
            actin_actin_status, actin_actin_lifetime);
        auto am_prev = pack_am_bond_pairs(previous_am_pairs, am_bonds_prev);
        auto am_current = pack_am_bond_pairs(active_am_pairs, am_bonds);
        auto aa_recovery = pack_recovery_state(
            active_recovery_pairs, actin_recovery_until, current_step);

        const hsize_t aa_prev_count = static_cast<hsize_t>(aa_prev.status.size());
        const hsize_t aa_current_count = static_cast<hsize_t>(aa_current.status.size());
        const hsize_t am_prev_count = static_cast<hsize_t>(am_prev.size() / 2);
        const hsize_t am_current_count = static_cast<hsize_t>(am_current.size() / 2);
        const hsize_t aa_recovery_count = static_cast<hsize_t>(aa_recovery.until.size());
        if (aa_prev_count > max_aa_pairs || aa_current_count > max_aa_pairs ||
            am_prev_count > max_am_pairs || am_current_count > max_am_pairs ||
            aa_recovery_count > max_recovery_pairs) {
            // Degrade gracefully: the compact /state bond snapshot is truncated to capacity
            // (see pad_int_vector) and its stored counts are clamped below. Trajectory and
            // /actin per-frame fields are unaffected; only a later resume would see the
            // over-capacity bonds dropped from this checkpoint.
            static bool warned_state_capacity = false;
            if (!warned_state_capacity) {
                fprintf(stderr,
                        "Warning: bond pair list exceeds compact /state capacity at step %llu "
                        "(aa<=%llu/%llu, am<=%llu/%llu, rec=%llu/%llu); truncating resume snapshot, "
                        "trajectory data unaffected.\n",
                        static_cast<unsigned long long>(current_step),
                        static_cast<unsigned long long>(std::max(aa_prev_count, aa_current_count)),
                        static_cast<unsigned long long>(max_aa_pairs),
                        static_cast<unsigned long long>(std::max(am_prev_count, am_current_count)),
                        static_cast<unsigned long long>(max_am_pairs),
                        static_cast<unsigned long long>(aa_recovery_count),
                        static_cast<unsigned long long>(max_recovery_pairs));
                warned_state_capacity = true;
            }
        }

        std::vector<int> aa_prev_pairs = aa_prev.pairs;
        std::vector<int> aa_prev_status = aa_prev.status;
        std::vector<int> aa_prev_lifetime = aa_prev.lifetime;
        std::vector<int> aa_current_pairs = aa_current.pairs;
        std::vector<int> aa_current_status = aa_current.status;
        std::vector<int> aa_current_lifetime = aa_current.lifetime;
        std::vector<int> aa_current_attach_step;
        std::vector<int> aa_current_status2_memory_until;
        aa_current_attach_step.reserve(aa_current.status.size());
        aa_current_status2_memory_until.reserve(aa_current.status.size());
        for (size_t idx = 0; idx + 1 < aa_current.pairs.size(); idx += 2) {
            const int a = aa_current.pairs[idx];
            const int b = aa_current.pairs[idx + 1];
            aa_current_attach_step.push_back(aa_attach_step[a][b]);
            aa_current_status2_memory_until.push_back(
                static_cast<int>(aa_status2_memory_until[a][b]));
        }
        std::vector<int> am_prev_pairs = am_prev;
        std::vector<int> am_current_pairs = am_current;
        std::vector<int> am_current_state;
        std::vector<double> am_current_target;
        std::vector<int> am_head_state_current;
        std::vector<int> am_head_bound_actin_current;
        std::vector<double> am_head_target_current;
        std::vector<double> am_head_binding_azimuth_current;
        std::vector<double> am_head_stroke_progress_current;
        std::vector<double> am_head_strain_energy_current;
        std::vector<double> am_head_detachment_rate_current;
        am_current_state.reserve(am_current_count);
        am_current_target.reserve(am_current_count);
        for (size_t idx = 0; idx + 1 < am_current.size(); idx += 2) {
            const int a = am_current[idx];
            const int m = am_current[idx + 1];
            am_current_state.push_back(dynamic_am ? static_cast<int>(am_state[a][m]) : 2);
            am_current_target.push_back(dynamic_am ? am_target_axial[a][m]
                                                   : _am_axial_coordinate(a, m));
        }
        const size_t total_head_count = static_cast<size_t>(myosin.n) *
                                        static_cast<size_t>(2 * am_heads_per_half);
        am_head_state_current.reserve(total_head_count);
        am_head_bound_actin_current.reserve(total_head_count);
        am_head_target_current.reserve(total_head_count);
        am_head_binding_azimuth_current.reserve(total_head_count);
        am_head_stroke_progress_current.reserve(total_head_count);
        am_head_strain_energy_current.reserve(total_head_count);
        am_head_detachment_rate_current.reserve(total_head_count);
        for (int m = 0; m < myosin.n; ++m) {
            const auto& myosin_halves =
                am_heads_by_myosin[static_cast<size_t>(m)];
            for (int half_index = 0; half_index < 2; ++half_index) {
                const auto& half_heads =
                    myosin_halves[static_cast<size_t>(half_index)];
                for (int h = 0; h < static_cast<int>(half_heads.size()); ++h) {
                    const auto& head = half_heads[static_cast<size_t>(h)];
                    am_head_state_current.push_back(static_cast<int>(head.state));
                    am_head_bound_actin_current.push_back(head.bound_actin);
                    am_head_target_current.push_back(head.target_axial);
                    am_head_binding_azimuth_current.push_back(
                        head.binding_azimuth);
                    am_head_stroke_progress_current.push_back(
                        head.stroke_progress);
                    const auto state = static_cast<am_motor_kinetics::State>(
                        head.state);
                    double strain_energy = 0.0;
                    double detachment_rate = 0.0;
                    if (head.bound_actin >= 0 &&
                        (state == am_motor_kinetics::State::weak ||
                         state == am_motor_kinetics::State::strong)) {
                        const int head_id = _am_head_id(half_index, h);
                        strain_energy = _am_head_strain_energy(
                            head.bound_actin, m, head_id, state);
                        double chemical_rate = am_k_weak_off;
                        if (state == am_motor_kinetics::State::strong) {
                            chemical_rate =
                                am_motor_kinetics::load_slowed_rate_with_floor(
                                am_k_strong_off,
                                _am_head_active_force(
                                    head.bound_actin, m, head_id),
                                am_force_scale,
                                am_strong_off_floor_fraction);
                        }
                        detachment_rate = _am_head_detachment_rate(
                            head.bound_actin, m, head_id, state,
                            chemical_rate);
                    }
                    am_head_strain_energy_current.push_back(strain_energy);
                    am_head_detachment_rate_current.push_back(detachment_rate);
                }
            }
        }
        std::vector<int> aa_recovery_pairs = aa_recovery.pairs;
        std::vector<int> aa_recovery_until_values = aa_recovery.until;
        std::vector<int> titin_bound_actin_current;
        std::vector<int> titin_last_bound_actin_current;
        std::vector<int> titin_attach_step_current;
        std::vector<int> titin_last_detach_step_current;
        std::vector<double> titin_distance_current;
        std::vector<double> titin_extension_current;
        std::vector<double> titin_force_current;
        std::vector<double> titin_attachment_distance_current;
        std::vector<double> titin_attachment_energy_current;
        std::vector<double> titin_attachment_force_current;
        if (titin_mode == TitinMode::explicit_dynamic) {
            for (const auto& titin : explicit_titins) {
                titin_bound_actin_current.push_back(
                    titin.bound ? titin.bound_actin : -1);
                titin_last_bound_actin_current.push_back(
                    titin.last_bound_actin);
                titin_attach_step_current.push_back(
                    static_cast<int>(titin.attach_step));
                titin_last_detach_step_current.push_back(
                    static_cast<int>(titin.last_detach_step));
                titin_distance_current.push_back(titin.distance);
                titin_extension_current.push_back(titin.extension);
                titin_force_current.push_back(titin.force);
                titin_attachment_distance_current.push_back(
                    titin.attachment_distance);
                titin_attachment_energy_current.push_back(
                    titin.attachment_energy);
                titin_attachment_force_current.push_back(
                    titin.attachment_force);
            }
        }

        pad_int_vector(aa_prev_pairs, static_cast<size_t>(max_aa_pairs * 2), -1);
        pad_int_vector(aa_prev_status, static_cast<size_t>(max_aa_pairs), 0);
        pad_int_vector(aa_prev_lifetime, static_cast<size_t>(max_aa_pairs), 0);
        pad_int_vector(aa_current_pairs, static_cast<size_t>(max_aa_pairs * 2), -1);
        pad_int_vector(aa_current_status, static_cast<size_t>(max_aa_pairs), 0);
        pad_int_vector(aa_current_lifetime, static_cast<size_t>(max_aa_pairs), 0);
        pad_int_vector(aa_current_attach_step, static_cast<size_t>(max_aa_pairs), -1);
        pad_int_vector(aa_current_status2_memory_until,
                       static_cast<size_t>(max_aa_pairs), 0);
        pad_int_vector(am_prev_pairs, static_cast<size_t>(max_am_pairs * 2), -1);
        pad_int_vector(am_current_pairs, static_cast<size_t>(max_am_pairs * 2), -1);
        pad_int_vector(am_current_state, static_cast<size_t>(max_am_pairs), 0);
        am_current_target.resize(static_cast<size_t>(max_am_pairs), 0.0);
        pad_int_vector(am_head_state_current, total_head_count, 0);
        pad_int_vector(am_head_bound_actin_current, total_head_count, -1);
        am_head_target_current.resize(total_head_count, 0.0);
        am_head_binding_azimuth_current.resize(total_head_count, 0.0);
        am_head_stroke_progress_current.resize(total_head_count, 0.0);
        am_head_strain_energy_current.resize(total_head_count, 0.0);
        am_head_detachment_rate_current.resize(total_head_count, 0.0);
        pad_int_vector(aa_recovery_pairs, static_cast<size_t>(max_recovery_pairs * 2), -1);
        pad_int_vector(aa_recovery_until_values, static_cast<size_t>(max_recovery_pairs), 0);

        append_to_dataset_int(group_state, "aa_pairs_prev_count", {static_cast<int>(std::min(aa_prev_count, max_aa_pairs))}, {1, 1});
        append_to_dataset_int(group_state, "aa_pairs_current_count", {static_cast<int>(std::min(aa_current_count, max_aa_pairs))}, {1, 1});
        append_to_dataset_int(group_state, "am_pairs_prev_count", {static_cast<int>(std::min(am_prev_count, max_am_pairs))}, {1, 1});
        append_to_dataset_int(group_state, "am_pairs_current_count", {static_cast<int>(std::min(am_current_count, max_am_pairs))}, {1, 1});
        append_to_dataset_int(group_state, "aa_recovery_count", {static_cast<int>(std::min(aa_recovery_count, max_recovery_pairs))}, {1, 1});
        append_to_dataset_int(group_state, "aa_pairs_prev", aa_prev_pairs, {1, max_aa_pairs, 2});
        append_to_dataset_int(group_state, "aa_status_prev", aa_prev_status, {1, max_aa_pairs, 1});
        append_to_dataset_int(group_state, "aa_lifetime_prev", aa_prev_lifetime, {1, max_aa_pairs, 1});
        append_to_dataset_int(group_state, "aa_pairs_current", aa_current_pairs, {1, max_aa_pairs, 2});
        append_to_dataset_int(group_state, "aa_status_current", aa_current_status, {1, max_aa_pairs, 1});
        append_to_dataset_int(group_state, "aa_lifetime_current", aa_current_lifetime, {1, max_aa_pairs, 1});
        append_to_dataset_int(group_state, "aa_attach_step_current", aa_current_attach_step, {1, max_aa_pairs, 1});
        append_to_dataset_int(group_state, "aa_status2_memory_until_current",
                              aa_current_status2_memory_until,
                              {1, max_aa_pairs, 1});
        append_to_dataset_int(group_state, "am_pairs_prev", am_prev_pairs, {1, max_am_pairs, 2});
        append_to_dataset_int(group_state, "am_pairs_current", am_current_pairs, {1, max_am_pairs, 2});
        append_to_dataset_int(group_state, "am_state_current", am_current_state, {1, max_am_pairs});
        append_to_dataset(group_state, "am_target_axial_current", am_current_target, {1, max_am_pairs});
        append_to_dataset_int(group_state, "am_head_state_current", am_head_state_current,
                              {1, myosin_width, heads_per_myosin});
        append_to_dataset_int(group_state, "am_head_bound_actin_current",
                              am_head_bound_actin_current,
                              {1, myosin_width, heads_per_myosin});
        append_to_dataset(group_state, "am_head_target_axial_current", am_head_target_current,
                          {1, myosin_width, heads_per_myosin});
        append_to_dataset(group_state, "am_head_binding_azimuth_current",
                          am_head_binding_azimuth_current,
                          {1, myosin_width, heads_per_myosin});
        append_to_dataset(group_state, "am_head_stroke_progress_current",
                          am_head_stroke_progress_current,
                          {1, myosin_width, heads_per_myosin});
        append_to_dataset(group_state, "am_head_strain_energy_current",
                          am_head_strain_energy_current,
                          {1, myosin_width, heads_per_myosin});
        append_to_dataset(group_state, "am_head_detachment_rate_current",
                          am_head_detachment_rate_current,
                          {1, myosin_width, heads_per_myosin});
        append_to_dataset(group_state, "myosin_material_e1_current",
                          flatten_vec_values(myosin_material_e1),
                          {1, myosin_width, 3});
        append_to_dataset(group_state, "myosin_material_e2_current",
                          flatten_vec_values(myosin_material_e2),
                          {1, myosin_width, 3});
        append_to_dataset_int(group_state, "aa_recovery_pairs", aa_recovery_pairs,
                              {1, max_recovery_pairs, 2});
        append_to_dataset_int(group_state, "aa_recovery_until_values", aa_recovery_until_values,
                              {1, max_recovery_pairs, 1});

        std::vector<double> neighbor_last_actin_x;
        std::vector<double> neighbor_last_actin_y;
        std::vector<double> neighbor_last_actin_z;
        std::vector<double> neighbor_last_myosin_x;
        std::vector<double> neighbor_last_myosin_y;
        std::vector<double> neighbor_last_myosin_z;
        neighbor_list.get_last_species_positions(
            neighbor_last_actin_x, neighbor_last_actin_y, neighbor_last_actin_z,
            neighbor_last_myosin_x, neighbor_last_myosin_y, neighbor_last_myosin_z);
        append_to_dataset(group_state, "neighbor_last_actin_x", neighbor_last_actin_x, {1, actin_width});
        append_to_dataset(group_state, "neighbor_last_actin_y", neighbor_last_actin_y, {1, actin_width});
        append_to_dataset(group_state, "neighbor_last_actin_z", neighbor_last_actin_z, {1, actin_width});
        append_to_dataset(group_state, "neighbor_last_myosin_x", neighbor_last_myosin_x, {1, myosin_width});
        append_to_dataset(group_state, "neighbor_last_myosin_y", neighbor_last_myosin_y, {1, myosin_width});
        append_to_dataset(group_state, "neighbor_last_myosin_z", neighbor_last_myosin_z, {1, myosin_width});
        append_to_dataset(group_state, "am_q_current", am_thin_activation,
                          {1, actin_width});
        std::vector<int> turnover_steps(actin_last_turnover_step.begin(),
                                        actin_last_turnover_step.end());
        append_to_dataset_int(
            group_state, "actin_last_turnover_step_current", turnover_steps,
            {1, actin_width});
        if (titin_mode == TitinMode::explicit_dynamic) {
            append_to_dataset_int(group_state, "titin_bound_actin_current",
                                  titin_bound_actin_current,
                                  {1, titin_width});
            append_to_dataset_int(group_state,
                                  "titin_last_bound_actin_current",
                                  titin_last_bound_actin_current,
                                  {1, titin_width});
            append_to_dataset_int(group_state, "titin_attach_step_current",
                                  titin_attach_step_current,
                                  {1, titin_width});
            append_to_dataset_int(group_state,
                                  "titin_last_detach_step_current",
                                  titin_last_detach_step_current,
                                  {1, titin_width});
            append_to_dataset(group_state, "titin_distance_current",
                              titin_distance_current, {1, titin_width});
            append_to_dataset(group_state, "titin_extension_current",
                              titin_extension_current, {1, titin_width});
            append_to_dataset(group_state, "titin_force_current",
                              titin_force_current, {1, titin_width});
            append_to_dataset(group_state,
                              "titin_attachment_distance_current",
                              titin_attachment_distance_current,
                              {1, titin_width});
            append_to_dataset(group_state,
                              "titin_attachment_energy_current",
                              titin_attachment_energy_current,
                              {1, titin_width});
            append_to_dataset(group_state,
                              "titin_attachment_force_current",
                              titin_attachment_force_current,
                              {1, titin_width});
            append_to_dataset_int(
                group_state, "titin_association_counts",
                {static_cast<int>(explicit_titin_supported_associations),
                 static_cast<int>(explicit_titin_unsupported_associations),
                 static_cast<int>(explicit_titin_rebind_associations),
                 static_cast<int>(explicit_titin_unsupported_releases)},
                {1, 4});
        }

        // Save RNG states so resumed trajectories can be bitwise reproducible.
        if (rng != nullptr) {
            const size_t main_state_size = gsl_rng_size(rng);
            if (main_state_size > 0) {
                ensure_state_dataset_int("rng_main_state", static_cast<hsize_t>(main_state_size));
                const auto* state_ptr = static_cast<const unsigned char*>(gsl_rng_state(rng));
                std::vector<int> main_state(main_state_size, 0);
                for (size_t idx = 0; idx < main_state_size; ++idx) {
                    main_state[idx] = static_cast<int>(state_ptr[idx]);
                }
                append_to_dataset_int(group_state, "rng_main_state", main_state,
                                      {1, static_cast<hsize_t>(main_state_size)});
            }
        }
        if (!rng_engines.empty() && rng_engines[0] != nullptr) {
            const int thread_count = static_cast<int>(rng_engines.size());
            const size_t thread_state_size = gsl_rng_size(rng_engines[0]);
            const hsize_t flat_width = static_cast<hsize_t>(thread_count) *
                                       static_cast<hsize_t>(thread_state_size);
            if (thread_state_size > 0) {
                ensure_state_dataset_int("rng_thread_state", flat_width);
                ensure_state_dataset_int("rng_thread_count", 1);
                ensure_state_dataset_int("rng_thread_state_size", 1);

                std::vector<int> thread_state(flat_width, 0);
                for (int t = 0; t < thread_count; ++t) {
                    if (rng_engines[t] == nullptr) {
                        continue;
                    }
                    const size_t local_size = gsl_rng_size(rng_engines[t]);
                    if (local_size != thread_state_size) {
                        continue;
                    }
                    const auto* local_ptr =
                        static_cast<const unsigned char*>(gsl_rng_state(rng_engines[t]));
                    size_t base = static_cast<size_t>(t) * thread_state_size;
                    for (size_t idx = 0; idx < thread_state_size; ++idx) {
                        thread_state[base + idx] = static_cast<int>(local_ptr[idx]);
                    }
                }
                append_to_dataset_int(group_state, "rng_thread_state", thread_state, {1, flat_width});
                append_to_dataset_int(group_state, "rng_thread_count", {thread_count}, {1, 1});
                append_to_dataset_int(group_state, "rng_thread_state_size",
                                      {static_cast<int>(thread_state_size)}, {1, 1});
            }
        }
        maybe_test_sleep("SARCOMERE_TEST_SLEEP_BEFORE_STATE_COMMIT", "state_commit");
        append_to_dataset_int(group_state, "state_commit_count", {1}, {1, 1});
        file.flush(H5F_SCOPE_GLOBAL);
    }

    // Sparse, frame-ragged local ACTN2 state.  The offset table maps each
    // saved step to a contiguous slice of identity and coordinate records,
    // avoiding an O(N_actin^2) trajectory representation.
    if (aa_local_multilink) {
        H5::H5File file(filename, H5F_ACC_RDWR);
        H5::Group group_cb(file.openGroup("/catch_bond"));
        if (!file.nameExists("/catch_bond/local_link_identity")) {
            create_empty_dataset_int(
                file, "/catch_bond", "local_link_identity",
                {0, 7}, {H5S_UNLIMITED, 7}, {256, 7});
            create_empty_dataset(
                file, "/catch_bond", "local_link_coordinates",
                {0, 2}, {H5S_UNLIMITED, 2}, {256, 2});
            create_empty_dataset_int(
                file, "/catch_bond", "local_link_frame_offsets",
                {0, 3}, {H5S_UNLIMITED, 3}, {64, 3});
        }
        // One row per identity record. The primary arm loads are the
        // motor-allocated, mechanically transmitted loads used by KMC. Raw
        // spring-arm forces and the allocated motor budgets are appended for
        // auditing and must never be confused with the lifetime input.
        if (!file.nameExists("/catch_bond/local_link_mechanics")) {
            create_empty_dataset(
                file, "/catch_bond", "local_link_mechanics",
                {0, 16}, {H5S_UNLIMITED, 16}, {256, 16});
        }
        if (!file.nameExists("/catch_bond/local_link_mechanics_v2")) {
            create_empty_dataset(
                file, "/catch_bond", "local_link_mechanics_v2",
                {0, 23}, {H5S_UNLIMITED, 23}, {256, 23});
        }
        if (!file.nameExists(
                "/catch_bond/local_link_mechanics_v2_frame_offsets")) {
            create_empty_dataset_int(
                file, "/catch_bond", "local_link_mechanics_v2_frame_offsets",
                {0, 4}, {H5S_UNLIMITED, 4}, {64, 4});
        }
        if (!file.nameExists("/catch_bond/local_link_arm_state")) {
            create_empty_dataset(
                file, "/catch_bond", "local_link_arm_state",
                {0, 9}, {H5S_UNLIMITED, 9}, {256, 9});
        }
        if (aa_status_mode ==
                actn2_load_kinetics::StatusMode::motor_transmitted_register &&
            !file.nameExists("/catch_bond/local_link_motor_support")) {
            // Rows align one-for-one with local_link_identity.
            create_empty_dataset(
                file, "/catch_bond", "local_link_motor_support",
                {0, 20}, {H5S_UNLIMITED, 20}, {256, 20});
        }
        H5::DataSet identity_dataset = group_cb.openDataSet("local_link_identity");
        H5::DataSpace identity_space = identity_dataset.getSpace();
        hsize_t identity_dims[2] = {0, 0};
        identity_space.getSimpleExtentDims(identity_dims);
        H5::DataSet mechanics_v2_dataset =
            group_cb.openDataSet("local_link_mechanics_v2");
        H5::DataSpace mechanics_v2_space = mechanics_v2_dataset.getSpace();
        hsize_t mechanics_v2_dims[2] = {0, 0};
        mechanics_v2_space.getSimpleExtentDims(mechanics_v2_dims);
        const size_t count = aa_active_link_slots.size();
        std::vector<int> identity;
        std::vector<double> coordinates;
        std::vector<double> mechanics_rows;
        std::vector<double> mechanics_rows_v2;
        std::vector<double> motor_support_rows;
        std::vector<double> arm_state_rows;
        identity.reserve(count * 7);
        coordinates.reserve(count * 2);
        mechanics_rows.reserve(count * 16);
        mechanics_rows_v2.reserve(count * 23);
        motor_support_rows.reserve(count * 20);
        arm_state_rows.reserve(count * 9);
        const double ramp = _current_interaction_ramp();
        const double single_link_alignment_k = ramp *
            actn2_load_kinetics::fixed_single_link_stiffness(
                kappa_aa, aa_stiffness_reference_links);
        std::vector<LocalAALinkResult> saved_mechanics(aa_links.size());
        std::vector<actn2_load_kinetics::ArmLoads> raw_arm_loads(aa_links.size());
        std::vector<double> total_raw_favorable_arm_load(
            static_cast<size_t>(actin.n), 0.0);
        for (size_t slot : aa_active_link_slots) {
            const auto& link = aa_links[slot];
            if (aa_status_mode ==
                actn2_load_kinetics::StatusMode::motor_transmitted_register) {
                const double age = link.candidate_myosin_i >= 0 &&
                        current_step >=
                            static_cast<size_t>(link.candidate_since)
                    ? (current_step -
                       static_cast<size_t>(link.candidate_since)) * dt
                    : 0.0;
                const double challenger_age = link.challenger_myosin_i >= 0 &&
                        current_step >=
                            static_cast<size_t>(link.challenger_since)
                    ? (current_step -
                       static_cast<size_t>(link.challenger_since)) * dt
                    : 0.0;
                motor_support_rows.insert(motor_support_rows.end(), {
                    static_cast<double>(link.candidate_myosin_i),
                    static_cast<double>(link.candidate_myosin_j),
                    static_cast<double>(link.candidate_half_i),
                    static_cast<double>(link.candidate_half_j),
                    link.candidate_support_ema,
                    link.candidate_support_instantaneous, age,
                    link.candidate_locked ? 1.0 : 0.0,
                    static_cast<double>(link.candidate_since),
                    static_cast<double>(link.candidate_last_seen),
                    link.candidate_raw_support_i,
                    link.candidate_raw_support_j,
                    link.candidate_allocated_support_i,
                    link.candidate_allocated_support_j,
                    static_cast<double>(link.challenger_myosin_i),
                    static_cast<double>(link.challenger_myosin_j),
                    static_cast<double>(link.challenger_half_i),
                    static_cast<double>(link.challenger_half_j),
                    link.challenger_support_ema, challenger_age});
            }
            if (aa_arm_state_model &&
                actn2_arm_kinetics::is_tether(link.arm_state)) {
                continue;
            }
            const double link_radial_k =
                ramp * aa_link_radial_stiffness(link.status);
            saved_mechanics[slot] = compute_local_aa_link(
                actin, link.actin_i, link.actin_j, link.s_i, link.s_j, box,
                link_radial_k, ramp * aa_link_perpendicular_stiffness(link.status),
                single_link_alignment_k, aa_optimal, aa_antiparallel_overlap,
                aa_link_register_ratio(link.status),
                aa_register_saturation_error, aa_axial_register_force_max,
            ramp * aa_radial_force_max);
            const auto& mech = saved_mechanics[slot];
            const std::array<double, 3> link_force{
                mech.force_on_first.x, mech.force_on_first.y,
                mech.force_on_first.z};
            const std::array<double, 3> direction_i{
                actin.direction[link.actin_i].x,
                actin.direction[link.actin_i].y,
                actin.direction[link.actin_i].z};
            const std::array<double, 3> direction_j{
                actin.direction[link.actin_j].x,
                actin.direction[link.actin_j].y,
                actin.direction[link.actin_j].z};
            raw_arm_loads[slot] = actn2_load_kinetics::signed_arm_loads(
                link_force, direction_i, direction_j);
            total_raw_favorable_arm_load[static_cast<size_t>(link.actin_i)] +=
                std::max(0.0, raw_arm_loads[slot].first);
            total_raw_favorable_arm_load[static_cast<size_t>(link.actin_j)] +=
                std::max(0.0, raw_arm_loads[slot].second);
        }
        std::vector<double> productive_motor_load(
            static_cast<size_t>(actin.n), 0.0);
        for (int i = 0; i < actin.n; ++i) {
            productive_motor_load[static_cast<size_t>(i)] =
                actn2_load_kinetics::norm(
                    _actin_productive_motor_force_vector(i));
        }
        for (size_t slot : aa_active_link_slots) {
            const auto& link = aa_links[slot];
            identity.insert(identity.end(), {
                static_cast<int>(link.id >> 32),
                static_cast<int>(link.id & 0xffffffffULL),
                link.actin_i, link.actin_j, static_cast<int>(link.status),
                static_cast<int>(link.attach_step),
                static_cast<int>(link.status2_memory_until)});
            coordinates.insert(coordinates.end(), {link.s_i, link.s_j});
            arm_state_rows.insert(arm_state_rows.end(), {
                static_cast<double>(static_cast<int>(link.arm_state)),
                static_cast<double>(link.arm_state_since),
                static_cast<double>(link.last_partner),
                actn2_arm_kinetics::arm_i_bound(link.arm_state) ? 1.0 : 0.0,
                actn2_arm_kinetics::arm_j_bound(link.arm_state) ? 1.0 : 0.0,
                link.arm_force_i, link.arm_force_j, link.contact_force,
                link.contact_penetration});
            if (aa_arm_state_model &&
                actn2_arm_kinetics::is_tether(link.arm_state)) {
                mechanics_rows.insert(mechanics_rows.end(), 16, 0.0);
                mechanics_rows_v2.insert(mechanics_rows_v2.end(), 23, 0.0);
                continue;
            }
            const auto& mech = saved_mechanics[slot];
            const auto raw_arms = raw_arm_loads[slot];
            const bool motor_allocated_mode = aa_status_mode ==
                actn2_load_kinetics::StatusMode::arm_load_productive;
            const bool productive_register_mode = aa_status_mode ==
                actn2_load_kinetics::StatusMode::productive_register;
            const double allocated_motor_i = motor_allocated_mode
                ? actn2_load_kinetics::allocated_motor_load(
                    raw_arms.first,
                    total_raw_favorable_arm_load[
                        static_cast<size_t>(link.actin_i)],
                    productive_motor_load[static_cast<size_t>(link.actin_i)])
                : 0.0;
            const double allocated_motor_j = motor_allocated_mode
                ? actn2_load_kinetics::allocated_motor_load(
                    raw_arms.second,
                    total_raw_favorable_arm_load[
                        static_cast<size_t>(link.actin_j)],
                    productive_motor_load[static_cast<size_t>(link.actin_j)])
                : 0.0;
            const actn2_load_kinetics::ArmLoads arms = motor_allocated_mode
                ? actn2_load_kinetics::ArmLoads{
                      std::min(std::max(0.0, raw_arms.first), allocated_motor_i),
                      std::min(std::max(0.0, raw_arms.second), allocated_motor_j)}
                : raw_arms;
            const double productive_support = _aa_productive_opposed_load(
                link.actin_i, link.actin_j);
            const double bilateral_load = arms.favorable_bilateral();
            const double registry_threshold = link.status == 2
                ? aa_status2_registry_off : aa_status2_registry_on;
            const double arm_load_threshold = link.status == 2
                ? aa_status2_force_off : aa_status2_force_on;
            const bool registry_ok = actn2_load_kinetics::registry_gate(
                static_cast<int>(link.status), mech.direction_dot,
                mech.registry_error, aa_status2_registry_on,
                aa_status2_registry_off);
            const bool arm_load_ok =
                arms.first > arm_load_threshold &&
                arms.second > arm_load_threshold;
            const bool productive_support_ok = productive_support > 0.0;
            const bool status2_gate = productive_register_mode
                ? actn2_load_kinetics::classify_productive_register(
                      static_cast<int>(link.status), productive_support,
                      aa_status2_force_on, mech.direction_dot,
                      mech.registry_error, aa_status2_registry_on,
                      aa_status2_registry_off) == 2
                : motor_allocated_mode
                ? arm_load_ok && productive_support_ok && registry_ok
                : actn2_load_kinetics::classify_arm_loaded_productive(
                      static_cast<int>(link.status), arms, productive_support,
                      aa_status2_force_on, aa_status2_force_off) == 2;
            const double lifetime_load = link.status == 2
                ? (productive_register_mode
                       ? productive_support : bilateral_load)
                : 0.0;
            const double lifetime = catch_bond_kinetics::lifetime_for_status(
                link.status, base_lifetime, directional_base_lifetime,
                lifetime_coeff, lifetime_load, stabilize_status1);
            const double legacy_tension = std::max(
                0.0, ramp * aa_link_radial_stiffness(link.status) *
                    (mech.length - aa_optimal));
            mechanics_rows.insert(mechanics_rows.end(), {
                mech.connector.x, mech.connector.y, mech.connector.z,
                mech.length, mech.length - aa_optimal,
                arms.first, arms.second, bilateral_load,
                productive_support, status2_gate ? 1.0 : 0.0,
                lifetime, legacy_tension,
                raw_arms.first, raw_arms.second,
                allocated_motor_i, allocated_motor_j});
            mechanics_rows_v2.insert(mechanics_rows_v2.end(), {
                mech.connector.x, mech.connector.y, mech.connector.z,
                mech.length, mech.length - aa_optimal,
                arms.first, arms.second, bilateral_load,
                productive_support, status2_gate ? 1.0 : 0.0,
                lifetime, legacy_tension,
                raw_arms.first, raw_arms.second,
                allocated_motor_i, allocated_motor_j,
                mech.direction_dot, mech.registry_error, registry_threshold,
                registry_ok ? 1.0 : 0.0,
                arm_load_ok ? 1.0 : 0.0,
                productive_support_ok ? 1.0 : 0.0,
                status2_gate ? 1.0 : 0.0});
        }
        if (count > 0) {
            append_to_dataset_int(group_cb, "local_link_identity", identity,
                                  {static_cast<hsize_t>(count), 7});
            append_to_dataset(group_cb, "local_link_coordinates", coordinates,
                              {static_cast<hsize_t>(count), 2});
            append_to_dataset(group_cb, "local_link_mechanics", mechanics_rows,
                              {static_cast<hsize_t>(count), 16});
            append_to_dataset(group_cb, "local_link_mechanics_v2",
                              mechanics_rows_v2,
                              {static_cast<hsize_t>(count), 23});
            append_to_dataset(group_cb, "local_link_arm_state",
                              arm_state_rows,
                              {static_cast<hsize_t>(count), 9});
            if (aa_status_mode ==
                actn2_load_kinetics::StatusMode::motor_transmitted_register) {
                append_to_dataset(group_cb, "local_link_motor_support",
                                  motor_support_rows,
                                  {static_cast<hsize_t>(count), 20});
            }
        }
        append_to_dataset_int(
            group_cb, "local_link_frame_offsets",
            {static_cast<int>(current_step), static_cast<int>(identity_dims[0]),
             static_cast<int>(count)}, {1, 3});
        append_to_dataset_int(
            group_cb, "local_link_mechanics_v2_frame_offsets",
            {static_cast<int>(current_step), static_cast<int>(identity_dims[0]),
             static_cast<int>(mechanics_v2_dims[0]), static_cast<int>(count)},
            {1, 4});
        file.flush(H5F_SCOPE_GLOBAL);
    }

    // Sparse frame-ragged myomesin occupancy.  Empty frames still receive an
    // offset row, so occupancy is unambiguous and arbitrary-frame resume can
    // restore the exact molecular graph.
    if (k_mm > 0.0) {
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/myomesin")) file.createGroup("/myomesin");
        if (!file.nameExists("/myomesin/bridge_state")) {
            create_empty_dataset(
                file, "/myomesin", "bridge_state",
                {0, 8}, {H5S_UNLIMITED, 8}, {128, 8});
            create_empty_dataset_int(
                file, "/myomesin", "frame_offsets",
                {0, 3}, {H5S_UNLIMITED, 3}, {128, 3});
        }
        if (!file.nameExists("/myomesin/next_id_by_frame")) {
            create_empty_dataset_int(
                file, "/myomesin", "next_id_by_frame",
                {0, 2}, {H5S_UNLIMITED, 2}, {128, 2});
        }
        H5::Group group(file.openGroup("/myomesin"));
        H5::DataSet state_dataset = group.openDataSet("bridge_state");
        H5::DataSpace state_space = state_dataset.getSpace();
        hsize_t state_dims[2] = {0, 0};
        state_space.getSimpleExtentDims(state_dims, nullptr);
        std::vector<double> rows;
        rows.reserve(myomesin_bridges.size() * 8);
        for (const auto& bridge : myomesin_bridges) {
            rows.insert(rows.end(), {
                static_cast<double>(bridge.id >> 32),
                static_cast<double>(bridge.id & 0xffffffffULL),
                static_cast<double>(bridge.myosin_i),
                static_cast<double>(bridge.myosin_j),
                static_cast<double>(bridge.attach_step),
                bridge.axial_offset, bridge.radial_distance, bridge.force});
        }
        if (!myomesin_bridges.empty()) {
            append_to_dataset(
                group, "bridge_state", rows,
                {static_cast<hsize_t>(myomesin_bridges.size()), 8});
        }
        append_to_dataset_int(
            group, "frame_offsets",
            {static_cast<int>(current_step), static_cast<int>(state_dims[0]),
             static_cast<int>(myomesin_bridges.size())},
            {1, 3});
        append_to_dataset_int(
            group, "next_id_by_frame",
            {static_cast<int>(myomesin_next_id >> 32),
             static_cast<int>(myomesin_next_id & 0xffffffffULL)},
            {1, 2});
    }

    // Flush any recorded catch-bond events to the HDF5 file
    if (!cb_breakage_events.empty() || !cb_limit_events.empty() || !aa_completed_lifetimes.empty()) {
        H5::H5File file(filename, H5F_ACC_RDWR);
        H5::Group group_cb(file.openGroup("/catch_bond"));

        if (!cb_breakage_events.empty()) {
            hsize_t event_width = 11 + 2 * max_myosin_bonds;
            hsize_t n_events = cb_breakage_events.size() / event_width;
            append_to_dataset(group_cb, "breakage", cb_breakage_events,
                               {n_events, event_width});
            cb_breakage_events.clear();
        }

        if (!cb_limit_events.empty()) {
            hsize_t limit_width = 5;
            hsize_t n_limit = cb_limit_events.size() / limit_width;
            append_to_dataset(group_cb, "limit_removal", cb_limit_events,
                               {n_limit, limit_width});
            cb_limit_events.clear();
        }

        if (!aa_completed_lifetimes.empty()) {
            hsize_t n = aa_completed_lifetimes.size();
            append_to_dataset(group_cb, "completed_lifetimes", aa_completed_lifetimes, {n});
            aa_completed_lifetimes.clear();
        }
    }

    if (!am_transition_events.empty()) {
        constexpr hsize_t transition_width = 11;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/actin_myo/head_transitions")) {
            create_empty_dataset(
                file, "/actin_myo", "head_transitions",
                {0, transition_width}, {H5S_UNLIMITED, transition_width},
                {100, transition_width});
        }
        H5::Group group_am(file.openGroup("/actin_myo"));
        const hsize_t n_events = am_transition_events.size() / transition_width;
        append_to_dataset(group_am, "head_transitions", am_transition_events,
                          {n_events, transition_width});
        am_transition_events.clear();
    }

    if (!am_forced_rupture_events.empty()) {
        constexpr hsize_t rupture_width = 15;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/actin_myo/forced_ruptures")) {
            create_empty_dataset(
                file, "/actin_myo", "forced_ruptures",
                {0, rupture_width}, {H5S_UNLIMITED, rupture_width},
                {100, rupture_width});
        }
        H5::Group group_am(file.openGroup("/actin_myo"));
        const hsize_t n_events =
            am_forced_rupture_events.size() / rupture_width;
        append_to_dataset(group_am, "forced_ruptures",
                          am_forced_rupture_events,
                          {n_events, rupture_width});
        am_forced_rupture_events.clear();
    }
    if (!am_forced_rupture_events_v2.empty()) {
        constexpr hsize_t rupture_width = 19;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/actin_myo/forced_ruptures_v2")) {
            create_empty_dataset(
                file, "/actin_myo", "forced_ruptures_v2",
                {0, rupture_width}, {H5S_UNLIMITED, rupture_width},
                {100, rupture_width});
        }
        H5::Group group_am(file.openGroup("/actin_myo"));
        const hsize_t n_events =
            am_forced_rupture_events_v2.size() / rupture_width;
        append_to_dataset(group_am, "forced_ruptures_v2",
                          am_forced_rupture_events_v2,
                          {n_events, rupture_width});
        am_forced_rupture_events_v2.clear();
    }

    if (!am_thick_transition_events.empty()) {
        constexpr hsize_t transition_width = 7;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/actin_myo/thick_transitions")) {
            create_empty_dataset(
                file, "/actin_myo", "thick_transitions",
                {0, transition_width}, {H5S_UNLIMITED, transition_width},
                {100, transition_width});
        }
        H5::Group group_am(file.openGroup("/actin_myo"));
        const hsize_t n_events =
            am_thick_transition_events.size() / transition_width;
        append_to_dataset(group_am, "thick_transitions",
                          am_thick_transition_events,
                          {n_events, transition_width});
        am_thick_transition_events.clear();
    }

    if (!actin_turnover_events.empty()) {
        constexpr hsize_t event_width = 20;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/actin/turnover_events")) {
            create_empty_dataset(
                file, "/actin", "turnover_events",
                {0, event_width}, {H5S_UNLIMITED, event_width},
                {100, event_width});
        }
        H5::Group group_actin(file.openGroup("/actin"));
        const hsize_t n_events = actin_turnover_events.size() / event_width;
        append_to_dataset(group_actin, "turnover_events",
                          actin_turnover_events, {n_events, event_width});
        actin_turnover_events.clear();
    }

    if (!aa_state_events.empty()) {
        constexpr hsize_t state_width = 13;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/catch_bond/state_transitions")) {
            create_empty_dataset(
                file, "/catch_bond", "state_transitions",
                {0, state_width}, {H5S_UNLIMITED, state_width},
                {100, state_width});
        }
        H5::Group group_cb(file.openGroup("/catch_bond"));
        const hsize_t n_events = aa_state_events.size() / state_width;
        append_to_dataset(group_cb, "state_transitions", aa_state_events,
                          {n_events, state_width});
        aa_state_events.clear();
    }

    if (!aa_load_events.empty()) {
        constexpr hsize_t load_width = 11;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/catch_bond/load_transitions")) {
            create_empty_dataset(
                file, "/catch_bond", "load_transitions",
                {0, load_width}, {H5S_UNLIMITED, load_width},
                {100, load_width});
        }
        H5::Group group_cb(file.openGroup("/catch_bond"));
        const hsize_t n_events = aa_load_events.size() / load_width;
        append_to_dataset(group_cb, "load_transitions", aa_load_events,
                          {n_events, load_width});
        aa_load_events.clear();
    }

    if (!aa_local_link_events.empty()) {
        constexpr hsize_t event_width = 11;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/catch_bond/local_link_events")) {
            create_empty_dataset(
                file, "/catch_bond", "local_link_events",
                {0, event_width}, {H5S_UNLIMITED, event_width},
                {100, event_width});
        }
        H5::Group group_cb(file.openGroup("/catch_bond"));
        const hsize_t n_events = aa_local_link_events.size() / event_width;
        append_to_dataset(group_cb, "local_link_events",
                          aa_local_link_events, {n_events, event_width});
        aa_local_link_events.clear();
    }
    if (!aa_local_link_events_v2.empty()) {
        constexpr hsize_t event_width = 21;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/catch_bond/local_link_events_v2")) {
            create_empty_dataset(
                file, "/catch_bond", "local_link_events_v2",
                {0, event_width}, {H5S_UNLIMITED, event_width},
                {100, event_width});
        }
        H5::Group group_cb(file.openGroup("/catch_bond"));
        const hsize_t n_events = aa_local_link_events_v2.size() / event_width;
        append_to_dataset(group_cb, "local_link_events_v2",
                          aa_local_link_events_v2, {n_events, event_width});
        aa_local_link_events_v2.clear();
    }
    if (!aa_local_link_events_v3.empty()) {
        constexpr hsize_t event_width = 35;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/catch_bond/local_link_events_v3")) {
            create_empty_dataset(
                file, "/catch_bond", "local_link_events_v3",
                {0, event_width}, {H5S_UNLIMITED, event_width},
                {100, event_width});
        }
        H5::Group group_cb(file.openGroup("/catch_bond"));
        const hsize_t n_events = aa_local_link_events_v3.size() / event_width;
        append_to_dataset(group_cb, "local_link_events_v3",
                          aa_local_link_events_v3, {n_events, event_width});
        aa_local_link_events_v3.clear();
    }
    if (!aa_tip_accessibility_events.empty()) {
        constexpr hsize_t event_width = 13;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/catch_bond/tip_accessibility_events")) {
            create_empty_dataset(
                file, "/catch_bond", "tip_accessibility_events",
                {0, event_width}, {H5S_UNLIMITED, event_width},
                {100, event_width});
        }
        H5::Group group_cb(file.openGroup("/catch_bond"));
        const hsize_t n_events =
            aa_tip_accessibility_events.size() / event_width;
        append_to_dataset(group_cb, "tip_accessibility_events",
                          aa_tip_accessibility_events,
                          {n_events, event_width});
        aa_tip_accessibility_events.clear();
    }
    if (!aa_arm_events.empty()) {
        constexpr hsize_t event_width = 14;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/catch_bond/arm_events")) {
            create_empty_dataset(
                file, "/catch_bond", "arm_events",
                {0, event_width}, {H5S_UNLIMITED, event_width},
                {100, event_width});
        }
        H5::Group group_cb(file.openGroup("/catch_bond"));
        const hsize_t n_events = aa_arm_events.size() / event_width;
        append_to_dataset(group_cb, "arm_events", aa_arm_events,
                          {n_events, event_width});
        aa_arm_events.clear();
    }
    if (!aa_motor_pair_events.empty()) {
        constexpr hsize_t event_width = 20;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/catch_bond/motor_pair_events")) {
            create_empty_dataset(
                file, "/catch_bond", "motor_pair_events",
                {0, event_width}, {H5S_UNLIMITED, event_width},
                {100, event_width});
        }
        H5::Group group_cb(file.openGroup("/catch_bond"));
        const hsize_t n_events = aa_motor_pair_events.size() / event_width;
        append_to_dataset(group_cb, "motor_pair_events",
                          aa_motor_pair_events, {n_events, event_width});
        aa_motor_pair_events.clear();
    }
    if (!explicit_titin_events.empty()) {
        constexpr hsize_t event_width = 10;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/titin")) file.createGroup("/titin");
        if (!file.nameExists("/titin/events")) {
            create_empty_dataset(
                file, "/titin", "events",
                {0, event_width}, {H5S_UNLIMITED, event_width},
                {100, event_width});
        }
        H5::Group group_titin(file.openGroup("/titin"));
        const hsize_t n_events = explicit_titin_events.size() / event_width;
        append_to_dataset(group_titin, "events", explicit_titin_events,
                          {n_events, event_width});
        explicit_titin_events.clear();
    }
    if (!myomesin_events.empty()) {
        constexpr hsize_t event_width = 9;
        H5::H5File file(filename, H5F_ACC_RDWR);
        if (!file.nameExists("/myomesin")) file.createGroup("/myomesin");
        if (!file.nameExists("/myomesin/events")) {
            create_empty_dataset(
                file, "/myomesin", "events",
                {0, event_width}, {H5S_UNLIMITED, event_width},
                {100, event_width});
        }
        H5::Group group(file.openGroup("/myomesin"));
        append_to_dataset(
            group, "events", myomesin_events,
            {static_cast<hsize_t>(myomesin_events.size() / event_width),
             event_width});
        myomesin_events.clear();
    }
}

void Sarcomere::save_resume_snapshot() {
    H5::H5File file(filename, H5F_ACC_RDWR);
    const std::string final_resume_path = "/resume";
    const std::string temp_resume_path = "/resume_new";
    const std::string backup_resume_path = "/resume_prev";
    delete_link_if_exists(file, temp_resume_path);
    H5::Group group_resume = file.createGroup(temp_resume_path);

    const hsize_t n_actins = static_cast<hsize_t>(actin.n);
    const hsize_t n_myosins = static_cast<hsize_t>(myosin.n);
    const hsize_t max_aa_pairs = static_cast<hsize_t>(actin.n) *
                                 static_cast<hsize_t>(std::max(10, 2 * max_myosin_bonds));
    const hsize_t max_am_pairs = static_cast<hsize_t>(myosin.n) *
                                 static_cast<hsize_t>(max_myosin_bonds);
    const hsize_t max_recovery_pairs = max_aa_pairs;

    auto replace_dataset_int = [&](const std::string& name,
                                   const std::vector<hsize_t>& dims,
                                   const std::vector<int>& data) {
        const std::string full_path = temp_resume_path + "/" + name;
        if (file.nameExists(full_path)) {
            H5Ldelete(file.getId(), full_path.c_str(), H5P_DEFAULT);
        }
        H5::DataSpace dataspace(dims.size(), dims.data());
        H5::IntType datatype(H5::PredType::STD_I32LE);
        H5::DataSet dataset = group_resume.createDataSet(name, datatype, dataspace);
        dataset.write(data.data(), H5::PredType::STD_I32LE);
    };

    auto replace_dataset_double = [&](const std::string& name,
                                      const std::vector<hsize_t>& dims,
                                      const std::vector<double>& data) {
        const std::string full_path = temp_resume_path + "/" + name;
        if (file.nameExists(full_path)) {
            H5Ldelete(file.getId(), full_path.c_str(), H5P_DEFAULT);
        }
        H5::DataSpace dataspace(dims.size(), dims.data());
        H5::FloatType datatype(H5::PredType::IEEE_F64LE);
        H5::DataSet dataset = group_resume.createDataSet(name, datatype, dataspace);
        dataset.write(data.data(), H5::PredType::IEEE_F64LE);
    };

    auto aa_prev = pack_aa_bond_state(
        previous_aa_pairs, actin_actin_bonds_prev,
        actin_actin_status_prev, actin_actin_lifetime_prev);
    auto aa_current = pack_aa_bond_state(
        active_aa_pairs, actin_actin_bonds,
        actin_actin_status, actin_actin_lifetime);
    auto am_prev = pack_am_bond_pairs(previous_am_pairs, am_bonds_prev);
    auto am_current = pack_am_bond_pairs(active_am_pairs, am_bonds);
    auto aa_recovery = pack_recovery_state(
        active_recovery_pairs, actin_recovery_until, current_step);

    if (aa_recovery.until.size() > static_cast<size_t>(max_recovery_pairs)) {
        // Degrade gracefully (truncate) instead of aborting; see pad_int_vector / save_state.
        static bool warned_resume_capacity = false;
        if (!warned_resume_capacity) {
            fprintf(stderr, "Warning: recovery pair list (%zu) exceeds compact resume capacity (%llu); "
                    "truncating resume snapshot.\n",
                    aa_recovery.until.size(), static_cast<unsigned long long>(max_recovery_pairs));
            warned_resume_capacity = true;
        }
    }

    std::vector<int> aa_prev_pairs = aa_prev.pairs;
    std::vector<int> aa_prev_status = aa_prev.status;
    std::vector<int> aa_prev_lifetime = aa_prev.lifetime;
    std::vector<int> aa_current_pairs = aa_current.pairs;
    std::vector<int> aa_current_status = aa_current.status;
    std::vector<int> aa_current_lifetime = aa_current.lifetime;
    std::vector<int> aa_current_attach_step;
    std::vector<int> aa_current_status2_memory_until;
    aa_current_attach_step.reserve(aa_current.status.size());
    aa_current_status2_memory_until.reserve(aa_current.status.size());
    for (size_t idx = 0; idx + 1 < aa_current.pairs.size(); idx += 2) {
        const int a = aa_current.pairs[idx];
        const int b = aa_current.pairs[idx + 1];
        aa_current_attach_step.push_back(aa_attach_step[a][b]);
        aa_current_status2_memory_until.push_back(
            static_cast<int>(aa_status2_memory_until[a][b]));
    }
    std::vector<int> am_prev_pairs = am_prev;
    std::vector<int> am_current_pairs = am_current;
    std::vector<int> am_current_state;
    std::vector<double> am_current_target;
    std::vector<int> am_head_state_current;
    std::vector<int> am_head_bound_actin_current;
    std::vector<double> am_head_target_current;
    std::vector<double> am_head_binding_azimuth_current;
    std::vector<double> am_head_stroke_progress_current;
    am_current_state.reserve(am_current.size() / 2);
    am_current_target.reserve(am_current.size() / 2);
    for (size_t idx = 0; idx + 1 < am_current.size(); idx += 2) {
        const int a = am_current[idx];
        const int m = am_current[idx + 1];
        am_current_state.push_back(dynamic_am ? static_cast<int>(am_state[a][m]) : 2);
        am_current_target.push_back(dynamic_am ? am_target_axial[a][m]
                                               : _am_axial_coordinate(a, m));
    }
    const size_t total_head_count = static_cast<size_t>(myosin.n) *
                                    static_cast<size_t>(2 * am_heads_per_half);
    am_head_state_current.reserve(total_head_count);
    am_head_bound_actin_current.reserve(total_head_count);
    am_head_target_current.reserve(total_head_count);
    am_head_binding_azimuth_current.reserve(total_head_count);
    am_head_stroke_progress_current.reserve(total_head_count);
    for (const auto& myosin_halves : am_heads_by_myosin) {
        for (const auto& half_heads : myosin_halves) {
            for (const auto& head : half_heads) {
                am_head_state_current.push_back(static_cast<int>(head.state));
                am_head_bound_actin_current.push_back(head.bound_actin);
                am_head_target_current.push_back(head.target_axial);
                am_head_binding_azimuth_current.push_back(
                    head.binding_azimuth);
                am_head_stroke_progress_current.push_back(
                    head.stroke_progress);
            }
        }
    }
    std::vector<int> aa_recovery_pairs = aa_recovery.pairs;
    std::vector<int> aa_recovery_until_values = aa_recovery.until;
    const hsize_t local_link_rows = static_cast<hsize_t>(
        std::max<size_t>(1, aa_active_link_slots.size()));
    std::vector<int> local_link_identity(
        static_cast<size_t>(local_link_rows) * 7, -1);
    std::vector<double> local_link_coordinates(
        static_cast<size_t>(local_link_rows) * 2, 0.0);
    std::vector<int> local_link_motor_identity(
        static_cast<size_t>(local_link_rows) * 12, -1);
    std::vector<double> local_link_motor_evidence(
        static_cast<size_t>(local_link_rows) * 7, 0.0);
    std::vector<int> local_link_arm_state(
        static_cast<size_t>(local_link_rows) * 3, -1);
    std::vector<double> local_link_arm_loads(
        static_cast<size_t>(local_link_rows) * 4, 0.0);
    for (size_t record = 0; record < aa_active_link_slots.size(); ++record) {
        const auto& link = aa_links[aa_active_link_slots[record]];
        local_link_identity[7 * record] = static_cast<int>(link.id >> 32);
        local_link_identity[7 * record + 1] =
            static_cast<int>(link.id & 0xffffffffULL);
        local_link_identity[7 * record + 2] = link.actin_i;
        local_link_identity[7 * record + 3] = link.actin_j;
        local_link_identity[7 * record + 4] = static_cast<int>(link.status);
        local_link_identity[7 * record + 5] = static_cast<int>(link.attach_step);
        local_link_identity[7 * record + 6] =
            static_cast<int>(link.status2_memory_until);
        local_link_coordinates[2 * record] = link.s_i;
        local_link_coordinates[2 * record + 1] = link.s_j;
        local_link_motor_identity[12 * record] = link.candidate_myosin_i;
        local_link_motor_identity[12 * record + 1] = link.candidate_myosin_j;
        local_link_motor_identity[12 * record + 2] = link.candidate_half_i;
        local_link_motor_identity[12 * record + 3] = link.candidate_half_j;
        local_link_motor_identity[12 * record + 4] =
            static_cast<int>(link.candidate_since);
        local_link_motor_identity[12 * record + 5] =
            static_cast<int>(link.candidate_last_seen);
        local_link_motor_identity[12 * record + 6] =
            link.candidate_locked ? 1 : 0;
        local_link_motor_identity[12 * record + 7] =
            link.challenger_myosin_i;
        local_link_motor_identity[12 * record + 8] =
            link.challenger_myosin_j;
        local_link_motor_identity[12 * record + 9] =
            link.challenger_half_i;
        local_link_motor_identity[12 * record + 10] =
            link.challenger_half_j;
        local_link_motor_identity[12 * record + 11] =
            static_cast<int>(link.challenger_since);
        local_link_motor_evidence[7 * record] = link.candidate_support_ema;
        local_link_motor_evidence[7 * record + 1] =
            link.candidate_support_instantaneous;
        local_link_motor_evidence[7 * record + 2] =
            link.candidate_raw_support_i;
        local_link_motor_evidence[7 * record + 3] =
            link.candidate_raw_support_j;
        local_link_motor_evidence[7 * record + 4] =
            link.candidate_allocated_support_i;
        local_link_motor_evidence[7 * record + 5] =
            link.candidate_allocated_support_j;
        local_link_motor_evidence[7 * record + 6] =
            link.challenger_support_ema;
        local_link_arm_state[3 * record] =
            static_cast<int>(link.arm_state);
        local_link_arm_state[3 * record + 1] =
            static_cast<int>(link.arm_state_since);
        local_link_arm_state[3 * record + 2] = link.last_partner;
        local_link_arm_loads[4 * record] = link.arm_force_i;
        local_link_arm_loads[4 * record + 1] = link.arm_force_j;
        local_link_arm_loads[4 * record + 2] = link.contact_force;
        local_link_arm_loads[4 * record + 3] = link.contact_penetration;
    }
    const hsize_t titin_rows = static_cast<hsize_t>(explicit_titins.size());
    std::vector<int> titin_state_int;
    std::vector<double> titin_state_double;
    if (titin_mode == TitinMode::explicit_dynamic) {
        titin_state_int.reserve(explicit_titins.size() * 5);
        titin_state_double.reserve(explicit_titins.size() * 6);
        for (const auto& titin : explicit_titins) {
            titin_state_int.insert(titin_state_int.end(), {
                titin.bound ? 1 : 0,
                titin.bound ? titin.bound_actin : -1,
                titin.last_bound_actin,
                static_cast<int>(titin.attach_step),
                static_cast<int>(titin.last_detach_step)});
            titin_state_double.insert(titin_state_double.end(), {
                titin.distance, titin.extension, titin.force,
                titin.attachment_distance, titin.attachment_energy,
                titin.attachment_force});
        }
    }
    const hsize_t myomesin_rows = static_cast<hsize_t>(
        std::max<size_t>(1, static_cast<size_t>(myosin.n) *
                               static_cast<size_t>(myomesin_sites_per_filament) /
                               2));
    std::vector<int> myomesin_state_int(
        static_cast<size_t>(myomesin_rows) * 5, -1);
    std::vector<double> myomesin_state_double(
        static_cast<size_t>(myomesin_rows) * 3, 0.0);
    for (size_t index = 0; index < myomesin_bridges.size(); ++index) {
        const auto& bridge = myomesin_bridges[index];
        myomesin_state_int[5 * index] =
            static_cast<int>(bridge.id >> 32);
        myomesin_state_int[5 * index + 1] =
            static_cast<int>(bridge.id & 0xffffffffULL);
        myomesin_state_int[5 * index + 2] = bridge.myosin_i;
        myomesin_state_int[5 * index + 3] = bridge.myosin_j;
        myomesin_state_int[5 * index + 4] =
            static_cast<int>(bridge.attach_step);
        myomesin_state_double[3 * index] = bridge.axial_offset;
        myomesin_state_double[3 * index + 1] = bridge.radial_distance;
        myomesin_state_double[3 * index + 2] = bridge.force;
    }

    pad_int_vector(aa_prev_pairs, static_cast<size_t>(max_aa_pairs * 2), -1);
    pad_int_vector(aa_prev_status, static_cast<size_t>(max_aa_pairs), 0);
    pad_int_vector(aa_prev_lifetime, static_cast<size_t>(max_aa_pairs), 0);
    pad_int_vector(aa_current_pairs, static_cast<size_t>(max_aa_pairs * 2), -1);
    pad_int_vector(aa_current_status, static_cast<size_t>(max_aa_pairs), 0);
    pad_int_vector(aa_current_lifetime, static_cast<size_t>(max_aa_pairs), 0);
    pad_int_vector(aa_current_attach_step, static_cast<size_t>(max_aa_pairs), -1);
    pad_int_vector(aa_current_status2_memory_until,
                   static_cast<size_t>(max_aa_pairs), 0);
    pad_int_vector(am_prev_pairs, static_cast<size_t>(max_am_pairs * 2), -1);
    pad_int_vector(am_current_pairs, static_cast<size_t>(max_am_pairs * 2), -1);
    pad_int_vector(am_current_state, static_cast<size_t>(max_am_pairs), 0);
    am_current_target.resize(static_cast<size_t>(max_am_pairs), 0.0);
    pad_int_vector(am_head_state_current, total_head_count, 0);
    pad_int_vector(am_head_bound_actin_current, total_head_count, -1);
    am_head_target_current.resize(total_head_count, 0.0);
    am_head_binding_azimuth_current.resize(total_head_count, 0.0);
    am_head_stroke_progress_current.resize(total_head_count, 0.0);
    pad_int_vector(aa_recovery_pairs, static_cast<size_t>(max_recovery_pairs * 2), -1);
    pad_int_vector(aa_recovery_until_values, static_cast<size_t>(max_recovery_pairs), 0);

    replace_dataset_int("current_step", {1}, {static_cast<int>(current_step)});
    replace_dataset_double("actin_center", {n_actins, 3}, flatten_3d_array(actin.center));
    replace_dataset_double("actin_direction", {n_actins, 3}, flatten_3d_array(actin.direction));
    replace_dataset_double("myosin_center", {n_myosins, 3}, flatten_3d_array(myosin.center));
    replace_dataset_double("myosin_direction", {n_myosins, 3}, flatten_3d_array(myosin.direction));
    replace_dataset_double("myosin_material_e1_current", {n_myosins, 3},
                           flatten_vec_values(myosin_material_e1));
    replace_dataset_double("myosin_material_e2_current", {n_myosins, 3},
                           flatten_vec_values(myosin_material_e2));
    replace_dataset_double("am_q_current", {n_actins}, am_thin_activation);
    replace_dataset_int(
        "actin_last_turnover_step_current", {n_actins},
        std::vector<int>(actin_last_turnover_step.begin(),
                         actin_last_turnover_step.end()));
    replace_dataset_int("aa_pairs_prev_count", {1}, {static_cast<int>(std::min(aa_prev.status.size(), static_cast<size_t>(max_aa_pairs)))});
    replace_dataset_int("aa_pairs_current_count", {1}, {static_cast<int>(std::min(aa_current.status.size(), static_cast<size_t>(max_aa_pairs)))});
    replace_dataset_int("am_pairs_prev_count", {1}, {static_cast<int>(std::min(am_prev.size() / 2, static_cast<size_t>(max_am_pairs)))});
    replace_dataset_int("am_pairs_current_count", {1}, {static_cast<int>(std::min(am_current.size() / 2, static_cast<size_t>(max_am_pairs)))});
    replace_dataset_int("aa_pairs_prev", {max_aa_pairs, 2}, aa_prev_pairs);
    replace_dataset_int("aa_status_prev", {max_aa_pairs, 1}, aa_prev_status);
    replace_dataset_int("aa_lifetime_prev", {max_aa_pairs, 1}, aa_prev_lifetime);
    replace_dataset_int("aa_pairs_current", {max_aa_pairs, 2}, aa_current_pairs);
    replace_dataset_int("aa_status_current", {max_aa_pairs, 1}, aa_current_status);
    replace_dataset_int("aa_lifetime_current", {max_aa_pairs, 1}, aa_current_lifetime);
    replace_dataset_int("aa_attach_step_current", {max_aa_pairs, 1}, aa_current_attach_step);
    replace_dataset_int("aa_status2_memory_until_current", {max_aa_pairs, 1},
                        aa_current_status2_memory_until);
    replace_dataset_int("am_pairs_prev", {max_am_pairs, 2}, am_prev_pairs);
    replace_dataset_int("am_pairs_current", {max_am_pairs, 2}, am_current_pairs);
    replace_dataset_int("am_state_current", {max_am_pairs}, am_current_state);
    replace_dataset_double("am_target_axial_current", {max_am_pairs}, am_current_target);
    const hsize_t heads_per_myosin = static_cast<hsize_t>(2 * am_heads_per_half);
    replace_dataset_int("am_head_state_current",
                        {n_myosins, heads_per_myosin},
                        am_head_state_current);
    replace_dataset_int("am_head_bound_actin_current",
                        {n_myosins, heads_per_myosin},
                        am_head_bound_actin_current);
    replace_dataset_double("am_head_target_axial_current",
                           {n_myosins, heads_per_myosin},
                           am_head_target_current);
    replace_dataset_double("am_head_binding_azimuth_current",
                           {n_myosins, heads_per_myosin},
                           am_head_binding_azimuth_current);
    replace_dataset_double("am_head_stroke_progress_current",
                           {n_myosins, heads_per_myosin},
                           am_head_stroke_progress_current);
    replace_dataset_int("aa_recovery_count", {1}, {static_cast<int>(std::min(aa_recovery.until.size(), static_cast<size_t>(max_recovery_pairs)))});
    replace_dataset_int("aa_recovery_pairs", {max_recovery_pairs, 2}, aa_recovery_pairs);
    replace_dataset_int("aa_recovery_until_values", {max_recovery_pairs, 1}, aa_recovery_until_values);
    replace_dataset_int("aa_local_link_count", {1},
                        {static_cast<int>(aa_active_link_slots.size())});
    replace_dataset_int("aa_local_link_identity", {local_link_rows, 7},
                        local_link_identity);
    replace_dataset_double("aa_local_link_coordinates", {local_link_rows, 2},
                           local_link_coordinates);
    replace_dataset_int("aa_local_link_motor_identity", {local_link_rows, 12},
                        local_link_motor_identity);
    replace_dataset_double("aa_local_link_motor_evidence", {local_link_rows, 7},
                           local_link_motor_evidence);
    replace_dataset_int("aa_local_link_arm_state", {local_link_rows, 3},
                        local_link_arm_state);
    replace_dataset_double("aa_local_link_arm_loads", {local_link_rows, 4},
                           local_link_arm_loads);
    replace_dataset_int("aa_local_next_id", {2},
                        {static_cast<int>(aa_next_link_id >> 32),
                         static_cast<int>(aa_next_link_id & 0xffffffffULL)});
    if (titin_mode == TitinMode::explicit_dynamic) {
        replace_dataset_int("titin_state_int", {titin_rows, 5},
                            titin_state_int);
        replace_dataset_double("titin_state_double", {titin_rows, 6},
                               titin_state_double);
        // Carried across resume so the per-frame provenance series stays a
        // monotone cumulative count instead of restarting at zero.
        replace_dataset_int(
            "titin_association_counts", {1, 4},
            {static_cast<int>(explicit_titin_supported_associations),
             static_cast<int>(explicit_titin_unsupported_associations),
             static_cast<int>(explicit_titin_rebind_associations),
             static_cast<int>(explicit_titin_unsupported_releases)});
    }
    replace_dataset_int(
        "myomesin_bridge_count", {1},
        {static_cast<int>(myomesin_bridges.size())});
    replace_dataset_int("myomesin_state_int", {myomesin_rows, 5},
                        myomesin_state_int);
    replace_dataset_double("myomesin_state_double", {myomesin_rows, 3},
                           myomesin_state_double);
    replace_dataset_int(
        "myomesin_next_id", {2},
        {static_cast<int>(myomesin_next_id >> 32),
         static_cast<int>(myomesin_next_id & 0xffffffffULL)});

    std::vector<double> neighbor_last_actin_x;
    std::vector<double> neighbor_last_actin_y;
    std::vector<double> neighbor_last_actin_z;
    std::vector<double> neighbor_last_myosin_x;
    std::vector<double> neighbor_last_myosin_y;
    std::vector<double> neighbor_last_myosin_z;
    neighbor_list.get_last_species_positions(
        neighbor_last_actin_x, neighbor_last_actin_y, neighbor_last_actin_z,
        neighbor_last_myosin_x, neighbor_last_myosin_y, neighbor_last_myosin_z);
    replace_dataset_double("neighbor_last_actin_x", {n_actins}, neighbor_last_actin_x);
    replace_dataset_double("neighbor_last_actin_y", {n_actins}, neighbor_last_actin_y);
    replace_dataset_double("neighbor_last_actin_z", {n_actins}, neighbor_last_actin_z);
    replace_dataset_double("neighbor_last_myosin_x", {n_myosins}, neighbor_last_myosin_x);
    replace_dataset_double("neighbor_last_myosin_y", {n_myosins}, neighbor_last_myosin_y);
    replace_dataset_double("neighbor_last_myosin_z", {n_myosins}, neighbor_last_myosin_z);
    replace_dataset_double("cb_breakage_pending",
                           {static_cast<hsize_t>(cb_breakage_events.size())},
                           cb_breakage_events);
    replace_dataset_double("cb_limit_pending",
                           {static_cast<hsize_t>(cb_limit_events.size())},
                           cb_limit_events);
    replace_dataset_double("aa_completed_lifetimes_pending",
                           {static_cast<hsize_t>(aa_completed_lifetimes.size())},
                           aa_completed_lifetimes);
    replace_dataset_double("am_transition_events_pending",
                           {static_cast<hsize_t>(am_transition_events.size())},
                           am_transition_events);
    replace_dataset_double("am_forced_rupture_events_pending",
                           {static_cast<hsize_t>(am_forced_rupture_events.size())},
                           am_forced_rupture_events);
    replace_dataset_double("am_forced_rupture_events_v2_pending",
                           {static_cast<hsize_t>(am_forced_rupture_events_v2.size())},
                           am_forced_rupture_events_v2);
    replace_dataset_double(
        "am_thick_transition_events_pending",
        {static_cast<hsize_t>(am_thick_transition_events.size())},
        am_thick_transition_events);
    replace_dataset_double(
        "actin_turnover_events_pending",
        {static_cast<hsize_t>(actin_turnover_events.size())},
        actin_turnover_events);
    replace_dataset_double("aa_state_events_pending",
                           {static_cast<hsize_t>(aa_state_events.size())},
                           aa_state_events);
    replace_dataset_double("aa_load_events_pending",
                           {static_cast<hsize_t>(aa_load_events.size())},
                           aa_load_events);
    replace_dataset_double("aa_local_link_events_pending",
                           {static_cast<hsize_t>(aa_local_link_events.size())},
                           aa_local_link_events);
    replace_dataset_double("aa_local_link_events_v2_pending",
                           {static_cast<hsize_t>(aa_local_link_events_v2.size())},
                           aa_local_link_events_v2);
    replace_dataset_double("aa_local_link_events_v3_pending",
                           {static_cast<hsize_t>(aa_local_link_events_v3.size())},
                           aa_local_link_events_v3);
    replace_dataset_double(
        "aa_tip_accessibility_events_pending",
        {static_cast<hsize_t>(aa_tip_accessibility_events.size())},
        aa_tip_accessibility_events);
    replace_dataset_double("aa_motor_pair_events_pending",
                           {static_cast<hsize_t>(aa_motor_pair_events.size())},
                           aa_motor_pair_events);
    replace_dataset_double("aa_arm_events_pending",
                           {static_cast<hsize_t>(aa_arm_events.size())},
                           aa_arm_events);
    replace_dataset_double("explicit_titin_events_pending",
                           {static_cast<hsize_t>(explicit_titin_events.size())},
                           explicit_titin_events);
    replace_dataset_double("myomesin_events_pending",
                           {static_cast<hsize_t>(myomesin_events.size())},
                           myomesin_events);

    if (rng != nullptr) {
        const size_t main_state_size = gsl_rng_size(rng);
        std::vector<int> main_state(main_state_size, 0);
        const auto* state_ptr = static_cast<const unsigned char*>(gsl_rng_state(rng));
        for (size_t idx = 0; idx < main_state_size; ++idx) {
            main_state[idx] = static_cast<int>(state_ptr[idx]);
        }
        replace_dataset_int("rng_main_state", {static_cast<hsize_t>(main_state_size)}, main_state);
    }
    if (!rng_engines.empty() && rng_engines[0] != nullptr) {
        const int thread_count = static_cast<int>(rng_engines.size());
        const size_t thread_state_size = gsl_rng_size(rng_engines[0]);
        const hsize_t flat_width = static_cast<hsize_t>(thread_count) *
                                   static_cast<hsize_t>(thread_state_size);
        std::vector<int> thread_state(flat_width, 0);
        for (int t = 0; t < thread_count; ++t) {
            if (rng_engines[t] == nullptr) {
                continue;
            }
            const auto* local_ptr =
                static_cast<const unsigned char*>(gsl_rng_state(rng_engines[t]));
            const size_t base = static_cast<size_t>(t) * thread_state_size;
            for (size_t idx = 0; idx < thread_state_size; ++idx) {
                thread_state[base + idx] = static_cast<int>(local_ptr[idx]);
            }
        }
        replace_dataset_int("rng_thread_state", {flat_width}, thread_state);
        replace_dataset_int("rng_thread_count", {1}, {thread_count});
        replace_dataset_int("rng_thread_state_size", {1}, {static_cast<int>(thread_state_size)});
    }
    file.flush(H5F_SCOPE_GLOBAL);
    maybe_test_sleep("SARCOMERE_TEST_SLEEP_BEFORE_RESUME_PROMOTE", "resume_promote");
    promote_temp_group(file, final_resume_path, temp_resume_path, backup_resume_path);
}

int Sarcomere::load_state(int& n_frames, int frame_index){
    {
        H5::H5File file(filename, H5F_ACC_RDONLY);
        H5::Group meta = file.openGroup("/meta");
        if (aa_zdisc_tip_contact) {
            constexpr const char* expected_scope =
                "all_nearby_myosin_tips_head_independent";
            if (!meta.attrExists("aa_tip_accessibility_scope")) {
                throw std::runtime_error(
                    "Cannot resume ACTN2 tip contact from a checkpoint "
                    "without head-independent accessibility metadata");
            }
            H5::Attribute scope_attr =
                meta.openAttribute("aa_tip_accessibility_scope");
            H5::StrType scope_type = scope_attr.getStrType();
            std::string stored_scope;
            scope_attr.read(scope_type, stored_scope);
            if (stored_scope != expected_scope) {
                throw std::runtime_error(
                    "Cannot resume ACTN2 tip contact with legacy "
                    "motor-scoped accessibility semantics");
            }
        }
        if (dynamic_am) {
            int stored_model_version = 0;
            if (!meta.attrExists("am_model_version")) {
                throw std::runtime_error(
                    "Cannot resume dynamic AM: trajectory predates AA/AM spacing schema v9");
            }
            meta.openAttribute("am_model_version").read(
                H5::PredType::NATIVE_INT, &stored_model_version);
            if (stored_model_version < 22 || !meta.attrExists("am_heads_per_half")) {
                throw std::runtime_error(
                    "Cannot resume dynamic AM checkpoint predating continuous "
                    "strain-dependent retention, crown mechanics, and the ATP-limited "
                    "strong-head off-rate, motor-allocated ACTN2 catch load, and "
                    "overlap-scaled neighbor head capacity, and fixed per-link "
                    "ACTN2 mechanics, finite-systole heartbeat gating, and "
                    "occupied-pair angular compliance and optional explicit "
                    "per-head transverse compliance with independently coupled "
                    "axial/transverse escape barriers, global actin-site "
                    "ownership, and nonproductive double-overlap binding v22");
            }
            int stored_heads = 1;
            meta.openAttribute("am_heads_per_half").read(
                H5::PredType::NATIVE_INT, &stored_heads);
            if (stored_heads != am_heads_per_half) {
                throw std::runtime_error(
                    "Cannot resume with a different am_heads_per_half than trajectory metadata");
            }
            auto require_int_attr = [&](const char* name, int expected) {
                if (!meta.attrExists(name)) {
                    throw std::runtime_error(
                        std::string("Cannot resume: missing AM metadata ") + name);
                }
                int stored = 0;
                meta.openAttribute(name).read(H5::PredType::NATIVE_INT, &stored);
                if (stored != expected) {
                    throw std::runtime_error(
                        std::string("Cannot resume with a different ") + name);
                }
            };
            auto require_double_attr = [&](const char* name, double expected) {
                if (!meta.attrExists(name)) {
                    throw std::runtime_error(
                        std::string("Cannot resume: missing AM metadata ") + name);
                }
                double stored = 0.0;
                meta.openAttribute(name).read(
                    H5::PredType::NATIVE_DOUBLE, &stored);
                if (std::abs(stored - expected) >
                    1e-12 * std::max({1.0, std::abs(stored), std::abs(expected)})) {
                    throw std::runtime_error(
                        std::string("Cannot resume with a different ") + name);
                }
            };
            require_int_attr("am_max_heads_per_actin_half",
                             am_max_heads_per_actin_half);
            require_int_attr("am_sites_per_actin", am_sites_per_actin);
            require_int_attr("am_max_heads_per_actin",
                             am_max_heads_per_actin);
            require_int_attr("am_thin_activation_enabled",
                             am_thin_activation_enabled ? 1 : 0);
            require_int_attr("am_spatial_heads", am_spatial_heads ? 1 : 0);
            require_int_attr("sterics_enabled", sterics_enabled ? 1 : 0);
            require_int_attr("am_tether_saturation_heads",
                             am_tether_saturation_heads);
            require_double_attr("am_head_transverse_k",
                                am_head_transverse_k);
            require_double_attr("am_alignment_k", am_alignment_k);
            require_double_attr("am_alignment_slack_rad",
                                am_alignment_slack);
            require_int_attr("am_strain_dependent_detachment",
                             am_strain_dependent_detachment ? 1 : 0);
            require_double_attr("am_strain_detachment_theta",
                                am_strain_detachment_theta);
            require_double_attr(
                "am_transverse_strain_detachment_theta",
                am_transverse_strain_detachment_theta);
            require_double_attr("am_strain_detachment_max_exponent",
                                am_strain_detachment_max_exponent);
            require_double_attr("am_strong_off_floor_fraction",
                                am_strong_off_floor_fraction);
            require_double_attr("am_cycle_period", am_cycle_period);
            require_double_attr("am_activation_rise", am_activation_rise);
            require_double_attr("am_activation_decay", am_activation_decay);
            require_double_attr("am_activation_floor", am_activation_floor);
            require_double_attr("am_activation_active_fraction",
                                am_activation_active_fraction);
            require_double_attr("actin_radius", actin_radius);
            require_double_attr("am_surface_gap", am_surface_gap);
            require_double_attr("am_steric_k", am_steric_k);
            require_double_attr("am_preferred_centerline", am_preferred_centerline);
            require_double_attr("aa_preferred_centerline", aa_optimal);
            require_double_attr("aa_spacing_ratio", aa_spacing_ratio);
            require_double_attr("aa_capture_margin", aa_capture_margin);
            require_int_attr("am_thin_feedback_enabled",
                             am_thin_feedback_enabled ? 1 : 0);
            require_double_attr("am_q_calcium_open_rate", am_q_calcium_open_rate);
            require_double_attr("am_q_calcium_half", am_q_calcium_half);
            require_double_attr("am_q_calcium_hill", am_q_calcium_hill);
            require_double_attr("am_q_strong_feedback_rate",
                                am_q_strong_feedback_rate);
            require_double_attr("am_q_strong_half", am_q_strong_half);
            require_double_attr("am_q_strong_hill", am_q_strong_hill);
            require_double_attr("am_q_close_rate", am_q_close_rate);
            if (meta.attrExists("am_q_initialized_at_calcium_equilibrium")) {
                require_int_attr(
                    "am_q_initialized_at_calcium_equilibrium",
                    am_q_initialized_at_calcium_equilibrium ? 1 : 0);
                require_double_attr("am_initial_q", am_initial_q);
            } else if (am_q_initialized_at_calcium_equilibrium) {
                throw std::runtime_error(
                    "Cannot resume calcium-equilibrium q initialization from "
                    "a trajectory without q-initialization metadata");
            }
            require_int_attr("am_thick_feedback_enabled",
                             am_thick_feedback_enabled ? 1 : 0);
            require_double_attr("am_thick_recruit_rate", am_thick_recruit_rate);
            require_double_attr("am_thick_fold_rate", am_thick_fold_rate);
            require_double_attr("am_thick_tension_scale", am_thick_tension_scale);
            require_double_attr("am_thick_initial_on_fraction",
                                am_thick_initial_on_fraction);
            require_double_attr("am_thick_force_clamp", am_thick_force_clamp);
        }
        auto require_promoted_int = [&](const char* name, int expected,
                                        int legacy_default) {
            if (!meta.attrExists(name)) {
                if (expected != legacy_default) {
                    throw std::runtime_error(
                        std::string("Cannot resume: trajectory lacks metadata for ") + name);
                }
                return;
            }
            int stored = 0;
            meta.openAttribute(name).read(H5::PredType::NATIVE_INT, &stored);
            if (stored != expected) {
                throw std::runtime_error(
                    std::string("Cannot resume with a different ") + name);
            }
        };
        auto require_promoted_double = [&](const char* name, double expected,
                                           double legacy_default) {
            if (!meta.attrExists(name)) {
                if (std::abs(expected - legacy_default) > 1e-12) {
                    throw std::runtime_error(
                        std::string("Cannot resume: trajectory lacks metadata for ") + name);
                }
                return;
            }
            double stored = 0.0;
            meta.openAttribute(name).read(
                H5::PredType::NATIVE_DOUBLE, &stored);
            if (std::abs(stored - expected) >
                1e-12 * std::max({1.0, std::abs(stored), std::abs(expected)})) {
                throw std::runtime_error(
                    std::string("Cannot resume with a different ") + name);
            }
        };
        require_promoted_int(
            "actin_turnover_enabled", actin_turnover_enabled ? 1 : 0, 0);
        require_promoted_double(
            "actin_turnover_rate", actin_turnover_rate, 0.0);
        require_promoted_double(
            "actin_turnover_start_time", actin_turnover_start_time, 0.0);
        require_promoted_int(
            "actin_turnover_max_placement_attempts",
            actin_turnover_max_placement_attempts, 1000);
        require_promoted_double(
            "actin_turnover_clearance_ratio",
            actin_turnover_clearance_ratio, 1.0);
        require_promoted_double(
            "aa_barbed_lateral_range", aa_barbed_lateral_range, 0.0);
        require_promoted_double(
            "aa_barbed_lateral_k", aa_barbed_lateral_k, 0.0);
        require_promoted_double(
            "aa_barbed_lateral_force_max",
            aa_barbed_lateral_force_max, 0.0);
        require_promoted_double(
            "aa_barbed_lateral_axial_window",
            aa_barbed_lateral_axial_window, 0.25);
        require_promoted_int(
            "aa_register_all_statuses", aa_register_all_statuses ? 1 : 0, 0);
        require_promoted_double(
            "aa_axial_register_force_max", aa_axial_register_force_max, 0.0);
        require_promoted_double("aa_barbed_zone", aa_barbed_zone, -1.0);
        require_promoted_int(
            "titin_mode", static_cast<int>(titin_mode),
            static_cast<int>(TitinMode::legacy_crosslinked));
        require_promoted_double("titin_k", titin_k, 0.0);
        require_promoted_double("k_mm", k_mm, 0.0);
        if (k_mm > 0.0) {
            require_promoted_int("myomesin_schema_version", 1, 0);
            require_promoted_int("myomesin_sites_per_filament",
                                 myomesin_sites_per_filament, 2);
            require_promoted_double("myomesin_k_on", myomesin_k_on, 5000.0);
            require_promoted_double("myomesin_k_off", myomesin_k_off, 1.0);
            require_promoted_double("myomesin_axial_capture",
                                    myomesin_axial_capture, 0.25);
            require_promoted_double("myomesin_radial_rest",
                                    myomesin_radial_rest, 0.11);
            require_promoted_double("myomesin_radial_capture",
                                    myomesin_radial_capture, 0.04);
            require_promoted_double("myomesin_radial_k",
                                    myomesin_radial_k, 0.0);
            require_promoted_double("myomesin_alignment_cos",
                                    myomesin_alignment_cos, 0.94);
            require_promoted_double("myomesin_bell_distance",
                                    myomesin_bell_distance, 0.0005);
            require_promoted_double("myomesin_force_max",
                                    myomesin_force_max, 10.0);
        }
        if (titin_k > 0.0 || meta.attrExists("titin_rest_length")) {
            require_promoted_double(
                "titin_rest_length", titin_rest_length, 0.0);
        }
        if (titin_mode == TitinMode::explicit_dynamic) {
            require_promoted_int("titin_schema_version", 4, 0);
            require_promoted_double(
                "titin_k_off_unsupported",
                explicit_titin_k_off_unsupported, 0.0);
            require_promoted_int(
                "titin_rebind_requires_actn2",
                explicit_titin_rebind_requires_actn2 ? 1 : 0, 0);
            require_promoted_int(
                "titin_require_live_actn2",
                explicit_titin_require_live_actn2 ? 1 : 0, 0);
            require_promoted_int(
                "titin_links_per_half", explicit_titin_per_half, 3);
            require_promoted_double(
                "titin_capture_distance",
                explicit_titin_capture_distance, 0.70);
            require_promoted_double("titin_k_on", explicit_titin_k_on, 5000.0);
            require_promoted_double("titin_k_off", explicit_titin_k_off, 0.5);
            require_promoted_int(
                "titin_strain_gated_kinetics",
                explicit_titin_strain_gated_kinetics ? 1 : 0, 1);
            require_promoted_double(
                "titin_transition_theta",
                explicit_titin_transition_theta, 0.5);
            require_promoted_double(
                "titin_rebind_memory",
                explicit_titin_rebind_memory, 0.0);
            require_promoted_double(
                "titin_rebind_rate_multiplier",
                explicit_titin_rebind_rate_multiplier, 1.0);
        }
        if (aa_local_multilink || meta.attrExists("aa_local_multilink")) {
            auto require_local_int = [&](const char* name, int expected) {
                if (!meta.attrExists(name)) {
                    throw std::runtime_error(
                        std::string("Cannot resume: missing local ACTN2 metadata ") + name);
                }
                int stored = 0;
                meta.openAttribute(name).read(H5::PredType::NATIVE_INT, &stored);
                if (stored != expected) {
                    throw std::runtime_error(
                        std::string("Cannot resume with a different ") + name);
                }
            };
            auto require_local_double = [&](const char* name, double expected) {
                if (!meta.attrExists(name)) {
                    throw std::runtime_error(
                        std::string("Cannot resume: missing local ACTN2 metadata ") + name);
                }
                double stored = 0.0;
                meta.openAttribute(name).read(
                    H5::PredType::NATIVE_DOUBLE, &stored);
                if (std::abs(stored - expected) >
                    1e-12 * std::max({1.0, std::abs(stored), std::abs(expected)})) {
                    throw std::runtime_error(
                        std::string("Cannot resume with a different ") + name);
                }
            };
            require_local_int("aa_local_multilink", aa_local_multilink ? 1 : 0);
            if (aa_local_multilink) {
                require_local_int("aa_local_link_schema",
                                  aa_arm_state_model ? 5 : 4);
                require_local_int("aa_prebundle_links_per_pair",
                                  aa_prebundle_links_per_pair);
                require_local_double("aa_site_exclusion", aa_site_exclusion);
                require_local_double("aa_overlap_capture", aa_overlap_capture);
                require_local_double(
                    "aa_axial_register_ratio", aa_axial_register_ratio);
                require_local_double(
                    "aa_status2_registry_on", aa_status2_registry_on);
                require_local_double(
                    "aa_status2_registry_off", aa_status2_registry_off);
                require_local_double(
                    "aa_register_saturation_error",
                    aa_register_saturation_error);
                require_local_double("aa_perpendicular_k", aa_perpendicular_k);
                if (meta.attrExists("aa_status1_radial_k")) {
                    require_local_double(
                        "aa_status1_radial_k", aa_status1_radial_k);
                } else if (aa_status1_radial_k > 0.0) {
                    throw std::runtime_error(
                        "Cannot resume an older trajectory with a non-legacy "
                        "aa_status1_radial_k");
                }
                require_promoted_double(
                    "aa_status1_site_diffusion",
                    aa_status1_site_diffusion, 0.0);
                require_promoted_int(
                    "aa_status2_locked_until_detach",
                    aa_status2_locked_until_detach ? 1 : 0, 0);
                require_promoted_double(
                    "aa_status2_unlock_rate", aa_status2_unlock_rate, 0.0);
                require_promoted_double(
                    "aa_status2_unlock_force_scale",
                    aa_status2_unlock_force_scale, 0.1);
                require_promoted_double(
                    "aa_status2_unlock_registry_scale",
                    aa_status2_unlock_registry_scale, 0.02);
                require_promoted_int(
                    "aa_arm_state_model", aa_arm_state_model ? 1 : 0, 0);
                require_promoted_double(
                    "aa_arm_open_rate_weak", aa_arm_open_rate_weak, 20.0);
                require_promoted_double(
                    "aa_arm_open_rate_mature", aa_arm_open_rate_mature, 1.0);
                require_promoted_double(
                    "aa_arm_bell_distance", aa_arm_bell_distance, 0.001);
                require_promoted_double(
                    "aa_tether_off_rate", aa_tether_off_rate, 5.0);
                require_promoted_double(
                    "aa_tether_rebind_rate", aa_tether_rebind_rate, 100.0);
                require_promoted_double(
                    "aa_tether_capture_distance",
                    aa_tether_capture_distance, 0.10);
                require_promoted_double(
                    "aa_arm_rebind_theta", aa_arm_rebind_theta, 0.5);
                require_promoted_double(
                    "aa_actn2_contact_radius", aa_actn2_contact_radius, 0.005);
                require_promoted_double(
                    "aa_myosin_contact_epsilon",
                    aa_myosin_contact_epsilon, 0.00428);
                require_promoted_double(
                    "aa_myosin_contact_force_cap",
                    aa_myosin_contact_force_cap, 40.0);
                require_promoted_double(
                    "aa_myosin_contact_shell_width",
                    aa_myosin_contact_shell_width, 0.0);
                require_promoted_double(
                    "aa_myosin_contact_shell_energy",
                    aa_myosin_contact_shell_energy, 0.0);
                require_promoted_int(
                    "aa_zdisc_tip_contact", aa_zdisc_tip_contact ? 1 : 0, 0);
                require_promoted_double(
                    "aa_zdisc_tip_onset", aa_zdisc_tip_onset, 0.050);
                require_promoted_double(
                    "aa_zdisc_tip_stiffness", aa_zdisc_tip_stiffness, 1000.0);
                require_promoted_double(
                    "aa_zdisc_tip_force_max", aa_zdisc_tip_force_max, 20.0);
                require_promoted_double(
                    "aa_zdisc_tip_transverse_flat",
                    aa_zdisc_tip_transverse_flat, 0.100);
                require_promoted_double(
                    "aa_zdisc_tip_transverse_cutoff",
                    aa_zdisc_tip_transverse_cutoff, 0.130);
                require_promoted_double(
                    "aa_occlusion_rupture_penetration",
                    aa_occlusion_rupture_penetration, 0.010);
                require_promoted_double(
                    "aa_association_max_direction_dot",
                    aa_association_max_direction_dot, 1.0);
                require_promoted_double(
                    "aa_link_rupture_length", aa_link_rupture_length, 0.0);
                require_local_int("aa_stiffness_reference_links",
                                  aa_stiffness_reference_links);
                require_local_int("aa_total_link_pool", aa_total_link_pool);
            }
        }
        if (meta.attrExists("aa_status_mode")) {
            int stored_mode = 0;
            meta.openAttribute("aa_status_mode").read(
                H5::PredType::NATIVE_INT, &stored_mode);
            if (stored_mode != static_cast<int>(aa_status_mode)) {
                throw std::runtime_error(
                    "Cannot resume with a different aa_status_mode than the trajectory metadata");
            }
            if (meta.attrExists("aa_status2_force_on") &&
                meta.attrExists("aa_status2_force_off")) {
                double stored_on = 0.0;
                double stored_off = 0.0;
                meta.openAttribute("aa_status2_force_on").read(
                    H5::PredType::NATIVE_DOUBLE, &stored_on);
                meta.openAttribute("aa_status2_force_off").read(
                    H5::PredType::NATIVE_DOUBLE, &stored_off);
                if (std::abs(stored_on - aa_status2_force_on) > 1e-12 ||
                    std::abs(stored_off - aa_status2_force_off) > 1e-12) {
                    throw std::runtime_error(
                        "Cannot resume with different AA status force thresholds than the trajectory metadata");
                }
            } else if (aa_status_mode !=
                       actn2_load_kinetics::StatusMode::legacy_motor_identity) {
                throw std::runtime_error(
                    "Cannot resume load-classified AA kinetics without stored force thresholds");
            }
        } else if (aa_status_mode !=
                   actn2_load_kinetics::StatusMode::legacy_motor_identity) {
            throw std::runtime_error(
                "Cannot resume load-classified AA kinetics from a trajectory without AA status-mode metadata");
        }
    }
    int target_frame = load_from_file(filename, actin, myosin, actin_actin_bonds, n_frames, frame_index);

    const size_t aa_stride = static_cast<size_t>(actin.n) * static_cast<size_t>(actin.n);
    const size_t am_stride = static_cast<size_t>(actin.n) * static_cast<size_t>(myosin.n);

    auto load_state_vector = [&](H5::Group& group, const std::string& dataset_name,
                                 size_t expected_width, std::vector<int>& out) -> bool {
        try {
            if (!group.nameExists(dataset_name)) {
                return false;
            }
            std::vector<hsize_t> dims;
            std::vector<double> raw = load_from_dataset(group, dataset_name, dims);
            if (dims.size() < 2 || target_frame < 0 || target_frame >= static_cast<int>(dims[0])) {
                return false;
            }
            if (static_cast<size_t>(dims[1]) != expected_width) {
                return false;
            }
            size_t offset = static_cast<size_t>(target_frame) * expected_width;
            out.resize(expected_width);
            for (size_t i = 0; i < expected_width; ++i) {
                out[i] = static_cast<int>(std::llround(raw[offset + i]));
            }
            return true;
        } catch (H5::Exception&) {
            return false;
        }
    };

    auto load_state_vector_double = [&](H5::Group& group, const std::string& dataset_name,
                                        size_t expected_width, std::vector<double>& out) -> bool {
        try {
            if (!group.nameExists(dataset_name)) {
                return false;
            }
            std::vector<hsize_t> dims;
            std::vector<double> raw = load_from_dataset(group, dataset_name, dims);
            if (dims.size() < 2 || target_frame < 0 || target_frame >= static_cast<int>(dims[0])) {
                return false;
            }
            if (static_cast<size_t>(dims[1]) != expected_width) {
                return false;
            }
            size_t offset = static_cast<size_t>(target_frame) * expected_width;
            out.resize(expected_width);
            std::copy_n(raw.begin() + static_cast<std::ptrdiff_t>(offset),
                        static_cast<std::ptrdiff_t>(expected_width),
                        out.begin());
            return true;
        } catch (H5::Exception&) {
            return false;
        }
    };

    auto assign_aa_matrix = [&](const std::vector<int>& flat, auto& matrix) {
        if (flat.size() != aa_stride) {
            return;
        }
        for (int i = 0; i < actin.n; ++i) {
            for (int j = 0; j < actin.n; ++j) {
                matrix[i][j] = flat[static_cast<size_t>(i) * actin.n + j];
            }
        }
    };

    auto assign_am_matrix = [&](const std::vector<int>& flat, auto& matrix) {
        if (flat.size() != am_stride) {
            return;
        }
        for (int i = 0; i < actin.n; ++i) {
            for (int j = 0; j < myosin.n; ++j) {
                matrix[i][j] = flat[static_cast<size_t>(i) * myosin.n + j];
            }
        }
    };

    auto load_state_tensor_int = [&](H5::Group& group, const std::string& dataset_name,
                                     size_t expected_rows, size_t expected_cols,
                                     std::vector<int>& out) -> bool {
        try {
            if (!group.nameExists(dataset_name)) {
                return false;
            }
            std::vector<hsize_t> dims;
            std::vector<double> raw = load_from_dataset(group, dataset_name, dims);
            if (dims.size() < 3 || target_frame < 0 || target_frame >= static_cast<int>(dims[0])) {
                return false;
            }
            if (static_cast<size_t>(dims[1]) != expected_rows ||
                static_cast<size_t>(dims[2]) != expected_cols) {
                return false;
            }
            const size_t frame_width = expected_rows * expected_cols;
            const size_t offset = static_cast<size_t>(target_frame) * frame_width;
            out.resize(frame_width);
            for (size_t i = 0; i < frame_width; ++i) {
                out[i] = static_cast<int>(std::llround(raw[offset + i]));
            }
            return true;
        } catch (H5::Exception&) {
            return false;
        }
    };

    auto load_state_tensor_double = [&](H5::Group& group, const std::string& dataset_name,
                                        size_t expected_rows, size_t expected_cols,
                                        std::vector<double>& out) -> bool {
        try {
            if (!group.nameExists(dataset_name)) {
                return false;
            }
            std::vector<hsize_t> dims;
            std::vector<double> raw = load_from_dataset(group, dataset_name, dims);
            if (dims.size() < 3 || target_frame < 0 ||
                target_frame >= static_cast<int>(dims[0])) {
                return false;
            }
            if (static_cast<size_t>(dims[1]) != expected_rows ||
                static_cast<size_t>(dims[2]) != expected_cols) {
                return false;
            }
            const size_t frame_width = expected_rows * expected_cols;
            const size_t offset = static_cast<size_t>(target_frame) * frame_width;
            out.resize(frame_width);
            std::copy_n(raw.begin() + static_cast<std::ptrdiff_t>(offset),
                        static_cast<std::ptrdiff_t>(frame_width), out.begin());
            return true;
        } catch (H5::Exception&) {
            return false;
        }
    };

    auto load_fixed_vector_int = [&](H5::Group& group, const std::string& dataset_name,
                                     size_t expected_width, std::vector<int>& out) -> bool {
        try {
            if (!group.nameExists(dataset_name)) {
                return false;
            }
            std::vector<hsize_t> dims;
            std::vector<double> raw = load_from_dataset(group, dataset_name, dims);
            if (dims.size() != 1 || static_cast<size_t>(dims[0]) != expected_width) {
                return false;
            }
            out.resize(expected_width);
            for (size_t i = 0; i < expected_width; ++i) {
                out[i] = static_cast<int>(std::llround(raw[i]));
            }
            return true;
        } catch (H5::Exception&) {
            return false;
        }
    };

    auto load_fixed_vector_double = [&](H5::Group& group, const std::string& dataset_name,
                                        size_t expected_width, std::vector<double>& out) -> bool {
        try {
            if (!group.nameExists(dataset_name)) {
                return false;
            }
            std::vector<hsize_t> dims;
            out = load_from_dataset(group, dataset_name, dims);
            return dims.size() == 1 && static_cast<size_t>(dims[0]) == expected_width;
        } catch (H5::Exception&) {
            return false;
        }
    };

    auto load_any_vector_double = [&](H5::Group& group, const std::string& dataset_name,
                                      std::vector<double>& out) -> bool {
        try {
            if (!group.nameExists(dataset_name)) {
                return false;
            }
            std::vector<hsize_t> dims;
            out = load_from_dataset(group, dataset_name, dims);
            return dims.size() == 1;
        } catch (H5::Exception&) {
            return false;
        }
    };

    auto load_fixed_matrix_int = [&](H5::Group& group, const std::string& dataset_name,
                                     size_t expected_rows, size_t expected_cols,
                                     std::vector<int>& out) -> bool {
        try {
            if (!group.nameExists(dataset_name)) {
                return false;
            }
            std::vector<hsize_t> dims;
            std::vector<double> raw = load_from_dataset(group, dataset_name, dims);
            if (dims.size() != 2 || static_cast<size_t>(dims[0]) != expected_rows ||
                static_cast<size_t>(dims[1]) != expected_cols) {
                return false;
            }
            out.resize(expected_rows * expected_cols);
            for (size_t i = 0; i < out.size(); ++i) {
                out[i] = static_cast<int>(std::llround(raw[i]));
            }
            return true;
        } catch (H5::Exception&) {
            return false;
        }
    };

    auto load_fixed_matrix_double = [&](H5::Group& group, const std::string& dataset_name,
                                        size_t expected_rows, size_t expected_cols,
                                        std::vector<double>& out) -> bool {
        try {
            if (!group.nameExists(dataset_name)) {
                return false;
            }
            std::vector<hsize_t> dims;
            out = load_from_dataset(group, dataset_name, dims);
            return dims.size() == 2 && static_cast<size_t>(dims[0]) == expected_rows &&
                   static_cast<size_t>(dims[1]) == expected_cols;
        } catch (H5::Exception&) {
            return false;
        }
    };

    auto clear_aa_state = [&](auto& bonds, auto& status, auto& lifetime) {
        for (int i = 0; i < actin.n; ++i) {
            std::fill(bonds[i].begin(), bonds[i].end(), 0);
            std::fill(status[i].begin(), status[i].end(), 0);
            std::fill(lifetime[i].begin(), lifetime[i].end(), 0);
        }
    };

    auto clear_am_state = [&](auto& bonds) {
        for (int i = 0; i < actin.n; ++i) {
            std::fill(bonds[i].begin(), bonds[i].end(), 0);
        }
    };

    auto clear_recovery_state = [&]() {
        for (int i = 0; i < actin.n; ++i) {
            std::fill(actin_recovery_until[i].begin(), actin_recovery_until[i].end(), 0);
        }
    };

    auto rebuild_explicit_titin_ownership = [&]() {
        explicit_titin_owner_by_actin.assign(static_cast<size_t>(actin.n), -1);
        for (size_t index = 0; index < explicit_titins.size(); ++index) {
            auto& titin = explicit_titins[index];
            if (!titin.bound) {
                titin.bound_actin = -1;
                continue;
            }
            if (titin.bound_actin < 0 || titin.bound_actin >= actin.n) {
                throw std::runtime_error(
                    "Explicit titin checkpoint contains an invalid bound-actin index");
            }
            int& owner = explicit_titin_owner_by_actin[
                static_cast<size_t>(titin.bound_actin)];
            if (owner >= 0) {
                throw std::runtime_error(
                    "Explicit titin checkpoint assigns one actin to multiple titins");
            }
            owner = static_cast<int>(index);
        }
    };

    auto restore_aa_state_from_pairs = [&](H5::Group& group,
                                           const std::string& count_name,
                                           const std::string& pair_name,
                                           const std::string& status_name,
                                           const std::string& lifetime_name,
                                           auto& bonds,
                                           auto& status,
                                           auto& lifetime) -> bool {
        std::vector<int> counts;
        std::vector<int> pairs;
        std::vector<int> flat_status;
        std::vector<int> flat_lifetime;
        const size_t max_aa_pairs = static_cast<size_t>(actin.n) *
                                    static_cast<size_t>(std::max(10, 2 * max_myosin_bonds));
        if (!load_state_vector(group, count_name, 1, counts) ||
            !load_state_tensor_int(group, pair_name, max_aa_pairs, 2, pairs) ||
            !load_state_tensor_int(group, status_name, max_aa_pairs, 1, flat_status) ||
            !load_state_tensor_int(group, lifetime_name, max_aa_pairs, 1, flat_lifetime)) {
            return false;
        }
        const size_t count = static_cast<size_t>(std::max(0, counts[0]));
        clear_aa_state(bonds, status, lifetime);
        for (size_t idx = 0; idx < count; ++idx) {
            const int a = pairs[2 * idx];
            const int b = pairs[2 * idx + 1];
            if (a < 0 || a >= actin.n || b < 0 || b >= actin.n || a == b) {
                continue;
            }
            bonds[a][b] = 1;
            bonds[b][a] = 1;
            status[a][b] = flat_status[idx];
            status[b][a] = flat_status[idx];
            lifetime[a][b] = flat_lifetime[idx];
            lifetime[b][a] = flat_lifetime[idx];
        }
        return true;
    };

    auto restore_am_state_from_pairs = [&](H5::Group& group,
                                           const std::string& count_name,
                                           const std::string& pair_name,
                                           auto& bonds) -> bool {
        std::vector<int> counts;
        std::vector<int> pairs;
        const size_t max_am_pairs = static_cast<size_t>(myosin.n) *
                                    static_cast<size_t>(max_myosin_bonds);
        if (!load_state_vector(group, count_name, 1, counts) ||
            !load_state_tensor_int(group, pair_name, max_am_pairs, 2, pairs)) {
            return false;
        }
        const size_t count = static_cast<size_t>(std::max(0, counts[0]));
        clear_am_state(bonds);
        for (size_t idx = 0; idx < count; ++idx) {
            const int a = pairs[2 * idx];
            const int m = pairs[2 * idx + 1];
            if (a < 0 || a >= actin.n || m < 0 || m >= myosin.n) {
                continue;
            }
            bonds[a][m] = 1;
        }
        return true;
    };

    auto restore_recovery_state_from_pairs = [&](H5::Group& group,
                                                 const std::string& count_name,
                                                 const std::string& pair_name,
                                                 const std::string& value_name) -> bool {
        std::vector<int> counts;
        std::vector<int> pairs;
        std::vector<int> until_values;
        const size_t max_recovery_pairs = static_cast<size_t>(actin.n) *
                                          static_cast<size_t>(std::max(10, 2 * max_myosin_bonds));
        if (!load_state_vector(group, count_name, 1, counts) ||
            !load_state_tensor_int(group, pair_name, max_recovery_pairs, 2, pairs) ||
            !load_state_tensor_int(group, value_name, max_recovery_pairs, 1, until_values)) {
            return false;
        }
        const size_t count = static_cast<size_t>(std::max(0, counts[0]));
        clear_recovery_state();
        for (size_t idx = 0; idx < count; ++idx) {
            const int a = pairs[2 * idx];
            const int b = pairs[2 * idx + 1];
            if (a < 0 || a >= actin.n || b < 0 || b >= actin.n || a == b) {
                continue;
            }
            const RecoveryStep until =
                static_cast<RecoveryStep>(std::max(0, until_values[idx]));
            actin_recovery_until[a][b] = until;
            actin_recovery_until[b][a] = until;
        }
        return true;
    };

    auto restore_aa_state_from_fixed_pairs = [&](H5::Group& group,
                                                 const std::string& count_name,
                                                 const std::string& pair_name,
                                                 const std::string& status_name,
                                                 const std::string& lifetime_name,
                                                 auto& bonds,
                                                 auto& status,
                                                 auto& lifetime) -> bool {
        std::vector<int> counts;
        std::vector<int> pairs;
        std::vector<int> flat_status;
        std::vector<int> flat_lifetime;
        const size_t max_aa_pairs = static_cast<size_t>(actin.n) *
                                    static_cast<size_t>(std::max(10, 2 * max_myosin_bonds));
        if (!load_fixed_vector_int(group, count_name, 1, counts) ||
            !load_fixed_matrix_int(group, pair_name, max_aa_pairs, 2, pairs) ||
            !load_fixed_matrix_int(group, status_name, max_aa_pairs, 1, flat_status) ||
            !load_fixed_matrix_int(group, lifetime_name, max_aa_pairs, 1, flat_lifetime)) {
            return false;
        }
        const size_t count = static_cast<size_t>(std::max(0, counts[0]));
        clear_aa_state(bonds, status, lifetime);
        for (size_t idx = 0; idx < count; ++idx) {
            const int a = pairs[2 * idx];
            const int b = pairs[2 * idx + 1];
            if (a < 0 || a >= actin.n || b < 0 || b >= actin.n || a == b) {
                continue;
            }
            bonds[a][b] = 1;
            bonds[b][a] = 1;
            status[a][b] = flat_status[idx];
            status[b][a] = flat_status[idx];
            lifetime[a][b] = flat_lifetime[idx];
            lifetime[b][a] = flat_lifetime[idx];
        }
        return true;
    };

    auto restore_am_state_from_fixed_pairs = [&](H5::Group& group,
                                                 const std::string& count_name,
                                                 const std::string& pair_name,
                                                 auto& bonds) -> bool {
        std::vector<int> counts;
        std::vector<int> pairs;
        const size_t max_am_pairs = static_cast<size_t>(myosin.n) *
                                    static_cast<size_t>(max_myosin_bonds);
        if (!load_fixed_vector_int(group, count_name, 1, counts) ||
            !load_fixed_matrix_int(group, pair_name, max_am_pairs, 2, pairs)) {
            return false;
        }
        const size_t count = static_cast<size_t>(std::max(0, counts[0]));
        clear_am_state(bonds);
        for (size_t idx = 0; idx < count; ++idx) {
            const int a = pairs[2 * idx];
            const int m = pairs[2 * idx + 1];
            if (a < 0 || a >= actin.n || m < 0 || m >= myosin.n) {
                continue;
            }
            bonds[a][m] = 1;
        }
        return true;
    };

    auto restore_recovery_state_from_fixed_pairs = [&](H5::Group& group,
                                                       const std::string& count_name,
                                                       const std::string& pair_name,
                                                       const std::string& value_name) -> bool {
        std::vector<int> counts;
        std::vector<int> pairs;
        std::vector<int> until_values;
        const size_t max_recovery_pairs = static_cast<size_t>(actin.n) *
                                          static_cast<size_t>(std::max(10, 2 * max_myosin_bonds));
        if (!load_fixed_vector_int(group, count_name, 1, counts) ||
            !load_fixed_matrix_int(group, pair_name, max_recovery_pairs, 2, pairs) ||
            !load_fixed_matrix_int(group, value_name, max_recovery_pairs, 1, until_values)) {
            return false;
        }
        const size_t count = static_cast<size_t>(std::max(0, counts[0]));
        clear_recovery_state();
        for (size_t idx = 0; idx < count; ++idx) {
            const int a = pairs[2 * idx];
            const int b = pairs[2 * idx + 1];
            if (a < 0 || a >= actin.n || b < 0 || b >= actin.n || a == b) {
                continue;
            }
            const RecoveryStep until =
                static_cast<RecoveryStep>(std::max(0, until_values[idx]));
            actin_recovery_until[a][b] = until;
            actin_recovery_until[b][a] = until;
        }
        return true;
    };

    bool loaded_state_group = false;
    bool restored_neighbor_cache = false;
    bool loaded_resume_snapshot = false;
    try {
        H5::H5File file(filename, H5F_ACC_RDWR);
        H5::Group group_state(file.openGroup("/state"));
        loaded_state_group = true;

        std::vector<int> flat;
        std::vector<hsize_t> dims;

        // Detect interrupted /state appends and rewind to the last complete checkpoint.
        // state_commit_count is appended last in save_state(), so its frame count marks
        // the number of fully committed state snapshots.
        const std::string state_commit_dataset =
            group_state.nameExists("state_commit_count") ? "state_commit_count" : "rng_thread_state_size";
        if (group_state.nameExists(state_commit_dataset)) {
            std::vector<double> dummy = load_from_dataset(group_state, state_commit_dataset, dims);
            (void)dummy;
            if (dims.size() >= 2) {
                const int safe_frames = static_cast<int>(dims[0]);
                const std::vector<std::string> state_datasets = {
                    "current_step",
                    "state_commit_count",
                    "aa_pairs_prev_count",
                    "aa_pairs_current_count",
                    "am_pairs_prev_count",
                    "am_pairs_current_count",
                    "aa_pairs_prev",
                    "aa_status_prev",
                    "aa_lifetime_prev",
                    "aa_pairs_current",
                    "aa_status_current",
                    "aa_lifetime_current",
                    "aa_attach_step_current",
                    "aa_status2_memory_until_current",
                    "am_pairs_prev",
                    "am_pairs_current",
                    "aa_recovery_count",
                    "aa_recovery_pairs",
                    "aa_recovery_until_values",
                    "neighbor_last_actin_x",
                    "neighbor_last_actin_y",
                    "neighbor_last_actin_z",
                    "neighbor_last_myosin_x",
                    "neighbor_last_myosin_y",
                    "neighbor_last_myosin_z",
                    "titin_bound_actin_current",
                    "titin_last_bound_actin_current",
                    "titin_attach_step_current",
                    "titin_last_detach_step_current",
                    "titin_distance_current",
                    "titin_extension_current",
                    "titin_force_current",
                    "titin_attachment_distance_current",
                    "titin_attachment_energy_current",
                    "titin_attachment_force_current",
                    "rng_main_state",
                    "rng_thread_state",
                    "rng_thread_count",
                    "rng_thread_state_size"
                };
                if (safe_frames >= 0 && dims[0] > static_cast<hsize_t>(safe_frames)) {
                    for (const auto& name : state_datasets) {
                        truncate_dataset_first_dim(group_state, name, static_cast<hsize_t>(safe_frames));
                    }
                    file.flush(H5F_SCOPE_GLOBAL);
                }
                if (target_frame >= safe_frames) {
                    target_frame = safe_frames - 1;
                    if (target_frame >= 0) {
                        load_from_file(filename, actin, myosin, actin_actin_bonds, n_frames, target_frame);
                    }
                }
            }
        }

        // Load current_step (if missing, fallback handled below).
        if (group_state.nameExists("current_step")) {
            std::vector<double> step_data = load_from_dataset(group_state, "current_step", dims);
            if (dims.size() >= 2 && target_frame >= 0 && target_frame < static_cast<int>(dims[0])) {
                current_step = static_cast<size_t>(std::llround(step_data[target_frame]));
            } else {
                current_step = 0;
            }
        } else {
            current_step = 0;
        }

        bool have_prev_aa_pairs = restore_aa_state_from_pairs(
            group_state,
            "aa_pairs_prev_count",
            "aa_pairs_prev",
            "aa_status_prev",
            "aa_lifetime_prev",
            actin_actin_bonds_prev,
            actin_actin_status_prev,
            actin_actin_lifetime_prev);
        if (!have_prev_aa_pairs) {
            if (load_state_vector(group_state, "actin_actin_bonds_prev", aa_stride, flat)) {
                assign_aa_matrix(flat, actin_actin_bonds_prev);
            }
            if (load_state_vector(group_state, "actin_actin_status_prev", aa_stride, flat)) {
                assign_aa_matrix(flat, actin_actin_status_prev);
            }
            if (load_state_vector(group_state, "actin_actin_lifetime_prev", aa_stride, flat)) {
                assign_aa_matrix(flat, actin_actin_lifetime_prev);
            }
        }

        bool have_prev_am_pairs = restore_am_state_from_pairs(
            group_state, "am_pairs_prev_count", "am_pairs_prev", am_bonds_prev);
        if (!have_prev_am_pairs && load_state_vector(group_state, "am_bonds_prev", am_stride, flat)) {
            assign_am_matrix(flat, am_bonds_prev);
        }

        if (restore_recovery_state_from_pairs(
                group_state, "aa_recovery_count", "aa_recovery_pairs", "aa_recovery_until_values")) {
            // loaded sparse recovery state
        } else if (load_state_tensor_int(group_state, "actin_recovery_until", aa_stride, 1, flat) ||
                   load_state_vector(group_state, "actin_recovery_until", aa_stride, flat)) {
            for (int i = 0; i < actin.n; ++i) {
                for (int j = 0; j < actin.n; ++j) {
                    actin_recovery_until[i][j] = static_cast<RecoveryStep>(std::max(
                        0, flat[static_cast<size_t>(i) * actin.n + j]));
                }
            }
        } else {
            for (int i = 0; i < actin.n; ++i) {
                std::fill(actin_recovery_until[i].begin(), actin_recovery_until[i].end(), 0);
            }
        }

        bool have_current_aa_pairs = restore_aa_state_from_pairs(
            group_state,
            "aa_pairs_current_count",
            "aa_pairs_current",
            "aa_status_current",
            "aa_lifetime_current",
            actin_actin_bonds,
            actin_actin_status,
            actin_actin_lifetime);
        bool have_current_bonds = false;
        bool have_current_status = false;
        bool have_current_lifetime = false;
        if (!have_current_aa_pairs) {
            have_current_bonds = load_state_vector(group_state, "actin_actin_bonds_current", aa_stride, flat);
            if (have_current_bonds) {
                assign_aa_matrix(flat, actin_actin_bonds);
            }

            have_current_status = load_state_vector(group_state, "actin_actin_status_current", aa_stride, flat);
            if (have_current_status) {
                assign_aa_matrix(flat, actin_actin_status);
            } else {
                // Legacy fallback: /actin/bonds stores only strong bonds.
                for (int i = 0; i < actin.n; ++i) {
                    for (int j = 0; j < actin.n; ++j) {
                        actin_actin_status[i][j] = (actin_actin_bonds[i][j] == 1) ? 2 : 0;
                    }
                }
            }

            have_current_lifetime = load_state_vector(group_state, "actin_actin_lifetime_current", aa_stride, flat);
            if (have_current_lifetime) {
                assign_aa_matrix(flat, actin_actin_lifetime);
            } else {
                for (int i = 0; i < actin.n; ++i) {
                    for (int j = 0; j < actin.n; ++j) {
                        if (actin_actin_bonds[i][j] == 1) {
                            actin_actin_lifetime[i][j] = std::max(1, actin_actin_lifetime_prev[i][j]);
                        } else {
                            actin_actin_lifetime[i][j] = 0;
                        }
                    }
                }
            }
        }
        for (int i = 0; i < actin.n; ++i) {
            std::fill(aa_attach_step[i].begin(), aa_attach_step[i].end(), -1);
        }
        {
            std::vector<int> counts_attach;
            std::vector<int> pairs_attach;
            std::vector<int> attach_values;
            const size_t max_aa_pairs = static_cast<size_t>(actin.n) *
                                        static_cast<size_t>(std::max(10, 2 * max_myosin_bonds));
            if (load_state_vector(group_state, "aa_pairs_current_count", 1, counts_attach) &&
                load_state_tensor_int(group_state, "aa_pairs_current", max_aa_pairs, 2, pairs_attach) &&
                load_state_tensor_int(group_state, "aa_attach_step_current", max_aa_pairs, 1, attach_values)) {
                const size_t count = static_cast<size_t>(std::max(0, counts_attach[0]));
                for (size_t idx = 0; idx < count; ++idx) {
                    const int a = pairs_attach[2 * idx];
                    const int b = pairs_attach[2 * idx + 1];
                    if (a < 0 || a >= actin.n || b < 0 || b >= actin.n || a == b) {
                        continue;
                    }
                    aa_attach_step[a][b] = attach_values[idx];
                    aa_attach_step[b][a] = attach_values[idx];
                }
            }
        }
        for (int i = 0; i < actin.n; ++i) {
            std::fill(aa_status2_memory_until[i].begin(),
                      aa_status2_memory_until[i].end(), 0);
        }
        {
            std::vector<int> counts_memory;
            std::vector<int> pairs_memory;
            std::vector<int> memory_values;
            const size_t max_aa_pairs = static_cast<size_t>(actin.n) *
                                        static_cast<size_t>(std::max(10, 2 * max_myosin_bonds));
            if (load_state_vector(group_state, "aa_pairs_current_count", 1, counts_memory) &&
                load_state_tensor_int(group_state, "aa_pairs_current", max_aa_pairs, 2, pairs_memory) &&
                load_state_tensor_int(group_state, "aa_status2_memory_until_current",
                                      max_aa_pairs, 1, memory_values)) {
                const size_t count = static_cast<size_t>(std::max(0, counts_memory[0]));
                for (size_t idx = 0; idx < count; ++idx) {
                    const int a = pairs_memory[2 * idx];
                    const int b = pairs_memory[2 * idx + 1];
                    if (a < 0 || a >= actin.n || b < 0 || b >= actin.n || a == b) continue;
                    const RecoveryStep until = static_cast<RecoveryStep>(
                        std::max(0, memory_values[idx]));
                    aa_status2_memory_until[a][b] = until;
                    aa_status2_memory_until[b][a] = until;
                }
            }
        }

        bool loaded_current_am_from_pairs = restore_am_state_from_pairs(
            group_state, "am_pairs_current_count", "am_pairs_current", am_bonds);
        bool loaded_current_am_from_dense = false;
        if (!loaded_current_am_from_pairs) {
            loaded_current_am_from_dense =
                load_state_vector(group_state, "am_bonds_current", am_stride, flat);
        }
        if (loaded_current_am_from_dense) {
            assign_am_matrix(flat, am_bonds);
        } else if (!loaded_current_am_from_pairs) {
            for (int i = 0; i < actin.n; ++i) {
                std::fill(am_bonds[i].begin(), am_bonds[i].end(), 0);
            }
            // Fallback for legacy files: reconstruct current am_bonds from /actin_myo/bonds.
            try {
                H5::Group group_am(file.openGroup("/actin_myo"));
                std::vector<double> am_bonds_all = load_from_dataset(group_am, "bonds", dims);
                if (dims.size() >= 3 && target_frame >= 0 && target_frame < static_cast<int>(dims[0])) {
                    size_t bonds_per_frame = static_cast<size_t>(dims[1]) * static_cast<size_t>(dims[2]);
                    size_t start = static_cast<size_t>(target_frame) * bonds_per_frame;
                    for (size_t idx = 0; idx + 1 < bonds_per_frame; idx += 2) {
                        int a = static_cast<int>(std::llround(am_bonds_all[start + idx]));
                        int m = static_cast<int>(std::llround(am_bonds_all[start + idx + 1]));
                        if (a >= 0 && a < actin.n && m >= 0 && m < myosin.n) {
                            am_bonds[a][m] = 1;
                        }
                    }
                }
            } catch (H5::Exception&) {
                // Keep zeros if legacy AM bonds are unavailable.
            }
        }

        for (int i = 0; i < actin.n; ++i) {
            std::fill(am_state[i].begin(), am_state[i].end(), 0);
            std::fill(am_seen[i].begin(), am_seen[i].end(), 0);
            std::fill(am_target_axial[i].begin(), am_target_axial[i].end(), 0.0);
        }
        if (dynamic_am) {
            _initialize_am_head_pools();
            std::vector<double> material_e1_values;
            std::vector<double> material_e2_values;
            if (!load_state_tensor_double(
                    group_state, "myosin_material_e1_current",
                    static_cast<size_t>(myosin.n), 3, material_e1_values) ||
                !load_state_tensor_double(
                    group_state, "myosin_material_e2_current",
                    static_cast<size_t>(myosin.n), 3, material_e2_values)) {
                throw std::runtime_error(
                    "Dynamic AM schema v8 checkpoint is missing material frames");
            }
            for (int m = 0; m < myosin.n; ++m) {
                myosin_material_e1[static_cast<size_t>(m)] = {
                    material_e1_values[3 * m], material_e1_values[3 * m + 1],
                    material_e1_values[3 * m + 2]};
                myosin_material_e2[static_cast<size_t>(m)] = {
                    material_e2_values[3 * m], material_e2_values[3 * m + 1],
                    material_e2_values[3 * m + 2]};
                const vec axis = myosin.direction[m];
                transport_myosin_material_frame(m, axis, axis, 0.0);
            }
            std::vector<int> am_head_states;
            std::vector<int> am_head_bound_actins;
            std::vector<double> am_head_targets;
            std::vector<double> am_head_binding_azimuths;
            std::vector<double> am_head_stroke_progress;
            const size_t heads_per_myosin = static_cast<size_t>(2 * am_heads_per_half);
            const bool have_head_state =
                load_state_tensor_int(group_state, "am_head_state_current",
                                      static_cast<size_t>(myosin.n), heads_per_myosin,
                                      am_head_states) &&
                load_state_tensor_int(group_state, "am_head_bound_actin_current",
                                      static_cast<size_t>(myosin.n), heads_per_myosin,
                                      am_head_bound_actins) &&
                load_state_tensor_double(group_state, "am_head_target_axial_current",
                                         static_cast<size_t>(myosin.n), heads_per_myosin,
                                         am_head_targets);
            const bool have_binding_azimuths = load_state_tensor_double(
                group_state, "am_head_binding_azimuth_current",
                static_cast<size_t>(myosin.n), heads_per_myosin,
                am_head_binding_azimuths);
            const bool have_stroke_progress = load_state_tensor_double(
                group_state, "am_head_stroke_progress_current",
                static_cast<size_t>(myosin.n), heads_per_myosin,
                am_head_stroke_progress);
            if (!have_head_state) {
                throw std::runtime_error(
                    "Dynamic AM schema v5 checkpoint is missing myosin-owned head state");
            }
            for (int a = 0; a < actin.n; ++a) {
                std::fill(am_bonds[a].begin(), am_bonds[a].end(), 0);
            }
            for (int m = 0; m < myosin.n; ++m) {
                for (int head_id = 0; head_id < static_cast<int>(heads_per_myosin); ++head_id) {
                    const size_t flat_index = static_cast<size_t>(m) * heads_per_myosin +
                                              static_cast<size_t>(head_id);
                    auto* head = _am_head(m, head_id);
                    if (head == nullptr) continue;
                    head->state = static_cast<BondFlag>(
                        std::clamp(am_head_states[flat_index], 0, 3));
                    head->bound_actin = am_head_bound_actins[flat_index];
                    head->target_axial = am_head_targets[flat_index];
                    if (have_binding_azimuths) {
                        head->binding_azimuth =
                            am_head_binding_azimuths[flat_index];
                    }
                    head->stroke_progress = have_stroke_progress
                        ? std::clamp(am_head_stroke_progress[flat_index],
                                     0.0, am_power_stroke)
                        : (head->state == static_cast<BondFlag>(
                               am_motor_kinetics::State::strong)
                               ? am_power_stroke : 0.0);
                    if (head->state == static_cast<BondFlag>(
                            am_motor_kinetics::State::detached) ||
                        head->state == static_cast<BondFlag>(
                            am_motor_kinetics::State::off_srx)) {
                        head->bound_actin = -1;
                    } else if (head->bound_actin < 0 || head->bound_actin >= actin.n) {
                        throw std::runtime_error(
                            "Dynamic AM checkpoint contains an invalid bound-actin index");
                    } else {
                        am_bonds[head->bound_actin][m] = 1;
                    }
                }
            }
            for (int a = 0; a < actin.n; ++a) {
                for (int m = 0; m < myosin.n; ++m) {
                    if (am_bonds[a][m] != 0) _sync_am_pair_summary(a, m);
                }
            }
            std::vector<double> restored_q;
            if (!load_state_vector_double(
                    group_state, "am_q_current",
                    static_cast<size_t>(actin.n), restored_q)) {
                throw std::runtime_error(
                    "Dynamic AM schema v4 checkpoint is missing thin activation q_i");
            }
            am_thin_activation = std::move(restored_q);
        }
        {
            std::vector<int> restored_turnover_steps;
            if (load_state_vector(
                    group_state, "actin_last_turnover_step_current",
                    static_cast<size_t>(actin.n), restored_turnover_steps)) {
                actin_last_turnover_step.resize(static_cast<size_t>(actin.n));
                for (int i = 0; i < actin.n; ++i) {
                    actin_last_turnover_step[static_cast<size_t>(i)] =
                        static_cast<RecoveryStep>(std::max(
                            0, restored_turnover_steps[static_cast<size_t>(i)]));
                }
            } else if (actin_turnover_enabled) {
                throw std::runtime_error(
                    "Turnover checkpoint is missing actin renewal ages");
            } else {
                actin_last_turnover_step.assign(
                    static_cast<size_t>(actin.n), 0);
            }
        }

        if (titin_mode == TitinMode::explicit_dynamic) {
            const size_t titin_count = explicit_titins.size();
            std::vector<int> bound_actin;
            std::vector<int> last_bound_actin;
            std::vector<int> attach_step;
            std::vector<int> last_detach_step;
            std::vector<double> distance;
            std::vector<double> extension;
            std::vector<double> force;
            std::vector<double> attachment_distance;
            std::vector<double> attachment_energy;
            std::vector<double> attachment_force;
            const bool have_titin_state =
                load_state_vector(group_state, "titin_bound_actin_current",
                                  titin_count, bound_actin) &&
                load_state_vector(group_state,
                                  "titin_last_bound_actin_current",
                                  titin_count, last_bound_actin) &&
                load_state_vector(group_state, "titin_attach_step_current",
                                  titin_count, attach_step) &&
                load_state_vector(group_state,
                                  "titin_last_detach_step_current",
                                  titin_count, last_detach_step) &&
                load_state_vector_double(group_state, "titin_distance_current",
                                         titin_count, distance) &&
                load_state_vector_double(group_state, "titin_extension_current",
                                         titin_count, extension) &&
                load_state_vector_double(group_state, "titin_force_current",
                                         titin_count, force) &&
                load_state_vector_double(
                    group_state, "titin_attachment_distance_current",
                    titin_count, attachment_distance) &&
                load_state_vector_double(
                    group_state, "titin_attachment_energy_current",
                    titin_count, attachment_energy) &&
                load_state_vector_double(
                    group_state, "titin_attachment_force_current",
                    titin_count, attachment_force);
            if (!have_titin_state) {
                throw std::runtime_error(
                    "Explicit titin checkpoint is missing molecular state");
            }
            for (size_t index = 0; index < titin_count; ++index) {
                auto& titin = explicit_titins[index];
                titin.bound_actin = bound_actin[index];
                titin.bound = titin.bound_actin >= 0;
                titin.last_bound_actin = last_bound_actin[index];
                titin.attach_step = static_cast<RecoveryStep>(
                    std::max(0, attach_step[index]));
                titin.last_detach_step = static_cast<RecoveryStep>(
                    std::max(0, last_detach_step[index]));
                titin.distance = distance[index];
                titin.extension = extension[index];
                titin.force = force[index];
                titin.attachment_distance = attachment_distance[index];
                titin.attachment_energy = attachment_energy[index];
                titin.attachment_force = attachment_force[index];
                titin.candidate_actin = -1;
                titin.candidate_distance = 0.0;
                titin.candidate_on_rate = 0.0;
            }
            rebuild_explicit_titin_ownership();
        }

        myomesin_bridges.clear();
        if (k_mm > 0.0) {
            if (!file.nameExists("/myomesin/frame_offsets") ||
                !file.nameExists("/myomesin/bridge_state")) {
                throw std::runtime_error(
                    "Sparse myomesin checkpoint is missing molecular state");
            }
            H5::Group myomesin_group(file.openGroup("/myomesin"));
            std::vector<hsize_t> offset_dims;
            const std::vector<double> offsets = load_from_dataset(
                myomesin_group, "frame_offsets", offset_dims);
            if (offset_dims.size() != 2 || offset_dims[1] != 3 ||
                target_frame < 0 ||
                target_frame >= static_cast<int>(offset_dims[0])) {
                throw std::runtime_error(
                    "Sparse myomesin frame offsets do not match trajectory");
            }
            const size_t offset_base =
                static_cast<size_t>(target_frame) * 3;
            const size_t start = static_cast<size_t>(std::max(
                0.0, offsets[offset_base + 1]));
            const size_t count = static_cast<size_t>(std::max(
                0.0, offsets[offset_base + 2]));
            std::vector<hsize_t> state_dims;
            const std::vector<double> states = load_from_dataset(
                myomesin_group, "bridge_state", state_dims);
            if (state_dims.size() != 2 || state_dims[1] != 8 ||
                start + count > static_cast<size_t>(state_dims[0])) {
                throw std::runtime_error(
                    "Sparse myomesin frame contains invalid state bounds");
            }
            std::vector<int> degree(static_cast<size_t>(myosin.n), 0);
            std::unordered_set<std::uint64_t> pairs;
            std::uint64_t largest_id = 0;
            for (size_t record = 0; record < count; ++record) {
                const size_t base = 8 * (start + record);
                MyomesinBridge bridge;
                bridge.id =
                    (static_cast<std::uint64_t>(static_cast<std::uint32_t>(
                         std::llround(states[base]))) << 32) |
                    static_cast<std::uint32_t>(
                        std::llround(states[base + 1]));
                bridge.myosin_i = static_cast<int>(
                    std::llround(states[base + 2]));
                bridge.myosin_j = static_cast<int>(
                    std::llround(states[base + 3]));
                bridge.attach_step = static_cast<RecoveryStep>(std::max(
                    0.0, states[base + 4]));
                bridge.axial_offset = states[base + 5];
                bridge.radial_distance = states[base + 6];
                bridge.force = states[base + 7];
                if (bridge.id == 0 || bridge.myosin_i < 0 ||
                    bridge.myosin_i >= myosin.n ||
                    bridge.myosin_j <= bridge.myosin_i ||
                    bridge.myosin_j >= myosin.n) {
                    throw std::runtime_error(
                        "Sparse myomesin checkpoint contains invalid bridge");
                }
                const auto key = counter_rng::canonical_pair_key(
                    bridge.myosin_i, bridge.myosin_j);
                if (!pairs.insert(key).second ||
                    ++degree[static_cast<size_t>(bridge.myosin_i)] >
                        myomesin_sites_per_filament ||
                    ++degree[static_cast<size_t>(bridge.myosin_j)] >
                        myomesin_sites_per_filament) {
                    throw std::runtime_error(
                        "Sparse myomesin checkpoint violates valency");
                }
                largest_id = std::max(largest_id, bridge.id);
                myomesin_bridges.push_back(bridge);
            }
            myomesin_next_id = largest_id + 1;
            if (myomesin_group.nameExists("next_id_by_frame")) {
                std::vector<hsize_t> next_dims;
                const std::vector<double> next_ids = load_from_dataset(
                    myomesin_group, "next_id_by_frame", next_dims);
                if (next_dims.size() != 2 || next_dims[1] != 2 ||
                    target_frame >= static_cast<int>(next_dims[0])) {
                    throw std::runtime_error(
                        "Sparse myomesin next-ID state does not match trajectory");
                }
                const size_t next_base =
                    static_cast<size_t>(target_frame) * 2;
                myomesin_next_id =
                    (static_cast<std::uint64_t>(
                         static_cast<std::uint32_t>(std::llround(
                             next_ids[next_base]))) << 32) |
                    static_cast<std::uint32_t>(std::llround(
                        next_ids[next_base + 1]));
            } else if (myomesin_group.nameExists("events")) {
                // Compatibility with schema-1 files written before the next
                // ID was stored per frame.  Events retain detached molecules
                // that are absent from the occupied-state rows.
                std::vector<hsize_t> event_dims;
                const std::vector<double> events = load_from_dataset(
                    myomesin_group, "events", event_dims);
                if (event_dims.size() == 2 && event_dims[1] == 9) {
                    for (size_t row = 0;
                         row < static_cast<size_t>(event_dims[0]); ++row) {
                        if (events[9 * row] >
                            static_cast<double>(current_step)) {
                            continue;
                        }
                        const auto event_id = static_cast<std::uint64_t>(
                            std::llround(events[9 * row + 2]));
                        myomesin_next_id = std::max(
                            myomesin_next_id, event_id + 1);
                    }
                }
            }
        }

        // Restore per-actin cb_status from the matrix so diagnostic output remains consistent.
        for (int i = 0; i < actin.n; ++i) {
            int max_status = 0;
            for (int j = 0; j < actin.n; ++j) {
                max_status = std::max(max_status, static_cast<int>(actin_actin_status[i][j]));
            }
            actin.cb_status[i] = max_status;
        }

        // Restore RNG states for exact reproducibility.
        if (rng != nullptr) {
            if (group_state.nameExists("rng_main_state")) {
                std::vector<double> rng_main_data = load_from_dataset(group_state, "rng_main_state", dims);
                const size_t expected = gsl_rng_size(rng);
                if (dims.size() >= 2 && target_frame >= 0 && target_frame < static_cast<int>(dims[0]) &&
                    static_cast<size_t>(dims[1]) == expected) {
                    size_t start = static_cast<size_t>(target_frame) * expected;
                    auto* state_ptr = static_cast<unsigned char*>(gsl_rng_state(rng));
                    for (size_t idx = 0; idx < expected; ++idx) {
                        int byte_value = static_cast<int>(std::llround(rng_main_data[start + idx]));
                        byte_value = std::clamp(byte_value, 0, 255);
                        state_ptr[idx] = static_cast<unsigned char>(byte_value);
                    }
                }
            }
        }
        if (!rng_engines.empty() && rng_engines[0] != nullptr &&
            group_state.nameExists("rng_thread_state")) {
                std::vector<double> thread_data = load_from_dataset(group_state, "rng_thread_state", dims);
                if (dims.size() >= 2 && target_frame >= 0 && target_frame < static_cast<int>(dims[0])) {
                    int file_thread_count = 0;
                    int file_state_size = 0;
                    if (group_state.nameExists("rng_thread_count") && group_state.nameExists("rng_thread_state_size")) {
                        std::vector<hsize_t> scalar_dims;
                        std::vector<double> tc = load_from_dataset(group_state, "rng_thread_count", scalar_dims);
                        if (scalar_dims.size() >= 2 && target_frame < static_cast<int>(scalar_dims[0])) {
                            file_thread_count = static_cast<int>(std::llround(tc[target_frame]));
                        }
                        std::vector<double> ts = load_from_dataset(group_state, "rng_thread_state_size", scalar_dims);
                        if (scalar_dims.size() >= 2 && target_frame < static_cast<int>(scalar_dims[0])) {
                            file_state_size = static_cast<int>(std::llround(ts[target_frame]));
                        }
                    }

                    if (file_thread_count <= 0) {
                        file_thread_count = static_cast<int>(rng_engines.size());
                    }
                    if (file_state_size <= 0) {
                        if (file_thread_count > 0) {
                            file_state_size = static_cast<int>(dims[1] / static_cast<hsize_t>(file_thread_count));
                        } else {
                            file_state_size = 0;
                        }
                    }

                    const size_t frame_width = static_cast<size_t>(dims[1]);
                    const size_t frame_start = static_cast<size_t>(target_frame) * frame_width;
                    const int local_thread_count = static_cast<int>(rng_engines.size());
                    const int restore_threads = std::min(local_thread_count, file_thread_count);
                    for (int t = 0; t < restore_threads; ++t) {
                        if (rng_engines[t] == nullptr) {
                            continue;
                        }
                        const size_t local_state_size = gsl_rng_size(rng_engines[t]);
                        const size_t copy_size = std::min(local_state_size, static_cast<size_t>(file_state_size));
                        auto* local_ptr = static_cast<unsigned char*>(gsl_rng_state(rng_engines[t]));
                        size_t thread_offset = frame_start + static_cast<size_t>(t) * static_cast<size_t>(file_state_size);
                        if (thread_offset + copy_size > thread_data.size()) {
                            break;
                        }
                        for (size_t idx = 0; idx < copy_size; ++idx) {
                            int byte_value = static_cast<int>(std::llround(thread_data[thread_offset + idx]));
                            byte_value = std::clamp(byte_value, 0, 255);
                            local_ptr[idx] = static_cast<unsigned char>(byte_value);
                        }
                    }
                    if (local_thread_count > file_thread_count) {
                        for (int t = file_thread_count; t < local_thread_count; ++t) {
                            if (rng_engines[t] == nullptr) {
                                continue;
                            }
                            gsl_rng_set(
                                rng_engines[t],
                                static_cast<unsigned long>(initial_seed) +
                                    static_cast<unsigned long>(t) +
                                    static_cast<unsigned long>(current_step));
                        }
                    }
                }
        }

        std::vector<double> neighbor_last_actin_x;
        std::vector<double> neighbor_last_actin_y;
        std::vector<double> neighbor_last_actin_z;
        std::vector<double> neighbor_last_myosin_x;
        std::vector<double> neighbor_last_myosin_y;
        std::vector<double> neighbor_last_myosin_z;
        bool have_neighbor_cache =
            load_state_vector_double(group_state, "neighbor_last_actin_x", static_cast<size_t>(actin.n), neighbor_last_actin_x) &&
            load_state_vector_double(group_state, "neighbor_last_actin_y", static_cast<size_t>(actin.n), neighbor_last_actin_y) &&
            load_state_vector_double(group_state, "neighbor_last_actin_z", static_cast<size_t>(actin.n), neighbor_last_actin_z) &&
            load_state_vector_double(group_state, "neighbor_last_myosin_x", static_cast<size_t>(myosin.n), neighbor_last_myosin_x) &&
            load_state_vector_double(group_state, "neighbor_last_myosin_y", static_cast<size_t>(myosin.n), neighbor_last_myosin_y) &&
            load_state_vector_double(group_state, "neighbor_last_myosin_z", static_cast<size_t>(myosin.n), neighbor_last_myosin_z);
        if (have_neighbor_cache) {
            // Reconstruct cached neighbor pairs from the original rebuild reference positions.
            neighbor_list.set_species_positions(neighbor_last_actin_x, neighbor_last_actin_y, neighbor_last_actin_z,
                                                neighbor_last_myosin_x, neighbor_last_myosin_y, neighbor_last_myosin_z);
            neighbor_list.rebuild_neighbor_list();
            // Then restore current particle positions while preserving the cached "last" positions.
            neighbor_list.set_species_positions(actin.center_x, actin.center_y, actin.center_z,
                                                myosin.center_x, myosin.center_y, myosin.center_z);
            restored_neighbor_cache = true;
        }

        if (frame_index < 0) {
            const std::string resume_group_path =
                file.nameExists("/resume") ? "/resume" :
                (file.nameExists("/resume_prev") ? "/resume_prev" : "");
            if (!resume_group_path.empty()) {
                try {
                    if (resume_group_path == "/resume_prev") {
                    }
                    H5::Group group_resume(file.openGroup(resume_group_path));
                    std::vector<int> ints;
                    std::vector<double> doubles;
                    size_t resume_step_value = current_step;

                    if (load_fixed_vector_int(group_resume, "current_step", 1, ints)) {
                        const size_t resume_step = static_cast<size_t>(std::max(0, ints[0]));
                        if (resume_step < current_step) {
                            throw H5::Exception("Sarcomere::load_state", "Stale resume snapshot");
                        }
                        current_step = resume_step;
                        resume_step_value = resume_step;
                    }
                    if (load_fixed_matrix_double(group_resume, "actin_center", static_cast<size_t>(actin.n), 3, doubles)) {
                        for (int i = 0; i < actin.n; ++i) {
                            actin.center[i].x = doubles[3 * i];
                            actin.center[i].y = doubles[3 * i + 1];
                            actin.center[i].z = doubles[3 * i + 2];
                        }
                    }
                    if (load_fixed_matrix_double(group_resume, "actin_direction", static_cast<size_t>(actin.n), 3, doubles)) {
                        for (int i = 0; i < actin.n; ++i) {
                            actin.direction[i].x = doubles[3 * i];
                            actin.direction[i].y = doubles[3 * i + 1];
                            actin.direction[i].z = doubles[3 * i + 2];
                        }
                    }
                    if (load_fixed_matrix_double(group_resume, "myosin_center", static_cast<size_t>(myosin.n), 3, doubles)) {
                        for (int i = 0; i < myosin.n; ++i) {
                            myosin.center[i].x = doubles[3 * i];
                            myosin.center[i].y = doubles[3 * i + 1];
                            myosin.center[i].z = doubles[3 * i + 2];
                        }
                    }
                    if (load_fixed_matrix_double(group_resume, "myosin_direction", static_cast<size_t>(myosin.n), 3, doubles)) {
                        for (int i = 0; i < myosin.n; ++i) {
                            myosin.direction[i].x = doubles[3 * i];
                            myosin.direction[i].y = doubles[3 * i + 1];
                            myosin.direction[i].z = doubles[3 * i + 2];
                        }
                    }
                    std::vector<double> material_e1_values;
                    std::vector<double> material_e2_values;
                    if (!load_fixed_matrix_double(
                            group_resume, "myosin_material_e1_current",
                            static_cast<size_t>(myosin.n), 3,
                            material_e1_values) ||
                        !load_fixed_matrix_double(
                            group_resume, "myosin_material_e2_current",
                            static_cast<size_t>(myosin.n), 3,
                            material_e2_values)) {
                        throw std::runtime_error(
                            "Resume snapshot is missing transported myosin frames");
                    }
                    for (int m = 0; m < myosin.n; ++m) {
                        myosin_material_e1[static_cast<size_t>(m)] = {
                            material_e1_values[3 * m],
                            material_e1_values[3 * m + 1],
                            material_e1_values[3 * m + 2]};
                        myosin_material_e2[static_cast<size_t>(m)] = {
                            material_e2_values[3 * m],
                            material_e2_values[3 * m + 1],
                            material_e2_values[3 * m + 2]};
                        const vec axis = myosin.direction[m];
                        transport_myosin_material_frame(m, axis, axis, 0.0);
                    }
                    actin.update_endpoints();
                    myosin.update_endpoints();

                    restore_aa_state_from_fixed_pairs(
                        group_resume,
                        "aa_pairs_prev_count",
                        "aa_pairs_prev",
                        "aa_status_prev",
                        "aa_lifetime_prev",
                        actin_actin_bonds_prev,
                        actin_actin_status_prev,
                        actin_actin_lifetime_prev);
                    restore_aa_state_from_fixed_pairs(
                        group_resume,
                        "aa_pairs_current_count",
                        "aa_pairs_current",
                        "aa_status_current",
                        "aa_lifetime_current",
                        actin_actin_bonds,
                        actin_actin_status,
                        actin_actin_lifetime);
                    for (int i = 0; i < actin.n; ++i) {
                        std::fill(aa_attach_step[i].begin(), aa_attach_step[i].end(), -1);
                    }
                    {
                        std::vector<int> counts_attach;
                        std::vector<int> pairs_attach;
                        std::vector<int> attach_values;
                        const size_t max_aa_pairs = static_cast<size_t>(actin.n) *
                                                    static_cast<size_t>(std::max(10, 2 * max_myosin_bonds));
                        if (load_fixed_vector_int(group_resume, "aa_pairs_current_count", 1, counts_attach) &&
                            load_fixed_matrix_int(group_resume, "aa_pairs_current", max_aa_pairs, 2, pairs_attach) &&
                            load_fixed_matrix_int(group_resume, "aa_attach_step_current", max_aa_pairs, 1, attach_values)) {
                            const size_t count = static_cast<size_t>(std::max(0, counts_attach[0]));
                            for (size_t idx = 0; idx < count; ++idx) {
                                const int a = pairs_attach[2 * idx];
                                const int b = pairs_attach[2 * idx + 1];
                                if (a < 0 || a >= actin.n || b < 0 || b >= actin.n || a == b) {
                                    continue;
                                }
                                aa_attach_step[a][b] = attach_values[idx];
                                aa_attach_step[b][a] = attach_values[idx];
                            }
                        }
                    }
                    for (int i = 0; i < actin.n; ++i) {
                        std::fill(aa_status2_memory_until[i].begin(),
                                  aa_status2_memory_until[i].end(), 0);
                    }
                    {
                        std::vector<int> counts_memory;
                        std::vector<int> pairs_memory;
                        std::vector<int> memory_values;
                        const size_t max_aa_pairs = static_cast<size_t>(actin.n) *
                                                    static_cast<size_t>(std::max(10, 2 * max_myosin_bonds));
                        if (load_fixed_vector_int(group_resume, "aa_pairs_current_count", 1, counts_memory) &&
                            load_fixed_matrix_int(group_resume, "aa_pairs_current", max_aa_pairs, 2, pairs_memory) &&
                            load_fixed_matrix_int(group_resume, "aa_status2_memory_until_current",
                                                  max_aa_pairs, 1, memory_values)) {
                            const size_t count = static_cast<size_t>(std::max(0, counts_memory[0]));
                            for (size_t idx = 0; idx < count; ++idx) {
                                const int a = pairs_memory[2 * idx];
                                const int b = pairs_memory[2 * idx + 1];
                                if (a < 0 || a >= actin.n || b < 0 || b >= actin.n || a == b) continue;
                                const RecoveryStep until = static_cast<RecoveryStep>(
                                    std::max(0, memory_values[idx]));
                                aa_status2_memory_until[a][b] = until;
                                aa_status2_memory_until[b][a] = until;
                            }
                        }
                    }

                    aa_links.clear();
                    aa_free_link_slots.clear();
                    aa_active_link_slots.clear();
                    aa_links_by_pair.clear();
                    aa_links_by_actin.assign(static_cast<size_t>(actin.n), {});
                    aa_next_link_id = 1;
                    if (aa_local_multilink) {
                        std::vector<int> local_counts;
                        std::vector<int> local_identity;
                        std::vector<double> local_coordinates;
                        std::vector<int> local_motor_identity;
                        std::vector<double> local_motor_evidence;
                        std::vector<int> local_arm_state;
                        std::vector<double> local_arm_loads;
                        std::vector<int> local_next_id;
                        if (!load_fixed_vector_int(
                                group_resume, "aa_local_link_count", 1,
                                local_counts)) {
                            throw std::runtime_error(
                                "Resume snapshot is missing local ACTN2 link count");
                        }
                        const size_t local_count = static_cast<size_t>(
                            std::max(0, local_counts[0]));
                        const size_t local_rows = std::max<size_t>(1, local_count);
                        if (!load_fixed_matrix_int(
                                group_resume, "aa_local_link_identity",
                                local_rows, 7, local_identity) ||
                            !load_fixed_matrix_double(
                                group_resume, "aa_local_link_coordinates",
                                local_rows, 2, local_coordinates) ||
                            !load_fixed_vector_int(
                                group_resume, "aa_local_next_id", 2,
                                local_next_id)) {
                            throw std::runtime_error(
                                "Resume snapshot is missing local ACTN2 link state");
                        }
                        const bool needs_motor_identity =
                            aa_status_mode == actn2_load_kinetics::StatusMode::
                                motor_transmitted_register;
                        const bool has_motor_identity =
                            load_fixed_matrix_int(
                                group_resume,
                                "aa_local_link_motor_identity",
                                local_rows, 12, local_motor_identity) &&
                            load_fixed_matrix_double(
                                group_resume,
                                "aa_local_link_motor_evidence",
                                local_rows, 7, local_motor_evidence);
                        if (needs_motor_identity && !has_motor_identity) {
                            throw std::runtime_error(
                                "Resume snapshot is missing identity-resolved "
                                "ACTN2 motor-pair evidence");
                        }
                        const bool has_arm_state =
                            load_fixed_matrix_int(
                                group_resume, "aa_local_link_arm_state",
                                local_rows, 3, local_arm_state) &&
                            load_fixed_matrix_double(
                                group_resume, "aa_local_link_arm_loads",
                                local_rows, 4, local_arm_loads);
                        if (aa_arm_state_model && !has_arm_state) {
                            throw std::runtime_error(
                                "Resume snapshot is missing ACTN2 two-arm state");
                        }
                        aa_links.reserve(local_count);
                        for (size_t record = 0; record < local_count; ++record) {
                            AACrosslink link;
                            link.id =
                                (static_cast<std::uint64_t>(static_cast<std::uint32_t>(
                                     local_identity[7 * record])) << 32) |
                                static_cast<std::uint32_t>(
                                    local_identity[7 * record + 1]);
                            link.actin_i = local_identity[7 * record + 2];
                            link.actin_j = local_identity[7 * record + 3];
                            const int status = local_identity[7 * record + 4];
                            if (link.actin_i < 0 || link.actin_i >= actin.n ||
                                link.actin_j < 0 || link.actin_j >= actin.n ||
                                link.actin_i == link.actin_j ||
                                (status != 1 && status != 2)) {
                                throw std::runtime_error(
                                    "Invalid local ACTN2 link in resume snapshot");
                            }
                            link.status = static_cast<BondFlag>(status);
                            link.attach_step = static_cast<RecoveryStep>(
                                static_cast<std::uint32_t>(
                                    local_identity[7 * record + 5]));
                            link.status2_memory_until = static_cast<RecoveryStep>(
                                static_cast<std::uint32_t>(
                                    local_identity[7 * record + 6]));
                            link.s_i = local_coordinates[2 * record];
                            link.s_j = local_coordinates[2 * record + 1];
                            if (has_arm_state) {
                                const int arm_state_value =
                                    local_arm_state[3 * record];
                                if (arm_state_value < static_cast<int>(
                                        actn2_arm_kinetics::State::weak_bridge) ||
                                    arm_state_value > static_cast<int>(
                                        actn2_arm_kinetics::State::tethered_to_j)) {
                                    throw std::runtime_error(
                                        "Resume snapshot contains invalid ACTN2 arm state");
                                }
                                link.arm_state = static_cast<
                                    actn2_arm_kinetics::State>(arm_state_value);
                                link.arm_state_since = static_cast<RecoveryStep>(
                                    std::max(0, local_arm_state[3 * record + 1]));
                                link.last_partner =
                                    local_arm_state[3 * record + 2];
                                link.arm_force_i = local_arm_loads[4 * record];
                                link.arm_force_j =
                                    local_arm_loads[4 * record + 1];
                                link.contact_force =
                                    local_arm_loads[4 * record + 2];
                                link.contact_penetration =
                                    local_arm_loads[4 * record + 3];
                            } else {
                                link.arm_state =
                                    actn2_arm_kinetics::bridge_state_for_status(
                                        status);
                                link.arm_state_since = link.attach_step;
                            }
                            if (has_motor_identity) {
                                link.candidate_myosin_i =
                                    local_motor_identity[12 * record];
                                link.candidate_myosin_j =
                                    local_motor_identity[12 * record + 1];
                                link.candidate_half_i = static_cast<int8_t>(
                                    local_motor_identity[12 * record + 2]);
                                link.candidate_half_j = static_cast<int8_t>(
                                    local_motor_identity[12 * record + 3]);
                                link.candidate_since =
                                    static_cast<RecoveryStep>(
                                        static_cast<std::uint32_t>(
                                            local_motor_identity[
                                                12 * record + 4]));
                                link.candidate_last_seen =
                                    static_cast<RecoveryStep>(
                                        static_cast<std::uint32_t>(
                                            local_motor_identity[
                                                12 * record + 5]));
                                link.candidate_locked =
                                    local_motor_identity[12 * record + 6] != 0;
                                link.candidate_support_ema =
                                    local_motor_evidence[7 * record];
                                link.candidate_support_instantaneous =
                                    local_motor_evidence[7 * record + 1];
                                link.candidate_raw_support_i =
                                    local_motor_evidence[7 * record + 2];
                                link.candidate_raw_support_j =
                                    local_motor_evidence[7 * record + 3];
                                link.candidate_allocated_support_i =
                                    local_motor_evidence[7 * record + 4];
                                link.candidate_allocated_support_j =
                                    local_motor_evidence[7 * record + 5];
                                link.challenger_myosin_i =
                                    local_motor_identity[12 * record + 7];
                                link.challenger_myosin_j =
                                    local_motor_identity[12 * record + 8];
                                link.challenger_half_i = static_cast<int8_t>(
                                    local_motor_identity[12 * record + 9]);
                                link.challenger_half_j = static_cast<int8_t>(
                                    local_motor_identity[12 * record + 10]);
                                link.challenger_since =
                                    static_cast<RecoveryStep>(
                                        static_cast<std::uint32_t>(
                                            local_motor_identity[
                                                12 * record + 11]));
                                link.challenger_support_ema =
                                    local_motor_evidence[7 * record + 6];
                            }
                            link.active = true;
                            aa_links.push_back(link);
                        }
                        aa_next_link_id =
                            (static_cast<std::uint64_t>(static_cast<std::uint32_t>(
                                 local_next_id[0])) << 32) |
                            static_cast<std::uint32_t>(local_next_id[1]);
                        _rebuild_local_aa_indices();
                    }

                    restore_am_state_from_fixed_pairs(
                        group_resume, "am_pairs_prev_count", "am_pairs_prev", am_bonds_prev);
                    restore_am_state_from_fixed_pairs(
                        group_resume, "am_pairs_current_count", "am_pairs_current", am_bonds);

                    if (dynamic_am) {
                        for (int a = 0; a < actin.n; ++a) {
                            std::fill(am_state[a].begin(), am_state[a].end(), 0);
                            std::fill(am_seen[a].begin(), am_seen[a].end(), 0);
                            std::fill(am_target_axial[a].begin(), am_target_axial[a].end(), 0.0);
                        }
                        _initialize_am_head_pools();
                        std::vector<int> am_head_states;
                        std::vector<int> am_head_bound_actins;
                        std::vector<double> am_head_targets;
                        std::vector<double> am_head_binding_azimuths;
                        std::vector<double> am_head_stroke_progress;
                        const size_t heads_per_myosin =
                            static_cast<size_t>(2 * am_heads_per_half);
                        const bool have_head_state =
                            load_fixed_matrix_int(
                                group_resume, "am_head_state_current",
                                static_cast<size_t>(myosin.n), heads_per_myosin,
                                am_head_states) &&
                            load_fixed_matrix_int(
                                group_resume, "am_head_bound_actin_current",
                                static_cast<size_t>(myosin.n), heads_per_myosin,
                                am_head_bound_actins) &&
                            load_fixed_matrix_double(
                                group_resume, "am_head_target_axial_current",
                                static_cast<size_t>(myosin.n), heads_per_myosin,
                                am_head_targets);
                        const bool have_binding_azimuths =
                            load_fixed_matrix_double(
                                group_resume,
                                "am_head_binding_azimuth_current",
                                static_cast<size_t>(myosin.n),
                                heads_per_myosin,
                                am_head_binding_azimuths);
                        const bool have_stroke_progress =
                            load_fixed_matrix_double(
                                group_resume,
                                "am_head_stroke_progress_current",
                                static_cast<size_t>(myosin.n),
                                heads_per_myosin,
                                am_head_stroke_progress);
                        if (!have_head_state) {
                            throw std::runtime_error(
                                "Resume snapshot is missing myosin-owned AM head state");
                        }
                        for (int a = 0; a < actin.n; ++a) {
                            std::fill(am_bonds[a].begin(), am_bonds[a].end(), 0);
                        }
                        for (int m = 0; m < myosin.n; ++m) {
                            for (int head_id = 0;
                                 head_id < static_cast<int>(heads_per_myosin); ++head_id) {
                                const size_t flat_index =
                                    static_cast<size_t>(m) * heads_per_myosin +
                                    static_cast<size_t>(head_id);
                                auto* head = _am_head(m, head_id);
                                if (head == nullptr) continue;
                                head->state = static_cast<BondFlag>(
                                    std::clamp(am_head_states[flat_index], 0, 3));
                                head->bound_actin = am_head_bound_actins[flat_index];
                                head->target_axial = am_head_targets[flat_index];
                                if (have_binding_azimuths) {
                                    head->binding_azimuth =
                                        am_head_binding_azimuths[flat_index];
                                }
                                head->stroke_progress = have_stroke_progress
                                    ? std::clamp(
                                          am_head_stroke_progress[flat_index],
                                          0.0, am_power_stroke)
                                    : (head->state == static_cast<BondFlag>(
                                           am_motor_kinetics::State::strong)
                                           ? am_power_stroke : 0.0);
                                if (head->state == static_cast<BondFlag>(
                                        am_motor_kinetics::State::detached) ||
                                    head->state == static_cast<BondFlag>(
                                        am_motor_kinetics::State::off_srx)) {
                                    head->bound_actin = -1;
                                } else if (head->bound_actin < 0 ||
                                           head->bound_actin >= actin.n) {
                                    throw std::runtime_error(
                                        "Resume snapshot contains invalid AM head ownership");
                                } else {
                                    am_bonds[head->bound_actin][m] = 1;
                                }
                            }
                        }
                        for (int a = 0; a < actin.n; ++a) {
                            for (int m = 0; m < myosin.n; ++m) {
                                if (am_bonds[a][m] != 0) {
                                    _sync_am_pair_summary(a, m);
                                }
                            }
                        }
                        std::vector<double> restored_q;
                        if (!load_fixed_vector_double(
                                group_resume, "am_q_current",
                                static_cast<size_t>(actin.n), restored_q)) {
                            throw std::runtime_error(
                                "Resume snapshot is missing thin activation q_i");
                        }
                        am_thin_activation = std::move(restored_q);
                    }
                    {
                        std::vector<int> restored_turnover_steps;
                        if (load_fixed_vector_int(
                                group_resume,
                                "actin_last_turnover_step_current",
                                static_cast<size_t>(actin.n),
                                restored_turnover_steps)) {
                            actin_last_turnover_step.resize(
                                static_cast<size_t>(actin.n));
                            for (int i = 0; i < actin.n; ++i) {
                                actin_last_turnover_step[static_cast<size_t>(i)] =
                                    static_cast<RecoveryStep>(std::max(
                                        0, restored_turnover_steps[
                                               static_cast<size_t>(i)]));
                            }
                        } else if (actin_turnover_enabled) {
                            throw std::runtime_error(
                                "Resume snapshot is missing actin renewal ages");
                        } else {
                            actin_last_turnover_step.assign(
                                static_cast<size_t>(actin.n), 0);
                        }
                    }

                    if (titin_mode == TitinMode::explicit_dynamic) {
                        const size_t titin_count = explicit_titins.size();
                        std::vector<int> titin_int;
                        std::vector<double> titin_double;
                        if (!load_fixed_matrix_int(
                                group_resume, "titin_state_int", titin_count,
                                5, titin_int) ||
                            !load_fixed_matrix_double(
                                group_resume, "titin_state_double", titin_count,
                                6, titin_double)) {
                            throw std::runtime_error(
                                "Resume snapshot is missing explicit titin molecular state");
                        }
                        for (size_t index = 0; index < titin_count; ++index) {
                            auto& titin = explicit_titins[index];
                            titin.bound = titin_int[5 * index] != 0;
                            titin.bound_actin = titin.bound
                                ? titin_int[5 * index + 1] : -1;
                            titin.last_bound_actin =
                                titin_int[5 * index + 2];
                            titin.attach_step = static_cast<RecoveryStep>(
                                std::max(0, titin_int[5 * index + 3]));
                            titin.last_detach_step = static_cast<RecoveryStep>(
                                std::max(0, titin_int[5 * index + 4]));
                            titin.distance = titin_double[6 * index];
                            titin.extension = titin_double[6 * index + 1];
                            titin.force = titin_double[6 * index + 2];
                            titin.attachment_distance =
                                titin_double[6 * index + 3];
                            titin.attachment_energy =
                                titin_double[6 * index + 4];
                            titin.attachment_force =
                                titin_double[6 * index + 5];
                            titin.candidate_actin = -1;
                            titin.candidate_distance = 0.0;
                            titin.candidate_on_rate = 0.0;
                        }
                        rebuild_explicit_titin_ownership();
                        // Association provenance is optional on resume: a
                        // snapshot written before schema 3 simply restarts the
                        // counters, which is visible in the series rather than
                        // silently wrong.
                        std::vector<int> association_counts;
                        if (load_fixed_matrix_int(
                                group_resume, "titin_association_counts", 1, 4,
                                association_counts) &&
                            association_counts.size() == 4) {
                            explicit_titin_supported_associations =
                                static_cast<std::uint64_t>(
                                    std::max(0, association_counts[0]));
                            explicit_titin_unsupported_associations =
                                static_cast<std::uint64_t>(
                                    std::max(0, association_counts[1]));
                            explicit_titin_rebind_associations =
                                static_cast<std::uint64_t>(
                                    std::max(0, association_counts[2]));
                            explicit_titin_unsupported_releases =
                                static_cast<std::uint64_t>(
                                    std::max(0, association_counts[3]));
                        }
                    }

                    if (k_mm > 0.0) {
                        const size_t bridge_rows = std::max<size_t>(
                            1, static_cast<size_t>(myosin.n) *
                                   static_cast<size_t>(
                                       myomesin_sites_per_filament) /
                                   2);
                        std::vector<int> bridge_count;
                        std::vector<int> bridge_int;
                        std::vector<double> bridge_double;
                        std::vector<int> next_id;
                        if (!load_fixed_vector_int(
                                group_resume, "myomesin_bridge_count", 1,
                                bridge_count) ||
                            !load_fixed_matrix_int(
                                group_resume, "myomesin_state_int",
                                bridge_rows, 5, bridge_int) ||
                            !load_fixed_matrix_double(
                                group_resume, "myomesin_state_double",
                                bridge_rows, 3, bridge_double) ||
                            !load_fixed_vector_int(
                                group_resume, "myomesin_next_id", 2,
                                next_id)) {
                            throw std::runtime_error(
                                "Resume snapshot is missing sparse myomesin state");
                        }
                        const size_t count = static_cast<size_t>(
                            std::max(0, bridge_count[0]));
                        if (count > bridge_rows) {
                            throw std::runtime_error(
                                "Resume snapshot has invalid myomesin count");
                        }
                        myomesin_bridges.clear();
                        std::vector<int> degree(
                            static_cast<size_t>(myosin.n), 0);
                        std::unordered_set<std::uint64_t> pairs;
                        for (size_t index = 0; index < count; ++index) {
                            MyomesinBridge bridge;
                            bridge.id =
                                (static_cast<std::uint64_t>(
                                     static_cast<std::uint32_t>(
                                         bridge_int[5 * index])) << 32) |
                                static_cast<std::uint32_t>(
                                    bridge_int[5 * index + 1]);
                            bridge.myosin_i = bridge_int[5 * index + 2];
                            bridge.myosin_j = bridge_int[5 * index + 3];
                            bridge.attach_step = static_cast<RecoveryStep>(
                                std::max(0, bridge_int[5 * index + 4]));
                            bridge.axial_offset = bridge_double[3 * index];
                            bridge.radial_distance =
                                bridge_double[3 * index + 1];
                            bridge.force = bridge_double[3 * index + 2];
                            if (bridge.id == 0 || bridge.myosin_i < 0 ||
                                bridge.myosin_i >= myosin.n ||
                                bridge.myosin_j <= bridge.myosin_i ||
                                bridge.myosin_j >= myosin.n) {
                                throw std::runtime_error(
                                    "Resume snapshot contains invalid myomesin bridge");
                            }
                            const auto key = counter_rng::canonical_pair_key(
                                bridge.myosin_i, bridge.myosin_j);
                            if (!pairs.insert(key).second ||
                                ++degree[static_cast<size_t>(bridge.myosin_i)] >
                                    myomesin_sites_per_filament ||
                                ++degree[static_cast<size_t>(bridge.myosin_j)] >
                                    myomesin_sites_per_filament) {
                                throw std::runtime_error(
                                    "Resume snapshot violates myomesin valency");
                            }
                            myomesin_bridges.push_back(bridge);
                        }
                        myomesin_next_id =
                            (static_cast<std::uint64_t>(
                                 static_cast<std::uint32_t>(next_id[0])) << 32) |
                            static_cast<std::uint32_t>(next_id[1]);
                    } else {
                        myomesin_bridges.clear();
                    }

                    if (restore_recovery_state_from_fixed_pairs(
                            group_resume,
                            "aa_recovery_count",
                            "aa_recovery_pairs",
                            "aa_recovery_until_values")) {
                        // loaded sparse recovery state
                    } else if (load_fixed_matrix_int(group_resume, "actin_recovery_until",
                                                     static_cast<size_t>(actin.n),
                                                     static_cast<size_t>(actin.n), ints)) {
                        for (int i = 0; i < actin.n; ++i) {
                            for (int j = 0; j < actin.n; ++j) {
                                actin_recovery_until[i][j] = static_cast<RecoveryStep>(std::max(
                                    0, ints[static_cast<size_t>(i) * actin.n + j]));
                            }
                        }
                    }

                    if (rng != nullptr &&
                        load_fixed_vector_int(group_resume, "rng_main_state", gsl_rng_size(rng), ints)) {
                        auto* state_ptr = static_cast<unsigned char*>(gsl_rng_state(rng));
                        for (size_t idx = 0; idx < ints.size(); ++idx) {
                            state_ptr[idx] = static_cast<unsigned char>(std::clamp(ints[idx], 0, 255));
                        }
                    }
                    if (!rng_engines.empty() && rng_engines[0] != nullptr &&
                        load_fixed_vector_int(group_resume, "rng_thread_state_size", 1, ints)) {
                        const int file_state_size = ints[0];
                        if (load_fixed_vector_int(group_resume, "rng_thread_count", 1, ints)) {
                            const int file_thread_count = ints[0];
                            const size_t flat_width =
                                static_cast<size_t>(std::max(0, file_thread_count)) *
                                static_cast<size_t>(std::max(0, file_state_size));
                            if (load_fixed_vector_int(group_resume, "rng_thread_state", flat_width, ints)) {
                                const int restore_threads =
                                    std::min(static_cast<int>(rng_engines.size()), file_thread_count);
                                for (int t = 0; t < restore_threads; ++t) {
                                    if (rng_engines[t] == nullptr) {
                                        continue;
                                    }
                                    auto* local_ptr =
                                        static_cast<unsigned char*>(gsl_rng_state(rng_engines[t]));
                                    const size_t local_size = gsl_rng_size(rng_engines[t]);
                                    const size_t copy_size =
                                        std::min(local_size, static_cast<size_t>(std::max(0, file_state_size)));
                                    const size_t base =
                                        static_cast<size_t>(t) * static_cast<size_t>(std::max(0, file_state_size));
                                    for (size_t idx = 0; idx < copy_size; ++idx) {
                                        local_ptr[idx] =
                                            static_cast<unsigned char>(std::clamp(ints[base + idx], 0, 255));
                                    }
                                }
                                const int local_thread_count = static_cast<int>(rng_engines.size());
                                if (local_thread_count > file_thread_count) {
                                    for (int t = file_thread_count; t < local_thread_count; ++t) {
                                        if (rng_engines[t] == nullptr) {
                                            continue;
                                        }
                                        gsl_rng_set(
                                            rng_engines[t],
                                            static_cast<unsigned long>(initial_seed) +
                                                static_cast<unsigned long>(t) +
                                                static_cast<unsigned long>(current_step));
                                    }
                                }
                            }
                        }
                    }
                    
                    if (load_fixed_vector_double(group_resume, "neighbor_last_actin_x", static_cast<size_t>(actin.n), neighbor_last_actin_x) &&
                        load_fixed_vector_double(group_resume, "neighbor_last_actin_y", static_cast<size_t>(actin.n), neighbor_last_actin_y) &&
                        load_fixed_vector_double(group_resume, "neighbor_last_actin_z", static_cast<size_t>(actin.n), neighbor_last_actin_z) &&
                        load_fixed_vector_double(group_resume, "neighbor_last_myosin_x", static_cast<size_t>(myosin.n), neighbor_last_myosin_x) &&
                        load_fixed_vector_double(group_resume, "neighbor_last_myosin_y", static_cast<size_t>(myosin.n), neighbor_last_myosin_y) &&
                        load_fixed_vector_double(group_resume, "neighbor_last_myosin_z", static_cast<size_t>(myosin.n), neighbor_last_myosin_z)) {
                        neighbor_list.set_species_positions(neighbor_last_actin_x, neighbor_last_actin_y, neighbor_last_actin_z,
                                                            neighbor_last_myosin_x, neighbor_last_myosin_y, neighbor_last_myosin_z);
                        neighbor_list.rebuild_neighbor_list();
                        neighbor_list.set_species_positions(actin.center_x, actin.center_y, actin.center_z,
                                                            myosin.center_x, myosin.center_y, myosin.center_z);
                        restored_neighbor_cache = true;
                    } else {
                        restored_neighbor_cache = false;
                    }

                    if (!load_any_vector_double(group_resume, "cb_breakage_pending", cb_breakage_events)) {
                        cb_breakage_events.clear();
                    }
                    if (!load_any_vector_double(group_resume, "cb_limit_pending", cb_limit_events)) {
                        cb_limit_events.clear();
                    }
                    if (!load_any_vector_double(group_resume, "aa_completed_lifetimes_pending",
                                                aa_completed_lifetimes)) {
                        aa_completed_lifetimes.clear();
                    }
                    if (!load_any_vector_double(group_resume, "am_transition_events_pending",
                                                am_transition_events)) {
                        am_transition_events.clear();
                    }
                    if (!load_any_vector_double(
                            group_resume, "am_forced_rupture_events_pending",
                            am_forced_rupture_events)) {
                        am_forced_rupture_events.clear();
                    }
                    if (!load_any_vector_double(
                            group_resume, "am_forced_rupture_events_v2_pending",
                            am_forced_rupture_events_v2)) {
                        am_forced_rupture_events_v2.clear();
                    }
                    if (!load_any_vector_double(
                            group_resume, "am_thick_transition_events_pending",
                            am_thick_transition_events)) {
                        am_thick_transition_events.clear();
                    }
                    if (!load_any_vector_double(
                            group_resume, "actin_turnover_events_pending",
                            actin_turnover_events)) {
                        actin_turnover_events.clear();
                    }
                    if (!load_any_vector_double(group_resume, "aa_state_events_pending",
                                                aa_state_events)) {
                        aa_state_events.clear();
                    }
                    if (!load_any_vector_double(group_resume, "aa_load_events_pending",
                                                aa_load_events)) {
                        aa_load_events.clear();
                    }
                    if (!load_any_vector_double(
                            group_resume, "aa_local_link_events_pending",
                            aa_local_link_events)) {
                        aa_local_link_events.clear();
                    }
                    if (!load_any_vector_double(
                            group_resume, "aa_local_link_events_v2_pending",
                            aa_local_link_events_v2)) {
                        aa_local_link_events_v2.clear();
                    }
                    if (!load_any_vector_double(
                            group_resume, "aa_local_link_events_v3_pending",
                            aa_local_link_events_v3)) {
                        aa_local_link_events_v3.clear();
                    }
                    if (!load_any_vector_double(
                            group_resume,
                            "aa_tip_accessibility_events_pending",
                            aa_tip_accessibility_events)) {
                        aa_tip_accessibility_events.clear();
                    }
                    if (!load_any_vector_double(
                            group_resume, "aa_motor_pair_events_pending",
                            aa_motor_pair_events)) {
                        aa_motor_pair_events.clear();
                    }
                    if (!load_any_vector_double(
                            group_resume, "aa_arm_events_pending",
                            aa_arm_events)) {
                        aa_arm_events.clear();
                    }
                    if (!load_any_vector_double(
                            group_resume, "explicit_titin_events_pending",
                            explicit_titin_events)) {
                        explicit_titin_events.clear();
                    }
                    if (!load_any_vector_double(
                            group_resume, "myomesin_events_pending",
                            myomesin_events)) {
                        myomesin_events.clear();
                    }
                    if (resume_step_value < current_step) {
                        throw H5::Exception("Sarcomere::load_state", "Resume snapshot older than committed /state");
                    }
                    loaded_resume_snapshot = true;
                } catch (H5::Exception&) {
                    loaded_resume_snapshot = loaded_state_group;
                }
            } else {
                loaded_resume_snapshot = loaded_state_group;
            }
        }

    } catch (H5::Exception&) {
        loaded_state_group = false;
    }

    if (!loaded_state_group) {
        throw std::runtime_error("Could not load /state group for resume. Refusing approximate fallback.");
    }

    if (!restored_neighbor_cache) {
        // Legacy checkpoint fallback: rebuild from current positions.
        neighbor_list.set_species_positions(actin.center_x, actin.center_y, actin.center_z,
                                            myosin.center_x, myosin.center_y, myosin.center_z);
        neighbor_list.rebuild_neighbor_list();
    }

    for (int i = 0; i < actin.n; ++i) {
        myosinIndicesPerActin.deleteAllConnections(i);
    }
    for (int m = 0; m < myosin.n; ++m) {
        actinIndicesPerMyosin.deleteAllConnections(m);
        for (auto& temp_conn : actinIndicesPerMyosin_temp) {
            temp_conn.deleteAllConnections(m);
        }
    }
    for (int i = 0; i < actin.n; ++i) {
        for (int m = 0; m < myosin.n; ++m) {
            if (am_bonds[i][m] != 1) {
                continue;
            }
            myosinIndicesPerActin.addConnection(i, m);
            actinIndicesPerMyosin.addConnection(m, i);
        }
    }
    has_myosin_bond_pairs = false;

    _rebuild_active_pair_lists_from_dense();
    return target_frame;
}
