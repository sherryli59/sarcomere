#include "sarcomere.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <gsl/gsl_randist.h>

namespace {

struct PackedAABondState {
    std::vector<int> pairs;
    std::vector<int> status;
    std::vector<int> lifetime;
};

struct PackedRecoveryState {
    std::vector<int> pairs;
    std::vector<int> until;
};

PackedAABondState pack_aa_bond_state(
    const std::vector<std::vector<int>>& bonds,
    const std::vector<std::vector<int>>& status,
    const std::vector<std::vector<int>>& lifetime) {
    PackedAABondState packed;
    const int n = static_cast<int>(bonds.size());
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (bonds[i][j] != 1) {
                continue;
            }
            packed.pairs.push_back(i);
            packed.pairs.push_back(j);
            packed.status.push_back(status[i][j]);
            packed.lifetime.push_back(lifetime[i][j]);
        }
    }
    return packed;
}

std::vector<int> pack_am_bond_pairs(const std::vector<std::vector<int>>& bonds) {
    std::vector<int> packed;
    const int n_actins = static_cast<int>(bonds.size());
    for (int i = 0; i < n_actins; ++i) {
        const int n_myosins = static_cast<int>(bonds[i].size());
        for (int j = 0; j < n_myosins; ++j) {
            if (bonds[i][j] != 1) {
                continue;
            }
            packed.push_back(i);
            packed.push_back(j);
        }
    }
    return packed;
}

PackedRecoveryState pack_recovery_state(
    const std::vector<std::vector<size_t>>& recovery_until,
    size_t current_step) {
    PackedRecoveryState packed;
    const int n = static_cast<int>(recovery_until.size());
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const size_t until = recovery_until[i][j];
            if (until <= current_step) {
                continue;
            }
            packed.pairs.push_back(i);
            packed.pairs.push_back(j);
            packed.until.push_back(
                static_cast<int>(std::min(until, static_cast<size_t>(std::numeric_limits<int>::max()))));
        }
    }
    return packed;
}

void pad_int_vector(std::vector<int>& values, size_t target_size, int pad_value) {
    if (values.size() < target_size) {
        values.resize(target_size, pad_value);
    }
}

}


// Parameterized Constructor
Sarcomere::Sarcomere(int& n_actins, int& n_myosins, vector box0, double& actin_length, double& myosin_length,
        double& myosin_radius, double& am_cutoff, double& am_optimal, double& aa_cutoff, double& aa_optimal,
        double& k_on, double& k_off,
        double& base_lifetime, double& lifetime_coeff, double& diff_coeff_ratio, double& k_aa, double& kappa_aa, double& k_am, double& kappa_am, double& k_mm, double& v_am,
        std::string& filename, gsl_rng* rng, int& seed, int& fix_myosin, double& dt, double tau_rec,
        double titin_k, double titin_rest_length, bool& directional, int max_myosin_bonds,
        double max_actin_force_param, double max_myosin_force_param,
        double max_actin_torque_param, double max_myosin_torque_param)
            : actin(n_actins, actin_length, box0, rng),
              myosin(n_myosins, myosin_length, myosin_radius, box0, rng),
              myosinIndicesPerActin(n_actins),
              actinIndicesPerMyosin(n_myosins),
              neighbor_list(0.0, box0, 0.0, is_periodic),
                actin_actin_bonds(n_actins, std::vector<int>(n_actins, 0)),
                actin_actin_status(n_actins, std::vector<int>(n_actins, 0)),
                actin_actin_status_prev(n_actins, std::vector<int>(n_actins, 0)),
                actin_actin_lifetime(n_actins, std::vector<int>(n_actins, 0)),
                actin_recovery_until(n_actins, std::vector<size_t>(n_actins, 0)),
                am_bonds(n_actins, std::vector<int>(n_myosins, 0)),
                actin_forces_temp(omp_get_max_threads(), std::vector<vec>(n_actins, {0, 0, 0})),
                myosin_forces_temp(omp_get_max_threads(), std::vector<vec>(n_myosins, {0, 0, 0})),
                myosin_velocities_temp(omp_get_max_threads(), std::vector<vec>(n_myosins, {0, 0, 0})),
                actin_torques_temp(omp_get_max_threads(), std::vector<vec>(n_actins, {0, 0, 0})),
                myosin_torques_temp(omp_get_max_threads(), std::vector<vec>(n_myosins, {0, 0, 0})),
                actin_cb_status_temp(omp_get_max_threads(), std::vector<int>(n_actins, 0)),
                myosin_f_load_temp(omp_get_max_threads(), std::vector<std::array<double, 2>>(n_myosins, {0.0, 0.0})),
                cb_breakage_events_temp(omp_get_max_threads()),
                aa_completed_lifetimes_temp(omp_get_max_threads()),
                myosin_f_load(n_myosins, {0.0, 0.0}),
                actinIndicesPerMyosin_temp(omp_get_max_threads(), utils::MoleculeConnection(n_myosins)),
                rng_engines(omp_get_max_threads(), nullptr),
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
            myomesin_cutoff = (k_mm > 0.0) ? myosin.length : 0.0;
            double neighbor_cutoff = std::max(am_cutoff, aa_cutoff);
            if (k_mm > 0.0) {
                neighbor_cutoff = std::max(neighbor_cutoff, myomesin_cutoff);
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
            n_myosins_per_actin.resize(n_actins);
            am_interaction.resize(n_actins);
            for (int i = 0; i < n_actins; i++) {
                am_interaction[i].resize(n_myosins);
            }
            actin.register_feature("f_load_cb");
            actin_f_load_cb = &actin["f_load_cb"];
            actin.register_feature("myosin_binding_ratio");
            actin.register_feature("crosslink_ratio");
            actin.register_feature("partial_binding_ratio");
            if (actin_f_load_cb) {
                actin_f_load_cb->assign(n_actins, 0.0);
            }
            // Initialize thread-local RNGs
            for (int t = 0; t < omp_get_max_threads(); ++t) {
                rng_engines[t] = gsl_rng_alloc(gsl_rng_mt19937); 
                gsl_rng_set(rng_engines[t], seed + t);           
            }
            set_periodicity(is_periodic);
            myosin_bond_matrix.assign(n_myosins, std::vector<int>(n_myosins, 0));
            has_myosin_bond_pairs = false;
        }

// Destructor
Sarcomere::~Sarcomere() {}

void Sarcomere::update_system() {
    // Advance global step counter each time the system is updated
    current_step++;
    const double bundle_strength = _current_bundle_strength();
    _update_neighbors();
    #pragma omp parallel
    {   
        _set_to_zero();  
        #pragma omp barrier  
        //Step 2: Compute actin-myosin binding
        #pragma omp for schedule(runtime)
        for (int i = 0; i < actin.n; i++) {
            _process_actin_myosin_binding(i);
        }
        #pragma omp barrier
        #pragma omp single
        { 
          _enforce_myosin_bond_limit(); 
        }
        #pragma omp barrier

        // Step 3: Concatenate actinIndicesPerMyosin connections
        #pragma omp for
        for (int i = 0; i < myosin.n; ++i) {
            for (int t = 0; t < omp_get_num_threads(); ++t) {
                auto indices = actinIndicesPerMyosin_temp[t].getConnections(i);
                for (int j = 0; j < indices.size(); j++) {
                    actinIndicesPerMyosin.addConnection(i, indices[j]);
                }
            }
        }

        #pragma omp barrier  
        if (base_lifetime > 0 || lifetime_coeff > 0) {
        // Step 4: Compute catch bonds
        #pragma omp for schedule(runtime)
        for (int i = 0; i < actin.n; i++) {
            _process_catch_bonds(i);
        }
        }
        #pragma omp barrier  

        // Step 5: Reduce actin catch-bond status using max over threads
        #pragma omp for
        for (int i = 0; i < actin.n; ++i) {
            for (int t = 0; t < omp_get_num_threads(); ++t) {
                actin.cb_status[i] = std::max(actin.cb_status[i], actin_cb_status_temp[t][i]);
            }
        }

        #pragma omp barrier

        #pragma omp single
        {
            _update_myosin_bond_matrix();
        }

        #pragma omp barrier

        // Step 6: Compute actin-myosin forces
        #pragma omp for schedule(runtime)
        for (int i = 0; i < actin.n; i++) {
            _calc_am_force_velocity(i);
        }

        if (titin_k > 0 && titin_rest_length > 0) {
            #pragma omp barrier
            #pragma omp for schedule(runtime)
            for (int i = 0; i < actin.n; ++i) {
                _apply_titin_forces(i);
            }
    }

    _volume_exclusion();
    // double k_theta = 1.0;
    // _apply_cb_alignment_bias(k_theta);

        #pragma omp barrier  

        // Step 8: Reduce actin forces and angular forces
        reduce_array(actin_forces_temp, actin.force);
        reduce_array(actin_torques_temp, actin.torque);

        // Step 9: Reduce myosin forces, velocities, and angular forces
        reduce_array(myosin_forces_temp, myosin.force);
        reduce_array(myosin_velocities_temp, myosin.velocity);
        reduce_array(myosin_torques_temp, myosin.torque);

        #pragma omp barrier
        if (bundle_strength > 0.0) {
            #pragma omp single
            {
                _apply_transverse_bundling(bundle_strength);
            }
        }
        #pragma omp barrier

        // Cap forces based on specified maxima
        #pragma omp for
        for (int i = 0; i < actin.n; ++i) {
            double cap = (actin.cb_status[i] < 2) ? max_actin_force : max_myosin_force;
            if (std::isfinite(cap) && cap > 0.0) {
                double mag = actin.force[i].norm();
                if (mag > cap && mag > 0.0) {
                    actin.force[i] *= (cap / mag);
                }
            }
        }

        #pragma omp for
        for (int i = 0; i < myosin.n; ++i) {
            double cap = max_myosin_force;
            if (std::isfinite(cap) && cap > 0.0) {
                double mag = myosin.force[i].norm();
                if (mag > cap && mag > 0.0) {
                    myosin.force[i] *= (cap / mag);
                }
            }
        }

        #pragma omp for
        for (int i = 0; i < actin.n; ++i) {
            double cap = max_actin_torque;
            if (std::isfinite(cap) && cap > 0.0) {
                double mag = actin.torque[i].norm();
                if (mag > cap && mag > 0.0) {
                    actin.torque[i] *= (cap / mag);
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
        }

        #pragma omp for
        for (int i = 0; i < myosin.n; ++i) {
            double dot = myosin.velocity[i].dot(myosin.direction[i]);
            size_t dir_idx;
            if (std::abs(dot) < EPS) {
                dir_idx = (myosin_f_load[i][0] >= myosin_f_load[i][1]) ? 0 : 1;
            } else {
                dir_idx = (dot >= 0.0) ? 0 : 1;
            }
            double load_factor = 1.0 - myosin_f_load[i][dir_idx];
            if (load_factor < 0.0) {
                load_factor = 0.0;
            }
            myosin.velocity[i].normalize();
            myosin.velocity[i] *= v_am * load_factor;
        }
    }

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
}


void Sarcomere::update_system_sterics_only() {
    // Advance step counter as well for sterics-only updates
    current_step++;
    const double bundle_strength = _current_bundle_strength();
    _update_neighbors();
    #pragma omp parallel
    {   
        _set_to_zero();  
        #pragma omp barrier  
        _myosin_exclusion();
        #pragma omp barrier  
        // Step 7: Reduce actin forces and angular forces
        reduce_array(actin_forces_temp, actin.force);
        reduce_array(actin_torques_temp, actin.torque);
        // Step 8: Reduce myosin forces, velocities, and angular forces
        reduce_array(myosin_forces_temp, myosin.force);
        reduce_array(myosin_torques_temp, myosin.torque);

        #pragma omp barrier
        if (bundle_strength > 0.0) {
            #pragma omp single
            {
                _apply_transverse_bundling(bundle_strength);
            }
        }
        #pragma omp barrier

        #pragma omp for
        for (int i = 0; i < actin.n; ++i) {
            double cap = max_actin_torque;
            if (std::isfinite(cap) && cap > 0.0) {
                double mag = actin.torque[i].norm();
                if (mag > cap && mag > 0.0) {
                    actin.torque[i] *= (cap / mag);
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
                }
            }
        }
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
    if (neighbor_list.needs_rebuild()) {
        printf("Rebuilding neighbor list\n");
        neighbor_list.rebuild_neighbor_list();
    }
}


void Sarcomere::_set_to_zero() {
    // Reset all temporary forces
    #pragma omp for
    for (int t = 0; t < omp_get_max_threads(); ++t) {
        cb_breakage_events_temp[t].clear();
        aa_completed_lifetimes_temp[t].clear();
    }

    #pragma omp single
    {
        if (kmc_break_flag.size() != static_cast<size_t>(actin.n)) {
            kmc_break_flag.assign(actin.n, std::vector<int>(actin.n, 0));
        }
    }

    #pragma omp for
    for (int t = 0; t < omp_get_max_threads(); ++t) {
        for (int i = 0; i < actin.n; ++i) {
            actin_forces_temp[t][i] = {0, 0, 0};
            actin_torques_temp[t][i] = {0, 0, 0};
            actin_cb_status_temp[t][i] = 0;
        }
    }

    #pragma omp for
    for (int t = 0; t < omp_get_max_threads(); ++t) {
        for (int i = 0; i < myosin.n; ++i) {
            myosin_forces_temp[t][i] = {0, 0, 0};
            myosin_velocities_temp[t][i] = {0, 0, 0};
            myosin_torques_temp[t][i] = {0, 0, 0};
            myosin_f_load_temp[t][i] = {0.0, 0.0};
            actinIndicesPerMyosin_temp[t].deleteAllConnections(i);
        }
    }

    #pragma omp for
    for (int i = 0; i < actin.n; i++) {
        actin.update_endpoints(i);
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
        actin["myosin_binding_ratio"][i] = 0;
        actin["crosslink_ratio"][i] = 1;
        actin["partial_binding_ratio"][i] = 0;
        actin_f_load_computed[i] = false;
        myosinIndicesPerActin.deleteAllConnections(i);
        // Ensure and reset per-step KMC break flags
        std::fill(kmc_break_flag[i].begin(), kmc_break_flag[i].end(), 0);
        for (int j = 0; j < actin.n; j++){
            actin_actin_bonds_prev[i][j] = actin_actin_bonds[i][j];
            actin_actin_status_prev[i][j] = actin_actin_status[i][j];
            actin_actin_bonds[i][j] = 0;
            actin_actin_status[i][j] = 0;
            actin_actin_lifetime_prev[i][j] = actin_actin_lifetime[i][j];
            actin_actin_lifetime[i][j] = 0;
        }
        for (int j = 0; j < myosin.n; j++){
            am_bonds_prev[i][j] = am_bonds[i][j];
            am_bonds[i][j] = 0;
        }
    }

    #pragma omp for
    for (int i = 0; i < myosin.n; i++) {
        myosin.update_endpoints(i);
        myosin.force[i] = {0, 0, 0};
        myosin.torque[i] = {0, 0, 0};
        myosin.velocity[i] = {0, 0, 0};
        actinIndicesPerMyosin.deleteAllConnections(i);
        myosin_f_load[i] = {0.0, 0.0};
    }
    #pragma omp for
    for (size_t i = 0; i < actin.center.size(); i++) {
        // Get actin and myosin neighbors for actin particle `i`
        actin_neighbors_by_species[i] = neighbor_list.get_neighbors_by_type(i);
    }
}

void Sarcomere::_process_actin_myosin_binding(int& i) {
    int thread_id = omp_get_thread_num();
    auto& local_actinIndicesPerMyosin = actinIndicesPerMyosin_temp[thread_id];
    auto myosin_neighbors = actin_neighbors_by_species[i].second;
    double f_load_contrib = 0;
    double f_load_contrib_cb = 0;
    for (int index = 0; index < myosin_neighbors.size(); index++) {
        int j = myosin_neighbors[index];
        am_interaction[i][j] = geometry::analyze_am(
            actin.left_end[i], actin.right_end[i], myosin.left_end[j], myosin.right_end[j],
            am_cutoff, box, is_periodic, directional);
        if (am_interaction[i][j].myosin_binding_ratio > EPS) {
            if (am_interaction[i][j].partial_binding_ratio > EPS || !directional) {
                double partial_ratio = am_interaction[i][j].partial_binding_ratio;
                double binding_ratio = am_interaction[i][j].myosin_binding_ratio;
                    if (actin_crosslink_ratio[i] > am_interaction[i][j].crosslinkable_ratio) {
                        actin_crosslink_ratio[i] = am_interaction[i][j].crosslinkable_ratio;
                        actin_crosslink_start[i] = am_interaction[i][j].crosslinkable_start;
                        actin_crosslink_end[i] = am_interaction[i][j].crosslinkable_end;
                        if (actin_crosslink_ratio[i] < EPS) {
                            auto actin_neighbors = actin_neighbors_by_species[i].first;
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
                                    10*k_aa);
                                double force_magnitude = repulsive_force.norm();
                                printf("Warning: Actin %d has zero crosslink ratio due to myosin %d,binding ratio is %f, partial binding ratio is %f, angle is %f, repulsion is %f\n", i, j, am_interaction[i][j].myosin_binding_ratio, 
                                    am_interaction[i][j].partial_binding_ratio, angle, force_magnitude);
                            }
                        }
                    }
                if (!directional || partial_ratio > EPS) {
                    // Record all actin–myosin attachments regardless of bond state
                    local_actinIndicesPerMyosin.addConnection(j, i);
                    myosinIndicesPerActin.addConnection(i, j);
                    n_myosins_per_actin[i]++;
                    double abs_cos = std::abs(actin.direction[i].dot(myosin.direction[j]));
                    if (abs_cos > actin_basic_tension[i]) {
                        actin_basic_tension[i] = abs_cos;
                    }
                    if (actin["partial_binding_ratio"][i] < partial_ratio) {
                        actin["partial_binding_ratio"][i] = partial_ratio;
                    }
                    double binding_ratio = am_interaction[i][j].myosin_binding_ratio;
                    if (actin["myosin_binding_ratio"][i] < binding_ratio) {
                        actin["myosin_binding_ratio"][i] = binding_ratio;
                    }

                    // Immediately register a bond when geometrically allowed
                    am_bonds[i][j] = 1;
                    // estimate load on actin due to this myosin
                    double partial = directional ? am_interaction[i][j].partial_binding_ratio
                                                 : am_interaction[i][j].myosin_binding_ratio;
                    double abs_cos_angle = std::abs(actin.direction[i].dot(myosin.direction[j]));
                    double contrib = 3.0 * std::min(partial, 1.0/3.0);
                    double contrib_cb = contrib * abs_cos_angle;
                    if (contrib > f_load_contrib){
                        f_load_contrib = contrib;
                    }
                    if (contrib_cb > f_load_contrib_cb) {
                        f_load_contrib_cb = contrib_cb;
                    }
                }
            }
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
    actin["crosslink_ratio"][i] = actin_crosslink_ratio[i];

}

void Sarcomere::_process_catch_bonds(int& i) {
    std::vector <int> actin_neighbors = actin_neighbors_by_species[i].first;
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


void Sarcomere::_enforce_myosin_bond_limit() {
    for (int j = 0; j < myosin.n; ++j) {
        struct BoundInfo {
            int actin_idx;
            bool cb;
            bool prev;
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
                bound.push_back({i, cb, prev, dot, true});
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
            myosinIndicesPerActin.deleteConnection(info.actin_idx, j);
            for (auto& temp_conn : actinIndicesPerMyosin_temp) {
                temp_conn.deleteConnection(j, info.actin_idx);
            }
            info.active = false;
        };

        auto cmp_priority = [](const BoundInfo* lhs, const BoundInfo* rhs) {
            auto score = [](const BoundInfo* entry) {
                int value = 0;
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

        if (per_group_limit > 0) {
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
    std::vector<int> myosin_indices = myosinIndicesPerActin.getConnections(i);
    const std::vector<int>& myosin_neighbors = actin_neighbors_by_species[i].second;
    vec velocity = v_am * actin.direction[i];
    for (int index = 0; index < myosin_indices.size(); index++) {
        int j = myosin_indices[index];
        if (am_bonds[i][j] != 1) {
            continue;
        }
        //scale by angle between actin and myosin
        double abs_cos_angle = std::abs(actin.direction[i].dot(myosin.direction[j]));
        double active_ratio = directional ? am_interaction[i][j].partial_binding_ratio
                                          : am_interaction[i][j].myosin_binding_ratio;
        double normalized_partial_ratio = 3.0 * std::min(active_ratio, 1.0/3.0);
        //apply exponential scaling to partial binding ratio
        normalized_partial_ratio = (1 - std::exp(-2 * normalized_partial_ratio)) / (1 - std::exp(-2));
        vector force_vec = compute_am_force_and_energy(
            actin, myosin, i, j, box, k_am * normalized_partial_ratio, kappa_am, am_cutoff, am_optimal);
        local_actin_forces[i].x += force_vec[0];
        local_actin_forces[i].y += force_vec[1];
        local_actin_forces[i].z += force_vec[2];
        local_myosin_forces[j].x -= force_vec[0];
        local_myosin_forces[j].y -= force_vec[1];
        local_myosin_forces[j].z -= force_vec[2];
        local_actin_torques[i].x += force_vec[3];
        local_actin_torques[i].y += force_vec[4];
        local_actin_torques[i].z += force_vec[5];
        local_myosin_torques[j].x += force_vec[6];
        local_myosin_torques[j].y += force_vec[7];
        local_myosin_torques[j].z += force_vec[8];

        if (actin.cb_status[i] == 2) {
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

    if (actin.cb_status[i] == 2) {
        if (myosin_indices.size() == 0) {
            printf("Warning: Actin %d has catch bond status but no myosin bound\n", i);
        }
        actin.velocity[i] = (1 - actin.f_load[i]) * velocity;
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
                    20*k_aa,
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


void Sarcomere::_apply_titin_forces(int& i) {
    int thread_id = omp_get_thread_num();
    auto& local_actin_forces = actin_forces_temp[thread_id];
    auto& local_myosin_forces = myosin_forces_temp[thread_id];

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
                printf("actin crosslink ratio %f, stretch %f, actin-myosin dot product %f, dot product %f, actin center (%f, %f, %f), actin direction (%f, %f, %f), myosin center (%f, %f, %f) \n",
                     actin_crosslink_ratio[act_b], stretch, am_dot_prod, dot_prod, 
                     actin.center[act_a].x, actin.center[act_a].y, actin.center[act_a].z,
                     actin.direction[act_a].x, actin.direction[act_a].y, actin.direction[act_a].z,
                    myo_center.x, myo_center.y, myo_center.z);
            }
        };

        apply_spring(i, j);
        apply_spring(j, i);
    }
}


void Sarcomere::_apply_myomesin_spring(int i, int j, std::vector<vec>& local_myosin_forces) {
    if (k_mm <= 0.0 || myomesin_cutoff <= 0.0) {
        return;
    }
    vec diff = utils::pbc_diff_masked(myosin.center[i], myosin.center[j], box, is_periodic);

        // get normalized directions
    vec di = myosin.direction[i];
    vec dj = myosin.direction[j];
    di.normalize();
    dj.normalize();
    // enforce apolar consistency: flip dj if it's anti-parallel to di
    if (di.dot(dj) < 0.0) {
        dj = -dj;
    }
    // average axis direction
    vec e_axis = di + dj;
    e_axis.normalize();
    // axial separation
    double axial = diff.dot(e_axis);
    double dist_axial = std::abs(axial);
    if (dist_axial <= 1e-12 || dist_axial > myomesin_cutoff) {
        return;
    }
    // unit force direction along ±e_axis
    vec unit_axial = (axial >= 0.0) ? e_axis : -e_axis;
    double stretch = dist_axial/myomesin_cutoff; //normalized stretch
    double force_scalar = -k_mm * stretch;
    vec force_on_i = force_scalar * unit_axial;
    vec force_on_j = -force_on_i;
    //print myosin centers, directions and force vectors
    // printf("Myomesin spring between myosin %d (center (%f, %f, %f), direction (%f, %f, %f)) and %d (center (%f, %f, %f), direction (%f, %f, %f)): stretch %f, force on i (%f, %f, %f)\n",
    //     i, myosin.center[i].x, myosin.center[i].y, myosin.center[i].z,
    //        myosin.direction[i].x, myosin.direction[i].y, myosin.direction[i].z,
    //     j, myosin.center[j].x, myosin.center[j].y, myosin.center[j].z,
    //        myosin.direction[j].x, myosin.direction[j].y, myosin.direction[j].z,
    //      stretch,
    //     force_on_i.x, force_on_i.y, force_on_i.z);
}



void Sarcomere::_volume_exclusion(){
    const double EPS_FORCE = 1e-9;
    double myosin_cutoff = 2.0 * myosin.radius;
    const double myomesin_distance_limit = 2.0 * am_cutoff;
    #pragma omp for schedule(runtime)
    for (int i = 0; i<myosin.n; i++){
        auto result = neighbor_list.get_neighbors_by_type(i+actin.n);
        const std::vector<int>& myosin_indices = result.second;
        auto& local_myosin_forces = myosin_forces_temp[omp_get_thread_num()];
        for (int index = 0; index < static_cast<int>(myosin_indices.size()); index++){
            int j = myosin_indices[index];
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
                    100 * k_aa,
                    local_myosin_forces[i],
                    local_myosin_forces[j],
                    seg_distance);
                if (seg_distance <= myomesin_distance_limit) {
                    _apply_myomesin_spring(i, j, local_myosin_forces);
                }
            }
        }
    }

    #pragma omp for schedule(runtime)
    for (int i = 0; i < actin.n; i++){
        auto result = neighbor_list.get_neighbors_by_type(i);
        const std::vector<int>& actin_indices = result.first;
        auto& local_actin_forces = actin_forces_temp[omp_get_thread_num()];
        if (actin.cb_status[i]==0) {
            continue;
        }
        for (int index = 0; index < static_cast<int>(actin_indices.size()); index++){
            int j = actin_indices[index];
            if (actin.cb_status[j]==0) {
                continue;
            }
            if (i<j){
                apply_actin_repulsion(
                    actin,
                    i,
                    j,
                    box,
                    aa_optimal,
                    10*k_aa,
                    local_actin_forces[i],
                    local_actin_forces[j]);
            }
        }
    }
}

void Sarcomere::_myosin_exclusion(){
    const double EPS_FORCE = 1e-9;
    double myosin_cutoff = 2.0 * myosin.radius;
    #pragma omp for schedule(runtime)
    for (int i = 0; i<myosin.n; i++){
        auto result = neighbor_list.get_neighbors_by_type(i+actin.n);
        const std::vector<int>& myosin_indices = result.second;
        auto& local_myosin_forces = myosin_forces_temp[omp_get_thread_num()];
        for (int index = 0; index < static_cast<int>(myosin_indices.size()); index++){
            int j = myosin_indices[index];
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
                    100*k_aa,
                    local_myosin_forces[i],
                    local_myosin_forces[j],
                    seg_distance);
            }
        }
    }
}

int Sarcomere::determine_cb_status(int& i, int& j){
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
    if (distance < aa_cutoff) {
        crosslink = true;
        printf("Actins %d and %d are crosslinked, distance: %f, cos_angle: %f, crosslink ratio: %f, %f\n",
               i, j, distance, cos_angle,crosslink_i, crosslink_j);
    }
    //}
    if (!crosslink){
        if (was_strong){
            int myosin_i = myosinIndicesPerActin.getConnections(i)[0];
            int myosin_j = myosinIndicesPerActin.getConnections(j)[0];
            printf("Actins %d and %d no longer crosslinked, distance: %f, cos_angle: %f, crosslink ratio: %f, %f myosin indices %d, %d\n",
                 i, j, distance, cos_angle,crosslink_i, crosslink_j, myosin_i, myosin_j);
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
    bool catch_bond =(actin_basic_tension[i]>EPS && actin_basic_tension[j]>EPS);
    if (directional){
        catch_bond = (catch_bond && cos_angle<0);
    }
    if (!catch_bond){
        if (was_strong){
            printf("Actins %d and %d no longer form catch bonds, distance: %f, cos_angle: %f, tensions: %f, %f\n",
                   i, j, distance, cos_angle, actin_basic_tension[i], actin_basic_tension[j]);
            record_break();
        }
        return 1;
    }
    auto& myosin_indices_i = myosinIndicesPerActin.getConnections(i);
    auto& myosin_indices_j = myosinIndicesPerActin.getConnections(j);
    if (myosin_indices_i.empty() || myosin_indices_j.empty()){
        if (was_strong){
            printf("Actins %d and %d no longer form catch bonds (no myosin), distance: %f, cos_angle: %f, tensions: %f, %f\n",
                   i, j, distance, cos_angle, actin_basic_tension[i], actin_basic_tension[j]);
            record_break();
        }
        return 1;
    }
    for (int mi : myosin_indices_i){
        if (am_bonds[i][mi] != 1) continue;
        for (int mj : myosin_indices_j){
            if (am_bonds[j][mj] != 1) continue;
            if (mi != mj){
                return 2;
            }
        }
    }
    if (was_strong){
        printf("Actins %d and %d no longer form catch bonds (no shared myosin), distance: %f, cos_angle: %f, tensions: %f, %f\n",
               i, j, distance, cos_angle, actin_basic_tension[i], actin_basic_tension[j]);
        record_break();
    }
    return 1;
}


bool Sarcomere::_cb_decide(int& i, int& j, int status){
    if (status == 0){
        return false;
    }
    bool was_bonded_prev = (actin_actin_bonds_prev[i][j] == 1 || actin_recovery_until[i][j] > current_step);
    if (!was_bonded_prev && bond_recovery_steps > 0 && actin_recovery_until[i][j] > current_step) {
        return false;
    }
    int thread_id = omp_get_thread_num();
    auto& local_actin_cb_status = actin_cb_status_temp[thread_id];
    auto& local_completed_lifetimes = aa_completed_lifetimes_temp[thread_id];

    double rand = gsl_rng_uniform(rng_engines[thread_id]);
    double abs_cos_angle = std::abs(actin.direction[i].dot(actin.direction[j]));
    double f_load_i = actin_f_load_cb ? (*actin_f_load_cb)[i] : actin.f_load[i];
    double f_load_j = actin_f_load_cb ? (*actin_f_load_cb)[j] : actin.f_load[j];
    double f_load = abs_cos_angle * std::min(f_load_i, f_load_j);
    if (actin_actin_bonds_prev[i][j] == 1) {
        double k_off_adjusted = dt /(base_lifetime + lifetime_coeff * f_load);
        if (rand < k_off_adjusted){
            if (f_load>0) {
            printf("k_off_adjusted: %f, rand: %f, f_load: %f, abs_cos_angle: %f, lifetime: %f\n",
                   k_off_adjusted, rand, f_load, abs_cos_angle,
                   base_lifetime + lifetime_coeff * f_load);
            printf("actual lifetime: %f\n", (current_step - aa_attach_step[i][j]) * dt);
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
            // Mark this pair as broken by KMC in this step
            kmc_break_flag[i][j] = 1;
            kmc_break_flag[j][i] = 1;
            if (bond_recovery_steps > 0) {
                size_t release = current_step + bond_recovery_steps;
                actin_recovery_until[i][j] = release;
                actin_recovery_until[j][i] = release;
            }
            return false;
        }
    }

    else{
        if (rand >= k_on * dt){
            //if (status==2){
            // printf("Actins %d and %d fail to form catch bond with status %d, rand %f, k_on*dt %f, f_load %f, abs_cos_angle %f\n",
            //     i, j, status, rand, k_on*dt, f_load, abs_cos_angle);
            //}
            return false;
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
 
    vector force_vec = compute_aa_force_and_energy(actin, i, j, box, k_aa, kappa_aa, aa_cutoff, aa_optimal);
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
    // If the bond is newly formed this step, stamp attach time
    if (actin_actin_bonds_prev[i][j] == 0) {
        if (aa_attach_step.size()) {
            aa_attach_step[i][j] = static_cast<int>(current_step);
            aa_attach_step[j][i] = static_cast<int>(current_step);
        }
        if (bond_recovery_steps > 0) {
            actin_recovery_until[i][j] = current_step;
            actin_recovery_until[j][i] = current_step;
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
    const std::vector<std::vector<int>>& actin_actin_bonds,
    const std::vector<std::vector<int>>& actin_actin_status,
    const utils::MoleculeConnection& myosinIndicesPerActin)
{
    // First, flatten the actin-actin bonds matrix (upper-triangle only)
    // into a vector of doubled indices.
    std::vector<double> flatActinBonds;
    for (int i = 0; i < actin.n; ++i) {
        for (int j = i + 1; j < actin.n; ++j) {
            if (actin_actin_bonds[i][j] == 1 &&
                actin_actin_status[i][j] == 2 &&
                actin.cb_status[i] == 2 &&
                actin.cb_status[j] == 2) {
                flatActinBonds.push_back(static_cast<double>(i));
                flatActinBonds.push_back(static_cast<double>(j));
            }
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

    // Extract actin–myosin bonds directly from the bond matrix.
    std::vector<double> flatActinMyosinBonds;
    for (int a = 0; a < actin.n; ++a) {
        for (int m = 0; m < myosin.n; ++m) {
            if (am_bonds[a][m] == 1) {
                flatActinMyosinBonds.push_back(static_cast<double>(a));
                flatActinMyosinBonds.push_back(static_cast<double>(m));
            }
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
        myosin_bond_matrix.assign(myosin.n, std::vector<int>(myosin.n, 0));
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

        ensure_state_dataset_int("current_step", 1);
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
        ensure_state_pair_dataset("am_pairs_prev", max_am_pairs, 2);
        ensure_state_pair_dataset("am_pairs_current", max_am_pairs, 2);
        ensure_state_dataset_int("aa_recovery_count", 1);
        ensure_state_pair_dataset("aa_recovery_pairs", max_recovery_pairs, 2);
        ensure_state_pair_dataset("aa_recovery_until_values", max_recovery_pairs, 1);
        ensure_state_dataset_double("neighbor_last_actin_x", actin_width);
        ensure_state_dataset_double("neighbor_last_actin_y", actin_width);
        ensure_state_dataset_double("neighbor_last_actin_z", actin_width);
        ensure_state_dataset_double("neighbor_last_myosin_x", myosin_width);
        ensure_state_dataset_double("neighbor_last_myosin_y", myosin_width);
        ensure_state_dataset_double("neighbor_last_myosin_z", myosin_width);

        // Save current_step (scalar value)
        std::vector<int> step_vec = {static_cast<int>(current_step)};
        append_to_dataset_int(group_state, "current_step", step_vec, {1, 1});

        auto aa_prev = pack_aa_bond_state(
            actin_actin_bonds_prev, actin_actin_status_prev, actin_actin_lifetime_prev);
        auto aa_current = pack_aa_bond_state(
            actin_actin_bonds, actin_actin_status, actin_actin_lifetime);
        auto am_prev = pack_am_bond_pairs(am_bonds_prev);
        auto am_current = pack_am_bond_pairs(am_bonds);
        auto aa_recovery = pack_recovery_state(actin_recovery_until, current_step);

        const hsize_t aa_prev_count = static_cast<hsize_t>(aa_prev.status.size());
        const hsize_t aa_current_count = static_cast<hsize_t>(aa_current.status.size());
        const hsize_t am_prev_count = static_cast<hsize_t>(am_prev.size() / 2);
        const hsize_t am_current_count = static_cast<hsize_t>(am_current.size() / 2);
        const hsize_t aa_recovery_count = static_cast<hsize_t>(aa_recovery.until.size());
        if (aa_prev_count > max_aa_pairs || aa_current_count > max_aa_pairs ||
            am_prev_count > max_am_pairs || am_current_count > max_am_pairs ||
            aa_recovery_count > max_recovery_pairs) {
            throw std::runtime_error("Bond pair list exceeds configured compact state capacity.");
        }

        std::vector<int> aa_prev_pairs = aa_prev.pairs;
        std::vector<int> aa_prev_status = aa_prev.status;
        std::vector<int> aa_prev_lifetime = aa_prev.lifetime;
        std::vector<int> aa_current_pairs = aa_current.pairs;
        std::vector<int> aa_current_status = aa_current.status;
        std::vector<int> aa_current_lifetime = aa_current.lifetime;
        std::vector<int> aa_current_attach_step;
        aa_current_attach_step.reserve(aa_current.status.size());
        for (size_t idx = 0; idx + 1 < aa_current.pairs.size(); idx += 2) {
            const int a = aa_current.pairs[idx];
            const int b = aa_current.pairs[idx + 1];
            aa_current_attach_step.push_back(aa_attach_step[a][b]);
        }
        std::vector<int> am_prev_pairs = am_prev;
        std::vector<int> am_current_pairs = am_current;
        std::vector<int> aa_recovery_pairs = aa_recovery.pairs;
        std::vector<int> aa_recovery_until_values = aa_recovery.until;

        pad_int_vector(aa_prev_pairs, static_cast<size_t>(max_aa_pairs * 2), -1);
        pad_int_vector(aa_prev_status, static_cast<size_t>(max_aa_pairs), 0);
        pad_int_vector(aa_prev_lifetime, static_cast<size_t>(max_aa_pairs), 0);
        pad_int_vector(aa_current_pairs, static_cast<size_t>(max_aa_pairs * 2), -1);
        pad_int_vector(aa_current_status, static_cast<size_t>(max_aa_pairs), 0);
        pad_int_vector(aa_current_lifetime, static_cast<size_t>(max_aa_pairs), 0);
        pad_int_vector(aa_current_attach_step, static_cast<size_t>(max_aa_pairs), -1);
        pad_int_vector(am_prev_pairs, static_cast<size_t>(max_am_pairs * 2), -1);
        pad_int_vector(am_current_pairs, static_cast<size_t>(max_am_pairs * 2), -1);
        pad_int_vector(aa_recovery_pairs, static_cast<size_t>(max_recovery_pairs * 2), -1);
        pad_int_vector(aa_recovery_until_values, static_cast<size_t>(max_recovery_pairs), 0);

        append_to_dataset_int(group_state, "aa_pairs_prev_count", {static_cast<int>(aa_prev_count)}, {1, 1});
        append_to_dataset_int(group_state, "aa_pairs_current_count", {static_cast<int>(aa_current_count)}, {1, 1});
        append_to_dataset_int(group_state, "am_pairs_prev_count", {static_cast<int>(am_prev_count)}, {1, 1});
        append_to_dataset_int(group_state, "am_pairs_current_count", {static_cast<int>(am_current_count)}, {1, 1});
        append_to_dataset_int(group_state, "aa_recovery_count", {static_cast<int>(aa_recovery_count)}, {1, 1});
        append_to_dataset_int(group_state, "aa_pairs_prev", aa_prev_pairs, {1, max_aa_pairs, 2});
        append_to_dataset_int(group_state, "aa_status_prev", aa_prev_status, {1, max_aa_pairs, 1});
        append_to_dataset_int(group_state, "aa_lifetime_prev", aa_prev_lifetime, {1, max_aa_pairs, 1});
        append_to_dataset_int(group_state, "aa_pairs_current", aa_current_pairs, {1, max_aa_pairs, 2});
        append_to_dataset_int(group_state, "aa_status_current", aa_current_status, {1, max_aa_pairs, 1});
        append_to_dataset_int(group_state, "aa_lifetime_current", aa_current_lifetime, {1, max_aa_pairs, 1});
        append_to_dataset_int(group_state, "aa_attach_step_current", aa_current_attach_step, {1, max_aa_pairs, 1});
        append_to_dataset_int(group_state, "am_pairs_prev", am_prev_pairs, {1, max_am_pairs, 2});
        append_to_dataset_int(group_state, "am_pairs_current", am_current_pairs, {1, max_am_pairs, 2});
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
}

void Sarcomere::save_resume_snapshot() {
    H5::H5File file(filename, H5F_ACC_RDWR);
    H5::Group group_resume;
    if (file.nameExists("/resume")) {
        group_resume = file.openGroup("/resume");
    } else {
        group_resume = file.createGroup("/resume");
    }

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
        const std::string full_path = "/resume/" + name;
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
        const std::string full_path = "/resume/" + name;
        if (file.nameExists(full_path)) {
            H5Ldelete(file.getId(), full_path.c_str(), H5P_DEFAULT);
        }
        H5::DataSpace dataspace(dims.size(), dims.data());
        H5::FloatType datatype(H5::PredType::IEEE_F64LE);
        H5::DataSet dataset = group_resume.createDataSet(name, datatype, dataspace);
        dataset.write(data.data(), H5::PredType::IEEE_F64LE);
    };

    auto aa_prev = pack_aa_bond_state(
        actin_actin_bonds_prev, actin_actin_status_prev, actin_actin_lifetime_prev);
    auto aa_current = pack_aa_bond_state(
        actin_actin_bonds, actin_actin_status, actin_actin_lifetime);
    auto am_prev = pack_am_bond_pairs(am_bonds_prev);
    auto am_current = pack_am_bond_pairs(am_bonds);
    auto aa_recovery = pack_recovery_state(actin_recovery_until, current_step);

    if (aa_recovery.until.size() > static_cast<size_t>(max_recovery_pairs)) {
        throw std::runtime_error("Recovery pair list exceeds configured compact resume capacity.");
    }

    std::vector<int> aa_prev_pairs = aa_prev.pairs;
    std::vector<int> aa_prev_status = aa_prev.status;
    std::vector<int> aa_prev_lifetime = aa_prev.lifetime;
    std::vector<int> aa_current_pairs = aa_current.pairs;
    std::vector<int> aa_current_status = aa_current.status;
    std::vector<int> aa_current_lifetime = aa_current.lifetime;
    std::vector<int> aa_current_attach_step;
    aa_current_attach_step.reserve(aa_current.status.size());
    for (size_t idx = 0; idx + 1 < aa_current.pairs.size(); idx += 2) {
        const int a = aa_current.pairs[idx];
        const int b = aa_current.pairs[idx + 1];
        aa_current_attach_step.push_back(aa_attach_step[a][b]);
    }
    std::vector<int> am_prev_pairs = am_prev;
    std::vector<int> am_current_pairs = am_current;
    std::vector<int> aa_recovery_pairs = aa_recovery.pairs;
    std::vector<int> aa_recovery_until_values = aa_recovery.until;

    pad_int_vector(aa_prev_pairs, static_cast<size_t>(max_aa_pairs * 2), -1);
    pad_int_vector(aa_prev_status, static_cast<size_t>(max_aa_pairs), 0);
    pad_int_vector(aa_prev_lifetime, static_cast<size_t>(max_aa_pairs), 0);
    pad_int_vector(aa_current_pairs, static_cast<size_t>(max_aa_pairs * 2), -1);
    pad_int_vector(aa_current_status, static_cast<size_t>(max_aa_pairs), 0);
    pad_int_vector(aa_current_lifetime, static_cast<size_t>(max_aa_pairs), 0);
    pad_int_vector(aa_current_attach_step, static_cast<size_t>(max_aa_pairs), -1);
    pad_int_vector(am_prev_pairs, static_cast<size_t>(max_am_pairs * 2), -1);
    pad_int_vector(am_current_pairs, static_cast<size_t>(max_am_pairs * 2), -1);
    pad_int_vector(aa_recovery_pairs, static_cast<size_t>(max_recovery_pairs * 2), -1);
    pad_int_vector(aa_recovery_until_values, static_cast<size_t>(max_recovery_pairs), 0);

    replace_dataset_int("current_step", {1}, {static_cast<int>(current_step)});
    replace_dataset_double("actin_center", {n_actins, 3}, flatten_3d_array(actin.center));
    replace_dataset_double("actin_direction", {n_actins, 3}, flatten_3d_array(actin.direction));
    replace_dataset_double("myosin_center", {n_myosins, 3}, flatten_3d_array(myosin.center));
    replace_dataset_double("myosin_direction", {n_myosins, 3}, flatten_3d_array(myosin.direction));
    replace_dataset_int("aa_pairs_prev_count", {1}, {static_cast<int>(aa_prev.status.size())});
    replace_dataset_int("aa_pairs_current_count", {1}, {static_cast<int>(aa_current.status.size())});
    replace_dataset_int("am_pairs_prev_count", {1}, {static_cast<int>(am_prev.size() / 2)});
    replace_dataset_int("am_pairs_current_count", {1}, {static_cast<int>(am_current.size() / 2)});
    replace_dataset_int("aa_pairs_prev", {max_aa_pairs, 2}, aa_prev_pairs);
    replace_dataset_int("aa_status_prev", {max_aa_pairs, 1}, aa_prev_status);
    replace_dataset_int("aa_lifetime_prev", {max_aa_pairs, 1}, aa_prev_lifetime);
    replace_dataset_int("aa_pairs_current", {max_aa_pairs, 2}, aa_current_pairs);
    replace_dataset_int("aa_status_current", {max_aa_pairs, 1}, aa_current_status);
    replace_dataset_int("aa_lifetime_current", {max_aa_pairs, 1}, aa_current_lifetime);
    replace_dataset_int("aa_attach_step_current", {max_aa_pairs, 1}, aa_current_attach_step);
    replace_dataset_int("am_pairs_prev", {max_am_pairs, 2}, am_prev_pairs);
    replace_dataset_int("am_pairs_current", {max_am_pairs, 2}, am_current_pairs);
    replace_dataset_int("aa_recovery_count", {1}, {static_cast<int>(aa_recovery.until.size())});
    replace_dataset_int("aa_recovery_pairs", {max_recovery_pairs, 2}, aa_recovery_pairs);
    replace_dataset_int("aa_recovery_until_values", {max_recovery_pairs, 1}, aa_recovery_until_values);
    if (file.nameExists("/resume/actin_recovery_until")) {
        H5Ldelete(file.getId(), "/resume/actin_recovery_until", H5P_DEFAULT);
    }

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
}

int Sarcomere::load_state(int& n_frames, int frame_index){
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

    auto assign_aa_matrix = [&](const std::vector<int>& flat, std::vector<std::vector<int>>& matrix) {
        if (flat.size() != aa_stride) {
            return;
        }
        for (int i = 0; i < actin.n; ++i) {
            for (int j = 0; j < actin.n; ++j) {
                matrix[i][j] = flat[static_cast<size_t>(i) * actin.n + j];
            }
        }
    };

    auto assign_am_matrix = [&](const std::vector<int>& flat, std::vector<std::vector<int>>& matrix) {
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

    auto clear_aa_state = [&](std::vector<std::vector<int>>& bonds,
                              std::vector<std::vector<int>>& status,
                              std::vector<std::vector<int>>& lifetime) {
        for (int i = 0; i < actin.n; ++i) {
            std::fill(bonds[i].begin(), bonds[i].end(), 0);
            std::fill(status[i].begin(), status[i].end(), 0);
            std::fill(lifetime[i].begin(), lifetime[i].end(), 0);
        }
    };

    auto clear_am_state = [&](std::vector<std::vector<int>>& bonds) {
        for (int i = 0; i < actin.n; ++i) {
            std::fill(bonds[i].begin(), bonds[i].end(), 0);
        }
    };

    auto clear_recovery_state = [&]() {
        for (int i = 0; i < actin.n; ++i) {
            std::fill(actin_recovery_until[i].begin(), actin_recovery_until[i].end(), 0);
        }
    };

    auto restore_aa_state_from_pairs = [&](H5::Group& group,
                                           const std::string& count_name,
                                           const std::string& pair_name,
                                           const std::string& status_name,
                                           const std::string& lifetime_name,
                                           std::vector<std::vector<int>>& bonds,
                                           std::vector<std::vector<int>>& status,
                                           std::vector<std::vector<int>>& lifetime) -> bool {
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
                                           std::vector<std::vector<int>>& bonds) -> bool {
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
            const size_t until = static_cast<size_t>(std::max(0, until_values[idx]));
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
                                                 std::vector<std::vector<int>>& bonds,
                                                 std::vector<std::vector<int>>& status,
                                                 std::vector<std::vector<int>>& lifetime) -> bool {
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
                                                 std::vector<std::vector<int>>& bonds) -> bool {
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
            const size_t until = static_cast<size_t>(std::max(0, until_values[idx]));
            actin_recovery_until[a][b] = until;
            actin_recovery_until[b][a] = until;
        }
        return true;
    };

    bool loaded_state_group = false;
    bool restored_neighbor_cache = false;
    bool loaded_resume_snapshot = false;
    try {
        H5::H5File file(filename, H5F_ACC_RDONLY);
        H5::Group group_state(file.openGroup("/state"));
        loaded_state_group = true;

        std::vector<int> flat;
        std::vector<hsize_t> dims;

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
                    actin_recovery_until[i][j] =
                        static_cast<size_t>(flat[static_cast<size_t>(i) * actin.n + j]);
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

        // Restore per-actin cb_status from the matrix so diagnostic output remains consistent.
        for (int i = 0; i < actin.n; ++i) {
            int max_status = 0;
            for (int j = 0; j < actin.n; ++j) {
                max_status = std::max(max_status, actin_actin_status[i][j]);
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
            if (file.nameExists("/resume")) {
                try {
                    H5::Group group_resume(file.openGroup("/resume"));
                    std::vector<int> ints;
                    std::vector<double> doubles;

                    if (load_fixed_vector_int(group_resume, "current_step", 1, ints)) {
                        const size_t resume_step = static_cast<size_t>(std::max(0, ints[0]));
                        if (resume_step < current_step) {
                            throw H5::Exception("Sarcomere::load_state", "Stale resume snapshot");
                        }
                        current_step = resume_step;
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

                    restore_am_state_from_fixed_pairs(
                        group_resume, "am_pairs_prev_count", "am_pairs_prev", am_bonds_prev);
                    restore_am_state_from_fixed_pairs(
                        group_resume, "am_pairs_current_count", "am_pairs_current", am_bonds);

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
                                actin_recovery_until[i][j] =
                                    static_cast<size_t>(ints[static_cast<size_t>(i) * actin.n + j]);
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
                    loaded_resume_snapshot = true;
                } catch (H5::Exception&) {
                    loaded_resume_snapshot = loaded_state_group;
                    printf("Warning: Could not load /resume snapshot from %s. Using /state as substitute.\n",
                           filename.c_str());
                }
            } else {
                loaded_resume_snapshot = loaded_state_group;
                printf("Warning: /resume snapshot missing in %s. Using /state as substitute.\n",
                       filename.c_str());
            }
        }

    } catch (H5::Exception&) {
        loaded_state_group = false;
    }

    if (!loaded_state_group) {
        printf("Warning: Could not load /state group. Falling back to approximate resume state.\n");
        current_step = 0;
        for (int i = 0; i < actin.n; ++i) {
            for (int j = 0; j < actin.n; ++j) {
                actin_actin_bonds_prev[i][j] = actin_actin_bonds[i][j];
                actin_actin_status_prev[i][j] = 0;
                actin_actin_lifetime_prev[i][j] = 0;
                actin_actin_status[i][j] = (actin_actin_bonds[i][j] == 1) ? 2 : 0;
                actin_actin_lifetime[i][j] = 0;
                actin_recovery_until[i][j] = 0;
            }
            for (int j = 0; j < myosin.n; ++j) {
                am_bonds_prev[i][j] = 0;
                am_bonds[i][j] = 0;
            }
        }
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

    return target_frame;
}
