#ifndef SARCOMERE_H
#define SARCOMERE_H

#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <tuple>
#include <array>
#include <omp.h>
#include <mutex>
#include <cstddef>

#include "components.h"
#include "utils.h"
#include "interaction.h"
#include "geometry.h"
#include "h5_utils.h"
#include "neighborlist.h"

// Alias for convenience.
using vector = std::vector<double>;
using interaction = geometry::am_interaction;
const double EPS = 1e-6;

class Sarcomere {
public:
    // Public Data Members
    Filament actin;
    Myosin myosin;
    std::array<bool,3> is_periodic{ {true, true, true} };
    NeighborList neighbor_list;
    utils::MoleculeConnection myosinIndicesPerActin;
    utils::MoleculeConnection actinIndicesPerMyosin;
    utils::MoleculeConnection actinIndicesPerActin;

    std::vector<std::vector<int>> actin_actin_bonds, actin_actin_bonds_prev;
    std::vector<std::vector<int>> actin_actin_status, actin_actin_status_prev;
    // Track lifetime (in steps) for each actin–actin catch bond
    std::vector<std::vector<int>> actin_actin_lifetime, actin_actin_lifetime_prev;
    std::vector<std::vector<int>> kmc_break_flag;
    // Track recovery windows and attachment timestamps for broken bonds
    std::vector<std::vector<size_t>> actin_recovery_until;
    std::vector<std::vector<int>> aa_attach_step;
    // Track actin–myosin bonds
    std::vector<std::vector<int>> am_bonds, am_bonds_prev;
    vector box;
    double k_aa, kappa_aa, k_on, k_off,
           kappa_am, k_am, v_am, myosin_radius_ratio,
           skin_distance, cutoff_radius, dt, base_lifetime, directional_base_lifetime,
           lifetime_coeff, diff_coeff_ratio;
    double am_cutoff, am_optimal;
    double aa_cutoff, aa_optimal;
    double k_mm = 0.0;
    double myomesin_optimal = 0.0;
    double myomesin_cutoff = 0.0;
    double titin_k, titin_rest_length;
    double k_bundle_max = 0.0;
    int bundle_ramp_steps = 0;
    double max_actin_force, max_myosin_force, max_actin_torque, max_myosin_torque;
    size_t bond_recovery_steps;
    bool directional;
    int fix_myosin;
    int max_myosin_bonds;
    int max_strong_actin_bonds;
    std::vector<std::vector<interaction>> am_interaction;
    vector actin_crosslink_ratio;
    std::vector<vec> actin_crosslink_start;
    std::vector<vec> actin_crosslink_end;
    std::vector<int> actin_n_bonds;
    std::vector<int> actin_strong_cb_count;
    std::vector<int> n_myosins_per_actin;
    std::vector<std::pair<std::vector<int>, std::vector<int>>> actin_neighbors_by_species;
    vector actin_basic_tension;
    std::vector<bool> actin_f_load_computed;
    std::vector<double>* actin_f_load_cb;
    std::vector<std::mutex> actin_f_load_mutex;
        
    gsl_rng* rng;
    std::string filename;
    int initial_seed = 0;

    std::vector<std::vector<vec>> actin_forces_temp,
                                    myosin_forces_temp, myosin_velocities_temp, actin_torques_temp, myosin_torques_temp;
    std::vector<std::vector<int>> actin_cb_status_temp;
    std::vector<std::vector<std::array<double, 2>>> myosin_f_load_temp;
    std::vector<std::vector<double>> cb_breakage_events_temp;
    std::vector<std::vector<double>> aa_completed_lifetimes_temp;
    std::vector<std::vector<double>> actin_chem_entropy_delta_temp;
    std::vector<std::vector<double>> actin_chem_input_delta_temp;
    std::vector<std::vector<double>> actin_chem_binding_delta_temp;
    std::vector<std::vector<double>> actin_chem_unbinding_delta_temp;
    std::vector<std::array<double, 2>> myosin_f_load;
    std::vector<utils::MoleculeConnection> actinIndicesPerMyosin_temp;
    std::vector<gsl_rng*> rng_engines;
    std::vector<std::vector<int>> myosin_bond_matrix;
    bool has_myosin_bond_pairs = false;

    // Record the global simulation step and catch-bond breakage events
    size_t current_step = 0;
    // Flat buffer storing (i, j, step, distance, cos_angle) for each breakage
    std::vector<double> cb_breakage_events;

    // Flat buffer storing (i, j, step, bond_count_i, bond_count_j) for
    // removals triggered by the max_strong_actin_bonds limit
    std::vector<double> cb_limit_events;
    // Completed lifetimes sampled at detachments
    std::vector<double> aa_completed_lifetimes;

    // Constructors & Destructor
    Sarcomere();
    Sarcomere(int& n_actins, int& n_myosins, vector box0, double& actin_length, double& myosin_length,
        double& myosin_radius, double& myosin_radius_ratio, double& aa_cutoff, double& aa_optimal,
        double& k_on,
        double& base_lifetime, double& directional_base_lifetime, double& lifetime_coeff,
        double& diff_coeff_ratio, double& k_aa, double& kappa_aa,
        double& k_am, double& kappa_am, double& k_mm, double& v_am, std::string& filename, gsl_rng* rng, int& seed,
        int& fix_myosin, double& dt, bool& directional, std::string& boundary_condition,
        int max_myosin_bonds, int max_strong_actin_bonds, double max_actin_force_param,
        double max_myosin_force_param, double max_actin_torque_param, double max_myosin_torque_param);
    Sarcomere(int& n_actins, int& n_myosins, vector box0, double& actin_length, double& myosin_length,
        double& myosin_radius, double& am_cutoff, double& am_optimal, double& aa_cutoff, double& aa_optimal,
         double& k_on, double& k_off,
        double& base_lifetime, double& directional_base_lifetime, double& lifetime_coeff,
        double& diff_coeff_ratio, double& k_aa, double& kappa_aa, double& k_am, double& kappa_am, double& k_mm, double& v_am,
        std::string& filename, gsl_rng* rng, int& seed, int& fix_myosin, double& dt, double tau_rec,
        double titin_k, double titin_rest_length, bool& directional, int max_myosin_bonds,
        double max_actin_force_param, double max_myosin_force_param,
        double max_actin_torque_param, double max_myosin_torque_param);
    ~Sarcomere();

    // Public Methods
    void myosin_on_a_lattice();
    void partial_fix(int& n_fixed_myosins);
    void cb();
    void bad_cb();
    void cb_off_angle();
    void am_off_angle();
    void set_myosin_direction_x_noise(double noise_std);
    void single_am();
    void sarcomeric_structure();
    void sarcomeric_structure_tight();
    void set_bundling_parameters(double max_strength, int ramp_steps);
    void update_system();
    void update_system_sterics_only();
    void set_periodicity(const std::array<bool,3>& periodic_axes);
    void new_file();
    void save_state();
    void save_resume_snapshot();
    int load_state(int& n_frames, int frame_index = -1);
    // Debug helper to compute catch-bond statistics for a single frame
    void debug_cb_stats();

private:
    // Private helper methods
    void _update_neighbors();
    void _set_to_zero();
    void _process_actin_myosin_binding(int& i);
    void _process_catch_bonds(int& i);
    void _calc_am_force_velocity(int& i);
    void _apply_titin_forces(int& i);
    void _apply_myomesin_spring(int i, int j, std::vector<vec>& local_myosin_forces);
    void _volume_exclusion();
    void _myosin_exclusion();
    void _myosin_repulsion(int& i, int& j);
    void _actin_repulsion(int& i, int& j);
    int determine_cb_status(int& i, int& j);
    bool _cb_decide(int& i, int& j, int status);
    void _record_chem_transition(int& i, int& j, double affinity, bool binding);
    void compute_actin_f_load(int& i);
    void _set_cb(int& i, int& j, int status);
    void _set_cb(int& i, std::vector<int> indices, std::vector<int> status);
    vec _alignment_torque(const vec& u, double k_bias);
    void _apply_cb_alignment_bias(double& k_theta_bias);
    void _update_myosin_bond_matrix();
    bool _myosin_pair_bonded(int mi, int mj) const;
    double _current_bundle_strength() const;
    void _apply_transverse_bundling(double k_bundle);
    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
        _extract_bonded_pairs(
        const std::vector<std::vector<int>>& actin_actin_bonds,
        const std::vector<std::vector<int>>& actin_actin_status,
        const utils::MoleculeConnection& myosinIndicesPerActin);
    void _enforce_actin_cb_limit();
    void _enforce_myosin_bond_limit();
};

#endif // SARCOMERE_H
