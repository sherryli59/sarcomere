#ifndef SARCOMERE_H
#define SARCOMERE_H

#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <omp.h>

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
    NeighborList neighbor_list;
    utils::MoleculeConnection myosinIndicesPerActin;
    utils::MoleculeConnection actinIndicesPerMyosin;
    utils::MoleculeConnection actinIndicesPerActin;

    std::vector<std::vector<int>> actin_actin_bonds, actin_actin_bonds_prev;
    // Track lifetime (in steps) for each actin–actin catch bond
    std::vector<std::vector<int>> actin_actin_lifetime, actin_actin_lifetime_prev;
    // Store normalized strength for each actin–actin bond
    std::vector<std::vector<double>> actin_actin_strength;
    // Track actin–myosin bonds for kinetic Monte Carlo updates
    std::vector<std::vector<int>> am_bonds, am_bonds_prev;
    // Store previous actin load used for k_off calculations
    std::vector<double> actin_f_load_prev;
    vector box;
    double k_aa, kappa_aa, cb_mult_factor, k_on, k_off,
           kappa_am, k_am, v_am, crosslinker_length, myosin_radius_ratio, skin_distance, cutoff_radius,
           dt, base_lifetime, lifetime_coeff, diff_coeff_ratio;
    bool directional;
    int fix_myosin;
    int max_myosin_bonds;
    int max_strong_actin_bonds;
    std::vector<std::vector<interaction>> am_interaction;
    vector actin_crosslink_ratio;
    std::vector<int> actin_n_bonds;
    std::vector<int> n_myosins_per_actin;
    std::vector<std::pair<std::vector<int>, std::vector<int>>> actin_neighbors_by_species;
    vector actin_basic_tension;
        
    gsl_rng* rng;
    std::string filename;

    std::vector<std::vector<vec>> actin_forces_temp, 
                                    myosin_forces_temp, myosin_velocities_temp, actin_torques_temp, myosin_torques_temp;
    std::vector<std::vector<double>> actin_cb_strengths_temp, myosin_f_load_temp;
    std::vector<utils::MoleculeConnection> actinIndicesPerMyosin_temp;
    std::vector<gsl_rng*> rng_engines;

    // Constructors & Destructor
    Sarcomere();
    Sarcomere(int& n_actins, int& n_myosins, vector box0, double& actin_length, double& myosin_length,
        double& myosin_radius, double& myosin_radius_ratio, double& crosslinker_length, double& k_on, double& k_off,
        double& base_lifetime, double& lifetime_coeff, double& diff_coeff_ratio, double& k_aa, double& kappa_aa, double& k_am, double& kappa_am, double& v_am,
        std::string& filename, gsl_rng* rng, int& seed, int& fix_myosin, double& dt, bool& directional, int max_myosin_bonds = 2, int max_strong_actin_bonds = 1);
    Sarcomere(int& n_actins, int& n_myosins, vector box0, double& actin_length, double& myosin_length,
        double& myosin_radius, double& myosin_radius_ratio, double& crosslinker_length, double& k_on,
        double& base_lifetime, double& lifetime_coeff, double& diff_coeff_ratio, double& k_aa, double& kappa_aa, double& k_am, double& kappa_am, double& v_am,
        std::string& filename, gsl_rng* rng, int& seed, int& fix_myosin, double& dt, bool& directional, std::string& boundary_condition, int max_myosin_bonds = 2, int max_strong_actin_bonds = 1);
    ~Sarcomere();

    // Public Methods
    void myosin_on_a_lattice();
    void partial_fix(int& n_fixed_myosins);
    void cb();
    void bad_cb();
    void cb_off_angle();
    void am_off_angle();
    void single_am();
    void sarcomeric_structure();
    void update_system();
    void update_system_sterics_only();
    void new_file();
    void save_state();
    void load_state(int& n_frames);
    // Debug helper to compute catch-bond statistics for a single frame
    void debug_cb_stats();

private:
    // Private helper methods
    void _update_neighbors();
    void _set_to_zero();
    void _process_actin_myosin_binding(int& i);
    void _process_catch_bonds(int& i);
    void _calc_am_force_velocity(int& i);
    void _volume_exclusion();
    void _myosin_exclusion();
    void _myosin_repulsion(int& i, int& j);
    void _actin_repulsion(int& i, int& j);
    double _get_cb_strength(int& i, int& j);
    void _get_f_load(int& i);
    void _set_cb(int& i, int& j, double& normalized_strength, bool& add_connection);
    void _set_cb(int& i, std::vector<int> indices, vector cb_strength);
    vec _alignment_torque(const vec& u, double k_bias);
    void _apply_cb_alignment_bias(double& k_theta_bias);
    std::pair<std::vector<double>, std::vector<double>>  _extract_bonded_pairs(
        const std::vector<std::vector<int>>& actin_actin_bonds,
        const utils::MoleculeConnection& myosinIndicesPerActin);
    void _enforce_myosin_bond_limit();
    void _enforce_actin_cb_limit();
};

#endif // SARCOMERE_H
