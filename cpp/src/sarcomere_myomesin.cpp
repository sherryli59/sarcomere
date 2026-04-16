#include "sarcomere.h"
#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <string>
#include <gsl/gsl_randist.h>


// Parameterized Constructor
Sarcomere::Sarcomere(int& n_actins, int& n_myosins, vector box0, double& actin_length, double& myosin_length,
        double& myosin_radius, double& am_cutoff, double& am_optimal, double& aa_cutoff, double& aa_optimal,
        double& k_on, double& k_off,
        double& base_lifetime, double& lifetime_coeff, double& diff_coeff_ratio, double& k_aa, double& kappa_aa, double& k_am, double& kappa_am, double& k_mm, double& v_am,
        std::string& filename, gsl_rng* rng, int& seed, int& fix_myosin, double& dt, double tau_rec,
        double titin_k, double titin_rest_length, bool& directional, int max_myosin_bonds,
        double max_actin_force_param, double max_myosin_force_param,
        double max_actin_torque_param, double max_myosin_torque_param,
        const std::array<bool,3>& periodic_axes,
        bool use_autodiff)
            : actin(n_actins, actin_length, box0, rng),
              myosin(n_myosins, myosin_length, myosin_radius, box0, rng),
              myosinIndicesPerActin(n_actins),
              actinIndicesPerMyosin(n_myosins),
              neighbor_list(0.0, box0, 0.0, periodic_axes),
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
                max_myosin_torque(max_myosin_torque_param),
                use_autodiff_forces(use_autodiff)

            {
            is_periodic = periodic_axes;
            actin.set_periodic_axes(is_periodic);
            myosin.set_periodic_axes(is_periodic);
            actin.initialize_within_box(rng);
            myosin.initialize_within_box(rng);
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
            myomesin_optimal = 2.0 * this->am_optimal;
            myomesin_cutoff = (k_mm > 0.0) ? myomesin_optimal * 1.1 : 0.0;
            if (tau_rec > 0) {
                bond_recovery_steps = static_cast<size_t>(std::ceil(tau_rec / dt));
            } else {
                bond_recovery_steps = 0;
            }
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
            actin_tension_history.resize(n_actins);
            actin_last_delta_pos.assign(n_actins, vec{0.0, 0.0, 0.0});
            actin_last_delta_rot.assign(n_actins, vec{0.0, 0.0, 0.0});
            myosin_last_delta_pos.assign(myosin.n, vec{0.0, 0.0, 0.0});
            myosin_last_delta_rot.assign(myosin.n, vec{0.0, 0.0, 0.0});
            actin_crosslink_ratio.resize(n_actins);
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

            const char* debug_env = std::getenv("SARCOMERE_DEBUG_ACTIN_FORCES");
            if (debug_env && *debug_env) {
                std::stringstream ss(debug_env);
                std::string token;
                while (std::getline(ss, token, ',')) {
                    token.erase(std::remove_if(token.begin(), token.end(), [](unsigned char ch){ return std::isspace(ch); }), token.end());
                    if (token.empty()) {
                        continue;
                    }
                    try {
                        int idx = std::stoi(token);
                        if (idx >= 0 && idx < n_actins) {
                            debug_force_actins.push_back(idx);
                        }
                    } catch (const std::exception&) {
                        continue;
                    }
                }
                std::sort(debug_force_actins.begin(), debug_force_actins.end());
                debug_force_actins.erase(std::unique(debug_force_actins.begin(), debug_force_actins.end()), debug_force_actins.end());
                if (!debug_force_actins.empty()) {
                    printf("[Sarcomere] Debug force tracing enabled for actins:");
                    for (int idx : debug_force_actins) {
                        printf(" %d", idx);
                    }
                    printf("\n");
                }
            }
        }

// Destructor
Sarcomere::~Sarcomere() {}


void Sarcomere::partial_fix(int& n_fixed){
    // Now each coordinate has three components: x, y, and z (with z = 0).
    std::vector<vector> myosin_positions;
    myosin_positions = {
        {0, -2, 0}, {0, -1, 0}, {0, 0, 0}, {0, 1, 0}, {0, 2, 0},
        {0, -2.5, 0}, {0, -1.5, 0}, {0, -0.5, 0}, {0, 0.5, 0}, {0, 1.5, 0}
    };
    for (int i = 0; i < n_fixed; i++){
        myosin.center[i].x = myosin_positions[i][0];
        myosin.center[i].y = myosin_positions[i][1];
        myosin.center[i].z = myosin_positions[i][2]; // set z coordinate to 0        
    }
    //set all myosin directions to x-axis
    for (int i = 0; i < myosin.n; i++){
        myosin.direction[i] = {1, 0, 0};
    }
    myosin.update_endpoints();
    update_system();
}

void Sarcomere::cb(){
    // For actin, now include a z coordinate equal to 0.
    std::vector<vector> actin_positions = {
        {0.5, -0.015, 0}, {-0.5, 0.015, 0}
    };
    for (int i = 0; i < actin_positions.size(); i++){
        actin.center[i].x = actin_positions[i][0];
        actin.center[i].y = actin_positions[i][1];
        actin.center[i].z = actin_positions[i][2]; // set z coordinate to 0
    }
    actin.direction[0] = {1, 0, 0};
    actin.direction[1] = {-1, 0, 0};
    
    std::vector<vector> myosin_positions = {
        {-1.33, 0.045, 0}, {-1.33, -0.015, 0}, {1.33, 0.015, 0}, {1.33, -0.045, 0}
    };
    for (int i = 0; i < myosin_positions.size(); i++){
        myosin.center[i].x = myosin_positions[i][0];
        myosin.center[i].y = myosin_positions[i][1];
        myosin.center[i].z = myosin_positions[i][2]; // set z coordinate to 0
        myosin.direction[i] = {1, 0, 0};
    }
}

void Sarcomere::set_myosin_direction_x_noise(double noise_std){
    const double sigma = (noise_std > 0.0) ? noise_std : 0.0;
    for (int i = 0; i < myosin.n; ++i) {
        double dx = 1.0;
        double dy = 0.0;
        double dz = 0.0;
        if (rng != nullptr && sigma > 0.0) {
            dx += gsl_ran_gaussian(rng, sigma);
            dy = gsl_ran_gaussian(rng, sigma);
            dz = gsl_ran_gaussian(rng, sigma);
        }
        vec dir{dx, dy, dz};
        double norm = dir.norm();
        if (norm < EPS) {
            dir = {1.0, 0.0, 0.0};
        } else {
            dir = dir / norm;
        }
        myosin.direction[i] = dir;
    }
    myosin.update_endpoints();
}

void Sarcomere::cb_off_angle(){
    if (actin.n < 2 || myosin.n < 4) {
        return;
    }

    constexpr double ANG2_DEG     = 150.0;
    constexpr double Y_ANCH_LEFT  = -0.015;
    constexpr double Y_ANCH_RIGHT =  0.015;

    const double angle_rad = ANG2_DEG * M_PI / 180.0;
    vec dir0{1.0, 0.0, 0.0};
    vec dir1{std::cos(angle_rad), std::sin(angle_rad), 0.0};
    dir1.normalize();
    printf("dir1: (%f, %f, %f)\n", dir1.x, dir1.y, dir1.z);

    actin.direction[0] = dir0;
    actin.direction[1] = dir1;

    const double half_len = 0.5 * actin.length;
    vec anchor0{0.0, Y_ANCH_LEFT, 0.0};   // treat as left endpoint for actin 0
    vec anchor1{0.0, Y_ANCH_RIGHT, 0.0};  // treat as right endpoint for actin 1

    vec center0 = anchor0 + dir0 * half_len;
    vec center1 = anchor1 + dir1 * half_len;
    printf("Actin 1 center: (%f, %f, %f)\n", center1.x, center1.y, center1.z);
    actin.center[0] = center0;
    actin.center[1] = center1;
    actin.update_endpoints();

    double y0 = actin.right_end_y[0];
    double y1 = actin.right_end_y[1];
    printf("Actin 0 right endpoint y: %f\n", y0);
    printf("Actin 1 right endpoint y: %f\n", y1);
    printf("Actin 0 left endpoint y: %f\n", actin.left_end_y[0]);
    printf("Actin 1 left endpoint y: %f\n", actin.left_end_y[1]);
    std::vector<vec> myosin_positions = {
        { 1.33, y0 - 0.03, 0.0 },
        { 1.33, y0 + 0.03, 0.0 },
        {-1.33, y1 - 0.03, 0.0 },
        {-1.33, y1 + 0.03, 0.0 }
    };


    for (size_t i = 0; i < myosin_positions.size() && i < static_cast<size_t>(myosin.n); ++i) {
        myosin.center[static_cast<int>(i)] = myosin_positions[i];
        myosin.direction[static_cast<int>(i)] = vec{1.0, 0.0, 0.0};
    }
    myosin.update_endpoints();
}

void Sarcomere::am_off_angle(){
    // For actin, now include a z coordinate equal to 0.
    std::vector<vector> actin_positions = {
       {-0.5, 0.015, 0}
    };
    for (int i = 0; i < actin_positions.size(); i++){
        actin.center[i].x = actin_positions[i][0];
        actin.center[i].y = actin_positions[i][1];
        actin.center[i].z = actin_positions[i][2];
    }

    // Actin directions: 0 deg and 160 deg apart
    actin.direction[0] = {-0.9397, 0.3420, 0};  // 160° from the first

    std::vector<vector> myosin_positions = {
        {-1.33, 0, 0}
    };
    for (int i = 0; i < myosin_positions.size(); i++){
        myosin.center[i].x = myosin_positions[i][0];
        myosin.center[i].y = myosin_positions[i][1];
        myosin.center[i].z = myosin_positions[i][2];
        myosin.direction[i] = {1, 0, 0};
    }
}

void Sarcomere::sarcomeric_structure_tight(){
    // set box to encompass all three dimensions
    box[0] = 5.32;
    box[1] = 5.32;
    box[2] = 5.32;

    // myosin heads: three rows at y = –0.32, 0, +0.32
    std::vector<vector> myosin_positions = {
        {-1, -0.03, 0.0}, {1, -0.06, 0.0},
        {-1,  0.03, 0.0}, {1,  0, 0.0},
        {-1,  0.09, 0.0}, {1,  0.06, 0.0}
    };
    for (int i = 0; i < myosin_positions.size(); i++){
        myosin.center[i].x   = myosin_positions[i][0];
        myosin.center[i].y   = myosin_positions[i][1];
        myosin.center[i].z   = myosin_positions[i][2];
        myosin.direction[i] = {1, 0, 0}; // set direction to x-axis
    }
    myosin.update_endpoints();

    std::vector<vector> actin_positions = {
        {-2.16, -0.06, 0.0}, 
        {-2.16, 0.0, 0.0}, 
        {-2.16,  0.06, 0.0}, 
        {0.1, -0.03, 0.0}, 
        {0.1, 0.03, 0.0},
        {0.1,  0.09, 0.0}
    };
    for (int i = 0; i < actin_positions.size(); i++){
        printf("Setting actin %d position to (%f, %f, %f)\n", i, actin_positions[i][0], actin_positions[i][1], actin_positions[i][2]);
        actin.center[i].x   = actin_positions[i][0];
        actin.center[i].y   = actin_positions[i][1];
        actin.center[i].z   = actin_positions[i][2];
        actin.direction[i] = {1, 0, 0}; // set direction to x-axis
    }
    int n = actin_positions.size();
    actin_positions = {
        {-0.1, -0.06, 0.0},
        {-0.1,  0.00, 0.0},
        {-0.1,  0.06, 0.0},
        { 2.16, -0.03, 0.0},
        { 2.16,  0.03, 0.0},
        { 2.16,  0.09, 0.0}
    };
    for (int i = 0; i < actin_positions.size(); i++){
        printf("Setting actin %d position to (%f, %f, %f)\n", i+n, actin_positions[i][0], actin_positions[i][1], actin_positions[i][2]);
        actin.center[i+n].x   = actin_positions[i][0];
        actin.center[i+n].y   = actin_positions[i][1];
        actin.center[i+n].z   = actin_positions[i][2];
        actin.direction[i+n] = {-1, 0, 0}; // set direction to negative x-axis
    }
    actin.update_endpoints();

    update_system();
}

void Sarcomere::sarcomeric_structure(){
    // set box to encompass all three dimensions
    box[0] = 5.32;
    box[1] = 5.32;
    box[2] = 5.32;

    // myosin heads: three rows at y = –0.32, 0, +0.32
    std::vector<vector> myosin_positions = {
        {-1.33, -0.03, 0.0}, {1.33, -0.06, 0.0},
        {-1.33,  0.03, 0.0}, {1.33,  0, 0.0},
        {-1.33,  0.09, 0.0}, {1.33,  0.06, 0.0}
    };
    for (int i = 0; i < myosin_positions.size(); i++){
        myosin.center[i].x   = myosin_positions[i][0];
        myosin.center[i].y   = myosin_positions[i][1];
        myosin.center[i].z   = myosin_positions[i][2];
        myosin.direction[i] = {1, 0, 0}; // set direction to x-axis
    }
    myosin.update_endpoints();

    std::vector<vector> actin_positions = {
        {-2.16, -0.06, 0.0}, 
        {-2.16, 0.0, 0.0}, 
        {-2.16,  0.06, 0.0}, 
        {0.5, -0.03, 0.0}, 
        {0.5, 0.03, 0.0},
        {0.5,  0.09, 0.0}
    };
    for (int i = 0; i < actin_positions.size(); i++){
        printf("Setting actin %d position to (%f, %f, %f)\n", i, actin_positions[i][0], actin_positions[i][1], actin_positions[i][2]);
        actin.center[i].x   = actin_positions[i][0];
        actin.center[i].y   = actin_positions[i][1];
        actin.center[i].z   = actin_positions[i][2];
        actin.direction[i] = {1, 0, 0}; // set direction to x-axis
    }

    int n = actin_positions.size();
    actin_positions = {
        {-0.5, -0.06, 0.0},
        {-0.5,  0.00, 0.0},
        {-0.5,  0.06, 0.0},
        { 2.16, -0.03, 0.0},
        { 2.16,  0.03, 0.0},
        { 2.16,  0.09, 0.0}
    };
    for (int i = 0; i < actin_positions.size(); i++){
        printf("Setting actin %d position to (%f, %f, %f)\n", i+n, actin_positions[i][0], actin_positions[i][1], actin_positions[i][2]);
        actin.center[i+n].x   = actin_positions[i][0];
        actin.center[i+n].y   = actin_positions[i][1];
        actin.center[i+n].z   = actin_positions[i][2];
        actin.direction[i+n] = {-1, 0, 0}; // set direction to negative x-axis
    }
    actin.update_endpoints();

    update_system();
}


void Sarcomere::update_system() {
    // Advance global step counter each time the system is updated
    current_step++;
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
            _enforce_myosin_bond_limit();
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
	_apply_wall_forces();
	// double k_theta = 1.0;
	// _apply_cb_alignment_bias(k_theta);

	    #pragma omp barrier

	    if (!debug_force_actins.empty()) {
	        #pragma omp single
	        {
	            const double FORCE_EPS = 1e-9;
	            for (int act_idx : debug_force_actins) {
	                if (act_idx < 0 || act_idx >= actin.n) {
	                    continue;
	                }
	                double total_fx = 0.0, total_fy = 0.0, total_fz = 0.0;
	                double total_tx = 0.0, total_ty = 0.0, total_tz = 0.0;
	                printf("[ForceTrace] step %zu actin %d cb_status=%d tension=%f\n",
	                       current_step, act_idx, actin.cb_status[act_idx],
	                       (act_idx < static_cast<int>(actin_basic_tension.size()) ? actin_basic_tension[act_idx] : 0.0));
	                for (size_t t = 0; t < actin_forces_temp.size(); ++t) {
	                    if (act_idx >= static_cast<int>(actin_forces_temp[t].size())) {
	                        continue;
	                    }
	                    const vec& f = actin_forces_temp[t][act_idx];
	                    const vec& torque = actin_torques_temp[t][act_idx];
	                    bool has_force = (std::fabs(f.x) > FORCE_EPS || std::fabs(f.y) > FORCE_EPS || std::fabs(f.z) > FORCE_EPS);
	                    bool has_torque = (std::fabs(torque.x) > FORCE_EPS || std::fabs(torque.y) > FORCE_EPS || std::fabs(torque.z) > FORCE_EPS);
	                    if (has_force || has_torque) {
	                        printf("  thread %zu actin force=(% .6e,% .6e,% .6e) torque=(% .6e,% .6e,% .6e)\n",
	                               t, f.x, f.y, f.z, torque.x, torque.y, torque.z);
	                    }
	                    total_fx += f.x;
	                    total_fy += f.y;
	                    total_fz += f.z;
	                    total_tx += torque.x;
	                    total_ty += torque.y;
	                    total_tz += torque.z;
	                }
                printf("  total actin force=(% .6e,% .6e,% .6e) torque=(% .6e,% .6e,% .6e)\n",
                       total_fx, total_fy, total_fz, total_tx, total_ty, total_tz);
                vec center = (static_cast<size_t>(act_idx) < actin.center.size()) ? actin.center[act_idx] : vec{0.0, 0.0, 0.0};
                vec direction = (static_cast<size_t>(act_idx) < actin.direction.size()) ? actin.direction[act_idx] : vec{0.0, 0.0, 0.0};
                vec delta_pos = (static_cast<size_t>(act_idx) < actin_last_delta_pos.size()) ? actin_last_delta_pos[act_idx] : vec{0.0, 0.0, 0.0};
                vec delta_rot = (static_cast<size_t>(act_idx) < actin_last_delta_rot.size()) ? actin_last_delta_rot[act_idx] : vec{0.0, 0.0, 0.0};
                printf("  center=(% .6e,% .6e,% .6e) direction=(% .6e,% .6e,% .6e)\n",
                       center.x, center.y, center.z, direction.x, direction.y, direction.z);
                printf("  delta_pos=(% .6e,% .6e,% .6e) delta_rot=(% .6e,% .6e,% .6e)\n",
                       delta_pos.x, delta_pos.y, delta_pos.z, delta_rot.x, delta_rot.y, delta_rot.z);
                const auto& attached_myosins = myosinIndicesPerActin.getConnections(act_idx);
                if (!attached_myosins.empty()) {
                    for (int myo_idx : attached_myosins) {
                        if (myo_idx < 0 || myo_idx >= myosin.n) {
                            continue;
	                        }
	                        double myo_fx = 0.0, myo_fy = 0.0, myo_fz = 0.0;
	                        double myo_tx = 0.0, myo_ty = 0.0, myo_tz = 0.0;
	                        printf("    attached myosin %d\n", myo_idx);
	                        for (size_t t = 0; t < myosin_forces_temp.size(); ++t) {
	                            if (myo_idx >= static_cast<int>(myosin_forces_temp[t].size())) {
	                                continue;
	                            }
	                            const vec& f = myosin_forces_temp[t][myo_idx];
	                            const vec& torque = myosin_torques_temp[t][myo_idx];
	                            bool has_force = (std::fabs(f.x) > FORCE_EPS || std::fabs(f.y) > FORCE_EPS || std::fabs(f.z) > FORCE_EPS);
	                            bool has_torque = (std::fabs(torque.x) > FORCE_EPS || std::fabs(torque.y) > FORCE_EPS || std::fabs(torque.z) > FORCE_EPS);
	                            if (has_force || has_torque) {
	                                printf("      thread %zu myosin force=(% .6e,% .6e,% .6e) torque=(% .6e,% .6e,% .6e)\n",
	                                       t, f.x, f.y, f.z, torque.x, torque.y, torque.z);
	                            }
	                            myo_fx += f.x;
	                            myo_fy += f.y;
	                            myo_fz += f.z;
                            myo_tx += torque.x;
                            myo_ty += torque.y;
                            myo_tz += torque.z;
                        }
                        printf("      total myosin force=(% .6e,% .6e,% .6e) torque=(% .6e,% .6e,% .6e)\n",
                               myo_fx, myo_fy, myo_fz, myo_tx, myo_ty, myo_tz);
                        vec my_center = (static_cast<size_t>(myo_idx) < myosin.center.size()) ? myosin.center[myo_idx] : vec{0.0, 0.0, 0.0};
                        vec my_direction = (static_cast<size_t>(myo_idx) < myosin.direction.size()) ? myosin.direction[myo_idx] : vec{0.0, 0.0, 0.0};
                        vec my_delta_pos = (static_cast<size_t>(myo_idx) < myosin_last_delta_pos.size()) ? myosin_last_delta_pos[myo_idx] : vec{0.0, 0.0, 0.0};
                        vec my_delta_rot = (static_cast<size_t>(myo_idx) < myosin_last_delta_rot.size()) ? myosin_last_delta_rot[myo_idx] : vec{0.0, 0.0, 0.0};
                        printf("      center=(% .6e,% .6e,% .6e) direction=(% .6e,% .6e,% .6e)\n",
                               my_center.x, my_center.y, my_center.z,
                               my_direction.x, my_direction.y, my_direction.z);
                        printf("      delta_pos=(% .6e,% .6e,% .6e) delta_rot=(% .6e,% .6e,% .6e)\n",
                               my_delta_pos.x, my_delta_pos.y, my_delta_pos.z,
                               my_delta_rot.x, my_delta_rot.y, my_delta_rot.z);
                        vec act_left = actin.left_end[act_idx];
                        vec act_right = actin.right_end[act_idx];
                        vec my_left = myosin.left_end[myo_idx];
                        vec my_right = myosin.right_end[myo_idx];
                        double am_distance = geometry::segment_segment_distance(
                            act_left, act_right, my_left, my_right, box, is_periodic);
                        printf("am_distance=% .6e\n", am_distance);
                    }
                } else {
                    printf("    no myosin attachments recorded for actin %d\n", act_idx);
                }
            }
	        }
	    }

        // Step 8: Reduce actin forces and angular forces
        reduce_array(actin_forces_temp, actin.force);
        reduce_array(actin_torques_temp, actin.torque);

        // Step 9: Reduce myosin forces, velocities, and angular forces
        reduce_array(myosin_forces_temp, myosin.force);
        reduce_array(myosin_velocities_temp, myosin.velocity);
        reduce_array(myosin_torques_temp, myosin.torque);

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
    _update_neighbors();
    #pragma omp parallel
    {   
        _set_to_zero();  
        #pragma omp barrier  
        _myosin_exclusion();
        _apply_wall_forces();
        #pragma omp barrier  
        // Step 7: Reduce actin forces and angular forces
        reduce_array(actin_forces_temp, actin.force);
        reduce_array(actin_torques_temp, actin.torque);
        // Step 8: Reduce myosin forces, velocities, and angular forces
        reduce_array(myosin_forces_temp, myosin.force);
        reduce_array(myosin_torques_temp, myosin.torque);
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

void Sarcomere::set_wall_parameters(double strength, double cutoff, double exponent) {
    wall_k = std::max(0.0, strength);
    wall_exponent = (exponent > 0.0) ? exponent : 2.0;
    double default_cutoff = (myosin.radius > 0.0) ? 0.2 * myosin.radius : 0.1;
    wall_cutoff = (cutoff > 0.0) ? cutoff : default_cutoff;
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
            am_cutoff, box, is_periodic);
        if (am_interaction[i][j].myosin_binding_ratio > 0) {
            double partial_ratio = am_interaction[i][j].partial_binding_ratio;
            double binding_ratio = am_interaction[i][j].myosin_binding_ratio;
            if (!directional || partial_ratio > EPS) {
                if (actin_crosslink_ratio[i] > am_interaction[i][j].crosslinkable_ratio) {
                    actin_crosslink_ratio[i] = am_interaction[i][j].crosslinkable_ratio;
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
                double partial = am_interaction[i][j].partial_binding_ratio;
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

    double current_tension = actin_basic_tension[i];
    auto& tension_hist = actin_tension_history[i];
    if (tension_hist.empty() || tension_hist.back().first != current_step) {
        tension_hist.emplace_back(current_step, current_tension);
    } else {
        tension_hist.back().second = current_tension;
    }

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
        std::vector<int> bound_actins;
        int cb2_count = 0;
        for (int i = 0; i < actin.n; ++i) {
            if (am_bonds[i][j] != 1) {
                continue;
            }
            if (actin.cb_status[i] == 2) {
                ++cb2_count;
            } else {
                bound_actins.push_back(i);
            }
        }

        int allowed_non_cb = std::max(0, max_myosin_bonds - cb2_count);
        if (static_cast<int>(bound_actins.size()) <= allowed_non_cb) {
            continue;
        }

        std::vector<int> priority(bound_actins.size(), 0);
        for (size_t idx = 0; idx < bound_actins.size(); ++idx) {
            int i = bound_actins[idx];
            bool cb = actin.cb_status[i] > 1;
            bool prev = am_bonds_prev[i][j] == 1;
            if (cb && prev) priority[idx] = 3;
            else if (cb) priority[idx] = 2;
            else if (prev) priority[idx] = 1;
        }

        auto order = utils::sort_indices(priority);
        for (size_t k = static_cast<size_t>(allowed_non_cb); k < order.size(); ++k) {
            int actin_idx = bound_actins[order[k]];
            am_bonds[actin_idx][j] = 0;
            myosinIndicesPerActin.deleteConnection(actin_idx, j);
            actinIndicesPerMyosin.deleteConnection(j, actin_idx);
            if (actin_idx >= 0 && actin_idx < static_cast<int>(n_myosins_per_actin.size())) {
                n_myosins_per_actin[actin_idx] = std::max(0, n_myosins_per_actin[actin_idx] - 1);
            }
            for (auto& temp_conn : actinIndicesPerMyosin_temp) {
                temp_conn.deleteConnection(j, actin_idx);
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
        // scale by angle between actin and myosin
        double abs_cos_angle = std::abs(actin.direction[i].dot(myosin.direction[j]));
        double normalized_partial_ratio = 3.0 * std::min(am_interaction[i][j].partial_binding_ratio, 1.0/3.0);
        if (abs_cos_angle < EPS || normalized_partial_ratio < EPS) {
            continue;
        }
        vector force_vec = use_autodiff_forces
            ? compute_am_force_and_energy_autodiff(
                  actin, myosin, i, j, box, k_am, kappa_am, am_cutoff, am_optimal)
            : compute_am_force_and_energy(
                  actin, myosin, i, j, box, k_am, kappa_am, am_cutoff, am_optimal);
        local_actin_forces[i].x += force_vec[0];
        local_actin_forces[i].y += force_vec[1];
        local_actin_forces[i].z += force_vec[2];
        local_myosin_forces[j].x -= force_vec[0];
        local_myosin_forces[j].y -= force_vec[1];
        local_myosin_forces[j].z -= force_vec[2];
        vec tau_i{force_vec[3], force_vec[4], force_vec[5]};
        vec tau_j{force_vec[6], force_vec[7], force_vec[8]};
        vec dir_i = actin.direction[i];
        dir_i.normalize();
        vec dir_j = myosin.direction[j];
        dir_j.normalize();
        tau_i -= dir_i * tau_i.dot(dir_i);
        tau_j -= dir_j * tau_j.dot(dir_j);
        local_actin_torques[i].x += tau_i.x;
        local_actin_torques[i].y += tau_i.y;
        local_actin_torques[i].z += tau_i.z;
        local_myosin_torques[j].x += tau_j.x;
        local_myosin_torques[j].y += tau_j.y;
        local_myosin_torques[j].z += tau_j.z;
        // double torque_i_mag = std::sqrt(tau_i.x*tau_i.x + tau_i.y*tau_i.y + tau_i.z*tau_i.z);
        // double torque_j_mag = std::sqrt(tau_j.x*tau_j.x + tau_j.y*tau_j.y + tau_j.z*tau_j.z);
        // printf("Actin %d and Myosin %d: torque magnitudes (%f, %f)\n", i, j, torque_i_mag, torque_j_mag);
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
        for (int j : myosin_neighbors) {
            double abs_cos_angle = std::abs(actin.direction[i].dot(myosin.direction[j]));
            if (am_bonds[i][j] == 1 | abs_cos_angle < 0.95) {
                continue;
            }
            apply_actin_myosin_repulsion(
                actin,
                myosin,
                i,
                j,
                box,
                am_cutoff*1.1,
                k_aa,
                local_actin_forces[i],
                local_myosin_forces[j]);
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
    double dist = diff.norm();
    if (dist <= 1e-12 || dist > myomesin_cutoff) {
        return;
    }
    vec unit = diff / dist;
    double stretch = dist - myomesin_optimal;
    double force_scalar = -k_mm * stretch;
    if (std::isfinite(max_myosin_force) && max_myosin_force > 0.0) {
        force_scalar = std::clamp(force_scalar, -max_myosin_force, max_myosin_force);
    }
    vec force_on_i = force_scalar * unit;
    vec force_on_j = -force_on_i;

    if (i < fix_myosin && j < fix_myosin) {
        return;
    }
    if (i < fix_myosin) {
        local_myosin_forces[j] += 2.0 * force_on_j;
        return;
    }
    if (j < fix_myosin) {
        local_myosin_forces[i] += 2.0 * force_on_i;
        return;
    }
    local_myosin_forces[i] += force_on_i;
    local_myosin_forces[j] += force_on_j;
}


void Sarcomere::_volume_exclusion(){
    const double EPS_FORCE = 1e-9;
    double myosin_cutoff = 2.0 * myosin.radius;
    const double myomesin_distance_limit = 2.0 * am_optimal * 1.1;
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
                    10 * k_aa,
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

void Sarcomere::_apply_wall_forces() {
    if (wall_k <= 0.0 || wall_cutoff <= 0.0) {
        return;
    }
    if (is_periodic[0] && is_periodic[1] && is_periodic[2]) {
        return;
    }

    auto apply_to_filament = [&](Filament& filament,
                                 std::vector<std::vector<vec>>& force_buffer,
                                 std::vector<std::vector<vec>>& torque_buffer,
                                 bool skip_fixed,
                                 int fixed_count) {
        #pragma omp for schedule(static)
        for (int idx = 0; idx < filament.n; ++idx) {
            if (skip_fixed && idx < fixed_count) {
                continue;
            }
            int thread_id = omp_get_thread_num();
            auto& local_forces = force_buffer[thread_id];
            auto& local_torques = torque_buffer[thread_id];

            vec left = filament.left_end[idx];
            vec right = filament.right_end[idx];
            vec center = filament.center[idx];
            vec segment = right - left;

            for (int axis = 0; axis < 3; ++axis) {
                if (is_periodic[axis]) {
                    continue;
                }
                if (axis >= static_cast<int>(box.size()) || box[axis] <= 0.0) {
                    continue;
                }
                double half_length = 0.5 * box[axis];
                double lower = -half_length;
                double upper = half_length;

                auto apply_plane = [&](double plane_coord, double normal_sign) {
                    vec normal{0.0, 0.0, 0.0};
                    if (axis == 0) normal.x = normal_sign;
                    else if (axis == 1) normal.y = normal_sign;
                    else normal.z = normal_sign;

                    double left_coord = (axis == 0) ? left.x : (axis == 1 ? left.y : left.z);
                    double right_coord = (axis == 0) ? right.x : (axis == 1 ? right.y : right.z);
                    double delta = right_coord - left_coord;
                    double clearance_left = (normal_sign > 0.0) ? (left_coord - plane_coord) : (plane_coord - left_coord);
                    double clearance_right = (normal_sign > 0.0) ? (right_coord - plane_coord) : (plane_coord - right_coord);
                    double clearance = 0.0;
                    vec q;

                    if (clearance_left <= 0.0 && clearance_right <= 0.0) {
                        double t = 0.0;
                        if (std::abs(delta) > 1e-12) {
                            t = (plane_coord - left_coord) / delta;
                        }
                        t = std::clamp(t, 0.0, 1.0);
                        q = left + segment * t;
                        clearance = 0.0;
                    } else if (clearance_left <= 0.0 || clearance_right <= 0.0) {
                        double t = 0.0;
                        if (std::abs(delta) > 1e-12) {
                            t = (plane_coord - left_coord) / delta;
                        }
                        t = std::clamp(t, 0.0, 1.0);
                        q = left + segment * t;
                        double q_coord = (axis == 0) ? q.x : (axis == 1 ? q.y : q.z);
                        clearance = (normal_sign > 0.0)
                            ? std::max(0.0, q_coord - plane_coord)
                            : std::max(0.0, plane_coord - q_coord);
                    } else {
                        if (clearance_left < clearance_right) {
                            clearance = clearance_left;
                            q = left;
                        } else {
                            clearance = clearance_right;
                            q = right;
                        }
                    }

                    clearance = std::max(0.0, clearance);
                    if (clearance >= wall_cutoff) {
                        return;
                    }
                    double ratio = (wall_cutoff - clearance) / wall_cutoff;
                    double magnitude = wall_k * std::pow(ratio, wall_exponent);
                    vec force = magnitude * normal;
                    local_forces[idx] += force;
                    vec arm = q - center;
                    local_torques[idx] += arm.cross(force);
                };

                apply_plane(lower, +1.0);
                apply_plane(upper, -1.0);
            }
        }
    };

    apply_to_filament(myosin, myosin_forces_temp, myosin_torques_temp, true, fix_myosin);
    apply_to_filament(actin, actin_forces_temp, actin_torques_temp, false, 0);
}

void Sarcomere::_myosin_exclusion(){
    const double EPS_FORCE = 1e-9;
    double myosin_cutoff = 2.0 * myosin.radius;
    const double myomesin_distance_limit = 2.0 * am_optimal * 1.1;
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
                    2*k_am,
                    local_myosin_forces[i],
                    local_myosin_forces[j],
                    seg_distance);
                if (seg_distance <= myomesin_distance_limit) {
                    _apply_myomesin_spring(i, j, local_myosin_forces);
                }
            }
        }
    }
}

int Sarcomere::determine_cb_status(int& i, int& j){
    double crosslink_i = std::clamp(actin["crosslink_ratio"][i], 0.0, 1.0);
    double crosslink_j = std::clamp(actin["crosslink_ratio"][j], 0.0, 1.0);
    if (crosslink_i <= EPS || crosslink_j <= EPS) {
        return 0;
    }

    vec crosslink_point_i = actin.left_end[i];
    crosslink_point_i += actin.direction[i] * (actin.length * crosslink_i);
    vec crosslink_point_j = actin.left_end[j];
    crosslink_point_j += actin.direction[j] * (actin.length * crosslink_j);

    // Compute geometric metrics using the first binding-zone points as endpoints
    double distance = geometry::segment_segment_distance(
        actin.left_end[i], crosslink_point_i, actin.left_end[j], crosslink_point_j, box, is_periodic);
    double cos_angle = actin.direction[i].dot(actin.direction[j]);

    bool was_strong = (actin_actin_status_prev[i][j] == 2);
    int thread_id = omp_get_thread_num();
    auto& local_breakage_events = cb_breakage_events_temp[thread_id];
    auto& local_completed_lifetimes = aa_completed_lifetimes_temp[thread_id];

    struct TensionSnapshot {
        double prev_value;
        double curr_value;
        size_t prev_step;
        size_t curr_step;
    };

    auto get_tension_snapshot = [&](int idx) -> TensionSnapshot {
        TensionSnapshot snap{actin_basic_tension[idx], actin_basic_tension[idx],
                             (current_step > 0) ? current_step - 1 : current_step,
                             current_step};
        const auto& history = actin_tension_history[idx];
        if (!history.empty()) {
            const auto& last = history.back();
            snap.curr_step = last.first;
            snap.curr_value = last.second;
            if (history.size() >= 2) {
                const auto& prev_entry = history[history.size() - 2];
                snap.prev_step = prev_entry.first;
                snap.prev_value = prev_entry.second;
            } else {
                snap.prev_step = snap.curr_step;
                snap.prev_value = snap.curr_value;
            }
        }
        return snap;
    };

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
    if ((crosslink_i > EPS) && (crosslink_j > EPS) || !directional) {
        if (distance < aa_cutoff) {
            crosslink = true;
        }
    }
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
            TensionSnapshot tension_i = get_tension_snapshot(i);
            TensionSnapshot tension_j = get_tension_snapshot(j);
            printf("Actins %d and %d no longer form catch bonds, distance: %f, cos_angle: %f, tensions: %f (prev %f @ step %zu), %f (prev %f @ step %zu)\n",
                   i, j, distance, cos_angle,
                   tension_i.curr_value, tension_i.prev_value, tension_i.prev_step,
                   tension_j.curr_value, tension_j.prev_value, tension_j.prev_step);
            record_break();
        }
        return 1;
    }
    auto& myosin_indices_i = myosinIndicesPerActin.getConnections(i);
    auto& myosin_indices_j = myosinIndicesPerActin.getConnections(j);
    if (myosin_indices_i.empty() || myosin_indices_j.empty()){
        if (was_strong){
            TensionSnapshot tension_i = get_tension_snapshot(i);
            TensionSnapshot tension_j = get_tension_snapshot(j);
            printf("Actins %d and %d no longer form catch bonds (no myosin), distance: %f, cos_angle: %f, tensions: %f (prev %f @ step %zu), %f (prev %f @ step %zu)\n",
                   i, j, distance, cos_angle,
                   tension_i.curr_value, tension_i.prev_value, tension_i.prev_step,
                   tension_j.curr_value, tension_j.prev_value, tension_j.prev_step);
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
        TensionSnapshot tension_i = get_tension_snapshot(i);
        TensionSnapshot tension_j = get_tension_snapshot(j);
        printf("Actins %d and %d no longer form catch bonds (no shared myosin), distance: %f, cos_angle: %f, tensions: %f (prev %f @ step %zu), %f (prev %f @ step %zu)\n",
               i, j, distance, cos_angle,
               tension_i.curr_value, tension_i.prev_value, tension_i.prev_step,
               tension_j.curr_value, tension_j.prev_value, tension_j.prev_step);
        record_break();
    }
    return 1;
}


bool Sarcomere::_cb_decide(int& i, int& j, int status){
    if (status == 0){
        return false;
    }
    bool was_bonded_prev = (actin_actin_bonds_prev[i][j] == 1);
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
            // if (f_load>0) {
            // printf("k_off_adjusted: %f, rand: %f, f_load: %f, abs_cos_angle: %f, lifetime: %f\n",
            //        k_off_adjusted, rand, f_load, abs_cos_angle,
            //        base_lifetime + lifetime_coeff * f_load);
            // printf("actual lifetime: %f\n", (current_step - aa_attach_step[i][j]) * dt);
            //     }
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
            return false;
        }
            // Bond is accepted or persists: record status locally
        if (status == 2){
            printf("Forming strong actin-actin bond between %d and %d with status %d (rand %f)\n", i, j, status, rand);
        }
    }
  
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
 
    vector force_vec = use_autodiff_forces
        ? compute_aa_force_and_energy_autodiff(actin, i, j, box, k_aa, kappa_aa, aa_cutoff, aa_optimal)
        : compute_aa_force_and_energy(actin, i, j, box, k_aa, kappa_aa, aa_cutoff, aa_optimal);
    local_actin_forces[i].x += force_vec[0];
    local_actin_forces[i].y += force_vec[1];
    local_actin_forces[i].z += force_vec[2];
    local_actin_forces[j].x -= force_vec[0];
    local_actin_forces[j].y -= force_vec[1];
    local_actin_forces[j].z -= force_vec[2];
    vec tau_i{force_vec[3], force_vec[4], force_vec[5]};
    vec tau_j{force_vec[6], force_vec[7], force_vec[8]};
    vec dir_i = actin.direction[i];
    dir_i.normalize();
    vec dir_j = actin.direction[j];
    dir_j.normalize();
    tau_i -= dir_i * tau_i.dot(dir_i);
    tau_j -= dir_j * tau_j.dot(dir_j);
    local_actin_torques[i].x += tau_i.x;
    local_actin_torques[i].y += tau_i.y;
    local_actin_torques[i].z += tau_i.z;
    local_actin_torques[j].x += tau_j.x;
    local_actin_torques[j].y += tau_j.y;
    local_actin_torques[j].z += tau_j.z;
    // // Print torque magnitudes for the two actins (thread-local accumulated torques)
    // double torque_i_mag = std::sqrt(
    //     local_actin_torques[i].x * local_actin_torques[i].x +
    //     local_actin_torques[i].y * local_actin_torques[i].y +
    //     local_actin_torques[i].z * local_actin_torques[i].z);
    // double torque_j_mag = std::sqrt(
    //     local_actin_torques[j].x * local_actin_torques[j].x +
    //     local_actin_torques[j].y * local_actin_torques[j].y +
    //     local_actin_torques[j].z * local_actin_torques[j].z);
    // double am_angle = actin.direction[i].dot(actin.direction[j]);
    // printf("Actin %d torque magnitude: %f, Actin %d torque magnitude: %f, dot prod %f\n",
    //        i, torque_i_mag, j, torque_j_mag, am_angle);
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
    // Align orientations toward the +/- x-axis using the same energy gradient as orientation_evolution.
    const vec target = {1.0, 0.0, 0.0};
    double dot = std::clamp(u.dot(target), -1.0, 1.0);
    vec cross_ut = u.cross(target);
    return k_bias * dot * cross_ut;
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

int Sarcomere::load_state(int& n_frames, int frame_index){
    int target_frame = load_from_file(filename, actin, myosin, actin_actin_bonds, n_frames, frame_index);
    update_system();
    return target_frame;
}
