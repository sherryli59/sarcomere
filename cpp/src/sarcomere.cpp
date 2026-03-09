#include "sarcomere.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>
#include <gsl/gsl_randist.h>


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

void Sarcomere::set_bundling_parameters(double max_strength, int ramp_steps) {
    k_bundle_max = (max_strength > 0.0) ? max_strength : 0.0;
    bundle_ramp_steps = (ramp_steps > 0) ? ramp_steps : 0;
}

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
        if (am_interaction[i][j].partial_binding_ratio > EPS || !directional) {
            double partial_ratio = am_interaction[i][j].partial_binding_ratio;
            double binding_ratio = am_interaction[i][j].myosin_binding_ratio;
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
        double normalized_partial_ratio = 3.0 * std::min(am_interaction[i][j].partial_binding_ratio, 1.0/3.0);
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
    vec crosslink_point_i = actin.left_end[i];
    crosslink_point_i += actin.direction[i] * (actin.length * crosslink_i);
    vec crosslink_point_j = actin.left_end[j];
    crosslink_point_j += actin.direction[j] * (actin.length * crosslink_j);

    // Compute geometric metrics using the first binding-zone points as endpoints
    double distance = geometry::segment_segment_distance(
        actin.left_end[i], crosslink_point_i, actin.left_end[j], crosslink_point_j, box, is_periodic);
    printf("distance %f\n",distance);
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
        printf("crosslink set to true\n");

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
