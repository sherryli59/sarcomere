#include "sarcomere.h"
#include <algorithm>

// Constructor
Sarcomere::Sarcomere() {}

// Parameterized Constructor
Sarcomere::Sarcomere(int& n_actins, int& n_myosins, vector box0, double& actin_length, double& myosin_length,
        double& myosin_radius, double& myosin_radius_ratio, double& crosslinker_length, double& k_on, double& k_off,
        double& base_lifetime, double& lifetime_coeff, double& diff_coeff_ratio, double& k_aa, double& kappa_aa, double& k_am, double& kappa_am, double& v_am,
        std::string& filename, gsl_rng* rng, int& seed, int& fix_myosin, double& dt, bool& directional)
            : actin(n_actins, actin_length, box0, rng),
              myosin(n_myosins, myosin_length, myosin_radius, box0, rng),
              myosinIndicesPerActin(n_actins),
              actinIndicesPerMyosin(n_myosins),
              neighbor_list(0.0, box0, 0.0),
                actin_actin_bonds(n_actins, std::vector<int>(n_actins, 0)),
                actin_actin_lifetime(n_actins, std::vector<int>(n_actins, 0)),
                actin_forces_temp(omp_get_max_threads(), std::vector<vec>(n_actins, {0, 0, 0})),
                myosin_forces_temp(omp_get_max_threads(), std::vector<vec>(n_myosins, {0, 0, 0})),
                myosin_velocities_temp(omp_get_max_threads(), std::vector<vec>(n_myosins, {0, 0, 0})),
                actin_torques_temp(omp_get_max_threads(), std::vector<vec>(n_actins, {0, 0, 0})),
                myosin_torques_temp(omp_get_max_threads(), std::vector<vec>(n_myosins, {0, 0, 0})),
                myosin_f_load_temp(omp_get_max_threads(), std::vector<double>(n_myosins, 0)),
                actin_cb_strengths_temp(omp_get_max_threads(), std::vector<double>(n_actins, 0)),
                actinIndicesPerMyosin_temp(omp_get_max_threads(), utils::MoleculeConnection(n_myosins)),
                rng_engines(omp_get_max_threads(), nullptr) 

            {
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
            this->crosslinker_length = crosslinker_length;
            this->skin_distance = skin_distance;
            this->filename = filename;
            this->rng = rng;
            this->fix_myosin = fix_myosin;
            this->dt = dt;
            this->myosin_radius_ratio = myosin_radius_ratio;
            this->base_lifetime = base_lifetime;
            this->lifetime_coeff = lifetime_coeff;
            this->diff_coeff_ratio = diff_coeff_ratio;
            this->directional = directional;
            cb_mult_factor = 1000;
            cutoff_radius = std::max(actin_length, myosin_length) +
                            std::max(2 * myosin_radius, crosslinker_length);
            double skin_distance = 0.15 * cutoff_radius;
            neighbor_list = NeighborList(cutoff_radius + skin_distance, box, skin_distance / 2);
            neighbor_list.initialize(actin.center_x, actin.center_y, actin.center_z,
                                    myosin.center_x, myosin.center_y, myosin.center_z);
            actin_actin_bonds_prev = actin_actin_bonds;
            actin_actin_lifetime_prev = actin_actin_lifetime;
            actin_basic_tension.resize(n_actins);
            actin_crosslink_ratio.resize(n_actins);
            actin_n_bonds.resize(n_actins);
            n_myosins_per_actin.resize(n_actins);
            am_interaction.resize(n_actins);
            for (int i = 0; i < n_actins; i++) {
                am_interaction[i].resize(n_myosins);
            }
            actin.register_feature("myosin_binding_ratio");
            actin.register_feature("crosslink_ratio");
            actin.register_feature("partial_binding_ratio");
            // Initialize thread-local RNGs
            for (int t = 0; t < omp_get_max_threads(); ++t) {
                rng_engines[t] = gsl_rng_alloc(gsl_rng_mt19937); 
                gsl_rng_set(rng_engines[t], seed + t);           
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
    // set all myosin directions to x-axis
    // for (int i = 0; i < myosin.n; i++){
    //     myosin.direction[i] = {1, 0, 0};
    // }
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
        {-1.33, 0, 0}, {1.33, 0, 0}
    };
    for (int i = 0; i < myosin_positions.size(); i++){
        myosin.center[i].x = myosin_positions[i][0];
        myosin.center[i].y = myosin_positions[i][1];
        myosin.center[i].z = myosin_positions[i][2]; // set z coordinate to 0
        myosin.direction[i] = {1, 0, 0};
    }
}

void Sarcomere::cb_off_angle(){
    // For actin, now include a z coordinate equal to 0.
    std::vector<vector> actin_positions = {
        {0.4, 0, 0}, {-0.4, 0.1, 0}
    };
    for (int i = 0; i < actin_positions.size(); i++){
        actin.center[i].x = actin_positions[i][0];
        actin.center[i].y = actin_positions[i][1];
        actin.center[i].z = actin_positions[i][2];
    }

    // Actin directions: 0 deg and 160 deg apart
    actin.direction[0] = {1, 0, 0};  // Reference direction
    actin.direction[1] = {-0.9397, 0.3420, 0};  // 160° from the first

    std::vector<vector> myosin_positions = {
        {-1.2, 0, 0}, {1.2, 0, 0}
    };
    for (int i = 0; i < myosin_positions.size(); i++){
        myosin.center[i].x = myosin_positions[i][0];
        myosin.center[i].y = myosin_positions[i][1];
        myosin.center[i].z = myosin_positions[i][2];
        myosin.direction[i] = {1, 0, 0};
    }
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


void Sarcomere::sarcomeric_structure(){
    // set box to encompass all three dimensions
    box[0] = 5.32;
    box[1] = 5.32;
    box[2] = 5.32;

    // myosin heads: three rows at y = –0.32, 0, +0.32
    std::vector<vector> myosin_positions = {
        {-1.33, -0.32, 0.0}, { 1.33, -0.32, 0.0},
        {-1.33,  0.00, 0.0}, { 1.33,  0.00, 0.0},
        {-1.33,  0.32, 0.0}, { 1.33,  0.32, 0.0}
    };
    for (int i = 0; i < myosin_positions.size(); i++){
        myosin.center[i].x   = myosin_positions[i][0];
        myosin.center[i].y   = myosin_positions[i][1];
        myosin.center[i].z   = myosin_positions[i][2];
        myosin.direction[i] = {1, 0, 0}; // set direction to x-axis
    }
    myosin.update_endpoints();

    // first half of actins: x = +0.5 and –2.16, nine per column, y ∈ {–0.47,…,+0.47}
    std::vector<vector> actin_positions = {
        { 0.5, -0.47, 0.0}, { 0.5, -0.32, 0.0}, { 0.5, -0.17, 0.0},
        { 0.5, -0.15, 0.0}, { 0.5,  0.00, 0.0}, { 0.5,  0.15, 0.0},
        { 0.5,  0.17, 0.0}, { 0.5,  0.32, 0.0}, { 0.5,  0.47, 0.0},
        {-2.16, -0.47, 0.0}, {-2.16, -0.32, 0.0}, {-2.16, -0.17, 0.0},
        {-2.16, -0.15, 0.0}, {-2.16,  0.00, 0.0}, {-2.16,  0.15, 0.0},
        {-2.16,  0.17, 0.0}, {-2.16,  0.32, 0.0}, {-2.16,  0.47, 0.0}
    };
    for (int i = 0; i < actin_positions.size(); i++){
        actin.center[i].x   = actin_positions[i][0];
        actin.center[i].y   = actin_positions[i][1];
        actin.center[i].z   = actin_positions[i][2];
        actin.direction[i] = {1, 0, 0}; // set direction to x-axis
    }

    // second half of actins: x = –0.5 and +2.16, same nine y’s
    int n = actin_positions.size();
    actin_positions = {
        {-0.5, -0.47, 0.0}, {-0.5, -0.32, 0.0}, {-0.5, -0.17, 0.0},
        {-0.5, -0.15, 0.0}, {-0.5,  0.00, 0.0}, {-0.5,  0.15, 0.0},
        {-0.5,  0.17, 0.0}, {-0.5,  0.32, 0.0}, {-0.5,  0.47, 0.0},
        { 2.16, -0.47, 0.0}, { 2.16, -0.32, 0.0}, { 2.16, -0.17, 0.0},
        { 2.16, -0.15, 0.0}, { 2.16,  0.00, 0.0}, { 2.16,  0.15, 0.0},
        { 2.16,  0.17, 0.0}, { 2.16,  0.32, 0.0}, { 2.16,  0.47, 0.0}
    };
    for (int i = 0; i < actin_positions.size(); i++){
        actin.center[i+n].x   = actin_positions[i][0];
        actin.center[i+n].y   = actin_positions[i][1];
        actin.center[i+n].z   = actin_positions[i][2];
        actin.direction[i+n] = {-1, 0, 0}; // set direction to negative x-axis
    }
    actin.update_endpoints();

    update_system();
}


void Sarcomere::update_system() {
    _update_neighbors();
    #pragma omp parallel
    {   
        _set_to_zero();  
        #pragma omp barrier  
        //Step 2: Compute actin-myosin binding
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < actin.n; i++) {
            _process_actin_myosin_binding(i);
        }

        #pragma omp barrier  

        // Step 3: Compute catch bonds
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < actin.n; i++) {
            _process_catch_bonds(i);
        }

        #pragma omp barrier  

        // Step 4: Reduce actin catch-bond strengths
        utils::reduce_array(actin_cb_strengths_temp, actin.cb_strength);


        #pragma omp barrier  

        // Step 5: Concatenate actinIndicesPerMyosin connections
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

        // Step 6: Compute actin-myosin forces
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < actin.n; i++) {
            _calc_am_force_velocity(i);
        }

        _volume_exclusion();
        double k_theta = 100.0;
        _apply_cb_alignment_bias(k_theta);

        #pragma omp barrier  

        // Step 7: Reduce actin forces and angular forces
        reduce_array(actin_forces_temp, actin.force);
        reduce_array(actin_torques_temp, actin.torque);

        // Step 8: Reduce myosin forces, velocities, and angular forces
        reduce_array(myosin_forces_temp, myosin.force);
        reduce_array(myosin_velocities_temp, myosin.velocity);
        reduce_array(myosin_torques_temp, myosin.torque);
        
    
        //set max velocity for myosin
        //std::vector<double> myosin_min_load(myosin.n, 0);
        // // Parallelize over myosins
        // #pragma omp for
        // for (int myosin_idx = 0; myosin_idx < myosin.n; ++myosin_idx) {
        //     double min_load = myosin_f_load_temp[0][myosin_idx];
        //     // Find max across threads
        //     for (int thread_idx = 1; thread_idx < omp_get_max_threads(); ++thread_idx) {
        //         min_load = std::min(min_load, myosin_f_load_temp[thread_idx][myosin_idx]);
        //     }
        //     myosin_min_load[myosin_idx] = min_load;
        //     printf("myosin %d, min_load: %f\n", myosin_idx, min_load);
        // }

        #pragma omp for
        for (int i = 0; i < myosin.n; i++){
            double v_max = v_am / (diff_coeff_ratio + 1);// * (1 - myosin_min_load[i]);
            //printf("myosin %d, v_max: %f\n", i, v_max);
            double v = myosin.velocity[i].norm();
            if (v>v_max){
                myosin.velocity[i] = myosin.velocity[i]/v*v_max;
            }
        }
    }
}


void Sarcomere::update_system_sterics_only() {
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
    }
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
        for (int i = 0; i < actin.n; ++i) {
            actin_forces_temp[t][i] = {0, 0, 0};
            actin_torques_temp[t][i] = {0, 0, 0};
            actin_cb_strengths_temp[t][i] = 0;
        }
    }

    #pragma omp for
    for (int t = 0; t < omp_get_max_threads(); ++t) {
        for (int i = 0; i < myosin.n; ++i) {
            myosin_forces_temp[t][i] = {0, 0, 0};
            myosin_velocities_temp[t][i] = {0, 0, 0};
            myosin_torques_temp[t][i] = {0, 0, 0};
            myosin_f_load_temp[t][i] = 1;
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
        actin.cb_strength[i] = 0;
        actin_basic_tension[i] = 0;
        actin_n_bonds[i] = 0;
        n_myosins_per_actin[i] = 0;
        actin_crosslink_ratio[i] = 1;
        actin["myosin_binding_ratio"][i] = 0;
        actin["crosslink_ratio"][i] = 1;
        actin["partial_binding_ratio"][i] = 0;
        myosinIndicesPerActin.deleteAllConnections(i);
        for (int j = 0; j < actin.n; j++){
            actin_actin_bonds_prev[i][j] = actin_actin_bonds[i][j];
            actin_actin_bonds[i][j] = 0;
            actin_actin_lifetime_prev[i][j] = actin_actin_lifetime[i][j];
            actin_actin_lifetime[i][j] = 0;
        }
    }

    #pragma omp for
    for (int i = 0; i < myosin.n; i++) {
        myosin.update_endpoints(i);
        myosin.force[i] = {0, 0, 0};
        myosin.torque[i] = {0, 0, 0};
        myosin.velocity[i] = {0, 0, 0};
        actinIndicesPerMyosin.deleteAllConnections(i);
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
    double cutoff = myosin.radius/ myosin_radius_ratio;
    for (int index = 0; index < myosin_neighbors.size(); index++) {
        int j = myosin_neighbors[index];
        
        am_interaction[i][j] = geometry::analyze_am(
            actin.left_end[i], actin.right_end[i], myosin.left_end[j], myosin.right_end[j],
            cutoff, box);

        if (am_interaction[i][j].myosin_binding_ratio>0){
            if (!directional || am_interaction[i][j].partial_binding_ratio > EPS){
                local_actinIndicesPerMyosin.addConnection(j, i);
                myosinIndicesPerActin.addConnection(i, j);
                if (actin_crosslink_ratio[i] > am_interaction[i][j].crosslinkable_ratio){
                    actin_crosslink_ratio[i] = am_interaction[i][j].crosslinkable_ratio;
                }
                if (am_interaction[i][j].partial_binding_ratio>EPS || !directional){
                    n_myosins_per_actin[i]++;
                    double abs_cos = std::abs(actin.direction[i].dot(myosin.direction[j]));
                    if (abs_cos > actin_basic_tension[i]) {
                        actin_basic_tension[i] = abs_cos;
                    }
                    if (actin["partial_binding_ratio"][i] < am_interaction[i][j].partial_binding_ratio) {
                        actin["partial_binding_ratio"][i] = am_interaction[i][j].partial_binding_ratio;
                    }
                    if (actin["myosin_binding_ratio"][i] < am_interaction[i][j].myosin_binding_ratio) {
                        actin["myosin_binding_ratio"][i] = am_interaction[i][j].myosin_binding_ratio;
                    }
                }
            }
        }
    }
    actin["crosslink_ratio"][i] = actin_crosslink_ratio[i];  
}

void Sarcomere::_process_catch_bonds(int& i) {
    int thread_id = omp_get_thread_num();
    std::vector <int> actin_neighbors = actin_neighbors_by_species[i].first;
    std::vector<double> cb_strengths;
    std::vector<int> cb_indices;
    for (int index = 0; index < actin_neighbors.size(); index++){
            int j = actin_neighbors[index];
            if (i>j){
                double cos_angle = actin.direction[i].dot(actin.direction[j]);
                double strength = _get_cb_strength(i,j);
                if (strength>EPS){
                    cb_indices.push_back(j);
                    cb_strengths.push_back(strength);
                }
            }
    }
    _set_cb(i, cb_indices,cb_strengths);
}

void Sarcomere::_calc_am_force_velocity(int& i) {
    int thread_id = omp_get_thread_num();
    // Thread-local temporary lists for forces and velocities
    auto& local_actin_forces = actin_forces_temp[thread_id];
    auto& local_actin_torques = actin_torques_temp[thread_id];
    auto& local_myosin_forces = myosin_forces_temp[thread_id];
    auto& local_myosin_velocities = myosin_velocities_temp[thread_id];
    auto& local_myosin_f_load = myosin_f_load_temp[thread_id];
    auto& local_myosin_torques = myosin_torques_temp[thread_id];
    std::vector<int> myosin_indices = myosinIndicesPerActin.getConnections(i);
    vec velocity = v_am * actin.direction[i];
    double f_load;
    double cutoff = myosin.radius / myosin_radius_ratio;
    for (int index = 0; index < myosin_indices.size(); index++) {
        int j = myosin_indices[index];
        double partial_binding_ratio = am_interaction[i][j].partial_binding_ratio;
        double k_am_adjusted = k_am * actin.cb_strength[i] * actin.f_load[i] * (partial_binding_ratio>EPS);
        double kappa_am_adjusted = kappa_am*std::min(partial_binding_ratio*3,1.);
        vector force_vec = compute_am_force_and_energy(
            actin, myosin, i, j, box, k_am_adjusted, kappa_am_adjusted, cutoff);
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

        if (actin.cb_strength[i]>1/cb_mult_factor && partial_binding_ratio>EPS){
            double binding_ratio_adjusted = std::min(partial_binding_ratio,1.0);
            double abs_cos_angle = std::abs(actin.direction[i].dot(myosin.direction[j]));
            binding_ratio_adjusted = binding_ratio_adjusted * abs_cos_angle;
            f_load = (1-std::exp(-2*binding_ratio_adjusted))/(1-std::exp(-2));
        }
        else{
            f_load = 0;
        }
        if (f_load<local_myosin_f_load[j]){
            local_myosin_f_load[j] = f_load;
        }
        velocity = velocity * (1-f_load);
        if (actin.cb_strength[i]>1/cb_mult_factor){
            actin.velocity[i] = velocity*diff_coeff_ratio/(diff_coeff_ratio+1);
            local_myosin_velocities[j] -= velocity/(diff_coeff_ratio+1);
        }
        else{
            actin.velocity[i] = velocity;
        }

    }
}



void Sarcomere::_volume_exclusion(){
    #pragma omp for schedule(dynamic)
    for (int i = 0; i<myosin.n; i++){
        auto result = neighbor_list.get_neighbors_by_type(i+actin.n);
        std::vector <int> myosin_indices = result.second;
        for (int index = 0; index < myosin_indices.size(); index++){
            int j = myosin_indices[index];
            if (i<j){
            _myosin_repulsion(i,myosin_indices[index]);}
        }
    }

    #pragma omp for schedule(dynamic)
    for (int i = 0; i < actin.n; i++){
        auto result = neighbor_list.get_neighbors_by_type(i);
        std::vector <int> actin_indices = result.first;
        for (int index = 0; index < actin_indices.size(); index++){
            int j = actin_indices[index];
            if (i<j){
            _actin_repulsion(i,actin_indices[index]);}
        }
    }
}

void Sarcomere::_myosin_exclusion(){
    #pragma omp for schedule(dynamic)
    for (int i = 0; i<myosin.n; i++){
        auto result = neighbor_list.get_neighbors_by_type(i+actin.n);
        std::vector <int> myosin_indices = result.second;
        for (int index = 0; index < myosin_indices.size(); index++){
            int j = myosin_indices[index];
            if (i<j){
            _myosin_repulsion(i,myosin_indices[index]);}
        }
    }
}

void Sarcomere::_myosin_repulsion(int& i, int& j){
    int thread_id = omp_get_thread_num();
    vec center_displacement = myosin.center[i] - myosin.center[j];
    center_displacement.pbc_wrap(box);
    double center_distance = center_displacement.norm();
    double cutoff = 2*myosin.radius;
    if (center_distance<=cutoff+myosin.length){
        auto& local_myosin_forces = myosin_forces_temp[thread_id]; 
        auto result = geometry::segment_segment_distance_w_normal(myosin.left_end[i], 
            myosin.right_end[i], myosin.left_end[j], myosin.right_end[j], box);
        double distance = result.first;
        if (distance<cutoff){
            vec normal_vector = result.second["normal"];
            double norm = normal_vector.norm();
            double factor = std::min((50*(cutoff - distance)/cutoff), 50.0);
            if (norm==0){
                normal_vector.x = center_displacement.x;
                normal_vector.y = center_displacement.y;
                normal_vector = normal_vector/center_distance;
            }
            else {
                normal_vector = normal_vector/norm;
            }
            // vec normal_vector = center_displacement/center_distance;
            if (i<fix_myosin){
                local_myosin_forces[j]-=2*factor*normal_vector;
                return;
            }
            else if (j<fix_myosin){
                local_myosin_forces[i]+=2*factor*normal_vector;
                return;
            }
            //check if any of the myosin is bound to cb-forming actin
            auto actin_indices_i = actinIndicesPerMyosin.getConnections(i);
            double cb_strength_i = 0;
            for (int index = 0; index < actin_indices_i.size(); index++){
                int k = actin_indices_i[index];
                cb_strength_i += actin.cb_strength[k];
            }
            auto actin_indices_j = actinIndicesPerMyosin.getConnections(j);
            double cb_strength_j = 0;
            for (int index = 0; index < actin_indices_j.size(); index++){
                int k = actin_indices_j[index];
                cb_strength_j += actin.cb_strength[k];
            }
            if (cb_strength_i>0.3 && cb_strength_j<0.3){
                //move myosin j only
                local_myosin_forces[j]= local_myosin_forces[j]- 2*factor*normal_vector;
            }
            else if (cb_strength_i<0.3 && cb_strength_j>0.3){
                //move myosin i only
                local_myosin_forces[i]+=2*factor*normal_vector;
            }
            else{
                local_myosin_forces[i]+=factor*normal_vector;
                local_myosin_forces[j]-=factor*normal_vector;
            }
        }
    }
}

void Sarcomere::_actin_repulsion(int& i, int& j){
    int thread_id = omp_get_thread_num();
    vec center_displacement = actin.center[i] - actin.center[j];
    center_displacement.pbc_wrap(box);
    double center_distance = center_displacement.norm();
    if (center_distance<=crosslinker_length+actin.length){
        double cos_angle = actin.direction[i].dot(actin.direction[j]);
        // Only apply repulsion if the angle difference is small enough
        auto& local_actin_forces = actin_forces_temp[thread_id];
        auto result = geometry::segment_segment_distance_w_normal(actin.left_end[i], 
            actin.right_end[i], actin.left_end[j], actin.right_end[j], box);
        double distance = result.first;
        double min_dist = 0.01;
        if (distance < min_dist){
            vec normal_vector = result.second["normal"];
            double norm = normal_vector.norm();
            double factor = 10*(min_dist - distance)/min_dist;
            if (norm==0){
                normal_vector.x = center_displacement.x;
                normal_vector.y = center_displacement.y;
                normal_vector = normal_vector/center_distance;
            }
            else {
                normal_vector = normal_vector/norm;
            }
            if (actin.cb_strength[i]>0.3 && actin.cb_strength[j]<0.3){
                //move actin j only
                local_actin_forces[j] = local_actin_forces[j]-2*factor*normal_vector;
            }
            else if (actin.cb_strength[i]<0.3 && actin.cb_strength[j]>0.3){
                //move actin i only
                local_actin_forces[i] = local_actin_forces[i]+2*factor*normal_vector;
            }
            else{
                local_actin_forces[i] = local_actin_forces[i]+factor*normal_vector;
                local_actin_forces[j] = local_actin_forces[j]-factor*normal_vector;
            }
        }
    }
}

double Sarcomere::_get_cb_strength(int& i, int& j){
    bool crosslink = false;
    if (actin_crosslink_ratio[i] > EPS && actin_crosslink_ratio[j] > EPS || ! directional){
        double distance = geometry::segment_segment_distance(actin.left_end[i], 
            actin.right_end[i], actin.left_end[j], actin.right_end[j], box);
        if (distance<crosslinker_length){
            crosslink = true;
        }
    }
    if (crosslink){ 
        double strength = 0;
        double cos_angle = actin.direction[i].dot(actin.direction[j]);
        double abs_cos_angle = std::abs(cos_angle);
        bool catch_bond =(actin_basic_tension[i]>EPS && actin_basic_tension[j]>EPS);
        if (directional){
            catch_bond = (catch_bond && cos_angle<0);
        }
        if (catch_bond){
            auto& myosin_indices_i = myosinIndicesPerActin.getConnections(i);
            auto& myosin_indices_j = myosinIndicesPerActin.getConnections(j);
            if (!utils::compare_indices(myosin_indices_i,myosin_indices_j)){
                strength = abs_cos_angle;
                if (actin.f_load[i] == 0){
                    _get_f_load(i);
                }
                if (actin.f_load[j] == 0){
                    _get_f_load(j);
                }
            }
        }
        else{
            strength = 1/cb_mult_factor*abs_cos_angle*abs_cos_angle;
        }
        return strength;
    }
    return 0.0;
}

void Sarcomere::_get_f_load(int& i){
    // i is the actin index
    actin.f_load[i] = 0;
    std::vector<int> myosin_indices = myosinIndicesPerActin.getConnections(i);
    double binding_ratio_adjusted_sum = 0;
    for (int index = 0; index < myosin_indices.size(); index++){
        int j = myosin_indices[index];
        double binding_ratio_adjusted;
        if (directional){
            binding_ratio_adjusted = std::min(am_interaction[i][j].partial_binding_ratio*3,1.0);
        }
        else{
            binding_ratio_adjusted = std::min(am_interaction[i][j].myosin_binding_ratio*3,1.0);
        }
        double abs_cos_angle = std::abs(actin.direction[i].dot(myosin.direction[j]));
        binding_ratio_adjusted_sum += binding_ratio_adjusted * abs_cos_angle * abs_cos_angle;
        if (binding_ratio_adjusted_sum>1){
            binding_ratio_adjusted_sum = 1;
            break;
        }
    }
    actin.f_load[i] = (1-std::exp(-2*binding_ratio_adjusted_sum))/(1-std::exp(-2));
}

void Sarcomere::_set_cb(int& i, int& j, double& normalized_strength, bool& add_connection){
    int thread_id = omp_get_thread_num();
    if (actin_actin_bonds[i][j] == 1) {
        std::cout << "Skipping duplicate bond: " << i << ", " << j << "\n";
    }
    auto& local_actin_forces = actin_forces_temp[thread_id];
    auto& local_actin_torques = actin_torques_temp[thread_id];
    auto& local_actin_cb_strengths = actin_cb_strengths_temp[thread_id];
    auto& local_myosin_forces = myosin_forces_temp[thread_id];
    auto& local_myosin_torques = myosin_torques_temp[thread_id];
    if (normalized_strength<EPS){
        return;
    }
    double rand = gsl_rng_uniform(rng_engines[thread_id]);
    double abs_cos_angle = std::abs(actin.direction[i].dot(actin.direction[j]));
    double f_load;
    if (normalized_strength<=1/cb_mult_factor){
        f_load = 0;
    }
    else {f_load = abs_cos_angle * std::min(actin.f_load[i],actin.f_load[j]);}
    if (actin_actin_bonds_prev[i][j] == 1) {
        double k_off_adjusted = dt /(base_lifetime+lifetime_coeff*f_load);
        if (rand < k_off_adjusted){ //k_off is actually k_off * dt
            //printf("breaking bond between %d and %d\n",i,j);
            return;
        }
    }
    else{
        if (rand >= k_on*dt){ //k_on is actually k_on * dt
            // printf("not forming bond between %d and %d\n",i,j);
            return;
        }
       // printf("forming bond between %d and %d\n",i,j);
    }
    vector force_vec = compute_aa_force_and_energy(actin,i, j, box,k_aa,kappa_aa); 
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
    local_actin_cb_strengths[i] += normalized_strength;
    local_actin_cb_strengths[j] += normalized_strength;
    actin_n_bonds[i] += 1;
    actin_n_bonds[j] += 1;
    actin_actin_bonds[i][j] = 1;
    actin_actin_bonds[j][i] = 1;
    // Update lifetime: increment if bond persisted, reset if new
    actin_actin_lifetime[i][j] = actin_actin_lifetime_prev[i][j] + 1;
    actin_actin_lifetime[j][i] = actin_actin_lifetime[i][j];
}

void Sarcomere::_set_cb(int& i, std::vector<int> indices, vector cb_strength){
    int thread_id = omp_get_thread_num();
    std::vector<size_t> sorted_indices = utils::sort_indices(cb_strength);
    //print all the sorted indices
    bool add_connection;
    for (int index = 0; index < sorted_indices.size(); index++){
        int j = indices[sorted_indices[index]];
        add_connection = true;
        _set_cb(i,j,cb_strength[sorted_indices[index]],add_connection);
    }
}

void Sarcomere::debug_cb_stats(){
    int close_opposite = 0;
    int diff_myosin = 0;
    for(int i=0;i<actin.n;i++){
        for(int j=i+1;j<actin.n;j++){
            double distance = geometry::segment_segment_distance(actin.left_end[i],
                actin.right_end[i], actin.left_end[j], actin.right_end[j], box);
            if(distance < crosslinker_length){
                double cos_angle = actin.direction[i].dot(actin.direction[j]);
                if(cos_angle < 0){
                    close_opposite++;
                    auto mi = myosinIndicesPerActin.getConnections(i);
                    auto mj = myosinIndicesPerActin.getConnections(j);
                    if(!mi.empty() && !mj.empty()){
                        bool shared = false;
                        for(int m: mi){
                            if(std::find(mj.begin(), mj.end(), m) != mj.end()){
                                shared = true;
                                break;
                            }
                        }
                        if(!shared){
                            diff_myosin++;
                            std::cout << "Pair " << i << "-" << j
                                      << " cb_strengths: " << actin.cb_strength[i]
                                      << ", " << actin.cb_strength[j] << std::endl;
                        }
                    }
                }
            }
        }
    }
    int active = 0;
    double total_life = 0.0;
    double max_life = 0.0;
    for(int i=0;i<actin.n;i++){
        for(int j=i+1;j<actin.n;j++){
            if(actin_actin_bonds[i][j]==1){
                active++;
                double lt = actin_actin_lifetime[i][j]*dt;
                total_life += lt;
                if(lt>max_life) max_life = lt;
            }
        }
    }
    double avg_life = active>0 ? total_life/active : 0.0;
    std::cout << "CB Debug: close_opposite_pairs=" << close_opposite
              << ", pairs_diff_myosin=" << diff_myosin
              << ", avg_lifetime=" << avg_life
              << ", max_lifetime=" << max_life << std::endl;
}


std::pair<std::vector<double>, std::vector<double>> Sarcomere::_extract_bonded_pairs(
    const std::vector<std::vector<int>>& actin_actin_bonds,
    const utils::MoleculeConnection& myosinIndicesPerActin)
{
    // First, flatten the actin-actin bonds matrix (upper-triangle only)
    // into a vector of doubled indices.
    std::vector<double> flatActinBonds;
    for (int i = 0; i < actin.n; ++i) {
        for (int j = i + 1; j < actin.n; ++j) {
            if (actin_actin_bonds[i][j] == 1) {
                flatActinBonds.push_back(static_cast<double>(i));
                flatActinBonds.push_back(static_cast<double>(j));
            }
        }
    }

    // Next, process the flattened actin bonds to extract unique myosin–myosin pairs.
    // For each actin–actin bond (i, j), retrieve the first myosin connected to actin i and j.
    std::set<std::pair<int,int>> uniqueMyosinBonds;
    for (size_t k = 0; k < flatActinBonds.size(); k += 2) {
        int i = static_cast<int>(flatActinBonds[k]);
        int j = static_cast<int>(flatActinBonds[k+1]);

        auto connections_i = myosinIndicesPerActin.getConnections(i);
        auto connections_j = myosinIndicesPerActin.getConnections(j);

        // Skip if either actin has no attached myosin.
        if (connections_i.empty() || connections_j.empty())
            continue;

        int myosin_i = connections_i[0]; // Use the first connection.
        int myosin_j = connections_j[0]; // Use the first connection.

        // Skip if both actins attach to the same myosin.
        if (myosin_i == myosin_j)
            continue;

        // Order the pair so that the smaller index comes first.
        if (myosin_i > myosin_j)
            std::swap(myosin_i, myosin_j);

        uniqueMyosinBonds.insert(std::make_pair(myosin_i, myosin_j));
    }

    // Flatten the set of unique myosin bonds into a vector<double>.
    std::vector<double> flattenedMyosinBonds;
    for (const auto& bond : uniqueMyosinBonds) {
        flattenedMyosinBonds.push_back(static_cast<double>(bond.first));
        flattenedMyosinBonds.push_back(static_cast<double>(bond.second));
    }

    // Return the pair: first element is actin bonds, second is myosin bonds.
    return {flatActinBonds, flattenedMyosinBonds };
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

            double acc_cb = 0.1; // accumulated cross-bridge strength
            for (int a : idxs) {
                double cb = actin.cb_strength[a];
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
            if (actin.cb_strength[i] < 0.1) continue; // skip if no active cross-bridge
            vec u = actin.direction[i];
            u.normalize();
            vec tau = _alignment_torque(u, k_theta_bias * actin.cb_strength[i]);
            actin_torques_temp[omp_get_thread_num()][i] += tau;
        }
}


void Sarcomere::new_file(){
    create_file(filename, actin, myosin);
}

void Sarcomere::save_state(){
    std::pair<std::vector<double>, std::vector<double>> bondPairData =
    _extract_bonded_pairs(actin_actin_bonds, myosinIndicesPerActin);
    std::vector<double> flatActinBonds = bondPairData.first;
    std::vector<double> flatMyosinBonds = bondPairData.second;
    append_to_file(filename, actin, myosin, flatActinBonds, flatMyosinBonds);
}

void Sarcomere::load_state(int& n_frames){
    load_from_file(filename, actin, myosin, actin_actin_bonds, n_frames);
    update_system();
}
