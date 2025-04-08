#include "sarcomere.h"

// Constructor
Sarcomere::Sarcomere() {}

// Parameterized Constructor
Sarcomere::Sarcomere(int& n_actins, int& n_myosins, vector box0, double& actin_length, double& myosin_length,
              double& myosin_radius, double& crosslinker_length, double& k_on, double& k_off,
              double& base_lifetime, double& lifetime_coeff, double& diff_coeff_ratio, double& k_aa, double& kappa_aa, double& k_am, double& kappa_am, double& v_am,
              std::string& filename, gsl_rng* rng, int& seed, int& fix_myosin, double& dt, bool& directional)
            : actin(n_actins, actin_length, box0, rng),
              myosin(n_myosins, myosin_length, myosin_radius, box0, rng),
              myosinIndicesPerActin(n_actins),
              actinIndicesPerMyosin(n_myosins),
              neighbor_list(0.0, box0, 0.0),
              actin_actin_bonds(n_actins, std::vector<int>(n_actins, 0)),
            actin_forces_temp(omp_get_max_threads(), std::vector<vec>(n_actins, {0, 0, 0})),
            myosin_forces_temp(omp_get_max_threads(), std::vector<vec>(n_myosins, {0, 0, 0})),
            myosin_velocities_temp(omp_get_max_threads(), std::vector<vec>(n_myosins, {0, 0, 0})),
            actin_angular_forces_temp(omp_get_max_threads(), std::vector<std::vector<double>>(n_actins, {0,0})),
            myosin_angular_forces_temp(omp_get_max_threads(), std::vector<std::vector<double>>(n_myosins, {0,0})),
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
            this->base_lifetime = base_lifetime;
            this->lifetime_coeff = lifetime_coeff;
            this->diff_coeff_ratio = diff_coeff_ratio;
            this->directional = directional;
            cb_mult_factor = 1000;
            cutoff_radius = std::max(actin_length, myosin_length) +
                            std::max(2 * myosin_radius, crosslinker_length);
            double skin_distance = 0.1 * cutoff_radius;
            neighbor_list = NeighborList(cutoff_radius + skin_distance, box, skin_distance / 2);
            neighbor_list.initialize(actin.center, myosin.center);
            actin_actin_bonds_prev = actin_actin_bonds;
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
            // for (int i = 0; i < n_myosins; i++) {
            //     myosin.theta[i] = 0;
            // }
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
    for (int i = 0; i < myosin.n; i++){
        myosin.theta[i] = 0;
        myosin.phi[i] = M_PI/2;
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
    actin.theta[0] = 0;
    actin.phi[0] = M_PI/2;
    actin.theta[1] = M_PI;
    actin.phi[1] = M_PI/2;
    
    std::vector<vector> myosin_positions = {
        {-1.33, 0, 0}, {1.33, 0, 0}
    };
    for (int i = 0; i < myosin_positions.size(); i++){
        myosin.center[i].x = myosin_positions[i][0];
        myosin.center[i].y = myosin_positions[i][1];
        myosin.center[i].z = myosin_positions[i][2]; // set z coordinate to 0
        myosin.theta[i] = 0;
        myosin.phi[i] = M_PI/2;
    }
}

void Sarcomere::sarcomeric_structure(){
    std::vector<vector> myosin_positions;
    box[0] = 5.32;
    box[1] = 5.32;
    box[2] = 5.32; // set box's third dimension
    myosin_positions = {
        {-1.33, -1.0, 0}, {1.33, -1.0, 0}, {-1.33, 0.0, 0}, {1.33, 0.0, 0}, {-1.33, 1.0, 0}, {1.33, 1.0, 0}
    };
    for (int i = 0; i < myosin_positions.size(); i++){
        myosin.center[i].x = myosin_positions[i][0];
        myosin.center[i].y = myosin_positions[i][1];
        myosin.center[i].z = myosin_positions[i][2]; // set z coordinate to 0
        myosin.theta[i] = 0;
        myosin.phi[i] = M_PI/2;
    }

    myosin.update_endpoints();
    
    std::vector<vector> actin_positions = {
        {0.5, -1.15, 0}, {0.5, -1.0, 0}, {0.5, -0.85, 0}, {0.5, -0.15, 0}, {0.5, 0.0, 0}, {0.5, 0.15, 0},
        {0.5, 0.85, 0}, {0.5, 1.0, 0}, {0.5, 1.15, 0}, 
        {-2.16, -1.15, 0}, {-2.16, -1.0, 0}, {-2.16, -0.85, 0},
        {-2.16, -0.15, 0}, {-2.16, 0.0, 0}, {-2.16, 0.15, 0},
        {-2.16, 0.85, 0}, {-2.16, 1.0, 0}, {-2.16, 1.15, 0}
    };
    for (int i = 0; i < actin_positions.size(); i++){
        actin.center[i].x = actin_positions[i][0];
        actin.center[i].y = actin_positions[i][1];
        actin.center[i].z = actin_positions[i][2]; // set z coordinate to 0
        actin.theta[i] = 0;
        actin.phi[i] = 0;
    }
    int n = actin_positions.size();
    actin_positions = {
        {-0.5, -1.15, 0}, {-0.5, -1.0, 0}, {-0.5, -0.85, 0}, {-0.5, -0.15, 0}, {-0.5, 0.0, 0}, {-0.5, 0.15, 0},
        {-0.5, 0.85, 0}, {-0.5, 1.0, 0}, {-0.5, 1.15, 0}, 
        {2.16, -1.15, 0}, {2.16, -1.0, 0}, {2.16, -0.85, 0},
        {2.16, -0.15, 0}, {2.16, 0.0, 0}, {2.16, 0.15, 0},
        {2.16, 0.85, 0}, {2.16, 1.0, 0}, {2.16, 1.15, 0}
    };
    for (int i = 0; i < actin_positions.size(); i++){
        actin.center[i+n].x = actin_positions[i][0];
        actin.center[i+n].y = actin_positions[i][1];
        actin.center[i+n].z = actin_positions[i][2]; // set z coordinate to 0
        actin.theta[i+n] = M_PI;
        actin.phi[i+n] = 0;
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
        // Step 2: Compute actin-myosin binding
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
            _actin_myosin_force(i);
        }

        _volume_exclusion();

        #pragma omp barrier  

        // Step 7: Reduce actin forces and angular forces
        utils::reduce_array(actin_forces_temp, actin.force);
        utils::reduce_array(actin_angular_forces_temp, actin.angular_force);

        // Step 8: Reduce myosin forces, velocities, and angular forces
        utils::reduce_array(myosin_forces_temp, myosin.force);
        utils::reduce_array(myosin_velocities_temp, myosin.velocity);
        utils::reduce_array(myosin_angular_forces_temp, myosin.angular_force);
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
        utils::reduce_array(actin_forces_temp, actin.force);
        utils::reduce_array(actin_angular_forces_temp, actin.angular_force);
        // Step 8: Reduce myosin forces, velocities, and angular forces
        utils::reduce_array(myosin_forces_temp, myosin.force);
        utils::reduce_array(myosin_angular_forces_temp, myosin.angular_force);
    }
}


void Sarcomere::_update_neighbors() {
    neighbor_list.set_species_positions(actin.center, myosin.center);
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
            actin_angular_forces_temp[t][i] = {0, 0};
            actin_cb_strengths_temp[t][i] = 0;
        }
    }

    #pragma omp for
    for (int t = 0; t < omp_get_max_threads(); ++t) {
        for (int i = 0; i < myosin.n; ++i) {
            myosin_forces_temp[t][i] = {0, 0, 0};
            myosin_velocities_temp[t][i] = {0, 0, 0};
            myosin_angular_forces_temp[t][i] = {0, 0};
            actinIndicesPerMyosin_temp[t].deleteAllConnections(i);
        }
    }

    #pragma omp for
    for (int i = 0; i < actin.n; i++) {
        actin.update_endpoints(i);
        actin.force[i] = {0, 0, 0};
        actin.angular_force[i] = {0, 0};
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
        }
    }

    #pragma omp for
    for (int i = 0; i < myosin.n; i++) {
        myosin.update_endpoints(i);
        myosin.force[i] = {0, 0, 0};
        myosin.angular_force[i] = {0, 0};
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
    double max_binding_ratio = -1;
    double myosin_binding_ratio_sum = 0;
    int myosin_index = -1;
    auto myosin_neighbors = actin_neighbors_by_species[i].second;
    for (int index = 0; index < myosin_neighbors.size(); index++) {
        int j = myosin_neighbors[index];
        am_interaction[i][j] = geometry::analyze_am(
            actin.left_end[i], actin.right_end[i], myosin.left_end[j], myosin.right_end[j],
            myosin.radius, box);
        if (am_interaction[i][j].myosin_binding_ratio>0){
            n_myosins_per_actin[i]++;
            myosin_binding_ratio_sum += am_interaction[i][j].myosin_binding_ratio;
            if (directional){
                if (am_interaction[i][j].partial_binding_ratio > max_binding_ratio) {
                    max_binding_ratio = am_interaction[i][j].partial_binding_ratio;
                    myosin_index = j;
                }
            }
            else{
                if (am_interaction[i][j].myosin_binding_ratio > max_binding_ratio) {
                    max_binding_ratio = am_interaction[i][j].myosin_binding_ratio;
                    myosin_index = j;
                }
            }
        }
    }
    if (myosin_index >= 0) {
        local_actinIndicesPerMyosin.addConnection(myosin_index, i);
        myosinIndicesPerActin.addConnection(i, myosin_index);
        double myosin_binding_ratio = am_interaction[i][myosin_index].myosin_binding_ratio;
        actin_crosslink_ratio[i] = am_interaction[i][myosin_index].crosslinkable_ratio-(myosin_binding_ratio_sum-myosin_binding_ratio);
        actin["myosin_binding_ratio"][i] = myosin_binding_ratio;
        actin["crosslink_ratio"][i] = actin_crosslink_ratio[i];
        actin["partial_binding_ratio"][i] = am_interaction[i][myosin_index].partial_binding_ratio;
        double a_m_angle = actin.theta[i] - myosin.theta[myosin_index];
        double abs_cos = std::abs(std::cos(a_m_angle));
        actin_basic_tension[i] = abs_cos;
    }
}

void Sarcomere::_process_catch_bonds(int& i) {
    int thread_id = omp_get_thread_num();
    std::vector <int> actin_neighbors = actin_neighbors_by_species[i].first;
    std::vector<double> cb_strengths;
    std::vector<int> cb_indices;
    for (int index = 0; index < actin_neighbors.size(); index++){
            int j = actin_neighbors[index];
            if (i>j){
                double cos_angle = std::cos(actin.theta[i] - actin.theta[j]);
                double strength = _get_cb_strength(i,j);
                if (strength>EPS){
                    cb_indices.push_back(j);
                    cb_strengths.push_back(strength);
                }
            }
    }
    _set_cb(i, cb_indices,cb_strengths);
}

void Sarcomere::_actin_myosin_force(int& i) {
    int thread_id = omp_get_thread_num();
    // Thread-local temporary lists for forces and velocities
    auto& local_actin_forces = actin_forces_temp[thread_id];
    auto& local_actin_angular_forces = actin_angular_forces_temp[thread_id];
    auto& local_myosin_forces = myosin_forces_temp[thread_id];
    auto& local_myosin_velocities = myosin_velocities_temp[thread_id];
    auto& local_myosin_angular_forces = myosin_angular_forces_temp[thread_id];
    std::vector<int> myosin_indices = myosinIndicesPerActin.getConnections(i);
    vec velocity = {v_am*cos(actin.theta[i]),v_am*sin(actin.theta[i])};
    for (int index = 0; index < myosin_indices.size(); index++) {
        int j = myosin_indices[index];
        double k_am_adjusted = k_am * actin.cb_strength[i] * actin.f_load[i] * (actin["partial_binding_ratio"][i]>EPS);
        vector force_vec = compute_am_force_and_energy(
            actin, myosin, i, j, box, k_am_adjusted, kappa_am, myosin.radius);
        local_actin_forces[i].x += force_vec[0];
        local_actin_forces[i].y += force_vec[1];
        local_actin_forces[i].z += force_vec[2];
        local_myosin_forces[j].x -= force_vec[0];
        local_myosin_forces[j].y -= force_vec[1];
        local_myosin_forces[j].z -= force_vec[2];
        local_actin_angular_forces[i][0] += force_vec[3];
        local_actin_angular_forces[i][1] += force_vec[4];
        local_myosin_angular_forces[j][0] += force_vec[5];
        local_myosin_angular_forces[j][1] += force_vec[6];
        velocity = velocity * (1-actin.f_load[i]*(actin.cb_strength[i]>1/cb_mult_factor)*(actin["partial_binding_ratio"][i]>EPS));
        actin.velocity[i] = velocity*diff_coeff_ratio/(diff_coeff_ratio+1);
        local_myosin_velocities[j] -= velocity/(diff_coeff_ratio+1);
    }
}

void Sarcomere::_volume_exclusion(){
    #pragma omp for schedule(dynamic)
    for (int i = 0; i<myosin.n; i++){
        auto result = neighbor_list.get_neighbors_by_type(i+actin.n);
        std::vector <int> myosin_indices = result.second;
        for (int index = 0; index < myosin_indices.size(); index++){
            //printf("i: %d, j: %d, thread_id: %d \n", i, myosin_indices[index], thread_id);
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
            //printf("i: %d, j: %d, thread_id: %d \n", i, myosin_indices[index], thread_id);
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
            vec normal_vector = result.second["vector"];
            double norm = normal_vector.norm();
            double factor = std::min((300*(cutoff - distance)/cutoff),100.);
            if (norm==0){
                normal_vector.x = center_displacement.x;
                normal_vector.y = center_displacement.y;
                normal_vector = normal_vector/center_distance;
            }
            else {
                normal_vector = normal_vector/norm;
            }
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
                local_myosin_forces[j]-=2*factor*normal_vector;
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
        auto& local_actin_forces = actin_forces_temp[thread_id];
        auto result = geometry::segment_segment_distance_w_normal(actin.left_end[i], 
            actin.right_end[i], actin.left_end[j], actin.right_end[j], box);
        double distance = result.first;
        double min_dist = 0.01;
        if (distance < min_dist){
            vec normal_vector = result.second["vector"];
            double norm = normal_vector.norm();
            double factor = std::min(20*(min_dist - distance)/min_dist,10.);
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
            // printf("actin %d and actin %d: (%f, %f), (%f, %f), distance: %f \n", i, j, actin.center[i].x, actin.center[i].y, actin.center[j].x, actin.center[j].y, distance);
            // printf("force for actin %d: %f, %f \n", i, local_actin_forces[i].x, local_actin_forces[i].y);
            // printf("force for actin %d: %f, %f \n", j, local_actin_forces[j].x, local_actin_forces[j].y);
        }

    }
}

double Sarcomere::_get_cb_strength(int& i, int& j){
    double angle = actin.theta[i] - actin.theta[j];
    double cos_angle = std::cos(angle);
    bool crosslink = false;
    if (actin_crosslink_ratio[i]>EPS && actin_crosslink_ratio[j]>EPS || ! directional){
        double distance = geometry::segment_segment_distance(actin.left_end[i], 
            actin.right_end[i], actin.left_end[j], actin.right_end[j], box);
        if (distance<crosslinker_length){
            crosslink = true;
        }
    }
    if (crosslink){ 
        double strength = 0;
        double abs_cos_angle = std::abs(cos_angle);
        bool catch_bond =(n_myosins_per_actin[i]<2 && n_myosins_per_actin[j]<2 &&
            actin_basic_tension[i]>EPS && actin_basic_tension[j]>EPS);
        if (directional){
            catch_bond = (catch_bond && cos_angle<0);
        }
        if (catch_bond){
            int myosin_index_i = myosinIndicesPerActin.getConnections(i)[0];
            int myosin_index_j = myosinIndicesPerActin.getConnections(j)[0];
            if (myosin_index_i != myosin_index_j){
                strength = abs_cos_angle;
                _get_f_load(i,myosin_index_i);
                _get_f_load(j,myosin_index_j);
            }
        }
        else{
            strength = 1/cb_mult_factor*abs_cos_angle;
        }
        return strength;
    }
    return 0.0;
}

void Sarcomere::_get_f_load(int& i, int& j){
    // i is the actin index and j is the myosin index
    double binding_ratio_adjusted = std::min(am_interaction[i][j].myosin_binding_ratio*3,1.0);
    binding_ratio_adjusted = 1-std::exp(-4*binding_ratio_adjusted);
    double abs_cos_angle = std::abs(std::cos(actin.theta[i] - myosin.theta[j]));
    actin.f_load[i] = binding_ratio_adjusted * abs_cos_angle;
}

void Sarcomere::_set_cb(int& i, int& j, double& normalized_strength, bool& add_connection){
    int thread_id = omp_get_thread_num();
    //printf("thread_id: %d, i: %d, j: %d, strength: %f \n", thread_id, i, j, normalized_strength);
    if (actin_actin_bonds[i][j] == 1) {
        std::cout << "Skipping duplicate bond: " << i << ", " << j << "\n";
    }
    auto& local_actin_forces = actin_forces_temp[thread_id];
    auto& local_actin_angular_forces = actin_angular_forces_temp[thread_id];
    auto& local_actin_cb_strengths = actin_cb_strengths_temp[thread_id];
    auto& local_myosin_forces = myosin_forces_temp[thread_id];
    auto& local_myosin_angular_forces = myosin_angular_forces_temp[thread_id];
    if (normalized_strength<=1/cb_mult_factor && (actin_n_bonds[i]>=2 || actin_n_bonds[j]>=2)){
        return;
    }
    if (normalized_strength<EPS){
        return;
    }
    double rand = gsl_rng_uniform(rng_engines[thread_id]);
    double abs_cos_angle = std::abs(std::cos(actin.theta[i] - actin.theta[j]));
    double f_load = std::min(actin.f_load[i],actin.f_load[j]);
    if (actin_actin_bonds_prev[i][j] == 1) {
        double k_off_adjusted = dt * abs_cos_angle/(base_lifetime+lifetime_coeff*f_load);
        if (rand < k_off_adjusted){ //k_off is actually k_off * dt
            printf("breaking bond between %d and %d\n",i,j);
            return;
        }
    }
    else{
        if (rand >= k_on*dt){ //k_on is actually k_on * dt
            // printf("not forming bond between %d and %d\n",i,j);
            return;
        }
        printf("forming bond between %d and %d\n",i,j);
    }
    vector force_vec = compute_aa_force_and_energy(actin,i, j, box,k_aa,kappa_aa); 
    local_actin_forces[i].x += force_vec[0];
    local_actin_forces[i].y += force_vec[1];
    local_actin_forces[i].z += force_vec[2];
    local_actin_forces[j].x -= force_vec[0];
    local_actin_forces[j].y -= force_vec[1];
    local_actin_forces[j].z -= force_vec[2];
    local_actin_angular_forces[i][0] += force_vec[3];
    local_actin_angular_forces[i][1] += force_vec[4];
    local_actin_angular_forces[j][0] += force_vec[5];
    local_actin_angular_forces[j][1] += force_vec[6];
    local_actin_cb_strengths[i] += normalized_strength;
    local_actin_cb_strengths[j] += normalized_strength;
    actin_n_bonds[i] += 1;
    actin_n_bonds[j] += 1;
    actin_actin_bonds[i][j] = 1;
    actin_actin_bonds[j][i] = 1;
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

void Sarcomere::new_file(){
    create_file(filename, actin, myosin);
}

void Sarcomere::save_state(){
    append_to_file(filename, actin, myosin);
}

void Sarcomere::load_state(){
    load_from_file(filename, actin, myosin);
    update_system();
}

