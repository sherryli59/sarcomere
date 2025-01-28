#ifndef SARCOMERE_H
#define SARCOMERE_H

#include <iostream>

#ifndef COMPONENTS_H
#include "components.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

#ifndef INTERACTION_H
#include "interaction.h"
#endif

#ifndef GEOMETRY_H
#include "geometry.h"
#endif

#ifndef H5_UTILS_H
#include "h5_utils.h"
#endif

#ifndef NEIGHBORLIST_H
#include "neighborlist.h"
#endif

#include <cmath>
#include <vector>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <omp.h>

using vector = std::vector<double>;
using interaction = geometry::am_interaction;
const double EPS = 1e-6;

class Sarcomere
{
    public:
        Filament actin;
        Myosin myosin;
        NeighborList neighbor_list; 
        utils::MoleculeConnection myosinIndicesPerActin;
        utils::MoleculeConnection actinIndicesPerMyosin;
        utils::MoleculeConnection actinIndicesPerActin;
        std::vector<std::vector<int>> actin_actin_bonds, actin_actin_bonds_prev;
        vector box;
        double k_aa, kappa_aa, cb_mult_factor, k_on, k_off,
            kappa_am, k_am, v_am, crosslinker_length, skin_distance, total_energy, cutoff_radius;
        std::vector<std::vector<interaction>> am_interaction;
        vector actin_crosslink_ratio;
        std::vector<int> actin_n_bonds;
        std::vector<std::pair<std::vector<int>, std::vector<int>>> actin_neighbors_by_species;
        vector actin_basic_tension; 
        gsl_rng * rng;
        std::string filename;
        // Lists to prevent simultaneous modifications
        std::vector<std::vector<vec>> actin_forces_temp;
        std::vector<std::vector<vec>> myosin_forces_temp, myosin_velocities_temp;
        std::vector<std::vector<double>> actin_angular_forces_temp, myosin_angular_forces_temp, actin_cb_strengths_temp;
        std::vector<utils::MoleculeConnection> actinIndicesPerMyosin_temp;
        std::vector<gsl_rng*> rng_engines;


        Sarcomere() {
            box.resize(2);
            box[0] = 0;
            box[1] = 0;
            k_aa = 0;
            cb_mult_factor = 0;
            total_energy = 0;
        }

        Sarcomere(int& n_actins, int& n_myosins, vector box0, double& actin_length, double& myosin_length,
                  double& myosin_radius, double& crosslinker_length, double& k_on, double& k_off,
                  double& cb_mult_factor, double& k_aa, double& kappa_aa, double& k_am, double& kappa_am, double& v_am,
                  std::string& filename, gsl_rng* rng, int& seed)
            : actin(n_actins, actin_length, box0, rng),
              myosin(n_myosins, myosin_length, myosin_radius, box0, rng),
              myosinIndicesPerActin(n_actins),
              actinIndicesPerMyosin(n_myosins),
              actinIndicesPerActin(n_actins),
              neighbor_list(0.0, box0, 0.0),
              actin_actin_bonds(n_actins, std::vector<int>(n_actins, 0)),
            actin_forces_temp(omp_get_max_threads(), std::vector<vec>(n_actins, {0, 0})),
            myosin_forces_temp(omp_get_max_threads(), std::vector<vec>(n_myosins, {0, 0})),
            myosin_velocities_temp(omp_get_max_threads(), std::vector<vec>(n_myosins, {0, 0})),
            actin_angular_forces_temp(omp_get_max_threads(), std::vector<double>(n_actins, 0)),
            myosin_angular_forces_temp(omp_get_max_threads(), std::vector<double>(n_myosins, 0)),
            actin_cb_strengths_temp(omp_get_max_threads(), std::vector<double>(n_actins, 0)),
            actinIndicesPerMyosin_temp(omp_get_max_threads(), utils::MoleculeConnection(n_myosins)),
            rng_engines(omp_get_max_threads(), nullptr) 

            {
            box.resize(2);
            box[0] = box0[0];
            box[1] = box0[1];
            this->k_aa = k_aa;
            this->kappa_aa = kappa_aa;
            this->k_on = k_on;
            this->k_off = k_off;
            this->cb_mult_factor = cb_mult_factor;
            this->k_am = k_am;
            this->kappa_am = kappa_am;
            this->v_am = v_am;
            this->crosslinker_length = crosslinker_length;
            this->skin_distance = skin_distance;
            this->filename = filename;
            this->rng = rng;
            cutoff_radius = std::max(actin_length, myosin_length) +
                            std::max(2 * myosin_radius, crosslinker_length);
            double skin_distance = 0.05 * cutoff_radius;
            neighbor_list = NeighborList(cutoff_radius + skin_distance, box, skin_distance / 2);
            neighbor_list.initialize(actin.center, myosin.center);
            actin_actin_bonds_prev = actin_actin_bonds;
            actin_basic_tension.resize(n_actins);
            actin_crosslink_ratio.resize(n_actins);
            actin_n_bonds.resize(n_actins);
            am_interaction.resize(n_actins);
            for (int i = 0; i < n_actins; i++) {
                am_interaction[i].resize(n_myosins);
            }
            actin.register_feature("myosin_binding_ratio");
            actin.register_feature("crosslink_ratio");
            actin.register_feature("partial_binding_ratio");
            for (int i = 0; i < n_myosins; i++) {
                myosin.theta[i] = 0;
            }
            // Initialize thread-local RNGs
            for (int t = 0; t < omp_get_max_threads(); ++t) {
                rng_engines[t] = gsl_rng_alloc(gsl_rng_mt19937); 
                gsl_rng_set(rng_engines[t], seed + t);           
            }
        }

        ~Sarcomere() {}

        void myosin_on_a_lattice(){
            std::vector<vector> myosin_positions;
            box[0] = 16;
            box[1] = 8;
            myosin_positions = {{-4,-3},{4,-3},{-4,0},{4,0}, {-4,3},{4,3}};
            for (int i = 0; i < myosin_positions.size(); i++){
                myosin.center[i].x = myosin_positions[i][0];
                myosin.center[i].y = myosin_positions[i][1];
                myosin.theta[i] = 0;
            }
            //myosin.n = myosin_positions.size();
            myosin.update_endpoints();
            update_system();
        }

        void partial_fix(){
            std::vector<vector> myosin_positions;
            myosin_positions = {{0,-0.5},{0,0.5}};
            for (int i = 0; i < myosin_positions.size(); i++){
                myosin.center[i].x = myosin_positions[i][0];
                myosin.center[i].y = myosin_positions[i][1];
                myosin.theta[i] = 0;
            }
            for (int i = 0; i <myosin.n; i++){
                myosin.theta[i] = 0;
            }
            //myosin.n = myosin_positions.size();
            myosin.update_endpoints();
            update_system();
        }

        void cb(){
            std::vector<vector> actin_positions = {{0.5,-0.015},{-0.5,0.015}};
            for (int i = 0; i < actin_positions.size(); i++){
                actin.center[i].x = actin_positions[i][0];
                actin.center[i].y = actin_positions[i][1];
            }
            actin.theta[0] = 0;
            actin.theta[1] = M_PI;
            std::vector<vector> myosin_positions = {{-1.33,0},{1.33,0}};
            for (int i = 0; i < myosin_positions.size(); i++){
                myosin.center[i].x = myosin_positions[i][0];
                myosin.center[i].y = myosin_positions[i][1];
                myosin.theta[i] = 0;
            }
        }

        void sarcomeric_structure(){
            std::vector<vector> myosin_positions;
            box[0] = 5.32;
            box[1] = 5.32;
            myosin_positions = {{-1.33, -1.0}, {1.33, -1.0}, {-1.33, 0.0}, {1.33, 0.0}, {-1.33, 1.0}, {1.33, 1.0}};
            for (int i = 0; i < myosin_positions.size(); i++){
                myosin.center[i].x = myosin_positions[i][0];
                myosin.center[i].y = myosin_positions[i][1];
                myosin.theta[i] = 0;
            }
            //myosin.n = myosin_positions.size();
            myosin.update_endpoints();
            //alpha_actinin.n = alpha_actinin_positions.size();
            std::vector<vector> actin_positions;
            actin_positions = {{0.5, -1.15}, {0.5, -1.0}, {0.5, -0.85}, {0.5, -0.15}, {0.5, 0.0}, {0.5, 0.15}, 
                       {0.5, 0.85}, {0.5, 1.0}, {0.5, 1.15}, {-2.16, -1.15}, {-2.16, -1.0}, {-2.16, -0.85}, 
                       {-2.16, -0.15}, {-2.16, 0.0}, {-2.16, 0.15}, {-2.16, 0.85}, {-2.16, 1.0}, {-2.16, 1.15}};
            for (int i = 0; i < actin_positions.size(); i++){
                actin.center[i].x = actin_positions[i][0];
                actin.center[i].y = actin_positions[i][1];
                actin.theta[i] = 0;
            }
            int n = actin_positions.size();
            actin_positions = {{-0.5, -1.15}, {-0.5, -1.0}, {-0.5, -0.85}, {-0.5, -0.15}, {-0.5, 0.0}, {-0.5, 0.15}, 
                       {-0.5, 0.85}, {-0.5, 1.0}, {-0.5, 1.15},{2.16, -1.15}, {2.16, -1.0}, {2.16, -0.85}, 
                       {2.16, -0.15}, {2.16, 0.0}, {2.16, 0.15}, {2.16, 0.85}, {2.16, 1.0}, {2.16, 1.15}};
            for (int i = 0; i < actin_positions.size(); i++){
                actin.center[i+n].x = actin_positions[i][0];
                actin.center[i+n].y = actin_positions[i][1];
                actin.theta[i+n] = M_PI;
            }
            //actin.n = actin_positions.size()*2;
            actin.update_endpoints();
            update_system();
            //save_state();
        }

        void update_system() {
            omp_set_num_threads(1);
            neighbor_list.set_species_positions(actin.center, myosin.center);
            if (neighbor_list.needs_rebuild()) {
                neighbor_list.rebuild_neighbor_list();
            }

            // Clear previous neighbors data and prepare to store new neighbors
            actin_neighbors_by_species.clear();
            actin_neighbors_by_species.resize(actin.center.size());
            #pragma omp parallel
            {
                // Step 2: _set_to_zero
                _set_to_zero();
                #pragma omp barrier // Ensure _set_to_zero completes for all threads
                // Step 3: _update1
                _update1();

                #pragma omp barrier // Ensure _update1 completes before proceeding
                // Step 4: _update2
                _update2();

                #pragma omp barrier // Ensure _update2 completes before reducing thread-local values
                // Step 5: Reduce `actin_cb_strengths_temp` into `actin.cb_strength`
                #pragma omp for
                for (int i = 0; i < actin.n; ++i) {
                    for (int t = 0; t < omp_get_num_threads(); ++t) {
                        actin.cb_strength[i] += actin_cb_strengths_temp[t][i];
                    }
                }
                // concatenate actinIndicesPerMyosin_temp
                #pragma omp for
                for (int i = 0; i < myosin.n; ++i){
                    for (int t = 0; t < omp_get_num_threads(); ++t){
                        auto indices = actinIndicesPerMyosin_temp[t].getConnections(i);
                        for (int j = 0; j < indices.size(); j++){
                            actinIndicesPerMyosin.addConnection(i, indices[j]);
                        }
                    }
                }

                #pragma omp barrier // Ensure the reduction completes before moving to _update3
                // Step 6: _update3
                _update3();

                #pragma omp barrier // Ensure _update3 completes before reducing thread-local forces

                // Step 7: Reduce `actin_forces_temp` and `actin_angular_forces_temp`
                #pragma omp for
                for (int i = 0; i < actin.n; ++i) {
                    for (int t = 0; t < omp_get_num_threads(); ++t) {
                        actin.force[i] += actin_forces_temp[t][i];
                        actin.angular_force[i] += actin_angular_forces_temp[t][i];
                    }
                }

                #pragma omp barrier // Ensure actin forces are reduced before processing myosin

                // Step 8: Reduce `myosin_forces_temp`, `myosin_velocities_temp`, and `myosin_angular_forces_temp`
                #pragma omp for
                for (int i = 0; i < myosin.n; ++i) {
                    for (int t = 0; t < omp_get_num_threads(); ++t) {
                        myosin.force[i] += myosin_forces_temp[t][i];
                        myosin.velocity[i] += myosin_velocities_temp[t][i];
                        myosin.angular_force[i] += myosin_angular_forces_temp[t][i];
                    }
                }
            }
        }


        void _set_to_zero() {
            int thread_id = omp_get_thread_num();
            // Thread-local temporary lists for forces and velocities
            auto& local_actin_forces = actin_forces_temp[thread_id];
            auto& local_actin_angular_forces = actin_angular_forces_temp[thread_id];
            auto& local_actin_cb_strengths = actin_cb_strengths_temp[thread_id];
            auto& local_myosin_forces = myosin_forces_temp[thread_id];
            auto& local_myosin_velocities = myosin_velocities_temp[thread_id];
            auto& local_myosin_angular_forces = myosin_angular_forces_temp[thread_id];
            auto& local_actinIndicesPerMyosin = actinIndicesPerMyosin_temp[thread_id];

            #pragma omp for
            for (int i = 0; i < actin.n; i++) {
                assert(thread_id == omp_get_thread_num());
                actin.update_endpoints(i);
                local_actin_forces[i] = {0, 0};
                local_actin_angular_forces[i] = 0;
                local_actin_cb_strengths[i] = 0;
                actin.force[i] = {0, 0};
                actin.angular_force[i] = 0;
                actin.velocity[i] = {0, 0};
                actin.tension[i] = 0;
                actin.cb_strength[i] = 0;
                actin_basic_tension[i] = 0;
                actin_n_bonds[i] = 0;
                actin_crosslink_ratio[i] = 1;
                actin["myosin_binding_ratio"][i] = 0;
                actin["crosslink_ratio"][i] = 1;
                actin["partial_binding_ratio"][i] = 0;
                myosinIndicesPerActin.deleteAllConnections(i);
                actinIndicesPerActin.deleteAllConnections(i);
                for (int j = 0; j < actin.n; j++){
                    actin_actin_bonds_prev[i][j] = actin_actin_bonds[i][j];
                    actin_actin_bonds[i][j] = 0;
                }
            }

            #pragma omp for
            for (int i = 0; i < myosin.n; i++) {
                myosin.update_endpoints(i);
                local_myosin_forces[i] = {0, 0};
                local_myosin_velocities[i] = {0, 0};
                local_myosin_angular_forces[i] = 0;
                myosin.force[i] = {0, 0};
                myosin.angular_force[i] = 0;
                myosin.velocity[i] = {0, 0};
                actinIndicesPerMyosin.deleteAllConnections(i);
                local_actinIndicesPerMyosin.deleteAllConnections(i);
            }
            #pragma omp for
            for (size_t i = 0; i < actin.center.size(); i++) {
                // Get actin and myosin neighbors for actin particle `i`
                actin_neighbors_by_species[i] = neighbor_list.get_neighbors_by_type(i);
            }
        }

        void _update1() {
            int thread_id = omp_get_thread_num();
            #pragma omp for
            for (int i = 0; i < actin.n; i++) {
                auto& local_actinIndicesPerMyosin = actinIndicesPerMyosin_temp[thread_id];
                double partial_binding_ratio = 0;
                double myosin_binding_ratio_sum = 0;
                int myosin_index = -1;
                auto myosin_neighbors = actin_neighbors_by_species[i].second;
                for (int index = 0; index < myosin_neighbors.size(); index++) {
                    int j = myosin_neighbors[index];
                    am_interaction[i][j] = geometry::analyze_am(
                        actin.left_end[i], actin.right_end[i], myosin.left_end[j], myosin.right_end[j],
                        myosin.radius, box);
                    myosin_binding_ratio_sum += am_interaction[i][j].myosin_binding_ratio;
                    if (am_interaction[i][j].partial_binding_ratio > partial_binding_ratio) {
                        partial_binding_ratio = am_interaction[i][j].partial_binding_ratio;
                        myosin_index = j;
                    }
                }
                if (myosin_index >= 0) {
                    local_actinIndicesPerMyosin.addConnection(myosin_index, i);
                    myosinIndicesPerActin.addConnection(i, myosin_index);
                    double myosin_binding_ratio = am_interaction[i][myosin_index].myosin_binding_ratio;
                    actin["myosin_binding_ratio"][i] = myosin_binding_ratio;
                    actin["crosslink_ratio"][i] = am_interaction[i][myosin_index].crosslinkable_ratio-(myosin_binding_ratio_sum-myosin_binding_ratio);
                    actin["partial_binding_ratio"][i] = am_interaction[i][myosin_index].partial_binding_ratio;
                }
            }
        }

        void _update2() {
            int thread_id = omp_get_thread_num();
            #pragma omp for
            for (int i = 0; i < actin.n; i++) {
                std::vector <int> actin_neighbors = actin_neighbors_by_species[i].first;
                std::vector<int> myosin_indices = myosinIndicesPerActin.getConnections(i);
                for (int index = 0; index < myosin_indices.size(); index++) {
                    int j = myosin_indices[index];
                    vec velocity = {v_am * cos(actin.theta[i]), v_am * sin(actin.theta[i])};
                    double a_m_angle = actin.theta[i] - myosin.theta[j];
                    double abs_cos = std::abs(std::cos(a_m_angle));
                    assert(i >= 0 && i < actin.n);
                    assert(j >= 0 && j < myosin.n);
                    double myosin_ratio = std::min(am_interaction[i][j].partial_binding_ratio * 3, 1.0);
                    actin_basic_tension[i] = myosin_ratio * abs_cos;
                    actin.velocity[i] = velocity;
                    std::vector<int> cb_indices;
                    std::vector<double> cb_strengths;
                    bool add_connection = true;
                    for (int index = 0; index < actin_neighbors.size(); index++){
                        int j = actin_neighbors[index];
                        if (i>j){
                            double strength = _get_cb_strength(i,j);
                            if (strength>EPS){
                                cb_indices.push_back(j);
                                cb_strengths.push_back(strength);
                            }
                        }
                    }
                    _set_cb(i, cb_indices,cb_strengths);
                }
            }
        }

        void _update3() {
            int thread_id = omp_get_thread_num();
            // Thread-local temporary lists for forces and velocities
            auto& local_actin_forces = actin_forces_temp[thread_id];
            auto& local_actin_angular_forces = actin_angular_forces_temp[thread_id];
            auto& local_myosin_forces = myosin_forces_temp[thread_id];
            auto& local_myosin_velocities = myosin_velocities_temp[thread_id];
            auto& local_myosin_angular_forces = myosin_angular_forces_temp[thread_id];
            #pragma omp for
            for (int i = 0; i < actin.n; i++) {
                std::vector<int> myosin_indices = myosinIndicesPerActin.getConnections(i);
                vec velocity = {v_am*cos(actin.theta[i]),v_am*sin(actin.theta[i])};
                for (int index = 0; index < myosin_indices.size(); index++) {
                    int j = myosin_indices[index];
                    double binding_ratio = am_interaction[i][j].myosin_binding_ratio;
                    double tension = actin_basic_tension[i] * (1 / cb_mult_factor + actin.cb_strength[i]);
                    vector force_vec = compute_am_force_and_energy(
                        actin, myosin, i, j, box, k_am * tension, kappa_am * tension, myosin.radius, true);
                    local_actin_forces[i].x += force_vec[0];
                    local_actin_forces[i].y += force_vec[1];
                    local_myosin_forces[j].x -= force_vec[0];
                    local_myosin_forces[j].y -= force_vec[1];
                    local_actin_angular_forces[i] += force_vec[2];
                    local_myosin_angular_forces[j] += force_vec[3];
                    double binding_ratio_adjusted = std::min(am_interaction[i][j].partial_binding_ratio*3,1.0);
                    double v_factor = 0.2*std::max(actin.cb_strength[i]-0.1,0.0)* std::max(1-binding_ratio_adjusted,0.0);
                    local_myosin_velocities[j] -= v_factor * velocity;
                }
                if (actin.cb_strength[i]>EPS){
                    double v_factor = std::max(0.1-actin.cb_strength[i],0.0);
                    actin.velocity[i].x *=v_factor;
                    actin.velocity[i].y *=v_factor;
                }
            }
            #pragma omp for
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

            #pragma omp for
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

        void _myosin_repulsion(int& i, int& j){
            int thread_id = omp_get_thread_num();
            vec center_displacement = myosin.center[i] - myosin.center[j];
            center_displacement.pbc_wrap(box);
            double center_distance = center_displacement.norm();
            if (center_distance<=2*myosin.radius+myosin.length){
                auto& local_myosin_forces = myosin_forces_temp[thread_id]; 
                auto result = geometry::segment_segment_distance_w_normal(myosin.left_end[i], 
                    myosin.right_end[i], myosin.left_end[j], myosin.right_end[j], box);
                double distance = result.first;
                if (distance<2*myosin.radius){
                    vec normal_vector = result.second["vector"];
                    double norm = normal_vector.norm();
                    double factor = std::min((100*(2*myosin.radius - distance)/(2*myosin.radius)),20.);
                    if (norm==0){
                        normal_vector.x = center_displacement.x;
                        normal_vector.y = center_displacement.y;
                        normal_vector.normalize();
                    }
                    else {
                        normal_vector = normal_vector/norm;
                    }
                    //check if any of the myosin is bound to cb-forming actin
                    auto actin_indices_i = actinIndicesPerMyosin.getConnections(i);
                    double cb_strength_i = 0;
                    //printf("cb_strength_i: %f \n", cb_strength_i);
                    for (int index = 0; index < actin_indices_i.size(); index++){
                        int k = actin_indices_i[index];
                        //printf("i: %d, k: %d, cb: %f, total cb %f \n", i, k, actin.cb_strength[k],cb_strength_i);
                        cb_strength_i += actin.cb_strength[k];
                    }
                    auto actin_indices_j = actinIndicesPerMyosin.getConnections(j);
                    double cb_strength_j = 0;
                    for (int index = 0; index < actin_indices_j.size(); index++){
                        int k = actin_indices_j[index];
                        //printf("j: %d, k: %d, cb: %f, total cb %f \n", j, k, actin.cb_strength[k], cb_strength_j);
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
                    double force_x = factor*normal_vector.x;
                    double force_y = factor*normal_vector.y;
                    //printf("thread_id: %d, i: %d, j: %d, force: %f, %f, total force: %f, %f, cb: %f, %f \n",
                     //thread_id, i, j, force_x, force_y, local_myosin_forces[i].x, local_myosin_forces[i].y, cb_strength_i, cb_strength_j);
                }

            }
        }

        void _actin_repulsion(int& i, int& j){
            int thread_id = omp_get_thread_num();
            vec center_displacement = actin.center[i] - actin.center[j];
            center_displacement.pbc_wrap(box);
            double center_distance = center_displacement.norm();
            if (center_distance<=crosslinker_length+actin.length){
                auto& local_actin_forces = actin_forces_temp[thread_id];
                auto result = geometry::segment_segment_distance_w_normal(actin.left_end[i], 
                    actin.right_end[i], actin.left_end[j], actin.right_end[j], box);
                double distance = result.first;
                double min_dist = crosslinker_length/2;
                if (distance < min_dist){
                    vec normal_vector = result.second["vector"];
                    double norm = normal_vector.norm();
                    double factor = std::min((4*(min_dist - distance)/min_dist),1.);
                    if (norm==0){
                        normal_vector.x = center_displacement.x;
                        normal_vector.y = center_displacement.y;
                        normal_vector.normalize();
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

        double _get_cb_strength(int& i, int& j){
            double angle = actin.theta[i] - actin.theta[j];
            double cos_angle = std::cos(angle);
            bool crosslink = false;
            vec right_i, right_j;
            if (cos_angle<0 && actin_crosslink_ratio[i]>EPS && actin_crosslink_ratio[j]>EPS){
                right_i = actin.left_end[i] + (actin.right_end[i]-actin.left_end[i])*actin_crosslink_ratio[i];
                right_j = actin.left_end[j] + (actin.right_end[j]-actin.left_end[j])*actin_crosslink_ratio[j];
                double distance = geometry::segment_segment_distance(actin.left_end[i], 
                    actin.right_end[i], actin.left_end[j], actin.right_end[j], box);
                if (distance<crosslinker_length){
                    crosslink = true;
                }
            }
            if (crosslink){ 
                double strength = 0;
                double abs_cos_angle = std::abs(cos_angle);
                if (actin_basic_tension[i]>EPS && actin_basic_tension[j]>EPS){
                    int myosin_index_i = myosinIndicesPerActin.getConnections(i)[0];
                    int myosin_index_j = myosinIndicesPerActin.getConnections(j)[0];
                    if (myosin_index_i != myosin_index_j){
                        double crosslink_ratio = std::sqrt(std::min(actin_crosslink_ratio[i]*3,1.0) * std::min(actin_crosslink_ratio[j]*3,1.0));
                        //strength = crosslink_ratio * actin_basic_tension[i]*actin_basic_tension[j]*abs_cos_angle;
                        strength =  actin_basic_tension[i]*actin_basic_tension[j]*(abs_cos_angle*abs_cos_angle);
                    }
                }
                else{
                    strength = 1/cb_mult_factor*abs_cos_angle;
                }
                return strength;
            }
            return 0.0;
        }
        
        void _set_cb(int& i, int& j, double& normalized_strength, bool& add_connection){
            int thread_id = omp_get_thread_num();
            if (actin_actin_bonds[i][j] == 1) {
                std::cout << "Skipping duplicate bond: " << i << ", " << j << "\n";
            }
            auto& local_actin_forces = actin_forces_temp[thread_id];
            auto& local_actin_angular_forces = actin_angular_forces_temp[thread_id];
            auto& local_actin_cb_strengths = actin_cb_strengths_temp[thread_id];
            auto& local_myosin_forces = myosin_forces_temp[thread_id];
            auto& local_myosin_angular_forces = myosin_angular_forces_temp[thread_id];

            if (normalized_strength<1/cb_mult_factor && (actin_n_bonds[i]>=2 || actin_n_bonds[j]>=2)){
                return;
            }
            vector force_vec;
            force_vec.resize(4);
            if (normalized_strength<EPS){
                return;
            }
            double rand = gsl_rng_uniform(rng_engines[thread_id]);
            if (actin_actin_bonds_prev[i][j] == 1) {
                double k_off_adjusted = k_off/normalized_strength;
                if (rand < k_off_adjusted){ //k_off is actually k_off * dt
                    //printf("breaking bond between %d and %d\n",i,j);
                    return;
                }
            }
            else{
                if (rand >= k_on){ //k_on is actually k_on * dt
                    // printf("not forming bond between %d and %d\n",i,j);
                    return;
                }
                else {
                    //printf("forming bond between %d and %d\n",i,j);
                }
            }
            double strength = normalized_strength*cb_mult_factor;
            //actin.tension[i] += actin_basic_tension[i]*strength;
            //actin.tension[j] += actin_basic_tension[j]*strength;
            force_vec = compute_aa_force_and_energy(actin,i, j, box,k_aa,kappa_aa); 
            local_actin_forces[i].x += force_vec[0];
            local_actin_forces[i].y += force_vec[1];
            local_actin_forces[j].x -= force_vec[0];
            local_actin_forces[j].y -= force_vec[1];
            local_actin_angular_forces[i] += force_vec[2];
            local_actin_angular_forces[j] += force_vec[3]; 
            local_actin_cb_strengths[i] += normalized_strength;
            local_actin_cb_strengths[j] += normalized_strength;
            actin_n_bonds[i] += 1;
            actin_n_bonds[j] += 1;
            if (add_connection){
                actinIndicesPerActin.addConnection(i,j);
                actinIndicesPerActin.addConnection(j,i);
                //printf("adding connection between %d and %d\n",i,j);
            }
            actin_actin_bonds[i][j] = 1;
            actin_actin_bonds[j][i] = 1;
        }

        void _set_cb(int& i, std::vector<int> indices, vector cb_strength){
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
     
        void new_file(){
            create_file(filename, actin, myosin);
        }

        void save_state(){
            append_to_file(filename, actin, myosin, total_energy, actinIndicesPerActin);
        }
        
        void load_state(){
            load_from_file(filename, actin, myosin, total_energy);
            update_system();
        }
        
};

#endif
