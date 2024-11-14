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


#include "h5_utils.h"

#include <cmath>
#include <vector>
#include <numeric>

class Sarcomere
{
    public:
        Filament actin;
        Myosin myosin;
        AlphaActinin alpha_actinin;
        //myosin indices that each actin is bound to
        utils::MoleculeConnection myosinIndicesPerActin;
        //actin indices that each myosin is bound to
        utils::MoleculeConnection actinIndicesPerMyosin;
        //actin indices that each actin forms a catch bond with
        utils::MoleculeConnection actinIndicesPerActin;
        std::vector <double> box;
        double k_aa, kappa_aa, cb_mult_factor,kappa_am, k_am, v_am;
        double total_energy;
        double ** a_m_energy;
        std::vector<double> ** a_m_normal_vec;
        std::vector<double> ** a_m_normal_start;
        std::vector<double> ** a_m_force;
        std::vector<double> ** a_m_velocity;

        double * actin_basic_tensions; 
        gsl_rng * rng;


        std::string filename;

        Sarcomere(){
            box.resize(2);
            box[0] = 0;
            box[1] = 0;
            k_aa = 0;
            cb_mult_factor = 0;
            total_energy = 0;
        }
        Sarcomere(int& n_actins, int& n_myosins, int& n_alpha_actinins, std::vector<double> box0, double& actin_length, double& myosin_length,
            double& cb_mult_factor, double& k_aa, double& kappa_aa, double& k_am, double& kappa_am,double& v_am,
            double& myosin_radius, double& alpha_actinin_radius, std::string& filename, gsl_rng * rng)
            : actin(n_actins, actin_length, box0, rng), 
                myosin(n_myosins, myosin_length, myosin_radius, box0, rng),
                alpha_actinin(n_alpha_actinins, alpha_actinin_radius,box0, rng),
                myosinIndicesPerActin(n_actins),
                actinIndicesPerMyosin(n_myosins),
                actinIndicesPerActin(n_actins)
            {
            box.resize(2);
            box[0] = box0[0];
            box[1] = box0[1];
            this->k_aa = k_aa;
            this->kappa_aa = kappa_aa;
            this->cb_mult_factor = cb_mult_factor;
            this->k_am = k_am;
            this->kappa_am = kappa_am;
            this->v_am = v_am;
            this->filename = filename;
            this->rng = rng;
            a_m_energy = new double*[actin.n];
            a_m_normal_vec = new std::vector<double>*[actin.n];
            a_m_normal_start = new std::vector<double>*[actin.n];
            a_m_force = new std::vector<double>*[actin.n];
            a_m_velocity = new std::vector<double>*[actin.n];
            for (int i = 0; i < actin.n; i++){
                a_m_energy[i] = new double[myosin.n];
                a_m_force[i] = new std::vector<double>[myosin.n];
                a_m_normal_vec[i] = new std::vector<double>[myosin.n];
                a_m_normal_start[i] = new std::vector<double>[myosin.n];
                a_m_velocity[i] = new std::vector<double>[myosin.n];
                for (int j = 0; j < myosin.n; j++){
                    a_m_energy[i][j] = 0;
                    a_m_force[i][j]= {0,0};
                    a_m_normal_vec[i][j] = {0,0};
                    a_m_normal_start[i][j] = {0,0};
                    a_m_velocity[i][j] = {0,0};
                }
            }
            actin_basic_tensions = new double[actin.n];

            
        }
        ~Sarcomere(){ 
            for (int i = 0; i < actin.n; i++){
                delete[] a_m_energy[i];
                //set to null ptr to avoid double free
                a_m_energy[i] = nullptr;
                delete[] a_m_normal_vec[i];
                a_m_normal_vec[i] = nullptr;
                delete[] a_m_normal_start[i];
                a_m_normal_start[i] = nullptr;
                delete[] a_m_force[i];
                a_m_force[i] = nullptr;
                delete[] a_m_velocity[i];
                a_m_velocity[i] = nullptr;
                
            }
            for (int i = 0; i < myosin.n; i++){
                a_m_energy[i] = nullptr;
            }
            printf("Sarcomere destructor called\n");
        }

        Sarcomere(const Sarcomere& other){
            actin = other.actin;
            myosin = other.myosin;
            alpha_actinin = other.alpha_actinin;
            box = other.box;
            k_aa = other.k_aa;
            k_am = other.k_am;
            cb_mult_factor = other.cb_mult_factor;
        }

        void myosin_on_a_lattice(){
            std::vector<std::vector<double>> myosin_positions;
            box[0] = 16;
            box[1] = 8;
            myosin_positions = {{-4,-3},{4,-3},{-4,0},{4,0}, {-4,3},{4,3}};
            for (int i = 0; i < myosin_positions.size(); i++){
                myosin.xs[i][0] = myosin_positions[i][0];
                myosin.xs[i][1] = myosin_positions[i][1];
                myosin.thetas[i] = 0;
            }
            //myosin.n = myosin_positions.size();
            myosin.update_endpoints();
            update_system();
        }
        
        void sarcomeric_structure(){
            std::vector<std::vector<double>> myosin_positions;
            box[0] = 16.2;
            box[1] = 12;
            myosin_positions = {{-4.05,-3},{4.05,-3},{-4.05,0},{4.05,0}, {-4.05,3},{4.05,3}};
            for (int i = 0; i < myosin_positions.size(); i++){
                myosin.xs[i][0] = myosin_positions[i][0];
                myosin.xs[i][1] = myosin_positions[i][1];
                myosin.thetas[i] = 0;
            }
            //myosin.n = myosin_positions.size();
            myosin.update_endpoints();
            std::vector<std::vector<double>> actin_positions;
            actin_positions = {{-6.6,-3.2},{1.5,-3.2},{-6.6,-3},{1.5,-3},{-6.6,-2.8},{1.5,-2.8},{-6.6,-0.2},{1.5,-0.2},{-6.6,0},{1.5,0},{-6.6,0.2},{1.5,0.2},{-6.6,2.8},{1.5,2.8},{-6.6,3},{1.5,3},{-6.6,3.2},{1.5,3.2}};
            for (int i = 0; i < actin_positions.size(); i++){
                actin.xs[i][0] = actin_positions[i][0];
                actin.xs[i][1] = actin_positions[i][1];
                actin.thetas[i] = 0;
            }
            int n = actin_positions.size();
            actin_positions = {{-1.5,-3.2},{6.6,-3.2},{-1.5,-3},{6.6,-3},{-1.5,-2.8},{6.6,-2.8},{-1.5,-0.2},{6.6,-0.2},{-1.5,0},{6.6,0},{-1.5,0.2},{6.6,0.2},{-1.5,2.8},{6.6,2.8},{-1.5,3},{6.6,3},{-1.5,3.2},{6.6,3.2}};
            for (int i = 0; i < actin_positions.size(); i++){
                actin.xs[i+n][0] = actin_positions[i][0];
                actin.xs[i+n][1] = actin_positions[i][1];
                actin.thetas[i+n] = M_PI;
            }
            //actin.n = actin_positions.size()*2;
            actin.update_endpoints();
            update_system();
            //save_state();
        }

        void partial_sarcomeric_structure(){
            std::vector<std::vector<double>> myosin_positions;
            box[0] = 19.2;
            box[1] = 12;
            myosin_positions = {{-4,-3},{4,-3},{-4,0},{4,0}, {-4,3},{4,3}};
            for (int i = 0; i < myosin_positions.size(); i++){
                myosin.xs[i][0] = myosin_positions[i][0];
                myosin.xs[i][1] = myosin_positions[i][1];
                myosin.thetas[i] = 0;
            }
            //myosin.n = myosin_positions.size();
            myosin.update_endpoints();
            //alpha_actinin.n = alpha_actinin_positions.size();
            std::vector<std::vector<double>> actin_positions;
            actin_positions = {{1.5,-3.2},{1.5,-3},{1.5,-2.8},{1.5,-0.2},{1.5,0},{1.5,0.2},{1.5,2.8},{1.5,3},{1.5,3.2}};
            for (int i = 0; i < actin_positions.size(); i++){
                actin.xs[i][0] = actin_positions[i][0];
                actin.xs[i][1] = actin_positions[i][1];
                actin.thetas[i] = 0;
            }
            int n = actin_positions.size();
            actin_positions = {{-1.5,-3.2},{-1.5,-3},{-1.5,-2.8},{-1.5,-0.2},{-1.5,0},{-1.5,0.2},{-1.5,2.8},{-1.5,3},{-1.5,3.2}};
            for (int i = 0; i < actin_positions.size(); i++){
                actin.xs[i+n][0] = actin_positions[i][0];
                actin.xs[i+n][1] = actin_positions[i][1];
                actin.thetas[i+n] = M_PI;
            }
            //actin.n = actin_positions.size()*2;
            actin.update_endpoints();
            update_system();
            //save_state();
        }

        void catch_bond(){
            std::vector<std::vector<double>> myosin_positions;
            box[0] = 19.2;
            box[1] = 12;
            myosin_positions = {{-4,0},{4,0}};
            for (int i = 0; i < myosin_positions.size(); i++){
                myosin.xs[i][0] = myosin_positions[i][0];
                myosin.xs[i][1] = myosin_positions[i][1];
                myosin.thetas[i] = 0;
            }
            //myosin.n = myosin_positions.size();
            myosin.update_endpoints();
            std::vector<std::vector<double>> alpha_actinin_positions;
            alpha_actinin_positions = {{0,0},{-9.6,0}};
            for (int i = 0; i < alpha_actinin_positions.size(); i++){
                alpha_actinin.xs[i][0] = alpha_actinin_positions[i][0];
                alpha_actinin.xs[i][1] = alpha_actinin_positions[i][1];
            }
            //alpha_actinin.n = alpha_actinin_positions.size();
            std::vector<std::vector<double>> actin_positions;
            actin_positions = {{1.5,0}};
            for (int i = 0; i < actin_positions.size(); i++){
                actin.xs[i][0] = actin_positions[i][0];
                actin.xs[i][1] = actin_positions[i][1];
                actin.thetas[i] = 0;
            }
            int n = actin_positions.size();
            actin_positions = {{-1.5,0}};
            for (int i = 0; i < actin_positions.size(); i++){
                actin.xs[i+n][0] = actin_positions[i][0];
                actin.xs[i+n][1] = actin_positions[i][1];
                actin.thetas[i+n] = M_PI;
            }
            //actin.n = actin_positions.size()*2;
            actin.update_endpoints();
            update_system();
            //save_state();
        }

        void update_system(){
            _set_to_zero();
            _update_connections();
            _update_velocities();
            _update_forces();
        }
        

        //initialize all interactions to zero
        void _set_to_zero(){
            actin.update_endpoints();
            myosin.update_endpoints();
            for (int i = 0; i < actin.n; i++){
                actin.forces[i][0] = 0;
                actin.forces[i][1] = 0;
                actin.angular_forces[i] = 0;
                 for (int j = 0; j < myosin.n; j++){
                    a_m_velocity[i][j] = {0,0};
                    a_m_force[i][j] = {0,0};
                }
                actin.velocities[i][0] = 0;
                actin.velocities[i][1] = 0;
                actin.tensions[i] = 0;
                actin_basic_tensions[i] = 0;
                actin.cb_strengths[i] = 0;
            }
            for (int i = 0; i < myosin.n; i++){
                myosin.forces[i][0] = 0;
                myosin.forces[i][1] = 0;
                myosin.angular_forces[i] = 0;
                myosin.velocities[i][0] = 0;
                myosin.velocities[i][1] = 0;
                for (int j = 0; j < i; j++){
                    _myosin_repulsion(i,j);
                }
            }
        }

        void _update_connections(){
            std::vector<int> n_actins_per_myosin;
            for (int j = 0; j < myosin.n; j++){
                n_actins_per_myosin.push_back(0);
            }
            //prioritize existing connections
            for (int i = 0; i < actin.n; i++){
                std::vector<int> myosin_indices = myosinIndicesPerActin.getConnections(i);
                for (int index = 0; index < myosin_indices.size(); index++){
                    int j = myosin_indices[index];
                    a_m_energy[i][j] = _actin_myosin_energy(i,j);
                    if (a_m_energy[i][j]!=0){
                        n_actins_per_myosin[j]++;
                    }
                    else{
                        myosinIndicesPerActin.deleteConnection(i,j);
                        actinIndicesPerMyosin.deleteConnection(j,i);
                    }
                }
            }
            for (int j = 0; j < myosin.n; j++){
                if (n_actins_per_myosin[j]>10){
                    continue;
                }
                for (int i = 0; i < actin.n; i++){
                    if (myosinIndicesPerActin.getConnections(i).size()==0){
                        a_m_energy[i][j] = _actin_myosin_energy(i,j);
                        if (a_m_energy[i][j]!=0){
                            n_actins_per_myosin[j]++;
                            myosinIndicesPerActin.addConnection(i,j);
                            actinIndicesPerMyosin.addConnection(j,i);
                        }
                    }
                    if (n_actins_per_myosin[j]>10){
                        break;
                    }
                }
            }

        }

        void _update_velocities(){
            //Assumes that the connections have been updated
            //To avoid double counting, we'll only consider catch-bond pairs (i,j) where i>j
            std::vector<int> myosin_bound_actin_indices;
            for (int i = 0; i < actin.n; i++){
                std::vector <int> myosin_indices = myosinIndicesPerActin.getConnections(i);
                for (int index = 0; index < myosin_indices.size(); index++){
                    int j = myosin_indices[index];
                    a_m_velocity[i][j][0] = v_am*cos(actin.thetas[i]);
                    a_m_velocity[i][j][1] = v_am*sin(actin.thetas[i]);
                    double a_m_angle = actin.thetas[i] - myosin.thetas[j];
                    double abs_cos = std::abs(std::cos(a_m_angle));
                    if (abs_cos>0.5){
                        double ratio = utils::find_boundary_points(actin.left_endpts[i],actin.right_endpts[i],
                        myosin.left_endpts[j],myosin.right_endpts[j],myosin.radius,box);
                        printf("ratio: %f\n",ratio);
                        ratio = std::min(1.0,ratio*3);
                        actin.tensions[i] = ratio*std::abs(std::cos(a_m_angle));
                        actin_basic_tensions[i] = actin.tensions[i];
                    }
                    actin.velocities[i][0] += a_m_velocity[i][j][0];
                    actin.velocities[i][1] += a_m_velocity[i][j][1];
                    myosin.velocities[j][0] += -a_m_velocity[i][j][0];
                    myosin.velocities[j][1] += -a_m_velocity[i][j][1];
                }
                
                std::vector <int> cb_actin_indices = actinIndicesPerActin.getConnections(i);
                if (cb_actin_indices.size()>0){
                    for (int index = 0; index < cb_actin_indices.size(); index++){
                        int j = cb_actin_indices[index];
                        if (i<=j){
                            continue;
                        }
                        if (actin.cb_strengths[i]>=1 || actin.cb_strengths[j]>=1){
                            actinIndicesPerActin.deleteConnection(i,j);
                            actinIndicesPerActin.deleteConnection(j,i);
                            continue;
                        }
                        auto result = utils::segment_segment_distance_w_normal(actin.left_endpts[i], 
                        actin.right_endpts[i], actin.left_endpts[j], actin.right_endpts[j], box);
                        double distance = result.first;
                        double angle = actin.thetas[i] - actin.thetas[j];
                        std::map<std::string, std::vector<double>> normal_vector_dict = result.second;
                        std::vector<double> normal_vector = normal_vector_dict["vector"];
                        if (distance<2*alpha_actinin.radius && std::cos(angle)<0){
                            double cb_interaction = _set_cb_interactions(i,j,normal_vector);
                            actin.cb_strengths[i] += cb_interaction;
                            actin.cb_strengths[j] += cb_interaction;
                        }
                        else {
                            actinIndicesPerActin.deleteConnection(i,j);
                            actinIndicesPerActin.deleteConnection(j,i);
                            //printf("deleted connection between actin %d and actin %d due to distance \n",i,j);
                        }   
                    }  
                }

                if (actin.cb_strengths[i]<1)
                {
                    bool add_connection = true;
                    _search_cb(i, myosin_bound_actin_indices,cb_actin_indices, add_connection);
                   
                    if (actin.cb_strengths[i]<1)
                    {
                        bool add_connection = false;
                        //indices to search over is all indices < i
                        std::vector<int> indices_to_search= std::vector<int>(i);
                        std::iota(indices_to_search.begin(), indices_to_search.end(), 0);
                        _search_cb(i,indices_to_search, myosin_bound_actin_indices, add_connection);
                    }

                }
                if (actin_basic_tensions[i]!=0){
                    myosin_bound_actin_indices.push_back(i);
                }
                if (actin.cb_strengths[i]>1){
                    actin.cb_strengths[i] = 1;
                }
            }
        }
        
        void _update_forces(){
            double angle;
            for (int i = 0; i < actin.n; i++){
                std::vector <int> myosin_indices = myosinIndicesPerActin.getConnections(i);
                for (int index = 0; index < myosin_indices.size(); index++){
                    int j = myosin_indices[index];
                    double a_m_distance = utils::norm(a_m_normal_vec[i][j]);
                    bool turn_on_spring = (a_m_distance>0.8*myosin.radius);
                    double tension = actin_basic_tensions[i]*(1+actin.cb_strengths[i]*cb_mult_factor);
                    std::vector<double> force_vec =  compute_am_force_and_energy(actin, myosin, i, j, box, k_am*tension,kappa_am*tension, turn_on_spring);
                    a_m_force[i][j][0] = force_vec[0];
                    a_m_force[i][j][1] = force_vec[1];
                    actin.angular_forces[i] +=force_vec[2];
                    myosin.angular_forces[j] +=force_vec[3];
                    if (actin.cb_strengths[i]<0 || actin.cb_strengths[i]>1){
                        printf("actin %d has cb strength %f \n",i,actin.cb_strengths[i]);
                    exit(1);
                    }
                    actin.forces[i][0] += a_m_force[i][j][0];
                    actin.forces[i][1] += a_m_force[i][j][1];
                    myosin.forces[j][0] += -a_m_force[i][j][0];
                    myosin.forces[j][1] += -a_m_force[i][j][1];
  
                }
                //check if the pulling force and velocity are in the same direction
                if (actin.cb_strengths[i]>0){
                    //double v_dot_f = actin.velocities[i][0]*actin.forces[i][0]+actin.velocities[i][1]*actin.forces[i][1];
                    double v_factor = std::max(0.9-actin.cb_strengths[i],0.0);
                    actin.velocities[i][0] *=v_factor;
                    actin.velocities[i][1] *=v_factor;  
                }
            }
        }

        std::vector<double> _get_titin_vector(int& myosin_index, int& actin_index){

            double optimal_len = 0.9*(0.5*myosin.length+actin.length);
            std::vector<double> titin_vector;
            std::vector<double> displacement_1;
            displacement_1.resize(2);
            displacement_1[0] = myosin.xs[myosin_index][0]-actin.left_endpts[actin_index][0];
            displacement_1[1] = myosin.xs[myosin_index][1]-actin.left_endpts[actin_index][1];
            utils::pbc_wrap(displacement_1, box);
            double length_1 = utils::norm(displacement_1);
            std::vector<double> displacement_2;
            displacement_2.resize(2);
            displacement_2[0] =actin.right_endpts[actin_index][0]-myosin.xs[myosin_index][0];
            displacement_2[1] =actin.right_endpts[actin_index][1]-myosin.xs[myosin_index][1];
            utils::pbc_wrap(displacement_2, box);
            double length_2 = utils::norm(displacement_2);
            if (length_1<length_2){
                titin_vector = displacement_1; 
                titin_vector[0]-=titin_vector[0]*optimal_len/length_1;
                titin_vector[1]-=titin_vector[1]*optimal_len/length_1;
            }
            else{
                titin_vector = displacement_2;
                titin_vector[0]-=titin_vector[0]*optimal_len/length_2;  
                titin_vector[1]-=titin_vector[1]*optimal_len/length_2;  
            }
            return titin_vector;
        }
        double _set_cb_interactions(int& i, int& j, std::vector<double> normal_vector){
            double k;
            double kappa;
            double strength=0;
            double angle = actin.thetas[i] - actin.thetas[j];
            double cos_angle =std::abs(std::cos(angle));
            std::vector<double> force_vec;
            force_vec.resize(4);
            if (actin_basic_tensions[i]!=0 && actin_basic_tensions[i]!=0){
                int myosin_index_i = myosinIndicesPerActin.getConnections(i)[0];
                int myosin_index_j = myosinIndicesPerActin.getConnections(j)[0];
                if (myosin_index_i != myosin_index_j 
                    && utils::check_cb(actin.left_endpts[i],actin.right_endpts[i],myosin.left_endpts[myosin_index_i],myosin.right_endpts[myosin_index_i],
                    actin.left_endpts[j],actin.right_endpts[j],myosin.left_endpts[myosin_index_j],myosin.right_endpts[myosin_index_j],box))
                 {
                    strength = std::sqrt(actin_basic_tensions[i]*actin_basic_tensions[j])*cos_angle;
                    strength = std::min(std::min(strength, 1-actin.cb_strengths[i]), 1-actin.cb_strengths[j]);
                    k = cb_mult_factor*k_aa*strength;
                    kappa = cb_mult_factor*kappa_aa*strength;
                    actin.tensions[i] += actin_basic_tensions[i]*strength*cb_mult_factor;
                    actin.tensions[j] += actin_basic_tensions[j]*strength*cb_mult_factor;
                    // std::vector<double> titin_vector = _get_titin_vector(myosin_index_j, i);
                    // actin.forces[i][0] += titin_vector[0]*5*k;
                    // actin.forces[i][1] += titin_vector[1]*5*k;
                    // titin_vector = _get_titin_vector(myosin_index_i, j);
                    // actin.forces[j][0] += titin_vector[0]*5*k;
                    // actin.forces[j][1] += titin_vector[1]*5*k;
                    if (strength>1){
                        printf("strength is greater than 1\n");
                    }
                    force_vec = compute_aa_force_and_energy(actin, myosin,i,myosin_index_i, j, myosin_index_j,box,k,kappa);
                }
                else{
                    strength = 1/cb_mult_factor*cos_angle;
                    strength = std::min(std::min(strength, 1-actin.cb_strengths[i]), 1-actin.cb_strengths[j]);
                    k = k_aa*cos_angle;
                    kappa = kappa_aa*cos_angle;
                    force_vec = compute_aa_force_and_energy(actin,i, j, box,k,kappa);
                }
            }
            else{
                strength = 1/cb_mult_factor*cos_angle;
                strength = std::min(std::min(strength, 1-actin.cb_strengths[i]), 1-actin.cb_strengths[j]);
                k = k_aa*cos_angle;
                kappa = kappa_aa*cos_angle;
                force_vec = compute_aa_force_and_energy(actin,i, j, box,k,kappa);
            }
            // actin.forces[i][0] += normal_vector[0]*k;
            // actin.forces[i][1] += normal_vector[1]*k;
            // actin.forces[j][0] += -normal_vector[0]*k;
            // actin.forces[j][1] += -normal_vector[1]*k;
            if (force_vec.size()!=4){
                printf("force vec size is not 4\n");
                exit(1);
            }
            actin.forces[i][0] += force_vec[0];
            actin.forces[i][1] += force_vec[1];
            actin.forces[j][0] -= force_vec[0];
            actin.forces[j][1] -= force_vec[1];
            actin.angular_forces[i] += force_vec[2];
            actin.angular_forces[j] += force_vec[3];
            // double dist_sq = normal_vector[0]*normal_vector[0]+normal_vector[1]*normal_vector[1];
            // actin.angular_forces[i] += std::sin(angle)*0.5*k/cos_angle*(dist_sq+0.1);
            // actin.angular_forces[j] += -std::sin(angle)*0.5*k/cos_angle*(dist_sq+0.1); 
            //check if the forces are nan
            if (std::isnan(actin.forces[i][0]) || std::isnan(actin.forces[i][1]) || std::isnan(actin.forces[j][0]) || std::isnan(actin.forces[j][1])){
                printf("k %f\n", k);
                printf("actin forces are nan\n");
                printf("cos_angle %f\n", cos_angle);
                exit(0);
            }
            //check if angular forces are nan
            if (std::isnan(actin.angular_forces[i]) || std::isnan(actin.angular_forces[j])){
                printf("k %f\n", k);
                printf("cos_angle %f\n", cos_angle);
                printf("actin angular forces are nan\n");
                exit(0);
            }
            //printf("strength between actin %d and actin %d is %f\n",i,j,strength);
            //printf("cb strength for actin %d and actin %d are %f, %f\n",i,j,actin.cb_strengths[i], actin.cb_strengths[j]);
            return strength;
        }

        double _get_cb_strengths(int& i, int& j, double * binding_site_i, double * binding_site_j){
            double strength;
            double cos_angle =std::abs(std::cos(actin.thetas[i] - actin.thetas[j]));
            if (actin_basic_tensions[i]!=0 && actin_basic_tensions[j]!=0){
                int myosin_index_i = myosinIndicesPerActin.getConnections(i)[0];
                int myosin_index_j = myosinIndicesPerActin.getConnections(j)[0];

                // double dist1 = utils::segment_segment_distance(binding_site_i, binding_site_j,
                //  myosin.left_endpts[myosin_index_i], myosin.right_endpts[myosin_index_i], box);
                // double dist2 = utils::segment_segment_distance(binding_site_i, binding_site_j,
                //  myosin.left_endpts[myosin_index_j], myosin.right_endpts[myosin_index_j], box);
                // printf("dist1 %f, dist2 %f\n", dist1, dist2);
                // if (dist1<myosin.radius || dist2<myosin.radius){
                //     return 0;
                // }

                if (myosin_index_i != myosin_index_j){
                    strength = std::sqrt(actin_basic_tensions[i]*actin_basic_tensions[j])*cos_angle;
                    strength = std::min(std::min(strength, 1-actin.cb_strengths[i]), 1-actin.cb_strengths[j]);
                    if (strength>1){
                        printf("strength is greater than 1\n");
                    }
                }
            }
            else{
                strength = 1/cb_mult_factor*cos_angle;
                strength = std::min(std::min(strength, 1-actin.cb_strengths[i]), 1-actin.cb_strengths[j]);
                strength = std::max(strength, 0.0);
            }
            return strength;
        }

        void _search_cb(int&i, std::vector<int>& indices_to_search, std::vector<int>& indices_to_exclude, bool& add_connection)
        {
            std::vector<int> indices;
            std::vector<std::vector<double>> normal_vecs;
            std::vector<double> cb_strengths;
            for (int index = 0; index < indices_to_search.size(); index++){
                int j = indices_to_search[index];
                if (std::find(indices_to_exclude.begin(), indices_to_exclude.end(), j) != indices_to_exclude.end()){
                    continue;
                }
                if (actin_basic_tensions[j]!=0 && actin.cb_strengths[j]<1){
                    auto result = utils::segment_segment_distance_w_normal(actin.left_endpts[i], 
                        actin.right_endpts[i], actin.left_endpts[j], actin.right_endpts[j], box);
                    double distance = result.first;
                    double angle = actin.thetas[i] - actin.thetas[j];
                    double binding_site_i [2] = {result.second["start"][0], result.second["start"][1]};
                    double binding_site_j [2] = {result.second["start"][0]+result.second["vector"][0],
                         result.second["start"][1]+result.second["vector"][1]};                   
                    if (distance<2*alpha_actinin.radius && std::cos(angle)<0){
                        double strength = _get_cb_strengths(i,j, binding_site_i, binding_site_j);
                        if (strength>0.0){
                            cb_strengths.push_back(strength);
                            indices.push_back(j);
                            normal_vecs.push_back(result.second["vector"]);
                        }
                    }
                }
            }
            //sort by strength
            std::vector<size_t> sorted_indices = utils::sort_indices(cb_strengths);
            for (int index = 0; index < sorted_indices.size(); index++){
                int j = indices[sorted_indices[index]];
                if (actin.cb_strengths[i]>=1){
                    actin.cb_strengths[i] = 1;
                    break;
                }
                if (actin.cb_strengths[j]>=1){
                    actin.cb_strengths[j] = 1;
                    continue;
                }
                double cb_interaction = _set_cb_interactions(i,j,normal_vecs[sorted_indices[index]]);
                actin.cb_strengths[i] += cb_interaction;
                actin.cb_strengths[j] += cb_interaction;

                if (add_connection){
                    actinIndicesPerActin.addConnection(i,j);
                    actinIndicesPerActin.addConnection(j,i);
                }
            }
        }

        double _actin_myosin_energy(int& actin_index, int& myosin_index){
            int i = actin_index;
            int j = myosin_index;
            double center_displacement[2];
            center_displacement[0] = actin.xs[i][0] - myosin.xs[j][0];
            center_displacement[1] = actin.xs[i][1] - myosin.xs[j][1];
            utils::pbc_wrap(center_displacement, box);
            double center_distance = utils::norm(center_displacement, 2);
            if (center_distance>myosin.radius+0.5*(actin.length+myosin.length)){
                return 0;
            }
            else{
                auto result = utils::segment_segment_distance_w_normal(actin.left_endpts[i], 
                    actin.right_endpts[i], myosin.left_endpts[j], myosin.right_endpts[j], box);
                double distance = result.first;
                std::map<std::string, std::vector<double>> normal_vector_dict = result.second;
                std::vector<double> normal_vector = normal_vector_dict["vector"];
                std::vector<double> normal_start = normal_vector_dict["start"];
                a_m_normal_vec[i][j] = normal_vector;
                a_m_normal_start[i][j] = normal_start;
                double angle = std::abs(actin.thetas[i] - myosin.thetas[j]);
                angle = std::fmod(angle, 2*M_PI);
                angle = std::min(angle, 2*M_PI-angle);
                if (distance<myosin.radius){
                    double magnitude = std::abs(std::cos(angle));
                    return magnitude;
                }
                else{
                    return 0;
                }
            }
        }

        void _myosin_repulsion(int& i, int& j){
            double center_displacement[2];
            center_displacement[0] = myosin.xs[i][0] - myosin.xs[j][0];
            center_displacement[1] = myosin.xs[i][1] - myosin.xs[j][1];
            utils::pbc_wrap(center_displacement, box);
            double center_distance = utils::norm(center_displacement, 2);
            if (center_distance<=myosin.radius+myosin.length){
                auto result = utils::segment_segment_distance_w_normal(myosin.left_endpts[i], 
                    myosin.right_endpts[i], myosin.left_endpts[j], myosin.right_endpts[j], box);
                double distance = result.first;
                if (distance<2*myosin.radius){
                    std::vector<double> normal_vector = result.second["vector"];
                    if (normal_vector[0]==0 && normal_vector[1]==0){
                        normal_vector[0] = center_displacement[0]/center_distance;
                        normal_vector[1] = center_displacement[1]/center_distance;
                    }
                    else {
                        normal_vector[0] = -normal_vector[0]/distance;
                        normal_vector[1] = -normal_vector[1]/distance;
                    }
                    myosin.forces[i][0] += 2000*k_aa*normal_vector[0];
                    myosin.forces[i][1] += 2000*k_aa*normal_vector[1];
                    myosin.forces[j][0] += -2000*k_aa*normal_vector[0];
                    myosin.forces[j][1] += -2000*k_aa*normal_vector[1];
                }

            }
        }


        void new_file(){
            create_file(filename, actin.n, myosin.n, alpha_actinin.n);
        }

        void save_state(){
            append_to_file(filename, actin, myosin, alpha_actinin, total_energy);
        }
        
        void load_state(){
            load_from_file(filename, actin, myosin, alpha_actinin, total_energy);
            update_system();
        }
        
};

#endif
