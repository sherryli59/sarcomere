#ifndef SARCOMERE_H
#define SARCOMERE_H

#include <iostream>

#ifndef COMPONENTS_H
#include "components.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

#include "h5_utils.h"

#include <cmath>
#include <vector>

class Sarcomere
{
    public:
        Filament actin;
        Myosin myosin;
        AlphaActinin alpha_actinin;
        std::vector <double> box;
        double e_am, e_al, e_barrier, e_barrier_al, e_catch_bond, f_myosin;
        double total_energy;
        double ** a_m_energy;
        std::vector<double> ** a_m_force;
        double ** m_al_energy;
        double * al_energy; // actin-alpha actinin binding energy (including catch bonds) for each alpha actinin 
        std::vector<std::vector<std::vector<int>>> al_cb; // indices of actins that each alpha actinin forms a catch bond with
        //store the number of catch-bond-forming alpha-actinins each actin is bound to
        std::vector<std::vector<double>> al_cb_energy;
        int * actin_cb_count;


        std::string filename;

        Sarcomere(){
            box.resize(2);
            box[0] = 0;
            box[1] = 0;
            e_am = 0;
            e_al = 0;
            e_barrier = 0;
            e_barrier_al = 0;
            e_catch_bond = 0;
            f_myosin = 0;
            total_energy = 0;
        }
        Sarcomere(int& n_actins, int& n_myosins, int& n_alpha_actinins, std::vector<double> box0, double& actin_length, double& myosin_length,
            double& e_am, double& e_al, double& e_barrier, double& e_barrier_al, double& e_catch_bond,
            double& f_myosin, double& myosin_radius, double& alpha_actinin_radius, std::string& filename, gsl_rng * rng)
            : actin(n_actins, actin_length, box0, rng), 
                myosin(n_myosins, myosin_length, myosin_radius, box0, rng),
                alpha_actinin(n_alpha_actinins, alpha_actinin_radius,box0, rng){
            box.resize(2);
            box[0] = box0[0];
            box[1] = box0[1];
            this->e_am = e_am;
            this->e_al = e_al;
            this->e_barrier = e_barrier;
            this->e_barrier_al = e_barrier_al;
            this->e_catch_bond = e_catch_bond*(f_myosin+2.5);
            this->f_myosin = f_myosin;
            printf("e_am = %f, e_al = %f, e_barrier = %f, e_barrier_al = %f, e_catch_bond = %f, f_myosin = %f\n", e_am, e_al, e_barrier, e_barrier_al, e_catch_bond, f_myosin);
            a_m_energy = new double*[actin.n];
            a_m_force = new std::vector<double>*[actin.n];
            actin_cb_count = new int[actin.n];
            for (int i = 0; i < actin.n; i++){
                a_m_energy[i] = new double[myosin.n];
                a_m_force[i] = new std::vector<double>[myosin.n];
                for (int j = 0; j < myosin.n; j++){
                    a_m_force[i][j].resize(2);
                }
                actin_cb_count[i] = 0;
            }
            m_al_energy = new double*[myosin.n];
            for (int i = 0; i < myosin.n; i++){
                m_al_energy[i] = new double[alpha_actinin.n];
            }
            al_energy = new double[alpha_actinin.n];
            al_cb.resize(alpha_actinin.n);
            al_cb_energy.resize(alpha_actinin.n);
            for (int i = 0; i < alpha_actinin.n; i++){
                al_cb[i].resize(0);
                al_cb_energy[i].resize(0);
            }
            total_energy = get_energy();
            this->filename = filename;
        }
        ~Sarcomere(){ 
            printf("start destructor\n");
            for (int i = 0; i < actin.n; i++){
                delete[] a_m_energy[i];
                delete[] a_m_force[i];
            }
            delete[] a_m_energy;
            delete[] a_m_force;
            delete[] al_energy;
            for (int i = 0; i < myosin.n; i++){
                delete[] m_al_energy[i];
            }
            delete[] m_al_energy;
            printf("Sarcomere destructor called\n");
        }

        Sarcomere(const Sarcomere& other){
            actin = other.actin;
            myosin = other.myosin;
            alpha_actinin = other.alpha_actinin;
            box = other.box;
            e_am = other.e_am;
            e_al = other.e_al;
            e_barrier = other.e_barrier;
            e_barrier_al = other.e_barrier_al;
            e_catch_bond = other.e_catch_bond;
            f_myosin = other.f_myosin;
        }

        void sarcomeric_structure(){
            std::vector<std::vector<double>> myosin_positions;
            box[0] = 19.2;
            box[1] = 12;
            myosin_positions = {{-4.8,-3},{4.8,-3},{-4.8,0},{4.8,0}, {-4.8,3},{4.8,3}};
            for (int i = 0; i < myosin_positions.size(); i++){
                myosin.xs[i][0] = myosin_positions[i][0];
                myosin.xs[i][1] = myosin_positions[i][1];
                myosin.thetas[i] = 0;
            }
            myosin.n = myosin_positions.size();
            myosin.update_endpoints();
            std::vector<std::vector<double>> alpha_actinin_positions;
            alpha_actinin_positions = {{-9.6,-3.2},{0,-3.2},{-9.6,-3},{0,-3},{-9.6,-2.8},{0,-2.8},{0,-0.2},{-9.6,-0.2},{0,0},{-9.6,0},{-9.6,0.2},{0,0.2},{-9.6,2.8},{0,2.8},{-9.6,3},{0,3},{-9.6,3.2},{0,3.2}};
            for (int i = 0; i < alpha_actinin_positions.size(); i++){
                alpha_actinin.xs[i][0] = alpha_actinin_positions[i][0];
                alpha_actinin.xs[i][1] = alpha_actinin_positions[i][1];
            }
            alpha_actinin.n = alpha_actinin_positions.size();
            std::vector<std::vector<double>> actin_positions;
            actin_positions = {{-8.1,-3.2},{1.5,-3.2},{-8.1,-3},{1.5,-3},{-8.1,-2.8},{1.5,-2.8},{-8.1,-0.2},{1.5,-0.2},{-8.1,0},{1.5,0},{-8.1,0.2},{1.5,0.2},{-8.1,2.8},{1.5,2.8},{-8.1,3},{1.5,3},{-8.1,3.2},{1.5,3.2}};
            for (int i = 0; i < actin_positions.size(); i++){
                actin.xs[i][0] = actin_positions[i][0];
                actin.xs[i][1] = actin_positions[i][1];
                actin.thetas[i] = M_PI;
            }
            actin.n = actin_positions.size();
            actin_positions = {{-1.5,-3.2},{8.1,-3.2},{-1.5,-3},{8.1,-3},{-1.5,-2.8},{8.1,-2.8},{-1.5,-0.2},{8.1,-0.2},{-1.5,0},{8.1,0},{-1.5,0.2},{8.1,0.2},{-1.5,2.8},{8.1,2.8},{-1.5,3},{8.1,3},{-1.5,3.2},{8.1,3.2}};
            for (int i = 0; i < actin_positions.size(); i++){
                actin.xs[i+actin.n][0] = actin_positions[i][0];
                actin.xs[i+actin.n][1] = actin_positions[i][1];
                actin.thetas[i+actin.n] = 0;
            }
            actin.n = actin_positions.size()*2;
            actin.update_endpoints();
            total_energy = get_energy();
            printf("total energy: %f\n", total_energy);
            save_state();

        }

        void partial_sarcomeric_structure(){
            std::vector<std::vector<double>> myosin_positions;
            box[0] = 19.2;
            box[1] = 12;
            myosin_positions = {{-4.8,-3},{4.8,-3},{-4.8,0},{4.8,0}, {-4.8,3},{4.8,3}};
            for (int i = 0; i < myosin_positions.size(); i++){
                myosin.xs[i][0] = myosin_positions[i][0];
                myosin.xs[i][1] = myosin_positions[i][1];
                myosin.thetas[i] = 0;
            }
            myosin.n = myosin_positions.size();
            myosin.update_endpoints();
            std::vector<std::vector<double>> alpha_actinin_positions;
            alpha_actinin_positions = {{-9.6,-3.2},{0,-3.2},{-9.6,-3},{0,-3},{-9.6,-2.8},{0,-2.8},{0,-0.2},{-9.6,-0.2},{0,0},{-9.6,0},{-9.6,0.2},{0,0.2},{-9.6,2.8},{0,2.8},{-9.6,3},{0,3},{-9.6,3.2},{0,3.2}};
            for (int i = 0; i < alpha_actinin_positions.size(); i++){
                alpha_actinin.xs[i][0] = alpha_actinin_positions[i][0];
                alpha_actinin.xs[i][1] = alpha_actinin_positions[i][1];
            }
            alpha_actinin.n = alpha_actinin_positions.size();
            std::vector<std::vector<double>> actin_positions;
            actin_positions = {{1.5,-3.2},{1.5,-3},{1.5,-2.8},{1.5,-0.2},{1.5,0},{1.5,0.2},{1.5,2.8},{1.5,3},{1.5,3.2}};
            for (int i = 0; i < actin_positions.size(); i++){
                actin.xs[i][0] = actin_positions[i][0];
                actin.xs[i][1] = actin_positions[i][1];
                actin.thetas[i] = 0;
            }
            actin.n = actin_positions.size();
            actin_positions = {{-1.5,-3.2},{-1.5,-3},{-1.5,-2.8},{-1.5,-0.2},{-1.5,0},{-1.5,0.2},{-1.5,2.8},{-1.5,3},{-1.5,3.2}};
            for (int i = 0; i < actin_positions.size(); i++){
                actin.xs[i+actin.n][0] = actin_positions[i][0];
                actin.xs[i+actin.n][1] = actin_positions[i][1];
                actin.thetas[i+actin.n] = M_PI;
            }
            actin.n = actin_positions.size()*2;
            actin.update_endpoints();
            total_energy = get_energy();
            printf("total energy: %f\n", total_energy);
            save_state();

        }

        void myosins_on_a_lattice(){
            std::vector<std::vector<double>> myosin_positions;
            box[0] = 19.2;
            box[1] = 8;
            myosin_positions = {{-4,0},{4,0},{-4,3},{4,3},{-4,-3},{4,-3}};
            for (int i = 0; i < myosin_positions.size(); i++){
                myosin.xs[i][0] = myosin_positions[i][0];
                myosin.xs[i][1] = myosin_positions[i][1];
                myosin.thetas[i] = 0;
            }
            myosin.n = myosin_positions.size();
            myosin.update_endpoints();
            //half of the actins have theta = 0, the other half have theta = pi
            for (int i = 0; i < actin.n; i++){
                if (i < actin.n/2){
                    actin.thetas[i] = 0;
                } else {
                    actin.thetas[i] = M_PI;
                }
            }
            actin.update_endpoints();
            total_energy = get_energy();
            printf("total energy: %f\n", total_energy);
            save_state();

        }

        void equilibrate_myosin(int i, double x, double y, double theta, double& delta_energy){
            double repulsion_old = e_barrier*myosin.self_repulsion(i);
            myosin.displace(i, x, y, theta);
            double repulsion_new = e_barrier*myosin.self_repulsion(i);
            delta_energy = repulsion_new - repulsion_old;
        }

        void equilibrate_alpha_actinin(int i, double x, double y, double& delta_energy){
            delta_energy = 0;
            double repulsion_old = e_barrier_al*alpha_actinin.self_repulsion(i);
            alpha_actinin.displace(i, x, y);
            for (int j = 0; j < myosin.n; j++){
                double e_new = myosin_alpha_actinin_energy(j , i);
                delta_energy += e_new - m_al_energy[j][i];
                m_al_energy[j][i] = e_new;
            }
            double repulsion_new = e_barrier_al*alpha_actinin.self_repulsion(i);
            delta_energy += repulsion_new - repulsion_old;
        }

        double get_energy(){
            double energy = 0;
            for (int i = 0; i < actin.n; i++){
                for (int j = 0; j < myosin.n; j++){
                    double a_m = actin_myosin_energy(i , j);
                    energy += a_m;
                    a_m_energy[i][j] = a_m;
                }
            }
            bool update_energy = false; //a_m energy already updated
            update_force(update_energy);
            double * a_al_binding_e = new double[alpha_actinin.n];
            for (int j = 0; j < actin.n; j++){
                actin_cb_count[j] = 0;
            }
            for (int i = 0; i < alpha_actinin.n; i++){
                a_al_binding_e[i] = actin_alpha_actinin_energy(i);
            }
            for (int i = 0; i < alpha_actinin.n; i++){
                int n_cb = al_cb[i].size();
                double al_cb_e_i = 0;
                for (int j = 0; j < n_cb; j++){
                    double al_cb_e_ij = al_cb_energy[i][j];
                    std::vector<int> actin_indices = al_cb[i][j];
                    int rescaling_factor = actin_cb_count[actin_indices[0]]*actin_cb_count[actin_indices[1]];
                    al_cb_e_ij = al_cb_e_ij/rescaling_factor;
                    al_cb_e_i += al_cb_e_ij;
                    if (al_cb_e_i < e_catch_bond){
                        break;
                    }
                }
                al_cb_e_i = std::max(al_cb_e_i, e_catch_bond);
                energy +=  a_al_binding_e[i]+al_cb_e_i;
                al_energy[i] = a_al_binding_e[i]+al_cb_e_i;
            }
            for (int i = 0; i < myosin.n; i++){
                for (int j = 0; j < alpha_actinin.n; j++){
                    double m_al = myosin_alpha_actinin_energy(i , j);
                    energy += m_al;
                    m_al_energy[i][j] = m_al;
                }
            }
            delete[] a_al_binding_e;
            energy += e_barrier*myosin.total_self_repulsion();
            energy += e_barrier_al*alpha_actinin.total_self_repulsion();
            return energy;
        }
        
        void update_force(bool& update_energy){
            for (int i = 0; i<myosin.n; i++){
                myosin.forces[i][0] = 0;
                myosin.forces[i][1] = 0;
            }
            for (int i = 0; i < actin.n; i++){
                actin.forces[i][0] = 0;
                actin.forces[i][1] = 0;
                for (int j = 0; j < myosin.n; j++){
                    std::vector<double> force = myosin_force_on_actin(i, j, update_energy);
                    actin.forces[i][0] += force[0];
                    actin.forces[i][1] += force[1];
                    myosin.forces[j][0] -= force[0];
                    myosin.forces[j][1] -= force[1];
                }
            }
        }

        double actin_myosin_energy(int& actin_index, int& myosin_index){
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
                double distance = utils::segment_segment_distance(actin.left_endpts[i], 
                    actin.right_endpts[i], myosin.left_endpts[j], myosin.right_endpts[j], box);
                double angle = std::abs(actin.thetas[i] - myosin.thetas[j]);
                angle = std::fmod(angle, 2*M_PI);
                angle = std::min(angle, 2*M_PI-angle);
                if (distance<myosin.radius && (angle<M_PI/4||angle>3*M_PI/4)){
                    // std::cout<<i<<" "<<j<<std::endl;
                    // std::cout<<actin.left_endpts[i][0]<<" "<<actin.left_endpts[i][1]<<" "<<actin.right_endpts[i][0]<<" "<<actin.right_endpts[i][1]<<std::endl;
                    // std::cout<<myosin.left_endpts[j][0]<<" "<<myosin.left_endpts[j][1]<<" "<<myosin.right_endpts[j][0]<<" "<<myosin.right_endpts[j][1]<<std::endl;
                    double magnitude = std::abs(std::cos(angle));
                    return magnitude*e_am;
                }
                else{
                    return 0;
                }
            }
        }

        std::vector<double> myosin_force_on_actin(int& actin_index, int& myosin_index, bool& update_energy){
            std::vector<double> force(2);
            int i = actin_index;
            int j = myosin_index;
            double energy;
            if (update_energy){
                energy = actin_myosin_energy(i, j);
                a_m_energy[i][j] = energy;
                
            }
            else{
                energy = a_m_energy[i][j];
            }
            if (energy == 0){
                force[0] = 0;
                force[1] = 0;
            }
            else{
                double orientation[2];
                orientation[0] = cos(actin.thetas[i]);
                orientation[1] = sin(actin.thetas[i]);
                force[0] = f_myosin*energy/e_am*orientation[0];
                force[1] = f_myosin*energy/e_am*orientation[1];
            }
            a_m_force[i][j][0] = force[0];
            a_m_force[i][j][1] = force[1];
            return force;
        }

        double actin_alpha_actinin_binding_energy(int& actin_index, int& alpha_actinin_index){
            int i = actin_index;
            int j = alpha_actinin_index;
            double energy = 0;
            double distance = utils::point_segment_distance(alpha_actinin.xs[j], actin.left_endpts[i], actin.right_endpts[i], box);
            if (distance<alpha_actinin.radius){
                energy = e_al;
            }
            return energy;
        }


        double actin_alpha_actinin_cb_energy(int& actin_index, int& alpha_actinin_index){
            int i = actin_index;
            int j = alpha_actinin_index;
            if (e_catch_bond!=0 && abs(actin.forces[i][0])+abs(actin.forces[i][1])>0){
                double catch_bond_strength = 0;
                for (int k=0; k<i; k++){ //actin k<i to avoid double counting
                    if (abs(actin.forces[k][0])+abs(actin.forces[k][1])>0){
                        //check if actin k is bound to alpha-actinin j
                        double kj_distance = utils::point_segment_distance(alpha_actinin.xs[j], 
                        actin.left_endpts[k], actin.right_endpts[k], box);
                        if (kj_distance<alpha_actinin.radius){
                            double force_mag_i = std::sqrt(actin.forces[i][0]*actin.forces[i][0]+actin.forces[i][1]*actin.forces[i][1]);
                            double force_mag_k = std::sqrt(actin.forces[k][0]*actin.forces[k][0]+actin.forces[k][1]*actin.forces[k][1]);
                            double pairwise_cosine = (actin.forces[i][0]*actin.forces[k][0]+actin.forces[i][1]*actin.forces[k][1])/(force_mag_i*force_mag_k);
                            double scaling_factor = std::min(force_mag_i/f_myosin,1.0)*std::min(force_mag_k/f_myosin,1.0);
                            pairwise_cosine = pairwise_cosine*scaling_factor;
                            pairwise_cosine = std::min(pairwise_cosine, 0.0);
                            if (abs(pairwise_cosine)>1e-4){
                                bool bound_to_different_myosins = false;
                                for (int l=0; l<myosin.n; l++){
                                    if ((a_m_energy[k][l]!=0 && a_m_energy[i][l]==0)||(a_m_energy[k][l]==0 && a_m_energy[i][l]!=0))
                                    {
                                        bound_to_different_myosins = true;
                                        break;
                                    }
                                }
                                if (bound_to_different_myosins){
                                    //check how many alpha-actinins are bound to actin k and actin i
                                    
                                    catch_bond_strength+=abs(pairwise_cosine);
                                    al_cb_energy[j].push_back(e_catch_bond*abs(pairwise_cosine));
                                    std::vector<int> actin_index = {k, i};
                                    al_cb[j].push_back(actin_index);
                                    actin_cb_count[k]++;
                                    actin_cb_count[i]++;
                                    // printf("pairwise_cosine = %f\n", pairwise_cosine);
                                    printf("actin %d and %d are bound to alpha-actinin %d\n", i, k, j);
                                    // std::cout<<alpha_actinin.xs[j][0]<<" "<<alpha_actinin.xs[j][1]<<std::endl;
                                    // std::cout<<actin.left_endpts[k][0]<<" "<<actin.left_endpts[k][1]<<std::endl;
                                    // std::cout<<actin.right_endpts[k][0]<<" "<<actin.right_endpts[k][1]<<std::endl;
                                    // std::cout<<actin.left_endpts[i][0]<<" "<<actin.left_endpts[i][1]<<std::endl;
                                    // std::cout<<actin.right_endpts[i][0]<<" "<<actin.right_endpts[i][1]<<std::endl;
                                    if (catch_bond_strength>1){
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }

                double catch_bond_energy = e_catch_bond*std::min(catch_bond_strength, 1.0);
                return catch_bond_energy;
            }
            else{
                return 0;
            }
        }

        double actin_alpha_actinin_energy(int& alpha_actinin_index){
            // The sun of all actin-alpha-actinin interactions for a given alpha-actinin
            int i = alpha_actinin_index;
            double binding_e = 0;
            double cb_e = 0;
            double binding_ej, cb_ej;
            al_cb[i].clear();
            for (int j = 0; j < actin.n; j++){
                binding_ej = actin_alpha_actinin_binding_energy(j , i);
                if (binding_ej!=0){
                    cb_ej = actin_alpha_actinin_cb_energy(j , i);
                }
                else{
                    cb_ej = 0;
                }
                binding_e += binding_ej;
                cb_e += cb_ej;
            }
            binding_e = std::max(binding_e, 2*e_al);
            cb_e = std::max(cb_e, e_catch_bond);
            return binding_e;
        }

        
        double myosin_alpha_actinin_energy(int& myosin_index, int& alpha_actinin_index){
            int i = myosin_index;
            int j = alpha_actinin_index;
            double distance = utils::point_segment_distance(alpha_actinin.xs[j], myosin.left_endpts[i], myosin.right_endpts[i], box);
            //printf("Distance between myosin %d and alpha-actinin %d: %f\n", i, j, distance);
            if (distance<alpha_actinin.radius+myosin.radius){
                return e_barrier_al;
            }
            else{
                return 0;
            }
        }

        void displace_actin(int i, double x, double y, double theta, double& delta_energy, std::vector<double>& avg_force){
            actin.displace(i, x, y, theta);
            delta_energy = 0;
            avg_force.resize(2);
            std::vector <double> force(2); 
            force[0] = 0;
            force[1] = 0;
            bool update_energy = false;
            for (int j = 0; j < myosin.n; j++){
                double e_new = actin_myosin_energy(i , j);
                delta_energy += e_new - a_m_energy[i][j];
                a_m_energy[i][j] = e_new;
                std::vector<double> f_old = a_m_force[i][j];
                std::vector<double> f = myosin_force_on_actin(i, j, update_energy);
                myosin.forces[j][0] += -(f[0] - f_old[0]);
                myosin.forces[j][1] += -(f[1] - f_old[1]);
                force[0] += f[0];
                force[1] += f[1];
                
            }
            avg_force[0] = (actin.forces[i][0]+force[0])/2;
            avg_force[1] = (actin.forces[i][1]+force[1])/2;
            actin.forces[i][0] = force[0];
            actin.forces[i][1] = force[1];
            double * a_al_binding_e = new double[alpha_actinin.n];
            for (int j = 0; j < actin.n; j++){
                actin_cb_count[j] = 0;
            }
            for (int j = 0; j < alpha_actinin.n; j++){
                a_al_binding_e[j] = actin_alpha_actinin_energy(j);
            }
            for (int j = 0; j < alpha_actinin.n; j++){
                int n_cb = al_cb[j].size();
                double al_cb_e_j = 0;
                for (int k = 0; k < n_cb; k++){
                    double al_cb_e_jk = al_cb_energy[j][k];
                    std::vector<int> actin_indices = al_cb[j][k];
                    int rescaling_factor = actin_cb_count[actin_indices[0]]*actin_cb_count[actin_indices[1]];
                    if (al_cb_e_jk!=0){
                        if (rescaling_factor == 0){
                            printf("Error: rescaling factor is zero\n");
                        }
                        al_cb_e_jk = al_cb_e_jk/rescaling_factor;
                        al_cb_e_j += al_cb_e_jk;
                        if (al_cb_e_j < e_catch_bond){
                            break;
                        }
                    }
                   
                }
                al_cb_e_j = std::max(al_cb_e_j, e_catch_bond);
                delta_energy +=  a_al_binding_e[j]+al_cb_e_j- al_energy[j];
                al_energy[j] = a_al_binding_e[j]+al_cb_e_j;
            }
            total_energy += delta_energy;
            delete[] a_al_binding_e;
        }

        void displace_myosin(int i, double x, double y, double theta, double& delta_energy, std::vector<double>& avg_force){
            delta_energy = 0;
            //MYOSIN DISPLACEMENT IS RELATIVELY EXPENSIVE RIGHT NOW
            double repulsion_old = e_barrier*myosin.self_repulsion(i);
            myosin.displace(i, x, y, theta);
            delta_energy = 0;
            std::vector <double> force(2); 
            force[0] = 0;
            force[1] = 0;
            bool update_energy = false;
            for (int j = 0; j < actin.n; j++){
                double e_new = actin_myosin_energy(j, i);
                delta_energy += e_new - a_m_energy[j][i];
                a_m_energy[j][i] = e_new;
                std::vector<double> f_old = a_m_force[j][i];
                std::vector<double> f = myosin_force_on_actin(j, i, update_energy);
                actin.forces[j][0] += f[0] - f_old[0];
                actin.forces[j][1] += f[1] - f_old[1];
                force[0] -= f[0];
                force[1] -= f[1];
            }
            avg_force[0] = (myosin.forces[i][0]+force[0])/2;
            avg_force[1] = (myosin.forces[i][1]+force[1])/2;
            myosin.forces[i][0] = force[0];
            myosin.forces[i][1] = force[1];
            if (e_catch_bond!=0){
                //need to update actin_alpha-actinin energy as well to account for changes in catch bonds
                double * a_al_binding_e = new double[alpha_actinin.n];
                for (int j = 0; j < actin.n; j++){
                    actin_cb_count[j] = 0;
                }
                for (int j = 0; j < alpha_actinin.n; j++){
                    a_al_binding_e[j] = actin_alpha_actinin_energy(j);
                }
                for (int j = 0; j < alpha_actinin.n; j++){
                    int n_cb = al_cb[j].size();
                    double al_cb_e_j = 0;
                    for (int k = 0; k < n_cb; k++){
                        double al_cb_e_jk = al_cb_energy[j][k];
                        std::vector<int> actin_indices = al_cb[j][k];
                        int rescaling_factor = actin_cb_count[actin_indices[0]]*actin_cb_count[actin_indices[1]];
                        al_cb_e_jk = al_cb_e_jk/rescaling_factor;
                        al_cb_e_j += al_cb_e_jk;
                        if (al_cb_e_j < e_catch_bond){
                            break;
                        }
                    }
                    al_cb_e_j = std::max(al_cb_e_j, e_catch_bond);
                    delta_energy +=  a_al_binding_e[j]+al_cb_e_j- al_energy[j];
                    al_energy[j] = a_al_binding_e[j]+al_cb_e_j;
                }
                delete[] a_al_binding_e;

            }

            for (int j = 0; j < alpha_actinin.n; j++){
                double e_new = myosin_alpha_actinin_energy(i , j);
                delta_energy += e_new - m_al_energy[i][j];
                m_al_energy[i][j] = e_new;
            }
            double repulsion_new = e_barrier*myosin.self_repulsion(i);
            delta_energy += repulsion_new - repulsion_old;
            total_energy += delta_energy;
        }

        void displace_alpha_actinin(int i, double x, double y, double& delta_energy){
            delta_energy = 0;
            double repulsion_old = e_barrier_al*alpha_actinin.self_repulsion(i);
            alpha_actinin.displace(i, x, y);
            // double a_al_e = actin_alpha_actinin_energy(i);
            // delta_energy +=  a_al_e- al_energy[i];
            // al_energy[i] = a_al_e;
            double * a_al_binding_e = new double[alpha_actinin.n];
            for (int j = 0; j < actin.n; j++){
                actin_cb_count[j] = 0;
            }
            for (int j = 0; j < alpha_actinin.n; j++){
                a_al_binding_e[j] = actin_alpha_actinin_energy(j);
            }
            for (int j = 0; j < alpha_actinin.n; j++){
                int n_cb = al_cb[j].size();
                double al_cb_e_j = 0;
                for (int k = 0; k < n_cb; k++){
                    double al_cb_e_jk = al_cb_energy[j][k];
                    std::vector<int> actin_indices = al_cb[j][k];
                    int rescaling_factor = actin_cb_count[actin_indices[0]]*actin_cb_count[actin_indices[1]];
                    al_cb_e_jk = al_cb_e_jk/rescaling_factor;
                    al_cb_e_j += al_cb_e_jk;
                    if (al_cb_e_j < e_catch_bond){
                        break;
                    }
                }
                al_cb_e_j = std::max(al_cb_e_j, e_catch_bond);
                delta_energy +=  a_al_binding_e[j]+al_cb_e_j- al_energy[j];
                al_energy[j] = a_al_binding_e[j]+al_cb_e_j;
            }
            for (int j = 0; j < myosin.n; j++){
                double e_new = myosin_alpha_actinin_energy(j , i);
                delta_energy += e_new - m_al_energy[j][i];
                m_al_energy[j][i] = e_new;
            }
            double repulsion_new = e_barrier_al*alpha_actinin.self_repulsion(i);
            delta_energy += repulsion_new - repulsion_old;
            total_energy += delta_energy;  
            delete[] a_al_binding_e;

        }

        void new_file(){
            create_file(filename, actin.n, myosin.n, alpha_actinin.n);
        }

        void save_state(){
            append_to_file(filename, actin, myosin, alpha_actinin, total_energy);
            // for (int i=0; i<myosin.n; i++){
            //     printf("Myosin %d: %f %f\n", i, myosin.xs[i][0], myosin.xs[i][1]);
            //     printf("Myosin endpoints: (%f %f), (%f %f)\n", myosin.left_endpts[i][0], myosin.left_endpts[i][1], myosin.right_endpts[i][0], myosin.right_endpts[i][1]);
            //     }
            //     //print positions of alpha-actinins
            // for (int i=0; i<alpha_actinin.n; i++){
            //     printf("Alpha-actinin %d: %f %f\n", i, alpha_actinin.xs[i][0], alpha_actinin.xs[i][1]);
            // }
            // for (int i=0; i<actin.n; i++){
            //     printf("Actin %d: %f %f\n", i, actin.xs[i][0], actin.xs[i][1]);
            //     printf("Actin endpoints: (%f %f), (%f %f)\n", actin.left_endpts[i][0], actin.left_endpts[i][1], actin.right_endpts[i][0], actin.right_endpts[i][1]);
            // }

        }
        
        void load_state(){
            load_from_file(filename, actin, myosin, alpha_actinin, total_energy);
            total_energy = get_energy();
            printf("Loaded state with energy %f\n", total_energy);
        }
        
        void check_energy(){
            double total_energy2 = get_energy();
            if (abs(total_energy2-total_energy)>1e-6){
                std::cout << "Total energy is not correct" << std::endl;
                std::cout << "Total energy: " <<total_energy << " Total energy2: " << total_energy2 << std::endl;
                exit(1);
            }
        }
        void check_forces(){
            bool update_energy = false; //calculating previous force
           
            for (int i = 0; i < actin.n; i++){
                std::vector<double> f = {0, 0};
                for (int j=0; j<myosin.n; j++){
                    std::vector<double> fij =  myosin_force_on_actin(i, j, update_energy);
                    f[0] += fij[0];
                    f[1] += fij[1];
                }
                //check if force is correct
                if (abs(f[0]-actin.forces[i][0])>1e-6 || abs(f[1]-actin.forces[i][1])>1e-6){
                    std::cout << "Force on actin " << i << " is not correct" << std::endl;
                    std::cout << "Force: " << f[0] << " " << f[1] << " Force2: " << actin.forces[i][0] << " " << actin.forces[i][1] << std::endl;
                    exit(1);
                }
            }
        }
};

#endif
