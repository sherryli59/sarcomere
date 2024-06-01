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
        double ** a_al_energy;
        double ** m_al_energy;
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
            this->e_catch_bond = e_catch_bond;
            this->f_myosin = f_myosin;
            printf("e_am = %f, e_al = %f, e_barrier = %f, e_barrier_al = %f, e_catch_bond = %f, f_myosin = %f\n", e_am, e_al, e_barrier, e_barrier_al, e_catch_bond, f_myosin);
            a_m_energy = new double*[actin.n];
            for (int i = 0; i < actin.n; i++){
                a_m_energy[i] = new double[myosin.n];
            }
            a_al_energy = new double*[actin.n];
            for (int i = 0; i < actin.n; i++){
                a_al_energy[i] = new double[alpha_actinin.n];
            }
            m_al_energy = new double*[myosin.n];
            for (int i = 0; i < myosin.n; i++){
                m_al_energy[i] = new double[alpha_actinin.n];
            }
            total_energy = get_energy();
            this->filename = filename;
        }
        ~Sarcomere(){ 
            printf("start destructor\n");
            for (int i = 0; i < actin.n; i++){
                delete[] a_m_energy[i];
                delete[] a_al_energy[i];
            }
            delete[] a_m_energy;
            delete[] a_al_energy;
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

        void initialize_myosins(){
            std::vector<std::vector<double>> myosin_positions;
            myosin_positions = {{0,0}, {3,0}, {0,3}, {3,3}};
            for (int i = 0; i < myosin_positions.size(); i++){
                myosin.xs[i][0] = myosin_positions[i][0];
                myosin.xs[i][1] = myosin_positions[i][1];
                myosin.thetas[i] = 0;
            }
            myosin.n = myosin_positions.size();
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
            for (int i = 0; i < actin.n; i++){
                for (int j = 0; j < alpha_actinin.n; j++){
                    double a_al = actin_alpha_actinin_energy(i , j);
                    energy += a_al;
                    a_al_energy[i][j] = a_al;
                }
            }
            
            for (int i = 0; i < myosin.n; i++){
                for (int j = 0; j < alpha_actinin.n; j++){
                    double m_al = myosin_alpha_actinin_energy(i , j);
                    energy += m_al;
                    m_al_energy[i][j] = m_al;
                }
            }
            energy += e_barrier*myosin.total_self_repulsion();
            energy += e_barrier_al*alpha_actinin.total_self_repulsion();
            return energy;
        }
        
        void update_force(bool& update_energy){
            for (int i = 0; i < actin.n; i++){
                actin.forces[i][0] = 0;
                actin.forces[i][1] = 0;
                for (int j = 0; j < myosin.n; j++){
                    std::vector<double> force = myosin_force_on_actin(i, j, update_energy);
                    actin.forces[i][0] += force[0];
                    actin.forces[i][1] += force[1];
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
                if (distance<myosin.radius){
                    return e_am;
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
                force[0] = f_myosin*orientation[0];
                force[1] = f_myosin*orientation[1];
            }
            return force;
        }

        double actin_alpha_actinin_energy(int& actin_index, int& alpha_actinin_index){
            // Assume all the actin-myosin interactions are up to date
            const int i = actin_index;
            const int j = alpha_actinin_index;
            double energy = 0;
            double distance = utils::point_segment_distance(alpha_actinin.xs[j], actin.left_endpts[i], actin.right_endpts[i], box);
            if (distance<alpha_actinin.radius){
                energy = e_al;
            }
            if (e_catch_bond!=0 && energy!=0){
                if (abs(actin.forces[i][0])+abs(actin.forces[i][1])>0){
                    double catch_bond_strength = 0;
                    for (int k=0; k<actin.n; k++){
                        if (k!=i){
                            
                            if (abs(actin.forces[k][0])+abs(actin.forces[k][1])>0){
                                //check if actin k is bound to alpha-actinin j
                                double kj_distance = utils::point_segment_distance(alpha_actinin.xs[j], 
                                actin.left_endpts[k], actin.right_endpts[k], box);
                                if (kj_distance<alpha_actinin.radius){
                                    double pairwise_cosine = actin.forces[i][0]*actin.forces[k][0]+actin.forces[i][1]*actin.forces[k][1];
                                    pairwise_cosine = pairwise_cosine/(f_myosin*f_myosin);
                                    pairwise_cosine = std::min(pairwise_cosine, 0.0);
                                    catch_bond_strength+=abs(pairwise_cosine);
                                    //print xs of actin i, k, and alpha-actinin j
                                    //if (abs(pairwise_cosine)>0){
                                    // printf("i=%d,j=%d,k=%d\n", i, j, k);
                                    // printf("catch_bond_strength=%f\n", abs(pairwise_cosine));
                                    // printf("Actin (%f, %f), Actin (%f, %f), Alpha-actinin (%f, %f)\n", actin.xs[i][0], actin.xs[i][1], actin.xs[k][0], actin.xs[k][1], alpha_actinin.xs[j][0], alpha_actinin.xs[j][1]);
                                    // }
                                    if (catch_bond_strength>2){
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    double catch_bond_energy = e_catch_bond*std::min(catch_bond_strength, 2.0);
                    // if (catch_bond_energy!=0){
                    //     printf("Catch bond energy: %f\n", catch_bond_energy);
                    // }

                    energy = energy+catch_bond_energy;
                }
            }
            // if (distance<alpha_actinin.radius){
            //     printf("Actin %d is bound to alpha-actinin %d\n", i, j);
            //}
            return energy;
        }
        
        double myosin_alpha_actinin_energy(int& myosin_index, int& alpha_actinin_index){
            int i = myosin_index;
            int j = alpha_actinin_index;
            double distance = utils::point_segment_distance(alpha_actinin.xs[j], myosin.left_endpts[i], myosin.right_endpts[i], box);
            //printf("Distance: %f\n", distance);
            if (distance<alpha_actinin.radius+myosin.radius){
                return e_barrier;
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
                std::vector<double> f = myosin_force_on_actin(i, j, update_energy);
                force[0] += f[0];
                force[1] += f[1];
            }
            avg_force[0] = (actin.forces[i][0]+force[0])/2;
            avg_force[1] = (actin.forces[i][1]+force[1])/2;
            actin.forces[i][0] = force[0];
            actin.forces[i][1] = force[1];
            if (e_catch_bond!=0){
                //TEMPORARY FIX: update all actin-myosin interactions; should only update affected actins
                for (int j = 0; j < actin.n; j++){
                    for (int k = 0; k < alpha_actinin.n; k++){
                        double e_new = actin_alpha_actinin_energy(j , k);
                        delta_energy += e_new - a_al_energy[j][k];
                        a_al_energy[j][k] = e_new;
                    }
                }
            }
            else{
                //THIS DOESN'T WORK WHEN CATCH BOND IS ON; NEED TO UPDATE OTHER ACTINS AS WELL
                for (int j = 0; j < alpha_actinin.n; j++){
                    double e_new = actin_alpha_actinin_energy(i , j);
                    delta_energy += e_new - a_al_energy[i][j];
                    a_al_energy[i][j] = e_new;
                }
            };
            total_energy += delta_energy;
        }

        void displace_myosin(int i, double x, double y, double theta, double& delta_energy){
            //MYOSIN DISPLACEMENT IS RELATIVELY EXPENSIVE RIGHT NOW
            double repulsion_old = e_barrier*myosin.self_repulsion(i);
            bool update_energy = false; //calculating previous force
            double ** forces_old = new double*[actin.n];
            for (int j = 0; j < actin.n; j++){
                forces_old[j] = new double[2];
                std::vector<double> f = myosin_force_on_actin(j, i, update_energy);
                forces_old[j][0] = f[0];
                forces_old[j][1] = f[1];
            }
            myosin.displace(i, x, y, theta);
            delta_energy = 0;
            for (int j = 0; j < actin.n; j++){
                double e_new = actin_myosin_energy(j, i);
                delta_energy += e_new - a_m_energy[j][i];
                a_m_energy[j][i] = e_new;
                std::vector<double> f = myosin_force_on_actin(j, i, update_energy);
                actin.forces[j][0] = actin.forces[j][0] - forces_old[j][0] + f[0];
                actin.forces[j][1] = actin.forces[j][1] - forces_old[j][1] + f[1];
            }

            if (e_catch_bond!=0){
                //need to update actin_alpha-actinin energy as well to account for changes in catch bonds
                for (int j = 0; j < actin.n; j++){
                    for (int k = 0; k < alpha_actinin.n; k++){
                        double e_new = actin_alpha_actinin_energy(j, k);
                        delta_energy += e_new - a_al_energy[j][k];
                        a_al_energy[j][k] = e_new;
                    }
                }
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
            double repulsion_old = e_barrier_al*alpha_actinin.self_repulsion(i);
            alpha_actinin.displace(i, x, y);
            delta_energy = 0;
            for (int j = 0; j < actin.n; j++){
                double e_new = actin_alpha_actinin_energy(j , i);
                delta_energy += e_new - a_al_energy[j][i];
                a_al_energy[j][i] = e_new;
            }
            for (int j = 0; j < myosin.n; j++){
                double e_new = myosin_alpha_actinin_energy(j , i);
                delta_energy += e_new - m_al_energy[j][i];
                m_al_energy[j][i] = e_new;
            }
            double repulsion_new = e_barrier_al*alpha_actinin.self_repulsion(i);
            delta_energy += repulsion_new - repulsion_old;
            total_energy += delta_energy;    
        }

        void new_file(){
            create_file(filename, actin.n, myosin.n, alpha_actinin.n);
        }

        void save_state(){
            append_to_file(filename, actin, myosin, alpha_actinin, total_energy);
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
                if (f[0] != actin.forces[i][0] || f[1] != actin.forces[i][1]){
                    std::cout << "Force on actin " << i << " is not correct" << std::endl;
                    std::cout << "Force: " << f[0] << " " << f[1] << " Force2: " << actin.forces[i][0] << " " << actin.forces[i][1] << std::endl;
                }
            }
        }
};

#endif
