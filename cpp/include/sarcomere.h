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


class Sarcomere
{
    public:
        Filament actin;
        Myosin myosin;
        AlphaActinin alpha_actinin;
        double * box;
        double e_am, e_al, e_barrier, f_myosin;
        bool catch_bond;
        double total_energy;
        double ** a_m_energy;
        double ** a_al_energy;
        double ** m_al_energy;
        std::string filename;

        Sarcomere(){
            box = new double[2];
            box[0] = 0;
            box[1] = 0;
            e_am = 0;
            e_al = 0;
            e_barrier = 0;
            f_myosin = 0;
            catch_bond = false;
            total_energy = 0;
        }
        Sarcomere(int& n_actins, int& n_myosins, int& n_alpha_actinins, double* box0, double& actin_length, double& myosin_length,
            double& e_am, double& e_al, double& f_myosin, double& myosin_radius, double& alpha_actinin_radius, bool& catch_bond, gsl_rng * rng)
            : actin(n_actins, actin_length, box0, rng), 
                myosin(n_myosins, myosin_length, myosin_radius, box0, rng),
                alpha_actinin(n_alpha_actinins, alpha_actinin_radius,box0, rng){
            box = new double[2];
            box[0] = box0[0];
            box[1] = box0[1];
            this->e_am = e_am;
            this->e_al = e_al;
            e_barrier = 100.0;
            this->f_myosin = f_myosin;
            this->catch_bond = catch_bond;
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
            filename = "traj.h5";
        }
        ~Sarcomere(){ 
            delete[] box;  
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
        }

        Sarcomere(const Sarcomere& other){
            actin = other.actin;
            myosin = other.myosin;
            alpha_actinin = other.alpha_actinin;
            box = new double[2];
            box[0] = other.box[0];
            box[1] = other.box[1];
            e_am = other.e_am;
            e_al = other.e_al;
            f_myosin = other.f_myosin;
            catch_bond = other.catch_bond;
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
            energy += e_barrier/100*alpha_actinin.total_self_repulsion();
            return energy;
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
        double actin_alpha_actinin_energy(int& actin_index, int& alpha_actinin_index){
            int i = actin_index;
            int j = alpha_actinin_index;
            double distance = utils::point_segment_distance(alpha_actinin.xs[j], actin.left_endpts[i], actin.right_endpts[i], box);
            if (distance<alpha_actinin.radius){
                return e_al;
            }
            else{
                return 0;
            }
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

        void displace_actin(int i, double x, double y, double theta, double& delta_energy){
            actin.displace(i, x, y, theta);
            delta_energy = 0;
            double total1 = 0;
            double total2 = 0;
            for (int j = 0; j < myosin.n; j++){
                double e_new = actin_myosin_energy(i , j);
                delta_energy += e_new - a_m_energy[i][j];
                total1 += e_new;
                total2 += a_m_energy[i][j];
                a_m_energy[i][j] = e_new;
            }
            for (int j = 0; j < alpha_actinin.n; j++){
                double e_new = actin_alpha_actinin_energy(i , j);
                delta_energy += e_new - a_al_energy[i][j];
                total1 += e_new;
                total2 += a_al_energy[i][j];
                a_al_energy[i][j] = e_new;
            }
            total_energy += delta_energy;
           
        }

        void displace_myosin(int i, double x, double y, double theta, double& delta_energy){
            double repulsion_old = e_barrier*myosin.self_repulsion(i);
            myosin.displace(i, x, y, theta);
            delta_energy = 0;
            for (int j = 0; j < actin.n; j++){
                double e_new = actin_myosin_energy(j, i);
                delta_energy += e_new - a_m_energy[j][i];
                a_m_energy[j][i] = e_new;
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
            double repulsion_old = e_barrier/100*alpha_actinin.self_repulsion(i);
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
            double repulsion_new = e_barrier/100*alpha_actinin.self_repulsion(i);
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
            printf("actin_xs[0][0] = %f\n", actin.xs[0][0]);
            load_from_file(filename, actin, myosin, alpha_actinin, total_energy);
            total_energy = get_energy();
            printf("Loaded state with energy %f\n", total_energy);
        }
        
};

#endif
