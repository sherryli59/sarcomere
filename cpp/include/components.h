#ifndef COMPONENTS_H
#define COMPONENTS_H

#ifndef UTILS_H
#include "utils.h"
#endif

#include <cmath>
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


class Filament
{
    public:
        int n;
        double * box;
        double length;
        double ** xs;
        double * thetas;
        double ** left_endpts;
        double ** right_endpts;
        double ** force;

        Filament(){
            n = 0;
            length = 0;
            xs = NULL;
            thetas = NULL;

        }
        Filament(int n0, double length0, double* box0, gsl_rng * rng){
            n = n0;
            length = length0;
            xs = new double*[n];
            thetas = new double[n];
            force = new double*[n];
            left_endpts = new double*[n];
            right_endpts = new double*[n];
            box = box0;
            for (int i = 0; i < n; i++){
                xs[i] = new double[2];
                xs[i][0] = gsl_ran_flat(rng, 0, box[0]);
                xs[i][1] = gsl_ran_flat(rng, 0, box[1]);
                thetas[i] = gsl_ran_flat(rng, 0, 2*M_PI);
                left_endpts[i] = new double[2];
                right_endpts[i] = new double[2];
                force[i] = new double[2];
                force[i][0] = 0;
                force[i][1] = 0;
            }
            update_endpoints();
        }

        virtual ~Filament(){
            for (int i = 0; i < n; i++){
                delete[] xs[i];
                delete[] force[i];
                delete[] left_endpts[i];
                delete[] right_endpts[i];

            }
            delete[] thetas;
            printf("Filament destructor called\n");
        }
        Filament(const Filament& other){
            n = other.n;
            box = other.box;
            length = other.length;
            xs = new double*[n];
            thetas = new double[n];
            for (int i = 0; i < n; i++){
                xs[i] = new double[2];
                xs[i][0] = other.xs[i][0];
                xs[i][1] = other.xs[i][1];
                thetas[i] = other.thetas[i];
                force[i] = new double[2];
                force[i][0] = other.force[i][0];
                force[i][1] = other.force[i][1];
            }
        }

        void displace(int& i, double& dx, double& dy){
            xs[i][0] += dx;
            xs[i][1] += dx;
            utils::pbc_wrap(xs[i], box);
            update_endpoints_i(i);
        }

        void displace(int& i, double& dx, double& dy, double& dtheta){
            xs[i][0] += dx;
            xs[i][1] += dy;
            thetas[i] += dtheta;
            utils::pbc_wrap(xs[i], box);
            utils::angle_wrap(thetas[i]);
            update_endpoints_i(i);
        }

        void update_endpoints_i(int& i){
            double * segments = new double[2];
            segments[0] = length*cos(thetas[i]);
            segments[1] = length*sin(thetas[i]);
            left_endpts[i][0] = xs[i][0] - 0.5*segments[0];
            left_endpts[i][1] = xs[i][1] - 0.5*segments[1];
            right_endpts[i][0] = xs[i][0] + 0.5*segments[0];
            right_endpts[i][1] = xs[i][1] + 0.5*segments[1];
            delete[] segments;
        }
        
        void update_endpoints(){
            for (int i = 0; i < n; i++){
                update_endpoints_i(i);
            }
        }
        void update_thetas(double * new_thetas){
            for (int i = 0; i < n; i++){
                thetas[i] = new_thetas[i];
            }
            update_endpoints();
        }
        void update_xs(double ** new_xs){
            for (int i = 0; i < n; i++){
                xs[i][0] = new_xs[i][0];
                xs[i][1] = new_xs[i][1];
            }
            update_endpoints();
        }
};


class Myosin: public Filament
{
    public:
        double radius;
        
        Myosin(){
            Filament();
        }

        Myosin(int n0, double length0, double radius0,  double* box0, gsl_rng * rng):
            Filament(n0, length0, box0,  rng)
            {
                radius = radius0;

            }   

        Myosin(const Myosin& other):
            Filament(other)
            {
                radius = other.radius;
            }

        double total_self_repulsion(){
            double e = 0;
            for (int i = 0; i < n; i++){
                for (int j=0; j < i; j++){
                    double r = utils::segment_segment_distance(left_endpts[i], right_endpts[i], left_endpts[j], right_endpts[j], box);
                    if (r < 2*radius){
                        e += 1;
                    }
                }
            }
            return e;
        }
        double self_repulsion(int& i){
            double e = 0;
            for (int j = 0; j < n; j++){
                if (j != i){
                    double r = utils::segment_segment_distance(left_endpts[i], right_endpts[i], left_endpts[j], right_endpts[j], box);
                    if (r < 2*radius){
                        e += 1;
                    }
                }
            }
            return e;
        }
};

class AlphaActinin
{
    public:
        int n;
        double radius;
        double ** xs;
        double * box;
        AlphaActinin(){
            n = 0;
            radius = 0;
            xs = NULL;
        }
        AlphaActinin(int& n0, double& radius0, double * box0, gsl_rng * rng){
            n = n0;
            radius = radius0;
            xs = new double*[n];
            box = box0;
            for (int i = 0; i < n; i++){
                xs[i] = new double[2];
                xs[i][0] = gsl_ran_flat(rng, 0, box[0]);
                xs[i][1] = gsl_ran_flat(rng, 0, box[1]);
            }
        }

        ~AlphaActinin(){
            for (int i = 0; i < n; i++){
                delete[] xs[i];
            }
            delete[] xs;
        }
        AlphaActinin(const AlphaActinin& other){
            n = other.n;
            radius = other.radius;
            xs = new double*[n];
            for (int i = 0; i < n; i++){
                xs[i] = new double[2];
                xs[i][0] = other.xs[i][0];
                xs[i][1] = other.xs[i][1];
            }
        }
        double total_self_repulsion(){
            double e = 0;
            for (int i = 0; i < n; i++){
                for (int j=0; j < i; j++){
                    double r = utils::point_point_distance(xs[i], xs[j], box);
                    if (r < 2*radius){
                        e += 1;
                    }
                }
            }
            return e;
        }
        double self_repulsion(int& i){
            double e = 0;
            for (int j = 0; j < n; j++){
                if (j != i){
                    double r = utils::point_point_distance(xs[i], xs[j], box);
                    if (r < 2*radius){
                        e += 1;
                    }
                }
            }
            return e;
        }
        void displace(int& i, double& dx, double& dy){
            xs[i][0] += dx;
            xs[i][1] += dx;
            utils::pbc_wrap(xs[i], box);
        }
};

#endif
