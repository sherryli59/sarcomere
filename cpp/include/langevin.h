#ifndef Langevin_H
#define Langevin_H

#ifndef SARCOMERE_H
#include "sarcomere.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

#ifndef H5_UTILS_H
#include "h5_utils.h"
#endif

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iostream>




class Langevin{
    public: 
        Langevin(Sarcomere& model0, double& beta0, double& dt0, double& D0, int& update_myosin_every0,
            int& update_dt_every0, int& save_every0, bool& resume):
            model(model0), beta(beta0), dt(dt0), D(D0), update_myosin_every(update_myosin_every0), update_dt_every(update_dt_every0), save_every(save_every0){
            acc_rate=0;
            if (resume){
                model.load_state();
                printf("Resuming from file %s\n", model.filename.c_str());
                //printf("Total energy: %f\n", model.get_energy());
                //print positions of myosins
                for (int i=0; i<model.myosin.n; i++){
                    printf("Myosin %d: %f %f\n", i, model.myosin.center[i].x, model.myosin.center[i].y);
                    printf("Myosin endpoints: (%f %f), (%f %f)\n", model.myosin.left_end[i].x, model.myosin.left_end[i].y, model.myosin.right_end[i].x, model.myosin.right_end[i].y);
                }
                //print positions of actins
                for (int i=0; i<model.actin.n; i++){
                    printf("Actin %d: %f %f\n", i, model.actin.center[i].x, model.actin.center[i].y);
                    printf("Actin endpoints: (%f %f), (%f %f)\n", model.actin.left_end[i].x, model.actin.left_end[i].y, model.actin.right_end[i].x, model.actin.right_end[i].y);
                }
            }
            else{
                model.new_file();
            }
        }
        double dt, beta, acc_rate, D;
        int update_myosin_every, update_dt_every, save_every;
        Sarcomere& model;
        ~Langevin(){
        }

      
        
        void run_langevin(int nsteps, gsl_rng * rng, int& fix_myosin){
            for (int i=0; i<nsteps; i++){
                //bool update_myosin = (update_myosin_every>0 && i%update_myosin_every==0);
                if (i%save_every==0){
                    std::cout << "Step " << i << std::endl;
                    model.save_state();
                    //model.check_energy();
                }
                sample_step(dt, D, rng, fix_myosin);
            }
        }
        void sample_step(double& dt, double& D, gsl_rng * rng, int& fix_myosin){
            model.update_system();
            //double gamma = 1/(beta*D);
            //get noise
            int n_randns = (model.myosin.n+model.actin.n)*3;
            double noise[n_randns];
            for (int i=0; i<n_randns; i++)
                noise[i]=gsl_rng_uniform(rng)*2-1;
            int n_acc_randns = model.myosin.n+model.actin.n;
            double acc_rand[n_acc_randns];
            for (int i=0; i<n_acc_randns; i++)
                acc_rand[i]=gsl_rng_uniform(rng);
            double dx, dy, dtheta;
            std::vector <double> force(2);
            for (int i=fix_myosin; i<model.myosin.n; i++){
                dx = model.myosin.force[i].x*beta*D*dt+model.myosin.velocity[i].x*dt+sqrt(2*D*dt)*noise[i*3];
                dy = model.myosin.force[i].y*beta*D*dt+model.myosin.velocity[i].y*dt+sqrt(2*D*dt)*noise[i*3+1];
                // dx = model.myosin.force[i].x*beta*D*dt+sqrt(2*D*dt)*noise[i*3];
                // dy = model.myosin.force[i].y*beta*D*dt+sqrt(2*D*dt)*noise[i*3+1];
                dtheta = model.myosin.angular_force[i]*beta*D*dt+sqrt(2*D*dt)*noise[i*3+2]*M_PI/5;
                model.myosin.displace(i, dx, dy, dtheta);
            }
            int offset = model.myosin.n*3;
            for (int i=0; i<model.actin.n; i++){
                // printf("actin %d\n", i);
                // printf("force*beta*D*dt: %f %f\n", model.actin.force[i].x*beta*D*dt, model.actin.force[i].y*beta*D*dt);
                // printf("velocity*dt: %f %f\n", model.actin.velocity[i].x*dt, model.actin.velocity[i].y*dt);
                // printf("sqrt(2*D*dt)*noise: %f %f\n", sqrt(2*D*dt)*noise[offset+i*3], sqrt(2*D*dt)*noise[offset+i*3+1]);
                dx = model.actin.force[i].x*beta*D*dt+model.actin.velocity[i].y*dt+sqrt(2*D*dt)*noise[offset+i*3];
                dy = model.actin.force[i].y*beta*D*dt+model.actin.velocity[i].y*dt+sqrt(2*D*dt)*noise[offset+i*3+1];
                double noise_x = sqrt(2*D*dt)*noise[offset+i*3];
                double noise_y = sqrt(2*D*dt)*noise[offset+i*3+1];
                double velocity_x = model.actin.velocity[i].x*dt;
                double velocity_y = model.actin.velocity[i].y*dt;
                double force_x = model.actin.force[i].x*beta*D*dt;
                double force_y = model.actin.force[i].y*beta*D*dt;
                dtheta = model.actin.angular_force[i]*beta*D*dt+sqrt(2*D*dt)*noise[offset+i*3+2]*M_PI/5;
                double noise_theta = sqrt(2*D*dt)*noise[offset+i*3+2]*M_PI/5;
                double force_theta = model.actin.angular_force[i]*beta*D*dt;
                // printf("angular force: %f\n", model.actin.angular_force[i]);
                // printf("dtheta: %f\n", dtheta);
                
                if (dx>0.04||dy>0.04){
                    printf("actin %d\n", i);
                    printf("force*beta*D*dt: %f %f\n", model.actin.force[i].x*beta*D*dt, model.actin.force[i].y*beta*D*dt);
                    printf("velocity*dt: %f %f\n", model.actin.velocity[i].x*dt, model.actin.velocity[i].y*dt);
                    printf("sqrt(2*D*dt)*noise: %f %f\n", sqrt(2*D*dt)*noise[offset+i*3], sqrt(2*D*dt)*noise[offset+i*3+1]);
                    printf("dx, dy: %f %f\n", dx, dy);
                    printf("dtheta: %f\n", dtheta);
                }
                model.actin.displace(i, dx, dy, dtheta);
            }
        }
};

#endif

