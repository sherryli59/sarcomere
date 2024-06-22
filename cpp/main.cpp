#include "include/mc.h"
#include "include/sarcomere.h"
#include "include/components.h"
#include "include/cxxopts.hpp"
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

int main(int argc, char* argv[]){
    int nsteps = 500000;
    int seed = 0;
    double dt = 0.01;
    double beta = 1;
    double D = 0.1;
    int update_dt_every = 500;
    int save_every = 500;
    double e_al = -1;
    double e_am = -3;
    double e_barrier = 100;
    double e_barrier_al = 3;
    double e_catch_bond = -2;
    double f_myosin = 500;
    int n_actins = 50;
    int n_myosins = 4;
    int n_alpha_actinins = 50;
    double Lx = 10;
    double Ly = 10;
    double actin_length = 3;
    double myosin_length = 3;
    double myosin_radius = 0.5;
    double alpha_actinin_radius = 0.1;
    bool resume = false;
    std::string filename = "data/traj.h5";
    try {
        cxxopts::Options options("sarcomere", "simulate sarcomere assembly using Monte Carlo");

        options.add_options()
            ("nsteps", "Number of steps", cxxopts::value<int>(nsteps)->default_value("500000"))
            ("seed", "Seed value", cxxopts::value<int>(seed)->default_value("0"))
            ("dt", "Time step", cxxopts::value<double>(dt)->default_value("0.01"))
            ("beta", "Beta value", cxxopts::value<double>(beta)->default_value("1.0"))
            ("D", "Diffusion coefficient", cxxopts::value<double>(D)->default_value("0.1"))
            ("update_dt_every", "Update dt every", cxxopts::value<int>(update_dt_every)->default_value("500"))
            ("save_every", "Save every", cxxopts::value<int>(save_every)->default_value("500"))
            ("e_al", "Energy AL", cxxopts::value<double>(e_al)->default_value("-1.0"))
            ("e_am", "Energy AM", cxxopts::value<double>(e_am)->default_value("-3.0"))
            ("e_barrier", "Energy barrier", cxxopts::value<double>(e_barrier)->default_value("100.0"))
            ("e_barrier_al", "Energy barrier AL", cxxopts::value<double>(e_barrier_al)->default_value("3.0"))
            ("e_catch_bond", "Energy catch bond", cxxopts::value<double>(e_catch_bond)->default_value("-2"))
            ("f_myosin", "Force Myosin", cxxopts::value<double>(f_myosin)->default_value("500.0"))
            ("n_actins", "Number of actins", cxxopts::value<int>(n_actins)->default_value("50"))
            ("n_myosins", "Number of myosins", cxxopts::value<int>(n_myosins)->default_value("4"))
            ("n_alpha_actinins", "Number of alpha actinins", cxxopts::value<int>(n_alpha_actinins)->default_value("50"))
            ("Lx", "Lx", cxxopts::value<double>(Lx)->default_value("10.0"))
            ("Ly", "Ly", cxxopts::value<double>(Ly)->default_value("10.0"))
            ("actin_length", "Actin length", cxxopts::value<double>(actin_length)->default_value("3.0"))
            ("myosin_length", "Myosin length", cxxopts::value<double>(myosin_length)->default_value("3.0"))
            ("myosin_radius", "Myosin radius", cxxopts::value<double>(myosin_radius)->default_value("0.5"))
            ("alpha_actinin_radius", "Alpha actinin radius", cxxopts::value<double>(alpha_actinin_radius)->default_value("0.1"))
            ("resume", "Resume", cxxopts::value<bool>(resume)->default_value("false"))
            ("filename", "Filename", cxxopts::value<std::string>(filename)->default_value("data/traj.h5"))
            ("h, help", "Print usage");

        auto result = options.parse(argc, argv);

        if (result.count("help")) {
            std::cout << options.help() << std::endl;
            return 0;
        }
    } catch (const std::exception& e)
	  {
		    std::cerr << "Error parsing options: " << e.what() << std::endl;
		    return 1;
	  }


    gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);
    std::vector <double> box(2);
    box[0] = Lx;
    box[1] = Ly;
    //Filament actin(n_actins, actin_length, box, rng);

    Sarcomere model(n_actins, n_myosins, n_alpha_actinins, box,
                      actin_length, myosin_length, e_am, e_al, e_barrier, e_barrier_al, e_catch_bond,
                       f_myosin, myosin_radius, alpha_actinin_radius,filename,rng);
    MC mc(model, beta, dt, D, update_dt_every, save_every, resume);
    int myosin_nsteps = 1000000;
    int aa_nsteps = 1000000;
    double dt0 = 0.5;
    mc.equilibrate(myosin_nsteps,aa_nsteps,dt0, rng);
    mc.run_mc(nsteps, rng);
    gsl_rng_free(rng);
    return 0;
}
