#include "include/mc.h"
#include "include/sarcomere.h"
#include "include/components.h"
#include "include/cxxopts.hpp"

int main(int argc, char* argv[]){
    int nsteps = 1000000;
    int seed = 0;
    double dt = 0.01;
    double beta = 1;
    double D = 0.1;
    int update_dt_every = 100;
    int save_every = 100;
    double e_al = -1;
    double e_am = -3;
    double f_myosin = 5;
    int n_actins = 50;
    int n_myosins = 4;
    int n_alpha_actinins = 50;
    double Lx0 = 10;
    double Ly0 = 10;
    double actin_length = 3;
    double myosin_length = 3;
    double myosin_radius = 0.5;
    double alpha_actinin_radius = 0.2;
    bool catch_bond = false;
    bool resume = false;
    try {
        cxxopts::Options options("sarcomere", "simulate sarcomere assembly using Monte Carlo");

        options.add_options()
            ("nsteps", "Number of steps", cxxopts::value<int>(nsteps)->default_value("100000"))
            ("seed", "Seed value", cxxopts::value<int>(seed)->default_value("0"))
            ("dt", "Time step", cxxopts::value<double>(dt)->default_value("0.01"))
            ("beta", "Beta value", cxxopts::value<double>(beta)->default_value("1.0"))
            ("D", "Diffusion coefficient", cxxopts::value<double>(D)->default_value("0.1"))
            ("update_dt_every", "Update dt every", cxxopts::value<int>(update_dt_every)->default_value("100"))
            ("save_every", "Save every", cxxopts::value<int>(save_every)->default_value("100"))
            ("e_al", "Energy AL", cxxopts::value<double>(e_al)->default_value("-1.0"))
            ("e_am", "Energy AM", cxxopts::value<double>(e_am)->default_value("-3.0"))
            ("f_myosin", "Force Myosin", cxxopts::value<double>(f_myosin)->default_value("5.0"))
            ("n_actins", "Number of actins", cxxopts::value<int>(n_actins)->default_value("30"))
            ("n_myosins", "Number of myosins", cxxopts::value<int>(n_myosins)->default_value("2"))
            ("n_alpha_actinins", "Number of alpha actinins", cxxopts::value<int>(n_alpha_actinins)->default_value("30"))
            ("Lx0", "Lx0", cxxopts::value<double>(Lx0)->default_value("10.0"))
            ("Ly0", "Ly0", cxxopts::value<double>(Ly0)->default_value("4.0"))
            ("actin_length", "Actin length", cxxopts::value<double>(actin_length)->default_value("3.0"))
            ("myosin_length", "Myosin length", cxxopts::value<double>(myosin_length)->default_value("3.0"))
            ("myosin_radius", "Myosin radius", cxxopts::value<double>(myosin_radius)->default_value("0.5"))
            ("alpha_actinin_radius", "Alpha actinin radius", cxxopts::value<double>(alpha_actinin_radius)->default_value("0.2"))
            ("catch_bond", "Catch bond", cxxopts::value<bool>(catch_bond)->default_value("false"))
            ("resume", "Resume", cxxopts::value<bool>(resume)->default_value("false"))
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
    double * box = new double[2];
    box[0] = Lx0;
    box[1] = Ly0;
    //Filament actin(n_actins, actin_length, box, rng);

    Sarcomere model(n_actins, n_myosins, n_alpha_actinins, box,
                      actin_length, myosin_length, e_am, e_al, f_myosin, myosin_radius, alpha_actinin_radius, catch_bond, rng);
    MC mc(model, beta, dt, D, update_dt_every, save_every, resume);
    mc.run_mc(nsteps, rng);
    return 0;
}
