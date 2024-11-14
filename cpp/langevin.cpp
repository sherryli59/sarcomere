#include "include/langevin.h"
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
    double diff_coeff = 0.01;
    int update_dt_every = 500;
    int update_myosin_every = 1;
    int save_every = 500;
    double cb_mult_factor = 200;
    double k_aa =0.4;
    double kappa_aa = 0.8;
    double k_am =0.3;
    double kappa_am = 0.8;
    double v_am =0.5;
    int n_actins = 50;
    int n_myosins = 4;
    double Lx = 10;
    double Ly = 10;
    double actin_length = 3;
    double myosin_length = 3;
    double myosin_radius = 0.5;
    double crosslinker_length = 0.3;
    bool resume = false;
    int fix_myosin = 0;
    std::string filename = "data/traj.h5";
    std::string init_struc = "random";
    try {
        cxxopts::Options options("sarcomere", "simulate sarcomere assembly using Monte Carlo");

        options.add_options()
            ("nsteps", "Number of steps", cxxopts::value<int>(nsteps)->default_value("500000"))
            ("seed", "Seed value", cxxopts::value<int>(seed)->default_value("0"))
            ("dt", "Time step", cxxopts::value<double>(dt)->default_value("0.01"))
            ("beta", "Beta value", cxxopts::value<double>(beta)->default_value("1.0"))
            ("diff_coeff", "Diffusion coefficient", cxxopts::value<double>(diff_coeff)->default_value("0.01"))
            ("update_dt_every", "Update dt every", cxxopts::value<int>(update_dt_every)->default_value("500"))
            ("update_myosin_every", "Update myosin every", cxxopts::value<int>(update_myosin_every)->default_value("1"))
            ("save_every", "Save every", cxxopts::value<int>(save_every)->default_value("500"))
            ("cb_mult_factor", "Catch bond multiplier factor", cxxopts::value<double>(cb_mult_factor)->default_value("200"))
            ("k_aa", "k_aa", cxxopts::value<double>(k_aa)->default_value("0.2"))
            ("kappa_aa", "kappa_aa", cxxopts::value<double>(kappa_aa)->default_value("0.6"))
            ("k_am", "k_am", cxxopts::value<double>(k_am)->default_value("0.2"))
            ("kappa_am", "kappa_am", cxxopts::value<double>(kappa_am)->default_value("0.6"))
            ("v_am", "v_am", cxxopts::value<double>(v_am)->default_value("0.5"))
            ("n_actins", "Number of actins", cxxopts::value<int>(n_actins)->default_value("50"))
            ("n_myosins", "Number of myosins", cxxopts::value<int>(n_myosins)->default_value("4"))
            ("Lx", "Lx", cxxopts::value<double>(Lx)->default_value("10.0"))
            ("Ly", "Ly", cxxopts::value<double>(Ly)->default_value("10.0"))
            ("actin_length", "Actin length", cxxopts::value<double>(actin_length)->default_value("3.0"))
            ("myosin_length", "Myosin length", cxxopts::value<double>(myosin_length)->default_value("3.0"))
            ("myosin_radius", "Myosin radius", cxxopts::value<double>(myosin_radius)->default_value("0.5"))
            ("crosslinker_length", "Crosslinker length", cxxopts::value<double>(crosslinker_length)->default_value("0.3"))
            ("resume", "Resume", cxxopts::value<bool>(resume)->default_value("false"))
            ("fix_myosin", "Fix myosin", cxxopts::value<int>(fix_myosin)->default_value("0"))
            ("filename", "Filename", cxxopts::value<std::string>(filename)->default_value("data/traj.h5"))
            ("initial_structure", "Type of initial structure", cxxopts::value<std::string>(init_struc)->default_value("random"))
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
    Sarcomere model(n_actins, n_myosins, box,
                      actin_length, myosin_length, cb_mult_factor,  k_aa, kappa_aa, k_am, kappa_am, v_am,
                       myosin_radius, crosslinker_length,filename,rng);
    Langevin sim(model, beta, dt, diff_coeff, update_myosin_every, update_dt_every, save_every, resume);
    int myosin_nsteps = 1000;
    int aa_nsteps = 1000;
    double dt0 = 0.5; 
    if (!resume){
        if (init_struc == "partial_sarcomere") {
        sim.model.partial_sarcomeric_structure();}
        //  else if (init_struc == "catch_bond") {
        //     sim.model.catch_bond();
        // }
        // else if (init_struc == "sarcomere"){
        //     sim.model.sarcomeric_structure();
        // }
        else if (init_struc == "fix_myosin"){
            sim.model.myosin_on_a_lattice();
        }
        else if (init_struc == "partial"){
            sim.model.partial_fix();
        }
    }
    sim.run_langevin(nsteps, rng, fix_myosin);
    gsl_rng_free(rng);
    return 0;
}
