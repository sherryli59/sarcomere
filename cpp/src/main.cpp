#include "langevin.h"
#include "sarcomere.h"
#include "components.h"
#include "cxxopts.hpp"
#include <array>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

int main(int argc, char* argv[]){
    int nsteps;
    int seed;
    double dt;
    double beta;
    double actin_diff_coeff_trans;
    double actin_diff_coeff_rot;
    double max_actin_displacement;
    double max_myosin_displacement;
    double myosin_diff_coeff_trans;
    double myosin_diff_coeff_rot;
    double max_displacement;
    double max_actin_rotation;
    double max_myosin_rotation;
    int save_every;
    double k_on;
    double k_off;
    double base_lifetime;
    double lifetime_coeff;
    double k_aa;
    double kappa_aa;
    double k_am;
    double kappa_am;
    double k_mm;
    double v_am;
    double wall_k;
    double wall_dcut;
    double wall_exponent;
    std::string periodic_mask_str;
    int n_actins;
    int n_myosins;
    double Lx;
    double Ly;
    double Lz;
    double actin_length;
    double myosin_length;
    double myosin_radius;
    double am_cutoff;
    double am_optimal;
    double aa_cutoff;
    double aa_optimal;
    double tau_rec;
    double titin_k;
    bool resume;
    int resume_frame;
    bool directional;
    bool deterministic;
    bool use_autodiff_forces;
    int n_fixed_myosins;
    int max_myosin_bonds;
    int dimension;

    std::string filename;
    std::string init_struc;
    try {
        cxxopts::Options options("sarcomere", "simulate sarcomere assembly using Monte Carlo");

        options.add_options()
            ("nsteps", "Number of steps", cxxopts::value<int>(nsteps)->default_value("500000"))
            ("seed", "Seed value", cxxopts::value<int>(seed)->default_value("0"))
            ("dt", "Time step", cxxopts::value<double>(dt)->default_value("0.000001"))
            ("beta", "Beta value", cxxopts::value<double>(beta)->default_value("241.0"))
            ("actin_diff_coeff_trans", "Actin translational diffusion coefficient", cxxopts::value<double>(actin_diff_coeff_trans)->default_value("0.2"))
            ("actin_diff_coeff_rot", "Actin rotational diffusion coefficient", cxxopts::value<double>(actin_diff_coeff_rot)->default_value("0.1"))
            ("myosin_diff_coeff_trans", "Myosin translational diffusion coefficient", cxxopts::value<double>(myosin_diff_coeff_trans)->default_value("0.05"))
            ("myosin_diff_coeff_rot", "Myosin rotational diffusion coefficient", cxxopts::value<double>(myosin_diff_coeff_rot)->default_value("0.05"))
            ("max_actin_displacement", "Maximum deterministic displacement per step from forces for actin",
             cxxopts::value<double>(max_actin_displacement)->default_value("0.005"))
            ("max_myosin_displacement", "Maximum deterministic displacement per step from forces for myosin",
             cxxopts::value<double>(max_myosin_displacement)->default_value("0.005"))
            ("max_actin_rotation", "Maximum rotational magnitude per step for actin",
             cxxopts::value<double>(max_actin_rotation)->default_value("0.005"))
            ("max_myosin_rotation", "Maximum rotational magnitude per step for myosin",
             cxxopts::value<double>(max_myosin_rotation)->default_value("0.005"))
            ("save_every", "Save every", cxxopts::value<int>(save_every)->default_value("200"))
            ("k_on", "k_on", cxxopts::value<double>(k_on)->default_value("5000"))
            ("k_off", "k_off", cxxopts::value<double>(k_off)->default_value("1"))
            ("base_lifetime", "Base lifetime", cxxopts::value<double>(base_lifetime)->default_value("0.001"))
            ("lifetime_coeff", "Lifetime coefficient", cxxopts::value<double>(lifetime_coeff)->default_value("0.4"))
            ("k_aa", "k_aa", cxxopts::value<double>(k_aa)->default_value("300"))
            ("kappa_aa", "kappa_aa", cxxopts::value<double>(kappa_aa)->default_value("100"))
            ("k_am", "k_am", cxxopts::value<double>(k_am)->default_value("300"))
            ("kappa_am", "kappa_am", cxxopts::value<double>(kappa_am)->default_value("100"))
            ("k_mm", "Myomesin spring constant", cxxopts::value<double>(k_mm)->default_value("0.0"))
            ("v_am", "v_am", cxxopts::value<double>(v_am)->default_value("5"))
            ("periodic_mask", "Periodicity mask as three characters (e.g. 110 => periodic in x,y)",
             cxxopts::value<std::string>(periodic_mask_str)->default_value("111"))
            ("wall_k", "Wall spring constant", cxxopts::value<double>(wall_k)->default_value("0.0"))
            ("wall_dcut", "Wall interaction cutoff distance", cxxopts::value<double>(wall_dcut)->default_value("0.0"))
            ("wall_exponent", "Wall force exponent", cxxopts::value<double>(wall_exponent)->default_value("2.0"))
            ("use_autodiff_forces", "Use autodiff-based AA/AM force calculations (slower)", cxxopts::value<bool>(use_autodiff_forces)->default_value("false"))
            ("n_actins", "Number of actins", cxxopts::value<int>(n_actins)->default_value("50"))
            ("n_myosins", "Number of myosins", cxxopts::value<int>(n_myosins)->default_value("4"))
            ("Lx", "Lx", cxxopts::value<double>(Lx)->default_value("10"))
            ("Ly", "Ly", cxxopts::value<double>(Ly)->default_value("10"))
            ("Lz", "Lz", cxxopts::value<double>(Lz)->default_value("10"))
            ("actin_length", "Actin length", cxxopts::value<double>(actin_length)->default_value("1"))
            ("myosin_length", "Myosin length", cxxopts::value<double>(myosin_length)->default_value("1.5"))
            ("myosin_radius", "Myosin radius", cxxopts::value<double>(myosin_radius)->default_value("0.025"))
            ("am_cutoff", "cutoff for am interaction range", cxxopts::value<double>(am_cutoff)->default_value("0.05"))
            ("am_optimal", "optimal distance for am interaction", cxxopts::value<double>(am_optimal)->default_value("0.03"))
            // ("aa_cutoff", "cutoff for aa interaction range", cxxopts::value<double>(aa_cutoff)->default_value("0.05"))
            // ("aa_optimal", "optimal distance for aa interaction", cxxopts::value<double>(aa_optimal)->default_value("0.03"))
            ("tau_rec", "Cooldown time after KMC break", cxxopts::value<double>(tau_rec)->default_value("0.002"))
            ("titin_k", "Titin effective spring constant", cxxopts::value<double>(titin_k)->default_value("0.0"))
            ("resume", "Resume", cxxopts::value<bool>(resume)->default_value("false"))
            ("resume_frame", "Frame index to load when resuming (0-based; default loads latest)",
             cxxopts::value<int>(resume_frame)->default_value("-1"))
            ("directional", "Directional", cxxopts::value<bool>(directional)->default_value("true"))
            ("deterministic", "Enable reproducible deterministic scheduling", cxxopts::value<bool>(deterministic)->default_value("false"))
            ("n_fixed_myosins", "Number of fixed myosins", cxxopts::value<int>(n_fixed_myosins)->default_value("0"))
            ("filename", "Filename", cxxopts::value<std::string>(filename)->default_value("data/traj.h5"))
            ("initial_structure", "Type of initial structure", cxxopts::value<std::string>(init_struc)->default_value("random"))
            ("max_myosin_bonds", "Maximum actin bonds per myosin",cxxopts::value<int>(max_myosin_bonds)->default_value("6"))
            ("dimension", "Simulation dimensionality (2 or 3)",
             cxxopts::value<int>(dimension)->default_value("3"))
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


    if (dimension != 2 && dimension != 3) {
        std::cerr << "dimension must be 2 or 3" << std::endl;
        return 1;
    }
    bool is3D = (dimension == 3);

    if (!resume && resume_frame >= 0) {
        std::cout << "Warning: --resume_frame ignored because --resume was not set.\n";
        resume_frame = -1;
    }
    if (resume_frame < -1) {
        std::cout << "Warning: --resume_frame must be -1 or non-negative. Using latest frame instead.\n";
        resume_frame = -1;
    }

#ifdef _OPENMP
    if (deterministic) {
        omp_set_dynamic(0);
        omp_set_schedule(omp_sched_static, 0);
    } else {
        omp_set_schedule(omp_sched_dynamic, 0);
    }
#endif

    gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);
    std::vector <double> box(3);
    box[0] = Lx;
    box[1] = Ly;
    box[2] = Lz;
    double diff_coeff_ratio = actin_diff_coeff_trans/myosin_diff_coeff_trans;
    double titin_rest_length = 0.5 * myosin_length + (2.0 / 3.0) * actin_length;
    auto compute_max_force = [&](double max_disp, double diff_coeff) {
        if (max_disp <= 0.0 || diff_coeff <= 0.0 || !std::isfinite(beta) || !std::isfinite(dt)) {
            return std::numeric_limits<double>::infinity();
        }
        double val = max_disp / (beta * diff_coeff * dt);
        if (!(val > 0.0) || !std::isfinite(val)) {
            return std::numeric_limits<double>::infinity();
        }
        return val;
    };
    auto compute_max_torque = [&](double max_rotation, double diff_coeff_rot) {
        if (max_rotation <= 0.0 || diff_coeff_rot <= 0.0 || !std::isfinite(beta) || !std::isfinite(dt)) {
            printf("max_rotation = %f, diff_coeff_rot = %f, beta = %f, dt = %f\n", max_rotation, diff_coeff_rot, beta, dt);
            return std::numeric_limits<double>::infinity();
        }
        double val = max_rotation / (beta * diff_coeff_rot * dt);
        if (!(val > 0.0) || !std::isfinite(val)) {
            return std::numeric_limits<double>::infinity();
        }
        return val;
    };
    double max_actin_force = compute_max_force(max_actin_displacement, actin_diff_coeff_trans);
    double max_myosin_force = compute_max_force(max_myosin_displacement, myosin_diff_coeff_trans);
    double max_actin_torque = compute_max_torque(max_actin_rotation, actin_diff_coeff_rot);
    double max_myosin_torque = compute_max_torque(max_myosin_rotation, myosin_diff_coeff_rot);
    printf("Max actin force: %f, max myosin force: %f\n", max_actin_force, max_myosin_force);
    if (std::isfinite(max_actin_torque)) {
        printf("Max actin torque: %f\n", max_actin_torque);
    } else {
        printf("Max actin torque: inf\n");
    }
    if (std::isfinite(max_myosin_torque)) {
        printf("Max myosin torque: %f\n", max_myosin_torque);
    } else {
        printf("Max myosin torque: inf\n");
    }
    aa_cutoff = am_cutoff;
    aa_optimal = am_optimal;
    auto parse_periodic_mask = [](const std::string& mask) {
        std::array<bool,3> periodic{true, true, true};
        if (!mask.empty()) {
            if (mask.size() != 3) {
                throw std::invalid_argument("periodic_mask must have length 3 (e.g. 110)");
            }
            for (size_t i = 0; i < 3; ++i) {
                char c = mask[i];
                if (c == '1' || c == 'T' || c == 't' || c == 'y' || c == 'Y') {
                    periodic[i] = true;
                } else if (c == '0' || c == 'F' || c == 'f' || c == 'n' || c == 'N') {
                    periodic[i] = false;
                } else {
                    throw std::invalid_argument("periodic_mask characters must be 0/1 (got '" + std::string(1, c) + "')");
                }
            }
        }
        return periodic;
    };

    std::array<bool,3> periodic_axes = parse_periodic_mask(periodic_mask_str);

    Sarcomere model(n_actins, n_myosins, box, actin_length, myosin_length,
                        myosin_radius, am_cutoff, am_optimal, aa_cutoff, aa_optimal,
                        k_on, k_off,
                        base_lifetime, lifetime_coeff, diff_coeff_ratio,
                          k_aa, kappa_aa, k_am, kappa_am, k_mm, v_am,
                        filename,rng, seed, n_fixed_myosins, dt, tau_rec,
                        titin_k, titin_rest_length,
                        directional, max_myosin_bonds, max_actin_force, max_myosin_force,
                        max_actin_torque, max_myosin_torque, periodic_axes, use_autodiff_forces);
    model.set_periodicity(periodic_axes);
    model.set_wall_parameters(wall_k, wall_dcut, wall_exponent);
    if (!is3D) {
        for (int i = 0; i < n_actins; ++i) {
            model.actin.center[i].z = 0;
            model.actin.direction[i].z = 0;
            model.actin.force[i].z = 0;
            model.actin.velocity[i].z = 0;
            model.actin.torque[i].z = 0;
        }
        for (int i = 0; i < n_myosins; ++i) {
            model.myosin.center[i].z = 0;
            model.myosin.direction[i].z = 0;
            model.myosin.force[i].z = 0;
            model.myosin.velocity[i].z = 0;
            model.myosin.torque[i].z = 0;
        }
        model.actin.update_endpoints();
        model.myosin.update_endpoints();
    }

    Langevin sim(model, beta, dt, actin_diff_coeff_trans,actin_diff_coeff_rot, myosin_diff_coeff_trans, myosin_diff_coeff_rot, save_every, resume, is3D,
                    max_actin_displacement, max_myosin_displacement,
                    max_actin_rotation, max_myosin_rotation, resume_frame);
    if (!resume){
        if (init_struc == "sarcomere") {
        sim.model.sarcomeric_structure_tight();}

        else if (init_struc == "partial"){
            sim.model.partial_fix(n_fixed_myosins);
        }
        else if (init_struc == "cb"){
            sim.model.cb();
        }
        else if (init_struc == "cb_off_angle"){
            sim.model.cb_off_angle();
        }
        int n_volume_exclusion = 0;
        sim.volume_exclusion(n_volume_exclusion, rng, n_fixed_myosins);
    }
    sim.run_langevin(nsteps, rng, n_fixed_myosins);
    gsl_rng_free(rng);
    return 0;
}
