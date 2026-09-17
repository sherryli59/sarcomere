#include "rod_mobility.h"

#include <cstdlib>
#include <iomanip>
#include <iostream>

int main(int argc, char** argv) {
    if (argc != 5 && argc != 9) {
        std::cerr << "usage: rod_mobility_calculator L_um d_um T_K eta_Pa_s "
                     "[chi_parallel chi_perp chi_tumble chi_roll]\n";
        return 2;
    }
    const double length = std::strtod(argv[1], nullptr);
    const double diameter = std::strtod(argv[2], nullptr);
    const double temperature = std::strtod(argv[3], nullptr);
    const double viscosity = std::strtod(argv[4], nullptr);
    const double chi_parallel = argc == 9 ? std::strtod(argv[5], nullptr) : 1.0;
    const double chi_perp = argc == 9 ? std::strtod(argv[6], nullptr) : 1.0;
    const double chi_tumble = argc == 9 ? std::strtod(argv[7], nullptr) : 1.0;
    const double chi_roll = argc == 9 ? std::strtod(argv[8], nullptr) : 1.0;
    try {
        const auto diffusion = rod_mobility::finite_cylinder(
            length, diameter, temperature, viscosity,
            chi_parallel, chi_perp, chi_tumble, chi_roll);
        std::cout << std::setprecision(10)
                  << "kBT_pN_um=" << rod_mobility::kbt(temperature) << '\n'
                  << "beta_per_pN_um=" << rod_mobility::beta(temperature) << '\n'
                  << "D_parallel_um2_s=" << diffusion.parallel << '\n'
                  << "D_perpendicular_um2_s=" << diffusion.perpendicular << '\n'
                  << "D_tumble_s-1=" << diffusion.tumble << '\n'
                  << "D_roll_s-1=" << diffusion.roll << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
