#include "gtest/gtest.h"
#include "include/utils.h"
#include "include/sarcomere.h"
#include <chrono>
#include <iostream>
#include <cmath>

TEST(NeighborList, RebuildTiming)
{
    int n_actins = 1500;
    int n_myosins = 130;
    std::vector<double> box{6.0, 4.0, 3.2};
    double actin_length = 1.0;
    double myosin_length = 1.5;
    double myosin_radius = 0.4;
    double myosin_radius_ratio = 0.4;
    double crosslinker_length = 0.06;
    double k_on = 1000.0;
    double k_off = 1.0;
    double base_lifetime = 0.001;
    double lifetime_coeff = 0.4;
    double diff_coeff_ratio = 10.0; // 0.5 / 0.05 from run_test.sh defaults
    double k_aa = 300.0;
    double kappa_aa = 50.0;
    double k_am = 50.0;
    double kappa_am = 200.0;
    double v_am = 5.0;
    std::string filename = "test.h5";
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    int seed = 0;
    int fix_myosin = 0;
    double dt = 0.00001;
    bool directional = true;

    Sarcomere model(n_actins, n_myosins, box, actin_length, myosin_length,
                    myosin_radius, myosin_radius_ratio, crosslinker_length,
                    k_on, k_off, base_lifetime, lifetime_coeff, diff_coeff_ratio,
                    k_aa, kappa_aa, k_am, kappa_am, v_am,
                    filename, rng, seed, fix_myosin, dt, directional);

    auto start = std::chrono::high_resolution_clock::now();
    model.neighbor_list.set_species_positions(model.actin.center_x, model.actin.center_y, model.actin.center_z,
                                              model.myosin.center_x, model.myosin.center_y, model.myosin.center_z);
    model.neighbor_list.rebuild_neighbor_list();
    auto end = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(end - start).count();
    std::cout << "Neighbor list rebuild took " << ms << " ms" << std::endl;
    gsl_rng_free(rng);
}

TEST(Vec, Distance)
{
    utils::vec a{0.0, 0.0, 0.0};
    utils::vec b{1.0, 2.0, 2.0};
    double expected = std::sqrt(1.0*1.0 + 2.0*2.0 + 2.0*2.0);
    EXPECT_NEAR(a.distance(b), expected, 1e-12);
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
