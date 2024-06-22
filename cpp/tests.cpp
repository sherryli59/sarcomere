//write an example google test
#include "gtest/gtest.h"
#include "include/mc.h"
#include "include/utils.h"
#include <vector>


// TEST(MC, Constructor){
//     double dt = 0.01;
//     double beta = 1;
//     double D = 0.1;
//     int update_dt_every = 100;
//     int save_every = 1000;
//     bool resume = false;
//     Sarcomere model;
//     MC mc(model, beta, dt, D, update_dt_every, save_every, resume);
//     EXPECT_EQ(mc.dt, dt);
//     EXPECT_EQ(mc.beta, beta);
//     EXPECT_EQ(mc.D, D);
//     EXPECT_EQ(mc.update_dt_every, update_dt_every);
//     EXPECT_EQ(mc.save_every, save_every);
//     EXPECT_EQ(mc.acc_rate, 0);
// }


// TEST(UTILS, point_segment_distance)
// {
//     using namespace utils;
//     double p[2] = {0.193261 -1.670936};
//     double p1[2] = {-6,-0.5};
//     double p2[2] = {-3,0.5};
//     std::vector<double> box = {10,10};
//     double d1 = point_segment_distance(p, p1, p2, box);
//     double d = point_segment_distance_new(p, p1, p2, box);
//     // distance should be smaller than 1  
//     printf("Distance: %f\n", d1);  
//     EXPECT_EQ(d1, d);
//     EXPECT_LT(d, 0.7);

// }

TEST(UTILS, segment_segment_distance)
{
    using namespace utils;
    // double p1[2] = {-0.160637,0.746029};
    // double p2[2] = {-2.945040,1.862764};
    // double q1[2] = {-0.375845,2.159613};
    // double q2[2] = {-2.720322,0.287869};
    double p1[2] = {5,5.5};
    double p2[2] = {5,4.5};
    double q1[2] = {5.5,5};
    double q2[2] = {4.5,5};
    // double p1[2] = {10,9.5};
    // double p2[2] = {10,10.5};
    // double q1[2] = {-0.5,0};
    // double q2[2] = {0.5,0};
    // double p1[2] = {1.34538, 1.93029};
    // double p2[2] = {-1.43645, 3.05342};
    // double q1[2] = {-1.13567, 0.0448536};
    // double q2[2] = {1.02296, -2.03849};
    std::vector<double> box = {10,10};
    double d = segment_segment_distance(p1, p2, q1,q2, box);
    // distance should be smaller than 1  
    printf("Distance: %f\n", d);  
    //EXPECT_GT(d, 0.5);
    EXPECT_EQ(d, 0);
}



int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
