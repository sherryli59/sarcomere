#include "gtest/gtest.h"
#include "include/utils.h"
#include <cmath>

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
