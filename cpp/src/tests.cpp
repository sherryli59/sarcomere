#include <gtest/gtest.h>
#include <tuple>
#include <vector>
#include <map>
#include <random>
#include <cmath>
#include "geometry.h"
#include "utils.h"
#include "interaction.h"
#include "autodiff/forward/real.hpp"  // Adjust the include path as needed.
using autodiff::real;
using vec = utils::vec;

// Define a small epsilon for floating point comparisons
constexpr double EPS = 1e-6;


// Test case for two parallel segments that are offset by 1 in the y direction.
TEST(SegmentSegmentDistanceTest, ParallelSegments) {
    // Define endpoints as arrays of autodiff::real.
    real A[3] = {0, 0, 0};
    real B[3] = {1, 0, 0};
    real C[3] = {0, 1, 0};
    real D[3] = {1, 1, 0};
    std::vector<double> box = {100, 100, 100}; // Large box to avoid PBC effects.
    
    // For parallel segments offset by 1, the minimum distance should be 1.
    real dist = segment_segment_distance(A, B, C, D, box);
    EXPECT_NEAR(dist.val(), 1.0, 1e-6);
}

// // Test case for intersecting segments (which should return a distance of zero).
// TEST(SegmentSegmentDistanceTest, IntersectingSegments) {
//     // These segments cross each other.
//     real A[3] = {0, 0, 0};
//     real B[3] = {1, 1, 0};
//     real C[3] = {0, 1, 0};
//     real D[3] = {1, 0, 0};
//     std::vector<double> box = {100, 100, 100};
    
//     real dist = segment_segment_distance(A, B, C, D, box);
//     EXPECT_NEAR(dist.val(), 0.0, 1e-6);
// }

// // Test case for skew segments in parallel planes.
// // Here segment1 is from (0,0,0) to (1,1,0) and segment2 is from (0,1,1) to (1,0,1).
// // Their xy projections intersect and the planes are separated by 1, so the minimal distance is 1.
// TEST(SegmentSegmentDistanceTest, SkewSegments) {
//     real A[3] = {0, 0, 0};
//     real B[3] = {1, 1, 0};
//     real C[3] = {0, 1, 1};
//     real D[3] = {1, 0, 1};
//     std::vector<double> box = {100, 100, 100};
    
//     real dist = segment_segment_distance(A, B, C, D, box);
//     EXPECT_NEAR(dist.val(), 1.0, 1e-6);
// }

TEST(SegmentDistanceTest, RandomizedConsistencyCheck) {
    const int num_trials = 1000;
    constexpr double delta = 1e-2;  // Small reduction for testing empty subsegment.
    //std::random_device rd;
    double seed = 42;
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dis_center(0.0, 1.0);
    std::uniform_real_distribution<double> dis_offset(0.0, 1.0);

    // Generate a random unit vector using spherical coordinates.
    auto random_unit_vector = [&]() -> vec {
        std::uniform_real_distribution<double> dis_phi(0, 2 * M_PI);
        std::uniform_real_distribution<double> dis_theta(0, M_PI);
        double phi = dis_phi(gen);
        double theta = dis_theta(gen);
        double x = std::sin(theta) * std::cos(phi);
        double y = std::sin(theta) * std::sin(phi);
        double z = std::cos(theta);
        return {x, y, z};
    };

    // Generate a random segment of length 1 centered at a given point.
    auto generate_random_segment = [&](const vec &center) -> std::pair<vec, vec> {
        vec dir = random_unit_vector();
        double half_length = 0.5;
        vec A = {center.x - half_length * dir.x,
                 center.y - half_length * dir.y,
                 center.z - half_length * dir.z};
        vec B = {center.x + half_length * dir.x,
                 center.y + half_length * dir.y,
                 center.z + half_length * dir.z};
        return std::make_pair(A, B);
    };

    // Use a large periodic box.
    std::vector<double> box = {10, 10, 10};

    for (int i = 0; i < num_trials; ++i) {
        // Generate two random centers.
        vec center1 = {dis_center(gen), dis_center(gen), dis_center(gen)};
        vec offset = {dis_offset(gen), dis_offset(gen), dis_offset(gen)};
        vec center2 = {center1.x + offset.x, center1.y + offset.y, center1.z + offset.z};

        // Generate segments for each center.
        auto seg1 = generate_random_segment(center1);
        auto seg2 = generate_random_segment(center2);
        vec A = seg1.first, B = seg1.second;
        vec C = seg2.first, D = seg2.second;
        geometry::apply_pbc(A, B, C, D, box);
        // Compute the minimum distance and associated info.
        auto [d_min, info] = geometry::segment_segment_distance_w_normal(A, B, C, D, box);
        // Compute the subsegment of AB within distance d_min of CD.
        auto [interval_length, start, end] = geometry::subsegment_within_distance(A, B, C, D, d_min+delta);

        // If the subsegment is empty, that's a failure.
        if (interval_length <= EPS) {
            std::cerr << "Trial " << i << " FAILED: Expected non-empty subsegment for d = " << d_min+delta << "\n";
            std::cerr << "Segment A: " << A << ", B: " << B << "\n";
            std::cerr << "Segment C: " << C << ", D: " << D << "\n";
            std::cerr << "subsegment_within_distance returned: (length: " << interval_length
                      << ", start: " << start << ", end: " << end << ")\n";
            std::cerr << "segment_segment_distance_w_normal returned: (d_min: " << d_min 
                      << ", start: " << info["start"] << ", end: " << info["end"]
                      << ", normal: " << info["normal"] << ")\n";
            FAIL() << "Non-empty subsegment expected but got empty.";
        }
        double d_lower = std::max(0.0, d_min - delta);
        // Now, if we reduce the threshold slightly, the subsegment should be empty.
        auto [interval_length_reduced, start_reduced, end_reduced] =
            geometry::subsegment_within_distance(A, B, C, D, d_lower);
        if (std::abs(interval_length_reduced) > EPS) {
            std::cerr << "Trial " << i << " FAILED: Expected empty subsegment for d = " << (d_lower) << "\n";
            std::cerr << "Segment A: " << A << ", B: " << B << "\n";
            std::cerr << "Segment C: " << C << ", D: " << D << "\n";
            std::cerr << "subsegment_within_distance returned: (length: " << interval_length_reduced
                      << ", start: " << start_reduced << ", end: " << end_reduced << ")\n";
            std::cerr << "segment_segment_distance_w_normal returned: (d_min: " << d_min 
                      << ", start: " << info["start"] << ", end: " << info["end"]
                      << ", normal: " << info["normal"] << ")\n";
            FAIL() << "Empty subsegment expected but got non-empty.";
        }
    }
}

TEST(SegmentDistanceTest, NoOverlap) {
    vec A = {0, 0, 0};
    vec B = {1, 0, 0};
    vec C = {0, 1, 0};
    vec D = {1, 1, 0};
    double d = 0.1;
    std::vector<double> box = {10, 10, 10};
    geometry::apply_pbc(A, B, C, D, box);
    auto [interval_length, start, end] = geometry::subsegment_within_distance(A, B, C, D, d);
    
    // Since min distance > d, expect interval length to be 0
    EXPECT_NEAR(interval_length, 0.0, EPS);
}

TEST(SegmentDistanceTest, ClosestPointComputation) {
    vec A = {0, 0, 0};
    vec B = {1, 0, 0};
    vec C = {0, 1, 0};
    vec D = {1, 1, 0};
    std::vector<double> box = {10, 10, 10};
    printf("no autodiff\n");
    auto [dist, info] = geometry::segment_segment_distance_w_normal(A, B, C, D, box);
    dist = geometry::segment_segment_distance(A, B, C, D, box);
    EXPECT_GT(dist, 0);
    EXPECT_NEAR(dist, 1.0, EPS);  // Expected min distance
}

TEST(SegmentDistanceTest, EdgeCaseParallel) {
    vec A = {0, 0, 0};
    vec B = {1, 0, 0};
    vec C = {0, 0.1, 0};
    vec D = {1, 0.1, 0};
    double d = 0.15;
    std::vector<double> box = {10, 10, 10};
    geometry::apply_pbc(A, B, C, D, box);
    auto [interval_length, start, end] = geometry::subsegment_within_distance(A, B, C, D, d);

    // Since the whole segment is within distance `d`, interval should be full length.
    EXPECT_NEAR(interval_length, 1.0, EPS);
}


TEST(SegmentDistanceTest, Sample1) {

    vec A = {-4.382204, -0.462398, -0.077329};
    vec B = {-5.218863, -0.174826, 0.388831};
    vec C = {-4.731318, -0.337676, 0.330925};
    vec D = {-4.450215, -0.907024, 1.689903};
    double d = 0.2;
    std::vector<double> box = {12, 5, 3};
    geometry::apply_pbc(A, B, C, D, box);
    auto [interval_length, start, end] = geometry::subsegment_within_distance(A, B, C, D, d);
    auto [dist, info] = geometry::segment_segment_distance_w_normal(A, B, C, D, box);
    printf("dist: %f, interval_length: %f\n", dist, interval_length);
}

TEST(SegmentDistanceTest, IntersectingSegments) {
    vec A = {0, 0, 0};
    vec B = {1, 1, 0};
    vec C = {0, 1, 0};
    vec D = {1, 0, 0};
    double d = 0.1;
    std::vector<double> box = {10, 10, 10};
    auto [dist, info] = geometry::segment_segment_distance_w_normal(A, B, C, D, box);
    
    // If segments intersect, min distance should be 0.
    EXPECT_NEAR(dist, 0.0, EPS);
}
