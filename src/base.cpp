#include <vector>
#include <random>
#include <cmath>
#include <opencv2/opencv.hpp>

using namespace std;

class Rod {
public:
    double x, y;    // Center of mass position
    double theta;   // Orientation angle
    double length;  // Rod length

    Rod(double x, double y, double theta, double length)
        : x(x), y(y), theta(theta), length(length) {}
};

class Simulation {
private:
    std::vector<Rod> rods;
    double dt;  // Time step
    double D_t; // Translational diffusion coefficient
    double D_r; // Rotational diffusion coefficient
    std::mt19937 gen;
    std::normal_distribution<> normal_dist;

public:
    Simulation(int num_rods, double box_size, double rod_length, double dt, double D_t, double D_r)
        : dt(dt), D_t(D_t), D_r(D_r), gen(std::random_device{}()), normal_dist(0.0, 1.0) {
        // Initialize rods with random positions and orientations
        std::uniform_real_distribution<> uniform_dist(0.0, box_size);
        std::uniform_real_distribution<> angle_dist(0.0, 2 * M_PI);
        
        rods.reserve(num_rods);
        for (int i = 0; i < num_rods; ++i) {
            rods.emplace_back(uniform_dist(gen), uniform_dist(gen), angle_dist(gen), rod_length);
        }
    }

    void step() {
        for (auto& rod : rods) {
            // Translational motion
            double dx = sqrt(2 * D_t * dt) * normal_dist(gen);
            double dy = sqrt(2 * D_t * dt) * normal_dist(gen);
            
            // Rotational motion
            double dtheta = sqrt(2 * D_r * dt) * normal_dist(gen);
            
            // Update rod state
            rod.x += dx;
            rod.y += dy;
            rod.theta += dtheta;
            
            // Normalize angle to [0, 2π)
            rod.theta = fmod(rod.theta, 2 * M_PI);
            if (rod.theta < 0) rod.theta += 2 * M_PI;
        }
    }


    void draw_state(const std::string& filename, int image_size = 800) {
        // Create a white image
        cv::Mat img(image_size, image_size, CV_8UC3, cv::Scalar(255, 255, 255));

        // Calculate scaling factor
        double scale = image_size / box_size;

        for (const auto& rod : rods) {
            // Calculate rod endpoints
            double x1 = rod.x + 0.5 * rod.length * std::cos(rod.theta);
            double y1 = rod.y + 0.5 * rod.length * std::sin(rod.theta);
            double x2 = rod.x - 0.5 * rod.length * std::cos(rod.theta);
            double y2 = rod.y - 0.5 * rod.length * std::sin(rod.theta);

            // Scale and flip y-coordinate (image origin is top-left)
            cv::Point p1(static_cast<int>(x1 * scale), image_size - static_cast<int>(y1 * scale));
            cv::Point p2(static_cast<int>(x2 * scale), image_size - static_cast<int>(y2 * scale));

            // Draw the rod
            cv::line(img, p1, p2, cv::Scalar(0, 0, 255), 2);  // Red color, thickness 2
        }

        // Save the image
        cv::imwrite(filename, img);
    }
};




int main() {
    int num_rods = 1000;
    double box_size = 100.0;
    double rod_length = 1.0;
    double dt = 0.01;
    double D_t = 1.0;
    double D_r = 0.1;

    Simulation sim(num_rods, box_size, rod_length, dt, D_t, D_r);

    for (int step = 0; step < 10000; ++step) {
        sim.step();
    }

    return 0;
}