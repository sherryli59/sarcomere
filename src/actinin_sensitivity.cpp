#include <vector>
#include <random>
#include <cmath>
#include <opencv2/opencv.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

class AlphaActinin
{
    public:
        int n;
        double radius;
        double ** xs;
        std::vector <double> box;

        AlphaActinin(){
            n = 0;
            radius = 0;
            xs = NULL;
        }

        AlphaActinin(int& n0, double& radius0, std::vector <double> box0, gsl_rng * rng){
            n = n0;
            radius = radius0;
            xs = new double*[n];
            box = box0;
            for (int i = 0; i < n; i++){
                xs[i] = new double[2];
                xs[i][0] = gsl_ran_flat(rng, -0.5*box[0],0.5*box[0]);
                xs[i][1] = gsl_ran_flat(rng, -0.5*box[1], 0.5*box[1]);
            }
        }

        ~AlphaActinin(){
            for (int i = 0; i < n; i++){
                delete[] xs[i];
            }
            delete[] xs;
        }

        AlphaActinin(const AlphaActinin& other){
            n = other.n;
            radius = other.radius;
            xs = new double*[n];
            for (int i = 0; i < n; i++){
                xs[i] = new double[2];
                xs[i][0] = other.xs[i][0];
                xs[i][1] = other.xs[i][1];
            }
        }
};


class Actin {
    public:
        int n;
        double ** xs;
        double * thetas;
        double ** left_endpts;
        double ** right_endpts;
        
        Actin(){
            n = 0;
            xs = NULL;
            thetas = NULL;
        }

        Actin(int n0, gsl_rng * rng, double *box){
            n = n0;
            xs = new double*[n];
            thetas = new double[n];
            left_endpts = new double*[n];
            right_endpts = new double*[n];
            for (int i = 0; i < n; i++){
                xs[i] = new double[2];
                xs[i][0] = box[0] * gsl_ran_flat(rng, -0.5, 0.5);
                xs[i][1] = box[1] * gsl_ran_flat(rng, -0.5, 0.5);
                thetas[i] = gsl_ran_flat(rng, 0, 2*M_PI);
                left_endpts[i] = new double[2];
                right_endpts[i] = new double[2];
            }
            update_endpoints();
        }

        void update_endpoints(){
            for (int i = 0; i < n; i++){
                std::vector <double> segments(2);
                segments[0] = cos(thetas[i]);
                segments[1] = sin(thetas[i]);
                left_endpts[i][0] = xs[i][0] - 0.5*segments[0];
                left_endpts[i][1] = xs[i][1] - 0.5*segments[1];
                right_endpts[i][0] = xs[i][0] + 0.5*segments[0];
                right_endpts[i][1] = xs[i][1] + 0.5*segments[1];
            }
        }
};

void step(Actin& actin, double dt, double D_t, double D_r, gsl_rng * rng){
    
    // compute the interaction between actin / alpha-actinin
    // get number of bound dimers
    // compute total force on each actin from dimers
    // update actin positions
    if (n_dimers==0){
        for (int i = 0; i < actin.n; i++){
            actin.xs[i][0] += gsl_ran_gaussian(rng, sqrt(2*D_t*dt));
            actin.xs[i][0] += f * cos(actin.thetas[i]) * dt;
        }
    }
    else{
        double dx = gsl_ran_gaussian(rng, sqrt(2*D_t*dt));
         for (int i = 0; i < actin.n; i++){
            actin.xs[i][0] += dx;
        }
    }
    // binding and unbinding of actinin
    // force dependent rate equation for the actinin binding / unbinding
    
    

}



int main() {
    int n_actin = 2;
    int n_actinin = 10;
    double box_size = 10.0;
    double dt = 0.01;
    double D_t = 1.0;
    double D_r = 0.1;


    for (int step = 0; step < 10000; ++step) {
        sim.step();
    }

    return 0;
}



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


