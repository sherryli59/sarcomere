#ifndef COMPONENTS_H
#define COMPONENTS_H

#ifndef UTILS_H
#include "utils.h"
#endif

#include <cmath>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
using vec = utils::vec;


class Filament {
public:
    int n;
    std::vector<double> box;
    double length;
    std::vector<vec> center;
    std::vector<double> theta;
    std::vector<vec> left_end;
    std::vector<vec> right_end;
    std::vector<vec> force;
    std::vector<double> angular_force;
    std::vector<vec> velocity;
    std::vector<double> tension;
    std::vector<double> cb_strength;

    // Dictionary to store 1D features
    std::unordered_map<std::string, std::vector<double>> custom_features;

    Filament() {
        n = 0;
        length = 0;
    }

    Filament(int n0, double length0, std::vector<double> box0, gsl_rng* rng) {
        n = n0;
        length = length0;
        box = box0;
        center.resize(n);
        theta.resize(n);
        left_end.resize(n);
        right_end.resize(n);
        force.resize(n);
        angular_force.resize(n);
        velocity.resize(n);
        tension.resize(n);
        cb_strength.resize(n);
        for (int i = 0; i < n; i++) {
            center[i].x = gsl_ran_flat(rng, -0.5 * box[0], 0.5 * box[0]);
            center[i].y = gsl_ran_flat(rng, -0.5 * box[1], 0.5 * box[1]);
            theta[i] = gsl_ran_flat(rng, 0, 2 * M_PI);
        }
        update_endpoints();
    }

    virtual ~Filament() {
        printf("Filament destructor called\n");
    }

    Filament(const Filament& other) {
        n = other.n;
        box = other.box;
        length = other.length;
        cb_strength = other.cb_strength;
        center = other.center;
        theta = other.theta;
        left_end = other.left_end;
        right_end = other.right_end;
        force = other.force;
        angular_force = other.angular_force;
        velocity = other.velocity;
        tension = other.tension;
        custom_features = other.custom_features;
    }

    void displace(int& i, double& dx, double& dy) {
        center[i].x += dx;
        center[i].y += dy;
        center[i].pbc_wrap(box);
        update_endpoints_i(i);
    }

    void displace(int& i, double& dx, double& dy, double& dtheta) {
        center[i].x += dx;
        center[i].y += dy;
        center[i].pbc_wrap(box);
        theta[i] += dtheta;
        utils::angle_wrap(theta[i]);
        update_endpoints_i(i);
    }

    void update_endpoints_i(int& i) {
        std::vector<double> segments(2);
        segments[0] = length * cos(theta[i]);
        segments[1] = length * sin(theta[i]);

        left_end[i].x = center[i].x - 0.5 * segments[0];
        left_end[i].y = center[i].y - 0.5 * segments[1];
        right_end[i].x = center[i].x + 0.5 * segments[0];
        right_end[i].y = center[i].y + 0.5 * segments[1];
    }

    void update_endpoints() {
        for (int i = 0; i < n; i++) {
            update_endpoints_i(i);
        }
    }

    void update_theta(std::vector<double> new_theta) {
        for (int i = 0; i < n; i++) {
            theta[i] = new_theta[i];
        }
        update_endpoints();
    }

    void update_center(std::vector<vec> new_center) {
        for (int i = 0; i < n; i++) {
            center[i] = new_center[i];
        }
        update_endpoints();
    }

    // Function to register a new 1D feature of length n
    void register_feature(const std::string& name) {
        if (custom_features.find(name) == custom_features.end()) {
            custom_features[name] = std::vector<double>(n, 0.0); // Initialize with zeros
            std::cout << "1D feature " << name << " registered successfully.\n";
        } else {
            std::cout << "1D feature " << name << " already exists.\n";
        }
    }

    // Access operator overloading for easy access to 1D features
    std::vector<double>& operator[](const std::string& name) {
        if (custom_features.find(name) == custom_features.end()) {
            throw std::runtime_error("Error: 1D feature " + name + " not found.");
        }
        return custom_features[name];
    }
};



class Myosin: public Filament
{
    public:
        double radius;
        
        Myosin(){
            Filament();
        }

        Myosin(int n0, double length0, double radius0,  std::vector <double> box0, gsl_rng * rng):
            Filament(n0, length0, box0,  rng)
            {
                radius = radius0;
            }   

        Myosin(const Myosin& other):
            Filament(other)
            {
                radius = other.radius;
            }
};



#endif
