#include "components.h"
#include <cstdio>
#include <stdexcept>
#include <cmath>
#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//===================
// Filament Methods
//===================

// Default constructor.
Filament::Filament() : n(0), length(0) {
    // Optionally initialize other members if needed.
}

// Parameterized constructor.
Filament::Filament(int n0, double length0, std::vector<double> box0, gsl_rng* rng)
    : n(n0), length(length0), box(box0)
{
    center.resize(n);
    direction.resize(n);
    left_end.resize(n);
    right_end.resize(n);
    force.resize(n);
    torque.resize(n);
    velocity.resize(n);
    f_load.resize(n);
    cb_strength.resize(n);
    
    // Randomly initialize the center positions and theta values.
    for (int i = 0; i < n; i++) {
        center[i].x = gsl_ran_flat(rng, -0.5 * box[0], 0.5 * box[0]);
        center[i].y = gsl_ran_flat(rng, -0.5 * box[1], 0.5 * box[1]);
        center[i].z = gsl_ran_flat(rng, -0.5 * box[2], 0.5 * box[2]);
        double x = gsl_ran_gaussian(rng, 1.0);
        double y = gsl_ran_gaussian(rng, 1.0);
        double z = gsl_ran_gaussian(rng, 1.0);
        double norm = sqrt(x*x + y*y + z*z);
        direction[i].x = x / norm;
        direction[i].y = y / norm;
        direction[i].z = z / norm;
    }
    update_endpoints();
}

// Destructor.
Filament::~Filament() {
    printf("Filament destructor called\n");
}

// Copy constructor.
Filament::Filament(const Filament& other) {
    n = other.n;
    box = other.box;
    length = other.length;
    cb_strength = other.cb_strength;
    center = other.center;
    direction = other.direction;
    left_end = other.left_end;
    right_end = other.right_end;
    force = other.force;
    torque = other.torque;
    velocity = other.velocity;
    f_load = other.f_load;
    custom_features = other.custom_features;
}

// Displace function (translation only).
void Filament::displace(int& i, double& dx, double& dy, double& dz) {
    center[i].x += dx;
    center[i].y += dy;
    center[i].z += dz;
    center[i].pbc_wrap(box);
    update_endpoints(i);
}


// Update endpoints for the i-th filament in 3D.
void Filament::update_endpoints(int& i) {
    left_end[i] = center[i] - 0.5 * length * direction[i];
    right_end[i] = center[i] + 0.5 * length * direction[i];
}


// Update endpoints for all filaments.
void Filament::update_endpoints() {
    for (int i = 0; i < n; i++) {
        update_endpoints(i);
    }
}



// Update center positions for all filaments.
void Filament::update_center(std::vector<vec> new_center) {
    for (int i = 0; i < n; i++) {
        center[i] = new_center[i];
    }
    update_endpoints();
}

// Register a new 1D feature.
void Filament::register_feature(const std::string& name) {
    if (custom_features.find(name) == custom_features.end()) {
        custom_features[name] = std::vector<double>(n, 0.0); // Initialize with zeros.
        std::cout << "1D feature " << name << " registered successfully.\n";
    } else {
        std::cout << "1D feature " << name << " already exists.\n";
    }
}

// Overload operator[] to access custom features.
std::vector<double>& Filament::operator[](const std::string& name) {
    if (custom_features.find(name) == custom_features.end()) {
        throw std::runtime_error("Error: 1D feature " + name + " not found.");
    }
    return custom_features[name];
}

//===================
// Myosin Methods
//===================

// Default constructor.
Myosin::Myosin() : Filament() {
    // The base default constructor is automatically called.
}

// Parameterized constructor.
Myosin::Myosin(int n0, double length0, double radius0, std::vector<double> box0, gsl_rng* rng)
    : Filament(n0, length0, box0, rng), radius(radius0)
{
    // Additional initialization if needed.
}

// Copy constructor.
Myosin::Myosin(const Myosin& other) : Filament(other) {
    radius = other.radius;
}

