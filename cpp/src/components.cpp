#include "components.h"
#include <cstdio>
#include <stdexcept>
#include <cmath>
#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <omp.h>

//===================
// Filament Methods
//===================

// Default constructor.
Filament::Filament()
    : n(0), length(0),
      center(center_x, center_y, center_z),
      direction(direction_x, direction_y, direction_z),
      left_end(left_end_x, left_end_y, left_end_z),
      right_end(right_end_x, right_end_y, right_end_z),
      force(force_x, force_y, force_z),
      torque(torque_x, torque_y, torque_z),
      velocity(velocity_x, velocity_y, velocity_z) {}

// Parameterized constructor.
Filament::Filament(int n0, double length0, std::vector<double> box0, gsl_rng* rng)
    : n(n0), box(box0), length(length0),
      center(center_x, center_y, center_z),
      direction(direction_x, direction_y, direction_z),
      left_end(left_end_x, left_end_y, left_end_z),
      right_end(right_end_x, right_end_y, right_end_z),
      force(force_x, force_y, force_z),
      torque(torque_x, torque_y, torque_z),
      velocity(velocity_x, velocity_y, velocity_z)
{
    center_x.resize(n);
    center_y.resize(n);
    center_z.resize(n);
    direction_x.resize(n);
    direction_y.resize(n);
    direction_z.resize(n);
    left_end_x.resize(n);
    left_end_y.resize(n);
    left_end_z.resize(n);
    right_end_x.resize(n);
    right_end_y.resize(n);
    right_end_z.resize(n);
    force_x.resize(n);
    force_y.resize(n);
    force_z.resize(n);
    torque_x.resize(n);
    torque_y.resize(n);
    torque_z.resize(n);
    velocity_x.resize(n);
    velocity_y.resize(n);
    velocity_z.resize(n);
    f_load.resize(n);
    cb_status.resize(n);

    // Randomly initialize the center positions and directions.
    for (int i = 0; i < n; i++) {
        center_x[i] = gsl_ran_flat(rng, -0.5 * box[0], 0.5 * box[0]);
        center_y[i] = gsl_ran_flat(rng, -0.5 * box[1], 0.5 * box[1]);
        center_z[i] = gsl_ran_flat(rng, -0.5 * box[2], 0.5 * box[2]);
        double x = gsl_ran_gaussian(rng, 1.0);
        double y = gsl_ran_gaussian(rng, 1.0);
        double z = gsl_ran_gaussian(rng, 1.0);
        double norm = sqrt(x*x + y*y + z*z);
        direction_x[i] = x / norm;
        direction_y[i] = y / norm;
        direction_z[i] = z / norm;
    }
    update_endpoints();
}

// Destructor.
Filament::~Filament() {
    printf("Filament destructor called\n");
}

// Copy constructor.
Filament::Filament(const Filament& other)
    : n(other.n), box(other.box), length(other.length),
      center(center_x, center_y, center_z),
      direction(direction_x, direction_y, direction_z),
      left_end(left_end_x, left_end_y, left_end_z),
      right_end(right_end_x, right_end_y, right_end_z),
      force(force_x, force_y, force_z),
      torque(torque_x, torque_y, torque_z),
      velocity(velocity_x, velocity_y, velocity_z),
      custom_features(other.custom_features)
{
    center_x = other.center_x;
    center_y = other.center_y;
    center_z = other.center_z;
    direction_x = other.direction_x;
    direction_y = other.direction_y;
    direction_z = other.direction_z;
    left_end_x = other.left_end_x;
    left_end_y = other.left_end_y;
    left_end_z = other.left_end_z;
    right_end_x = other.right_end_x;
    right_end_y = other.right_end_y;
    right_end_z = other.right_end_z;
    force_x = other.force_x;
    force_y = other.force_y;
    force_z = other.force_z;
    torque_x = other.torque_x;
    torque_y = other.torque_y;
    torque_z = other.torque_z;
    velocity_x = other.velocity_x;
    velocity_y = other.velocity_y;
    velocity_z = other.velocity_z;
    f_load = other.f_load;
    cb_status = other.cb_status;
}

// Displace function (translation only).
void Filament::displace(int& i, double& dx, double& dy, double& dz) {
    center_x[i] += dx;
    center_y[i] += dy;
    center_z[i] += dz;
    vec temp{center_x[i], center_y[i], center_z[i]};
    temp.pbc_wrap(box);
    center_x[i] = temp.x;
    center_y[i] = temp.y;
    center_z[i] = temp.z;
    update_endpoints(i);
}

// Update endpoints for the i-th filament in 3D.
void Filament::update_endpoints(int& i) {
    left_end_x[i] = center_x[i] - 0.5 * length * direction_x[i];
    left_end_y[i] = center_y[i] - 0.5 * length * direction_y[i];
    left_end_z[i] = center_z[i] - 0.5 * length * direction_z[i];
    right_end_x[i] = center_x[i] + 0.5 * length * direction_x[i];
    right_end_y[i] = center_y[i] + 0.5 * length * direction_y[i];
    right_end_z[i] = center_z[i] + 0.5 * length * direction_z[i];
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
        center_x[i] = new_center[i].x;
        center_y[i] = new_center[i].y;
        center_z[i] = new_center[i].z;
    }
    update_endpoints();
}

// Reduce thread-local vec arrays into a VecArray target
void reduce_array(std::vector<std::vector<vec>>& temp_array, Filament::VecArray& target_array) {
    #pragma omp for
    for (size_t i = 0; i < target_array.size(); ++i) {
        for (int t = 0; t < omp_get_num_threads(); ++t) {
            target_array[i] += temp_array[t][i];
        }
    }
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

