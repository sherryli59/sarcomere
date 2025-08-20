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

// Define an alias for vec using the type from utils.
using vec = utils::vec;

class Filament {
public:
    int n;
    std::vector<double> box;
    double length;

    // Structure-of-arrays storage for vector quantities
    std::vector<double> center_x, center_y, center_z;
    std::vector<double> direction_x, direction_y, direction_z;
    std::vector<double> left_end_x, left_end_y, left_end_z;
    std::vector<double> right_end_x, right_end_y, right_end_z;
    std::vector<double> force_x, force_y, force_z;
    std::vector<double> torque_x, torque_y, torque_z;
    std::vector<double> velocity_x, velocity_y, velocity_z;
    std::vector<double> f_load;
    std::vector<int> cb_status;

    struct VecRef {
        double &x, &y, &z;
        operator vec() const { return {x, y, z}; }
        VecRef& operator=(const vec& v) { x = v.x; y = v.y; z = v.z; return *this; }
        VecRef& operator+=(const vec& v) { x += v.x; y += v.y; z += v.z; return *this; }
        VecRef& operator-=(const vec& v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
        VecRef& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }
        vec operator-(const VecRef& other) const { return {x - other.x, y - other.y, z - other.z}; }
        vec operator+(const VecRef& other) const { return {x + other.x, y + other.y, z + other.z}; }
        vec operator*(double s) const { return {x * s, y * s, z * s}; }
        friend vec operator*(double s, const VecRef& v) { return {v.x * s, v.y * s, v.z * s}; }
        vec operator/(double s) const { return {x / s, y / s, z / s}; }
        void normalize() {
            double norm = std::sqrt(x * x + y * y + z * z);
            if (norm > 0) { x /= norm; y /= norm; z /= norm; }
        }
        double dot(const VecRef& other) const {
            return x * other.x + y * other.y + z * other.z;
        }
        double norm() const { return std::sqrt(x * x + y * y + z * z); }
    };

    class VecArray {
        std::vector<double>* x; std::vector<double>* y; std::vector<double>* z;
    public:
        VecArray(std::vector<double>& x_, std::vector<double>& y_, std::vector<double>& z_)
            : x(&x_), y(&y_), z(&z_) {}
        VecRef operator[](size_t i) { return VecRef{(*x)[i], (*y)[i], (*z)[i]}; }
        const VecRef operator[](size_t i) const {
            return VecRef{const_cast<double&>((*x)[i]), const_cast<double&>((*y)[i]), const_cast<double&>((*z)[i])};
        }
        size_t size() const { return x->size(); }
    };

    VecArray center;
    VecArray direction;
    VecArray left_end;
    VecArray right_end;
    VecArray force;
    VecArray torque;
    VecArray velocity;

    // Dictionary to store 1D features.
    std::unordered_map<std::string, std::vector<double>> custom_features;

    // Constructors and destructor.
    Filament();
    Filament(int n0, double length0, std::vector<double> box0, gsl_rng* rng);
    virtual ~Filament();
    Filament(const Filament& other);

    // Member functions.
    void displace(int& i, double& dx, double& dy, double& dz);
    void update_endpoints(int& i);
    void update_endpoints();
    void update_center(std::vector<vec> new_center);

    // Register a new 1D feature of length n.
    void register_feature(const std::string& name);

    // Operator overload for easy access to 1D features.
    std::vector<double>& operator[](const std::string& name);
};

// Helper to reduce temporary arrays into VecArray targets
void reduce_array(std::vector<std::vector<vec>>& temp_array, Filament::VecArray& target_array);

class Myosin : public Filament {
public:
    double radius;
    // Constructors.
    Myosin();
    Myosin(int n0, double length0, double radius0, std::vector<double> box0, gsl_rng* rng);
    Myosin(const Myosin& other);
};

#endif // COMPONENTS_H
