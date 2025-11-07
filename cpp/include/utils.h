#ifndef UTILS_H
#define UTILS_H

#include <array>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <utility>
#include <omp.h>
#include <iostream>
#include <stdexcept>


namespace utils {

inline double wrap_axis(double x, double L, bool periodic) {
    if (!periodic || L == 0.0) {
        return x;
    }
    return x - L * std::round(x / L);
}

//------------------------------------------------------------------------------
// The vec structure: represents a 3D vector/point and provides common operations.
//------------------------------------------------------------------------------
struct vec {
    double x;
    double y;
    double z;

    // Addition operator.
    vec operator+(const vec& p) const {
        return vec { x + p.x, y + p.y,  z + p.z };
    }
    vec& operator+=(const vec& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    // Subtraction operator.
    vec operator-(const vec& p) const {
        return vec { x - p.x, y - p.y, z - p.z };
    }
    vec& operator-=(const vec& other) {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }

    // Multiplication by scalar.
    vec operator*(double c) const {
        return vec { x * c, y * c , z * c};
    }
    friend vec operator*(double c, const vec& v) {
        return v * c;
    }

    // Division by scalar.
    vec operator/(double c) const {
        return vec { x / c, y / c , z / c};
    }

    // Euclidean norm.
    double norm() const {
        return std::sqrt(x*x + y*y + z*z);
    }

    // Unary minus.
    vec operator-() const {
        return vec { -x, -y, -z };
    }

    // Wrap the coordinates according to periodic boundary conditions.
    // The box vector is assumed to contain the periodic lengths in x and y.
    void pbc_wrap(const std::vector<double>& box, const std::array<bool,3>& periodic) {
        double Lx = box.size() > 0 ? box[0] : 0.0;
        double Ly = box.size() > 1 ? box[1] : 0.0;
        double Lz = box.size() > 2 ? box[2] : 0.0;
        x = wrap_axis(x, Lx, periodic[0]);
        y = wrap_axis(y, Ly, periodic[1]);
        z = wrap_axis(z, Lz, periodic[2]);
    }

    void pbc_wrap(const std::vector<double>& box) {
        static constexpr std::array<bool,3> default_mask{true, true, true};
        pbc_wrap(box, default_mask);
    }

    // Euclidean distance (no periodic boundaries).
    double distance(const vec& p) const {
        return std::sqrt((x - p.x) * (x - p.x) +
                         (y - p.y) * (y - p.y) +
                         (z - p.z) * (z - p.z));
    }

    // Distance squared with periodic boundary conditions.
    double distance_squared(const vec& p, const std::vector<double>& box, const std::array<bool,3>& periodic) const {
        double Lx = box.size() > 0 ? box[0] : 0.0;
        double Ly = box.size() > 1 ? box[1] : 0.0;
        double Lz = box.size() > 2 ? box[2] : 0.0;
        double dx = wrap_axis(x - p.x, Lx, periodic[0]);
        double dy = wrap_axis(y - p.y, Ly, periodic[1]);
        double dz = wrap_axis(z - p.z, Lz, periodic[2]);
        return dx * dx + dy * dy + dz * dz;
    }
    double distance_squared(const vec& p, const std::vector<double>& box) const {
        static constexpr std::array<bool,3> default_mask{true, true, true};
        return distance_squared(p, box, default_mask);
    }
    
    // Distance computed with periodic boundary conditions.
    double distance(const vec& p, const std::vector<double>& box, const std::array<bool,3>& periodic) const {
        return std::sqrt(distance_squared(p, box, periodic));
    }
    double distance(const vec& p, const std::vector<double>& box) const {
        static constexpr std::array<bool,3> default_mask{true, true, true};
        return distance(p, box, default_mask);
    }

    // Dot product.
    double dot(const vec& p) const {
        return x * p.x + y * p.y + z * p.z;
    }
    


    // Squared norm.
    double norm_squared() const {
        return x * x + y * y + z * z;
    }

    void normalize() {
        double n = norm();
        if (n > 1e-12) {
            x /= n;
            y /= n;
            z /= n;
        } else {
            x = 1.0; y = 0.0; z = 0.0;
        }
    }
    vec normalized() const {
        vec v = *this;
        v.normalize();
        return v;
    }

    vec cross(const vec& p) const {
        return vec {
            y * p.z - z * p.y,
            z * p.x - x * p.z,
            x * p.y - y * p.x
        };
    }  
    
};

inline vec rodrigues_rotate(const vec& x, const vec& rotvec) {
    double theta = rotvec.norm();
    if (theta < 1e-12) {
        vec wx = rotvec.cross(x);
        vec wwx = rotvec.cross(wx);
        return x + wx + 0.5 * wwx;
    }
    vec axis = rotvec / theta;
    vec axx = axis.cross(x);
    double c = std::cos(theta);
    double s = std::sin(theta);
    double axis_dot = axis.dot(x);
    return x * c + axx * s + axis * (axis_dot * (1.0 - c));
}

//------------------------------------------------------------------------------
// Free function declarations (definitions are provided in utils.cpp)
//------------------------------------------------------------------------------

// Overload operator<< for vec
std::ostream& operator<<(std::ostream &os, const vec &v);


bool compare_indices(const std::vector<int>& a, const std::vector<int>& b);

// Wrap a coordinate value using periodic boundary conditions.
double pbc_wrap(double x, double& box);

double angle_between(const vec& u1, const vec& u2);

vec pbc_diff_masked(const vec& a, const vec& b, const std::vector<double>& box, const std::array<bool,3>& periodic);

//------------------------------------------------------------------------------
// MoleculeConnection class declaration
//------------------------------------------------------------------------------
class MoleculeConnection {
public:
    MoleculeConnection();
    MoleculeConnection(int numA);

    // Add a connection from molecule A (index aIndex) to molecule B (index bIndex).
    void addConnection(int aIndex, int bIndex);
    // Delete a connection from molecule A to molecule B.
    void deleteConnection(int aIndex, int bIndex);
    // Delete all connections for molecule A.
    void deleteAllConnections(int aIndex);
    // Retrieve the connections for molecule A.
    const std::vector<int>& getConnections(int aIndex) const;

private:
    std::vector<std::vector<int>> connections;
};

//------------------------------------------------------------------------------
// Template functions (must be in the header)
//------------------------------------------------------------------------------

// Returns a vector of indices that sort the input vector in descending order.
template <typename T>
std::vector<size_t> sort_indices(const std::vector<T>& vec) {
    std::vector<size_t> indices(vec.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        indices[i] = i;
    }
    std::sort(indices.begin(), indices.end(),
              [&vec](size_t a, size_t b) { return vec[a] > vec[b]; });
    return indices;
}

// Overload += operator for std::vector<double>
std::vector<double>& operator+=(std::vector<double>& a, const std::vector<double>& b);

// Overload + operator for std::vector<double>
std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b);

// Reduces a temporary array (with one sub-array per thread) into a target array.
template <typename T>
void reduce_array(std::vector<std::vector<T>>& temp_array, std::vector<T>& target_array) {
    #pragma omp for
    for (size_t i = 0; i < target_array.size(); ++i) {
        for (int t = 0; t < omp_get_num_threads(); ++t) {
            target_array[i] += temp_array[t][i];
        }
    }
}

} // namespace utils

#endif // UTILS_H
