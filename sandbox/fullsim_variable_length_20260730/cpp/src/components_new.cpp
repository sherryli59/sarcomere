#include "components.h"
#include <cstdio>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <omp.h>

namespace {

constexpr double ZERO_TOL = 1e-8;

inline double get_component(const vec& v, int axis) {
    return (axis == 0) ? v.x : (axis == 1 ? v.y : v.z);
}

inline void set_component(vec& v, int axis, double value) {
    if (axis == 0) {
        v.x = value;
    } else if (axis == 1) {
        v.y = value;
    } else {
        v.z = value;
    }
}

// New: sample a unit direction with per-axis asymmetric bounds
// u_min[u] <= u[u] <= u_max[u] for u in {x,y,z}
vec sample_direction_with_asym_limits(const std::array<double,3>& u_min,
                                      const std::array<double,3>& u_max,
                                      gsl_rng* rng)
{
    // Fast rejection sampler on S^2
    const int max_attempts = 5000;
    for (int attempt = 0; attempt < max_attempts; ++attempt) {
        double x = gsl_ran_gaussian(rng, 1.0);
        double y = gsl_ran_gaussian(rng, 1.0);
        double z = gsl_ran_gaussian(rng, 1.0);
        double n = std::sqrt(x*x + y*y + z*z);
        if (n <= 1e-15) continue;
        x /= n; y /= n; z /= n;
        if (x >= u_min[0] && x <= u_max[0] &&
            y >= u_min[1] && y <= u_max[1] &&
            z >= u_min[2] && z <= u_max[2]) {
            return {x,y,z};
        }
    }

    // Fallback: set the tightest axis first, sample the rest, renormalize
    std::array<int,3> order{0,1,2};
    std::sort(order.begin(), order.end(), [&](int a,int b){
        return (u_max[a]-u_min[a]) < (u_max[b]-u_min[b]);
    });

    // Fix the tightest axis component uniformly within its range
    double comp[3]{0,0,0};
    int k = order[0];
    comp[k] = gsl_ran_flat(rng, u_min[k], u_max[k]);

    // Sample the remaining two from a circle of radius sqrt(1-comp[k]^2)
    double R = std::sqrt(std::max(0.0, 1.0 - comp[k]*comp[k]));
    double theta = gsl_ran_flat(rng, 0.0, 2.0*M_PI);
    int a = order[1], b = order[2];
    comp[a] = R * std::cos(theta);
    comp[b] = R * std::sin(theta);

    // If any component is outside its bounds, clamp and renormalize a couple times
    for (int it = 0; it < 3; ++it) {
        for (int ax = 0; ax < 3; ++ax) {
            comp[ax] = std::min(std::max(comp[ax], u_min[ax]), u_max[ax]);
        }
        double nn = std::sqrt(comp[0]*comp[0]+comp[1]*comp[1]+comp[2]*comp[2]);
        if (nn > 1e-15) { comp[0]/=nn; comp[1]/=nn; comp[2]/=nn; }
    }
    return {comp[0],comp[1],comp[2]};
}

} // namespace

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
      periodic_axes(other.periodic_axes),
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
    temp.pbc_wrap(box, periodic_axes);
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

void Filament::set_periodic_axes(const std::array<bool,3>& periodic) {
    periodic_axes = periodic;
    for (int i = 0; i < n; ++i) {
        vec wrapped{center_x[i], center_y[i], center_z[i]};
        wrapped.pbc_wrap(box, periodic_axes);
        center_x[i] = wrapped.x;
        center_y[i] = wrapped.y;
        center_z[i] = wrapped.z;
    }
    update_endpoints();
}

void Filament::initialize_within_box(gsl_rng* rng) {
    if (box.size() < 3) box.resize(3, 0.0);
    const double H[3] = {0.5*box[0], 0.5*box[1], 0.5*box[2]};

    for (int i = 0; i < n; ++i) {
        // (A) sample any center uniformly inside the box on each axis
        vec c{
            (H[0] > 0 ? gsl_ran_flat(rng, -H[0], H[0]) : 0.0),
            (H[1] > 0 ? gsl_ran_flat(rng, -H[1], H[1]) : 0.0),
            (H[2] > 0 ? gsl_ran_flat(rng, -H[2], H[2]) : 0.0)
        };

        // (B) compute per-axis asymmetric bounds for u: [u_min, u_max]
        // Non-periodic: keep both filament tips inside [-H, H] along that axis.
        //   c ± (L/2) * u ∈ [-H, H] ⇒
        //     lower = max(-2(H + c)/L, -2(H - c)/L)
        //     upper = min( 2(H - c)/L,  2(H + c)/L)
        // Periodic: [-1, 1]
        std::array<double,3> umin{-1.0,-1.0,-1.0}, umax{1.0,1.0,1.0};
        for (int ax = 0; ax < 3; ++ax) {
            if (!periodic_axes[ax]) {
                const double ca = (ax==0? c.x : ax==1? c.y : c.z);
                const double denom = std::max(length, ZERO_TOL);
                const double lower = std::max(
                    -2.0 * (H[ax] + ca) / denom,
                    -2.0 * (H[ax] - ca) / denom);
                const double upper = std::min(
                    2.0 * (H[ax] - ca) / denom,
                    2.0 * (H[ax] + ca) / denom);
                double u_min_ax = lower;
                double u_max_ax = upper;
                // clamp to unit-vector range
                umin[ax] = std::max(u_min_ax, -1.0);
                umax[ax] = std::min(u_max_ax,  1.0);
                // fallback if the interval vanished due to numerical issues
                if (umin[ax] > umax[ax]) {
                    const double mid = std::clamp(0.5 * (umin[ax] + umax[ax]), -1.0, 1.0);
                    umin[ax] = umax[ax] = mid;
                }
            }
        }

        // (C) sample a unit orientation respecting those asymmetric bounds
        vec u = sample_direction_with_asym_limits(umin, umax, rng);

        // (D) write state
        center_x[i] = c.x; center_y[i] = c.y; center_z[i] = c.z;
        direction_x[i] = u.x; direction_y[i] = u.y; direction_z[i] = u.z;

        update_endpoints(i);

        // zero aux
        force_x[i] = force_y[i] = force_z[i] = 0.0;
        torque_x[i] = torque_y[i] = torque_z[i] = 0.0;
        velocity_x[i] = velocity_y[i] = velocity_z[i] = 0.0;
        f_load[i] = 0.0; cb_status[i] = 0;
    }
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
