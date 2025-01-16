
#ifndef INTERACTION_H
#define INTERACTION_H


#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#ifndef COMPONENTS_H
#include "components.h"
#endif

using Eigen::VectorXd;
using autodiff::ArrayXreal;
using autodiff::real;


// smooth_round function
template<typename T>
T smooth_round(T x) {
    return round(val(x));
}

// Function to calculate the norm of an array
template<typename T>
T norm(T* x, int dim) {
    T norm_val = 0;
    for (int i = 0; i < dim; i++) {
        norm_val += x[i] * x[i];
    }
    return sqrt(norm_val);
}

// Function to calculate the norm of a vector
template<typename T>
T norm(const std::vector<T>& x) {
    T norm_val = 0;
    for (int i = 0; i < x.size(); i++) {
        norm_val += x[i] * x[i];
    }
    return sqrt(norm_val);
}

// Periodic boundary conditions (PBC) wrap for arrays
template<typename T>
void pbc_wrap(T* x, const std::vector<double>& box) {
    x[0] = x[0] - box[0] * smooth_round(x[0] / box[0]);
    x[1] = x[1] - box[1] * smooth_round(x[1] / box[1]);
}

// Periodic boundary conditions (PBC) wrap for vectors
template<typename T>
void pbc_wrap(std::vector<T>& x, const std::vector<double>& box) {
    x[0] = x[0] - box[0] * smooth_round(x[0] / box[0]);
    x[1] = x[1] - box[1] * smooth_round(x[1] / box[1]);
}

// Periodic boundary condition for a single coordinate
template<typename T>
T pbc_wrap(T x, const T& box) {
    return x - box * smooth_round(x / box);
}

// Function to wrap an angle between -π and π
template<typename T>
void angle_wrap(T& theta) {
    theta = theta - 2 * M_PI * smooth_round(theta / (2 * M_PI));
}


// Template function for computing the shortest distance between a point and a segment (with PBC)
template<typename T>
T point_segment_distance(const T* x, const T* a, const T* b, const std::vector<double>& box, std::vector<T>& norm_vec) {
    //print x a b
    T ab[2] = {b[0] - a[0], b[1] - a[1]};
    T ab_norm = sqrt(ab[0] * ab[0] + ab[1] * ab[1]);
    T ab_normalized[2] = {ab[0] / ab_norm, ab[1] / ab_norm};

    T ap[2] = {x[0] - a[0], x[1] - a[1]};
    pbc_wrap(ap, box);
    T ap_ab = ap[0] * ab_normalized[0] + ap[1] * ab_normalized[1];
    ap_ab = std::max(T(0), std::min(ap_ab, ab_norm));
    norm_vec = {ap[0] - ap_ab * ab_normalized[0], ap[1] - ap_ab * ab_normalized[1]};
    return sqrt(norm_vec[0] * norm_vec[0] + norm_vec[1] * norm_vec[1]);
}

// Orientation function with PBC
template<typename T>
int orientation(const T* p, const T* q, const T* r, const std::vector<double>& box) {
    T pq[2] = {q[0] - p[0], q[1] - p[1]};
    T qr[2] = {r[0] - q[0], r[1] - q[1]};
    T val = pq[1] * qr[0] - pq[0] * qr[1];
    if (val == T(0)) return 0;  // Collinear
    return (val > 0) ? 1 : 2;  // Clockwise or Counterclockwise
}


// Template function for computing the shortest distance between two segments with PBC
template<typename T>
T segment_segment_distance(const T* a, const T* b, const T* c, const T* d, const std::vector<double>& box) {
    std::vector<T> norm_vec;
    T a_cd = point_segment_distance(a, c, d, box, norm_vec);
    T b_cd = point_segment_distance(b, c, d, box, norm_vec);
    T c_ab = point_segment_distance(c, a, b, box, norm_vec);
    T d_ab = point_segment_distance(d, a, b, box, norm_vec);
    T dist = std::min({a_cd, b_cd, c_ab, d_ab});
    // Check if segments intersect
    T ab_center[2] = {(a[0] + b[0]) / 2, (a[1] + b[1]) / 2};
    T cd_center[2] = {(c[0] + d[0]) / 2, (c[1] + d[1]) / 2};

    T displacement[2] = {
        box[0] * smooth_round((ab_center[0] - cd_center[0]) / box[0]),
        box[1] * smooth_round((ab_center[1] - cd_center[1]) / box[1])
    };

    // Adjust endpoints for displacement
    T a0[2] = {a[0] - displacement[0], a[1] - displacement[1]};
    T b0[2] = {b[0] - displacement[0], b[1] - displacement[1]};

    // Calculate orientation values
    int o1 = orientation(a0, b0, c, box);
    int o2 = orientation(a0, b0, d, box);
    int o3 = orientation(c, d, a0, box);
    int o4 = orientation(c, d, b0, box);
    // Check for intersection
    if ((o1 != o2) && (o3 != o4)) {
        dist = T(0);  // Segments intersect
    }

    return dist;
}

// Function to compute the energy between two segments
real aa_energy(const ArrayXreal& center1, const double& length1, const real& theta1,
             const ArrayXreal& center2, const double& length2, const real& theta2, const std::vector<double>& box,
             const double k_aa, const double kappa_aa) {
    // Calculate the endpoints of the first segment
    real x1_start = center1[0] - (length1 / 2) * cos(theta1);
    real y1_start = center1[1] - (length1 / 2) * sin(theta1);
    real x1_end = center1[0] + (length1 / 2) * cos(theta1);
    real y1_end = center1[1] + (length1 / 2) * sin(theta1);

    // Calculate the endpoints of the second segment
    real x2_start = center2[0] - (length2 / 2) * cos(theta2);
    real y2_start = center2[1] - (length2 / 2) * sin(theta2);
    real x2_end = center2[0] + (length2 / 2) * cos(theta2);
    real y2_end = center2[1] + (length2 / 2) * sin(theta2);

    real a[2] = {x1_start, y1_start};
    real b[2] = {x1_end, y1_end};
    real c[2] = {x2_start, y2_start};
    real d[2] = {x2_end, y2_end};

    // Compute distance between segments
    real dist = segment_segment_distance(a, b, c, d, box);
    dist = dist - 0.03;
    real angle = theta1 - theta2;
    angle_wrap(angle);
    angle = M_PI-abs(angle);
    real approx_force = k_aa * dist;
    return 0.5 * (k_aa * dist * dist + kappa_aa * angle * angle); // Energy
}


// Function to compute the energy between two segments
real am_energy1(const ArrayXreal& center1, const double& length1, const real& theta1,
             const ArrayXreal& center2, const double& length2, const real& theta2, const std::vector<double>& box,
             const double k_am, const double kappa_am, const double myosin_radius) {
    // Calculate the endpoints of the first segment
    real x1_start = center1[0] - (length1 / 2) * cos(theta1);
    real y1_start = center1[1] - (length1 / 2) * sin(theta1);
    real x1_end = center1[0] + (length1 / 2) * cos(theta1);
    real y1_end = center1[1] + (length1 / 2) * sin(theta1);

    // Calculate the endpoints of the second segment
    real x2_start = center2[0] - (length2 / 2) * cos(theta2);
    real y2_start = center2[1] - (length2 / 2) * sin(theta2);
    real x2_end = center2[0] + (length2 / 2) * cos(theta2);
    real y2_end = center2[1] + (length2 / 2) * sin(theta2);

    real a[2] = {x1_start, y1_start};
    real b[2] = {x1_end, y1_end};
    real c[2] = {x2_start, y2_start};
    real d[2] = {x2_end, y2_end};
    // Compute distance between segments
    real dist = segment_segment_distance(a, b, c, d, box);
    real strength = abs(cos(theta1 - theta2));
    real angle = theta1 - theta2;
    angle_wrap(angle);
    angle = abs(angle);
    angle = (angle < M_PI - angle) ? angle : (M_PI - angle);
    if (dist.val()>0.8*myosin_radius){
        if (dist.val()>myosin_radius){
            printf("dist: %f\n",dist.val());
            printf("endpoints after: %f %f %f %f %f %f %f %f\n",val(a[0]),val(a[1]),val(b[0]),val(b[1]),val(c[0]),val(c[1]),val(d[0]),val(d[1]));
            printf("something is wrong\n");
            exit(1);
        }
        real energy1 = 0.5 * k_am* val(strength) * dist * dist;
        real energy2 = 0.5 * kappa_am * angle * angle;

        return 0.5 * k_am* val(strength) * dist * dist + 0.5 * kappa_am * angle * angle; // Energy
    }
    else {
        real energy = 0.5 * kappa_am * angle * angle;
        return 0.5 * kappa_am * angle * angle;
    }
}


// Function to compute the energy between two segments
real am_energy(const real& theta1,const real& theta2, const double kappa_am) {
    real angle = theta1 - theta2;
    angle_wrap(angle);
    angle = abs(angle);
    angle = (angle < M_PI - angle) ? angle : (M_PI - angle);
    return 0.5 * kappa_am * angle * angle; // Energy
}


// Function to compute forces and energy using autodiff
std::vector<double> compute_aa_force_and_energy(Filament& actin,
                                int& actin1_index,  int& actin2_index, 
                               const std::vector<double>& box, const double k_aa, const double kappa_aa) {
    ArrayXreal center1(2);   
    center1 << actin.center[actin1_index].x, actin.center[actin1_index].y;      
    real theta1 = actin.theta[actin1_index];
    ArrayXreal center2(2);
    center2 << actin.center[actin2_index].x, actin.center[actin2_index].y;
    real theta2 = actin.theta[actin2_index];  
    real u;  
    std::vector<double> forces; 
    VectorXd forces_2 = -gradient(aa_energy, wrt(center1, theta1, theta2), 
            at(center1, actin.length, theta1, center2, actin.length, theta2, box, k_aa, kappa_aa), u);
    forces.resize(forces_2.size());
    VectorXd::Map(&forces[0], forces_2.size()) = forces_2;
    return forces;
}

// Function to compute forces and energy using autodiff
std::vector<double> compute_am_force_and_energy(Filament& actin, Myosin& myosin,
                                int& actin_index, int& myosin_index,
                               const std::vector<double>& box, const double k_am, const double kappa_am, const double myosin_radius, const bool turn_on_spring) {
    ArrayXreal center1(2);   
    center1 << actin.center[actin_index].x, actin.center[actin_index].y;      
    real theta1 = actin.theta[actin_index];
    ArrayXreal center2(2);
    center2 << myosin.center[myosin_index].x, myosin.center[myosin_index].y;
    real theta2 = myosin.theta[myosin_index];  
    real u;
    std::vector<double> forces; 
    VectorXd forces_2;    
    if (turn_on_spring){
        forces_2 = -gradient(am_energy1, wrt(center1, theta1, theta2), 
                at(center1, actin.length, theta1, center2, myosin.length, theta2, box, k_am,kappa_am, myosin_radius), u);
        forces.resize(forces_2.size());
        double energy = am_energy1(center1, actin.length, theta1, center2, myosin.length, theta2, box, k_am,kappa_am, myosin_radius).val();
        VectorXd::Map(&forces[0], forces_2.size()) = forces_2;
    }
    else{
        forces_2 = -gradient(am_energy, wrt(theta1, theta2), 
                at(theta1,theta2, kappa_am), u);
        //append zeros for center1 in the front
        forces.resize(forces_2.size()+2);
        forces[0] = 0;
        forces[1] = 0;
        VectorXd::Map(&forces[2], forces_2.size()) = forces_2;    
    }

    return forces;
}
#endif
