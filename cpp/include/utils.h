#ifndef UTILS_H
#define UTILS_H
#include <cmath>

namespace utils {

double norm(double * x, int dim){
    double norm = 0;
    for(int i = 0; i < dim; i++){
        norm += x[i]*x[i];
    }
    return sqrt(norm);
}

void pbc_wrap(double * x, double * box){
    x[0] = x[0] - box[0]*round(x[0]/box[0]);
    x[1] = x[1] - box[1]*round(x[1]/box[1]);
}

void angle_wrap(double& theta){
    theta = theta - 2*M_PI*round(theta/(2*M_PI));
}

double point_point_distance(double * x, double * y, double * box){
    double r[2];
    r[0] = x[0] - y[0];
    r[1] = x[1] - y[1];
    pbc_wrap(r, box);
    return norm(r, 2);
}

double point_segment_distance(double * x, double * a, double * b, double * box){
    // compute the shortest distance between a point x and a segment ab
    double ab[2];
    ab[0] = b[0] - a[0];
    ab[1] = b[1] - a[1];
    double ab_norm = norm(ab, 2);
    ab[0] = ab[0]/ab_norm;
    ab[1] = ab[1]/ab_norm;
    double ap[2];
    ap[0] = x[0] - a[0];
    ap[1] = x[1] - a[1];
    pbc_wrap(ap, box);
    double ap_ab = ap[0]*ab[0] + ap[1]*ab[1];
    ap_ab = fmax(0, fmin(ap_ab, ab_norm));
    double vec[2];
    vec[0] = ap[0] - ap_ab*ab[0];
    vec[1] = ap[1] - ap_ab*ab[1];
    double dist = norm(vec, 2);
    return dist;
}

double distance(double* p1, double* p2) {
    return std::hypot(p2[0] - p1[0], p2[1] - p1[1]);
}

// Function to calculate the distance between a point and a line segment 
double point_segment_distance_new(double* p, double* p1, double* p2, double* box) {
    double l2 = (p2[0] - p1[0]) * (p2[0] - p1[0]) + (p2[1] - p1[1]) * (p2[1] - p1[1]);
    if (l2 == 0) return distance(p, p1);  // Segment is actually a point
    double t = ((p[0] - p1[0]) * (p2[0] - p1[0]) + (p[1] - p1[1]) * (p2[1] - p1[1])) / l2;
    if (t < 0) return distance(p, p1);  // Beyond the 'p1' end of the segment
    else if (t > 1) return distance(p, p2);  // Beyond the 'p2' end of the segment
    double projection_x = p1[0] + t * (p2[0] - p1[0]);
    double projection_y = p1[1] + t * (p2[1] - p1[1]);
    return distance(p, new double[2]{projection_x, projection_y});
}

int orientation(double * p, double * q, double * r){
    double val = (q[1] - p[1]) * (r[0] - q[0]) - (q[0] - p[0]) * (r[1] - q[1]);
    if(val == 0){
        return 0;
    }else if(val > 0){
        return 1;
    }else{
        return 2;
    }
}

double segment_segment_distance(double * a, double * b, double * c, double * d, double * box){
    double a_cd = point_segment_distance(a, c, d, box);
    double b_cd = point_segment_distance(b, c, d, box);
    double c_ab = point_segment_distance(c, a, b, box);
    double d_ab = point_segment_distance(d, a, b, box);
    double dist = fmin(fmin(a_cd, b_cd), fmin(c_ab, d_ab));
    // check if segments intersect
    int o1 = orientation(a, c, b);
    int o2 = orientation(a, c, d);
    int o3 = orientation(b, d, a);
    int o4 = orientation(b, d, c);
    if((o1 != o2) && (o3 != o4)){
        dist = 0;
    }
    return dist;
}
}


#endif 