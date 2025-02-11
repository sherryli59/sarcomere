#ifndef UTILS_H
#define UTILS_H
#include <cmath>
#include <vector>
#include <map>
#include <algorithm> 
#include <numeric>
#include <utility>



namespace utils {

struct vec
{
    double x;
    double y;
    //defining the + operator for the point struct
    vec operator+(const vec& p) const
    {
        vec sum;
        sum.x = x + p.x;
        sum.y = y + p.y;
        return sum;
    }
    //defining the - operator for the point struct
    vec operator-(const vec& p) const
    {
        vec diff;
        diff.x = x - p.x;
        diff.y = y - p.y;
        return diff;
    }
    //define the norm of the point
    double norm() const
    {
        return sqrt(x*x + y*y);
    }

    void pbc_wrap(std::vector <double> box){
        x = x - box[0]*round(x/box[0]);
        y = y - box[1]*round(y/box[1]);
    }
    
};

double norm(double * x, int dim){
    double norm = 0;
    for(int i = 0; i < dim; i++){
        norm += x[i]*x[i];
    }
    return sqrt(norm);
}

double norm(std::vector <double> x){
    double norm = 0;
    for(int i = 0; i < x.size(); i++){
        norm += x[i]*x[i];
    }
    return sqrt(norm);
}

void pbc_wrap(double * x, std::vector <double> box){
    x[0] = x[0] - box[0]*round(x[0]/box[0]);
    x[1] = x[1] - box[1]*round(x[1]/box[1]);
}

void pbc_wrap(std::vector <double> & x, std::vector <double>& box){
    x[0] = x[0] - box[0]*round(x[0]/box[0]);
    x[1] = x[1] - box[1]*round(x[1]/box[1]);
}


double pbc_wrap(double x, double& box){
    return x - box*round(x/box);
}

void angle_wrap(double& theta){
    theta = theta - 2*M_PI*round(theta/(2*M_PI));
}

double point_point_distance(double * x, double * y, std::vector <double> box){
    double r[2];
    r[0] = x[0] - y[0];
    r[1] = x[1] - y[1];
    pbc_wrap(r, box);
    return norm(r, 2);
}

double point_segment_distance(double * x, double * a, double * b, std::vector <double> box){
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
    std::vector<double> vec(2);
    vec[0] = ap[0] - ap_ab*ab[0];
    vec[1] = ap[1] - ap_ab*ab[1];
    double dist = norm(vec);
    return dist;
}

std::pair<double,std::vector<double>> point_segment_distance_w_normal(double * x, double * a, double * b, std::vector <double> box){
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
    std::vector<double> vec(2);
    vec[0] = ap_ab*ab[0] - ap[0];
    vec[1] = ap_ab*ab[1] - ap[1];
    double dist = norm(vec);
    return std::make_pair(dist, vec);
}

double distance(double* p1, double* p2) {
    return std::hypot(p2[0] - p1[0], p2[1] - p1[1]);
}

// Function to calculate the distance between a point and a line segment 
double point_segment_distance_new(double* p, double* p1, double* p2, std::vector <double> box) {
    double l2 = (p2[0] - p1[0]) * (p2[0] - p1[0]) + (p2[1] - p1[1]) * (p2[1] - p1[1]);
    if (l2 == 0) return distance(p, p1);  // Segment is actually a point
    double t = ((p[0] - p1[0]) * (p2[0] - p1[0]) + (p[1] - p1[1]) * (p2[1] - p1[1])) / l2;
    if (t < 0) return distance(p, p1);  // Beyond the 'p1' end of the segment
    else if (t > 1) return distance(p, p2);  // Beyond the 'p2' end of the segment
    double projection_x = p1[0] + t * (p2[0] - p1[0]);
    double projection_y = p1[1] + t * (p2[1] - p1[1]);
    return distance(p, new double[2]{projection_x, projection_y});
}

int orientation(double * p, double * q, double * r, std::vector <double> box){
    
    double val = pbc_wrap(q[1] - p[1],box[1]) * pbc_wrap(r[0] - q[0], box[0]) - pbc_wrap(q[0] - p[0],box[0]) * pbc_wrap(r[1] - q[1],box[1]);
    if(val == 0){
        return 0;
    }else if(val > 0){
        return 1;
    }else{
        return 2;
    }
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


double segment_segment_distance(double * a, double * b, double * c, double * d, std::vector <double> box){
    //segments are ab and cd
    double a_cd = point_segment_distance(a, c, d, box);
    double b_cd = point_segment_distance(b, c, d, box);
    double c_ab = point_segment_distance(c, a, b, box);
    double d_ab = point_segment_distance(d, a, b, box);
    double dist = fmin(fmin(a_cd, b_cd), fmin(c_ab, d_ab));
    // check if segments intersect
    double ab_center[2];
    ab_center[0] = (a[0] + b[0])/2;
    ab_center[1] = (a[1] + b[1])/2;
    double cd_center[2];
    cd_center[0] = (c[0] + d[0])/2;
    cd_center[1] = (c[1] + d[1])/2;
    double displacement[2];
    displacement[0] = box[0]*round((ab_center[0] - cd_center[0])/box[0]);
    displacement[1] = box[1]*round((ab_center[1] - cd_center[1])/box[1]);
    double a0[2], b0[2];
    a0[0]=a[0]-displacement[0];
    a0[1]=a[1]-displacement[1];
    b0[0]=b[0]-displacement[0];
    b0[1]=b[1]-displacement[1];
    int o1 = orientation(a0, b0, c);
    int o2 = orientation(a0, b0, d);
    int o3 = orientation(c, d, a0);
    int o4 = orientation(c, d, b0);
    if((o1 != o2) && (o3 != o4)){
        dist = 0;
    }
    return dist;
}

std::pair<double,std::map<std::string, std::vector<double>>> segment_segment_distance_w_normal(double * a, double * b, double * c, double * d, std::vector <double> box){
    //segments are ab and cd
    std::pair<double, std::vector<double>> a_cd = point_segment_distance_w_normal(a, c, d, box);
    double a_cd_dist = a_cd.first;
    std::pair<double, std::vector<double>> b_cd = point_segment_distance_w_normal(b, c, d, box);
    double b_cd_dist = b_cd.first;
    std::pair<double, std::vector<double>> c_ab = point_segment_distance_w_normal(c, a, b, box);
    double c_ab_dist = c_ab.first;
    std::pair<double, std::vector<double>> d_ab = point_segment_distance_w_normal(d, a, b, box);
    double d_ab_dist = d_ab.first;
    double dist = fmin(fmin(a_cd_dist, b_cd_dist), fmin(c_ab_dist, d_ab_dist));
    std::map<std::string, std::vector<double>> normal_dict;
    std::vector<double> vec(2);
    if(dist == a_cd_dist){
        vec[0] = a[0];
        vec[1] = a[1];
        normal_dict["start"] = vec;
        normal_dict["vector"] = a_cd.second;
        vec[0] += a_cd.second[0];
        vec[1] += a_cd.second[1];
        normal_dict["end"] = vec;
    }else if(dist == b_cd_dist){
        vec[0] = b[0];
        vec[1] = b[1];
        normal_dict["start"] = vec;
        normal_dict["vector"] = b_cd.second;
        vec[0] += b_cd.second[0];
        vec[1] += b_cd.second[1];
        normal_dict["end"] = vec;

    }else if(dist == c_ab_dist){
        vec[0] = c[0];
        vec[1] = c[1];
        normal_dict["end"] = vec;
        normal_dict["vector"] = c_ab.second;
        normal_dict["vector"][0] = -normal_dict["vector"][0];
        normal_dict["vector"][1] = -normal_dict["vector"][1];
        vec[0] += c_ab.second[0];
        vec[1] += c_ab.second[1];
        normal_dict["start"] = vec;
    }else{
        vec[0] = d[0];
        vec[1] = d[1];
        normal_dict["end"] = vec;
        normal_dict["vector"] = d_ab.second;
        normal_dict["vector"][0] = -normal_dict["vector"][0];
        normal_dict["vector"][1] = -normal_dict["vector"][1];
        vec[0] += d_ab.second[0];
        vec[1] += d_ab.second[1];
        normal_dict["start"] = vec;
    }
    // check if segments intersect
    double ab_center[2];
    ab_center[0] = (a[0] + b[0])/2;
    ab_center[1] = (a[1] + b[1])/2;
    double cd_center[2];
    cd_center[0] = (c[0] + d[0])/2;
    cd_center[1] = (c[1] + d[1])/2;
    double displacement[2];
    displacement[0] = box[0]*round((ab_center[0] - cd_center[0])/box[0]);
    displacement[1] = box[1]*round((ab_center[1] - cd_center[1])/box[1]);
    double a0[2], b0[2];
    a0[0]=a[0]-displacement[0];
    a0[1]=a[1]-displacement[1];
    b0[0]=b[0]-displacement[0];
    b0[1]=b[1]-displacement[1];
    int o1 = orientation(a0, b0, c);
    int o2 = orientation(a0, b0, d);
    int o3 = orientation(c, d, a0);
    int o4 = orientation(c, d, b0);
    if((o1 != o2) && (o3 != o4)){
        dist = 0;
        normal_dict["vector"] = {0,0};
    }
    return std::make_pair(dist, normal_dict);
}

class MoleculeConnection {
public:
    std::vector<std::vector<int>> connections;
    MoleculeConnection(){
    }
    MoleculeConnection(int numA) : connections(numA) {
        // Initialize the vector for each A molecule
        for (int i = 0; i < numA; i++) {
            connections[i] = std::vector<int>();
        }
    }

    // Add a co nnection from molecule A to molecule B
    void addConnection(int aIndex, int bIndex) {
        if (aIndex >= 0 && aIndex < connections.size()) {
            // Avoid adding duplicates (optional)
            if (std::find(connections[aIndex].begin(), connections[aIndex].end(), bIndex) == connections[aIndex].end()) {
                connections[aIndex].push_back(bIndex);
            }
        }
    }

    // Delete a connection from molecule A to molecule B
    void deleteConnection(int aIndex, int bIndex) {
        if (aIndex >= 0 && aIndex < connections.size()) {
            auto& connList = connections[aIndex];
            auto it = std::find(connList.begin(), connList.end(), bIndex);
            if (it != connList.end()) {
                connList.erase(it); // Remove the connection
            } 
        }
    }

    void deleteAllConnections(int aIndex) {
        if (aIndex >= 0 && aIndex < connections.size()) {
            connections[aIndex].clear();
        }
    }

    // Query connections for a specific molecule A
    const std::vector<int>& getConnections(int aIndex) const {
        if (aIndex >= 0 && aIndex < connections.size()) {
            return connections[aIndex];
        }
        else {
            static const std::vector<int> emptyList;
            return emptyList;
        }
    }

};

template <typename T>
std::vector<size_t> sort_indices(const std::vector<T>& vec) {
    // Initialize indices with 0, 1, ..., n-1
    std::vector<size_t> indices(vec.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        indices[i] = i;
    }

    // Sort indices based on values in the vector
    std::sort(indices.begin(), indices.end(),
              [&vec](size_t a, size_t b) { return vec[a] < vec[b]; });

    return indices;
}

// Function to check for intersection and find the intersection point of two line segments, considering PBC
bool getIntersection(const double* A, const double* B, const double* C, const double* D, double* intersection, const std::vector<double>& box) {
    // Adjust the start point A based on PBC (shifted toward B)
    double AB[2] = {A[0] - B[0], A[1] - B[1]};
    pbc_wrap(AB, box);
    double A_adj[2] = {B[0] + AB[0], B[1] + AB[1]};  // Adjusted A with PBC

    // Line A_adjB represented as a1*x + b1*y = c1
    double a1 = B[1] - A_adj[1];
    double b1 = A_adj[0] - B[0];
    double c1 = a1 * A_adj[0] + b1 * A_adj[1];

    // Adjust the start point C based on PBC (shifted toward D)
    double CD[2] = {C[0] - D[0], C[1] - D[1]};
    pbc_wrap(CD, box);
    double C_adj[2] = {D[0] + CD[0], D[1] + CD[1]};  // Adjusted C with PBC

    // Line C_adjD represented as a2*x + b2*y = c2
    double a2 = D[1] - C_adj[1];
    double b2 = C_adj[0] - D[0];
    double c2 = a2 * C_adj[0] + b2 * C_adj[1];

    // Compute the determinant to check if the lines are parallel
    double determinant = a1 * b2 - a2 * b1;
    printf("A,B,C,D: (%f, %f), (%f, %f), (%f, %f), (%f, %f)\n", A_adj[0], A_adj[1], B[0], B[1], C_adj[0], C_adj[1], D[0], D[1]);
    if (determinant == 0) {
        // The lines are parallel, no intersection
        return false;
    } else {
        // Intersection point exists
        intersection[0] = (b2 * c1 - b1 * c2) / determinant;
        intersection[1] = (a1 * c2 - a2 * c1) / determinant;
        printf("Intersection point: (%f, %f)\n", intersection[0], intersection[1]);
        // Check if the intersection point is within both segments
        if (intersection[0] >= std::min(A_adj[0], B[0]) && intersection[0] <= std::max(A_adj[0], B[0]) &&
            intersection[1] >= std::min(A_adj[1], B[1]) && intersection[1] <= std::max(A_adj[1], B[1]) &&
            intersection[0] >= std::min(C_adj[0], D[0]) && intersection[0] <= std::max(C_adj[0], D[0]) &&
            intersection[1] >= std::min(C_adj[1], D[1]) && intersection[1] <= std::max(C_adj[1], D[1])) {
            return true;  // The intersection is within both segments
        }
    }
    return false;  // The intersection is outside the segments
}

// Function to rotate a vector (x, y) by an angle theta (in radians)
void rotate_vector(double* vec, double cos_theta, double sin_theta) {
    double x_new = cos_theta * vec[0] - sin_theta * vec[1];
    double y_new = sin_theta * vec[0] + cos_theta * vec[1];
    vec[0] = x_new;
    vec[1] = y_new;
}

bool check_overlap(double * start1, double * end1, double * start2, double * end2, std::vector <double> box){
    // adjust end1 and end2 based on PBC
    double vec1[2] = {end1[0] - start1[0], end1[1] - start1[1]};
    pbc_wrap(vec1, box);
    double end1_adj[2] = {start1[0] + vec1[0], start1[1] + vec1[1]};  // Adjusted end1 with PBC
    double vec2[2] = {end2[0] - start2[0], end2[1] - start2[1]};
    pbc_wrap(vec2, box);
    double end2_adj[2] = {start2[0] + vec2[0], start2[1] + vec2[1]};  // Adjusted end2 with PBC    //rotate and shift the whole system so that start1 is at the origin and end1 is on the positive x-axis
    double start2_adj[2] = {start2[0] - start1[0], start2[1] - start1[1]};
    end2_adj[0] = end2_adj[0] - start1[0];
    end2_adj[1] = end2_adj[1] - start1[1];
    double cos_theta = vec1[0]/sqrt(vec1[0]*vec1[0] + vec1[1]*vec1[1]);
    double sin_theta = vec1[1]/sqrt(vec1[0]*vec1[0] + vec1[1]*vec1[1]);
    rotate_vector(start2_adj, cos_theta, -sin_theta);
    rotate_vector(end2_adj, cos_theta, -sin_theta);
    //min image convention    
    end2_adj[0] = end2_adj[0] - round(start2_adj[0]/box[0])*box[0];
    end2_adj[1] = end2_adj[1] - round(start2_adj[1]/box[1])*box[1];
    start2_adj[0] = start2_adj[0] - round(start2_adj[0]/box[0])*box[0];
    start2_adj[1] = start2_adj[1] - round(start2_adj[1]/box[1])*box[1];
    if (end2_adj[0] < 0){
        return true;
    }
    return false;
}

bool check_cb(double * actin_left1, double * actin_right1, double * myosin_left1, double * myosin_right1,
            double * actin_left2, double * actin_right2, double * myosin_left2, double * myosin_right2,
     std::vector <double> box){
    double dist1 = point_segment_distance(actin_left1, myosin_left1, myosin_right1, box);
    double dist2 = point_segment_distance(actin_right1, myosin_left1, myosin_right1, box);
    if (dist1<dist2){
        return false;
    }
    dist1 = point_segment_distance(actin_left2, myosin_left2, myosin_right2, box);
    dist2 = point_segment_distance(actin_right2, myosin_left2, myosin_right2, box);
    if (dist1<dist2){
        return false;
    }
    return check_overlap(actin_left1, actin_right1, actin_left2, actin_right2, box);
}

// Function to compute the dot product of two vectors
double dotProduct(double x1, double y1, double x2, double y2) {
    return x1 * x2 + y1 * y2;
}

// Function to calculate the points on segment 1 that are distance d away from segment 2
double analyze_overlap(double* P1_left, double* P1_right, 
                        double* P2_left, double* P2_right, double d, std::vector<double> box) {
                                           
    // Direction vector of segment 1
    double dir1_x = P1_right[0] - P1_left[0];
    double dir1_y = P1_right[1] - P1_left[1];

    // Direction vector of segment 2
    double dir2_x = P2_right[0] - P2_left[0];
    double dir2_y = P2_right[1] - P2_left[1];
    pbc_wrap(dir1_x, box[0]);
    pbc_wrap(dir1_y, box[1]);
    pbc_wrap(dir2_x, box[0]);
    pbc_wrap(dir2_y, box[1]);

    // Normal vector of segment 2 (perpendicular to the direction vector)
    double nx2 = -dir2_y;
    double ny2 = dir2_x;

    // Normalize the normal vector of segment 2
    double norm_factor = std::sqrt(nx2 * nx2 + ny2 * ny2);
    nx2 /= norm_factor;
    ny2 /= norm_factor;
    std::vector<double> left_vec = {P1_left[0] - P2_left[0], P1_left[1] - P2_left[1]};
    pbc_wrap(left_vec, box);
    double dot1 = dotProduct(left_vec[0], left_vec[1], nx2, ny2);
    double dot2 = dotProduct(dir1_x, dir1_y, nx2, ny2);
    double t1 = (d - dot1) / dot2;
    double t2 = (-d - dot1) / dot2;
    // Clamp t1 and t2 between 0 and 1 to ensure the points are within the segment
    t1 = std::max(0.0, std::min(1.0, t1));
    t2 = std::max(0.0, std::min(1.0, t2));
    
    double ratio = std::abs(t1 - t2);
    // // Points on segment 1
    // double* point1 = new double[2];
    // double* point2 = new double[2];

    // point1[0] = P1_left[0] + t1 * dir1_x;
    // point1[1] = P1_left[1] + t1 * dir1_y;

    // point2[0] = P1_left[0] + t2 * dir1_x;
    // point2[1] = P1_left[1] + t2 * dir1_y;
    return ratio;
}

}



#endif 