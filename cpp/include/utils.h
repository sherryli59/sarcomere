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
    //defining the * operator for the point struct
    vec operator*(double c) const
    {
        vec scaled;
        scaled.x = c * x;
        scaled.y = c * y;
        return scaled;
    }
    friend vec operator*(double c, const vec& v) {
        return v * c;  // Reuse the vec * double multiplication defined above
    }
    //defining the / operator for the point struct
    vec operator/(double c) const
    {
        vec scaled;
        scaled.x = x / c;
        scaled.y = y / c;
        return scaled;
    }
    //define the norm of the point
    double norm() const
    {
        return sqrt(x*x + y*y);
    }

    vec operator-() const
    {
        vec neg;
        neg.x = -x;
        neg.y = -y;
        return neg;
    }

    void pbc_wrap(std::vector <double> box){
        x = x - box[0]*round(x/box[0]);
        y = y - box[1]*round(y/box[1]);
    }

    double distance(const vec& p) const
    {
        return std::sqrt((x - p.x) * (x - p.x) + (y - p.y) * (y - p.y));
    }

    double distance(const vec& p, std::vector<double> box) const
    {
        double dx = x - p.x;
        double dy = y - p.y;
        dx = dx - box[0]*round(dx/box[0]);
        dy = dy - box[1]*round(dy/box[1]);
        return std::sqrt(dx*dx + dy*dy);
    }

    double dot(const vec& p) const
    {
        return x * p.x + y * p.y;
    }
    vec normalize() const
    {
        double n = norm();
        if (n == 0)
        {
            return *this;
        }
        return *this / norm();
    }

    double norm_squared() const
    {
        return x*x + y*y;
    }
    


};




double pbc_wrap(double x, double& box){
    return x - box*round(x/box);
}

void angle_wrap(double& theta){
    theta = theta - 2*M_PI*round(theta/(2*M_PI));
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
              [&vec](size_t a, size_t b) { return vec[a] > vec[b]; });

    return indices;
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
    pbc_wrap(left_vec[0], box[0]);
    pbc_wrap(left_vec[1], box[1]);
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