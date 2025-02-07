#include "utils.h"
#include <algorithm>  // For std::find
#include <cmath>      // For std::sqrt, std::round, and M_PI

namespace utils {

//------------------------------------------------------------------------------
// Definitions of free functions
//------------------------------------------------------------------------------

double pbc_wrap(double x, double& box) {
    return x - box * std::round(x / box);
}

void angle_wrap(double& theta) {
    theta = theta - 2 * M_PI * std::round(theta / (2 * M_PI));
}

double dotProduct(double x1, double y1, double x2, double y2) {
    return x1 * x2 + y1 * y2;
}

double analyze_overlap(double* P1_left, double* P1_right, 
                         double* P2_left, double* P2_right, double d, 
                         std::vector<double> box) {
    // Direction vector of segment 1.
    double dir1_x = P1_right[0] - P1_left[0];
    double dir1_y = P1_right[1] - P1_left[1];

    // Direction vector of segment 2.
    double dir2_x = P2_right[0] - P2_left[0];
    double dir2_y = P2_right[1] - P2_left[1];

    // Wrap the directional components.
    pbc_wrap(dir1_x, box[0]);
    pbc_wrap(dir1_y, box[1]);
    pbc_wrap(dir2_x, box[0]);
    pbc_wrap(dir2_y, box[1]);

    // Compute the normal vector (perpendicular) of segment 2.
    double nx2 = -dir2_y;
    double ny2 = dir2_x;

    // Normalize the normal vector.
    double norm_factor = std::sqrt(nx2 * nx2 + ny2 * ny2);
    nx2 /= norm_factor;
    ny2 /= norm_factor;

    // Compute the vector from the left endpoint of segment 2 to that of segment 1.
    std::vector<double> left_vec = { P1_left[0] - P2_left[0], P1_left[1] - P2_left[1] };
    pbc_wrap(left_vec[0], box[0]);
    pbc_wrap(left_vec[1], box[1]);

    // Project the left_vec and segment 1 direction onto the normal.
    double dot1 = dotProduct(left_vec[0], left_vec[1], nx2, ny2);
    double dot2 = dotProduct(dir1_x, dir1_y, nx2, ny2);

    // Calculate the parameter values along segment 1.
    double t1 = (d - dot1) / dot2;
    double t2 = (-d - dot1) / dot2;

    // Clamp t1 and t2 to the [0, 1] range.
    t1 = std::max(0.0, std::min(1.0, t1));
    t2 = std::max(0.0, std::min(1.0, t2));

    double ratio = std::abs(t1 - t2);
    return ratio;
}

//------------------------------------------------------------------------------
// Definitions of MoleculeConnection member functions
//------------------------------------------------------------------------------

MoleculeConnection::MoleculeConnection() {
    // Default constructor.
}

MoleculeConnection::MoleculeConnection(int numA) : connections(numA) {
    // Initialize each molecule A's connection vector.
    for (int i = 0; i < numA; i++) {
        connections[i] = std::vector<int>();
    }
}

void MoleculeConnection::addConnection(int aIndex, int bIndex) {
    if (aIndex >= 0 && aIndex < static_cast<int>(connections.size())) {
        // Avoid adding duplicate connections.
        if (std::find(connections[aIndex].begin(), connections[aIndex].end(), bIndex) == connections[aIndex].end()) {
            connections[aIndex].push_back(bIndex);
        }
    }
}

void MoleculeConnection::deleteConnection(int aIndex, int bIndex) {
    if (aIndex >= 0 && aIndex < static_cast<int>(connections.size())) {
        auto& connList = connections[aIndex];
        auto it = std::find(connList.begin(), connList.end(), bIndex);
        if (it != connList.end()) {
            connList.erase(it);
        }
    }
}

void MoleculeConnection::deleteAllConnections(int aIndex) {
    if (aIndex >= 0 && aIndex < static_cast<int>(connections.size())) {
        connections[aIndex].clear();
    }
}

const std::vector<int>& MoleculeConnection::getConnections(int aIndex) const {
    if (aIndex >= 0 && aIndex < static_cast<int>(connections.size())) {
        return connections[aIndex];
    } else {
        static const std::vector<int> emptyList;
        return emptyList;
    }
}

} // namespace utils
