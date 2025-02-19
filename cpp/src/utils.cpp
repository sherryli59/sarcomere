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

// Wrap an angle into the interval [0, π) if is_range_pi is true, or [–π, π) otherwise.
void angle_wrap(double& theta, bool& is_range_pi){
    if (is_range_pi) {
        theta = theta - M_PI * std::round(theta / M_PI);
    } else {
        theta = theta - 2 * M_PI * std::round(theta / (2 * M_PI));
    }
}

// Overload += operator for std::vector<double>
std::vector<double>& operator+=(std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vectors must be of the same size for element-wise addition.");
    }

    for (size_t i = 0; i < a.size(); ++i) {
        a[i] += b[i];  // Element-wise addition
    }
    return a;
}

// Overload + operator for std::vector<double>
std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> result = a; // Make a copy
    result += b; // Use the overloaded += operator
    return result;
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
