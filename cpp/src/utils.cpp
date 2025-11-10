#include "utils.h"
#include <algorithm>  // For std::find
#include <cmath>      // For std::sqrt, std::round, and M_PI

namespace utils {

//------------------------------------------------------------------------------
// Definitions of free functions
//------------------------------------------------------------------------------


bool compare_indices(const std::vector<int>& a, const std::vector<int>& b) {
    if (a.size() != b.size()) {
        return false;
    }
    std::vector<int> sorted_a = a;
    std::vector<int> sorted_b = b;
    std::sort(sorted_a.begin(), sorted_a.end());
    std::sort(sorted_b.begin(), sorted_b.end());
    return std::equal(sorted_a.begin(), sorted_a.end(), sorted_b.begin());
}

double pbc_wrap(double x, double& box) {
    return wrap_axis(x, box, true);
}

vec pbc_diff_masked(const vec& a, const vec& b, const std::vector<double>& box, const std::array<bool,3>& periodic) {
    double Lx = box.size() > 0 ? box[0] : 0.0;
    double Ly = box.size() > 1 ? box[1] : 0.0;
    double Lz = box.size() > 2 ? box[2] : 0.0;
    return vec {
        wrap_axis(a.x - b.x, Lx, periodic[0]),
        wrap_axis(a.y - b.y, Ly, periodic[1]),
        wrap_axis(a.z - b.z, Lz, periodic[2])
    };
}

void pbc_wrap_centered(vec& x, const std::vector<double>& box, const std::array<bool,3>& periodic) {
    auto wrap_centered_axis = [](double coord, double Lk) {
        double q = coord;
        double H = 0.5 * Lk;
        q -= Lk * std::round(q / Lk);
        if (q <= -H) {
            q += Lk;
        }
        return q;
    };

    const double Lx = box.size() > 0 ? box[0] : 0.0;
    const double Ly = box.size() > 1 ? box[1] : 0.0;
    const double Lz = box.size() > 2 ? box[2] : 0.0;

    if (periodic[0] && Lx > 0.0) {
        x.x = wrap_centered_axis(x.x, Lx);
    }
    if (periodic[1] && Ly > 0.0) {
        x.y = wrap_centered_axis(x.y, Ly);
    }
    if (periodic[2] && Lz > 0.0) {
        x.z = wrap_centered_axis(x.z, Lz);
    }
}

void clamp_orientation_for_box_centered(vec& u, double h, double r,
                                        const std::vector<double>& box,
                                        const std::array<bool,3>& periodic) {
    auto box_len = [&](int idx) -> double {
        return (static_cast<size_t>(idx) < box.size()) ? box[idx] : 0.0;
    };

    std::array<double,3> cap{1.0, 1.0, 1.0};
    for (int k = 0; k < 3; ++k) {
        const double Lk = box_len(k);
        if (!periodic[k] && Lk > 0.0) {
            const double H = 0.5 * Lk;
            const double ck = (H - r) / std::max(h, 1e-12);
            cap[k] = std::clamp(ck, 0.0, 1.0);
        }
    }

    if (std::abs(u.x) <= cap[0] && std::abs(u.y) <= cap[1] && std::abs(u.z) <= cap[2]) {
        return;
    }

    auto clip = [](double value, double limit) {
        return std::copysign(std::min(std::abs(value), limit), value);
    };

    u.x = clip(u.x, cap[0]);
    u.y = clip(u.y, cap[1]);
    u.z = clip(u.z, cap[2]);

    double n2 = u.x * u.x + u.y * u.y + u.z * u.z;
    if (n2 < 1e-16) {
        int kbest = 0;
        double best_len = box_len(kbest);
        for (int k = 1; k < 3; ++k) {
            const double Lk = box_len(k);
            if (Lk > best_len) {
                kbest = k;
                best_len = Lk;
            }
        }
        u = vec{
            (kbest == 0) ? 1.0 : 0.0,
            (kbest == 1) ? 1.0 : 0.0,
            (kbest == 2) ? 1.0 : 0.0
        };
    } else {
        double n = std::sqrt(n2);
        u.x /= n;
        u.y /= n;
        u.z /= n;
    }
}

void no_flux_slide_capsule_centered(vec& c, vec& v, const vec& u,
                                    double h, double r,
                                    const std::vector<double>& box,
                                    const std::array<bool,3>& periodic) {
    auto box_len = [&](int idx) -> double {
        return (static_cast<size_t>(idx) < box.size()) ? box[idx] : 0.0;
    };

    for (int k = 0; k < 3; ++k) {
        const double Lk = box_len(k);
        if (periodic[k] || Lk <= 0.0) {
            continue;
        }

        const double H = 0.5 * Lk;
        const double uk = (k == 0) ? u.x : (k == 1 ? u.y : u.z);
        const double pad = h * std::abs(uk) + r;
        const double lo = -(H - pad);
        const double hi = (H - pad);

        double* ck = (k == 0) ? &c.x : (k == 1 ? &c.y : &c.z);
        double* vk = (k == 0) ? &v.x : (k == 1 ? &v.y : &v.z);

        if (lo > hi) {
            *ck = 0.0;
            *vk = 0.0;
            continue;
        }

        double before = *ck;
        if (*ck < lo) {
            *ck = lo;
        } else if (*ck > hi) {
            *ck = hi;
        }

        if (*ck != before) {
            *vk = 0.0;
        }
    }
}


double angle_between(const vec& u1, const vec& u2) {
    double dot = u1.dot(u2);
    dot = std::clamp(dot, -1.0, 1.0);
    return std::acos(dot);
}


std::ostream& operator<<(std::ostream &os, const vec &v) {
    os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    return os;
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
