#include "components.h"
#include "interaction.h"

#include <Eigen/Eigenvalues>
#include <cassert>
#include <cmath>
#include <gsl/gsl_rng.h>

namespace {
constexpr int count = 6;

Eigen::VectorXd transverse_force(
    Filament& actin, Myosin& myosin, const std::vector<double>& box,
    double k_am, double k_aa, double d_am, double d_aa) {
    Eigen::VectorXd force = Eigen::VectorXd::Zero(2 * count);
    actin.update_endpoints();
    myosin.update_endpoints();
    for (int i = 0; i < count; ++i) {
        int m = 0;
        const auto am = compute_am_force_and_energy(
            actin, myosin, i, m, box, k_am, 1.0, d_am + 0.05, d_am);
        force[2 * i] += am[1];
        force[2 * i + 1] += am[2];
    }
    for (int i = 0; i < count; ++i) {
        int j = (i + 1) % count;
        const auto aa = compute_aa_force_and_energy(
            actin, i, j, box, k_aa, 1.0, d_aa + 0.02, d_aa);
        force[2 * i] += aa[0 + 1];
        force[2 * i + 1] += aa[0 + 2];
        force[2 * j] -= aa[0 + 1];
        force[2 * j + 1] -= aa[0 + 2];
    }
    return force;
}
}

int main() {
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    const std::vector<double> box{4.0, 4.0, 4.0};
    Filament actin(count, 1.0, box, rng);
    Myosin myosin(1, 1.5, 0.025, box, rng);
    actin.set_periodic_axes({false, false, false});
    myosin.set_periodic_axes({false, false, false});
    myosin.center[0] = {0.0, 0.0, 0.0};
    myosin.direction[0] = {1.0, 0.0, 0.0};
    constexpr double k_am = 10.0;
    constexpr double k_aa = 10.0;
    constexpr double d_am = 0.055;
    constexpr double d_aa = 0.0506;
    constexpr double radius =
        (k_am * d_am + k_aa * d_aa) / (k_am + k_aa);
    for (int i = 0; i < count; ++i) {
        const double phi = 2.0 * M_PI * i / count;
        actin.center[i] = {-0.30, radius * std::cos(phi),
                           radius * std::sin(phi)};
        actin.direction[i] = {1.0, 0.0, 0.0};
    }
    Eigen::VectorXd equilibrium = transverse_force(
        actin, myosin, box, k_am, k_aa, d_am, d_aa);
    assert(equilibrium.norm() < 1e-8);

    // Numerical Hessian of the actual AM+AA force implementation.  Fixing the
    // central myosin leaves one global ring-rotation zero mode; every other
    // transverse mode must have non-negative curvature.
    constexpr double epsilon = 1e-6;
    Eigen::MatrixXd hessian(2 * count, 2 * count);
    for (int coordinate = 0; coordinate < 2 * count; ++coordinate) {
        const int index = coordinate / 2;
        double* component = coordinate % 2 == 0
            ? &actin.center_y[index] : &actin.center_z[index];
        *component += epsilon;
        const Eigen::VectorXd plus = transverse_force(
            actin, myosin, box, k_am, k_aa, d_am, d_aa);
        *component -= 2.0 * epsilon;
        const Eigen::VectorXd minus = transverse_force(
            actin, myosin, box, k_am, k_aa, d_am, d_aa);
        *component += epsilon;
        hessian.col(coordinate) = -(plus - minus) / (2.0 * epsilon);
    }
    hessian = 0.5 * (hessian + hessian.transpose());
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(hessian);
    assert(solver.info() == Eigen::Success);
    const auto eigenvalues = solver.eigenvalues();
    assert(eigenvalues.minCoeff() > -1e-5);
    int positive = 0;
    for (int i = 0; i < eigenvalues.size(); ++i) {
        if (eigenvalues[i] > 1e-4) ++positive;
    }
    assert(positive >= 11);
    gsl_rng_free(rng);
    return 0;
}
