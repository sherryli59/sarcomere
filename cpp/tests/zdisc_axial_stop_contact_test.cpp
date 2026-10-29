#include "interaction.h"
#include "zdisc_axial_stop_reference.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include <gsl/gsl_rng.h>

// compute_zdisc_axial_stop_contact against the standalone axial-stop law
// (analysis/standalone_actn2_axial_contact.py) and its defining properties.

namespace {
constexpr double kOnset = 0.050, kStiffness = 1000.0, kForceMax = 30.0;
constexpr double kGateRadius = 0.130, kGateCutoff = 0.140, kEpsilon = 0.005;
constexpr double kHalf = 0.750;

struct Both {
    double energy = 0.0;
    vec force_i{0, 0, 0}, force_j{0, 0, 0}, force_m{0, 0, 0};
    vec torque_i{0, 0, 0}, torque_j{0, 0, 0}, torque_m{0, 0, 0};
    bool applied = false;
};

// Place attachment sites exactly at p and q: actin 0 points +x with s=0.1,
// actin 1 points -x with s=0.1 (both 0.1 um from their centers).
void place(Filament& actin, const vec& p, const vec& q) {
    actin.direction[0] = {1, 0, 0};
    actin.direction[1] = {-1, 0, 0};
    actin.center[0] = p - 0.1 * vec{1, 0, 0};
    actin.center[1] = q - 0.1 * vec{-1, 0, 0};
    actin.update_endpoints();
}

Both both_halves(const Filament& actin, const Myosin& myosin,
                 const std::vector<double>& box, int nodes = 16) {
    Both total;
    for (int end_sign : {-1, 1}) {
        const auto c = compute_zdisc_axial_stop_contact(
            actin, myosin, 0, 1, 0.1, 0.1, 0, end_sign, box, kOnset,
            kStiffness, kForceMax, kGateRadius, kGateCutoff, kEpsilon, nodes);
        if (!c.applied) continue;
        total.applied = true;
        total.energy += c.energy;
        total.force_i += c.force_on_actin_i;
        total.force_j += c.force_on_actin_j;
        total.force_m += c.force_on_myosin;
        total.torque_i += c.physical_torque_on_actin_i;
        total.torque_j += c.physical_torque_on_actin_j;
        total.torque_m += c.physical_torque_on_myosin;
    }
    return total;
}

double lateral(const vec& f, const vec& axis) {
    return (f - f.dot(axis) * axis).norm();
}

double phi(double depth) {
    if (depth <= 0.0) return 0.0;
    const double saturation = kForceMax / kStiffness;
    return depth <= saturation ? 0.5 * kStiffness * depth * depth
                               : kForceMax * (depth - 0.5 * saturation);
}
}  // namespace

int main() {
    const std::vector<double> box{8.0, 8.0, 8.0};
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    Filament actin(2, 1.0, box, rng);
    Myosin myosin(1, 2.0 * kHalf, 0.05, box, rng);
    actin.set_periodic_axes({false, false, false});
    myosin.set_periodic_axes({false, false, false});
    auto set_myosin = [&](const vec& center, vec axis) {
        myosin.center[0] = center;
        myosin.direction[0] = axis;
        myosin.update_endpoints();
    };

    // 1. Same energy, endpoint forces and myosin torque as the standalone law.
    double worst_energy = 0, worst_force = 0, worst_refined = 0;
    for (const auto& ref : kAxialStopReference) {
        set_myosin({ref.center[0], ref.center[1], ref.center[2]},
                   {ref.axis[0], ref.axis[1], ref.axis[2]});
        place(actin, {ref.p0[0], ref.p0[1], ref.p0[2]},
              {ref.p1[0], ref.p1[1], ref.p1[2]});
        const vec fi{ref.force_i[0], ref.force_i[1], ref.force_i[2]};
        const vec fj{ref.force_j[0], ref.force_j[1], ref.force_j[2]};
        const vec tm{ref.torque_m[0], ref.torque_m[1], ref.torque_m[2]};
        for (int nodes : {16, 128}) {
            const auto c = both_halves(actin, myosin, box, nodes);
            if (c.applied != (ref.energy > 0.0) ||
                std::abs(c.energy - ref.energy) > 1e-4) {
                std::cerr << "case " << (&ref - kAxialStopReference) << " nodes "
                          << nodes << " energy " << c.energy << " ref "
                          << ref.energy << " Fi " << c.force_i.x << ","
                          << c.force_i.y << "," << c.force_i.z << " ref "
                          << fi.x << "," << fi.y << "," << fi.z << "\n";
            }
            assert(c.applied == (ref.energy > 0.0));
            const double e = std::abs(c.energy - ref.energy);
            const double f = std::max({(c.force_i - fi).norm(),
                                       (c.force_j - fj).norm(),
                                       (c.torque_m - tm).norm()});
            if (nodes == 16) {
                worst_energy = std::max(worst_energy, e);
                worst_force = std::max(worst_force, f);
            } else {
                worst_refined = std::max(worst_refined, std::max(e, f));
            }
        }
    }
    // Production 16-node rule on the split support matches the reference to
    // the reference's own accuracy (4096 nodes across the force kink).
    std::cerr << "reference: 16-node |dE| " << worst_energy << " |dF| "
              << worst_force << "; 128-node " << worst_refined << "\n";
    assert(worst_energy < 1e-8 && worst_force < 1e-6);
    assert(worst_refined < 1e-6);

    // The probe that exposed the weighted law: 40 nm outward connector,
    // 20 nm burial.  Weighted alpha=3 gave 19.9 pN axial and 30.3 pN radial.
    set_myosin({0, 0, 0}, {1, 0, 0});
    place(actin, {.730, .080, 0}, {.730, .120, 0});
    {
        const auto c = both_halves(actin, myosin, box);
        const vec net = c.force_i + c.force_j;
        assert(std::abs(net.x - 29.999296328004895) < 1e-6);
        assert(std::hypot(net.y, net.z) == 0.0);
    }

    // 2. No sideways force on either endpoint or on the myosin, and the
    //    myosin reaction is equal and opposite.  Angular momentum balances.
    std::mt19937_64 random(20261029);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    auto uniform = [&](double lo, double hi) { return lo + (hi - lo) * unit(random); };
    vec axis{1, .03, -.02};
    axis.normalize();
    const vec center{.008, -.003, .004};
    set_myosin(center, axis);
    int contacts = 0;
    for (int trial = 0; trial < 2000; ++trial) {
        const vec p{uniform(.62, .84), uniform(-.16, .16), uniform(-.16, .16)};
        const vec q = p + vec{uniform(-.05, .05), uniform(-.07, .07), uniform(-.07, .07)};
        place(actin, p, q);
        const auto c = both_halves(actin, myosin, box);
        if (!c.applied) continue;
        ++contacts;
        assert(lateral(c.force_i, axis) < 1e-12);
        assert(lateral(c.force_j, axis) < 1e-12);
        assert((c.force_m + c.force_i + c.force_j).norm() < 1e-12);
        assert(std::abs((c.force_i + c.force_j).dot(axis)) <= kForceMax + 1e-9);
        const vec balance = (p - center).cross(c.force_i) +
            (q - center).cross(c.force_j) + c.torque_m;
        assert(balance.norm() < 1e-12);
    }
    assert(contacts > 500);

    // 3. Axial components are exact derivatives of the energy (fixed radial
    //    positions): rigid connector shift, one endpoint, and the myosin.
    const double h = 1e-7;
    double worst_fd = 0.0;
    for (int trial = 0; trial < 300; ++trial) {
        const vec p{uniform(.66, .82), uniform(-.12, .12), uniform(-.12, .12)};
        const vec q = p + vec{uniform(-.03, .03), uniform(-.05, .05), uniform(-.05, .05)};
        place(actin, p, q);
        const auto c = both_halves(actin, myosin, box, 128);
        if (!c.applied) continue;
        auto energy_at = [&](const vec& dp, const vec& dq, const vec& dc) {
            place(actin, p + dp, q + dq);
            set_myosin(center + dc, axis);
            const double e = both_halves(actin, myosin, box, 128).energy;
            place(actin, p, q);
            set_myosin(center, axis);
            return e;
        };
        const vec zero{0, 0, 0}, step = h * axis;
        const double rigid = -(energy_at(step, step, zero) -
                               energy_at(-1.0 * step, -1.0 * step, zero)) / (2 * h);
        const double first = -(energy_at(step, zero, zero) -
                               energy_at(-1.0 * step, zero, zero)) / (2 * h);
        const double body = -(energy_at(zero, zero, step) -
                              energy_at(zero, zero, -1.0 * step)) / (2 * h);
        worst_fd = std::max({worst_fd,
            std::abs(rigid - (c.force_i + c.force_j).dot(axis)),
            std::abs(first - c.force_i.dot(axis)),
            std::abs(body - c.force_m.dot(axis))});
    }
    assert(worst_fd < 1e-3);

    // 4. Gate, onset and saturation edges for a point connector.
    set_myosin({0, 0, 0}, {1, 0, 0});
    auto point_force = [&](double x, double rho) {
        place(actin, {x, rho, 0}, {x, rho, 0});
        const auto c = both_halves(actin, myosin, box);
        return (c.force_i + c.force_j).x;
    };
    const double z_full = .760, z_eps = std::hypot(z_full, kEpsilon);
    const double full = std::min(kStiffness * (kHalf + kOnset - z_eps), kForceMax) *
        z_full / z_eps;
    assert(std::abs(point_force(z_full, 0.0) - full) < 1e-9);
    assert(std::abs(point_force(z_full, kGateRadius) - full) < 1e-9);
    assert(std::abs(point_force(z_full, 0.135) - 0.5 * full) < 1e-9);
    assert(point_force(z_full, 0.13999) > 0.0);
    assert(point_force(z_full, 0.14001) == 0.0);
    assert(point_force(.7996, .05) > 0.0);
    assert(point_force(.8001, .05) == 0.0);
    assert(std::abs(point_force(.60, .05) - kForceMax * .60 / std::hypot(.60, kEpsilon)) < 1e-9);
    // The opposite tip pushes the other way, and each half sees only its own.
    assert(std::abs(point_force(-z_full, 0.05) + full) < 1e-9);
    place(actin, {-z_full, .05, 0}, {-z_full, .05, 0});
    assert(!compute_zdisc_axial_stop_contact(
        actin, myosin, 0, 1, .1, .1, 0, 1, box, kOnset, kStiffness,
        kForceMax, kGateRadius, kGateCutoff, kEpsilon).applied);

    // 5. Documented non-conservation: a point loop inside the gate does no
    //    work; one crossing the gate edge does Phi*(g_in - g_out).
    auto loop_work = [&](double za, double zb, double ra, double rb) {
        const vec corners[5]{{za, ra, 0}, {zb, ra, 0}, {zb, rb, 0},
                             {za, rb, 0}, {za, ra, 0}};
        double work = 0.0;
        const int steps = 4096;
        for (int side = 0; side < 4; ++side) {
            const vec delta = (1.0 / steps) * (corners[side + 1] - corners[side]);
            for (int k = 0; k < steps; ++k) {
                const vec mid = corners[side] + (k + 0.5) * delta;
                place(actin, mid, mid);
                const auto c = both_halves(actin, myosin, box);
                work += (c.force_i + c.force_j).dot(delta);
            }
        }
        return work;
    };
    const auto depth = [](double z) { return kHalf + kOnset - std::hypot(z, kEpsilon); };
    assert(std::abs(loop_work(.770, .790, .080, .120)) < 1e-9);
    const double analytic = phi(depth(.770)) - phi(depth(.790));
    assert(std::abs(loop_work(.770, .790, .100, .150) - analytic) < 1e-7);

    // 6. Periodic images: translation by a box vector changes nothing, and a
    //    tilted long myosin reaches a non-nearest image.
    actin.set_periodic_axes({true, true, true});
    myosin.set_periodic_axes({true, true, true});
    set_myosin({0, 3.9, 0}, {1, 0, 0});
    place(actin, {.74, -3.98, 0}, {.74, -3.94, 0});
    const auto wrapped = both_halves(actin, myosin, box);
    assert(wrapped.applied);
    place(actin, {8.74, 4.02, 0}, {8.74, 4.06, 0});
    const auto shifted = both_halves(actin, myosin, box);
    assert(std::abs(shifted.energy - wrapped.energy) < 1e-12);
    assert((shifted.force_m - wrapped.force_m).norm() < 1e-9);
    {
        const double s = std::sqrt(.5);
        Myosin long_tilted(1, 1.5, .05, box, rng);
        long_tilted.set_periodic_axes({true, true, true});
        long_tilted.center[0] = {0, 0, 0};
        long_tilted.direction[0] = {s, s, 0};
        long_tilted.update_endpoints();
        // Nearest image: rho ~0.68 um, no contact.  The image shifted by
        // -1 in x sits 0.03 um from the axis on the -tip half.
        place(actin, {.48, -.48, 0}, {.49, -.48, 0});
        const auto image = compute_zdisc_axial_stop_contact(
            actin, long_tilted, 0, 1, .1, .1, 0, -1, {1, 2, 2}, kOnset,
            kStiffness, kForceMax, kGateRadius, kGateCutoff, kEpsilon);
        assert(image.applied);
    }

    // 7. Invalid parameters apply nothing.
    set_myosin({0, 0, 0}, {1, 0, 0});
    actin.set_periodic_axes({false, false, false});
    myosin.set_periodic_axes({false, false, false});
    place(actin, {.73, .08, 0}, {.73, .12, 0});
    assert(!compute_zdisc_axial_stop_contact(
        actin, myosin, 0, 1, .1, .1, 0, 1, box, kOnset, kStiffness,
        kForceMax, kGateRadius, kGateRadius, kEpsilon).applied);
    assert(!compute_zdisc_axial_stop_contact(
        actin, myosin, 0, 1, .1, .1, 0, 1, box, kOnset, kStiffness,
        kForceMax, kGateRadius, kGateCutoff, 0.0).applied);

    gsl_rng_free(rng);
    std::cout << "zdisc axial-stop contact: standalone reference (16-node max |dE| "
              << worst_energy << " pN um, |dF| " << worst_force
              << " pN; 128-node " << worst_refined << "), " << contacts
              << " axial-only contacts, FD " << worst_fd << " pN, loop work passed\n";
    return 0;
}
