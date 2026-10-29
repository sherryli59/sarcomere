#include "actn2_maturation_geometry.h"
#include "reference_geometry.h"
#include <cassert>
#include <chrono>
#include <fstream>
#include <random>

using utils::vec;
const std::array<bool, 3> pbc{true, true, true};

int main(int argc, char** argv) {
    assert(argc == 2);
    const std::string input = argv[1];
    std::ifstream fl(input + "/links.txt"), fm(input + "/myosins.txt"),
                  fb(input + "/box.txt");
    int nl = 0, nm = 0;
    fl >> nl; fm >> nm;
    assert(fl && fm && nl > 0 && nm > 0);
    std::vector<vec> p(nl), q(nl), c(nm), u(nm);
    for (int k = 0; k < nl; ++k)
        fl >> p[k].x >> p[k].y >> p[k].z >> q[k].x >> q[k].y >> q[k].z;
    for (int k = 0; k < nm; ++k)
        fm >> c[k].x >> c[k].y >> c[k].z >> u[k].x >> u[k].y >> u[k].z;
    std::vector<double> box(3);
    fb >> box[0] >> box[1] >> box[2];
    assert(fl && fm && fb && box[0] > 0 && box[1] > 0 && box[2] > 0);
    constexpr double radius = actn2_maturation::deletion_radius;
    constexpr double buffer = actn2_maturation::deletion_buffer;
    auto baseline = [&](int k) {
        for (int m = 0; m < nm; ++m)
            if (actn2_reference::connector_overlaps_myosin(
                p[k], q[k], c[m], u[m], 1.5, box, pbc, radius, buffer)) return true;
        return false;
    };
    auto optimized = [&](int k) {
        vec d = q[k] - p[k]; d.pbc_wrap(box, pbc);
        const vec midpoint = p[k] + .5 * d;
        const double reach = std::hypot(.75-buffer, radius) + .5*d.norm();
        for (int m = 0; m < nm; ++m) {
            vec rel = midpoint - c[m]; rel.pbc_wrap(box, pbc);
            if (rel.norm_squared() > reach*reach) continue;
            if (actn2_maturation::connector_overlaps_myosin(
                p[k], q[k], c[m], u[m], 1.5, box, pbc, radius, buffer)) return true;
        }
        return false;
    };
    int hits = 0;
    for (int k = 0; k < nl; ++k) {
        assert(baseline(k) == optimized(k));
        hits += baseline(k);
    }
    std::mt19937_64 random(20261028);
    std::uniform_real_distribution<double> value(-3, 3), positive(.01, 3);
    for (int k = 0; k < 100000; ++k) {
        vec a{value(random),value(random),value(random)};
        vec b{value(random),value(random),value(random)};
        vec center{value(random),value(random),value(random)};
        vec axis{value(random),value(random),value(random)};
        double length=positive(random), r=positive(random)*.1;
        std::vector<double> domain{positive(random),positive(random),positive(random)};
        std::array<bool,3> periodic{bool(k&1),bool(k&2),bool(k&4)};
        assert(actn2_reference::connector_overlaps_myosin(a,b,center,axis,length,domain,periodic,r) ==
               actn2_maturation::connector_overlaps_myosin(a,b,center,axis,length,domain,periodic,r));
    }
    auto timing = [&](auto function) {
        volatile int count = 0;
        const auto start = std::chrono::steady_clock::now();
        for (int rep = 0; rep < 20; ++rep)
            for (int k = 0; k < nl; ++k) count += function(k);
        assert(count == 20*hits);
        return std::chrono::duration<double,std::milli>(
            std::chrono::steady_clock::now()-start).count()/20;
    };
    const double old_ms = timing(baseline), new_ms = timing(optimized);
    std::cout << "{\"links\":" << nl << ",\"myosins\":" << nm
              << ",\"deleted\":" << hits << ",\"random_equivalence_cases\":100000"
              << ",\"baseline_ms\":" << old_ms << ",\"optimized_ms\":" << new_ms
              << ",\"speedup\":" << old_ms/new_ms << "}\n";
}
