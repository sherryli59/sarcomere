#include "neighborlist.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <vector>

namespace {

void require_cross_boundary_pair(const std::array<bool, 3>& periodic) {
    const std::vector<double> box{4.0, 2.0, 2.0};
    // cutoff=1.8 creates exactly two x cells. Actin is global particle 0 in
    // the upper x cell; myosin is global particle 1 in the lower x cell. The
    // historical +offset-only rule missed this ordering in nonperiodic boxes.
    NeighborList neighbors(1.8, box, 0.05, periodic);
    neighbors.initialize(
        {0.10}, {0.0}, {0.0},
        {-0.10}, {0.0}, {0.0});

    const auto split = neighbors.get_neighbors_by_type(0);
    assert(split.first.empty());
    assert(split.second.size() == 1);
    assert(split.second.front() == 0);

    const auto reverse = neighbors.get_neighbors_by_type(1);
    assert(reverse.first.size() == 1);
    assert(reverse.first.front() == 0);
    assert(reverse.second.empty());
}

}  // namespace

int main() {
    require_cross_boundary_pair({false, false, false});
    // The duplicate-cell guard must also prevent duplicate pairs when -1 and
    // +1 wrap to the same cell under two-cell periodicity.
    require_cross_boundary_pair({true, true, true});
    std::cout << "two-cell cross-boundary neighbor test passed\n";
    return 0;
}
