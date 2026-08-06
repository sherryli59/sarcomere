#include "aa_neighbor_identity.h"

#include <cassert>

int main() {
    using aa_neighbor_identity::Registry;
    Registry slots;

    assert(slots.reserve(0, 1, 0, 2));
    assert(slots.reserve(0, 2, 0, 2));
    assert(!slots.reserve(0, 3, 0, 2));
    assert(slots.degree(0) == 2);

    // A live status-1 edge survives molecular turnover for the short memory.
    slots.observe(0, 1, 1, 10, 100);
    slots.expire_unseen(11, 5);
    assert(slots.contains(0, 1));
    slots.expire_unseen(16, 5);
    assert(!slots.contains(0, 1));

    // Catch-supported topology receives the longer status-2 protection.
    slots.observe(0, 2, 2, 20, 100);
    slots.expire_unseen(26, 5);
    assert(slots.contains(0, 2));
    slots.expire_unseen(121, 5);
    assert(!slots.contains(0, 2));

    // A contour rupture releases the molecular neighbor identity immediately.
    assert(slots.reserve(1, 2, 200, 2));
    assert(slots.degree(1) == 1);
    assert(slots.degree(2) == 1);
    slots.erase(1, 2);
    assert(!slots.contains(1, 2));
    assert(slots.degree(1) == 0);
    assert(slots.degree(2) == 0);

    // A pinned intended lattice edge is not annealed away.
    assert(slots.reserve(4, 5, 0, 2, true));
    slots.expire_unseen(100000, 0);
    assert(slots.contains(4, 5));
    return 0;
}
