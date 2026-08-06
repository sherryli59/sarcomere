#include "aa_refractory.h"

#include <cassert>

int main() {
    using aa_refractory::can_attempt_association;
    assert(can_attempt_association(true, 10, 100));
    assert(!can_attempt_association(false, 10, 11));
    assert(can_attempt_association(false, 11, 11));
    assert(can_attempt_association(false, 12, 11));
    assert(can_attempt_association(false, 0, 0));
    assert(aa_refractory::mechanically_active(true, false));
    assert(aa_refractory::mechanically_active(false, true));
    assert(!aa_refractory::mechanically_active(false, false));
    return 0;
}
