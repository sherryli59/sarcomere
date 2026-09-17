#ifndef AA_REFRACTORY_H
#define AA_REFRACTORY_H

#include <cstddef>

namespace aa_refractory {

// A persistent bond is always evaluated for detachment.  An unbound pair may
// attempt association only after its post-detachment refractory interval.
inline bool can_attempt_association(
    bool bonded_previous_step,
    std::size_t current_step,
    std::size_t refractory_until) {
    return bonded_previous_step || current_step >= refractory_until;
}

// Capture distance gates formation only. A previously occupied bond remains a
// mechanical spring until its explicit KMC detachment event.
inline bool mechanically_active(bool bonded_previous_step, bool within_capture) {
    return bonded_previous_step || within_capture;
}

}  // namespace aa_refractory

#endif
