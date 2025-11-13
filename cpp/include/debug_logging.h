#ifndef DEBUG_LOGGING_H
#define DEBUG_LOGGING_H

#include "utils.h"
#include <cstdio>

namespace debug_logging {

inline void log_myosin_force(int myosin_index, const vec& contribution, const char* source) {
    constexpr double kThreshold = 100.0;
    double magnitude = contribution.norm();
    if (magnitude > kThreshold) {
        printf("[MyosinForce] idx=%d source=%s delta=(% .6e,% .6e,% .6e) |delta|=% .6e\n",
               myosin_index,
               source ? source : "unknown",
               contribution.x,
               contribution.y,
               contribution.z,
               magnitude);
    }
}

}  // namespace debug_logging

#endif  // DEBUG_LOGGING_H
