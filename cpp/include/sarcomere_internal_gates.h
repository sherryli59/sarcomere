#pragma once

// Shared between sarcomere.cpp and sarcomere_actn2.cpp, split out on
// 2026-09-17.  Moved verbatim out of the anonymous namespace that used to
// sit at sarcomere.cpp:1727.  Anonymous namespace in a header gives each
// translation unit its own copy, which is identical by construction here:
// the value is a pure function of two environment variables.

#include <cstdlib>
#include <string>

namespace {
// Escape hatch for A/B timing only.  Every gate it controls is exact -- each
// rejects solely cases where the kernel it guards provably contributes zero --
// so switching it off changes runtime and nothing else.  It exists because
// node-to-node speed varies ~10%, which makes a cumulative speedup stitched
// across jobs meaningless; with this, one job measures the whole set.
const bool opt_gates_enabled = [] {
    for (const char* name : {"SARC_OPT_GATES", "SARC_AM_PREFILTER"}) {
        const char* v = std::getenv(name);
        if (v && std::string(v) == "0") return false;
    }
    return true;
}();
}  // namespace
