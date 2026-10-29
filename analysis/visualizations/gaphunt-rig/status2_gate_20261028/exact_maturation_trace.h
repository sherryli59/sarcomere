#ifndef EXACT_MATURATION_TRACE_H
#define EXACT_MATURATION_TRACE_H
// Diagnostic-only logger. Included solely by the isolated replay translation
// unit; the production simulator and gate implementation remain untouched.
#include <cstdio>
#include <cstdlib>

inline void record_exact_maturation(
    const Sarcomere& model, const Sarcomere::AACrosslink& link, bool admitted) {
    static FILE* file = [] {
        const char* path = std::getenv("ACTN2_GATE_AUDIT_PATH");
        if (!path) std::abort();
        FILE* result = std::fopen(path, "w");
        if (!result) std::abort();
        return result;
    }();
    #pragma omp critical(exact_maturation_trace)
    {
        static size_t veto_geometry_records = 0;
        const bool geometry = admitted || veto_geometry_records++ < 20;
        std::fprintf(file,
            "{\"step\":%zu,\"link_id\":%llu,\"admitted\":%s,\"actins\":[%d,%d]",
            model.current_step, static_cast<unsigned long long>(link.id),
            admitted ? "true" : "false", link.actin_i, link.actin_j);
        if (geometry) {
            const auto patch = compute_actn2_patch_geometry(
                model.actin, link.actin_i, link.actin_j,
                link.s_i, link.s_j, model.box);
            std::fprintf(file,
                ",\"first\":[%.17g,%.17g,%.17g],\"second\":[%.17g,%.17g,%.17g],\"myosins\":[",
                patch.first.x, patch.first.y, patch.first.z,
                patch.second.x, patch.second.y, patch.second.z);
            for (int m = 0; m < model.myosin.n; ++m) {
                const vec c = model.myosin.center[m];
                const vec u = model.myosin.direction[m];
                std::fprintf(file, "%s[%.17g,%.17g,%.17g,%.17g,%.17g,%.17g]",
                    m ? "," : "", c.x,c.y,c.z,u.x,u.y,u.z);
            }
            std::fprintf(file, "]");
        }
        std::fprintf(file, "}\n");
        std::fflush(file);
    }
}
#endif
