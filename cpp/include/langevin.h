#ifndef LANGEVIN_H
#define LANGEVIN_H

#include "sarcomere.h"
#include "utils.h"
#include "h5_utils.h"
#include "rod_mobility.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <omp.h>
#include <mutex>
#include <cstdint>

// Declare the global mutex for saving (defined in the source file)
extern std::mutex save_mutex;

enum class Actn2TipProjectionTraversal {
    spatial_incidence,
    sparse_incidence,
    full_scan,
};

struct Actn2TipProjectionProfile {
    std::uint64_t calls = 0;
    std::uint64_t active_links = 0;
    std::uint64_t unique_actin_pairs = 0;
    std::uint64_t constraints = 0;
    std::uint64_t broadphase_tests = 0;
    std::uint64_t broadphase_candidates = 0;
    std::uint64_t solve_evaluations = 0;
    std::uint64_t audit_evaluations = 0;
    std::uint64_t corrections = 0;
    std::uint64_t sweeps = 0;
    std::uint64_t exhausted_calls = 0;
    std::uint64_t audited_calls = 0;
    std::uint64_t certified_without_audit_calls = 0;
    double build_seconds = 0.0;
    double solve_seconds = 0.0;
    double audit_seconds = 0.0;
    double total_seconds = 0.0;
};

struct LangevinRuntimeProfile {
    std::uint64_t steps = 0;
    double update_system_seconds = 0.0;
    double noise_seconds = 0.0;
    double myosin_motion_seconds = 0.0;
    double actin_motion_seconds = 0.0;
    double myosin_contact_seconds = 0.0;
    double actn2_contact_seconds = 0.0;
    double sample_total_seconds = 0.0;
    double save_state_seconds = 0.0;
    double save_resume_seconds = 0.0;
    double loop_total_seconds = 0.0;
};

class Langevin {
public:
    // Constructor.
    // Note: Parameters are passed by reference as in the original code.
    Langevin(Sarcomere& model0, double& beta0, double& dt0, double& D0_actin_trans,
        double& D0_actin_rot, double& D0_myosin_trans, double& D0_myosin_rot,
        int& save_every0, bool& resume, bool is3D,
        double max_actin_displacement, double max_myosin_displacement,
        double max_actin_rotation, double max_myosin_rotation,
        int resume_frame = -1);
    // Destructor.
    ~Langevin();

    // Run the Langevin simulation for a given number of steps.
    void run_langevin(int nsteps, gsl_rng* rng, int& fix_myosin);
    void volume_exclusion(int nsteps, gsl_rng* rng, int& fix_myosin,
                          bool save_frames = true);
    // Take a single simulation step.
    void sample_step(double& dt, gsl_rng* rng, int& fix_myosin);
    void resolve_myosin_hard_contacts(int fix_myosin, int max_sweeps = 4);
    void resolve_actn2_tip_hard_contacts(int fix_myosin,
                                         int max_sweeps = 64,
                                         double tolerance = 1e-9,
                                         Actn2TipProjectionTraversal traversal =
                                             Actn2TipProjectionTraversal::spatial_incidence);
    const Actn2TipProjectionProfile& actn2_tip_projection_profile() const {
        return actn2_tip_profile;
    }
    void reset_actn2_tip_projection_profile() {
        actn2_tip_profile = {};
    }
    const LangevinRuntimeProfile& runtime_profile() const {
        return run_profile;
    }
    void reset_runtime_profile() {
        run_profile = {};
    }
    void set_rod_mobility(const rod_mobility::Diffusion& actin_mobility,
                          const rod_mobility::Diffusion& myosin_mobility);
    void set_isometric_actin_axial(bool enabled) {
        constrain_actin_axial = enabled;
    }
    void set_thermal_noise_scale(double scale) {
        thermal_noise_scale = scale;
    }

    // Data members.
    double dt, beta, D_actin_trans, D_actin_rot, D_myosin_trans, D_myosin_rot;
    double max_actin_displacement, max_myosin_displacement;
    double max_actin_rotation, max_myosin_rotation;
    rod_mobility::Diffusion actin_rod_mobility, myosin_rod_mobility;
    bool physical_rod_mobility = false;
    bool constrain_actin_axial = false;
    double thermal_noise_scale = 1.0;
    bool is3D;
    bool skip_initial_save;
    int save_every, start_step, loaded_frame_index;
    // When true, interpret `nsteps` argument to run_* as the total target
    // number of steps (i.e., inclusive of previously-run steps) while
    // resuming. Set in constructor when resume mode is active.
    bool interpret_nsteps_as_total_when_resuming = false;
    Actn2TipProjectionProfile actn2_tip_profile;
    // Per-link survivor buckets for the parallel discovery pass,
    // kept across calls so the sweep costs no allocations.
    std::vector<std::vector<int>> tip_discovery_scratch;
    // Reused per-pair myosin candidate lists, and the stamp array used to
    // build them without a sort.  The token only ever increments, so the
    // stamps never need clearing.
    std::vector<std::vector<int>> tip_pair_candidates;
    std::vector<std::uint64_t> tip_pair_union_stamp;
    std::uint64_t tip_pair_union_token = 0;
    LangevinRuntimeProfile run_profile;
    Sarcomere& model;
};

#endif // LANGEVIN_H
