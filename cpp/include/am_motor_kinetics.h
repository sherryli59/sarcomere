#ifndef AM_MOTOR_KINETICS_H
#define AM_MOTOR_KINETICS_H

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace am_motor_kinetics {

enum class State : std::uint8_t {
    detached = 0,
    weak = 1,
    strong = 2,
    off_srx = 3,
};

inline double transition_probability(double rate, double dt) {
    if (!(rate > 0.0) || !(dt > 0.0)) return 0.0;
    const double exponent = std::min(rate * dt, 700.0);
    return -std::expm1(-exponent);
}

// Select the requested rounded number of initially ON heads while spreading
// them deterministically across the ordered half-filament pool.  Rounding is
// important for small coarse-grained pools: nine heads at a 0.10 equilibrium
// fraction should contain one ON head, rather than being truncated to zero.
inline bool initially_on_head(int head_index, int head_count,
                              double on_fraction) {
    if (head_index < 0 || head_count <= 0 || head_index >= head_count) {
        return false;
    }
    const int target_on = std::clamp(
        static_cast<int>(std::lround(static_cast<double>(head_count) *
                                    std::clamp(on_fraction, 0.0, 1.0))),
        0, head_count);
    return (head_index + 1) * target_on / head_count >
           head_index * target_on / head_count;
}

// A normalized calcium-like transient. A non-positive period selects tonic
// activation, which is useful for matching the legacy time-independent model.
inline double activation(double time, double period, double rise, double decay,
                         double floor, double active_fraction = 0.5) {
    floor = std::clamp(floor, 0.0, 1.0);
    if (!(period > 0.0)) return 1.0;
    rise = std::max(rise, 1e-12);
    decay = std::max(decay, rise + 1e-12);
    double phase = std::fmod(std::max(time, 0.0), period);
    active_fraction = std::clamp(active_fraction, 0.0, 1.0);
    if (!(active_fraction > 0.0) ||
        phase >= active_fraction * period) {
        return floor;
    }
    const double transient = std::exp(-phase / decay) - std::exp(-phase / rise);
    const double peak_time = rise * decay * std::log(decay / rise) / (decay - rise);
    const double peak = std::exp(-peak_time / decay) - std::exp(-peak_time / rise);
    const double normalized = peak > 0.0 ? std::clamp(transient / peak, 0.0, 1.0) : 0.0;
    return floor + (1.0 - floor) * normalized;
}

inline double hill_gate(double signal, double half_activation, double exponent) {
    signal = std::max(0.0, signal);
    half_activation = std::max(half_activation, 1e-12);
    exponent = std::max(exponent, 1e-12);
    const double signal_power = std::pow(signal, exponent);
    const double half_power = std::pow(half_activation, exponent);
    const double denominator = signal_power + half_power;
    return denominator > 0.0 ? signal_power / denominator : 0.0;
}

// Weak, non-force-producing encounters remain possible on a nominally closed
// thin filament, whereas recruitment into the strong working state requires
// thin-filament activation.  The basal fraction represents low-affinity
// collision/weak-binding access; it must not create active force by itself.
inline double weak_attachment_gate(double thin_activation,
                                   double basal_fraction) {
    const double q = std::clamp(thin_activation, 0.0, 1.0);
    const double basal = std::clamp(basal_fraction, 0.0, 1.0);
    return basal + (1.0 - basal) * q;
}

// Exact update for dq/dt = k_open(1-q) - k_close*q, with calcium and
// currently strong heads contributing independently to k_open. This form is
// bounded and remains stable when rates are fast relative to the integrator.
inline double thin_activation_step(
    double q, double calcium, double strong_fraction, double dt,
    double calcium_open_rate, double calcium_half, double calcium_hill,
    double strong_feedback_rate, double strong_half, double strong_hill,
    double close_rate) {
    q = std::clamp(q, 0.0, 1.0);
    const double calcium_gate =
        hill_gate(calcium, calcium_half, calcium_hill);
    // Strong binding potentiates a calcium-permissive thin filament. It does
    // not create a calcium-independent open state, which would latch the model
    // through diastole when residual strong heads are slow to detach.
    const double open_rate =
        std::max(0.0, calcium_open_rate) * calcium_gate +
        std::max(0.0, strong_feedback_rate) * calcium_gate *
            hill_gate(strong_fraction, strong_half, strong_hill);
    close_rate = std::max(0.0, close_rate);
    const double total_rate = open_rate + close_rate;
    if (!(dt > 0.0) || !(total_rate > 0.0)) return q;
    const double equilibrium = open_rate / total_rate;
    const double updated = equilibrium +
        (q - equilibrium) * std::exp(-std::min(total_rate * dt, 700.0));
    return std::clamp(updated, 0.0, 1.0);
}

// Calcium-only stationary state used to prepare a tonically activated thin
// filament without pre-running any binding chemistry at fixed geometry.
// Strong-head feedback is deliberately absent because all heads are initially
// detached in this preparation.
inline double thin_activation_calcium_equilibrium(
    double calcium, double calcium_open_rate, double calcium_half,
    double calcium_hill, double close_rate) {
    const double open_rate = std::max(0.0, calcium_open_rate) *
        hill_gate(calcium, calcium_half, calcium_hill);
    close_rate = std::max(0.0, close_rate);
    const double total_rate = open_rate + close_rate;
    return total_rate > 0.0 ? std::clamp(open_rate / total_rate, 0.0, 1.0)
                            : 0.0;
}

// Resisting force slows both the power stroke and strong-state detachment.
// force_scale is an explicit model-unit force scale, avoiding an implicit
// conversion between the simulation's coarse-grained force and pN.
inline double load_slowed_rate(double zero_load_rate, double resisting_force,
                               double force_scale) {
    if (!(zero_load_rate > 0.0)) return 0.0;
    if (!(force_scale > 0.0)) return zero_load_rate;
    const double exponent = -std::max(0.0, resisting_force) / force_scale;
    return zero_load_rate * std::exp(std::max(exponent, -700.0));
}

// A parallel load-insensitive chemical pathway prevents a resisting load from
// freezing an attached motor indefinitely. `floor_fraction` is the fraction
// of the zero-load rate retained at asymptotically large load.
inline double load_slowed_rate_with_floor(double zero_load_rate,
                                          double resisting_force,
                                          double force_scale,
                                          double floor_fraction) {
    if (!(zero_load_rate > 0.0)) return 0.0;
    const double floor = std::clamp(floor_fraction, 0.0, 1.0);
    if (!(force_scale > 0.0)) return zero_load_rate;
    const double exponent = -std::max(0.0, resisting_force) / force_scale;
    const double catch_branch = std::exp(std::max(exponent, -700.0));
    return zero_load_rate * (floor + (1.0 - floor) * catch_branch);
}

// Elastic distortion lowers the activation barrier for leaving an attached
// state.  `beta_energy` is dimensionless (E/kBT) and `barrier_fraction`
// specifies how much of that stored energy couples to the transition state.
// The exponent cap prevents a badly distorted diagnostic configuration from
// overflowing before the KMC event can remove the head.
inline double strain_accelerated_rate(double base_rate, double beta_energy,
                                      double barrier_fraction,
                                      double maximum_exponent = 20.0) {
    if (!(base_rate > 0.0)) return 0.0;
    if (!(beta_energy > 0.0) || !(barrier_fraction > 0.0)) {
        return base_rate;
    }
    const double exponent = std::clamp(
        barrier_fraction * beta_energy, 0.0,
        std::max(0.0, maximum_exponent));
    return base_rate * std::exp(exponent);
}

// Axial and transverse distortion need not couple equally to the chemical
// escape coordinate.  In particular, a coarse-grained transverse head spring
// can provide mechanical reach without automatically becoming a strong slip-
// bond penalty.  Energies are in the simulation energy unit (pN um), beta has
// the reciprocal unit, and both barrier fractions are dimensionless.
inline double anisotropic_strain_accelerated_rate(
    double base_rate, double beta,
    double axial_energy, double axial_barrier_fraction,
    double transverse_energy, double transverse_barrier_fraction,
    double maximum_exponent = 20.0) {
    if (!(base_rate > 0.0)) return 0.0;
    const double exponent = std::clamp(
        std::max(0.0, beta) * (
            std::max(0.0, axial_barrier_fraction) *
                std::max(0.0, axial_energy) +
            std::max(0.0, transverse_barrier_fraction) *
                std::max(0.0, transverse_energy)),
        0.0, std::max(0.0, maximum_exponent));
    return base_rate * std::exp(exponent);
}

// Convert a saturating pair-level transverse stiffness into an equal
// per-head contribution. The total is continuous in occupancy and bounded by
// `pair_stiffness`; in particular, the final head carries only 1/saturation
// of the mature pair constraint.
inline double occupancy_scaled_per_head_stiffness(double pair_stiffness,
                                                  int occupied_heads,
                                                  int saturation_heads) {
    if (!(pair_stiffness > 0.0) || occupied_heads <= 0) return 0.0;
    const int saturation = std::max(1, saturation_heads);
    const double total = pair_stiffness * std::min(
        1.0, static_cast<double>(occupied_heads) /
                 static_cast<double>(saturation));
    return total / static_cast<double>(occupied_heads);
}

// Number of simultaneous AM attachments available to one thin-filament
// neighbor of a myosin half.  The conserved half-filament pool is shared over
// the approximately six thin filaments in the cardiac lattice, and only the
// axially overlapping fraction of those slots is available.  A positive but
// sub-slot overlap is rounded up to one physical head; zero overlap has zero
// capacity.  A configured full-overlap cap <= 0 disables neighbor sharing.
inline int overlap_limited_neighbor_capacity(int heads_per_half,
                                             int full_overlap_cap,
                                             double productive_overlap,
                                             double head_bearing_length) {
    const int heads = std::max(0, heads_per_half);
    if (heads == 0 || !(productive_overlap > 0.0) ||
        !(head_bearing_length > 0.0)) return 0;
    const int full_capacity = full_overlap_cap > 0
        ? std::min(heads, full_overlap_cap) : heads;
    const double fraction = std::clamp(
        productive_overlap / head_bearing_length, 0.0, 1.0);
    return std::min(full_capacity, std::max(
        1, static_cast<int>(std::ceil(full_capacity * fraction - 1e-12))));
}

inline int weak_slots_after_strong(int neighbor_capacity,
                                   int strong_occupancy) {
    return std::max(0, neighbor_capacity - std::max(0, strong_occupancy));
}

inline double thick_recruitment_rate(double zero_tension_rate,
                                     double half_tension,
                                     double tension_scale) {
    if (!(zero_tension_rate > 0.0)) return 0.0;
    if (!(tension_scale > 0.0)) return zero_tension_rate;
    const double exponent = std::min(
        std::max(0.0, half_tension) / tension_scale, 20.0);
    return zero_tension_rate * std::exp(exponent);
}

inline double thick_folding_rate(double zero_tension_rate,
                                 double half_tension,
                                 double tension_scale) {
    return load_slowed_rate(
        zero_tension_rate, std::max(0.0, half_tension), tension_scale);
}

// A strongly bound cross-bridge is a tension-only active element. A power
// stroke advances its rest coordinate in the productive direction; subsequent
// overshoot unloads the head instead of making the active motor pull backward.
inline double active_crossbridge_force(double stiffness,
                                       double productive_strain,
                                       double stall_force,
                                       bool bidirectional = false) {
    if (!(stiffness > 0.0) || !std::isfinite(productive_strain)) return 0.0;
    if (!bidirectional) {
        if (!(productive_strain > 0.0)) return 0.0;
        const double force = stiffness * productive_strain;
        return stall_force > 0.0 ? std::min(force, stall_force) : force;
    }
    // The attached head is an elastic element, not a ratchet.  The power
    // stroke advances its rest coordinate one-way (chemistry), but once the
    // thin filament overshoots that rest point the bound head is strained the
    // other way and pulls BACK.  Rectifying this to zero leaves nothing
    // opposing over-penetration past the M-line, which is the restoring force
    // real cross-bridges provide during lengthening.
    const double force = stiffness * productive_strain;
    if (stall_force > 0.0) {
        return std::clamp(force, -stall_force, stall_force);
    }
    return force;
}

// A weakly attached head is a passive elastic bridge.  Unlike the active
// working stroke it resists displacement in either axial direction; rapid KMC
// turnover makes this a transient friction rather than a permanent clamp.
inline double passive_crossbridge_force(double stiffness, double strain) {
    if (!(stiffness > 0.0) || !std::isfinite(strain)) return 0.0;
    return stiffness * strain;
}

// Linear Hill-like force-velocity branch for the finite working stroke. The
// internal motor coordinate advances at unloaded_velocity at zero load and
// stops at stall. Motion is bounded separately by the remaining stroke.
// Free-energy change of advancing the head's rest coordinate by the full
// working stroke `d` while the thin-filament attachment sits at axial
// coordinate `z`:  dG = 1/2 k [(d-z)^2 - z^2] = k d (d/2 - z).
inline double stroke_mechanical_free_energy(double stiffness, double stroke,
                                            double z) {
    return stiffness * stroke * (0.5 * stroke - z);
}

// Forward / reverse rates for the discrete working stroke, split by a Bell
// coefficient and satisfying k_+/k_- = exp(-beta (dG_mech - dG_chem)).
struct StrokeRates { double forward; double reverse; };
inline StrokeRates stroke_transition_rates(double base_rate, double beta,
                                           double dG_mech, double dG_chem,
                                           double theta, double max_exponent) {
    StrokeRates r{0.0, 0.0};
    if (!(base_rate > 0.0)) return r;
    const double dG = dG_mech - dG_chem;
    const double th = std::clamp(theta, 0.0, 1.0);
    const double fwd = std::clamp(-beta * th * dG, -max_exponent, max_exponent);
    const double rev = std::clamp(beta * (1.0 - th) * dG, -max_exponent, max_exponent);
    r.forward = base_rate * std::exp(fwd);
    r.reverse = base_rate * std::exp(rev);
    return r;
}

inline double working_stroke_velocity(double unloaded_velocity,
                                      double resisting_force,
                                      double stall_force) {
    if (!(unloaded_velocity > 0.0)) return 0.0;
    if (!(stall_force > 0.0)) return unloaded_velocity;
    const double load_fraction = std::clamp(
        std::max(0.0, resisting_force) / stall_force, 0.0, 1.0);
    return unloaded_velocity * (1.0 - load_fraction);
}

}  // namespace am_motor_kinetics

#endif
