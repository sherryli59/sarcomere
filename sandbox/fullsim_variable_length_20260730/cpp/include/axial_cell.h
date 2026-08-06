#ifndef AXIAL_CELL_H
#define AXIAL_CELL_H

#include "sarcomere.h"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace axial_cell {

enum class Mode {
    fixed,
    stress_clamp,
    series_spring,
};

// One coarse-grained titin molecule. Its thick-filament binding domain is
// permanently assigned to a material site on one myosin half in this minimal
// assembly model. Its Z-terminal domain binds one individual actin barbed end;
// there is deliberately no Z-node, ACTN2-component, or centroid anchor.
struct ExplicitTitin {
    int myosin = -1;
    int half = -1;
    int track = -1;
    double anchor_azimuth = 0.0;
    int bound_actin = -1;
    int candidate_actin = -1;
    // Previous molecular partner retained only for finite geminate rebinding.
    // This is not a bond: no force is applied while the titin is unbound.
    int last_bound_actin = -1;
    std::size_t attach_step = 0;
    std::size_t last_detach_step = 0;
    double distance = 0.0;
    double extension = 0.0;
    double force = 0.0;
    double candidate_distance = 0.0;
    double candidate_energy = 0.0;
    double candidate_on_rate = 0.0;
    double attachment_distance = 0.0;
    double attachment_energy = 0.0;
    double attachment_force = 0.0;
    bool bound = false;
};

struct Config {
    Mode mode = Mode::fixed;
    // Fixed length normally bypasses the entire axial-cell module. Diagnostic
    // positive controls may keep Lx frozen while still applying configured
    // reciprocal titin/register forces.
    bool apply_internal_forces_when_fixed = false;
    double mobility = 0.0;               // um / (pN s)
    double applied_tension = 0.0;        // pN; positive expands Lx
    double series_k = 0.0;               // pN / um
    double series_rest_length = 0.0;     // um
    double min_length = 3.0;             // um
    double max_length = 6.0;             // um

    // Diagnostic always-connected barbed-end/Z--M springs.  These are kept
    // separate from the production bond-state-dependent titin implementation.
    double titin_k = 0.0;                // pN / um
    double titin_rest_length = 0.0;      // um

    // Explicit local titin tethers. Capture depends only on the distance and
    // polarity of an individual actin barbed-end binding site. It is independent
    // of myosin-head occupancy and ACTN2 connectivity. Bound tethers exert a
    // reciprocal 3D endpoint force and are slack in compression.
    bool explicit_titin = false;
    int explicit_titin_per_half = 1;
    // Maximum simultaneously bound titins per coarse actin. Keep this at one
    // for resolved hexagonal actins; a minimal rod may represent a bundle and
    // accept several independent molecular tracks.
    int explicit_titin_max_per_actin = 1;
    double explicit_titin_k = 0.0;              // pN / um
    double explicit_titin_slack_length = 0.0;   // um
    double explicit_titin_capture_distance = 0.0;  // um
    double explicit_titin_k_on = 0.0;           // 1 / s
    double explicit_titin_k_off = 0.0;          // 1 / s
    double explicit_titin_force_scale = 0.0;    // pN; <=0 is force independent
    // When enabled, elastic strain enters both transition rates so that
    // k_on(d) / k_off(d) carries the Boltzmann factor exp[-beta U(d)].
    bool explicit_titin_strain_gated_kinetics = false;
    double explicit_titin_transition_theta = 0.5;
    // Seconds for which a just-detached titin may rebind its previous actin
    // by molecular identity and distance before returning to global search.
    double explicit_titin_rebind_memory = 0.0;
    // Intramolecular effective-concentration factor during that memory window.
    // The elastic Boltzmann penalty is still applied after this multiplier.
    double explicit_titin_rebind_rate_multiplier = 1.0;
    std::vector<ExplicitTitin> explicit_titins;

    // Positive-control Z-disc register, applied to the two intended AA pairs.
    double register_k = 0.0;             // pN / um
    double register_target = 0.1;        // um
    double register_force_cap = 0.0;     // pN; <=0 means no cap
    // The register coordinate is the signed separation of the two left
    // endpoints.  When enabled, apply the matching orientational gradient as
    // well as the center force so this diagnostic spring is conservative for
    // rotating rods.
    bool register_endpoint_torque = false;

    // Phenomenological cross-midpoint feedback.  It turns on only after the
    // AA overlap exceeds `midpoint_target`. Positive sign restores toward the
    // target; negative sign assists further sliding.
    double midpoint_feedback_k = 0.0;     // pN / um
    double midpoint_target = 0.1;         // um
    double midpoint_force_cap = 0.0;      // pN
    double midpoint_feedback_sign = 1.0;  // +1 restoring, -1 assisting
};

struct Diagnostics {
    double length_before = 0.0;
    double length_after = 0.0;
    double q_am = 0.0;
    double q_aa = 0.0;
    double q_titin = 0.0;
    double q_register = 0.0;
    double q_midpoint = 0.0;
    double q_external = 0.0;
    double q_total = 0.0;
    double dlength_dt = 0.0;
    double net_added_force = 0.0;
    int titin_bound_links = 0;
    int titin_candidate_links = 0;
    int titin_tension_links = 0;
    int titin_compressed_links = 0;
    int titin_released_links = 0;
    int titin_newly_bound_links = 0;
    double titin_mean_distance = 0.0;
    double titin_total_force = 0.0;
    double titin_attachment_energy_total = 0.0;
    double titin_attachment_energy_max = 0.0;
    double titin_attachment_force_max = 0.0;
};

Mode parse_mode(const std::string& name);
const char* mode_name(Mode mode);

double wrap_delta(double value, double length);
double pair_generalized_force(double force_on_first_x,
                              double first_x,
                              double second_x,
                              double length);
double integrate_length(double length, double generalized_force,
                        double mobility, double dt,
                        double min_length, double max_length);
double tension_only_extension(double distance, double rest_length);
double tension_only_force(double distance, double rest_length,
                          double stiffness);
double tension_only_energy(double distance, double rest_length,
                           double stiffness);
double strain_gated_association_rate(double base_rate, double beta,
                                     double theta, double elastic_energy);
double strain_gated_dissociation_rate(double base_rate, double beta,
                                      double theta, double elastic_energy);

// Construct material thick-filament sites. Calling this does not bind actin
// and does not create or query a Z-disc-like object.
void initialize_explicit_titins(const Sarcomere& model, Config& config);
vec explicit_titin_thick_anchor(const Sarcomere& model,
                                const ExplicitTitin& titin);

double aa_overlap(const Sarcomere& model, int first, int second);
double mean_aa_overlap(const Sarcomere& model);
double mean_linked_aa_overlap(const Sarcomere& model);
double mean_am_overlap(const Sarcomere& model);

// Adds only the opt-in diagnostic forces and measures all cell-force
// components at the current state. Fixed mode is a strict no-op.
Diagnostics prepare_step(Sarcomere& model, Config& config);

// Advances and affinely rescales the cell using a previously prepared force.
// Fixed mode is a strict no-op.
void advance_cell(Sarcomere& model, const Config& config,
                  double dt, Diagnostics& diagnostics);

// Public for unit tests and for controlled initialization.
void set_axial_length(Sarcomere& model, double new_length);

}  // namespace axial_cell

#endif
