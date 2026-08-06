#ifndef SARCOMERE_H
#define SARCOMERE_H

#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <tuple>
#include <array>
#include <omp.h>
#include <mutex>
#include <cstddef>
#include <cstdint>
#include <unordered_map>

#include "components.h"
#include "utils.h"
#include "interaction.h"
#include "geometry.h"
#include "h5_utils.h"
#include "neighborlist.h"
#include "am_motor_kinetics.h"
#include "actn2_load_kinetics.h"
#include "aa_neighbor_identity.h"
#include "rod_mobility.h"

// Alias for convenience.
using vector = std::vector<double>;
const double EPS = 1e-6;

class Sarcomere {
public:
    using BondFlag = std::uint8_t;
    using RecoveryStep = std::uint32_t;
    struct AMInteractionState {
        double myosin_binding_ratio = 0.0;
        double partial_binding_ratio = 0.0;
    };
    struct AMHeadState {
        BondFlag state = static_cast<BondFlag>(am_motor_kinetics::State::detached);
        int bound_actin = -1;
        int8_t half = -1;
        double anchor_axial = 0.0;
        double anchor_azimuth = 0.0;
        double binding_azimuth = 0.0;
        // Fixed material coordinate on the bound actin (um).  Despite the
        // legacy field name, this is not a continually advancing motor target.
        double target_axial = 0.0;
        double stroke_progress = 0.0;
    };
    struct AACrosslink {
        std::uint64_t id = 0;
        int actin_i = -1;
        int actin_j = -1;
        double s_i = 0.0;
        double s_j = 0.0;
        BondFlag status = 1;
        RecoveryStep attach_step = 0;
        RecoveryStep status2_memory_until = 0;
        double support_ema = 0.0;   // motor-timescale average of productive opposed load
        bool active = true;
    };
    enum class AMTransitionReason : int {
        attachment = 0,
        weak_detachment = 1,
        powerstroke = 2,
        strong_detachment = 3,
        geometric_rupture = 4,
        neighbor_loss = 5,
        capacity_limit = 6,
    };
    enum AMRetentionFailure : int {
        am_retention_ok = 0,
        am_not_in_neighbor_list = 1 << 0,
        am_beyond_rupture_distance = 1 << 1,
        am_wrong_productive_half = 1 << 2,
        am_no_productive_axial_overlap = 1 << 3,
        am_association_gate_closed = 1 << 4,
        am_unexpected_candidate_loss = 1 << 5,
        am_head_beyond_retention_reach = 1 << 6,
    };
    struct AMRetentionDiagnostic {
        int failure_mask = am_retention_ok;
        double segment_distance = -1.0;
        double rupture_distance = -1.0;
        double actin_myosin_dot = 0.0;
        int expected_half = 0;
        int head_half = 0;
        double productive_axial_overlap = 0.0;
        double binding_ratio = 0.0;
        double partial_ratio = 0.0;
        double head_connector_distance = -1.0;
        double head_surface_distance = -1.0;
        double head_excess_extension = -1.0;
        double strain_energy = 0.0;
        double capture_reach = -1.0;
        double retention_reach = -1.0;
    };

    // Public Data Members
    Filament actin;
    Myosin myosin;
    std::array<bool,3> is_periodic{ {true, true, true} };
    NeighborList neighbor_list;
    utils::MoleculeConnection myosinIndicesPerActin;
    utils::MoleculeConnection actinIndicesPerMyosin;
    utils::MoleculeConnection actinIndicesPerActin;

    std::vector<std::vector<BondFlag>> actin_actin_bonds, actin_actin_bonds_prev;
    std::vector<std::vector<BondFlag>> actin_actin_status, actin_actin_status_prev;
    // Track lifetime (in steps) for each actin–actin catch bond
    std::vector<std::vector<int>> actin_actin_lifetime, actin_actin_lifetime_prev;
    // Track recovery windows and attachment timestamps for broken bonds
    std::vector<std::vector<RecoveryStep>> actin_recovery_until;
    // End step of the optional status-2 retention window. Mechanical support
    // is still required to establish status 2; this only bridges short losses
    // for an already-bonded status-2 AA pair and does not create residual load.
    std::vector<std::vector<RecoveryStep>> aa_status2_memory_until;
    std::vector<std::vector<int>> aa_attach_step;
    // Track actin–myosin bonds
    std::vector<std::vector<BondFlag>> am_bonds, am_bonds_prev;
    vector box;
    double k_aa, kappa_aa, k_on, k_off,
           kappa_am, k_am, v_am, myosin_radius_ratio,
           skin_distance, cutoff_radius, dt, base_lifetime, directional_base_lifetime,
           lifetime_coeff, diff_coeff_ratio;
    double am_cutoff, am_optimal;
    double actin_radius = 0.015;
    double am_surface_gap = 0.015;
    double am_steric_k = 1000.0;
    double actin_steric_k = 10000.0;
    double myosin_steric_k = 1000.0;
    double six_ring_angle_noise_deg = 20.0;
    double six_ring_radial_noise = 0.01;
    double six_ring_axis_noise_deg = 0.0;
    double aa_cutoff, aa_optimal;
    double aa_antiparallel_overlap = 0.1;
    // Initial Z-disc interdigitation used only by the minimal rigs; <=0
    // reproduces the perfectly registered start.
    double init_aa_overlap = -1.0;
    // Z-disc width: ACTN2 binds only within this distance of a filament's
    // barbed end. <=0 restores the old rule (bind anywhere in the mutual
    // overlap). Cardiac Z-discs are 100-140 nm, so ~0.1 um.
    double aa_barbed_zone = -1.0;
    // Optional per-molecule radial stiffness for weak/status-1 ACTN2 links.
    // A non-positive value preserves the legacy behavior in which status 1
    // and status 2 use the same fixed per-link stiffness.
    double aa_status1_radial_k = -1.0;
    double aa_axial_register_ratio = 1.0;
    // Exploratory switch: make the passive barbed-end register mechanical in
    // every attached ACTN2 state.  Status still controls catch-bond lifetime.
    // Default false preserves the production status-2-only behavior.
    bool aa_register_all_statuses = false;
    // >0 uses the saturating barbed-end register force with this ceiling
    // (pN) instead of the Gaussian form that decays to zero far from the
    // target.  Small-error stiffness is unchanged either way.
    double aa_axial_register_force_max = 0.0;
    double aa_status2_support_window = 0.0;  // s; <=0 keeps the instantaneous rule
    double aa_status2_registry_on = 0.1;
    double aa_status2_registry_off = 0.15;
    double aa_register_saturation_error = 0.05;

    double aa_link_register_ratio(BondFlag status) const {
        return (aa_register_all_statuses || status == 2)
            ? aa_axial_register_ratio : 0.0;
    }
    double aa_link_radial_stiffness(BondFlag status) const {
        if (status == 1 && aa_status1_radial_k > 0.0) {
            return aa_status1_radial_k;
        }
        return actn2_load_kinetics::fixed_single_link_stiffness(
            k_aa, aa_stiffness_reference_links);
    }
    bool am_bidirectional_crossbridge = false;
    double am_overlap_target_ratio = 0.0;  // legacy sliding throttle; <=0 disables
    // Sparse, spatially resolved ACTN2 links are the production model. The
    // old one-edge-per-pair model remains available through the CLI control.
    bool aa_local_multilink = true;
    int aa_prebundle_links_per_pair = 2;
    double aa_site_exclusion = 0.0192;
    double aa_overlap_capture = 0.1;
    double aa_perpendicular_k = 50.0;
    // Optional molecular polarity gate for new local ACTN2 associations.
    // A value >= 1 disables the gate; otherwise the two actin directors must
    // have dot product below this threshold. This constrains binding chemistry,
    // not neighbor identity or any preassembled lattice topology.
    double aa_association_max_direction_dot = 1.0;
    // Finite ACTN2 contour length. A positive value forces an overstretched
    // local link to detach and immediately releases its remembered neighbor
    // slot so another nearby antiparallel partner can be sampled.
    double aa_link_rupture_length = 0.0;
    int aa_stiffness_reference_links = 3;
    int aa_total_link_pool = 0;
    std::uint64_t aa_next_link_id = 1;
    std::vector<AACrosslink> aa_links;
    std::vector<size_t> aa_free_link_slots;
    std::vector<size_t> aa_active_link_slots;
    std::vector<std::vector<size_t>> aa_links_by_actin;
    std::unordered_map<std::uint64_t, std::vector<size_t>> aa_links_by_pair;
    // Scratch buffers for _update_local_aa_links().  They were function locals
    // while that routine ran inside an `omp single`; the orphaned worksharing
    // version needs every thread to see the SAME buffer, so they live on the
    // object and are reused (capacity is retained) across steps.
    std::vector<LocalAALinkResult> aa_scratch_mechanics;
    std::vector<actn2_load_kinetics::ArmLoads> aa_scratch_arm_loads;
    std::unordered_map<std::uint64_t, double> aa_scratch_pair_tension;
    std::vector<double> aa_scratch_total_arm_load;
    std::vector<double> aa_scratch_productive_load;
    // Per-actin motor-force vectors are invariant across the whole routine
    // (nothing in it touches AM state), so they are evaluated once per actin
    // instead of twice per incident link.
    std::vector<std::array<double, 3>> aa_scratch_motor_force;
    std::vector<std::array<double, 3>> aa_scratch_productive_force;
    // Association candidates, evaluated in parallel against the pre-commit
    // state and then committed serially in canonical order.
    struct AAAssociationProposal {
        double s_i = 0.0;
        double s_j = 0.0;
        bool fire = false;
    };
    std::vector<std::pair<int, int>> aa_scratch_pairs;
    std::vector<AAAssociationProposal> aa_scratch_proposals;
    std::vector<char> aa_scratch_dirty_actin;
    std::vector<std::vector<int>> aa_scratch_sorted_neighbors;
    std::vector<size_t> aa_scratch_pair_offset;
    size_t aa_scratch_pair_count = 0;
    // Neighbor identity is a topological layer separate from molecular ACTN2
    // occupancy.  Short detachments therefore do not instantly replace a
    // lattice edge; status-2 support protects the vacated slot for longer.
    bool aa_neighbor_identity_enabled = true;
    size_t aa_neighbor_slot_memory_steps = 0;
    size_t aa_status2_slot_memory_steps = 0;
    aa_neighbor_identity::Registry aa_neighbor_slots;
    double aa_spacing_ratio = 0.92;
    double am_preferred_centerline = 0.0;
    double aa_capture_margin = 0.02;
    double k_mm = 0.0;
    double myomesin_optimal = 0.0;
    double myomesin_cutoff = 0.0;
    double titin_k, titin_rest_length;
    double k_bundle_max = 0.0;
    int bundle_ramp_steps = 0;
    size_t interaction_ramp_start_step = 0;
    size_t interaction_ramp_steps = 0;
    double max_actin_force, max_myosin_force, max_actin_torque, max_myosin_torque;
    bool physical_rod_mobility = false;
    bool isometric_actin_axial = false;
    double temperature_kelvin = 310.0;
    double thermal_beta = 0.0;
    // Integrator settings copied into trajectory metadata.  Tensor mobility
    // lives below; these retain the legacy fallback and Brownian-noise scale
    // so a diagnostic file fully describes how motion was generated.
    double legacy_actin_D_trans = 0.0;
    double legacy_actin_D_rot = 0.0;
    double legacy_myosin_D_trans = 0.0;
    double legacy_myosin_D_rot = 0.0;
    double thermal_noise_scale_metadata = 1.0;
    rod_mobility::Diffusion actin_rod_mobility;
    rod_mobility::Diffusion myosin_rod_mobility;
    size_t bond_recovery_steps;
    bool directional;
    bool stabilize_status1 = false;
    actn2_load_kinetics::StatusMode aa_status_mode =
        actn2_load_kinetics::StatusMode::legacy_motor_identity;
    double aa_status2_force_on = 1.0;
    double aa_status2_force_off = 0.5;
    // Explicit cross-bridge chemomechanics. main.cpp enables this by default;
    // the false member initializer keeps direct library construction benign
    // until set_am_kinetics() supplies the run configuration.
    bool dynamic_am = false;
    bool am_offcenter_torque = true;
    bool am_spatial_heads = true;
    bool sterics_enabled = true;
    int am_heads_per_half = 1;
    // Maximum weak+strong heads that one fully overlapping actin may occupy
    // on one myosin half. Cardiac thick filaments have about six surrounding
    // thin-filament neighbors, so 24 shares the default 141-head pool without
    // reducing the conserved half-filament total. Zero disables the cap.
    int am_max_heads_per_actin_half = 24;
    double am_k_on = 100.0;
    double am_weak_binding_floor = 0.1;
    double am_k_weak_off = 200.0;
    double am_k_powerstroke = 200.0;
    double am_k_strong_off = 160.0;
    double am_strong_off_floor_fraction = 1.0;
    double am_force_scale = 1.0;
    // Cardiac cross-bridge stiffness in model units (pN/um). 1000 pN/um is
    // 1 pN/nm, within the single-molecule beta-cardiac range.
    double am_motor_k = 1000.0;
    double am_power_stroke = 0.005;
    double am_unloaded_velocity = 5.0;
    double am_stall_force = 4.5;
    double am_weak_stiffness_ratio = 0.1;
    // Pair-level flat-bottom containment: force-free near the preferred AM
    // spacing, stiff only outside the capture band so transverse drift cannot
    // masquerade as chemical or axial-overlap dissociation.
    double am_tether_k = 3000.0;
    double am_tether_slack = 0.002;
    // Passive aggregate angular compliance of an occupied AM pair. It is
    // force-free inside the slack angle and does not replace active,
    // off-center head torque.
    double am_alignment_k = 30.0;
    double am_alignment_slack = 2.0 * M_PI / 180.0;
    // Pair-level transverse stiffness grows with occupied heads and saturates
    // at this population instead of appearing/disappearing with the last head.
    int am_tether_saturation_heads = 6;
    bool am_binding_azimuth_selection = true;
    double am_binding_azimuth_half_width = M_PI / 6.0;
    // Emergency contour-extension guard. Ordinary loss is sampled through
    // strain-dependent KMC rather than imposed at this boundary.
    double am_rupture_distance = 0.1;
    bool am_strain_dependent_detachment = true;
    double am_strain_detachment_theta = 0.1;
    double am_strain_detachment_max_exponent = 20.0;
    double am_cycle_period = 0.0;
    double am_activation_rise = 0.02;
    double am_activation_decay = 0.2;
    double am_activation_floor = 0.0;
    double am_activation_active_fraction = 0.5;
    double am_activation = 1.0;
    bool am_thin_activation_enabled = true;
    bool am_thin_feedback_enabled = true;
    double am_q_calcium_open_rate = 20.0;
    double am_q_calcium_half = 0.5;
    double am_q_calcium_hill = 3.0;
    double am_q_strong_feedback_rate = 100.0;
    double am_q_strong_half = 0.02;
    double am_q_strong_hill = 2.0;
    double am_q_close_rate = 20.0;
    std::vector<double> am_thin_activation;
    bool am_thick_feedback_enabled = false;
    double am_thick_recruit_rate = 2.0;
    double am_thick_fold_rate = 20.0;
    double am_thick_tension_scale = 0.1;
    double am_thick_initial_on_fraction = 0.1;
    double am_thick_force_clamp = -1.0;
    std::vector<std::array<double, 2>> am_half_tension;
    std::vector<vec> myosin_material_e1;
    std::vector<vec> myosin_material_e2;
    size_t am_status2_memory_steps = 0;
    std::vector<std::vector<BondFlag>> am_state;
    std::vector<std::vector<BondFlag>> am_seen;
    std::vector<std::vector<double>> am_target_axial;
    // Physical heads are owned by a thick-filament half, not duplicated for
    // every geometrically neighboring actin. Index 0 is the negative axial
    // half and index 1 is the positive axial half.
    std::vector<std::array<std::vector<AMHeadState>, 2>> am_heads_by_myosin;
    int fix_myosin;
    int max_myosin_bonds;
    int max_strong_actin_bonds;
    int max_aa_bonds_per_actin = 0;
    bool aa_energy_gated_association = false;
    // Provisional finite-time kinetic compromise from the passive six-ring
    // sweep. This has no effect unless energy gating or the passive DB control
    // is enabled.
    double aa_transition_state_theta = 0.1;
    bool aa_passive_detailed_balance = false;
    // The unbound encounter complex is much softer and less orientationally
    // selective than a mature two-arm ACTN2 bridge.
    double aa_encounter_k = 10.0;
    double aa_encounter_kappa = 0.0;
    std::vector<std::vector<AMInteractionState>> am_interaction;
    vector actin_crosslink_ratio;
    std::vector<vec> actin_crosslink_start;
    std::vector<vec> actin_crosslink_end;
    std::vector<int> actin_n_bonds;
    std::vector<int> actin_strong_cb_count;
    std::vector<int> n_myosins_per_actin;
    std::vector<std::pair<std::vector<int>, std::vector<int>>> actin_neighbors_by_species;
    vector actin_basic_tension;
    std::vector<bool> actin_f_load_computed;
    std::vector<double>* actin_f_load_cb;
    std::vector<std::mutex> actin_f_load_mutex;
        
    gsl_rng* rng;
    std::string filename;
    int initial_seed = 0;

    std::vector<std::vector<vec>> actin_forces_temp,
                                    myosin_forces_temp, myosin_velocities_temp, actin_torques_temp, myosin_torques_temp;
    std::vector<std::vector<int>> actin_cb_status_temp;
    std::vector<std::vector<std::array<double, 2>>> myosin_f_load_temp;
    std::vector<std::vector<double>> cb_breakage_events_temp;
    std::vector<std::vector<double>> aa_completed_lifetimes_temp;
    std::vector<std::vector<double>> am_transition_events_temp;
    std::vector<std::vector<double>> aa_state_events_temp;
    std::vector<std::vector<double>> actin_chem_entropy_delta_temp;
    std::vector<std::vector<double>> actin_chem_input_delta_temp;
    std::vector<std::vector<double>> actin_chem_binding_delta_temp;
    std::vector<std::vector<double>> actin_chem_unbinding_delta_temp;
    std::vector<std::array<double, 2>> myosin_f_load;
    struct AMHeadCandidate {
        int actin_index = -1;
        int half_index = 0;
        double geometry_gate = 0.0;
        double productive_overlap = 0.0;
        double binding_ratio = 0.0;
        double partial_ratio = 0.0;
    };
    std::vector<std::vector<AMHeadCandidate>> am_head_candidates_;
    std::vector<utils::MoleculeConnection> actinIndicesPerMyosin_temp;
    std::vector<gsl_rng*> rng_engines;
    std::vector<std::vector<BondFlag>> myosin_bond_matrix;
    bool has_myosin_bond_pairs = false;

    // Compact indexes into the dense compatibility matrices. These keep the
    // per-step reset and checkpoint paths proportional to active pairs.
    using BondPair = std::pair<int, int>;
    std::vector<BondPair> active_aa_pairs;
    std::vector<BondPair> previous_aa_pairs;
    // Persistent buffer for _record_aa_state_changes(); see there.
    std::vector<BondPair> aa_state_change_candidates;
    std::vector<BondPair> active_am_pairs;
    std::vector<BondPair> previous_am_pairs;
    std::vector<BondPair> active_recovery_pairs;
    std::vector<std::vector<BondPair>> active_aa_pairs_temp;
    std::vector<std::vector<BondPair>> active_am_pairs_temp;
    std::vector<std::vector<BondPair>> recovery_pairs_temp;

    // Record the global simulation step and catch-bond breakage events
    size_t current_step = 0;
    // Flat buffer storing (i, j, step, distance, cos_angle) for each breakage
    std::vector<double> cb_breakage_events;

    // Flat buffer storing (i, j, step, bond_count_i, bond_count_j) for
    // removals triggered by the max_strong_actin_bonds limit
    std::vector<double> cb_limit_events;
    // Completed lifetimes sampled at detachments
    std::vector<double> aa_completed_lifetimes;
    // Rows: actin, myosin, step, old_state, new_state, reason, activation,
    // axial_strain, active_motor_force, preceding_actin_cb_status, head_id.
    std::vector<double> am_transition_events;
    // Rows: actin, myosin, step, head_id, old_state, broad_reason,
    // retention_failure_mask, segment_distance, rupture_distance,
    // actin_myosin_dot, expected_half, head_half, productive_axial_overlap,
    // binding_ratio, partial_ratio. One row is emitted for every head removed
    // by a pair-level forced rupture rather than an ordinary KMC transition.
    std::vector<double> am_forced_rupture_events;
    // The 15 legacy fields above followed by exact material-head connector
    // distance, head-to-actin surface distance, association reach, retention
    // reach. This distinguishes pair-neighbor loss from true per-head rupture.
    std::vector<double> am_forced_rupture_events_v2;
    // Rows: myosin, half_sign, head_id, step, old_state, new_state, tension.
    std::vector<double> am_thick_transition_events;
    // Rows: actin_i, actin_j, step, old_bond, new_bond, old_status,
    // new_status, motor_force_i/j, strong_count_i/j, load_i/j.
    std::vector<double> aa_state_events;
    // Rows: actin_i, actin_j, step, old_status, new_status,
    // motor_only_load, opposed_drive, transmitted_load, combined_candidate,
    // selected_classifier_load, polarity_dot. Schema version is stored in /meta.
    std::vector<double> aa_load_events;
    // Rows: link_id, actin_i, actin_j, step, event_code (1 attach,
    // 2 detach, 3 status), old_status, new_status, tensile_load,
    // kinetic_lifetime, s_i, s_j.
    std::vector<double> aa_local_link_events;
    // Backward-compatible extended rows: the 11 fields above followed by
    // signed arm_i, signed arm_j, bilateral effective load, productive motor
    // support, instantaneous gate, and legacy scalar connector tension.
    std::vector<double> aa_local_link_events_v2;
    std::vector<int> previous_actin_cb_status;
    std::vector<int> aa_reserved_degree;

    // Constructors & Destructor
    Sarcomere();
    Sarcomere(int& n_actins, int& n_myosins, vector box0, double& actin_length, double& myosin_length,
        double& myosin_radius, double& myosin_radius_ratio, double& aa_cutoff, double& aa_optimal,
        double& k_on,
        double& base_lifetime, double& directional_base_lifetime, double& lifetime_coeff,
        double& diff_coeff_ratio, double& k_aa, double& kappa_aa,
        double& k_am, double& kappa_am, double& k_mm, double& v_am, std::string& filename, gsl_rng* rng, int& seed,
        int& fix_myosin, double& dt, bool& directional, std::string& boundary_condition,
        int max_myosin_bonds, int max_strong_actin_bonds, double max_actin_force_param,
        double max_myosin_force_param, double max_actin_torque_param, double max_myosin_torque_param);
    Sarcomere(int& n_actins, int& n_myosins, vector box0, double& actin_length, double& myosin_length,
        double& myosin_radius, double& am_cutoff, double& am_optimal, double& aa_cutoff, double& aa_optimal,
         double& k_on, double& k_off,
        double& base_lifetime, double& directional_base_lifetime, double& lifetime_coeff,
        double& diff_coeff_ratio, double& k_aa, double& kappa_aa, double& k_am, double& kappa_am, double& k_mm, double& v_am,
        std::string& filename, gsl_rng* rng, int& seed, int& fix_myosin, double& dt, double tau_rec,
        double titin_k, double titin_rest_length, bool& directional, int max_myosin_bonds,
        double max_actin_force_param, double max_myosin_force_param,
        double max_actin_torque_param, double max_myosin_torque_param);
    ~Sarcomere();

    // Public Methods
    void myosin_on_a_lattice();
    void partial_fix(int& n_fixed_myosins);
    void cb();
    void minimal_catch();
    void minimal_catch_variant(const std::string& variant);
    void aa_bundle_common();
    void aa_bundle_opposed();
    void bad_cb();
    void cb_off_angle();
    void am_off_angle();
    void set_actin_direction_x_noise(double noise_std);
    void set_myosin_direction_x_noise(double noise_std,
                                       double tilt_deg = 0.0);
    void single_am();
    void single_am_tilted();
    void six_actin_ring();
    void prebind_six_actin_ring();
    void bilateral_six_actin_ring();
    void prebind_bilateral_six_actin_ring();
    void seed_balanced_six_actin_heads(int heads_per_actin = 6);
    void force_balanced_six_actin_cell();
    void prebind_force_balanced_six_actin_cell();
    void seed_force_balanced_six_actin_heads(int heads_per_actin = 6);
    void shared_thick_thin_lattice();
    void prebind_shared_thick_thin_lattice(int links_per_edge = 1);
    void seed_shared_thick_thin_heads(int heads_per_pair = 6);
    void prebind_local_aa_pair(int first, int second, int link_count);
    void prebind_assembly_actin_links();
    void minimal_sarcomere();
    void two_sarcomere_unit();
    void three_membered_lattice(bool prebind_actn2 = true,
                                double transverse_jitter = 0.0,
                                double axial_jitter = 0.0,
                                double tilt_jitter = 0.0);
    void seed_three_membered_heads(int heads_per_actin = 6);
    void sarcomeric_structure();
    void sarcomeric_structure_tight();
    void set_bundling_parameters(double max_strength, int ramp_steps);
    void set_interaction_ramp(size_t start_step, size_t ramp_steps);
    void set_steric_stiffness(double actin_stiffness,
                              double myosin_stiffness);
    void set_am_kinetics(bool enabled, double k_on_am, double k_weak_off,
                         double k_powerstroke, double k_strong_off,
                         double force_scale, double motor_k, double power_stroke,
                         double stall_force, double unloaded_velocity,
                         double weak_stiffness_ratio, double rupture_distance,
                         double cycle_period, double activation_rise,
                         double activation_decay, double activation_floor,
                         double status2_memory_time,
                         double activation_active_fraction = 0.5);
    void set_am_head_count(int heads_per_half);
    void set_am_neighbor_head_cap(int full_overlap_cap);
    void set_am_weak_binding_floor(double basal_fraction);
    void set_am_proximity_tether(double stiffness, double slack,
                                 bool select_binding_azimuth,
                                 double azimuth_half_width_deg,
                                 double alignment_stiffness = 0.0,
                                 double alignment_slack_deg = 2.0);
    void set_am_retention_kinetics(bool strain_dependent,
                                   double barrier_fraction,
                                   double maximum_exponent,
                                   int tether_saturation_heads,
                                   double strong_off_floor_fraction = 1.0);
    void set_am_thin_activation(bool enabled, bool feedback_enabled,
                                double calcium_open_rate, double calcium_half,
                                double calcium_hill, double strong_feedback_rate,
                                double strong_half, double strong_hill,
                                double close_rate);
    void set_am_thick_regulation(bool enabled, double recruit_rate,
                                 double fold_rate, double tension_scale,
                                 double initial_on_fraction,
                                 double force_clamp);
    void set_am_radial_geometry(double actin_radius_value,
                                double preferred_surface_gap,
                                double steric_stiffness);
    void set_aa_spacing_geometry(double ratio,
                                 double am_centerline,
                                 double capture_margin);
    void set_aa_axial_overlap(double overlap, double register_ratio = 1.0,
                              double registry_on = 0.1,
                              double registry_off = 0.15,
                              double saturation_error = 0.05);
    void set_local_aa_model(bool enabled, int links_per_pair,
                            double site_exclusion, double overlap_capture,
                            double perpendicular_stiffness,
                            int total_link_pool,
                            int stiffness_reference_links = 3);
    void set_aa_neighbor_identity(bool enabled, double status1_memory,
                                  double status2_memory);
    bool reserve_aa_neighbor_pair(int first, int second,
                                  bool pinned = false);
    int aa_neighbor_degree(int actin_index) const;
    void set_physical_mobility(double temperature,
                               double beta_value,
                               const rod_mobility::Diffusion& actin_mobility,
                               const rod_mobility::Diffusion& myosin_mobility);
    void transport_myosin_material_frame(int myosin_index,
                                         const vec& old_axis,
                                         const vec& new_axis,
                                         double axial_rotation);
    void reset_myosin_material_frames();
    void set_aa_status_kinetics(actn2_load_kinetics::StatusMode mode,
                                double force_on, double force_off);
    void update_system();
    void update_system_sterics_only();
    void set_periodicity(const std::array<bool,3>& periodic_axes);
    void new_file();
    void save_state();
    void save_resume_snapshot();
    int load_state(int& n_frames, int frame_index = -1);
    // Debug helper to compute catch-bond statistics for a single frame
    void debug_cb_stats();

private:
    double _current_interaction_ramp() const;
    // Private helper methods
    void _update_neighbors();
    void _set_to_zero();
    void _finalize_active_pair_lists();
    void _rebuild_active_pair_lists_from_dense();
    void _process_actin_myosin_binding(int& i);
    void _process_dynamic_am_heads();
    void _update_thin_filament_activation();
    void _drop_unseen_am_bonds();
    void _initialize_am_head_pools();
    double _am_head_association_gate(int actin_index, int myosin_index,
                                     const AMHeadState& head,
                                     double binding_azimuth) const;
    bool _am_head_retention_reachable(int actin_index, int myosin_index,
                                      const AMHeadState& head) const;
    double _am_head_strain_energy(int actin_index, int myosin_index,
                                  int head_index,
                                  am_motor_kinetics::State state) const;
    double _am_head_detachment_rate(int actin_index, int myosin_index,
                                    int head_index,
                                    am_motor_kinetics::State state,
                                    double chemical_rate) const;
    int _am_bound_head_count(int actin_index, int myosin_index) const;
    double _am_per_head_tether_stiffness(int bound_head_count) const;
    AMHeadState* _am_head(int myosin_index, int head_id);
    const AMHeadState* _am_head(int myosin_index, int head_id) const;
    int _am_head_id(int half_index, int head_index) const;
    double _am_axial_coordinate(int actin_index, int myosin_index) const;
    double _am_productive_strain(int actin_index, int myosin_index,
                                 int head_index = -1) const;
    double _am_head_active_force(int actin_index, int myosin_index,
                                 int head_index) const;
    double _am_active_force(int actin_index, int myosin_index) const;
    double _am_resisting_force(int actin_index, int myosin_index) const;
    std::array<double, 3> _actin_motor_force_vector(int actin_index) const;
    std::array<double, 3> _actin_productive_motor_force_vector(
        int actin_index) const;
    double _aa_productive_opposed_load(int i, int j) const;
    double _aa_transmitted_load(int i, int j);
    void _record_am_transition(int actin_index, int myosin_index, int head_index,
                               am_motor_kinetics::State old_state,
                               am_motor_kinetics::State new_state,
                               AMTransitionReason reason);
    AMRetentionDiagnostic _am_retention_diagnostic(
        int actin_index, int myosin_index, int head_index,
        bool candidate_missing) const;
    void _record_am_forced_rupture(
        int actin_index, int myosin_index, int head_index,
        am_motor_kinetics::State old_state, AMTransitionReason reason,
        bool candidate_missing);
    void _record_am_thick_transition(int myosin_index, int half_index,
                                     int head_id,
                                     am_motor_kinetics::State old_state,
                                     am_motor_kinetics::State new_state,
                                     double half_tension);
    void _sync_am_pair_summary(int actin_index, int myosin_index);
    void _detach_am_heads_for_pair(int actin_index, int myosin_index,
                                   AMTransitionReason reason);
    void _process_catch_bonds(int& i);
    void _update_local_aa_links();
    std::uint64_t _aa_pair_key(int i, int j) const;
    std::pair<double, double> _local_aa_attachment_candidate(
        int i, int j, int ordinal, double& overlap) const;
    // Per-pair axial geometry shared by every candidate ordinal of a pair.
    struct AAPairAxis {
        vec axis;
        double pi = 0.0;
        double pj = 0.0;
        double lower = 0.0;
        double overlap = 0.0;
        double denominator_i = 0.0;
        double denominator_j = 0.0;
    };
    AAPairAxis _local_aa_pair_axis(int i, int j) const;
    std::pair<double, double> _local_aa_attachment_site(
        const AAPairAxis& geometry, int ordinal, double& overlap) const;
    int _local_aa_distinct_neighbors(int actin_index) const;
    bool _local_aa_site_available(int actin_index, int partner,
                                  double coordinate) const;
    AAAssociationProposal _evaluate_aa_association(int i, int j) const;
    bool _local_aa_site_outside_am_overlap(int actin_index,
                                           double axial_coordinate) const;
    void _rebuild_local_aa_indices();
    void _record_aa_state_changes();
    void _calc_am_force_velocity(int& i);
    void _apply_titin_forces(int& i);
    void _apply_myomesin_spring(int i, int j, std::vector<vec>& local_myosin_forces);
    void _volume_exclusion();
    void _myosin_exclusion();
    void _myosin_repulsion(int& i, int& j);
    void _actin_repulsion(int& i, int& j);
    int determine_cb_status(int& i, int& j);
    actn2_load_kinetics::ArmLoads _aa_signed_arm_loads(int i, int j);
    actn2_load_kinetics::ArmLoads _aa_effective_arm_loads(int i, int j);
    double _actin_motor_load(int actin_index) const;
    bool _cb_decide(int& i, int& j, int status);
    void _record_chem_transition(int& i, int& j, double affinity, bool binding);
    void compute_actin_f_load(int& i);
    void _set_cb(int& i, int& j, int status);
    void _set_cb(int& i, std::vector<int> indices, std::vector<int> status);
    vec _alignment_torque(const vec& u, double k_bias);
    void _apply_cb_alignment_bias(double& k_theta_bias);
    void _update_myosin_bond_matrix();
    bool _myosin_pair_bonded(int mi, int mj) const;
    double _current_bundle_strength() const;
    void _apply_transverse_bundling(double k_bundle);
    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
        _extract_bonded_pairs(
        const std::vector<std::vector<BondFlag>>& actin_actin_bonds,
        const std::vector<std::vector<BondFlag>>& actin_actin_status,
        const utils::MoleculeConnection& myosinIndicesPerActin);
    void _enforce_actin_cb_limit();
    void _enforce_myosin_bond_limit();
};

#endif // SARCOMERE_H
