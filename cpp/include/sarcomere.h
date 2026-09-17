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
#include "actn2_arm_kinetics.h"
#include "aa_neighbor_identity.h"
#include "rod_mobility.h"
#include "titin_kinetics.h"

// Alias for convenience.
using vector = std::vector<double>;
const double EPS = 1e-6;

class Sarcomere {
public:
    using BondFlag = std::uint8_t;
    using RecoveryStep = std::uint32_t;
    enum class TitinMode : int {
        // Historical occupancy-gated spring through a status-2 AA partner.
        legacy_crosslinked = 0,
        // Balanced four-actin/two-myosin cell: each actin barbed end is
        // tethered reciprocally to the center of its owning thick filament.
        balanced_owning = 1,
        // Explicit material domains on both thick-filament halves bind
        // individual, geometrically compatible actin barbed ends.  Links
        // turn over and are exactly slack under compression.
        explicit_dynamic = 2,
    };
    struct ExplicitTitin {
        int myosin = -1;
        int half = -1;
        int track = -1;
        double anchor_azimuth = 0.0;
        int bound_actin = -1;
        int candidate_actin = -1;
        int last_bound_actin = -1;
        RecoveryStep attach_step = 0;
        RecoveryStep last_detach_step = 0;
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
    // One coarse-grained myomesin/M-band module between two thick-filament
    // bare zones.  A module may represent several parallel molecules, but it
    // remains one local graph edge and consumes one finite-valency site on
    // each myosin.  At most one module is allowed for a given myosin pair.
    struct MyomesinBridge {
        std::uint64_t id = 0;
        int myosin_i = -1;
        int myosin_j = -1;
        RecoveryStep attach_step = 0;
        double axial_offset = 0.0;
        double radial_distance = 0.0;
        double force = 0.0;
    };
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
        // Opt-in explicit two-arm kinetics.  Legacy trajectories omit this
        // field and map status 1/2 to weak/mature bridges on resume.
        actn2_arm_kinetics::State arm_state =
            actn2_arm_kinetics::State::weak_bridge;
        RecoveryStep arm_state_since = 0;
        int last_partner = -1;
        double arm_force_i = 0.0;
        double arm_force_j = 0.0;
        double contact_force = 0.0;
        double contact_penetration = 0.0;
        RecoveryStep attach_step = 0;
        RecoveryStep status2_memory_until = 0;
        double support_ema = 0.0;   // motor-timescale average of productive opposed load
        // Candidate motor identities are exchangeable only in status 1. Once
        // status 2 forms they are frozen while the mature conformation lasts;
        // optional mature-to-weak relaxation releases them for annealing.
        int candidate_myosin_i = -1;
        int candidate_myosin_j = -1;
        int8_t candidate_half_i = 0;
        int8_t candidate_half_j = 0;
        RecoveryStep candidate_since = 0;
        RecoveryStep candidate_last_seen = 0;
        double candidate_support_ema = 0.0;
        double candidate_support_instantaneous = 0.0;
        double candidate_raw_support_i = 0.0;
        double candidate_raw_support_j = 0.0;
        double candidate_allocated_support_i = 0.0;
        double candidate_allocated_support_j = 0.0;
        int challenger_myosin_i = -1;
        int challenger_myosin_j = -1;
        int8_t challenger_half_i = 0;
        int8_t challenger_half_j = 0;
        RecoveryStep challenger_since = 0;
        double challenger_support_ema = 0.0;
        bool candidate_locked = false;
        bool active = true;
    };
    struct ProductiveMotorContribution {
        int myosin = -1;
        int8_t half = 0;
        std::array<double, 3> force{0.0, 0.0, 0.0};
    };
    enum class AMTransitionReason : int {
        attachment = 0,
        weak_detachment = 1,
        powerstroke = 2,
        strong_detachment = 3,
        geometric_rupture = 4,
        neighbor_loss = 5,
        capacity_limit = 6,
        site_occupied = 7,
        actin_turnover = 8,
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
    // Range of the myosin-myosin repulsion.  0 => 2*myosin.radius (hard core
    // only).  ~0.0953 holds the thick-filament lattice at the spacing that
    // seats an actin in the trigonal interstice.
    double myosin_exclusion_distance = 0.0;
    // Optional scaffold-free short-sarcomere contact.  It acts on actual
    // antiparallel pointed-end cap geometry, independent of ACTN2/AM status.
    // It is default-off until the compressed-cell validation is complete.
    bool aa_double_overlap_contact_enabled = false;
    double aa_double_overlap_contact_diameter = 0.012;
    double aa_double_overlap_contact_epsilon = 0.00428;
    double aa_double_overlap_contact_force_cap = 0.0;
    std::uint64_t aa_double_overlap_contact_count = 0;
    double aa_double_overlap_max_penetration = 0.0;
    double aa_double_overlap_max_force = 0.0;
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
    double aa_axial_register_ratio = 1.0;
    // When enabled, every attached local ACTN2 molecule carries axial
    // register mechanics. Status 1/2 still controls molecular lifetime.
    bool aa_register_all_statuses = false;
    // Status 1 is a *sliding* bond: radial + alignment only, so a
    // crosslinked pair can still slide axially to search for register.
    // The perpendicular term penalises the connector's component along
    // each filament axis, i.e. it resists exactly that sliding, so it
    // belongs to status 2 together with the axial register term.
    // Default true preserves historical behaviour. Set false to make
    // status 1 a sliding bond. Incompatible with
    // aa_status_mode=arm_load_productive, which gates status-2
    // promotion on exactly this load.
    bool aa_perpendicular_all_statuses = true;
    // >0 uses the saturating barbed-end register force with this ceiling
    // (pN) instead of the Gaussian form that decays to zero far from the
    // target.  Small-error stiffness is unchanged either way.
    double aa_axial_register_force_max = 0.0;
    // Saturating radial laws. <=0 keeps the historical unbounded
    // quadratic. When positive, force is f_max*tanh(ext*k/f_max): the
    // small-displacement stiffness is unchanged, but the force (and so
    // the displacement per step) is bounded, which is what lets dt rise
    // without resorting to a displacement clamp. A clamp breaks detailed
    // balance; this is a real potential and does not.
    double aa_radial_force_max = 0.0;
    double am_tether_force_max = 0.0;
    // Floor on the radial pair tether's bound-head occupancy multiplier.
    // 0 reproduces the occupancy-gated behaviour exactly.
    double am_tether_occupancy_floor = 0.0;
    // Reach window past am_preferred_centerline within which the floor
    // applies.  Mandatory: the tether is uncapped by default.
    double am_tether_occupancy_floor_reach = 0.04;
    double aa_status2_support_window = 0.0;  // s; <=0 keeps the instantaneous rule
    // Identity-resolved weak-state annealing. A candidate pair may exchange
    // while status 1, must persist for commit_time, then is frozen in status 2.
    double aa_motor_pair_support_window = 0.010;
    double aa_motor_pair_commit_time = 0.010;
    double aa_motor_pair_switch_ratio = 1.5;
    double aa_motor_pair_exchange_time = 0.001;
    double aa_motor_pair_dropout_time = 0.010;
    double aa_status2_registry_on = 0.1;
    double aa_status2_registry_off = 0.15;
    double aa_register_saturation_error = 0.05;
    bool am_bidirectional_crossbridge = false;
    // Crossover pull-through: when an actin slides past the central bare
    // zone into the opposite head-bearing half, let that half power stroke
    // too. Stroke direction is set by actin polarity in
    // compute_am_two_spring_head (productive_strain = stroke_progress - z,
    // z measured along the actin axis), so those heads pull the filament
    // the same way it was already moving and carry it out, instead of
    // trapping it as nonproductive weak encounters that occupy sites.
    bool am_crossover_pull_through = false;
    double am_overlap_target_ratio = 0.0;  // legacy sliding throttle; <=0 disables
    // Sparse, spatially resolved ACTN2 links are the production model. The
    // old one-edge-per-pair model remains available through the CLI control.
    bool aa_local_multilink = true;
    int aa_prebundle_links_per_pair = 1;
    double aa_site_exclusion = 0.0192;
    double aa_overlap_capture = 0.1;
    double aa_perpendicular_k = 50.0;
    // Weak/status-1 ACTN2 is a soft encounter complex. A negative value
    // preserves the historical behavior and uses the locked-link stiffness.
    double aa_status1_radial_k = -1.0;
    // One-dimensional diffusion coefficient (um^2/s) of each status-1
    // attachment domain along its actin barbed-end association zone. The
    // connector force supplies the Einstein drift. Status-2 sites remain
    // fixed material coordinates until molecular detachment.
    double aa_status1_site_diffusion = 0.0;
    // Physical two-arm ACTN2 turnover.  When enabled, material coordinates do
    // not continuously diffuse: annealing occurs when one arm opens while the
    // other remains tethered and the free arm rebinds locally.
    bool aa_arm_state_model = false;
    double aa_arm_open_rate_weak = 20.0;
    double aa_arm_open_rate_mature = 1.0;
    double aa_arm_bell_distance = 0.001;
    double aa_tether_off_rate = 5.0;
    double aa_tether_rebind_rate = 100.0;
    double aa_tether_capture_distance = 0.10;
    double aa_arm_rebind_theta = 0.5;
    // Finite ACTN2-connector/thick-filament contact.  The excluded distance is
    // myosin.radius + aa_actn2_contact_radius. epsilon is in pN um.
    double aa_actn2_contact_radius = 0.005;
    double aa_myosin_contact_epsilon = 0.00428;
    double aa_myosin_contact_force_cap = 40.0;
    // Optional conservative shoulder outside the physical hard core.  It
    // coarse-grains the crowded myosin-head corona.  Width=0 preserves the
    // original short-range WCA contact and disables association gating.
    double aa_myosin_contact_shell_width = 0.0;
    double aa_myosin_contact_shell_energy = 0.0;
    // Optional axial end-stop between a myosin tip and the local bilateral
    // ACTN2 patch. At most one patch is selected per myosin tip per step.
    bool aa_zdisc_tip_contact = false;
    double aa_zdisc_tip_onset = 0.050;
    double aa_zdisc_tip_stiffness = 1000.0;
    double aa_zdisc_tip_force_max = 20.0;
    // Lateral pressure holding ALL barbed-end pairs near a target spacing.
    // Zero range keeps the term inert so existing runs are unchanged.
    double aa_barbed_lateral_range = 0.0;
    double aa_barbed_lateral_k = 0.0;
    double aa_barbed_lateral_force_max = 0.0;
    double aa_barbed_lateral_axial_window = 0.25;
    double aa_zdisc_tip_transverse_flat = 0.100;
    double aa_zdisc_tip_transverse_cutoff = 0.130;
    double aa_occlusion_rupture_penetration = 0.010;
    // Optional molecular polarity gate for new local ACTN2 associations.
    // Values >= 1 disable it; lower values require u_i dot u_j below the
    // threshold without prescribing a neighbor graph.
    double aa_association_max_direction_dot = 1.0;
    // A positive contour length causes an overstretched explicit ACTN2 link
    // to rupture and releases its remembered neighbor slot for reannealing.
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
    std::vector<std::array<double, 2>> aa_scratch_contact_arm_loads;
    std::vector<double> aa_scratch_contact_penetration;
    std::vector<std::vector<std::pair<int, ACTN2MyosinContactResult>>>
        aa_scratch_contact_results;
    std::vector<std::uint8_t> aa_scratch_contact_evaluated;
    std::unordered_map<std::uint64_t, double> aa_scratch_pair_tension;
    std::vector<double> aa_scratch_total_arm_load;
    std::vector<double> aa_scratch_productive_load;
    // Per-actin motor-force vectors are invariant across the whole routine
    // (nothing in it touches AM state), so they are evaluated once per actin
    // instead of twice per incident link.
    std::vector<std::array<double, 3>> aa_scratch_motor_force;
    std::vector<std::array<double, 3>> aa_scratch_productive_force;
    std::vector<std::vector<ProductiveMotorContribution>>
        aa_scratch_productive_by_myosin;
    std::vector<std::vector<double>> aa_scratch_productive_source_demand;
    // Cumulative opt-in profiling. Pair-evaluation time is summed over
    // worker threads; local-update time is wall clock.
    double aa_motor_support_build_seconds = 0.0;
    double aa_motor_candidate_seconds = 0.0;
    double aa_motor_local_update_seconds = 0.0;
    std::uint64_t aa_motor_candidate_pairs_evaluated = 0;
    std::uint64_t aa_motor_profile_steps = 0;
    double aa_motor_update_started = 0.0;
    // Wall-clock attribution for the three top-level update_system phases.
    std::uint64_t system_profile_steps = 0;
    double system_neighbor_seconds = 0.0;
    double system_parallel_seconds = 0.0;
    // Per-section timing inside the omp parallel region of update_system().
    // Accumulated PER THREAD with no added barriers: inserting barriers to
    // get a clean team-wide split would change the execution graph, which
    // this code deliberately preserves for reproducibility.  Each thread's
    // entry therefore includes its own work plus any wait at the implicit
    // barrier that ends the section, so the team mean is the honest
    // attribution and the max approximates the critical path.
    enum AaPhase {
    kAaRebuildIndices,   // rebuild_indices (single)
    kAaArmKinetics,   // arm_kinetics (for)
    kAaSiteDiffusion,   // site_diffusion (for)
    kAaMotorForces,   // motor_forces (for)
    kAaArmLoads,   // arm_loads (single)
    kAaLinkLoop,   // link_loop_ordered (for ordered)
    kAaReindex,   // reindex_expire (single)
    kAaSortNeighbors,   // sort_neighbors (for)
    kAaPrefixSum,   // pair_prefix_sum (single)
    kAaScatterPairs,   // scatter_pairs (for)
    kAaEvalProposals,   // evaluate_proposals (for)
    kAaCommitLinks,   // commit_links (single)
    kAaTipContact,      // tip_contact only (single)
    kAaTipAndScatter,   // mechanics scatter (single)
    kAaPhaseCount
};

enum ProfileSection {
        kSecSetToZero = 0, kSecAmBinding, kSecThinActivation, kSecBondLimit,
        kSecConcatConnections, kSecAaChemistry, kSecChemReduce,
        kSecBondMatrix, kSecAmForce, kSecAaStateTitin, kSecVolumeExclusion,
        kSecReduce, kSecTail, kProfileSectionCount
    };
    std::vector<std::array<double, kProfileSectionCount>>
        parallel_section_seconds;
    static const char* profile_section_name(int section);
    double system_finalize_seconds = 0.0;
    double system_total_seconds = 0.0;
    // Association candidates, evaluated in parallel against the pre-commit
    // state and then committed serially in canonical order.
    struct AAAssociationProposal {
        double s_i = 0.0;
        double s_j = 0.0;
        int nearest_tip_myosin = -1;
        int nearest_tip_sign = 0;
        double nearest_tip_gap = std::numeric_limits<double>::infinity();
        double nearest_tip_rho = std::numeric_limits<double>::infinity();
        bool nearest_tip_motor_engaged = false;
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
    // Sparse, local, finite-valency M-band bridge model. k_mm is the axial
    // stiffness of one occupied bridge, rather than an always-on pair field.
    int myomesin_sites_per_filament = 2;
    double myomesin_k_on = 5000.0;
    double myomesin_k_off = 1.0;
    double myomesin_axial_capture = 0.25;
    double myomesin_radial_rest = 0.11;
    double myomesin_radial_capture = 0.04;
    double myomesin_radial_k = 0.0;
    double myomesin_alignment_cos = 0.94;
    double myomesin_bell_distance = 0.0005;
    double myomesin_force_max = 10.0;
    std::uint64_t myomesin_next_id = 1;
    std::vector<MyomesinBridge> myomesin_bridges;
    std::vector<double> myomesin_events;
    double titin_k, titin_rest_length;
    TitinMode titin_mode = TitinMode::legacy_crosslinked;
    int explicit_titin_per_half = 3;
    double explicit_titin_capture_distance = 0.70;
    double explicit_titin_k_on = 5000.0;
    double explicit_titin_k_off = 0.5;
    double explicit_titin_transition_theta = 0.5;
    bool explicit_titin_strain_gated_kinetics = true;
    // Finite local recapture of the same material barbed-end domain after
    // dissociation. The association barrier remains strain-gated.
    double explicit_titin_rebind_memory = 0.0;
    double explicit_titin_rebind_rate_multiplier = 1.0;
    // Capture specificity. The tension-only force law is exactly zero over the
    // whole slack region, so elastic gating cannot distinguish candidates there
    // and every barbed end inside the capture shell receives the full k_on.
    // When this gate is enabled, *forming* a new titin-actin association also
    // requires the individual actin barbed end to already carry at least one
    // live ACTN2 A-A link: local evidence of a nascent Z-disc complex, with no
    // centroid, Z node, or imposed lattice. Status 1 is sufficient; requiring
    // status 2 would make titin depend on mature motor load. The gate applies
    // to formation only. An already bound titin is never released by ACTN2
    // turnover, and the same-partner rebinding memory window bypasses it so the
    // complex can anneal rather than couple brittlely.
    bool explicit_titin_require_live_actn2 = false;
    // Serial-linkage dissociation. Titin does not reach actin on its own: its
    // Z-repeat region is held by ACTN2, which is what binds the filament. The
    // anchor is therefore actin -- ACTN2 -- titin in series, and its composite
    // lifetime is set by whichever element fails first. A single global k_off
    // implicitly assumes support is permanent, which is true of the
    // preassembled minimal cell it was calibrated on (4/4 status-2 pairs) and
    // false of a random assembly, where 97% of ACTN2 links die geometrically
    // with a median age near 6 ms.
    //
    // When this rate is positive it replaces explicit_titin_k_off for a bound
    // titin whose actin currently carries no live ACTN2 link. Nothing is force
    // released: dissociation stays stochastic and strain gated, and the
    // same-partner rebinding memory still lets a titin ride out a momentary
    // gap. Retention then follows ACTN2 maturation on its own, because a
    // status-2 link persists and a drift-by encounter does not.
    //
    // Zero or negative keeps the single-rate behaviour.
    double explicit_titin_k_off_unsupported = 0.0;
    // Whether same-partner recapture is also subject to the ACTN2 requirement.
    //
    // The rebinding memory was specified for a parallel anchor, where titin
    // holds actin directly and identity plus proximity are enough to
    // re-engage. Under the serial anchor there is nothing to rebind *through*
    // while the ACTN2 is absent, and leaving the exemption in place cancels
    // the serial law outright: a release is followed by an immediate
    // recapture, because recapture carries the x100 memory multiplier while
    // release runs at the unsupported rate.
    //
    // With this enabled the memory keeps its rate advantage but not its
    // exemption, so re-engagement is fast exactly when support returns. That
    // is the annealing the memory was introduced for, without the permanence.
    bool explicit_titin_rebind_requires_actn2 = false;
    std::vector<ExplicitTitin> explicit_titins;
    std::vector<int> explicit_titin_owner_by_actin;
    // Rows: step, event code (1 bind, 2 release), titin index, myosin,
    // half, track, actin, distance, elastic energy, force.
    std::vector<double> explicit_titin_events;
    // Cumulative association provenance, saved per frame as
    // /state/titin_association_counts = {supported, unsupported, rebind}.
    // "Supported" means the actin held a live ACTN2 link at the instant the
    // bond formed; "rebind" counts same-partner recaptures inside the memory
    // window, which are exempt from the gate and therefore counted separately
    // rather than being scored as unsupported.
    std::uint64_t explicit_titin_supported_associations = 0;
    std::uint64_t explicit_titin_unsupported_associations = 0;
    std::uint64_t explicit_titin_rebind_associations = 0;
    // Releases that occurred while the actin held no live ACTN2 link, i.e.
    // attributable to loss of the serial anchor rather than to titin's own
    // chemistry. Saved per frame alongside the association counters.
    std::uint64_t explicit_titin_unsupported_releases = 0;
    // Step at which a fixed-geometry chemistry pre-equilibration ended.
    // Stored as provenance; production frames retain the absolute KMC step.
    size_t release_step_origin = 0;
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
    // A stereospecific status-2 ACTN2 complex retains its material attachment
    // coordinates until detachment, contour rupture, or an enabled reversible
    // mature-to-weak transition. False preserves the historical frame-by-frame
    // classifier when the unlock rate is also zero.
    bool aa_status2_locked_until_detach = false;
    // Optional reversible mature-to-weak conformational transition.  A zero
    // rate preserves historical behavior.  Positive support suppresses the
    // unloaded rate; registry error outside the retention band accelerates it.
    double aa_status2_unlock_rate = 0.0;
    double aa_status2_unlock_force_scale = 0.1;
    double aa_status2_unlock_registry_scale = 0.02;
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
    // One conserved set of coarse material binding sites is shared by every
    // myosin touching an actin. This prevents each A-M pair from independently
    // recreating the same thin-filament binding capacity.
    int am_sites_per_actin = 182;
    // Global simultaneous head occupancy of one actin, independent of how
    // many thick filaments overlap it. This is distinct from the much denser
    // molecular binding-site lattice above.
    int am_max_heads_per_actin = 24;
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
    // Working-stroke treatment.  false = the historical prescribed
    // force-velocity ramp; true = discrete strain-dependent transition
    // obeying detailed balance (force-velocity becomes emergent).
    bool am_stroke_detailed_balance = false;
    double am_stroke_base_rate = 1000.0;
    double am_stroke_theta = 0.5;
    double am_stroke_free_energy = -1.0;  // <0 => derive from am_stall_force

    double am_weak_stiffness_ratio = 0.1;
    // Pair-level flat-bottom containment: force-free near the preferred AM
    // spacing, stiff only outside the capture band so transverse drift cannot
    // masquerade as chemical or axial-overlap dissociation.
    double am_tether_k = 3000.0;
    // Initial misalignment used by single_am_tilted (degrees).
    double probe_tilt_deg = 20.0;
    // Centre-trap stiffness (pN/um); 0 disables the trap entirely.
    double probe_confine_k = 0.0;
    // Anchors, captured lazily at the first force evaluation.
    std::vector<vec> probe_confine_anchor;
    double am_tether_slack = 0.002;
    // Optional passive compliance of each explicit A-M head connector in the
    // plane transverse to the actin axis. Zero is the historical control.
    // Mechanical stiffness and coupling to the chemical escape barrier are
    // independent; am_rupture_distance remains the emergency contour guard.
    double am_head_transverse_k = 0.0;
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
    // Applies only to negative axial strain in a strong head. Weak heads use
    // their chemical off-rate rather than a welded axial-strain penalty.
    double am_strain_detachment_theta = 0.1;
    // Applies to transverse elastic energy in either attached state. The
    // coarse-grained default is zero so the reach spring does not defeat its
    // own retention function.
    double am_transverse_strain_detachment_theta = 0.0;
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
    bool am_q_initialized_at_calcium_equilibrium = false;
    double am_initial_q = 0.0;
    std::vector<double> am_thin_activation;
    // Fixed-population whole-filament renewal.  This is an accelerated
    // topology-remodeling surrogate: a selected actin is detached from every
    // molecular partner and reinserted without changing N or filament length.
    bool actin_turnover_enabled = false;
    double actin_turnover_rate = 0.0;
    double actin_turnover_start_time = 0.0;
    int actin_turnover_max_placement_attempts = 1000;
    // Fraction of the nominal hard-core contact distance a renewed actin
    // must clear.  The steric potential is soft and force-capped, so the
    // equilibrium population itself sits inside the nominal contact
    // distance: at A800/M160 the median incumbent actin-actin separation is
    // 0.89 of 2*r_actin and only ~9% of incumbents clear the full distance.
    // Demanding ratio 1 therefore rejects essentially every proposal and
    // silently disables turnover.  The default admits candidates that are no
    // more overlapped than the incumbent 10th-25th percentile.
    double actin_turnover_clearance_ratio = 0.8;
    std::vector<RecoveryStep> actin_last_turnover_step;
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
    // Short-range actin-actin candidates for _volume_exclusion.  The shared
    // broadphase list is indexed on filament CENTRES, so its cutoff is set by
    // filament length (1.6 um) rather than by any interaction range: each
    // actin carries ~857 candidates of which ~4 are ever in range.  This is a
    // segment-distance-filtered subset, rebuilt only when the cell list is,
    // and used ONLY by the sterics loop whose interactions are all short
    // range.  Every other consumer of actin_neighbors_by_species keeps the
    // full list, because their ranges differ and a shared narrow list would
    // silently drop pairs.
    std::vector<std::vector<int>> actin_close_neighbors;
    // actin.direction is not guaranteed unit length, and three steric kernels
    // normalize it per pair.  Cache the normalized copy once per actin per step
    // so the pair gates can reproduce the kernels' own barbed ends bit for bit
    // without paying 2 normalizations on every one of the ~431k pairs.
    std::vector<vec> actin_unit_direction;
    // Visit stamps for the Z-disc tip-contact candidate scan.
    std::vector<long long> aa_tip_candidate_stamp;
    // Per-link survivors from the parallel tip-contact scan, plus one stamp
    // array per thread so the per-link myosin dedup needs no synchronisation.
    struct AaTipScanHit {
        size_t slot;
        int myosin_index;
        int end_sign;
        ZDiscMyosinTipContactResult contact;
    };
    std::vector<std::vector<AaTipScanHit>> aa_tip_scan_buckets;
    std::vector<std::vector<long long>> aa_tip_stamp_pool;
    // Per-thread phase timers for _update_local_aa_links.  Mean over
    // threads charges each thread its share of the barrier it waits on,
    // which is exactly what makes a serial region visible.
    std::vector<std::array<double, kAaPhaseCount>> aa_phase_seconds;
    static const char* aa_phase_name(int phase);
    double actin_close_cutoff = 0.0;
    bool actin_close_valid = false;
    vector actin_basic_tension;
    std::vector<bool> actin_f_load_computed;
    std::vector<double>* actin_f_load_cb;
    std::vector<std::mutex> actin_f_load_mutex;
        
    gsl_rng* rng;
    std::string filename;
    int initial_seed = 0;
    // Optional provenance for a geometry-only initialization. Molecular
    // occupancies are deliberately not imported from the source trajectory.
    std::string initial_geometry_source;
    int initial_geometry_source_frame_requested = -1;
    int initial_geometry_source_frame_loaded = -1;

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
        bool productive = true;
        double axial_overlap = 0.0;
        double binding_ratio = 0.0;
        double partial_ratio = 0.0;
    };
    std::vector<std::vector<AMHeadCandidate>> am_head_candidates_;
    // Per-step shared material-site claims. Existing attachments populate the
    // table before the parallel myosin KMC loop; successful new associations
    // claim a free site atomically so two myosins cannot bind the same actin
    // monomer during one timestep.
    std::vector<int> am_site_claims_;
    std::vector<int> am_site_claim_count_;
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
    // Width 20: actin,step,age_steps,success,attempts,released_am,
    // released_actn2,released_titin,old_center_xyz,new_center_xyz,
    // old_direction_xyz,new_direction_xyz.
    std::vector<double> actin_turnover_events;
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
    // Versioned audit rows retain the 21 v2 fields and append exact geometry,
    // component gate decisions, and geometric-versus-stochastic detachment
    // provenance. The schema is recorded in /meta.
    std::vector<double> aa_local_link_events_v3;
    // Exact accepted association/rebinding accessibility audit. Width 13:
    // link_id,step,event_code,actin_i,actin_j,s_i,s_j,nearest_myosin,
    // tip_sign,axial_gap,transverse_distance,tip_zone_occupied,
    // motor_engaged. Event codes: 1 association, 2 tether rebinding.
    std::vector<double> aa_tip_accessibility_events;
    // Exact identity-state events for motor_transmitted_register. Width 20;
    // schema is recorded in HDF5 metadata.
    std::vector<double> aa_motor_pair_events;
    // Rows: link_id,step,event_code,old_state,new_state,bound_actin,
    // old_partner,new_partner,s_bound,s_free,arm_force,contact_force,
    // contact_penetration,rate. Event codes are recorded in metadata.
    std::vector<double> aa_arm_events;
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
    void sarcomeric_structure();
    void sarcomeric_structure_tight();
    void set_bundling_parameters(double max_strength, int ramp_steps);
    void set_interaction_ramp(size_t start_step, size_t ramp_steps);
    void set_steric_stiffness(double actin_stiffness,
                              double myosin_stiffness,
                              double myosin_exclusion = 0.0);
    void configure_double_overlap_contact(bool enabled,
                                          double contact_diameter,
                                          double epsilon,
                                          double force_cap = 0.0);
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
    void set_am_site_count(int sites_per_actin);
    void set_am_global_head_cap(int maximum_heads);
    void set_am_weak_binding_floor(double basal_fraction);
    void set_am_proximity_tether(double stiffness, double slack,
                                 bool select_binding_azimuth,
                                 double azimuth_half_width_deg,
                                 double alignment_stiffness = 0.0,
                                 double alignment_slack_deg = 2.0);
    void set_am_head_transverse_tether(double stiffness);
    void set_am_retention_kinetics(bool strain_dependent,
                                   double axial_barrier_fraction,
                                   double transverse_barrier_fraction,
                                   double maximum_exponent,
                                   int tether_saturation_heads,
                                   double strong_off_floor_fraction = 1.0);
    void set_am_thin_activation(bool enabled, bool feedback_enabled,
                                double calcium_open_rate, double calcium_half,
                                double calcium_hill, double strong_feedback_rate,
                                double strong_half, double strong_hill,
                                double close_rate);
    void initialize_am_thin_activation_calcium_equilibrium();
    void set_am_thick_regulation(bool enabled, double recruit_rate,
                                 double fold_rate, double tension_scale,
                                 double initial_on_fraction,
                                 double force_clamp);
    void configure_actin_turnover(bool enabled, double rate,
                                  double start_time,
                                  int max_placement_attempts,
                                  double clearance_ratio);
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
    // Import only filament centers and directions from an HDF5 trajectory.
    // A negative frame selects the latest committed source frame. Returns the
    // resolved zero-based source frame index.
    int load_initial_geometry(const std::string& source_filename,
                              int frame_index = -1);
    void reset_release_diagnostics();
    void configure_explicit_titin(int links_per_half,
                                  double capture_distance,
                                  double association_rate,
                                  double dissociation_rate,
                                  bool strain_gated,
                                  double transition_theta,
                                  double rebind_memory = 0.0,
                                  double rebind_rate_multiplier = 1.0,
                                  bool require_live_actn2 = false,
                                  double k_off_unsupported = 0.0,
                                  bool rebind_requires_actn2 = false);
    void initialize_explicit_titins();
    void configure_myomesin_bridges(int sites_per_filament,
                                    double association_rate,
                                    double dissociation_rate,
                                    double axial_capture,
                                    double radial_rest,
                                    double radial_capture,
                                    double radial_stiffness,
                                    double alignment_cos,
                                    double bell_distance,
                                    double force_max);
    bool prebind_myomesin_bridge(int first, int second);
    int myomesin_degree(int myosin_index) const;
    // True when actin_index carries at least one live ACTN2 A-A link of any
    // status. Public so that occupancy diagnostics and tests can score titin
    // support with exactly the predicate the capture gate uses.
    bool actin_has_live_actn2(int actin_index) const;
    vec explicit_titin_thick_anchor(const ExplicitTitin& titin) const;
    void set_aa_status_kinetics(actn2_load_kinetics::StatusMode mode,
                                double force_on, double force_off);
    void update_system();
    void update_system_sterics_only();
    void set_periodicity(const std::array<bool,3>& periodic_axes);
    void new_file();
    void save_state();
    void save_resume_snapshot();
    int load_state(int& n_frames, int frame_index = -1);
    double aa_link_register_ratio(BondFlag status) const {
        return (aa_register_all_statuses || status == 2)
            ? aa_axial_register_ratio : 0.0;
    }
    double aa_link_perpendicular_stiffness(BondFlag status) const {
        return (aa_perpendicular_all_statuses || status == 2)
            ? aa_perpendicular_k : 0.0;
    }
    double aa_link_radial_stiffness(BondFlag status) const {
        if (status == 1 && aa_status1_radial_k > 0.0) {
            return aa_status1_radial_k;
        }
        return actn2_load_kinetics::fixed_single_link_stiffness(
            k_aa, aa_stiffness_reference_links);
    }
    struct AATipExclusionContact {
        int myosin = -1;
        int end_sign = 0;
        double axial_gap = std::numeric_limits<double>::infinity();
        double transverse_distance = std::numeric_limits<double>::infinity();
        double penetration = 0.0;
        bool motor_engaged = false;
    };
    struct AAMyosinExclusionEvaluation {
        double energy = 0.0;
        bool hard_core_overlap = false;
        // The bilateral bridge midpoint lies inside the axial Z-patch corona
        // of a nearby myosin tip, independently of motor-head occupancy.
        bool axial_tip_contact = false;
        // Maximum g0-g among eligible contacts. Diffusive weak sites already
        // in the contact zone may move only toward smaller penetration.
        double axial_tip_penetration = 0.0;
        int nearest_tip_myosin = -1;
        int nearest_tip_sign = 0;
        double nearest_tip_gap = std::numeric_limits<double>::infinity();
        double nearest_tip_rho = std::numeric_limits<double>::infinity();
        bool nearest_tip_motor_engaged = false;
        std::vector<AATipExclusionContact> axial_tip_contacts;
    };
    // Read-only diagnostic used by focused tests and trajectory probes. It
    // applies exactly the same local exclusion logic as association and weak
    // site diffusion without mutating molecular state.
    AAMyosinExclusionEvaluation evaluate_local_aa_myosin_exclusion(
        int actin_i, int actin_j, double s_i, double s_j) const;
    // Pure no-flux rule used by status-1 contour diffusion. Contacts are
    // identified by the physical (myosin, tip) pair, not just by depth.
    static bool aa_tip_diffusion_move_allowed(
        const AAMyosinExclusionEvaluation& old_evaluation,
        const AAMyosinExclusionEvaluation& proposed_evaluation,
        double tolerance = 1e-12);
    // Debug helper to compute catch-bond statistics for a single frame
    void debug_cb_stats();

private:
    double _current_interaction_ramp() const;
    // Private helper methods
    void _update_neighbors();
    void _refresh_actin_close_neighbors();
    void _process_actin_turnover();
    bool _renew_actin(int actin_index, int& placement_attempts,
                      int& released_am, int& released_actn2,
                      int& released_titin, vec& new_center,
                      vec& new_direction);
    bool _actin_renewal_position_is_clear(int actin_index,
                                          const vec& center,
                                          const vec& direction) const;
    void _set_to_zero();
    void _finalize_active_pair_lists();
    void _rebuild_active_pair_lists_from_dense();
    void _process_actin_myosin_binding(int& i);
    void _process_dynamic_am_heads();
    void _update_thin_filament_activation();
    void _drop_unseen_am_bonds();
    void _initialize_am_head_pools();
    int _am_site_index(double actin_coordinate) const;
    double _am_site_coordinate(int site_index) const;
    double _nearest_free_am_site_coordinate(
        int actin_index, int myosin_index, const AMHeadState& head,
        double binding_azimuth) const;
    double _claim_free_am_site_coordinate(
        int actin_index, int myosin_index, const AMHeadState& head,
        double binding_azimuth);
    void _release_am_site_claim(int actin_index, double actin_coordinate);
    void _rebuild_am_site_claims();
    bool _am_head_is_productive(int actin_index, int myosin_index,
                                const AMHeadState& head) const;
    // True when the actin overlaps the head-bearing band of `half`
    // (-1 or +1), i.e. it has reached past the bare zone into that half.
    bool _am_actin_reaches_half(int actin_index, int myosin_index,
                                int half) const;
    // Productivity for a half rather than a specific head, so candidate
    // generation and the per-head test cannot drift apart.
    bool _am_half_is_productive(int actin_index, int myosin_index,
                                int half) const;
    void _enforce_global_am_site_ownership();
    double _am_head_association_gate(int actin_index, int myosin_index,
                                     const AMHeadState& head,
                                     double binding_azimuth) const;
    bool _am_head_retention_reachable(int actin_index, int myosin_index,
                                      const AMHeadState& head) const;
    double _am_head_strain_energy(int actin_index, int myosin_index,
                                  int head_index,
                                  am_motor_kinetics::State state) const;
    std::array<double, 2> _am_head_detachment_energies(
        int actin_index, int myosin_index, int head_index,
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
    std::vector<ProductiveMotorContribution>
    _actin_productive_motor_contributions(int actin_index) const;
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
    std::vector<std::pair<double, double>>
    _local_aa_tether_capture_intervals(
        int bound_actin, int candidate_actin, double bound_coordinate) const;
    bool _local_aa_neighbor_pair_available(int first, int second) const;
    AAMyosinExclusionEvaluation _local_aa_myosin_exclusion(
        int actin_i, int actin_j, double s_i, double s_j) const;
    AAAssociationProposal _evaluate_aa_association(int i, int j) const;
    bool _local_aa_site_outside_am_overlap(int actin_index,
                                           double axial_coordinate) const;
    void _rebuild_local_aa_indices();
    void _record_aa_state_changes();
    void _calc_am_force_velocity(int& i);
    void _apply_titin_forces(int& i);
    void _update_explicit_titin();
    int _best_explicit_titin_candidate(size_t titin_index,
                                       const vec& thick_anchor) const;
    void _update_myomesin_bridges();
    void _apply_myomesin_bridge(MyomesinBridge& bridge,
                                std::vector<vec>& local_myosin_forces);
    bool _myomesin_pair_geometry(int i, int j, vec& axis,
                                 double& axial, vec& radial,
                                 double& radial_distance) const;
    bool _barbed_ends_within_reach(int i, int j) const;
    bool _aa_pair_may_associate(int i, int j) const;
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
    // Harmonic trap on each actin's centre toward a fixed anchor.  Applied
    // at the centre => zero torque, so orientation stays free.
    void _apply_probe_confinement();
    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
        _extract_bonded_pairs(
        const std::vector<std::vector<BondFlag>>& actin_actin_bonds,
        const std::vector<std::vector<BondFlag>>& actin_actin_status,
        const utils::MoleculeConnection& myosinIndicesPerActin);
    void _enforce_actin_cb_limit();
    void _enforce_myosin_bond_limit();
};

#endif // SARCOMERE_H
