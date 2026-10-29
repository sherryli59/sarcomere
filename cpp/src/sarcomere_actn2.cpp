// ACTN2 subsystem, split out of sarcomere.cpp on 2026-09-17.
//
// One contiguous 2,470-line run (former lines 3770-6239), moved
// byte-for-byte with no logic change.  _update_local_aa_links alone was
// 1,761 lines of it.

#include "sarcomere.h"
#include "aa_refractory.h"
#include "aa_association_kinetics.h"
#include "catch_bond_kinetics.h"
#include "actn2_maturation_geometry.h"
#include "counter_rng.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <thread>
#include <tuple>
#include <unordered_set>
#include <gsl/gsl_randist.h>

#include "sarcomere_internal_gates.h"

namespace {
// Used only by this unit (former sarcomere.cpp:21, referenced at 4340
// and 4526, both inside the moved block).

double boltzmann_exclusion_factor(double beta, double energy) {
    if (!(energy > 0.0) || !std::isfinite(energy)) {
        return std::isfinite(energy) ? 1.0 : 0.0;
    }
    const double exponent = std::clamp(beta * energy, 0.0, 50.0);
    return std::exp(-exponent);
}

}  // namespace

std::uint64_t Sarcomere::_aa_pair_key(int i, int j) const {
    const auto first = static_cast<std::uint32_t>(std::min(i, j));
    const auto second = static_cast<std::uint32_t>(std::max(i, j));
    return (static_cast<std::uint64_t>(first) << 32) |
           static_cast<std::uint64_t>(second);
}

void Sarcomere::_rebuild_local_aa_indices() {
    aa_active_link_slots.clear();
    aa_links_by_pair.clear();
    // assign() would destroy and reallocate every per-actin bucket twice per
    // step.  Clearing in place leaves the same empty buckets but keeps their
    // capacity, which is what the incremental add_link() path relies on too.
    if (aa_links_by_actin.size() != static_cast<size_t>(actin.n)) {
        aa_links_by_actin.assign(static_cast<size_t>(actin.n), {});
    } else {
        for (auto& bucket : aa_links_by_actin) bucket.clear();
    }
    for (size_t slot = 0; slot < aa_links.size(); ++slot) {
        const auto& link = aa_links[slot];
        if (!link.active) continue;
        aa_active_link_slots.push_back(slot);
        const bool arm_bridge = !aa_arm_state_model ||
            actn2_arm_kinetics::is_bridge(link.arm_state);
        if (arm_bridge) {
            aa_links_by_pair[_aa_pair_key(link.actin_i, link.actin_j)].push_back(slot);
            aa_links_by_actin[static_cast<size_t>(link.actin_i)].push_back(slot);
            aa_links_by_actin[static_cast<size_t>(link.actin_j)].push_back(slot);
        } else if (link.arm_state ==
                   actn2_arm_kinetics::State::tethered_to_i) {
            aa_links_by_actin[static_cast<size_t>(link.actin_i)].push_back(slot);
        } else if (link.arm_state ==
                   actn2_arm_kinetics::State::tethered_to_j) {
            aa_links_by_actin[static_cast<size_t>(link.actin_j)].push_back(slot);
        }
    }
}

void Sarcomere::set_aa_neighbor_identity(
    bool enabled, double status1_memory, double status2_memory) {
    aa_neighbor_identity_enabled = enabled;
    aa_neighbor_slot_memory_steps = status1_memory > 0.0
        ? static_cast<size_t>(std::ceil(status1_memory / dt)) : 0;
    aa_status2_slot_memory_steps = status2_memory > 0.0
        ? static_cast<size_t>(std::ceil(status2_memory / dt)) : 0;
    if (!enabled) aa_neighbor_slots.clear();
}

bool Sarcomere::reserve_aa_neighbor_pair(
    int first, int second, bool pinned) {
    if (!aa_neighbor_identity_enabled) return true;
    return aa_neighbor_slots.reserve(
        first, second, static_cast<aa_neighbor_identity::Step>(current_step),
        max_aa_bonds_per_actin, pinned);
}

int Sarcomere::aa_neighbor_degree(int actin_index) const {
    return aa_neighbor_identity_enabled
        ? aa_neighbor_slots.degree(actin_index) : 0;
}

// The axial frame of a filament pair does not depend on which ACTN2 site of
// that pair is being probed, so it is factored out of the ordinal loop in
// _update_local_aa_links().  Every expression below is copied verbatim from
// the original _local_aa_attachment_candidate() so the split is bitwise
// neutral.
Sarcomere::AAPairAxis Sarcomere::_local_aa_pair_axis(int i, int j) const {
    AAPairAxis geometry;
    vec ui = actin.direction[i];
    vec uj = actin.direction[j];
    ui.normalize();
    uj.normalize();
    if (ui.dot(uj) < 0.0) uj = -uj;
    vec axis = ui + uj;
    if (axis.norm() <= 1e-10) axis = ui;
    axis.normalize();
    const double pi = actin.center[i].dot(axis);
    vec center_displacement = actin.center[j] - actin.center[i];
    center_displacement.pbc_wrap(box, is_periodic);
    const double pj = pi + center_displacement.dot(axis);
    const double di = actin.direction[i].dot(axis);
    const double dj = actin.direction[j].dot(axis);
    const double hi = 0.5 * actin.length * std::abs(di);
    const double hj = 0.5 * actin.length * std::abs(dj);
    double lower = std::max(pi - hi, pj - hj);
    double upper = std::min(pi + hi, pj + hj);
    // ACTN2 is confined to the Z-disc: a zone of width aa_barbed_zone measured
    // from each filament's BARBED end, which is filament coordinate s = -L/2.
    // A crosslink can only form where BOTH filaments' zones overlap, so the
    // mutual overlap is intersected with each.  The surviving length is the
    // accessible binding length and scales the association rate through
    // overlap_extensive_factor; side-by-side pairs, whose zones never meet,
    // are left with zero accessible length and cannot bind at all.
    if (aa_barbed_zone > 0.0) {
        const double half = 0.5 * actin.length;
        const double ai = pi - half * di;
        const double bi = pi + (aa_barbed_zone - half) * di;
        const double aj = pj - half * dj;
        const double bj = pj + (aa_barbed_zone - half) * dj;
        lower = std::max(lower, std::max(std::min(ai, bi), std::min(aj, bj)));
        upper = std::min(upper, std::min(std::max(ai, bi), std::max(aj, bj)));
    }
    geometry.axis = axis;
    geometry.pi = pi;
    geometry.pj = pj;
    geometry.lower = lower;
    geometry.overlap = std::max(0.0, upper - lower);
    geometry.denominator_i = di;
    geometry.denominator_j = dj;
    return geometry;
}

std::pair<double, double> Sarcomere::_local_aa_attachment_site(
    const AAPairAxis& geometry, int ordinal, double& overlap) const {
    overlap = geometry.overlap;
    if (overlap <= 0.0) return {0.0, 0.0};
    // Place every coarse-grained ACTN2 molecule at its own material site.
    // The former two-slot expression repeated 0.35/0.65 for ordinal >= 2,
    // so increasing the molecular redundancy created duplicate springs rather
    // than additional binding sites.  Even spacing gives N distinct sites and
    // reduces to the same central pair (1/3, 2/3) when N == 2.
    const int site_count = std::max(1, aa_prebundle_links_per_pair);
    const int site_index = std::clamp(ordinal, 0, site_count - 1);
    const double fraction = static_cast<double>(site_index + 1) /
                            static_cast<double>(site_count + 1);
    const double target = geometry.lower + fraction * overlap;
    const double denominator_i = geometry.denominator_i;
    const double denominator_j = geometry.denominator_j;
    if (std::abs(denominator_i) <= 1e-8 || std::abs(denominator_j) <= 1e-8) {
        overlap = 0.0;
        return {0.0, 0.0};
    }
    const double s_i = std::clamp(
        (target - geometry.pi) / denominator_i,
        -0.5 * actin.length, 0.5 * actin.length);
    const double s_j = std::clamp(
        (target - geometry.pj) / denominator_j,
        -0.5 * actin.length, 0.5 * actin.length);
    return {s_i, s_j};
}

std::pair<double, double> Sarcomere::_local_aa_attachment_candidate(
    int i, int j, int ordinal, double& overlap) const {
    return _local_aa_attachment_site(
        _local_aa_pair_axis(i, j), ordinal, overlap);
}

bool Sarcomere::_local_aa_site_outside_am_overlap(
    int actin_index, double axial_coordinate) const {
    if (actin_index < 0 || actin_index >= actin.n) return false;
    // Legacy directional geometry defined the crosslinkable segment from the
    // barbed end (left_end, s=-L/2) up to the beginning of AM overlap.  With no
    // engaged myosin the ratio is one and the full actin remains available,
    // preserving prebundling before thick-filament incorporation.
    const double ratio = std::clamp(
        actin_crosslink_ratio[actin_index], 0.0, 1.0);
    if (ratio <= EPS) return false;
    const double allowed_end =
        -0.5 * actin.length + ratio * actin.length;
    return axial_coordinate <= allowed_end + 1e-12;
}

// Local evidence that this individual barbed end is part of a nascent Z-disc
// complex: at least one live ACTN2 molecule of any status holds it.  Status 1
// counts, because status 2 is gated on mature motor load and would recreate the
// original low-capture failure.  The per-actin index is rebuilt at the end of
// _update_local_aa_links(), which runs earlier in the same step than
// _update_explicit_titin(), so the buckets are current here.  A size mismatch
// means the local ACTN2 model is not active at all; that is reported as "no
// support" so the gate cannot silently pass on a stale index.
bool Sarcomere::actin_has_live_actn2(int actin_index) const {
    if (actin_index < 0 || actin_index >= actin.n) return false;
    if (aa_links_by_actin.size() != static_cast<size_t>(actin.n)) return false;
    for (size_t slot : aa_links_by_actin[static_cast<size_t>(actin_index)]) {
        if (slot >= aa_links.size() || !aa_links[slot].active) continue;
        const auto& link = aa_links[slot];
        if (!aa_arm_state_model ||
            (link.actin_i == actin_index &&
             actn2_arm_kinetics::arm_i_bound(link.arm_state)) ||
            (link.actin_j == actin_index &&
             actn2_arm_kinetics::arm_j_bound(link.arm_state))) {
            return true;
        }
    }
    return false;
}

int Sarcomere::_local_aa_distinct_neighbors(int actin_index) const {
    std::vector<int> neighbors;
    for (size_t slot : aa_links_by_actin[static_cast<size_t>(actin_index)]) {
        const auto& link = aa_links[slot];
        if (aa_arm_state_model &&
            !actn2_arm_kinetics::is_bridge(link.arm_state)) continue;
        neighbors.push_back(link.actin_i == actin_index
            ? link.actin_j : link.actin_i);
    }
    std::sort(neighbors.begin(), neighbors.end());
    neighbors.erase(std::unique(neighbors.begin(), neighbors.end()),
                    neighbors.end());
    return static_cast<int>(neighbors.size());
}

bool Sarcomere::_local_aa_site_available(
    int actin_index, int partner, double coordinate) const {
    for (size_t slot : aa_links_by_actin[static_cast<size_t>(actin_index)]) {
        const auto& link = aa_links[slot];
        if (!link.active) continue;
        // Tethers are indexed only on their still-bound actin, so their
        // surviving material site remains occupied while the free arm is
        // available to close.  Treating the entire tether as absent would let
        // a second molecule occupy the same pair-facing site during annealing.
        const int other = link.actin_i == actin_index
            ? link.actin_j : link.actin_i;
        // The reduced model has no ACTN2 azimuth.  Links to distinct
        // neighboring filaments may therefore share an axial level and
        // occupy different azimuthal faces; contour exclusion applies
        // only to duplicate links within the same filament pair.
        if (other != partner) continue;
        const double occupied = link.actin_i == actin_index ? link.s_i : link.s_j;
        if (std::abs(occupied - coordinate) < aa_site_exclusion) return false;
    }
    return true;
}

bool Sarcomere::_local_aa_neighbor_pair_available(
    int first, int second) const {
    return !aa_neighbor_identity_enabled ||
        aa_neighbor_slots.can_reserve(
            first, second, max_aa_bonds_per_actin);
}

std::vector<std::pair<double, double>>
Sarcomere::_local_aa_tether_capture_intervals(
    int bound_actin, int candidate_actin, double bound_coordinate) const {
    std::vector<std::pair<double, double>> intervals;
    if (bound_actin < 0 || bound_actin >= actin.n ||
        candidate_actin < 0 || candidate_actin >= actin.n ||
        bound_actin == candidate_actin) {
        return intervals;
    }

    // The free domain searches the candidate's continuous exposed contour,
    // not the canonical sites used to seed a new two-filament encounter.  Its
    // remaining arm stays at its material coordinate, making this a genuinely
    // local intramolecular capture problem.
    const double contour_lower = -0.5 * actin.length;
    double contour_upper = 0.5 * actin.length;
    if (aa_barbed_zone > 0.0) {
        contour_upper = std::min(
            contour_upper, contour_lower + aa_barbed_zone);
    }
    const double crosslink_ratio = std::clamp(
        actin_crosslink_ratio[candidate_actin], 0.0, 1.0);
    contour_upper = std::min(
        contour_upper, contour_lower + crosslink_ratio * actin.length);
    if (!(contour_upper > contour_lower)) return intervals;

    vec candidate_direction = actin.direction[candidate_actin];
    if (!(candidate_direction.norm() > 1e-12)) return intervals;
    candidate_direction.normalize();
    const vec bound_center = actin.center[bound_actin];
    const vec candidate_center = actin.center[candidate_actin];
    const vec bound_point = bound_center +
        bound_coordinate * actin.direction[bound_actin];
    vec center_from_bound = candidate_center - bound_point;
    center_from_bound.pbc_wrap(box, is_periodic);
    const double closest_coordinate =
        -center_from_bound.dot(candidate_direction);
    const vec transverse = center_from_bound +
        closest_coordinate * candidate_direction;
    const double radial_squared = transverse.dot(transverse);
    const double capture_squared =
        aa_tether_capture_distance * aa_tether_capture_distance;
    if (radial_squared >= capture_squared) return intervals;
    const double axial_reach = std::sqrt(std::max(
        0.0, capture_squared - radial_squared));
    const double reachable_lower = std::max(
        contour_lower, closest_coordinate - axial_reach);
    const double reachable_upper = std::min(
        contour_upper, closest_coordinate + axial_reach);
    if (!(reachable_upper > reachable_lower)) return intervals;
    intervals.push_back({reachable_lower, reachable_upper});

    // Remove pair-facing contour already occupied by another ACTN2 molecule.
    // Distinct neighboring filaments may still bind at the same axial level,
    // matching _local_aa_site_available().
    if (!(aa_site_exclusion > 0.0)) return intervals;
    for (size_t slot :
         aa_links_by_actin[static_cast<size_t>(candidate_actin)]) {
        const auto& link = aa_links[slot];
        if (!link.active) continue;
        const int other = link.actin_i == candidate_actin
            ? link.actin_j : link.actin_i;
        if (other != bound_actin) continue;
        const double occupied = link.actin_i == candidate_actin
            ? link.s_i : link.s_j;
        const double excluded_lower = occupied - aa_site_exclusion;
        const double excluded_upper = occupied + aa_site_exclusion;
        std::vector<std::pair<double, double>> remaining;
        remaining.reserve(intervals.size() + 1);
        for (const auto& interval : intervals) {
            if (excluded_upper <= interval.first ||
                excluded_lower >= interval.second) {
                remaining.push_back(interval);
                continue;
            }
            if (excluded_lower > interval.first) {
                remaining.push_back({
                    interval.first,
                    std::min(interval.second, excluded_lower)});
            }
            if (excluded_upper < interval.second) {
                remaining.push_back({
                    std::max(interval.first, excluded_upper),
                    interval.second});
            }
        }
        intervals.swap(remaining);
        if (intervals.empty()) break;
    }
    return intervals;
}

bool Sarcomere::_aa_status2_connector_excluded(
    int actin_i, int actin_j, double s_i, double s_j,
    double penetration_buffer) const {
    const auto patch = compute_actn2_patch_geometry(
        actin, actin_i, actin_j, s_i, s_j, box);
    // Compute the segment sphere once per link, then reject distant myosins
    // before invoking the exact periodic-image kernel. Keep all myosins:
    // endpoint force-neighbor lists cannot conservatively bound a connector.
    const double reach = std::hypot(
        0.5 * myosin.length - penetration_buffer,
        actn2_maturation::exclusion_radius - penetration_buffer) +
        patch.half_length;
    const double reach_sq = reach * reach;
    for (int m = 0; m < myosin.n; ++m) {
        vec relative = patch.midpoint - myosin.center[m];
        relative.pbc_wrap(box, myosin.periodic_axes);
        if (relative.norm_squared() > reach_sq) continue;
        if (actn2_maturation::connector_overlaps_myosin(
                patch.first, patch.second, myosin.center[m],
                myosin.direction[m], myosin.length, box, myosin.periodic_axes,
                actn2_maturation::exclusion_radius - penetration_buffer,
                penetration_buffer))
            return true;
    }
    return false;
}

Sarcomere::AAMyosinExclusionEvaluation
Sarcomere::_local_aa_myosin_exclusion(
    int actin_i, int actin_j, double s_i, double s_j) const {
    AAMyosinExclusionEvaluation evaluation;
    const bool radial_contact_enabled =
        aa_actn2_contact_radius > 0.0 &&
        aa_myosin_contact_shell_width > 0.0 &&
        aa_myosin_contact_shell_energy > 0.0;
    if (!radial_contact_enabled && !aa_zdisc_tip_contact) {
        return evaluation;
    }
    const auto& first_neighbors =
        actin_neighbors_by_species[static_cast<size_t>(actin_i)].second;
    const auto& second_neighbors =
        actin_neighbors_by_species[static_cast<size_t>(actin_j)].second;
    const double hard_core_distance =
        myosin.radius + aa_actn2_contact_radius;
    const double radial_outer_distance = hard_core_distance +
        (radial_contact_enabled ? aa_myosin_contact_shell_width : 0.0);
    const auto patch = compute_actn2_patch_geometry(
        actin, actin_i, actin_j, s_i, s_j, box);
    auto evaluate_myosin = [&](int myosin_index) {
        const bool may_contact_body = radial_contact_enabled &&
            actn2_patch_may_contact_myosin_body(
                patch, myosin, myosin_index, box, radial_outer_distance);
        const bool may_contact_tip = aa_zdisc_tip_contact &&
            actn2_patch_may_contact_myosin_tip(
                patch, myosin, myosin_index, box, aa_zdisc_tip_onset,
                aa_zdisc_tip_transverse_cutoff);
        if (!may_contact_body && !may_contact_tip) return;
        if (may_contact_body) {
            const auto contact = compute_actn2_myosin_contact(
                actin, myosin, actin_i, actin_j, s_i, s_j, myosin_index,
                box, hard_core_distance, aa_myosin_contact_epsilon,
                aa_myosin_contact_force_cap, aa_myosin_contact_shell_width,
                aa_myosin_contact_shell_energy);
            evaluation.energy += contact.energy;
            evaluation.hard_core_overlap = evaluation.hard_core_overlap ||
                contact.penetration > 0.0;
        }
        // Steric accessibility is a property of the thick-filament body, not
        // of motor occupancy.  Motor engagement is retained only as an audit
        // field; it never decides whether ACTN2 may occupy this patch.
        bool motor_engages_pair = false;
        if (may_contact_tip && myosin_index >= 0 &&
            myosin_index < static_cast<int>(am_heads_by_myosin.size())) {
            for (const auto& half_heads :
                 am_heads_by_myosin[static_cast<size_t>(myosin_index)]) {
                for (const auto& head : half_heads) {
                    if (head.state == static_cast<BondFlag>(
                            am_motor_kinetics::State::detached)) continue;
                    if (head.bound_actin == actin_i ||
                        head.bound_actin == actin_j) {
                        motor_engages_pair = true;
                        break;
                    }
                }
                if (motor_engages_pair) break;
            }
        }
        if (may_contact_tip) {
            for (int end_sign : {-1, 1}) {
                const auto tip = compute_zdisc_myosin_tip_contact(
                    actin, myosin, actin_i, actin_j, s_i, s_j,
                    myosin_index, end_sign, box, aa_zdisc_tip_onset,
                    aa_zdisc_tip_stiffness, aa_zdisc_tip_force_max,
                    aa_zdisc_tip_transverse_flat,
                    aa_zdisc_tip_transverse_cutoff);
                if (tip.eligible_half &&
                    tip.transverse_distance < aa_zdisc_tip_transverse_cutoff &&
                    tip.axial_gap < evaluation.nearest_tip_gap) {
                    evaluation.nearest_tip_myosin = myosin_index;
                    evaluation.nearest_tip_sign = end_sign;
                    evaluation.nearest_tip_gap = tip.axial_gap;
                    evaluation.nearest_tip_rho = tip.transverse_distance;
                    evaluation.nearest_tip_motor_engaged =
                        motor_engages_pair;
                }
                if (!tip.applied) continue;
                evaluation.axial_tip_contact = true;
                evaluation.axial_tip_penetration = std::max(
                    evaluation.axial_tip_penetration, tip.compression);
                evaluation.axial_tip_contacts.push_back({
                    myosin_index, end_sign, tip.axial_gap,
                    tip.transverse_distance, tip.compression,
                    motor_engages_pair});
            }
        }
    };
    for (int myosin_index : first_neighbors) {
        evaluate_myosin(myosin_index);
    }
    for (int myosin_index : second_neighbors) {
        if (std::find(first_neighbors.begin(), first_neighbors.end(),
                      myosin_index) == first_neighbors.end()) {
            evaluate_myosin(myosin_index);
        }
    }
    return evaluation;
}

Sarcomere::AAMyosinExclusionEvaluation
Sarcomere::evaluate_local_aa_myosin_exclusion(
    int actin_i, int actin_j, double s_i, double s_j) const {
    return _local_aa_myosin_exclusion(actin_i, actin_j, s_i, s_j);
}

bool Sarcomere::aa_tip_diffusion_move_allowed(
    const AAMyosinExclusionEvaluation& old_evaluation,
    const AAMyosinExclusionEvaluation& proposed_evaluation,
    double tolerance) {
    tolerance = std::max(0.0, tolerance);
    for (const auto& proposed_contact :
         proposed_evaluation.axial_tip_contacts) {
        const auto old_contact = std::find_if(
            old_evaluation.axial_tip_contacts.begin(),
            old_evaluation.axial_tip_contacts.end(),
            [&](const AATipExclusionContact& contact) {
                return contact.myosin == proposed_contact.myosin &&
                       contact.end_sign == proposed_contact.end_sign;
            });
        if (old_contact == old_evaluation.axial_tip_contacts.end() ||
            proposed_contact.penetration >
                old_contact->penetration + tolerance) {
            return false;
        }
    }
    return true;
}

// Pure evaluation of one association proposal against the current link state.
// It performs no mutation, so it may run concurrently for many pairs; the
// caller commits the surviving proposals serially in canonical order and
// re-evaluates any pair whose actins were touched by an earlier commit.
// Necessary conditions for _evaluate_aa_association to return a firing
// proposal.  Both are geometry-only, so they cannot change as links commit
// during the step -- which is what makes it safe to drop a pair from the
// candidate list rather than re-test it in the serial loop.
//
//   1. it opens with
//        if (dir_i . dir_j >= aa_association_max_direction_dot) return {};
//   2. it ends the geometric section with
//        if (candidate.length >= aa_cutoff) return {};
//      where length is |connector| between two attachment points
//      centre + s*direction with |s| <= L/2, i.e. points ON the two segments.
//      So segment_segment_distance <= length, and a pair whose segments are
//      further apart than aa_cutoff can never satisfy it.
bool Sarcomere::_aa_pair_may_associate(int i, int j) const {
    if (actin.direction[i].dot(actin.direction[j]) >=
        aa_association_max_direction_dot) {
        return false;
    }
    return geometry::segments_may_interact(
        actin.left_end[i], actin.right_end[i],
        actin.left_end[j], actin.right_end[j],
        aa_cutoff, box, is_periodic);
}

Sarcomere::AAAssociationProposal Sarcomere::_evaluate_aa_association(
    int i, int j) const {
    AAAssociationProposal proposal;
    if (actin.direction[i].dot(actin.direction[j]) >=
        aa_association_max_direction_dot) {
        return proposal;
    }
    const auto key = _aa_pair_key(i, j);
    const auto pair_it = aa_links_by_pair.find(key);
    const int occupied = pair_it == aa_links_by_pair.end()
        ? 0
        : static_cast<int>(pair_it->second.size());
    if (occupied >= aa_prebundle_links_per_pair) return proposal;
    if (aa_total_link_pool > 0 &&
        static_cast<int>(aa_active_link_slots.size()) >= aa_total_link_pool) {
        return proposal;
    }
    const bool has_neighbor_slot =
        aa_neighbor_identity_enabled && aa_neighbor_slots.contains(i, j);
    if (max_aa_bonds_per_actin > 0 && occupied == 0 && !has_neighbor_slot) {
        const int degree_i = aa_neighbor_identity_enabled
            ? aa_neighbor_slots.degree(i) : _local_aa_distinct_neighbors(i);
        const int degree_j = aa_neighbor_identity_enabled
            ? aa_neighbor_slots.degree(j) : _local_aa_distinct_neighbors(j);
        if (degree_i >= max_aa_bonds_per_actin ||
            degree_j >= max_aa_bonds_per_actin) return proposal;
    }
    double overlap = 0.0;
    std::pair<double, double> attachment{0.0, 0.0};
    bool found_site = false;
    // Occupancy alone does not identify which of the two axial slots
    // remains after a detachment.  Probe canonical slots explicitly
    // so either missing site can repopulate.  The pair axis is shared by
    // every ordinal, so it is built once.
    const AAPairAxis geometry = _local_aa_pair_axis(i, j);
    for (int ordinal = 0; ordinal < aa_prebundle_links_per_pair; ++ordinal) {
        double candidate_overlap = 0.0;
        const auto candidate_attachment =
            _local_aa_attachment_site(geometry, ordinal, candidate_overlap);
        if (candidate_overlap <= 0.0) continue;
        if (!_local_aa_site_available(i, j, candidate_attachment.first) ||
            !_local_aa_site_available(j, i, candidate_attachment.second)) {
            continue;
        }
        if (!_local_aa_site_outside_am_overlap(
                i, candidate_attachment.first) ||
            !_local_aa_site_outside_am_overlap(
                j, candidate_attachment.second)) {
            continue;
        }
        overlap = candidate_overlap;
        attachment = candidate_attachment;
        found_site = true;
        break;
    }
    if (!found_site) return proposal;
    const auto candidate = compute_local_aa_link(
        actin, i, j, attachment.first, attachment.second, box,
        aa_encounter_k, aa_encounter_kappa, 0.0, aa_optimal);
    if (candidate.length >= aa_cutoff) return proposal;
    const auto exclusion = _local_aa_myosin_exclusion(
        i, j, attachment.first, attachment.second);
    proposal.nearest_tip_myosin = exclusion.nearest_tip_myosin;
    proposal.nearest_tip_sign = exclusion.nearest_tip_sign;
    proposal.nearest_tip_gap = exclusion.nearest_tip_gap;
    proposal.nearest_tip_rho = exclusion.nearest_tip_rho;
    proposal.nearest_tip_motor_engaged =
        exclusion.nearest_tip_motor_engaged;
    // A molecular bridge cannot be born through the thick-filament hard core.
    // Configurations in the soft corona remain possible but pay their actual
    // conservative Boltzmann weight below.
    if (exclusion.hard_core_overlap || exclusion.axial_tip_contact) {
        return proposal;
    }
    const double overlap_factor =
        aa_association_kinetics::overlap_extensive_factor(
            overlap, aa_overlap_capture);
    double association_rate =
        (aa_prebundle_links_per_pair - occupied) * k_on * overlap_factor;
    if (aa_total_link_pool > 0) {
        association_rate *= std::max(0.0,
            1.0 - static_cast<double>(aa_active_link_slots.size()) /
                      static_cast<double>(aa_total_link_pool));
    }
    if (aa_energy_gated_association || aa_passive_detailed_balance) {
        association_rate = aa_association_kinetics::transition_state_on_rate(
            association_rate, thermal_beta, candidate.energy,
            aa_transition_state_theta);
    }
    association_rate *= boltzmann_exclusion_factor(
        thermal_beta, exclusion.energy);
    counter_rng::CounterRng association_rng(
        static_cast<std::uint64_t>(initial_seed),
        static_cast<std::uint64_t>(current_step),
        counter_rng::Stream::aa_pair_association,
        counter_rng::canonical_pair_key(i, j));
    proposal.s_i = attachment.first;
    proposal.s_j = attachment.second;
    proposal.fire = association_rng.uniform() <
        aa_association_kinetics::event_probability(association_rate, dt);
    return proposal;
}

const char* Sarcomere::aa_phase_name(int phase) {
    switch (phase) {
        case kAaRebuildIndices: return "rebuild_indices";
        case kAaArmKinetics: return "arm_kinetics";
        case kAaSiteDiffusion: return "site_diffusion";
        case kAaMotorForces: return "motor_forces";
        case kAaArmLoads: return "arm_loads";
        case kAaLinkLoop: return "link_loop_ordered";
        case kAaReindex: return "reindex_expire";
        case kAaSortNeighbors: return "sort_neighbors";
        case kAaPrefixSum: return "pair_prefix_sum";
        case kAaScatterPairs: return "scatter_pairs";
        case kAaEvalProposals: return "evaluate_proposals";
        case kAaCommitLinks: return "commit_links";
        case kAaTipContact: return "tip_contact";
        case kAaTipAndScatter: return "mechanics_scatter";
        default: return "unknown";
    }
}

void Sarcomere::_update_local_aa_links() {
    auto& phase_seconds = aa_phase_seconds[omp_get_thread_num()];
    auto phase_mark = std::chrono::steady_clock::now();
    int phase_index = 0;
    const auto stamp_aa_phase = [&](int next) {
        const auto now = std::chrono::steady_clock::now();
        phase_seconds[phase_index] +=
            std::chrono::duration<double>(now - phase_mark).count();
        phase_mark = now;
        phase_index = next;
    };
    // Called from inside the `omp parallel` region of update_system() by every
    // thread of the team.  The routine is a sequence of orphaned worksharing
    // constructs separated by `single` blocks that own the shared mutations
    // (index rebuilds, link creation/destruction, event emission).  All
    // cross-phase state therefore lives in aa_scratch_* members rather than
    // function locals, which would otherwise be private per thread.
    const double ramp = _current_interaction_ramp();
    const double single_link_alignment_k = ramp *
        actn2_load_kinetics::fixed_single_link_stiffness(
            kappa_aa, aa_stiffness_reference_links);
    const bool motor_transmitted_mode =
        aa_status_mode ==
        actn2_load_kinetics::StatusMode::motor_transmitted_register;

    stamp_aa_phase(kAaRebuildIndices);
    #pragma omp single
    {
    if (motor_transmitted_mode) aa_motor_update_started = omp_get_wtime();
    _rebuild_local_aa_indices();
    // Existing single-arm tethers turn over before bridge mechanics.  A free
    // arm closes by a local KMC event; the bound material coordinate never
    // slides or teleports.  Bridges that open later in this step first become
    // eligible to rebind on the next step, giving every molecular transition a
    // well-defined order.
    if (aa_arm_state_model) {
        const auto tether_slots = aa_active_link_slots;
        for (size_t slot : tether_slots) {
            auto& link = aa_links[slot];
            if (!link.active ||
                !actn2_arm_kinetics::is_tether(link.arm_state)) continue;
            const auto old_state = link.arm_state;
            const bool bound_to_i = old_state ==
                actn2_arm_kinetics::State::tethered_to_i;
            const int bound_actin = bound_to_i ? link.actin_i : link.actin_j;
            const int old_partner = bound_to_i ? link.actin_j : link.actin_i;
            const double bound_s = bound_to_i ? link.s_i : link.s_j;
            counter_rng::CounterRng detach_rng(
                static_cast<std::uint64_t>(initial_seed),
                static_cast<std::uint64_t>(current_step),
                counter_rng::Stream::aa_tether_detachment, link.id);
            const double off_probability =
                actn2_arm_kinetics::event_probability(aa_tether_off_rate, dt);
            const bool stochastic_detachment =
                detach_rng.uniform() < off_probability;
            const bool no_crosslinkable_contour =
                actn2_arm_kinetics::contour_fully_occluded(
                    actin_crosslink_ratio[bound_actin], EPS);
            if (stochastic_detachment || no_crosslinkable_contour) {
                aa_arm_events.insert(aa_arm_events.end(), {
                    static_cast<double>(link.id),
                    static_cast<double>(current_step), 3.0,
                    static_cast<double>(static_cast<int>(old_state)), 0.0,
                    static_cast<double>(bound_actin),
                    static_cast<double>(old_partner), -1.0, bound_s,
                    std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0, 0.0,
                    no_crosslinkable_contour
                        ? std::numeric_limits<double>::infinity()
                        : aa_tether_off_rate});
                link.active = false;
                aa_free_link_slots.push_back(slot);
                continue;
            }

            struct RebindCandidate {
                int actin = -1;
                double lower_s = 0.0;
                double upper_s = 0.0;
                double rate = 0.0;
                double energy = 0.0;
            };
            std::vector<RebindCandidate> candidates;
            double total_rate = 0.0;
            for (int candidate_actin :
                 actin_neighbors_by_species[static_cast<size_t>(bound_actin)].first) {
                if (candidate_actin == bound_actin) continue;
                if (actin.direction[bound_actin].dot(
                        actin.direction[candidate_actin]) >=
                    aa_association_max_direction_dot) continue;
                if (!_local_aa_neighbor_pair_available(
                        bound_actin, candidate_actin)) continue;
                const auto pair_it = aa_links_by_pair.find(
                    _aa_pair_key(bound_actin, candidate_actin));
                if (pair_it != aa_links_by_pair.end() &&
                    static_cast<int>(pair_it->second.size()) >=
                        aa_prebundle_links_per_pair) continue;
                const auto capture_intervals =
                    _local_aa_tether_capture_intervals(
                        bound_actin, candidate_actin, bound_s);
                for (const auto& interval : capture_intervals) {
                    const double interval_width =
                        interval.second - interval.first;
                    if (!(interval_width > 0.0)) continue;
                    // Midpoint quadrature resolves the continuous contour at
                    // no coarser than 10 nm (and at half the site-exclusion
                    // scale when that is smaller).  Each bin carries its
                    // physical contour measure, so refining the quadrature
                    // cannot multiply the molecular attempt rate.
                    const double bin_target = std::max(
                        1e-6, std::min(
                            0.010, aa_site_exclusion > 0.0
                                ? 0.5 * aa_site_exclusion : 0.010));
                    const int bin_count = std::max(
                        1, static_cast<int>(std::ceil(
                            interval_width / bin_target)));
                    const double bin_width = interval_width / bin_count;
                    for (int bin = 0; bin < bin_count; ++bin) {
                    const double lower_s = interval.first + bin * bin_width;
                    const double upper_s = lower_s + bin_width;
                    const double free_s = 0.5 * (lower_s + upper_s);
                    const auto candidate_mechanics = compute_local_aa_link(
                        actin,
                        bound_to_i ? bound_actin : candidate_actin,
                        bound_to_i ? candidate_actin : bound_actin,
                        bound_to_i ? bound_s : free_s,
                        bound_to_i ? free_s : bound_s,
                        box, aa_link_radial_stiffness(static_cast<BondFlag>(1)),
                        aa_link_perpendicular_stiffness(static_cast<BondFlag>(1)),
                        0.0, aa_optimal, 0.0, 0.0,
                        aa_register_saturation_error,
                        aa_axial_register_force_max, aa_radial_force_max);
                    if (candidate_mechanics.length >
                        aa_tether_capture_distance) continue;
                    const int candidate_i = bound_to_i
                        ? bound_actin : candidate_actin;
                    const int candidate_j = bound_to_i
                        ? candidate_actin : bound_actin;
                    const double candidate_s_i = bound_to_i
                        ? bound_s : free_s;
                    const double candidate_s_j = bound_to_i
                        ? free_s : bound_s;
                    const auto exclusion = _local_aa_myosin_exclusion(
                        candidate_i, candidate_j,
                        candidate_s_i, candidate_s_j);
                    if (exclusion.hard_core_overlap ||
                        exclusion.axial_tip_contact) continue;
                    const double local_closing_rate =
                        actn2_arm_kinetics::elastic_rebinding_rate(
                            aa_tether_rebind_rate,
                            candidate_mechanics.energy, thermal_beta,
                            aa_arm_rebind_theta);
                    const double rate = local_closing_rate *
                        boltzmann_exclusion_factor(
                            thermal_beta, exclusion.energy) *
                        bin_width / std::max(1e-12, aa_overlap_capture);
                    if (!(rate > 0.0)) continue;
                    candidates.push_back({candidate_actin, lower_s, upper_s, rate,
                                          candidate_mechanics.energy});
                    total_rate += rate;
                    }
                }
            }
            if (!(total_rate > 0.0)) continue;
            counter_rng::CounterRng rebind_rng(
                static_cast<std::uint64_t>(initial_seed),
                static_cast<std::uint64_t>(current_step),
                counter_rng::Stream::aa_arm_rebinding, link.id);
            if (rebind_rng.uniform() >=
                actn2_arm_kinetics::event_probability(total_rate, dt)) continue;
            double target = rebind_rng.uniform() * total_rate;
            const RebindCandidate* selected = &candidates.back();
            for (const auto& candidate : candidates) {
                target -= candidate.rate;
                if (target <= 0.0) {
                    selected = &candidate;
                    break;
                }
            }
            const double selected_free_s = selected->lower_s +
                (selected->upper_s - selected->lower_s) * rebind_rng.uniform();
            const int selected_i = bound_to_i
                ? bound_actin : selected->actin;
            const int selected_j = bound_to_i
                ? selected->actin : bound_actin;
            const double selected_s_i = bound_to_i
                ? bound_s : selected_free_s;
            const double selected_s_j = bound_to_i
                ? selected_free_s : bound_s;
            const auto selected_exclusion = _local_aa_myosin_exclusion(
                selected_i, selected_j, selected_s_i, selected_s_j);
            if (selected_exclusion.hard_core_overlap ||
                selected_exclusion.axial_tip_contact) continue;
            if (aa_zdisc_tip_contact &&
                selected_exclusion.nearest_tip_myosin >= 0 &&
                selected_exclusion.nearest_tip_gap <
                    aa_zdisc_tip_onset - 1e-12) {
                throw std::runtime_error(
                    "ACTN2 tether rebinding attempted to commit inside a "
                    "myosin-tip exclusion zone");
            }
            if (aa_neighbor_identity_enabled &&
                !aa_neighbor_slots.reserve(
                    bound_actin, selected->actin,
                    static_cast<aa_neighbor_identity::Step>(current_step),
                    max_aa_bonds_per_actin)) {
                continue;
            }
            link.last_partner = old_partner;
            if (bound_to_i) {
                link.actin_j = selected->actin;
                link.s_j = selected_free_s;
            } else {
                link.actin_i = selected->actin;
                link.s_i = selected_free_s;
            }
            link.status = 1;
            link.arm_state = actn2_arm_kinetics::State::weak_bridge;
            link.arm_state_since = static_cast<RecoveryStep>(std::min(
                current_step,
                static_cast<size_t>(std::numeric_limits<RecoveryStep>::max())));
            link.attach_step = link.arm_state_since;
            // Make this closure visible to later tether events in the same
            // serial KMC sweep.  Without the incremental index update, two
            // free arms could select the same pair-facing contour before the
            // end-of-phase rebuild.
            aa_links_by_pair[_aa_pair_key(
                bound_actin, selected->actin)].push_back(slot);
            aa_links_by_actin[static_cast<size_t>(selected->actin)]
                .push_back(slot);
            aa_tip_accessibility_events.insert(
                aa_tip_accessibility_events.end(), {
                    static_cast<double>(link.id),
                    static_cast<double>(current_step), 2.0,
                    static_cast<double>(link.actin_i),
                    static_cast<double>(link.actin_j), link.s_i, link.s_j,
                    static_cast<double>(
                        selected_exclusion.nearest_tip_myosin),
                    static_cast<double>(
                        selected_exclusion.nearest_tip_sign),
                    selected_exclusion.nearest_tip_gap,
                    selected_exclusion.nearest_tip_rho, 0.0,
                    selected_exclusion.nearest_tip_motor_engaged ? 1.0
                                                                  : 0.0});
            aa_arm_events.insert(aa_arm_events.end(), {
                static_cast<double>(link.id),
                static_cast<double>(current_step), 2.0,
                static_cast<double>(static_cast<int>(old_state)),
                static_cast<double>(static_cast<int>(link.arm_state)),
                static_cast<double>(bound_actin),
                static_cast<double>(old_partner),
                static_cast<double>(selected->actin), bound_s,
                selected_free_s, 0.0, 0.0, 0.0, selected->rate});
        }
        _rebuild_local_aa_indices();
    }
    aa_scratch_pair_tension.clear();
    aa_scratch_mechanics.assign(aa_links.size(), LocalAALinkResult{});
    aa_scratch_arm_loads.assign(aa_links.size(),
                                actn2_load_kinetics::ArmLoads{});
    aa_scratch_contact_arm_loads.assign(
        aa_links.size(), std::array<double, 2>{0.0, 0.0});
    aa_scratch_contact_penetration.assign(aa_links.size(), 0.0);
    aa_scratch_contact_results.assign(aa_links.size(), {});
    aa_scratch_contact_evaluated.assign(aa_links.size(), 0);
    aa_scratch_mature_connector_excluded.assign(aa_links.size(), 0);
    aa_scratch_total_arm_load.assign(static_cast<size_t>(actin.n), 0.0);
    aa_scratch_productive_load.assign(static_cast<size_t>(actin.n), 0.0);
    aa_scratch_motor_force.assign(static_cast<size_t>(actin.n),
                                  std::array<double, 3>{0.0, 0.0, 0.0});
    aa_scratch_productive_force.assign(static_cast<size_t>(actin.n),
                                       std::array<double, 3>{0.0, 0.0, 0.0});
    if (aa_scratch_productive_by_myosin.size() !=
        static_cast<size_t>(actin.n)) {
        aa_scratch_productive_by_myosin.resize(static_cast<size_t>(actin.n));
    }
    if (aa_scratch_productive_source_demand.size() !=
        static_cast<size_t>(actin.n)) {
        aa_scratch_productive_source_demand.resize(
            static_cast<size_t>(actin.n));
    }
    }

    auto& mechanics = aa_scratch_mechanics;
    auto& raw_arm_loads = aa_scratch_arm_loads;
    auto& contact_arm_loads = aa_scratch_contact_arm_loads;
    auto& contact_penetration = aa_scratch_contact_penetration;
    auto& contact_results = aa_scratch_contact_results;
    auto& contact_evaluated = aa_scratch_contact_evaluated;
    auto& pair_tension_sum = aa_scratch_pair_tension;
    auto& total_raw_favorable_arm_load = aa_scratch_total_arm_load;
    auto& productive_motor_load = aa_scratch_productive_load;
    const auto& active_slots = aa_active_link_slots;
    const long long n_active = static_cast<long long>(active_slots.size());

    // A weak ACTN2 encounter complex is not stereospecifically locked: each
    // binding domain diffuses in one dimension along its actin's barbed-end
    // association zone. The connector force supplies the Einstein drift.
    // Counter-based noise keeps each physical molecule's trajectory invariant
    // to OpenMP scheduling. Status-2 coordinates are deliberately untouched.
    if (!aa_arm_state_model && aa_status1_site_diffusion > 0.0 &&
        thermal_beta > 0.0 && dt > 0.0) {
        const double lower = -0.5 * actin.length;
        const double zone = aa_barbed_zone > 0.0
            ? std::min(actin.length, aa_barbed_zone) : actin.length;
        const double upper = lower + zone;
        const double sigma = std::sqrt(2.0 * aa_status1_site_diffusion * dt);
        auto reflect = [lower, upper](double value) {
            if (!(upper > lower)) return lower;
            const double width = upper - lower;
            double folded = std::fmod(value - lower, 2.0 * width);
            if (folded < 0.0) folded += 2.0 * width;
            return lower + (folded <= width ? folded
                                             : 2.0 * width - folded);
        };
        stamp_aa_phase(kAaArmKinetics);
        #pragma omp for schedule(static)
        for (long long index = 0; index < n_active; ++index) {
            const size_t slot = active_slots[static_cast<size_t>(index)];
            auto& link = aa_links[slot];
            if (!link.active || link.status != 1) continue;
            const auto weak_mechanics = compute_local_aa_link(
                actin, link.actin_i, link.actin_j, link.s_i, link.s_j, box,
                ramp * aa_link_radial_stiffness(link.status),
                ramp * aa_link_perpendicular_stiffness(link.status), single_link_alignment_k,
                aa_optimal, aa_antiparallel_overlap,
                aa_link_register_ratio(link.status),
                aa_register_saturation_error, aa_axial_register_force_max,
            ramp * aa_radial_force_max);
            counter_rng::CounterRng site_rng(
                static_cast<std::uint64_t>(initial_seed),
                static_cast<std::uint64_t>(current_step),
                counter_rng::Stream::aa_status1_site_diffusion,
                static_cast<std::uint64_t>(link.id));
            const double radius = std::sqrt(
                -2.0 * std::log(site_rng.uniform_pos()));
            const double angle = 6.2831853071795864769 * site_rng.uniform();
            const double gaussian_i = radius * std::cos(angle);
            const double gaussian_j = radius * std::sin(angle);
            const vec force_i = weak_mechanics.force_on_first;
            const vec force_j = -force_i;
            const double drift_i = aa_status1_site_diffusion * thermal_beta *
                force_i.dot(actin.direction[link.actin_i]) * dt;
            const double drift_j = aa_status1_site_diffusion * thermal_beta *
                force_j.dot(actin.direction[link.actin_j]) * dt;
            const double proposed_i = reflect(
                link.s_i + drift_i + sigma * gaussian_i);
            const double proposed_j = reflect(
                link.s_j + drift_j + sigma * gaussian_j);
            const auto old_exclusion = _local_aa_myosin_exclusion(
                link.actin_i, link.actin_j, link.s_i, link.s_j);
            const auto proposed_exclusion = _local_aa_myosin_exclusion(
                link.actin_i, link.actin_j, proposed_i, proposed_j);
            // A no-flux Brownian boundary is enforced independently for each
            // physical myosin tip. A move cannot enter a new tip zone, and a
            // site already inside a zone may only reduce that same tip's
            // penetration. Keying by (myosin,end_sign) prevents a move from
            // leaving one excluded zone by silently entering another.
            const bool exits_or_reduces_axial_contact =
                aa_tip_diffusion_move_allowed(
                    old_exclusion, proposed_exclusion);
            if (!proposed_exclusion.hard_core_overlap &&
                exits_or_reduces_axial_contact) {
                link.s_i = proposed_i;
                link.s_j = proposed_j;
            }
        }
    }

    // Phase 0.  The motor-force vectors of an actin are read-only for the whole
    // routine, so caching them per actin replaces the former two evaluations
    // per incident link (each of which rescanned every head of every bound
    // myosin) with one evaluation per actin.
    double motor_support_build_thread_seconds = 0.0;
    stamp_aa_phase(kAaSiteDiffusion);
    #pragma omp for schedule(static)
    for (int i = 0; i < actin.n; ++i) {
        aa_scratch_motor_force[static_cast<size_t>(i)] =
            _actin_motor_force_vector(i);
        if (motor_transmitted_mode) {
            const double started = omp_get_wtime();
            auto& contributions =
                aa_scratch_productive_by_myosin[static_cast<size_t>(i)];
            contributions = _actin_productive_motor_contributions(i);
            aa_scratch_productive_source_demand[
                static_cast<size_t>(i)].assign(contributions.size(), 0.0);
            std::array<double, 3> aggregate{0.0, 0.0, 0.0};
            for (const auto& contribution : contributions) {
                aggregate[0] += contribution.force[0];
                aggregate[1] += contribution.force[1];
                aggregate[2] += contribution.force[2];
            }
            aa_scratch_productive_force[static_cast<size_t>(i)] = aggregate;
            motor_support_build_thread_seconds += omp_get_wtime() - started;
        } else {
            aa_scratch_productive_by_myosin[static_cast<size_t>(i)].clear();
            aa_scratch_productive_source_demand[
                static_cast<size_t>(i)].clear();
            aa_scratch_productive_force[static_cast<size_t>(i)] =
                _actin_productive_motor_force_vector(i);
        }
        productive_motor_load[static_cast<size_t>(i)] =
            actn2_load_kinetics::norm(
                aa_scratch_productive_force[static_cast<size_t>(i)]);
    }
    if (motor_transmitted_mode) {
        #pragma omp atomic update
        aa_motor_support_build_seconds +=
            motor_support_build_thread_seconds;
    }

    // Phase 1a.  Unloaded mechanics and signed arm loads are independent per
    // link, so they are computed by the whole team.
    stamp_aa_phase(kAaMotorForces);
    #pragma omp for schedule(static)
    for (long long index = 0; index < n_active; ++index) {
        const size_t slot = active_slots[static_cast<size_t>(index)];
        const auto& link = aa_links[slot];
        if (aa_arm_state_model &&
            !actn2_arm_kinetics::is_bridge(link.arm_state)) {
            mechanics[slot] = LocalAALinkResult{};
            raw_arm_loads[slot] = actn2_load_kinetics::ArmLoads{};
            continue;
        }
        const double link_radial_k =
            ramp * aa_link_radial_stiffness(link.status);
        // Geometry is read-only after site diffusion. Evaluate established
        // mature links with the team, then consume the flags in canonical
        // slot order below. No RNG, event or floating-point reduction moves.
        if (link.status == 2) {
            aa_scratch_mature_connector_excluded[slot] =
                _aa_status2_connector_excluded(
                    link.actin_i, link.actin_j, link.s_i, link.s_j,
                    actn2_maturation::deletion_buffer);
        }
        mechanics[slot] = compute_local_aa_link(
            actin, link.actin_i, link.actin_j, link.s_i, link.s_j, box,
            link_radial_k, ramp * aa_link_perpendicular_stiffness(link.status),
            single_link_alignment_k, aa_optimal, aa_antiparallel_overlap,
            aa_link_register_ratio(link.status),
            aa_register_saturation_error, aa_axial_register_force_max,
            ramp * aa_radial_force_max);
        const std::array<double, 3> link_force{
            mechanics[slot].force_on_first.x,
            mechanics[slot].force_on_first.y,
            mechanics[slot].force_on_first.z};
        const std::array<double, 3> direction_i{
            actin.direction[link.actin_i].x,
            actin.direction[link.actin_i].y,
            actin.direction[link.actin_i].z};
        const std::array<double, 3> direction_j{
            actin.direction[link.actin_j].x,
            actin.direction[link.actin_j].y,
            actin.direction[link.actin_j].z};
        raw_arm_loads[slot] = actn2_load_kinetics::signed_arm_loads(
            link_force, direction_i, direction_j);
        if ((aa_myosin_contact_epsilon > 0.0 ||
             (aa_myosin_contact_shell_width > 0.0 &&
              aa_myosin_contact_shell_energy > 0.0)) &&
            aa_actn2_contact_radius > 0.0) {
            contact_evaluated[slot] = 1;
            std::vector<int> contact_myosins =
                actin_neighbors_by_species[static_cast<size_t>(link.actin_i)].second;
            const auto& second_neighbors =
                actin_neighbors_by_species[static_cast<size_t>(link.actin_j)].second;
            contact_myosins.insert(contact_myosins.end(),
                                    second_neighbors.begin(), second_neighbors.end());
            std::sort(contact_myosins.begin(), contact_myosins.end());
            contact_myosins.erase(
                std::unique(contact_myosins.begin(), contact_myosins.end()),
                contact_myosins.end());
            const auto patch = compute_actn2_patch_geometry(
                actin, link.actin_i, link.actin_j, link.s_i, link.s_j, box);
            const double outer_contact_distance =
                myosin.radius + aa_actn2_contact_radius +
                aa_myosin_contact_shell_width;
            for (int myosin_index : contact_myosins) {
                if (!actn2_patch_may_contact_myosin_body(
                        patch, myosin, myosin_index, box,
                        outer_contact_distance)) continue;
                const auto contact = compute_actn2_myosin_contact(
                    actin, myosin, link.actin_i, link.actin_j,
                    link.s_i, link.s_j, myosin_index, box,
                    myosin.radius + aa_actn2_contact_radius,
                    aa_myosin_contact_epsilon,
                    aa_myosin_contact_force_cap,
                    aa_myosin_contact_shell_width,
                    aa_myosin_contact_shell_energy);
                if (!contact.applied) continue;
                contact_results[slot].push_back({myosin_index, contact});
                contact_arm_loads[slot][0] +=
                    contact.arm_weight_i * contact.force;
                contact_arm_loads[slot][1] +=
                    contact.arm_weight_j * contact.force;
                contact_penetration[slot] = std::max(
                    contact_penetration[slot], contact.penetration);
            }
        }
    }

    // Phase 1b.  The two accumulators are floating-point reductions over
    // overlapping actin/pair keys.  They stay serial and in slot order so the
    // summation order is exactly the original one; the work is a handful of
    // adds per link and does not scale with head count.
    stamp_aa_phase(kAaArmLoads);
    #pragma omp single
    {
    for (size_t slot : active_slots) {
        const auto& link = aa_links[slot];
        if (aa_arm_state_model &&
            !actn2_arm_kinetics::is_bridge(link.arm_state)) continue;
        pair_tension_sum[_aa_pair_key(link.actin_i, link.actin_j)] +=
            std::max(0.0, ramp * aa_link_radial_stiffness(link.status) *
                (mechanics[slot].length - aa_optimal));
        total_raw_favorable_arm_load[static_cast<size_t>(link.actin_i)] +=
            std::max(0.0, raw_arm_loads[slot].first);
        total_raw_favorable_arm_load[static_cast<size_t>(link.actin_j)] +=
            std::max(0.0, raw_arm_loads[slot].second);
        if (motor_transmitted_mode) {
            const std::array<double, 3> connector_i{
                mechanics[slot].force_on_first.x,
                mechanics[slot].force_on_first.y,
                mechanics[slot].force_on_first.z};
            const std::array<double, 3> connector_j{
                -connector_i[0], -connector_i[1], -connector_i[2]};
            const auto& sources_i = aa_scratch_productive_by_myosin[
                static_cast<size_t>(link.actin_i)];
            const auto& sources_j = aa_scratch_productive_by_myosin[
                static_cast<size_t>(link.actin_j)];
            auto& demand_i = aa_scratch_productive_source_demand[
                static_cast<size_t>(link.actin_i)];
            auto& demand_j = aa_scratch_productive_source_demand[
                static_cast<size_t>(link.actin_j)];
            for (size_t source = 0; source < sources_i.size(); ++source) {
                demand_i[source] += actn2_load_kinetics::opposing_support(
                    sources_i[source].force, connector_i);
            }
            for (size_t source = 0; source < sources_j.size(); ++source) {
                demand_j[source] += actn2_load_kinetics::opposing_support(
                    sources_j[source].force, connector_j);
            }
        }
    }
    }

    // Phase 2. Catch-state classification and deletion bookkeeping commit
    // in canonical slot order. Read-only mechanics and mature-connector
    // geometry were evaluated by the team in phase 1a.
    double motor_candidate_thread_seconds = 0.0;
    std::uint64_t motor_candidate_thread_pairs = 0;
    stamp_aa_phase(kAaLinkLoop);
    // Serial by measurement, not by necessity: the loop holds only ~0.1 s of
    // work per 100 steps at production density, and `ordered` lock-step cost
    // 0.27 s on top of that at 8-32 threads.  Too cheap to parallelise.
    #pragma omp single
    {
    for (long long index = 0; index < n_active; ++index) {
        const size_t slot = active_slots[static_cast<size_t>(index)];
        auto& link = aa_links[slot];
        if (aa_arm_state_model &&
            !actn2_arm_kinetics::is_bridge(link.arm_state)) {
            continue;
        }
        const auto& mech = mechanics[slot];
        // This scalar connector tension is retained only for the legacy
        // shared-load model.  The production arm_load_productive classifier
        // below uses the signed force delivered to each actin arm together
        // with independently measured opposed productive motor support.
        const double link_tension = std::max(
            0.0, ramp * aa_link_radial_stiffness(link.status) *
                (mech.length - aa_optimal));
        const auto& motor_i =
            aa_scratch_motor_force[static_cast<size_t>(link.actin_i)];
        const auto& motor_j =
            aa_scratch_motor_force[static_cast<size_t>(link.actin_j)];
        const double dx = motor_i[0] - motor_j[0];
        const double dy = motor_i[1] - motor_j[1];
        const double dz = motor_i[2] - motor_j[2];
        const double opposed = 0.5 * std::sqrt(dx * dx + dy * dy + dz * dz);
        // Read-only lookup: phase 1b already inserted every active pair key,
        // and `find` (unlike `operator[]`) cannot mutate the shared map from
        // several threads at once.
        const auto tension_it =
            pair_tension_sum.find(_aa_pair_key(link.actin_i, link.actin_j));
        const double denominator = tension_it == pair_tension_sum.end()
            ? 0.0 : tension_it->second;
        const double weight = denominator > 0.0
            ? link_tension / denominator : 0.0;
        const auto arm_loads = raw_arm_loads[slot];
        double motor_transmitted_support = 0.0;
        bool motor_pair_committed = false;
        const actn2_load_kinetics::MotorPairIdentity old_candidate_identity{
            link.candidate_myosin_i, link.candidate_myosin_j,
            link.candidate_half_i, link.candidate_half_j};
        const bool candidate_was_locked = link.candidate_locked;
        actn2_load_kinetics::MotorPairEvent motor_pair_event =
            actn2_load_kinetics::MotorPairEvent::none;
        double motor_pair_best_support = 0.0;
        double motor_pair_current_support = 0.0;
        double motor_pair_candidate_age = 0.0;
        if (motor_transmitted_mode) {
            const double started = omp_get_wtime();
            const std::array<double, 3> connector_i{
                mech.force_on_first.x, mech.force_on_first.y,
                mech.force_on_first.z};
            const std::array<double, 3> connector_j{
                -connector_i[0], -connector_i[1], -connector_i[2]};
            const auto& sources_i = aa_scratch_productive_by_myosin[
                static_cast<size_t>(link.actin_i)];
            const auto& sources_j = aa_scratch_productive_by_myosin[
                static_cast<size_t>(link.actin_j)];
            const auto& demand_i = aa_scratch_productive_source_demand[
                static_cast<size_t>(link.actin_i)];
            const auto& demand_j = aa_scratch_productive_source_demand[
                static_cast<size_t>(link.actin_j)];
            double best_support = 0.0;
            double current_support = 0.0;
            double best_raw_i = 0.0;
            double best_raw_j = 0.0;
            double best_allocated_i = 0.0;
            double best_allocated_j = 0.0;
            double current_raw_i = 0.0;
            double current_raw_j = 0.0;
            double current_allocated_i = 0.0;
            double current_allocated_j = 0.0;
            int best_myosin_i = -1;
            int best_myosin_j = -1;
            int8_t best_half_i = 0;
            int8_t best_half_j = 0;
            auto identity_less = [](int mi, int8_t hi, int mj, int8_t hj,
                                    int other_mi, int8_t other_hi,
                                    int other_mj, int8_t other_hj) {
                return std::tie(mi, hi, mj, hj) <
                       std::tie(other_mi, other_hi, other_mj, other_hj);
            };
            for (size_t source_i_index = 0;
                 source_i_index < sources_i.size(); ++source_i_index) {
                const auto& source_i = sources_i[source_i_index];
                const double raw_support_i =
                    actn2_load_kinetics::opposing_support(
                        source_i.force, connector_i);
                const double support_i =
                    actn2_load_kinetics::allocated_source_support(
                        raw_support_i, demand_i[source_i_index],
                        actn2_load_kinetics::norm(source_i.force));
                for (size_t source_j_index = 0;
                     source_j_index < sources_j.size(); ++source_j_index) {
                    const auto& source_j = sources_j[source_j_index];
                    // The same thick filament pulling both thin filaments is
                    // not an M-band registration signal.
                    if (!actn2_load_kinetics::distinct_motor_sources(
                            source_i.myosin, source_j.myosin)) {
                        continue;
                    }
                    ++motor_candidate_thread_pairs;
                    const double raw_support_j =
                        actn2_load_kinetics::opposing_support(
                            source_j.force, connector_j);
                    const double support_j =
                        actn2_load_kinetics::allocated_source_support(
                            raw_support_j, demand_j[source_j_index],
                            actn2_load_kinetics::norm(source_j.force));
                    const double support =
                        actn2_load_kinetics::motor_transmitted_load(
                            arm_loads, support_i, support_j);
                    const bool is_current =
                        source_i.myosin == link.candidate_myosin_i &&
                        source_j.myosin == link.candidate_myosin_j &&
                        source_i.half == link.candidate_half_i &&
                        source_j.half == link.candidate_half_j;
                    if (is_current) {
                        current_support = support;
                        current_raw_i = raw_support_i;
                        current_raw_j = raw_support_j;
                        current_allocated_i = support_i;
                        current_allocated_j = support_j;
                    }
                    if (support > best_support ||
                        (support == best_support && support > 0.0 &&
                         (best_myosin_i < 0 || identity_less(
                              source_i.myosin, source_i.half,
                              source_j.myosin, source_j.half,
                              best_myosin_i, best_half_i,
                              best_myosin_j, best_half_j)))) {
                        best_support = support;
                        best_myosin_i = source_i.myosin;
                        best_myosin_j = source_j.myosin;
                        best_half_i = source_i.half;
                        best_half_j = source_j.half;
                        best_raw_i = raw_support_i;
                        best_raw_j = raw_support_j;
                        best_allocated_i = support_i;
                        best_allocated_j = support_j;
                    }
                }
            }

            actn2_load_kinetics::MotorPairState candidate_state;
            candidate_state.identity = old_candidate_identity;
            candidate_state.since_step = link.candidate_since;
            candidate_state.last_seen_step = link.candidate_last_seen;
            candidate_state.support_ema = link.candidate_support_ema;
            candidate_state.support_instantaneous =
                link.candidate_support_instantaneous;
            candidate_state.challenger_identity = {
                link.challenger_myosin_i, link.challenger_myosin_j,
                link.challenger_half_i, link.challenger_half_j};
            candidate_state.challenger_since_step = link.challenger_since;
            candidate_state.challenger_support_ema =
                link.challenger_support_ema;
            candidate_state.locked = link.candidate_locked;
            const actn2_load_kinetics::MotorPairIdentity best_identity{
                best_myosin_i, best_myosin_j, best_half_i, best_half_j};
            const auto candidate_update =
                actn2_load_kinetics::update_motor_pair_candidate(
                    candidate_state, best_identity, best_support,
                    current_support, static_cast<std::uint64_t>(current_step),
                    dt, aa_motor_pair_support_window,
                    aa_motor_pair_commit_time, aa_motor_pair_switch_ratio,
                    aa_motor_pair_exchange_time, aa_motor_pair_dropout_time);
            link.candidate_myosin_i = candidate_state.identity.myosin_i;
            link.candidate_myosin_j = candidate_state.identity.myosin_j;
            link.candidate_half_i = candidate_state.identity.half_i;
            link.candidate_half_j = candidate_state.identity.half_j;
            link.candidate_since = static_cast<RecoveryStep>(std::min(
                candidate_state.since_step,
                static_cast<std::uint64_t>(
                    std::numeric_limits<RecoveryStep>::max())));
            link.candidate_last_seen = static_cast<RecoveryStep>(std::min(
                candidate_state.last_seen_step,
                static_cast<std::uint64_t>(
                    std::numeric_limits<RecoveryStep>::max())));
            link.candidate_support_ema = candidate_state.support_ema;
            link.candidate_support_instantaneous =
                candidate_state.support_instantaneous;
            link.challenger_myosin_i =
                candidate_state.challenger_identity.myosin_i;
            link.challenger_myosin_j =
                candidate_state.challenger_identity.myosin_j;
            link.challenger_half_i =
                candidate_state.challenger_identity.half_i;
            link.challenger_half_j =
                candidate_state.challenger_identity.half_j;
            link.challenger_since = static_cast<RecoveryStep>(std::min(
                candidate_state.challenger_since_step,
                static_cast<std::uint64_t>(
                    std::numeric_limits<RecoveryStep>::max())));
            link.challenger_support_ema =
                candidate_state.challenger_support_ema;
            motor_pair_event = candidate_update.event;
            motor_pair_candidate_age = candidate_update.candidate_age;
            motor_pair_committed = candidate_update.committed;
            motor_pair_best_support = best_support;
            motor_pair_current_support = current_support;
            if (motor_pair_event ==
                    actn2_load_kinetics::MotorPairEvent::selected ||
                motor_pair_event ==
                    actn2_load_kinetics::MotorPairEvent::switched) {
                link.candidate_raw_support_i = best_raw_i;
                link.candidate_raw_support_j = best_raw_j;
                link.candidate_allocated_support_i = best_allocated_i;
                link.candidate_allocated_support_j = best_allocated_j;
            } else if (candidate_state.identity.valid() &&
                       !candidate_state.locked) {
                link.candidate_raw_support_i = current_raw_i;
                link.candidate_raw_support_j = current_raw_j;
                link.candidate_allocated_support_i = current_allocated_i;
                link.candidate_allocated_support_j = current_allocated_j;
            } else if (!candidate_state.identity.valid()) {
                link.candidate_raw_support_i = 0.0;
                link.candidate_raw_support_j = 0.0;
                link.candidate_allocated_support_i = 0.0;
                link.candidate_allocated_support_j = 0.0;
            }
            motor_transmitted_support = link.candidate_support_ema;
            motor_candidate_thread_seconds += omp_get_wtime() - started;
        }
        double allocated_motor_i = 0.0;
        double allocated_motor_j = 0.0;
        auto kinetic_arm_loads = arm_loads;
        if (aa_status_mode ==
            actn2_load_kinetics::StatusMode::arm_load_productive) {
            allocated_motor_i = actn2_load_kinetics::allocated_motor_load(
                arm_loads.first,
                total_raw_favorable_arm_load[static_cast<size_t>(link.actin_i)],
                productive_motor_load[static_cast<size_t>(link.actin_i)]);
            allocated_motor_j = actn2_load_kinetics::allocated_motor_load(
                arm_loads.second,
                total_raw_favorable_arm_load[static_cast<size_t>(link.actin_j)],
                productive_motor_load[static_cast<size_t>(link.actin_j)]);
            kinetic_arm_loads = {
                std::min(std::max(0.0, arm_loads.first), allocated_motor_i),
                std::min(std::max(0.0, arm_loads.second), allocated_motor_j)};
        }
        // Identical to _aa_productive_opposed_load(i, j) with the two
        // productive force vectors taken from the phase-0 cache.
        const double instantaneous_support =
            actn2_load_kinetics::productive_opposed_load(
                std::array<double, 3>{actin.direction[link.actin_i].x,
                                      actin.direction[link.actin_i].y,
                                      actin.direction[link.actin_i].z},
                std::array<double, 3>{actin.direction[link.actin_j].x,
                                      actin.direction[link.actin_j].y,
                                      actin.direction[link.actin_j].z},
                aa_scratch_productive_force[
                    static_cast<size_t>(link.actin_i)],
                aa_scratch_productive_force[
                    static_cast<size_t>(link.actin_j)]);
        // Status-2 formation previously required registry AND bilateral opposed
        // motor load in the SAME timestep.  With ~6 ms strong-head lifetimes that
        // coincidence is rare, so links stalled at status 1 indefinitely.  Average
        // the productive support over a motor timescale instead: a registry-
        // qualified pair may mature on load accumulated across a short window.
        // The catch-bond load still uses the instantaneous value.  The EMA may
        // additionally suppress the opt-in mature-to-weak relaxation below,
        // allowing it to bridge brief head-turnover gaps without storing load
        // indefinitely.
        if (aa_status2_support_window > 0.0) {
            const double alpha = std::min(
                1.0, dt / std::max(aa_status2_support_window, 1e-12));
            link.support_ema += alpha * (instantaneous_support - link.support_ema);
        } else {
            link.support_ema = instantaneous_support;
        }
        const double productive_support =
            (link.status == 2) ? instantaneous_support
                               : std::max(instantaneous_support, link.support_ema);
        double effective_load = std::min(link_tension, weight * opposed);
        if (aa_status_mode ==
                actn2_load_kinetics::StatusMode::productive_opposed ||
            aa_status_mode ==
                actn2_load_kinetics::StatusMode::productive_register) {
            effective_load = productive_support;
        } else if (aa_status_mode ==
                   actn2_load_kinetics::StatusMode::arm_load_productive) {
            effective_load = kinetic_arm_loads.favorable_bilateral();
        }
        const BondFlag old_status = link.status;
        const bool registry_ok = actn2_load_kinetics::registry_gate(
            static_cast<int>(old_status), mech.direction_dot,
            mech.registry_error, aa_status2_registry_on,
            aa_status2_registry_off);
        const bool registry_retention_ok = actn2_load_kinetics::registry_gate(
            2, mech.direction_dot, mech.registry_error,
            aa_status2_registry_on, aa_status2_registry_off);
        bool instantaneous_status2_gate = false;
        bool lock_motor_candidate = false;
        if (aa_status_mode ==
            actn2_load_kinetics::StatusMode::productive_register) {
            const int classified =
                actn2_load_kinetics::classify_productive_register(
                    static_cast<int>(old_status), productive_support,
                    aa_status2_force_on, mech.direction_dot,
                    mech.registry_error, aa_status2_registry_on,
                    aa_status2_registry_off);
            instantaneous_status2_gate = old_status != 2 && classified == 2;
            link.status = static_cast<BondFlag>(classified);
            if (link.status != 2) effective_load = 0.0;
        } else if (motor_transmitted_mode) {
            const int classified =
                actn2_load_kinetics::classify_motor_transmitted_register(
                    static_cast<int>(old_status),
                    motor_transmitted_support, motor_pair_committed,
                    aa_status2_force_on, mech.direction_dot,
                    mech.registry_error, aa_status2_registry_on,
                    aa_status2_registry_off);
            instantaneous_status2_gate = old_status != 2 && classified == 2;
            link.status = static_cast<BondFlag>(classified);
            lock_motor_candidate = link.status == 2;
            effective_load = link.status == 2
                ? arm_loads.favorable_bilateral() : 0.0;
        } else if (aa_status_mode ==
                   actn2_load_kinetics::StatusMode::arm_load_productive) {
            const int classified =
                actn2_load_kinetics::classify_arm_loaded_productive(
                    static_cast<int>(link.status), kinetic_arm_loads,
                    productive_support, aa_status2_force_on,
                    aa_status2_force_off);
            instantaneous_status2_gate = classified == 2 && registry_ok;
            link.status = static_cast<BondFlag>(
                classified == 2 && registry_ok ? 2 : 1);
            if (link.status != 2) effective_load = 0.0;
        } else {
            const double threshold = link.status == 2
                ? aa_status2_force_off : aa_status2_force_on;
            link.status = effective_load > threshold && registry_ok ? 2 : 1;
        }
        bool status2_unlocked = false;
        if (old_status == 2 && aa_status2_unlock_rate > 0.0) {
            // Status 2 is a reversible ACTN2 conformation.  It survives brief
            // head turnover through the existing motor-timescale filter, but
            // sustained loss of source-resolved support allows a stochastic
            // return to the weak, annealable state without detaching ACTN2.
            const double filtered_support = motor_transmitted_mode
                ? link.candidate_support_ema : link.support_ema;
            const double unlock_rate = actn2_load_kinetics::mature_unlock_rate(
                aa_status2_unlock_rate, filtered_support,
                aa_status2_unlock_force_scale, mech.registry_error,
                aa_status2_registry_off, aa_status2_unlock_registry_scale);
            const double unlock_probability =
                actn2_load_kinetics::mature_unlock_probability(unlock_rate, dt);
            counter_rng::CounterRng unlock_rng(
                static_cast<std::uint64_t>(initial_seed),
                static_cast<std::uint64_t>(current_step),
                counter_rng::Stream::aa_status2_unlock, link.id);
            status2_unlocked = unlock_rng.uniform() < unlock_probability;
            if (status2_unlocked) {
                link.status = 1;
                link.status2_memory_until = 0;
                link.candidate_locked = false;
                effective_load = 0.0;
            } else {
                link.status = 2;
                if (aa_status_mode ==
                        actn2_load_kinetics::StatusMode::productive_opposed ||
                    aa_status_mode ==
                        actn2_load_kinetics::StatusMode::productive_register) {
                    effective_load = instantaneous_support;
                } else if (motor_transmitted_mode) {
                    effective_load = arm_loads.favorable_bilateral();
                } else if (aa_status_mode ==
                           actn2_load_kinetics::StatusMode::arm_load_productive) {
                    effective_load = kinetic_arm_loads.favorable_bilateral();
                }
            }
        }
        if (aa_status2_unlock_rate <= 0.0 &&
            (aa_status2_locked_until_detach || aa_arm_state_model) &&
            old_status == 2) {
            // Stereospecific locking is a molecular state, not an
            // instantaneous force threshold.  Head turnover may remove the
            // load that created the lock, but it cannot make the binding
            // domains diffuse again without an actual ACTN2 unbinding event.
            // The status-2 lifetime below still responds to the current load,
            // and contour rupture remains an immediate escape path.
            link.status = 2;
        }
        const bool register_retention_mode =
            aa_status_mode ==
                actn2_load_kinetics::StatusMode::productive_register ||
            motor_transmitted_mode;
        if (!register_retention_mode &&
            link.status == 2 && am_status2_memory_steps > 0) {
            link.status2_memory_until = static_cast<RecoveryStep>(std::min(
                static_cast<size_t>(std::numeric_limits<RecoveryStep>::max()),
                current_step + am_status2_memory_steps));
        } else if (!register_retention_mode &&
                   link.status == 1 && registry_retention_ok &&
                   static_cast<size_t>(link.status2_memory_until) >= current_step) {
            if (aa_status_mode !=
                    actn2_load_kinetics::StatusMode::arm_load_productive ||
                actn2_load_kinetics::classify_existing_bond(
                    2, kinetic_arm_loads, aa_status2_force_on,
                    aa_status2_force_off) == 2) {
                // In the arm-resolved mode memory may bridge a brief loss of
                // strong-head support only while the physical link continues
                // to carry bilateral load. It cannot stabilize an unloaded
                // link or create status 2 from status 1.
                link.status = 2;
                effective_load = kinetic_arm_loads.favorable_bilateral();
            }
        }
        // All routes into status 2 (including status memory) pass through
        // the same material-connector gate.
        if (old_status != 2 && link.status == 2 &&
            _aa_status2_connector_excluded(
                link.actin_i, link.actin_j, link.s_i, link.s_j)) {
            link.status = 1;
            instantaneous_status2_gate = false;
            link.status2_memory_until = 0;
            if (aa_status_mode ==
                    actn2_load_kinetics::StatusMode::productive_register ||
                motor_transmitted_mode ||
                aa_status_mode ==
                    actn2_load_kinetics::StatusMode::arm_load_productive)
                effective_load = 0.0;
        }
        // Commit motor identity locking only after geometric admission.
        if (lock_motor_candidate && link.status == 2) {
            link.candidate_locked = true;
            link.challenger_myosin_i = -1;
            link.challenger_myosin_j = -1;
            link.challenger_half_i = 0;
            link.challenger_half_j = 0;
            link.challenger_since = 0;
            link.challenger_support_ema = 0.0;
        }
        if (aa_arm_state_model) {
            if (link.status == 2) {
                if (old_status != 2) {
                    link.arm_state_since = static_cast<RecoveryStep>(std::min(
                        current_step,
                        static_cast<size_t>(
                            std::numeric_limits<RecoveryStep>::max())));
                }
                link.arm_state = actn2_arm_kinetics::State::mature_bridge;
            } else {
                if (old_status == 2) {
                    link.arm_state_since = static_cast<RecoveryStep>(std::min(
                        current_step,
                        static_cast<size_t>(
                            std::numeric_limits<RecoveryStep>::max())));
                }
                link.arm_state = actn2_arm_kinetics::State::weak_bridge;
            }
        }
        mechanics[slot] = compute_local_aa_link(
            actin, link.actin_i, link.actin_j, link.s_i, link.s_j, box,
            ramp * aa_link_radial_stiffness(link.status),
            ramp * aa_link_perpendicular_stiffness(link.status),
            single_link_alignment_k, aa_optimal, aa_antiparallel_overlap,
            aa_link_register_ratio(link.status),
            aa_register_saturation_error, aa_axial_register_force_max,
            ramp * aa_radial_force_max);
        const double lifetime = catch_bond_kinetics::lifetime_for_status(
            link.status, base_lifetime, directional_base_lifetime,
            lifetime_coeff, effective_load, stabilize_status1);
        counter_rng::CounterRng detachment_rng(
            static_cast<std::uint64_t>(initial_seed),
            static_cast<std::uint64_t>(current_step),
            counter_rng::Stream::aa_link_detachment,
            static_cast<std::uint64_t>(link.id));
        const double registry_threshold = old_status == 2
            ? aa_status2_registry_off : aa_status2_registry_on;
        const double arm_load_threshold = old_status == 2
            ? aa_status2_force_off : aa_status2_force_on;
        const bool arm_load_ok =
            kinetic_arm_loads.first > arm_load_threshold &&
            kinetic_arm_loads.second > arm_load_threshold;
        const bool productive_support_ok = productive_support > 0.0;
        // Mature links use the same full-connector cylinder as admission,
        // eroded by a small penetration tolerance radially and axially.
        // Test old status too: loss of support must not evade geometric
        // deletion on the very step a mature connector enters the cylinder.
        // This is immediate reason-5 deletion, independent of the legacy
        // barbed-tip hazard, dwell, force switch and arm-state model.
        catch_bond_kinetics::DetachmentDecision detachment;
        const bool mature_connector = old_status == 2 || link.status == 2;
        const bool connector_excluded = old_status == 2
            ? aa_scratch_mature_connector_excluded[slot] != 0
            : (link.status == 2 && _aa_status2_connector_excluded(
                link.actin_i, link.actin_j, link.s_i, link.s_j,
                actn2_maturation::deletion_buffer));
        if (mature_connector && connector_excluded) {
            detachment.mature_connector_excluded = true;
            detachment.probability = 1.0;
        }
        if (!mature_connector && !aa_arm_state_model &&
            aa_overlap_occlusion_veto) {
            if (aa_occlusion_rupture_rate > 0.0) {
                // Burial is a hazard, not a guillotine.  A dedicated RNG
                // stream keeps aa_link_detachment's draw sequence intact,
                // so rate0 = 0 remains bitwise inert.
                counter_rng::CounterRng occlusion_rng(
                    static_cast<std::uint64_t>(initial_seed),
                    static_cast<std::uint64_t>(current_step),
                    counter_rng::Stream::aa_occlusion_rupture,
                    static_cast<std::uint64_t>(link.id));
                detachment = catch_bond_kinetics::
                    decide_barbed_tip_capsule_penetration_rate(
                        actin["aa_barbed_tip_kill_penetration"][link.actin_i],
                        actin["aa_barbed_tip_kill_penetration"][link.actin_j],
                        aa_occlusion_rupture_penetration,
                        aa_occlusion_rupture_rate, dt,
                        [&]() { return occlusion_rng.uniform(); });
            } else {
            detachment =
                catch_bond_kinetics::decide_barbed_tip_capsule_penetration(
                    actin["aa_barbed_tip_kill_penetration"][link.actin_i],
                    actin["aa_barbed_tip_kill_penetration"][link.actin_j],
                    aa_occlusion_rupture_penetration);
            }
            // Require the burial to PERSIST.  A tip that pokes in and bounces
            // straight back out must not execute a crosslink whose natural
            // lifetime is three orders of magnitude longer than the excursion.
            if (aa_occlusion_rupture_rate > 0.0) {
                // The hazard already spares transients; the dwell counter
                // would only advance on a firing draw, which is not a
                // dwell.  The two mechanisms are alternatives.
                link.occlusion_dwell = 0;
            } else if (detachment.am_overlap_occluded) {
                if (link.occlusion_dwell <
                        std::numeric_limits<RecoveryStep>::max()) {
                    link.occlusion_dwell += 1;
                }
                if (static_cast<size_t>(link.occlusion_dwell) <=
                        aa_occlusion_rupture_dwell_steps) {
                    detachment.am_overlap_occluded = false;
                    detachment.probability = 0.0;
                }
            } else {
                link.occlusion_dwell = 0;
            }
        }
        actn2_arm_kinetics::OpeningDecision arm_opening;
        const auto old_arm_state = link.arm_state;
        const double tensile_arm_i = link_tension + contact_arm_loads[slot][0];
        const double tensile_arm_j = link_tension + contact_arm_loads[slot][1];
        link.arm_force_i = tensile_arm_i;
        link.arm_force_j = tensile_arm_j;
        link.contact_force = contact_arm_loads[slot][0] +
                             contact_arm_loads[slot][1];
        link.contact_penetration = contact_penetration[slot];
        if (!detachment.detaches() &&
            aa_occlusion_rupture_penetration > 0.0 &&
            contact_penetration[slot] >= aa_occlusion_rupture_penetration) {
            detachment.actn2_myosin_collision = true;
            detachment.probability = 1.0;
        }
        if (aa_arm_state_model && !detachment.detaches() &&
            aa_link_rupture_length > 0.0 &&
            mechanics[slot].length > aa_link_rupture_length) {
            detachment.geometric = true;
            detachment.probability = 1.0;
        }
        const bool attachment_occluded_i = aa_arm_state_model &&
            actn2_arm_kinetics::contour_fully_occluded(
                actin_crosslink_ratio[link.actin_i], EPS);
        const bool attachment_occluded_j = aa_arm_state_model &&
            actn2_arm_kinetics::contour_fully_occluded(
                actin_crosslink_ratio[link.actin_j], EPS);
        if (aa_arm_state_model && !detachment.detaches() &&
            (attachment_occluded_i || attachment_occluded_j)) {
            // A scalar partial-overlap boundary is not itself a molecular
            // collision.  Partial overlap is handled by the explicit local
            // ACTN2-myosin contact above, whose force enters the Bell arm
            // hazard.  Deterministic opening remains only when an actin has
            // no crosslinkable contour at all, as required by the geometric
            // model.  If both contours vanish, deterministic link-id ordering
            // chooses i and the remaining tether is removed next step.
            arm_opening.opens_i = attachment_occluded_i;
            arm_opening.opens_j = !attachment_occluded_i && attachment_occluded_j;
            arm_opening.rate_i = attachment_occluded_i
                ? std::numeric_limits<double>::infinity() : 0.0;
            arm_opening.rate_j = attachment_occluded_j
                ? std::numeric_limits<double>::infinity() : 0.0;
            arm_opening.probability = 1.0;
        } else if (aa_arm_state_model && !detachment.detaches()) {
            const double zero_force_rate = link.status == 2
                ? aa_arm_open_rate_mature : aa_arm_open_rate_weak;
            const double rate_i = actn2_arm_kinetics::opening_rate(
                zero_force_rate, tensile_arm_i, aa_arm_bell_distance,
                thermal_beta);
            const double rate_j = actn2_arm_kinetics::opening_rate(
                zero_force_rate, tensile_arm_j, aa_arm_bell_distance,
                thermal_beta);
            counter_rng::CounterRng opening_rng(
                static_cast<std::uint64_t>(initial_seed),
                static_cast<std::uint64_t>(current_step),
                counter_rng::Stream::aa_arm_opening, link.id);
            arm_opening = actn2_arm_kinetics::decide_bridge_opening(
                rate_i, rate_j, dt, [&]() { return opening_rng.uniform(); });
        } else if (!detachment.detaches()) {
            detachment = catch_bond_kinetics::decide_detachment(
                mechanics[slot].length, aa_link_rupture_length, dt, lifetime,
                [&]() { return detachment_rng.uniform(); });
        }
        const bool contour_rupture = detachment.geometric;
        const bool overlap_occlusion = detachment.am_overlap_occluded;
        const bool collision_rupture = detachment.actn2_myosin_collision;
        const bool geometric_rupture =
            contour_rupture || overlap_occlusion || collision_rupture ||
            detachment.mature_connector_excluded;
        const bool arm_opens = arm_opening.opens_i || arm_opening.opens_j;
        const bool detaches = detachment.detaches();
        // Everything from here mutates shared, order-sensitive state: the
        // event streams are appended in slot order, aa_neighbor_slots is a
        // single hash map, and the free-slot stack decides which slot the next
        // association reuses.  `ordered` replays this block in the original
        // sequential order regardless of how the iterations were scheduled.
        {
        auto record_motor_pair_event = [&actn2_load_kinetics_unused =
                                             motor_transmitted_mode,
                                         this, &link,
                                         motor_pair_current_support,
                                         motor_pair_best_support,
                                         motor_pair_candidate_age](
            int event_code,
            const actn2_load_kinetics::MotorPairIdentity& old_identity,
            const actn2_load_kinetics::MotorPairIdentity& new_identity) {
            if (!actn2_load_kinetics_unused) return;
            aa_motor_pair_events.insert(aa_motor_pair_events.end(), {
                static_cast<double>(link.id >> 32),
                static_cast<double>(link.id & 0xffffffffULL),
                static_cast<double>(current_step),
                static_cast<double>(event_code),
                static_cast<double>(old_identity.myosin_i),
                static_cast<double>(old_identity.myosin_j),
                static_cast<double>(old_identity.half_i),
                static_cast<double>(old_identity.half_j),
                static_cast<double>(new_identity.myosin_i),
                static_cast<double>(new_identity.myosin_j),
                static_cast<double>(new_identity.half_i),
                static_cast<double>(new_identity.half_j),
                motor_pair_current_support, motor_pair_best_support,
                link.candidate_support_ema, motor_pair_candidate_age,
                link.candidate_raw_support_i,
                link.candidate_raw_support_j,
                link.candidate_allocated_support_i,
                link.candidate_allocated_support_j});
        };
        const actn2_load_kinetics::MotorPairIdentity new_candidate_identity{
            link.candidate_myosin_i, link.candidate_myosin_j,
            link.candidate_half_i, link.candidate_half_j};
        if (motor_pair_event !=
            actn2_load_kinetics::MotorPairEvent::none) {
            record_motor_pair_event(
                static_cast<int>(motor_pair_event), old_candidate_identity,
                new_candidate_identity);
        }
        if (!candidate_was_locked && link.candidate_locked) {
            record_motor_pair_event(4, new_candidate_identity,
                                    new_candidate_identity);
        }
        if (status2_unlocked && candidate_was_locked &&
            new_candidate_identity.valid()) {
            record_motor_pair_event(6, new_candidate_identity,
                                    new_candidate_identity);
        }
        if (link.status != old_status) {
            aa_local_link_events.insert(aa_local_link_events.end(), {
                static_cast<double>(link.id),
                static_cast<double>(link.actin_i),
                static_cast<double>(link.actin_j),
                static_cast<double>(current_step), 3.0,
                static_cast<double>(old_status),
                static_cast<double>(link.status), effective_load, lifetime,
                link.s_i, link.s_j});
            aa_local_link_events_v2.insert(aa_local_link_events_v2.end(), {
                static_cast<double>(link.id), static_cast<double>(link.actin_i),
                static_cast<double>(link.actin_j), static_cast<double>(current_step),
                3.0, static_cast<double>(old_status),
                static_cast<double>(link.status), effective_load, lifetime,
                link.s_i, link.s_j, kinetic_arm_loads.first,
                kinetic_arm_loads.second,
                kinetic_arm_loads.favorable_bilateral(), productive_support,
                instantaneous_status2_gate ? 1.0 : 0.0, link_tension,
                arm_loads.first, arm_loads.second,
                allocated_motor_i, allocated_motor_j});
            aa_local_link_events_v3.insert(aa_local_link_events_v3.end(), {
                static_cast<double>(link.id), static_cast<double>(link.actin_i),
                static_cast<double>(link.actin_j), static_cast<double>(current_step),
                3.0, static_cast<double>(old_status),
                static_cast<double>(link.status), effective_load, lifetime,
                link.s_i, link.s_j, kinetic_arm_loads.first,
                kinetic_arm_loads.second,
                kinetic_arm_loads.favorable_bilateral(), productive_support,
                instantaneous_status2_gate ? 1.0 : 0.0, link_tension,
                arm_loads.first, arm_loads.second,
                allocated_motor_i, allocated_motor_j,
                mechanics[slot].length, aa_link_rupture_length,
                mechanics[slot].direction_dot, mechanics[slot].registry_error,
                registry_threshold, registry_ok ? 1.0 : 0.0,
                arm_load_ok ? 1.0 : 0.0,
                productive_support_ok ? 1.0 : 0.0,
                instantaneous_status2_gate ? 1.0 : 0.0,
                0.0, detachment.probability, detachment.draw, 0.0, 0.0});
        }
        if (aa_neighbor_identity_enabled && !arm_opens) {
            aa_neighbor_slots.observe(
                link.actin_i, link.actin_j, static_cast<int>(link.status),
                static_cast<aa_neighbor_identity::Step>(current_step),
                static_cast<aa_neighbor_identity::Step>(
                    aa_status2_slot_memory_steps));
        }
        if (arm_opens && !detaches) {
            const int old_i = link.actin_i;
            const int old_j = link.actin_j;
            const int bound_actin = arm_opening.opens_i ? old_j : old_i;
            const int opened_partner = arm_opening.opens_i ? old_i : old_j;
            const double bound_s = arm_opening.opens_i ? link.s_j : link.s_i;
            const double opened_s = arm_opening.opens_i ? link.s_i : link.s_j;
            link.last_partner = opened_partner;
            link.status = 1;
            link.arm_state = arm_opening.opens_i
                ? actn2_arm_kinetics::State::tethered_to_j
                : actn2_arm_kinetics::State::tethered_to_i;
            link.arm_state_since = static_cast<RecoveryStep>(std::min(
                current_step,
                static_cast<size_t>(std::numeric_limits<RecoveryStep>::max())));
            if (new_candidate_identity.valid()) {
                record_motor_pair_event(
                    5, new_candidate_identity,
                    actn2_load_kinetics::MotorPairIdentity{});
            }
            link.candidate_locked = false;
            link.candidate_myosin_i = -1;
            link.candidate_myosin_j = -1;
            link.challenger_myosin_i = -1;
            link.challenger_myosin_j = -1;
            link.candidate_support_ema = 0.0;
            link.challenger_support_ema = 0.0;
            if (aa_neighbor_identity_enabled) {
                aa_neighbor_slots.erase(old_i, old_j);
            }
            aa_arm_events.insert(aa_arm_events.end(), {
                static_cast<double>(link.id),
                static_cast<double>(current_step), 1.0,
                static_cast<double>(static_cast<int>(old_arm_state)),
                static_cast<double>(static_cast<int>(link.arm_state)),
                static_cast<double>(bound_actin),
                static_cast<double>(opened_partner), -1.0, bound_s, opened_s,
                arm_opening.opens_i ? tensile_arm_i : tensile_arm_j,
                link.contact_force, link.contact_penetration,
                arm_opening.opens_i ? arm_opening.rate_i : arm_opening.rate_j});
        }
        if (detaches) {
            // Code 5 means that an ACTN2 detachment actually released a held
            // motor-pair identity.  Do not emit a release row for links that
            // never selected a candidate; those are ordinary ACTN2
            // detachments and otherwise inflate the partner-release count.
            if (new_candidate_identity.valid()) {
                record_motor_pair_event(
                    5, new_candidate_identity,
                    actn2_load_kinetics::MotorPairIdentity{});
            }
            if (geometric_rupture && aa_neighbor_identity_enabled) {
                aa_neighbor_slots.erase(link.actin_i, link.actin_j);
            }
            if (link.status == 2) {
                // Slot 0 always exists and, unlike omp_get_thread_num() inside
                // an `ordered` region, does not depend on which thread happened
                // to reach this iteration.  At one thread this is the index the
                // routine already used.
                aa_completed_lifetimes_temp[0].push_back(
                    (current_step - static_cast<size_t>(link.attach_step)) * dt);
            }
            aa_local_link_events.insert(aa_local_link_events.end(), {
                static_cast<double>(link.id),
                static_cast<double>(link.actin_i),
                static_cast<double>(link.actin_j),
                static_cast<double>(current_step), 2.0,
                static_cast<double>(link.status), 0.0, effective_load,
                lifetime, link.s_i, link.s_j});
            aa_local_link_events_v2.insert(aa_local_link_events_v2.end(), {
                static_cast<double>(link.id), static_cast<double>(link.actin_i),
                static_cast<double>(link.actin_j), static_cast<double>(current_step),
                2.0, static_cast<double>(link.status), 0.0, effective_load,
                lifetime, link.s_i, link.s_j, kinetic_arm_loads.first,
                kinetic_arm_loads.second,
                kinetic_arm_loads.favorable_bilateral(), productive_support,
                instantaneous_status2_gate ? 1.0 : 0.0, link_tension,
                arm_loads.first, arm_loads.second,
                allocated_motor_i, allocated_motor_j});
            aa_local_link_events_v3.insert(aa_local_link_events_v3.end(), {
                static_cast<double>(link.id), static_cast<double>(link.actin_i),
                static_cast<double>(link.actin_j), static_cast<double>(current_step),
                2.0, static_cast<double>(link.status), 0.0, effective_load,
                lifetime, link.s_i, link.s_j, kinetic_arm_loads.first,
                kinetic_arm_loads.second,
                kinetic_arm_loads.favorable_bilateral(), productive_support,
                instantaneous_status2_gate ? 1.0 : 0.0, link_tension,
                arm_loads.first, arm_loads.second,
                allocated_motor_i, allocated_motor_j,
                mechanics[slot].length, aa_link_rupture_length,
                mechanics[slot].direction_dot, mechanics[slot].registry_error,
                registry_threshold, registry_ok ? 1.0 : 0.0,
                arm_load_ok ? 1.0 : 0.0,
                productive_support_ok ? 1.0 : 0.0,
                instantaneous_status2_gate ? 1.0 : 0.0,
                static_cast<double>(static_cast<int>(detachment.reason())),
                detachment.probability, detachment.draw,
                geometric_rupture ? 1.0 : 0.0,
                detachment.stochastic ? 1.0 : 0.0});
            link.active = false;
            aa_free_link_slots.push_back(slot);
        }
        }  // formerly omp ordered
    }
    }  // omp single
    if (motor_transmitted_mode) {
        #pragma omp atomic update
        aa_motor_candidate_seconds += motor_candidate_thread_seconds;
        #pragma omp atomic update
        aa_motor_candidate_pairs_evaluated += motor_candidate_thread_pairs;
    }

    stamp_aa_phase(kAaReindex);
    #pragma omp single
    {
    _rebuild_local_aa_indices();
    if (aa_neighbor_identity_enabled) {
        aa_neighbor_slots.expire_unseen(
            static_cast<aa_neighbor_identity::Step>(current_step),
            static_cast<aa_neighbor_identity::Step>(
                aa_neighbor_slot_memory_steps));
    }
    if (aa_scratch_sorted_neighbors.size() != static_cast<size_t>(actin.n)) {
        aa_scratch_sorted_neighbors.assign(static_cast<size_t>(actin.n), {});
    }
    aa_scratch_pair_offset.assign(static_cast<size_t>(actin.n) + 1, 0);
    aa_scratch_dirty_actin.assign(static_cast<size_t>(actin.n), 0);
    }

    // Canonical proposal order: outer actin ascending, partner ascending,
    // restricted to i > j.  Materializing that list lets the (pure) evaluation
    // of every pair run across the team before the serial commit.  Building it
    // is itself O(N log k) and was the largest remaining serial block -- worse,
    // one thread rebuilding a list the whole team had just dirtied got *slower*
    // as the team grew.  It is now a per-actin sort in parallel, an O(N) prefix
    // sum, and a parallel scatter into exact offsets, which reproduces the same
    // ordering the serial loop produced.
    stamp_aa_phase(kAaSortNeighbors);
    #pragma omp for schedule(static)
    for (int i = 0; i < actin.n; ++i) {
        auto& sorted = aa_scratch_sorted_neighbors[static_cast<size_t>(i)];
        // Filter before sorting.  The surviving list is a subsequence of the
        // old one in the same ascending order, so the pair list, its ordering,
        // and therefore the greedy commit order below are all unchanged; only
        // pairs that provably cannot fire are absent.
        sorted.clear();
        // Source from the segment-distance-filtered list (~539 entries) rather
        // than the raw centre-based one (~890).  Its cutoff is 0.77 um and a
        // viable ACTN2 pair needs segment distance below aa_cutoff = 0.12 um,
        // so nothing that could associate is missing from it.
        if (opt_gates_enabled) {
            for (int j : actin_close_neighbors[static_cast<size_t>(i)]) {
                if (_aa_pair_may_associate(i, j)) sorted.push_back(j);
            }
        } else {
            sorted.assign(actin_neighbors_by_species[i].first.begin(),
                          actin_neighbors_by_species[i].first.end());
        }
        std::sort(sorted.begin(), sorted.end());
        size_t count = 0;
        for (int j : sorted) if (i > j) ++count;
        aa_scratch_pair_offset[static_cast<size_t>(i) + 1] = count;
    }

    stamp_aa_phase(kAaPrefixSum);
    #pragma omp single
    {
    for (int i = 0; i < actin.n; ++i) {
        aa_scratch_pair_offset[static_cast<size_t>(i) + 1] +=
            aa_scratch_pair_offset[static_cast<size_t>(i)];
    }
    const size_t total_pairs = aa_scratch_pair_offset[
        static_cast<size_t>(actin.n)];
    // resize(), not assign(): every entry is overwritten by the evaluation
    // loop below, so re-initializing the whole buffer each step is wasted.
    if (aa_scratch_pairs.size() < total_pairs) {
        aa_scratch_pairs.resize(total_pairs);
        aa_scratch_proposals.resize(total_pairs);
    }
    aa_scratch_pair_count = total_pairs;
    }

    stamp_aa_phase(kAaScatterPairs);
    #pragma omp for schedule(static)
    for (int i = 0; i < actin.n; ++i) {
        size_t at = aa_scratch_pair_offset[static_cast<size_t>(i)];
        for (int j : aa_scratch_sorted_neighbors[static_cast<size_t>(i)]) {
            if (i <= j) continue;
            aa_scratch_pairs[at++] = {i, j};
        }
    }

    // Phase 3a.  Evaluate every candidate against the post-detachment state.
    // _evaluate_aa_association() is const, so the whole team can run it.
    const long long n_pairs = static_cast<long long>(aa_scratch_pair_count);
    stamp_aa_phase(kAaEvalProposals);
    #pragma omp for schedule(static)
    for (long long k = 0; k < n_pairs; ++k) {
        aa_scratch_proposals[static_cast<size_t>(k)] = _evaluate_aa_association(
            aa_scratch_pairs[static_cast<size_t>(k)].first,
            aa_scratch_pairs[static_cast<size_t>(k)].second);
    }

    stamp_aa_phase(kAaCommitLinks);
    #pragma omp single
    {
    auto add_link = [&](int i, int j, const AAAssociationProposal& proposal) {
        const double s_i = proposal.s_i;
        const double s_j = proposal.s_j;
        if (aa_zdisc_tip_contact && proposal.nearest_tip_myosin >= 0 &&
            proposal.nearest_tip_gap < aa_zdisc_tip_onset - 1e-12) {
            throw std::runtime_error(
                "ACTN2 association attempted to commit inside a myosin-tip "
                "exclusion zone");
        }
        AACrosslink link;
        link.id = aa_next_link_id++;
        link.actin_i = i;
        link.actin_j = j;
        link.s_i = s_i;
        link.s_j = s_j;
        link.status = 1;
        link.arm_state = actn2_arm_kinetics::State::weak_bridge;
        link.attach_step = static_cast<RecoveryStep>(std::min(
            current_step,
            static_cast<size_t>(std::numeric_limits<RecoveryStep>::max())));
        link.arm_state_since = link.attach_step;
        link.active = true;
        aa_tip_accessibility_events.insert(
            aa_tip_accessibility_events.end(), {
                static_cast<double>(link.id),
                static_cast<double>(current_step), 1.0,
                static_cast<double>(i), static_cast<double>(j), s_i, s_j,
                static_cast<double>(proposal.nearest_tip_myosin),
                static_cast<double>(proposal.nearest_tip_sign),
                proposal.nearest_tip_gap, proposal.nearest_tip_rho, 0.0,
                proposal.nearest_tip_motor_engaged ? 1.0 : 0.0});
        aa_local_link_events.insert(aa_local_link_events.end(), {
            static_cast<double>(link.id), static_cast<double>(i),
            static_cast<double>(j), static_cast<double>(current_step),
            1.0, 0.0, 1.0, 0.0, base_lifetime, s_i, s_j});
        aa_local_link_events_v2.insert(aa_local_link_events_v2.end(), {
            static_cast<double>(link.id), static_cast<double>(i),
            static_cast<double>(j), static_cast<double>(current_step),
            1.0, 0.0, 1.0, 0.0, base_lifetime, s_i, s_j,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0});
        size_t slot = 0;
        if (!aa_free_link_slots.empty()) {
            slot = aa_free_link_slots.back();
            aa_free_link_slots.pop_back();
            aa_links[slot] = link;
        } else {
            slot = aa_links.size();
            aa_links.push_back(link);
            mechanics.emplace_back();
            contact_results.emplace_back();
            contact_evaluated.push_back(0);
        }
        contact_results[slot].clear();
        contact_evaluated[slot] = 0;
        mechanics[slot] = compute_local_aa_link(
            actin, i, j, s_i, s_j, box,
            ramp * aa_link_radial_stiffness(link.status),
            ramp * aa_link_perpendicular_stiffness(link.status),
            single_link_alignment_k, aa_optimal, aa_antiparallel_overlap,
            aa_link_register_ratio(link.status),
            aa_register_saturation_error, aa_axial_register_force_max,
            ramp * aa_radial_force_max);
        const auto& initial_mechanics = mechanics[slot];
        const bool initial_registry_ok = actn2_load_kinetics::registry_gate(
            1, initial_mechanics.direction_dot, initial_mechanics.registry_error,
            aa_status2_registry_on, aa_status2_registry_off);
        const double diagnostic_nan =
            std::numeric_limits<double>::quiet_NaN();
        aa_local_link_events_v3.insert(aa_local_link_events_v3.end(), {
            static_cast<double>(link.id), static_cast<double>(i),
            static_cast<double>(j), static_cast<double>(current_step),
            1.0, 0.0, 1.0, 0.0, base_lifetime, s_i, s_j,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            initial_mechanics.length, aa_link_rupture_length,
            initial_mechanics.direction_dot, initial_mechanics.registry_error,
            aa_status2_registry_on, initial_registry_ok ? 1.0 : 0.0,
            0.0, 0.0, 0.0,
            0.0, diagnostic_nan, diagnostic_nan, 0.0, 0.0});
        // Maintain sparse indices incrementally. Rebuilding every index after
        // one local KMC association made a successful event O(Nlinks).
        aa_active_link_slots.push_back(slot);
        aa_links_by_pair[_aa_pair_key(i, j)].push_back(slot);
        aa_links_by_actin[static_cast<size_t>(i)].push_back(slot);
        aa_links_by_actin[static_cast<size_t>(j)].push_back(slot);
    };

    // Association proposals are committed in canonical actin order, exactly as
    // the original single serial loop did.  A proposal only depends on state
    // that add_link()/reserve() touch at the two actins of the committed pair
    // (pair occupancy, occupied axial sites, neighbor degree), so a cached
    // verdict is valid unless one of its actins has already been committed
    // this step; those few pairs are re-evaluated here.  The global link pool,
    // when enabled, makes every rate depend on the running link count, so that
    // configuration re-evaluates unconditionally.
    const bool pool_limited = aa_total_link_pool > 0;
    for (size_t k = 0; k < aa_scratch_pair_count; ++k) {
        const int i = aa_scratch_pairs[k].first;
        const int j = aa_scratch_pairs[k].second;
        const bool stale = pool_limited ||
            aa_scratch_dirty_actin[static_cast<size_t>(i)] ||
            aa_scratch_dirty_actin[static_cast<size_t>(j)];
        const AAAssociationProposal proposal = stale
            ? _evaluate_aa_association(i, j)
            : aa_scratch_proposals[k];
        if (!proposal.fire) continue;
        if (aa_neighbor_identity_enabled &&
            !aa_neighbor_slots.reserve(
                i, j,
                static_cast<aa_neighbor_identity::Step>(current_step),
                max_aa_bonds_per_actin)) {
            continue;
        }
        add_link(i, j, proposal);
        aa_scratch_dirty_actin[static_cast<size_t>(i)] = 1;
        aa_scratch_dirty_actin[static_cast<size_t>(j)] = 1;
    }
    }  // omp single

    // Apply mechanics and create the old pair summaries as derived state.
    // This phase does no filament mechanics; it is a scatter of already
    // computed quantities into shared per-actin accumulators and the dense
    // pair matrices, so it stays serialized rather than paying for privatized
    // copies of those matrices.
    stamp_aa_phase(kAaTipContact);
    // ---- parallel tip-contact scan -------------------------------------
    // Read-only: patch geometry, the conservative sphere test and the contact
    // kernel.  Survivors go to per-link buckets; the order-sensitive reduction
    // runs serially below over the same link sequence.
    if (aa_zdisc_tip_contact) {
        const size_t n_slots = aa_active_link_slots.size();
        #pragma omp single
        {
        if (aa_tip_scan_buckets.size() < n_slots) {
            aa_tip_scan_buckets.resize(n_slots);
        }
        }
        #pragma omp for schedule(dynamic, 8)
        for (long long index = 0; index < static_cast<long long>(n_slots);
             ++index) {
            auto& bucket = aa_tip_scan_buckets[static_cast<size_t>(index)];
            bucket.clear();
            const size_t slot = aa_active_link_slots[static_cast<size_t>(index)];
            const auto& link = aa_links[slot];
            if (aa_arm_state_model &&
                !actn2_arm_kinetics::is_bridge(link.arm_state)) continue;
            const auto patch = compute_actn2_patch_geometry(
                actin, link.actin_i, link.actin_j, link.s_i, link.s_j, box);
            // A connector interior may contact a myosin absent from both
            // endpoint actins' neighbor lists. Use a conservative segment
            // sphere against every myosin before the exact contact query.
            for (int m = 0; m < myosin.n; ++m) {
                if (!actn2_patch_may_contact_myosin_tip(
                        patch, myosin, m, box, aa_zdisc_tip_onset,
                        aa_zdisc_tip_stop_gate_cutoff)) continue;
                for (int end_sign : {-1, 1}) {
                    auto contact = compute_zdisc_axial_stop_contact(
                        actin, myosin, link.actin_i, link.actin_j,
                        link.s_i, link.s_j, m, end_sign, box,
                        aa_zdisc_tip_onset, aa_zdisc_tip_stop_stiffness,
                        aa_zdisc_tip_stop_force_max,
                        aa_zdisc_tip_stop_gate_radius,
                        aa_zdisc_tip_stop_gate_cutoff,
                        aa_zdisc_tip_stop_epsilon);
                    if (!contact.applied) continue;
                    bucket.push_back({slot, m, end_sign, contact});
                }
            }
        }
    }
    #pragma omp single
    {
    const int thread_id = omp_get_thread_num();
    auto& local_forces = actin_forces_temp[thread_id];
    auto& local_torques = actin_torques_temp[thread_id];
    auto& local_myosin_forces = myosin_forces_temp[thread_id];
    auto& local_myosin_torques = myosin_torques_temp[thread_id];
    auto& local_status = actin_cb_status_temp[thread_id];
    auto& aa_local_link_count_feature = actin["aa_local_link_count"];
    auto& aa_local_status1_count_feature = actin["aa_local_status1_count"];
    auto& aa_local_status2_count_feature = actin["aa_local_status2_count"];
    auto& aa_local_tension_max_feature = actin["aa_local_tension_max"];
    auto& aa_zdisc_tip_force_feature = actin["aa_zdisc_tip_force"];
    auto& aa_zdisc_tip_contact_count_feature =
        actin["aa_zdisc_tip_contact_count"];
    if (aa_zdisc_tip_contact) {
        // One energy budget per (unordered actin pair, myosin): redundant
        // links on a pair contribute only the link with the greatest energy
        // summed over both tips, ties going to the lower link id.  Every
        // distinct pair contributes.  The scan is replayed in canonical link
        // order, so selection and application order are deterministic.
        struct PairMyosinSelection {
            size_t bucket = 0;      // index into aa_tip_scan_buckets
            size_t first_hit = 0;   // this myosin's hits are contiguous
            size_t hit_count = 0;
            double energy = 0.0;
            std::uint64_t link_id = 0;
        };
        struct PairMyosinKey {
            std::uint64_t pair;
            int myosin_index;
            bool operator==(const PairMyosinKey& other) const {
                return pair == other.pair && myosin_index == other.myosin_index;
            }
        };
        struct PairMyosinHash {
            size_t operator()(const PairMyosinKey& key) const {
                return std::hash<std::uint64_t>{}(
                    key.pair * 0x9e3779b97f4a7c15ULL ^
                    static_cast<std::uint64_t>(key.myosin_index));
            }
        };
        std::vector<PairMyosinSelection> selections;
        std::unordered_map<PairMyosinKey, size_t, PairMyosinHash> selection_of;
        for (size_t index = 0; index < aa_active_link_slots.size(); ++index) {
            const auto& hits = aa_tip_scan_buckets[index];
            for (size_t begin = 0; begin < hits.size();) {
                const int m = hits[begin].myosin_index;
                size_t end = begin;
                double energy = 0.0;
                while (end < hits.size() && hits[end].myosin_index == m) {
                    energy += hits[end].contact.energy;
                    ++end;
                }
                const auto& link = aa_links[hits[begin].slot];
                const PairMyosinSelection candidate{
                    index, begin, end - begin, energy, link.id};
                const PairMyosinKey key{
                    _aa_pair_key(link.actin_i, link.actin_j), m};
                const auto found = selection_of.find(key);
                if (found == selection_of.end()) {
                    selection_of.emplace(key, selections.size());
                    selections.push_back(candidate);
                } else {
                    auto& best = selections[found->second];
                    if (energy > best.energy + 1e-12 ||
                        (std::abs(energy - best.energy) <= 1e-12 &&
                         link.id < best.link_id)) {
                        best = candidate;
                    }
                }
                begin = end;
            }
        }
        for (const auto& selection : selections) {
            const auto& hits = aa_tip_scan_buckets[selection.bucket];
            for (size_t h = selection.first_hit;
                 h < selection.first_hit + selection.hit_count; ++h) {
                const auto& hit = hits[h];
                const auto& link = aa_links[hit.slot];
                const auto& contact = hit.contact;
                local_forces[link.actin_i] += contact.force_on_actin_i;
                local_forces[link.actin_j] += contact.force_on_actin_j;
                local_myosin_forces[hit.myosin_index] += contact.force_on_myosin;
                local_torques[link.actin_i] += contact.torque_on_actin_i;
                local_torques[link.actin_j] += contact.torque_on_actin_j;
                local_myosin_torques[hit.myosin_index] += contact.torque_on_myosin;
                _tb_add(kTsActn2Contact, thread_id, link.actin_i,
                        contact.torque_on_actin_i);
                _tb_add(kTsActn2Contact, thread_id, link.actin_j,
                        contact.torque_on_actin_j);
                aa_zdisc_tip_force_feature[link.actin_i] +=
                    contact.force_on_actin_i.norm();
                aa_zdisc_tip_force_feature[link.actin_j] +=
                    contact.force_on_actin_j.norm();
                aa_zdisc_tip_contact_count_feature[link.actin_i] += 1.0;
                aa_zdisc_tip_contact_count_feature[link.actin_j] += 1.0;
                // Per tip: total axial force, and the deepest contact point.
                const bool plus = hit.end_sign > 0;
                myosin[plus ? "aa_zdisc_tip_force_plus" :
                               "aa_zdisc_tip_force_minus"][hit.myosin_index] +=
                    contact.axial_force;
                auto& gap_feature = myosin[plus ? "aa_zdisc_tip_gap_plus" :
                                                  "aa_zdisc_tip_gap_minus"];
                const double deepest = gap_feature[hit.myosin_index];
                if (std::isnan(deepest) || contact.deepest_gap < deepest) {
                    gap_feature[hit.myosin_index] = contact.deepest_gap;
                    myosin[plus ? "aa_zdisc_tip_rho_plus" :
                                   "aa_zdisc_tip_rho_minus"][hit.myosin_index] =
                        contact.deepest_rho;
                    myosin[plus ? "aa_zdisc_tip_gate_plus" :
                                   "aa_zdisc_tip_gate_minus"][hit.myosin_index] =
                        contact.deepest_gate;
                }
            }
        }
    }
    stamp_aa_phase(kAaTipAndScatter);
    for (size_t slot : aa_active_link_slots) {
        const auto& link = aa_links[slot];
        if (aa_arm_state_model &&
            !actn2_arm_kinetics::is_bridge(link.arm_state)) continue;
        const auto& mech = mechanics[slot];
        local_forces[link.actin_i] += mech.force_on_first;
        local_forces[link.actin_j] -= mech.force_on_first;
        local_torques[link.actin_i] += mech.torque_on_first;
        local_torques[link.actin_j] += mech.torque_on_second;
        // Split by which ACTN2 stiffness supplies it.  The four terms sum to
        // torque_on_*, which the closure check asserts.
        _tb_add(kTsActn2Lever, thread_id, link.actin_i, mech.torque_first_lever);
        _tb_add(kTsActn2Perp,  thread_id, link.actin_i, mech.torque_first_perp);
        _tb_add(kTsActn2Align, thread_id, link.actin_i, mech.torque_first_align);
        _tb_add(kTsActn2Axial, thread_id, link.actin_i, mech.torque_first_axial);
        _tb_add(kTsActn2Lever, thread_id, link.actin_j, mech.torque_second_lever);
        _tb_add(kTsActn2Perp,  thread_id, link.actin_j, mech.torque_second_perp);
        _tb_add(kTsActn2Align, thread_id, link.actin_j, mech.torque_second_align);
        _tb_add(kTsActn2Axial, thread_id, link.actin_j, mech.torque_second_axial);
        if ((aa_myosin_contact_epsilon > 0.0 ||
             (aa_myosin_contact_shell_width > 0.0 &&
              aa_myosin_contact_shell_energy > 0.0)) &&
            aa_actn2_contact_radius > 0.0) {
            auto apply_contact = [&](
                int myosin_index,
                const ACTN2MyosinContactResult& contact) {
                local_forces[link.actin_i] += contact.force_on_actin_i;
                local_forces[link.actin_j] += contact.force_on_actin_j;
                local_myosin_forces[myosin_index] += contact.force_on_myosin;
                local_torques[link.actin_i] += contact.torque_on_actin_i;
                local_torques[link.actin_j] += contact.torque_on_actin_j;
                local_myosin_torques[myosin_index] += contact.torque_on_myosin;
                _tb_add(kTsActn2Contact, thread_id, link.actin_i,
                        contact.torque_on_actin_i);
                _tb_add(kTsActn2Contact, thread_id, link.actin_j,
                        contact.torque_on_actin_j);
            };
            if (slot < contact_evaluated.size() &&
                contact_evaluated[slot]) {
                for (const auto& entry : contact_results[slot]) {
                    apply_contact(entry.first, entry.second);
                }
            } else {
                // A bridge associated after phase 1 has no cached contact yet.
                // This is the only path that evaluates contact during force
                // application; pre-existing bridges reuse the kinetic pass.
                std::vector<int> contact_myosins =
                    actin_neighbors_by_species[
                        static_cast<size_t>(link.actin_i)].second;
                const auto& second_neighbors =
                    actin_neighbors_by_species[
                        static_cast<size_t>(link.actin_j)].second;
                contact_myosins.insert(
                    contact_myosins.end(), second_neighbors.begin(),
                    second_neighbors.end());
                std::sort(contact_myosins.begin(), contact_myosins.end());
                contact_myosins.erase(
                    std::unique(contact_myosins.begin(), contact_myosins.end()),
                    contact_myosins.end());
                const auto patch = compute_actn2_patch_geometry(
                    actin, link.actin_i, link.actin_j,
                    link.s_i, link.s_j, box);
                const double outer_contact_distance =
                    myosin.radius + aa_actn2_contact_radius +
                    aa_myosin_contact_shell_width;
                for (int myosin_index : contact_myosins) {
                    if (!actn2_patch_may_contact_myosin_body(
                            patch, myosin, myosin_index, box,
                            outer_contact_distance)) continue;
                    const auto contact = compute_actn2_myosin_contact(
                        actin, myosin, link.actin_i, link.actin_j,
                        link.s_i, link.s_j, myosin_index, box,
                        myosin.radius + aa_actn2_contact_radius,
                        aa_myosin_contact_epsilon,
                        aa_myosin_contact_force_cap,
                        aa_myosin_contact_shell_width,
                        aa_myosin_contact_shell_energy);
                    if (contact.applied) apply_contact(myosin_index, contact);
                }
            }
        }
        local_status[link.actin_i] = std::max(local_status[link.actin_i],
                                              static_cast<int>(link.status));
        local_status[link.actin_j] = std::max(local_status[link.actin_j],
                                              static_cast<int>(link.status));
        ++actin_n_bonds[link.actin_i];
        ++actin_n_bonds[link.actin_j];
        aa_local_link_count_feature[link.actin_i] += 1.0;
        aa_local_link_count_feature[link.actin_j] += 1.0;
        auto& status_feature = link.status == 2
            ? aa_local_status2_count_feature : aa_local_status1_count_feature;
        status_feature[link.actin_i] += 1.0;
        status_feature[link.actin_j] += 1.0;
        const double link_tension = std::max(
            0.0, ramp * aa_link_radial_stiffness(link.status) *
                (mech.length - aa_optimal));
        aa_local_tension_max_feature[link.actin_i] = std::max(
            aa_local_tension_max_feature[link.actin_i], link_tension);
        aa_local_tension_max_feature[link.actin_j] = std::max(
            aa_local_tension_max_feature[link.actin_j], link_tension);
        if (link.status == 2) {
            ++actin_strong_cb_count[link.actin_i];
            ++actin_strong_cb_count[link.actin_j];
        }
        const int first = std::min(link.actin_i, link.actin_j);
        const int second = std::max(link.actin_i, link.actin_j);
        actin_actin_bonds[first][second] = actin_actin_bonds[second][first] = 1;
        actin_actin_status[first][second] = actin_actin_status[second][first] =
            std::max(actin_actin_status[first][second], link.status);
        actin_actin_lifetime[first][second] = actin_actin_lifetime[second][first] =
            std::max(actin_actin_lifetime[first][second],
                static_cast<int>(current_step - link.attach_step + 1));
    active_aa_pairs_temp[thread_id].emplace_back(first, second);
    }
    if (motor_transmitted_mode) {
        aa_motor_local_update_seconds += omp_get_wtime() -
            aa_motor_update_started;
        ++aa_motor_profile_steps;
    }
    }  // omp single
    stamp_aa_phase(0);   // close the final phase
}

void Sarcomere::_process_catch_bonds(int& i) {
    std::vector <int> actin_neighbors = actin_neighbors_by_species[i].first;
    // A capture neighbor list is sufficient for new associations, but an
    // occupied spring must continue to be evaluated after it stretches beyond
    // capture range.  KMC, not neighbor-list loss, owns detachment.
    for (const auto& pair : previous_aa_pairs) {
        if (pair.first == i) actin_neighbors.push_back(pair.second);
        if (pair.second == i) actin_neighbors.push_back(pair.first);
    }
    std::sort(actin_neighbors.begin(), actin_neighbors.end());
    actin_neighbors.erase(
        std::unique(actin_neighbors.begin(), actin_neighbors.end()),
        actin_neighbors.end());
    std::vector<int> statuses;
    std::vector<int> cb_indices;
    for (int index = 0; index < actin_neighbors.size(); index++){
            int j = actin_neighbors[index];
            if (i>j){
                int status = determine_cb_status(i,j);
                if (status > 0){
                    bool acc = _cb_decide(i,j,status);
                    if (acc){
                        cb_indices.push_back(j);
                        statuses.push_back(status);
                    }
                }
            }
    }
    _set_cb(i, cb_indices,statuses);
}

void Sarcomere::_record_aa_state_changes() {
    // The candidate list is rebuilt from scratch every step but its size barely
    // moves, so a persistent buffer keeps the capacity instead of reallocating
    // (and copying previous_aa_pairs into a fresh allocation) each time.
    std::vector<BondPair>& candidates = aa_state_change_candidates;
    candidates.assign(previous_aa_pairs.begin(), previous_aa_pairs.end());
    for (const auto& local_pairs : active_aa_pairs_temp) {
        for (auto pair : local_pairs) {
            if (pair.first > pair.second) std::swap(pair.first, pair.second);
            candidates.push_back(pair);
        }
    }
    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());

    const int thread_id = omp_in_parallel() ? omp_get_thread_num() : 0;
    auto& events = aa_state_events_temp[thread_id];
    auto& am_motor_force_feature = actin["am_motor_force"];
    auto& am_strong_count_feature = actin["am_strong_count"];
    for (const auto& [i, j] : candidates) {
        const int old_bond = static_cast<int>(actin_actin_bonds_prev[i][j]);
        const int new_bond = static_cast<int>(actin_actin_bonds[i][j]);
        const int old_status = old_bond
            ? static_cast<int>(actin_actin_status_prev[i][j]) : 0;
        const int new_status = new_bond
            ? static_cast<int>(actin_actin_status[i][j]) : 0;
        if (old_bond == new_bond && old_status == new_status) continue;

        const double load_i = actin_f_load_cb ? (*actin_f_load_cb)[i] : actin.f_load[i];
        const double load_j = actin_f_load_cb ? (*actin_f_load_cb)[j] : actin.f_load[j];
        events.insert(events.end(), {
            static_cast<double>(i),
            static_cast<double>(j),
            static_cast<double>(current_step),
            static_cast<double>(old_bond),
            static_cast<double>(new_bond),
            static_cast<double>(old_status),
            static_cast<double>(new_status),
            am_motor_force_feature[i],
            am_motor_force_feature[j],
            am_strong_count_feature[i],
            am_strong_count_feature[j],
            load_i,
            load_j});

        const auto force_i = _actin_motor_force_vector(i);
        const auto force_j = _actin_motor_force_vector(j);
        const double motor_only = actn2_load_kinetics::catch_load(
            actn2_load_kinetics::StatusMode::motor_only,
            force_i, force_j, 0.0);
        const double opposed_drive =
            actn2_load_kinetics::opposed_motor_drive(force_i, force_j);
        const double transmitted = _aa_transmitted_load(i, j);
        const double combined_candidate = std::min(opposed_drive, transmitted);
        const double selected =
            aa_status_mode == actn2_load_kinetics::StatusMode::force_projection
            ? _aa_effective_arm_loads(i, j).favorable_bilateral()
            : actn2_load_kinetics::catch_load(
                  aa_status_mode, force_i, force_j, transmitted);
        aa_load_events.insert(aa_load_events.end(), {
            static_cast<double>(i),
            static_cast<double>(j),
            static_cast<double>(current_step),
            static_cast<double>(old_status),
            static_cast<double>(new_status),
            motor_only,
            opposed_drive,
            transmitted,
            combined_candidate,
            selected,
            actin.direction[i].dot(actin.direction[j])});
    }
}
