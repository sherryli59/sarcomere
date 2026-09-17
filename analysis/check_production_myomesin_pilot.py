#!/usr/bin/env python3
"""Summarize structural order and sparse-module topology in a production pilot."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np


def percentile(values, q):
    return float(np.percentile(values, q)) if len(values) else None


def band_order(x, lx, harmonics=(1, 2, 3, 4)):
    x = np.asarray(x) % lx
    return max(
        float(abs(np.exp(2j * np.pi * harmonic * x / lx).mean()))
        for harmonic in harmonics
    )


def component_sizes(n_vertices, edges):
    parent = list(range(n_vertices))
    size = [1] * n_vertices

    def find(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    for first, second in edges:
        root_first, root_second = find(first), find(second)
        if root_first == root_second:
            continue
        if size[root_first] < size[root_second]:
            root_first, root_second = root_second, root_first
        parent[root_second] = root_first
        size[root_first] += size[root_second]
    active = {vertex for edge in edges for vertex in edge}
    counts = {}
    for vertex in active:
        root = find(vertex)
        counts[root] = counts.get(root, 0) + 1
    return sorted(counts.values(), reverse=True)


def motor_pair_event_counts(events: np.ndarray) -> dict[str, int]:
    """Count identity events, excluding code-5 rows with no held identity."""
    counts = {
        name: int(np.sum(events[:, 3].astype(int) == code)) if len(events) else 0
        for code, name in (
            (1, "selected"),
            (2, "switched"),
            (3, "dropped_after_support_timeout"),
            (4, "locked_at_maturation"),
        )
    }
    if len(events):
        detachment_rows = events[events[:, 3].astype(int) == 5]
        valid_release = (
            (detachment_rows[:, 4] >= 0) & (detachment_rows[:, 5] >= 0)
        )
        counts["released_by_actn2_detachment"] = int(np.sum(valid_release))
        counts["actn2_detachment_without_held_identity"] = int(
            len(detachment_rows) - np.sum(valid_release)
        )
    else:
        counts["released_by_actn2_detachment"] = 0
        counts["actn2_detachment_without_held_identity"] = 0
    return counts


def summarize(path: Path, tail_fraction: float = 0.2) -> dict:
    with h5py.File(path, "r", locking=False) as handle:
        meta = handle["meta"].attrs
        dt = float(meta["dt"])
        lx = float(meta.get("box_x", meta.get("Lx", 4.825)))
        centers_a = handle["actin/center"][:]
        directions_a = handle["actin/direction"][:]
        centers_m = handle["myosin/center"][:]
        steps = handle["state/current_step"][:, 0].astype(int)
        n_myosins = centers_m.shape[1]
        status2 = (
            handle["actin/aa_local_status2_count"][:]
            if "actin/aa_local_status2_count" in handle else None
        )
        motor_pair_events = (
            handle["catch_bond/motor_pair_events"][:]
            if "catch_bond/motor_pair_events" in handle else np.empty((0, 20))
        )
        if (
            "catch_bond/local_link_motor_support" in handle
            and "catch_bond/local_link_mechanics_v2_frame_offsets" in handle
            and len(handle["catch_bond/local_link_mechanics_v2_frame_offsets"])
        ):
            all_motor_support = handle[
                "catch_bond/local_link_motor_support"
            ][:]
            support_offsets = handle[
                "catch_bond/local_link_mechanics_v2_frame_offsets"
            ][-1]
            support_start = int(support_offsets[2])
            support_count = int(support_offsets[3])
            final_motor_support = handle[
                "catch_bond/local_link_motor_support"
            ][support_start:support_start + support_count]
        else:
            all_motor_support = np.empty((0, 20))
            final_motor_support = np.empty((0, 20))
        force_caps = {
            species: float(np.sum(
                handle[f"{species}/cumulative_force_cap_events"][-1]
            ))
            if f"{species}/cumulative_force_cap_events" in handle else None
            for species in ("actin", "myosin")
        }
        module_enabled = "myomesin/frame_offsets" in handle
        if module_enabled:
            offsets = handle["myomesin/frame_offsets"][:]
            state = handle["myomesin/bridge_state"][:]
            events = handle["myomesin/events"][:] \
                if "myomesin/events" in handle else np.empty((0, 9))
        else:
            offsets = np.empty((0, 3), dtype=int)
            state = np.empty((0, 8))
            events = np.empty((0, 9))
        parameters = {
            "k_mm": float(meta.get("k_mm", 0.0)),
            "myomesin_radial_k": float(meta.get("myomesin_radial_k", 0.0)),
            "myomesin_force_max": float(meta.get("myomesin_force_max", 0.0)),
            "myomesin_sites_per_filament": int(
                meta.get("myomesin_sites_per_filament", 0)
            ),
            "myomesin_axial_capture": float(
                meta.get("myomesin_axial_capture", 0.0)
            ),
            "aa_motor_pair_commit_time": float(
                meta.get("aa_motor_pair_commit_time", 0.0)
            ),
            "aa_status2_force_on": float(meta.get("aa_status2_force_on", 0.0)),
        }

    final = len(steps) - 1
    actin_axes = directions_a[final] / np.linalg.norm(
        directions_a[final], axis=1, keepdims=True
    )
    plus = centers_a[final, :, 0] + 0.5 * actin_axes[:, 0]
    minus = centers_a[final, :, 0] - 0.5 * actin_axes[:, 0]
    # Cell A is a persistent Z-band phenotype, not a favorable single frame.
    # Score the last 20% of the trajectory as well as the endpoint so a noisy
    # transient cannot satisfy the replication gate.
    tail_start = max(0, int(np.floor(0.8 * len(steps))))
    z_band_tail = []
    for frame in range(tail_start, len(steps)):
        axes = directions_a[frame] / np.linalg.norm(
            directions_a[frame], axis=1, keepdims=True
        )
        plus_frame = centers_a[frame, :, 0] + 0.5 * axes[:, 0]
        minus_frame = centers_a[frame, :, 0] - 0.5 * axes[:, 0]
        z_band_tail.append(max(
            band_order(plus_frame, lx), band_order(minus_frame, lx)
        ))
    structural = {
        "time_ms": float(steps[final] * dt * 1e3),
        "m_band_order": band_order(centers_m[final, :, 0], lx),
        "z_band_order": max(band_order(plus, lx), band_order(minus, lx)),
        "z_band_order_tail_median": float(np.median(z_band_tail)),
        "z_band_order_tail_min": float(np.min(z_band_tail)),
        "a_band_order": band_order(centers_a[final, :, 0], lx),
        "status2_links_final": (
            float(np.sum(status2[final, :, 0]) / 2.0)
            if status2 is not None else None
        ),
        "force_cap_events": force_caps,
    }
    if status2 is not None:
        status2_links = np.sum(status2[:, :, 0], axis=1) / 2.0
        nonzero_status2 = np.flatnonzero(status2_links > 0)
        structural["status2_links_max"] = float(np.max(status2_links))
        structural["status2_first_time_ms"] = (
            float(steps[nonzero_status2[0]] * dt * 1e3)
            if len(nonzero_status2) else None
        )

    event_counts = motor_pair_event_counts(motor_pair_events)
    motor_pair_annealing = {"event_counts": event_counts}
    if len(all_motor_support):
        saved_selected = (
            (all_motor_support[:, 0] >= 0) & (all_motor_support[:, 1] >= 0)
        )
        saved_challenger = (
            (all_motor_support[:, 14] >= 0) & (all_motor_support[:, 15] >= 0)
        )
        motor_pair_annealing.update({
            "saved_support_rows": int(len(all_motor_support)),
            "saved_selected_identity_rows": int(np.sum(saved_selected)),
            "saved_challenger_identity_rows": int(np.sum(saved_challenger)),
            "saved_challenger_support_pn_max": float(
                np.max(all_motor_support[:, 18], initial=0.0)
            ),
            "saved_challenger_age_ms_max": float(
                1e3 * np.max(all_motor_support[:, 19], initial=0.0)
            ),
        })
        selected_rows = all_motor_support[saved_selected]
        commit_time = parameters["aa_motor_pair_commit_time"]
        force_on = parameters["aa_status2_force_on"]
        age_ready = selected_rows[:, 6] + 1e-15 >= commit_time
        force_ready = selected_rows[:, 4] > force_on
        motor_pair_annealing.update({
            "saved_candidate_age_ms_max": float(
                1e3 * np.max(selected_rows[:, 6], initial=0.0)
            ),
            "saved_candidate_support_ema_pn_max": float(
                np.max(selected_rows[:, 4], initial=0.0)
            ),
            "saved_candidate_rows_age_ready": int(np.sum(age_ready)),
            "saved_candidate_rows_force_ready": int(np.sum(force_ready)),
            "saved_candidate_rows_age_and_force_ready": int(
                np.sum(age_ready & force_ready)
            ),
        })
    else:
        motor_pair_annealing.update({
            "saved_support_rows": 0,
            "saved_selected_identity_rows": 0,
            "saved_challenger_identity_rows": 0,
            "saved_challenger_support_pn_max": None,
            "saved_challenger_age_ms_max": None,
            "saved_candidate_age_ms_max": None,
            "saved_candidate_support_ema_pn_max": None,
            "saved_candidate_rows_age_ready": 0,
            "saved_candidate_rows_force_ready": 0,
            "saved_candidate_rows_age_and_force_ready": 0,
        })
    if len(final_motor_support):
        valid = final_motor_support[:, 0] >= 0
        valid_challenger = (
            (final_motor_support[:, 14] >= 0)
            & (final_motor_support[:, 15] >= 0)
        )
        motor_pair_annealing.update({
            "final_actn2_links_evaluated": int(len(final_motor_support)),
            "final_selected_identities": int(np.sum(valid)),
            "final_locked_identities": int(np.sum(final_motor_support[:, 7] > 0.5)),
            "final_challenger_identities": int(np.sum(valid_challenger)),
            "final_candidate_age_ms_max": float(
                1e3 * np.max(final_motor_support[:, 6], initial=0.0)
            ),
            "final_support_ema_pn_max": float(
                np.max(final_motor_support[:, 4], initial=0.0)
            ),
            "final_support_instantaneous_pn_max": float(
                np.max(final_motor_support[:, 5], initial=0.0)
            ),
            "final_challenger_support_pn_max": float(
                np.max(final_motor_support[:, 18], initial=0.0)
            ),
            "final_challenger_age_ms_max": float(
                1e3 * np.max(final_motor_support[:, 19], initial=0.0)
            ),
        })
    else:
        motor_pair_annealing.update({
            "final_actn2_links_evaluated": 0,
            "final_selected_identities": 0,
            "final_locked_identities": 0,
            "final_challenger_identities": 0,
            "final_candidate_age_ms_max": None,
            "final_support_ema_pn_max": None,
            "final_support_instantaneous_pn_max": None,
            "final_challenger_support_pn_max": None,
            "final_challenger_age_ms_max": None,
        })
    result = {
        "trajectory": str(path),
        "parameters": parameters,
        "structural": structural,
        "motor_pair_annealing": motor_pair_annealing,
        "module_enabled": module_enabled,
    }
    if not module_enabled:
        return result

    frame_rows = []
    observations = {}
    final_ids = set()
    duplicate_edges = 0
    max_degree_seen = 0
    for frame, (_, start, count) in enumerate(offsets):
        rows = state[int(start):int(start + count)]
        edges = []
        ids = []
        degree = np.zeros(n_myosins, dtype=int)
        for row in rows:
            module_id = (int(row[0]) << 32) | int(row[1])
            pair = tuple(sorted((int(row[2]), int(row[3]))))
            ids.append(module_id)
            edges.append(pair)
            degree[pair[0]] += 1
            degree[pair[1]] += 1
            observations.setdefault(module_id, []).append(
                (int(offsets[frame, 0]), 1e3 * abs(float(row[5])))
            )
        duplicate_edges += len(edges) - len(set(edges))
        max_degree_seen = max(max_degree_seen, int(degree.max(initial=0)))
        components = component_sizes(n_myosins, edges)
        frame_rows.append({
            "step": int(offsets[frame, 0]),
            "count": len(edges),
            "largest_component": components[0] if components else 0,
            "max_degree": int(degree.max(initial=0)),
        })
        if frame == len(offsets) - 1:
            final_ids = set(ids)

    formation = {}
    association_count = release_count = invalid_release_count = 0
    for row in events:
        code = int(row[1])
        module_id = int(row[2])
        if code == 1:
            association_count += 1
            formation[module_id] = 1e3 * abs(float(row[5]))
        elif code == 2:
            release_count += 1
        elif code == 3:
            invalid_release_count += 1

    reductions = []
    persistent_reductions = []
    for module_id, samples in observations.items():
        if module_id not in formation or formation[module_id] < 20.0:
            continue
        initial = formation[module_id]
        final_mismatch = float(np.median([value for _, value in samples[-5:]]))
        reduction = (initial - final_mismatch) / initial
        reductions.append(reduction)
        if module_id in final_ids:
            persistent_reductions.append(reduction)

    tail_start = max(0, int((1.0 - tail_fraction) * len(frame_rows)))
    tail = frame_rows[tail_start:]
    counts = np.asarray([row["count"] for row in tail], dtype=float)
    largest = np.asarray(
        [row["largest_component"] for row in tail], dtype=float
    )
    axial = np.abs(state[:, 5]) * 1e3 if len(state) else np.asarray([])
    radial = state[:, 6] * 1e3 if len(state) else np.asarray([])
    forces = state[:, 7] if len(state) else np.asarray([])
    topology = {
        "association_events": association_count,
        "stochastic_release_events": release_count,
        "invalid_geometry_release_events": invalid_release_count,
        "duplicate_pair_rows": duplicate_edges,
        "max_degree": max_degree_seen,
        "module_count_tail_median": percentile(counts, 50),
        "module_count_tail_max": percentile(counts, 100),
        "largest_component_tail_median": percentile(largest, 50),
        "largest_component_tail_max": percentile(largest, 100),
        "largest_component_fraction_tail_max": (
            float(largest.max() / n_myosins) if len(largest) else 0.0
        ),
        "axial_mismatch_nm_median": percentile(axial, 50),
        "axial_mismatch_nm_p90": percentile(axial, 90),
        "radial_spacing_nm_median": percentile(radial, 50),
        "radial_spacing_nm_p10": percentile(radial, 10),
        "radial_spacing_nm_p90": percentile(radial, 90),
        "module_force_pn_max": percentile(forces, 100),
        "annealing_modules_evaluated": len(reductions),
        "annealing_reduction_fraction_median": percentile(reductions, 50),
        "annealing_fraction_positive": (
            float(np.mean(np.asarray(reductions) > 0.0)) if reductions else None
        ),
        "annealing_fraction_at_least_10pct": (
            float(np.mean(np.asarray(reductions) >= 0.10)) if reductions else None
        ),
        "persistent_modules_evaluated": len(persistent_reductions),
        "persistent_reduction_fraction_median": percentile(
            persistent_reductions, 50
        ),
    }
    gates = {
        "finite_valency": max_degree_seen <= parameters["myomesin_sites_per_filament"],
        "unique_pair_edges": duplicate_edges == 0,
        "no_system_spanning_component": (
            topology["largest_component_tail_max"]
            <= max(2, int(np.floor(0.45 * n_myosins)))
        ),
        "population_annealing": (
            len(reductions) > 0
            and topology["annealing_reduction_fraction_median"] > 0.0
            and topology["annealing_fraction_positive"] >= 0.5
        ),
    }
    result["topology"] = topology
    result["gates"] = gates
    result["pass"] = all(gates.values())
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trajectory", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--tail-fraction", type=float, default=0.2)
    parser.add_argument("--require-pass", action="store_true")
    args = parser.parse_args()
    result = summarize(args.trajectory, args.tail_fraction)
    rendered = json.dumps(result, indent=2)
    print(rendered)
    if args.output:
        args.output.write_text(rendered + "\n")
    return 1 if args.require_pass and not result.get("pass", True) else 0


if __name__ == "__main__":
    raise SystemExit(main())
