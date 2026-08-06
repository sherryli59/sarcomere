#!/usr/bin/env python3
"""Diagnose ACTN2 contour rupture and plot A-M overlap for final cb2 actins."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/sarcomere-mplconfig")

import h5py
import matplotlib.pyplot as plt
import numpy as np


def percentile_summary(values: np.ndarray) -> dict[str, float]:
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return {key: float("nan") for key in ("min", "p10", "p25", "median", "p75", "p90", "p99", "max")}
    quantiles = np.percentile(values, [0, 10, 25, 50, 75, 90, 99, 100])
    return dict(zip(("min", "p10", "p25", "median", "p75", "p90", "p99", "max"), map(float, quantiles)))


def minimum_image(vector: np.ndarray, box: np.ndarray, periodic: np.ndarray) -> np.ndarray:
    result = np.asarray(vector, dtype=float).copy()
    for axis in range(3):
        if periodic[axis] and box[axis] > 0:
            result[axis] -= box[axis] * np.round(result[axis] / box[axis])
    return result


def productive_overlap(
    actin_center: np.ndarray,
    actin_direction: np.ndarray,
    myosin_center: np.ndarray,
    myosin_direction: np.ndarray,
    box: np.ndarray,
    periodic: np.ndarray,
    actin_length: float,
    myosin_length: float,
    bare_zone_fraction: float,
) -> float:
    axis = np.asarray(myosin_direction, dtype=float)
    axis /= np.linalg.norm(axis)
    actin_axis = np.asarray(actin_direction, dtype=float)
    actin_axis /= np.linalg.norm(actin_axis)
    displacement = minimum_image(actin_center - myosin_center, box, periodic)
    center_projection = float(displacement @ axis)
    half_projection = 0.5 * actin_length * abs(float(actin_axis @ axis))
    actin_min = center_projection - half_projection
    actin_max = center_projection + half_projection
    productive_half = -1 if float(actin_axis @ axis) > 0.0 else 1
    bare_half = 0.5 * bare_zone_fraction * myosin_length
    head_min = -0.5 * myosin_length if productive_half < 0 else bare_half
    head_max = -bare_half if productive_half < 0 else 0.5 * myosin_length
    return max(0.0, min(actin_max, head_max) - max(actin_min, head_min))


def preceding_frame_indices(event_steps: np.ndarray, frame_steps: np.ndarray) -> np.ndarray:
    return np.maximum(0, np.searchsorted(frame_steps, event_steps, side="right") - 1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("trajectory", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    with h5py.File(args.trajectory, "r") as handle:
        meta = handle["meta"].attrs
        columns = meta["aa_local_link_events_v3_columns"]
        if isinstance(columns, bytes):
            columns = columns.decode()
        column_names = columns.split(",")
        col = {name: index for index, name in enumerate(column_names)}
        events = handle["catch_bond/local_link_events_v3"][:]
        detachments = events[events[:, col["event_code"]] == 2]
        geometric = detachments[detachments[:, col["geometric_detachment"]] == 1]
        stochastic = detachments[detachments[:, col["stochastic_detachment"]] == 1]
        associations = events[events[:, col["event_code"]] == 1]
        association_by_id = {int(row[col["link_id"]]): row for row in associations}
        formation = np.array([association_by_id[int(row[col["link_id"]])] for row in geometric])

        frame_steps = handle["state/current_step"][:, 0].astype(int)
        frame_index = preceding_frame_indices(geometric[:, col["step"]].astype(int), frame_steps)
        dt = float(meta["dt"])

        # Use saved per-link connector vectors at the latest output frame before
        # rupture. This is an approximate direction decomposition (<=1 ms old),
        # while the event row supplies the exact rupture length.
        offsets = handle["catch_bond/local_link_mechanics_v2_frame_offsets"][:]
        identities = handle["catch_bond/local_link_identity"][:]
        mechanics = handle["catch_bond/local_link_mechanics_v2"][:]
        actin_direction = handle["actin/direction"][:]
        axial_components = []
        transverse_components = []
        decomposition_age_ms = []
        decomposition_found = []
        for event, fi in zip(geometric, frame_index):
            step, identity_start, mechanics_start, count = map(int, offsets[fi])
            identity_slice = identities[identity_start : identity_start + count]
            matches = np.flatnonzero(identity_slice[:, 1] == int(event[col["link_id"]]))
            if matches.size == 0:
                axial_components.append(np.nan)
                transverse_components.append(np.nan)
                decomposition_age_ms.append(np.nan)
                decomposition_found.append(False)
                continue
            local = int(matches[0])
            connector = mechanics[mechanics_start + local, :3]
            connector_length = float(np.linalg.norm(connector))
            direction = actin_direction[fi, int(event[col["actin_i"]])]
            direction /= np.linalg.norm(direction)
            axial = abs(float(connector @ direction))
            transverse = np.sqrt(max(0.0, connector_length**2 - axial**2))
            # Preserve the latest saved connector's direction, but scale its
            # components to the exact event-time contour length.
            scale = float(event[col["connector_length"]]) / max(connector_length, 1e-12)
            axial_components.append(axial * scale)
            transverse_components.append(transverse * scale)
            decomposition_age_ms.append((float(event[col["step"]]) - step) * dt * 1e3)
            decomposition_found.append(True)

        axial_components = np.asarray(axial_components)
        transverse_components = np.asarray(transverse_components)
        decomposition_found = np.asarray(decomposition_found, dtype=bool)

        # Approximate external state at the preceding saved frame. These are
        # correlations, not an exact additive per-link force decomposition.
        pre_event = {}
        for dataset in (
            "am_motor_force",
            "am_tether_force",
            "titin_force",
            "am_strong_count",
            "am_weak_count",
            "titin_bound",
        ):
            data = handle[f"actin/{dataset}"][:, :, 0]
            values_i = data[frame_index, geometric[:, col["actin_i"]].astype(int)]
            values_j = data[frame_index, geometric[:, col["actin_j"]].astype(int)]
            pre_event[dataset] = np.maximum(values_i, values_j)

        # Final cb2 actins: maximum productive axial overlap over every current
        # A-M neighbor. A cb2 actin with no current A-M pair is retained at zero.
        final_cb = handle["actin/cb_status"][-1, :, 0]
        cb2_actins = np.flatnonzero(final_cb == 2)
        pair_count = int(handle["state/am_pairs_current_count"][-1, 0])
        am_pairs = handle["state/am_pairs_current"][-1, :pair_count].astype(int)
        final_actin_center = handle["actin/center"][-1]
        final_actin_direction = handle["actin/direction"][-1]
        final_myosin_center = handle["myosin/center"][-1]
        final_myosin_direction = handle["myosin/direction"][-1]
        box = np.array([meta["box_x"], meta["box_y"], meta["box_z"]], dtype=float)
        periodic = np.array([meta["periodic_x"], meta["periodic_y"], meta["periodic_z"]], dtype=bool)
        actin_length = float(meta["actin_length"])
        myosin_length = float(meta["myosin_length"])
        bare_fraction = float(meta["am_bare_zone_fraction"])
        head_bearing_length = 0.5 * myosin_length * (1.0 - bare_fraction)
        cb2_overlap = []
        cb2_partner_count = []
        for actin_id in cb2_actins:
            partners = am_pairs[am_pairs[:, 0] == actin_id, 1]
            overlaps = [
                productive_overlap(
                    final_actin_center[actin_id], final_actin_direction[actin_id],
                    final_myosin_center[myosin_id], final_myosin_direction[myosin_id],
                    box, periodic, actin_length, myosin_length, bare_fraction,
                )
                for myosin_id in partners
            ]
            cb2_overlap.append(max(overlaps, default=0.0))
            cb2_partner_count.append(len(partners))
        cb2_overlap = np.asarray(cb2_overlap)

        exact_length = geometric[:, col["connector_length"]]
        formation_length = formation[:, col["connector_length"]]
        registry_at_break = geometric[:, col["registry_error"]]
        registry_at_formation = formation[:, col["registry_error"]]
        status_at_break = geometric[:, col["old_status"]].astype(int)
        tension = geometric[:, col["single_link_tension"]]
        event_age_ms = (geometric[:, col["step"]] - formation[:, col["step"]]) * dt * 1e3

        summary = {
            "trajectory": str(args.trajectory),
            "final_time_ms": float(frame_steps[-1] * dt * 1e3),
            "association_count": int(len(associations)),
            "geometric_detachment_count": int(len(geometric)),
            "stochastic_detachment_count": int(len(stochastic)),
            "geometric_status1_count": int(np.sum(status_at_break == 1)),
            "geometric_status2_count": int(np.sum(status_at_break == 2)),
            "geometric_status1_fraction": float(np.mean(status_at_break == 1)),
            "formation_length_nm": percentile_summary(formation_length * 1e3),
            "rupture_length_nm": percentile_summary(exact_length * 1e3),
            "length_growth_nm": percentile_summary((exact_length - formation_length) * 1e3),
            "link_age_at_geometric_rupture_ms": percentile_summary(event_age_ms),
            "registry_at_formation_nm": percentile_summary(registry_at_formation * 1e3),
            "registry_at_rupture_nm": percentile_summary(registry_at_break * 1e3),
            "registry_change_nm": percentile_summary((registry_at_break - registry_at_formation) * 1e3),
            "single_link_tension_at_rupture_pN": percentile_summary(tension),
            "fraction_born_above_90nm": float(np.mean(formation_length > 0.09)),
            "fraction_born_above_95nm": float(np.mean(formation_length > 0.095)),
            "fraction_registry_over_40nm_at_rupture": float(np.mean(registry_at_break > 0.04)),
            "fraction_productive_support_positive_at_rupture": float(np.mean(geometric[:, col["productive_support"]] > 0)),
            "decomposition_coverage": float(np.mean(decomposition_found)),
            "decomposition_snapshot_age_ms": percentile_summary(np.asarray(decomposition_age_ms)[decomposition_found]),
            "estimated_axial_component_at_rupture_nm": percentile_summary(axial_components[decomposition_found] * 1e3),
            "estimated_transverse_component_at_rupture_nm": percentile_summary(transverse_components[decomposition_found] * 1e3),
            "estimated_transverse_fraction_at_rupture": percentile_summary((transverse_components / exact_length)[decomposition_found]),
            "pre_event_max_arm_motor_force_pN": percentile_summary(pre_event["am_motor_force"]),
            "pre_event_max_arm_head_tether_force_pN": percentile_summary(pre_event["am_tether_force"]),
            "pre_event_max_arm_titin_force_pN": percentile_summary(pre_event["titin_force"]),
            "pre_event_either_titin_bound_fraction": float(np.mean(pre_event["titin_bound"] > 0)),
            "pre_event_max_arm_strong_heads": percentile_summary(pre_event["am_strong_count"]),
            "pre_event_max_arm_weak_heads": percentile_summary(pre_event["am_weak_count"]),
            "cb2_actin_count": int(len(cb2_actins)),
            "cb2_zero_am_pair_count": int(np.sum(np.asarray(cb2_partner_count) == 0)),
            "cb2_max_productive_am_overlap_um": percentile_summary(cb2_overlap),
            "productive_half_length_um": head_bearing_length,
        }

    with (args.output_dir / "summary.json").open("w") as stream:
        json.dump(summary, stream, indent=2)

    with (args.output_dir / "cb2_am_axial_overlap.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["actin_id", "cb_status", "current_am_pair_count", "max_productive_am_axial_overlap_um"])
        for actin_id, partners, overlap in zip(cb2_actins, cb2_partner_count, cb2_overlap):
            writer.writerow([int(actin_id), 2, int(partners), float(overlap)])

    # Plot requested by the user. A histogram plus an individual-value strip
    # avoids hiding a small-n distribution behind a smooth density estimate.
    fig, ax = plt.subplots(figsize=(8.2, 4.8), constrained_layout=True)
    bins = np.linspace(0.0, head_bearing_length, 10)
    counts, edges, patches = ax.hist(cb2_overlap, bins=bins, color="#2b6cb0", edgecolor="white", linewidth=1.2)
    ax.scatter(cb2_overlap, np.full_like(cb2_overlap, -0.22), marker="|", s=180, linewidth=1.8, color="#153e75", clip_on=False)
    median = float(np.median(cb2_overlap)) if len(cb2_overlap) else 0.0
    ax.axvline(median, color="#d97706", linewidth=2.2, label=f"median = {median:.3f} µm")
    ax.axvline(head_bearing_length, color="#4a5568", linestyle="--", linewidth=1.6, label=f"full head-bearing half = {head_bearing_length:.3f} µm")
    ax.set_xlim(0, head_bearing_length * 1.04)
    ax.set_ylim(bottom=-0.5)
    ax.set_xlabel("Maximum productive A–M axial overlap per cb2 actin (µm)")
    ax.set_ylabel("Number of cb2 actins")
    ax.set_title(f"Final 40 ms snapshot: cb2 actins (n = {len(cb2_overlap)})")
    ax.grid(axis="y", alpha=0.22)
    ax.legend(frameon=False, loc="upper right")
    fig.savefig(args.output_dir / "cb2_am_axial_overlap_distribution.png", dpi=220)
    fig.savefig(args.output_dir / "cb2_am_axial_overlap_distribution.pdf")
    plt.close(fig)

    # Compact mechanism figure: entry geometry, contour threshold, and the
    # latest-frame direction decomposition for every geometric rupture.
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.0), constrained_layout=True)
    axes[0].hist(formation_length * 1e3, bins=np.linspace(0, 125, 26), color="#63b3ed", label="at formation")
    axes[0].hist(exact_length * 1e3, bins=np.linspace(0, 125, 26), color="#e53e3e", alpha=0.78, label="at rupture")
    axes[0].axvline(100, color="#2b6cb0", linestyle="--", label="100 nm capture")
    axes[0].axvline(120, color="#9b2c2c", linestyle="--", label="120 nm rupture")
    axes[0].set(xlabel="ACTN2 connector length (nm)", ylabel="Events", title="Narrow geometric margin")
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].scatter(registry_at_formation * 1e3, (registry_at_break - registry_at_formation) * 1e3, s=8, alpha=0.28, color="#805ad5")
    axes[1].axhline(0, color="#4a5568", linewidth=1)
    axes[1].set(xlabel="Registry error at formation (nm)", ylabel="Registry change before rupture (nm)", title="Axial registry barely changes")
    axes[2].scatter(axial_components[decomposition_found] * 1e3, transverse_components[decomposition_found] * 1e3, s=9, alpha=0.32, color="#2f855a")
    axes[2].set(xlabel="Estimated axial connector component (nm)", ylabel="Estimated transverse component (nm)", title="Rupture is mostly transverse/3D")
    for axis in axes:
        axis.grid(alpha=0.18)
    fig.savefig(args.output_dir / "geometric_rupture_mechanism.png", dpi=220)
    fig.savefig(args.output_dir / "geometric_rupture_mechanism.pdf")
    plt.close(fig)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
