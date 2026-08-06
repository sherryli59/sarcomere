#!/usr/bin/env python3
"""Measure endogenous ACTN2 force–velocity samples from a production trajectory."""

from __future__ import annotations

import argparse
import csv
import json
import os
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/sarcomere-mplconfig")

import h5py
import matplotlib.pyplot as plt
import numpy as np


def minimum_image(displacement: np.ndarray, box: np.ndarray) -> np.ndarray:
    return displacement - box * np.round(displacement / box)


def linear_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    slope, intercept = np.polyfit(x, y, 1)
    correlation = float(np.corrcoef(x, y)[0, 1])
    return float(slope), float(intercept), correlation


def clustered_slope_interval(
    x: np.ndarray,
    y: np.ndarray,
    pair_ids: np.ndarray,
    draws: int = 3000,
    seed: int = 20260805,
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    unique = np.unique(pair_ids)
    by_pair = {pair: np.flatnonzero(pair_ids == pair) for pair in unique}
    slopes = []
    for _ in range(draws):
        selected = rng.choice(unique, size=len(unique), replace=True)
        indices = np.concatenate([by_pair[pair] for pair in selected])
        if np.ptp(x[indices]) <= 1e-12:
            continue
        slopes.append(np.polyfit(x[indices], y[indices], 1)[0])
    return tuple(map(float, np.percentile(slopes, [2.5, 97.5])))


def clustered_mean_interval(
    values: np.ndarray,
    pair_ids: np.ndarray,
    draws: int = 3000,
    seed: int = 20260806,
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    unique = np.unique(pair_ids)
    by_pair = {pair: np.flatnonzero(pair_ids == pair) for pair in unique}
    means = []
    for _ in range(draws):
        selected = rng.choice(unique, size=len(unique), replace=True)
        indices = np.concatenate([by_pair[pair] for pair in selected])
        means.append(np.mean(values[indices]))
    return tuple(map(float, np.percentile(means, [2.5, 97.5])))


def bootstrap_mean_interval(values: np.ndarray, draws: int = 3000, seed: int = 1) -> tuple[float, float]:
    if len(values) == 1:
        return float(values[0]), float(values[0])
    rng = np.random.default_rng(seed)
    means = np.mean(rng.choice(values, size=(draws, len(values)), replace=True), axis=1)
    return tuple(map(float, np.percentile(means, [2.5, 97.5])))


def summarize_group(x: np.ndarray, y: np.ndarray, pairs: np.ndarray) -> dict[str, float | int | list[float]]:
    slope, intercept, correlation = linear_fit(x, y)
    low, high = clustered_slope_interval(x, y, pairs)
    mean_low, mean_high = clustered_mean_interval(y, pairs)
    return {
        "sample_count": int(len(x)),
        "unique_pair_count": int(len(np.unique(pairs))),
        "mean_load_pN_per_actin": float(np.mean(x)),
        "median_load_pN_per_actin": float(np.median(x)),
        "mean_extension_velocity_um_per_s": float(np.mean(y)),
        "mean_velocity_pair_clustered_95ci": [mean_low, mean_high],
        "slope_um_per_s_per_pN": slope,
        "slope_pair_clustered_95ci": [low, high],
        "intercept_um_per_s": intercept,
        "linear_zero_crossing_pN": float(-intercept / slope) if slope != 0 else float("nan"),
        "pearson_r": correlation,
        "conventional_decreasing_fv_supported": bool(high < 0.0),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("trajectory", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--start-ms", type=float, default=20.0)
    parser.add_argument("--window-ms", type=float, default=5.0)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    with h5py.File(args.trajectory, "r") as handle:
        meta = handle["meta"].attrs
        dt = float(meta["dt"])
        frame_steps = handle["state/current_step"][:, 0].astype(int)
        frame_time = frame_steps * dt
        frame_spacing = float(np.median(np.diff(frame_time)))
        width = int(round((args.window_ms * 1e-3) / frame_spacing))
        if width < 1 or not np.isclose(width * frame_spacing, args.window_ms * 1e-3):
            raise ValueError("window must be an integer multiple of the saved-frame spacing")
        start = int(np.searchsorted(frame_time, args.start_ms * 1e-3 - 1e-12))

        offsets = handle["catch_bond/local_link_mechanics_v2_frame_offsets"][:]
        identities = handle["catch_bond/local_link_identity"][:]
        mechanics = handle["catch_bond/local_link_mechanics_v2"][:]
        mechanics_columns = meta["aa_local_link_mechanics_v2_columns"]
        if isinstance(mechanics_columns, bytes):
            mechanics_columns = mechanics_columns.decode()
        mech_col = {name: index for index, name in enumerate(mechanics_columns.split(","))}

        # Each frame maps an antiparallel actin pair to total favorable
        # bilateral ACTN2 load and the number of mature (status-2) molecules.
        pair_frames: list[dict[tuple[int, int], dict[str, float]]] = []
        for row in offsets:
            _, identity_start, mechanics_start, count = map(int, row)
            frame_identity = identities[identity_start : identity_start + count]
            frame_mechanics = mechanics[mechanics_start : mechanics_start + count]
            accumulator: dict[tuple[int, int], dict[str, float]] = defaultdict(
                lambda: {"load": 0.0, "links": 0.0, "status2": 0.0}
            )
            for identity, mechanical in zip(frame_identity, frame_mechanics):
                actin_i, actin_j = sorted((int(identity[2]), int(identity[3])))
                pair = (actin_i, actin_j)
                accumulator[pair]["load"] += max(0.0, float(mechanical[mech_col["bilateral_load"]]))
                accumulator[pair]["links"] += 1.0
                accumulator[pair]["status2"] += float(identity[4] == 2)
            pair_frames.append(dict(accumulator))

        centers = handle["actin/center"][:]
        directions = handle["actin/direction"][:]
        actin_length = float(meta["actin_length"])
        head_count = (
            handle["actin/am_weak_count"][:, :, 0]
            + handle["actin/am_strong_count"][:, :, 0]
        )
        box = np.array([meta["box_x"], meta["box_y"], meta["box_z"]], dtype=float)

        samples: list[dict[str, float | int | str]] = []
        for first in range(start, len(frame_time) - width, width):
            last = first + width
            common_pairs = set(pair_frames[first])
            for frame in range(first + 1, last + 1):
                common_pairs.intersection_update(pair_frames[frame])

            for actin_i, actin_j in common_pairs:
                direction_dot = float(np.mean(np.sum(
                    directions[first : last + 1, actin_i]
                    * directions[first : last + 1, actin_j], axis=1
                )))
                if direction_dot > -0.5:
                    continue
                load = float(np.mean([
                    pair_frames[frame][(actin_i, actin_j)]["load"]
                    for frame in range(first, last + 1)
                ]))
                status2 = float(np.mean([
                    pair_frames[frame][(actin_i, actin_j)]["status2"]
                    for frame in range(first, last + 1)
                ]))
                mean_heads = float(np.mean(
                    head_count[first : last + 1, [actin_i, actin_j]]
                ))
                if mean_heads <= 0.5:
                    continue
                # The minimal balanced diagnostic defines extension from the
                # material barbed ends, not lab-frame center velocity. This
                # cancels pair-centroid drift and remains valid under rotation.
                def aa_overlap(frame: int) -> float:
                    barbed_i = centers[frame, actin_i] - 0.5 * actin_length * directions[frame, actin_i]
                    barbed_j = centers[frame, actin_j] - 0.5 * actin_length * directions[frame, actin_j]
                    connector = minimum_image(barbed_j - barbed_i, box)
                    return 0.5 * (
                        connector @ directions[frame, actin_i]
                        - connector @ directions[frame, actin_j]
                    )

                extension_velocity = (aa_overlap(first) - aa_overlap(last)) / (
                    frame_time[last] - frame_time[first]
                )
                state = "cb2" if status2 > 0.5 else "status1"
                samples.append({
                    "window_start_ms": float(frame_time[first] * 1e3),
                    "window_end_ms": float(frame_time[last] * 1e3),
                    "actin_i": actin_i,
                    "actin_j": actin_j,
                    "pair_id": f"{actin_i}:{actin_j}",
                    "state": state,
                    "mean_status2_links": status2,
                    "mean_am_heads_per_actin": mean_heads,
                    "actn2_load_pN_per_actin": load,
                    "aa_extension_velocity_um_per_s": float(extension_velocity),
                })

    with (args.output_dir / "force_velocity_samples.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(samples[0]))
        writer.writeheader()
        writer.writerows(samples)

    load = np.asarray([row["actn2_load_pN_per_actin"] for row in samples], dtype=float)
    velocity = np.asarray([row["aa_extension_velocity_um_per_s"] for row in samples], dtype=float)
    state = np.asarray([row["state"] for row in samples])
    pair_id = np.asarray([row["pair_id"] for row in samples])

    group_masks = {
        "all_crosslinked": np.ones(len(samples), dtype=bool),
        "status1": state == "status1",
        "cb2": state == "cb2",
    }
    summary = {
        "trajectory": str(args.trajectory),
        "analysis_start_ms": args.start_ms,
        "analysis_end_ms": float(frame_time[-1] * 1e3),
        "window_ms": args.window_ms,
        "positive_velocity_semantics": "A-A overlap shrinking (motor-driven extension)",
        "negative_velocity_semantics": "A-A overlap increasing (re-annealing)",
        "load_semantics": "sum of favorable bilateral ACTN2 axial load for the pair, pN per actin",
        "selection": "antiparallel pair continuously ACTN2-linked for full window and mean >0.5 A-M heads per actin",
        "groups": {
            name: summarize_group(load[mask], velocity[mask], pair_id[mask])
            for name, mask in group_masks.items()
        },
    }

    fixed_edges = np.array([0.0, 0.5, 1.0, 2.0, 4.0, np.inf])
    binned_rows = []
    for group, mask in group_masks.items():
        for bin_index, (left, right) in enumerate(zip(fixed_edges[:-1], fixed_edges[1:])):
            in_bin = mask & (load >= left) & (load < right)
            values = velocity[in_bin]
            if len(values) == 0:
                continue
            low, high = bootstrap_mean_interval(values, seed=100 + bin_index)
            binned_rows.append({
                "group": group,
                "load_bin_left_pN": left,
                "load_bin_right_pN": right,
                "sample_count": int(len(values)),
                "mean_load_pN_per_actin": float(np.mean(load[in_bin])),
                "mean_velocity_um_per_s": float(np.mean(values)),
                "velocity_ci_low": low,
                "velocity_ci_high": high,
            })
    summary["bins"] = binned_rows
    with (args.output_dir / "summary.json").open("w") as stream:
        json.dump(summary, stream, indent=2)
    with (args.output_dir / "force_velocity_bins.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(binned_rows[0]))
        writer.writeheader()
        writer.writerows(binned_rows)

    # Aligned facets expose the pooled relationship and the distinct weak and
    # mature states; combining them alone would obscure the maturation step.
    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.5), sharey=True, constrained_layout=True)
    panels = [
        ("all_crosslinked", "All continuously crosslinked pairs"),
        ("status1", "Weak ACTN2 pairs (status 1)"),
        ("cb2", "Pairs carrying mature ACTN2 (cb2)"),
    ]
    colors = {"all_crosslinked": "#2563a6", "status1": "#0f8a72", "cb2": "#7c3aed"}
    for ax, (group, title) in zip(axes, panels):
        mask = group_masks[group]
        ax.scatter(load[mask], velocity[mask], s=14, alpha=0.18, color=colors[group], edgecolor="none")
        rows = [row for row in binned_rows if row["group"] == group]
        bx = np.asarray([row["mean_load_pN_per_actin"] for row in rows])
        by = np.asarray([row["mean_velocity_um_per_s"] for row in rows])
        lo = np.asarray([row["velocity_ci_low"] for row in rows])
        hi = np.asarray([row["velocity_ci_high"] for row in rows])
        ax.errorbar(bx, by, yerr=np.vstack((by - lo, hi - by)), color=colors[group],
                    marker="o", markersize=6, linewidth=2, capsize=3, label="bin mean ± 95% CI")
        slope = summary["groups"][group]["slope_um_per_s_per_pN"]
        intercept = summary["groups"][group]["intercept_um_per_s"]
        fit_x = np.linspace(0, min(8.0, np.percentile(load[mask], 97.5)), 100)
        ax.plot(fit_x, intercept + slope * fit_x, linestyle="--", linewidth=1.5,
                color=colors[group], label=f"linear slope {slope:+.2f} µm s⁻¹ pN⁻¹")
        ax.axhline(0, color="#4b5563", linewidth=1)
        ax.set_title(title)
        ax.set_xlabel("Endogenous ACTN2 axial load (pN per actin)")
        ax.grid(alpha=0.18)
        ax.legend(frameon=False, fontsize=8)
    axes[0].set_ylabel("A–A extension velocity (µm/s)\npositive = overlap shrinking; zero = stall")
    axes[0].set_xlim(left=0)
    axes[1].set_xlim(left=0)
    axes[2].set_xlim(left=0)
    fig.savefig(args.output_dir / "crosslinked_force_velocity.png", dpi=220)
    fig.savefig(args.output_dir / "crosslinked_force_velocity.pdf")
    plt.close(fig)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
