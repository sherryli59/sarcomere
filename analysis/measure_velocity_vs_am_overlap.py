#!/usr/bin/env python3
"""Measure crosslinked-pair extension velocity versus productive A-M overlap."""

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


def productive_overlap(
    actin_center: np.ndarray,
    actin_direction: np.ndarray,
    myosin_center: np.ndarray,
    myosin_direction: np.ndarray,
    box: np.ndarray,
    actin_length: float,
    myosin_length: float,
    bare_fraction: float,
) -> float:
    axis = myosin_direction / np.linalg.norm(myosin_direction)
    thin_axis = actin_direction / np.linalg.norm(actin_direction)
    displacement = minimum_image(actin_center - myosin_center, box)
    projected_center = float(displacement @ axis)
    projected_half_length = 0.5 * actin_length * abs(float(thin_axis @ axis))
    actin_min = projected_center - projected_half_length
    actin_max = projected_center + projected_half_length
    productive_half = -1 if float(thin_axis @ axis) > 0.0 else 1
    bare_half = 0.5 * bare_fraction * myosin_length
    head_min = -0.5 * myosin_length if productive_half < 0 else bare_half
    head_max = -bare_half if productive_half < 0 else 0.5 * myosin_length
    return max(0.0, min(actin_max, head_max) - max(actin_min, head_min))


def pair_cluster_interval(
    values: np.ndarray,
    pairs: np.ndarray,
    statistic,
    draws: int = 3000,
    seed: int = 20260805,
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    unique = np.unique(pairs)
    by_pair = {pair: np.flatnonzero(pairs == pair) for pair in unique}
    estimates = []
    for _ in range(draws):
        sampled = rng.choice(unique, size=len(unique), replace=True)
        indices = np.concatenate([by_pair[pair] for pair in sampled])
        estimates.append(statistic(indices))
    return tuple(map(float, np.percentile(estimates, [2.5, 97.5])))


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
        spacing = float(np.median(np.diff(frame_time)))
        width = int(round(args.window_ms * 1e-3 / spacing))
        start = int(np.searchsorted(frame_time, args.start_ms * 1e-3 - 1e-12))
        if width < 1 or not np.isclose(width * spacing, args.window_ms * 1e-3):
            raise ValueError("window must be an integer multiple of frame spacing")

        actin_center = handle["actin/center"][:]
        actin_direction = handle["actin/direction"][:]
        myosin_center = handle["myosin/center"][:]
        myosin_direction = handle["myosin/direction"][:]
        head_state = handle["state/am_head_state_current"][:]
        bound_actin = handle["state/am_head_bound_actin_current"][:]
        head_count = (
            handle["actin/am_weak_count"][:, :, 0]
            + handle["actin/am_strong_count"][:, :, 0]
        )
        box = np.array([meta["box_x"], meta["box_y"], meta["box_z"]], dtype=float)
        actin_length = float(meta["actin_length"])
        myosin_length = float(meta["myosin_length"])
        bare_fraction = float(meta["am_bare_zone_fraction"])
        n_frames, n_actins = actin_center.shape[:2]

        # The overlap observable is restricted to myosins carrying at least one
        # material weak/strong head on that actin in the saved frame. This
        # avoids assigning an unrelated nearby thick filament as its A-M owner.
        attached_am_overlap = np.zeros((n_frames, n_actins), dtype=float)
        for frame in range(n_frames):
            attached = np.argwhere(head_state[frame] > 0)
            partners: dict[int, set[int]] = defaultdict(set)
            for myosin_id, head_id in attached:
                actin_id = int(bound_actin[frame, myosin_id, head_id])
                if actin_id >= 0:
                    partners[actin_id].add(int(myosin_id))
            for actin_id, myosin_ids in partners.items():
                attached_am_overlap[frame, actin_id] = max(
                    productive_overlap(
                        actin_center[frame, actin_id],
                        actin_direction[frame, actin_id],
                        myosin_center[frame, myosin_id],
                        myosin_direction[frame, myosin_id],
                        box,
                        actin_length,
                        myosin_length,
                        bare_fraction,
                    )
                    for myosin_id in myosin_ids
                )

        offsets = handle["catch_bond/local_link_mechanics_v2_frame_offsets"][:]
        identities = handle["catch_bond/local_link_identity"][:]
        pair_frames: list[dict[tuple[int, int], dict[str, float]]] = []
        for row in offsets:
            _, identity_start, _, count = map(int, row)
            frame_identity = identities[identity_start : identity_start + count]
            accumulator: dict[tuple[int, int], dict[str, float]] = defaultdict(
                lambda: {"links": 0.0, "status2": 0.0}
            )
            for identity in frame_identity:
                pair = tuple(sorted((int(identity[2]), int(identity[3]))))
                accumulator[pair]["links"] += 1.0
                accumulator[pair]["status2"] += float(identity[4] == 2)
            pair_frames.append(dict(accumulator))

        samples: list[dict[str, float | int | str]] = []
        for first in range(start, n_frames - width, width):
            last = first + width
            common_pairs = set(pair_frames[first])
            for frame in range(first + 1, last + 1):
                common_pairs.intersection_update(pair_frames[frame])

            for actin_i, actin_j in common_pairs:
                direction_dot = float(np.mean(np.sum(
                    actin_direction[first : last + 1, actin_i]
                    * actin_direction[first : last + 1, actin_j], axis=1
                )))
                if direction_dot > -0.5:
                    continue
                heads_i = float(np.mean(head_count[first : last + 1, actin_i]))
                heads_j = float(np.mean(head_count[first : last + 1, actin_j]))
                if min(heads_i, heads_j) <= 0.5:
                    continue

                def aa_overlap(frame: int) -> float:
                    barbed_i = (
                        actin_center[frame, actin_i]
                        - 0.5 * actin_length * actin_direction[frame, actin_i]
                    )
                    barbed_j = (
                        actin_center[frame, actin_j]
                        - 0.5 * actin_length * actin_direction[frame, actin_j]
                    )
                    connector = minimum_image(barbed_j - barbed_i, box)
                    return 0.5 * (
                        connector @ actin_direction[frame, actin_i]
                        - connector @ actin_direction[frame, actin_j]
                    )

                velocity = (aa_overlap(first) - aa_overlap(last)) / (
                    frame_time[last] - frame_time[first]
                )
                overlap_i = float(np.mean(attached_am_overlap[first : last + 1, actin_i]))
                overlap_j = float(np.mean(attached_am_overlap[first : last + 1, actin_j]))
                mean_overlap = 0.5 * (overlap_i + overlap_j)
                status2 = float(np.mean([
                    pair_frames[frame][(actin_i, actin_j)]["status2"]
                    for frame in range(first, last + 1)
                ]))
                state = "cb2" if status2 > 0.5 else "status1"
                samples.append({
                    "window_start_ms": float(frame_time[first] * 1e3),
                    "window_end_ms": float(frame_time[last] * 1e3),
                    "actin_i": actin_i,
                    "actin_j": actin_j,
                    "pair_id": f"{actin_i}:{actin_j}",
                    "state": state,
                    "mean_status2_links": status2,
                    "mean_attached_am_overlap_i_um": overlap_i,
                    "mean_attached_am_overlap_j_um": overlap_j,
                    "pair_mean_attached_am_overlap_um": mean_overlap,
                    "aa_extension_velocity_um_per_s": float(velocity),
                })

    with (args.output_dir / "velocity_vs_am_overlap_samples.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(samples[0]))
        writer.writeheader()
        writer.writerows(samples)

    overlap = np.asarray([row["pair_mean_attached_am_overlap_um"] for row in samples])
    velocity = np.asarray([row["aa_extension_velocity_um_per_s"] for row in samples])
    state = np.asarray([row["state"] for row in samples])
    pair_id = np.asarray([row["pair_id"] for row in samples])
    bins = np.array([0.0, 0.15, 0.30, 0.45, 0.60, 0.676])
    summary: dict[str, object] = {
        "trajectory": str(args.trajectory),
        "analysis_window_ms": [args.start_ms, float(frame_time[-1] * 1e3)],
        "averaging_window_ms": args.window_ms,
        "overlap_semantics": "mean of the two actins' maximum productive overlap with currently head-attached myosins",
        "velocity_semantics": "positive means A-A overlap shrinking; zero means stall; negative means re-annealing",
        "selection": "antiparallel pair continuously ACTN2-linked and both actins mean >0.5 A-M heads over the window",
        "groups": {},
        "bins": [],
    }
    binned_rows = []
    group_colors = {"status1": "#0f8a72", "cb2": "#7c3aed"}
    for group in ("status1", "cb2"):
        mask = state == group
        x, y, p = overlap[mask], velocity[mask], pair_id[mask]
        slope, intercept = np.polyfit(x, y, 1)
        slope_ci = pair_cluster_interval(
            x, p, lambda index: np.polyfit(x[index], y[index], 1)[0], seed=31 if group == "status1" else 47
        )
        mean_ci = pair_cluster_interval(y, p, lambda index: np.mean(y[index]), seed=61 if group == "status1" else 79)
        summary["groups"][group] = {
            "sample_count": int(len(x)),
            "unique_pair_count": int(len(np.unique(p))),
            "mean_am_overlap_um": float(np.mean(x)),
            "median_am_overlap_um": float(np.median(x)),
            "mean_velocity_um_per_s": float(np.mean(y)),
            "mean_velocity_pair_clustered_95ci": list(mean_ci),
            "slope_um_per_s_per_um_overlap": float(slope),
            "slope_pair_clustered_95ci": list(slope_ci),
            "intercept_um_per_s": float(intercept),
            "pearson_r": float(np.corrcoef(x, y)[0, 1]),
        }
        for index, (left, right) in enumerate(zip(bins[:-1], bins[1:])):
            selected = mask & (overlap >= left) & (overlap < right)
            if not np.any(selected):
                continue
            values = velocity[selected]
            rng = np.random.default_rng(1000 + index + (0 if group == "status1" else 100))
            bootstrap = np.mean(rng.choice(values, size=(3000, len(values)), replace=True), axis=1)
            row = {
                "group": group,
                "overlap_bin_left_um": left,
                "overlap_bin_right_um": right,
                "sample_count": int(len(values)),
                "mean_am_overlap_um": float(np.mean(overlap[selected])),
                "mean_velocity_um_per_s": float(np.mean(values)),
                "velocity_ci_low": float(np.percentile(bootstrap, 2.5)),
                "velocity_ci_high": float(np.percentile(bootstrap, 97.5)),
            }
            binned_rows.append(row)
    summary["bins"] = binned_rows
    with (args.output_dir / "summary.json").open("w") as stream:
        json.dump(summary, stream, indent=2)
    with (args.output_dir / "velocity_vs_am_overlap_bins.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(binned_rows[0]))
        writer.writeheader()
        writer.writerows(binned_rows)

    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.5), sharex=True, sharey=True, constrained_layout=True)
    for ax, group, title in zip(axes, ("status1", "cb2"), ("Weak ACTN2 pairs (status 1)", "Mature ACTN2 pairs (cb2)")):
        mask = state == group
        color = group_colors[group]
        ax.scatter(overlap[mask], velocity[mask], s=18, alpha=0.22, color=color, edgecolor="none")
        rows = [row for row in binned_rows if row["group"] == group]
        bx = np.asarray([row["mean_am_overlap_um"] for row in rows])
        by = np.asarray([row["mean_velocity_um_per_s"] for row in rows])
        lo = np.asarray([row["velocity_ci_low"] for row in rows])
        hi = np.asarray([row["velocity_ci_high"] for row in rows])
        ax.errorbar(bx, by, yerr=np.vstack((by - lo, hi - by)), marker="o", markersize=6,
                    color=color, linewidth=2, capsize=3, label="bin mean ± 95% CI")
        slope = summary["groups"][group]["slope_um_per_s_per_um_overlap"]
        intercept = summary["groups"][group]["intercept_um_per_s"]
        fit_x = np.linspace(0.0, 0.675, 100)
        ax.plot(fit_x, intercept + slope * fit_x, color=color, linestyle="--", linewidth=1.5,
                label=f"linear slope {slope:+.1f} s⁻¹")
        ax.axhline(0.0, color="#4b5563", linewidth=1)
        ax.set_title(title)
        ax.set_xlabel("Mean attached productive A–M overlap (µm)")
        ax.set_xlim(0.0, 0.69)
        ax.grid(alpha=0.18)
        ax.legend(frameon=False, fontsize=8)
    axes[0].set_ylabel("A–A extension velocity (µm/s)\npositive = overlap shrinking; zero = stall")
    fig.savefig(args.output_dir / "velocity_vs_am_overlap.png", dpi=220)
    fig.savefig(args.output_dir / "velocity_vs_am_overlap.pdf")
    plt.close(fig)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
