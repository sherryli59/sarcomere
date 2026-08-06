#!/usr/bin/env python3
"""Measure scaffold-free transverse organization around one thick filament."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


PATTERN = re.compile(r"(?P<condition>[a-zA-Z0-9_.-]+)_rep-(?P<rep>\d+)")


def read_vector(handle: h5py.File, group: str, prefix: str) -> np.ndarray:
    return np.stack([
        np.asarray(handle[f"{group}/{prefix}_{axis}"])[:, 0, 0]
        for axis in "xyz"
    ], axis=1)


def one_run(path: Path) -> tuple[dict[str, float | str], np.ndarray]:
    match = PATTERN.fullmatch(path.stem)
    if match is None:
        raise ValueError(path)
    with h5py.File(path, "r") as handle:
        actin_center = np.asarray(handle["actin/center"], dtype=float)
        actin_direction = np.asarray(handle["actin/direction"], dtype=float)
        myosin_center = np.asarray(handle["myosin/center"], dtype=float)[:, 0, :]
        myosin_direction = np.asarray(handle["myosin/direction"], dtype=float)[:, 0, :]
        e1 = read_vector(handle, "myosin", "material_e1")
        e2 = read_vector(handle, "myosin", "material_e2")
        bound = np.asarray(handle["state/am_head_bound_actin_current"], dtype=int)[:, 0, :]
        am_pair_count = np.asarray(
            handle["state/am_pairs_current_count"], dtype=int)[:, 0]
        aa_count = np.asarray(
            handle["state/aa_pairs_current_count"], dtype=int)[:, 0]
        aa_status_padded = np.asarray(
            handle["state/aa_status_current"], dtype=int)[:, :, 0]
        aa_pairs_padded = np.asarray(
            handle["state/aa_pairs_current"], dtype=int)
        steps = np.asarray(handle["state/current_step"], dtype=int)[:, 0]
        aa_status1 = np.array([
            np.count_nonzero(aa_status_padded[t, :aa_count[t]] == 1)
            for t in range(len(aa_count))
        ])
        aa_status2 = np.array([
            np.count_nonzero(aa_status_padded[t, :aa_count[t]] == 2)
            for t in range(len(aa_count))
        ])
        metadata = handle["meta"].attrs
        rho = float(metadata.get("aa_spacing_ratio", np.nan))
        d_am = float(metadata.get("am_preferred_centerline", np.nan))
        d_aa = float(metadata.get("aa_preferred_centerline", np.nan))
        dt = float(metadata.get("dt", np.nan))
        refractory_time = float(metadata.get("aa_refractory_time", 0.0))
        am_model_version = int(metadata.get("am_model_version", 0))
        actin_force_caps = (
            np.asarray(handle["actin/cumulative_force_cap_events"], dtype=float)
            if "actin/cumulative_force_cap_events" in handle else None)
        actin_torque_caps = (
            np.asarray(handle["actin/cumulative_torque_cap_events"], dtype=float)
            if "actin/cumulative_torque_cap_events" in handle else None)
    displacement = actin_center - myosin_center[:, None, :]
    axial = np.sum(displacement * myosin_direction[:, None, :], axis=2)
    transverse = displacement - axial[:, :, None] * myosin_direction[:, None, :]
    y = np.sum(transverse * e1[:, None, :], axis=2)
    z = np.sum(transverse * e2[:, None, :], axis=2)
    radius = np.sqrt(y * y + z * z)
    theta = np.arctan2(z, y)
    psi6 = np.abs(np.mean(np.exp(6j * theta), axis=1))
    gaps = np.diff(np.concatenate([
        np.sort(np.mod(theta, 2 * np.pi), axis=1),
        np.sort(np.mod(theta, 2 * np.pi), axis=1)[:, :1] + 2 * np.pi,
    ], axis=1), axis=1)
    gap_cv = np.std(gaps, axis=1) / np.mean(gaps, axis=1)
    radius_cv = np.std(radius, axis=1) / np.maximum(np.mean(radius, axis=1), 1e-12)
    alignment = np.abs(np.sum(
        actin_direction * myosin_direction[:, None, :], axis=2))
    nematic = np.mean(0.5 * (3.0 * alignment * alignment - 1.0), axis=1)
    distinct = (
        np.array([len(set(frame[frame >= 0].tolist())) for frame in bound])
        if am_model_version > 0 else am_pair_count.astype(float))
    intended_edges = {
        tuple(sorted((index, (index + 1) % 6))) for index in range(6)
    }
    observed_edges: list[set[tuple[int, int]]] = []
    for frame, count in zip(aa_pairs_padded, aa_count):
        observed_edges.append({
            tuple(sorted((int(pair[0]), int(pair[1]))))
            for pair in frame[:count] if pair[0] >= 0 and pair[1] >= 0
        })
    correct_count = np.array([
        len(edges & intended_edges) for edges in observed_edges], dtype=float)
    precision = np.array([
        len(edges & intended_edges) / len(edges) if edges else 1.0
        for edges in observed_edges], dtype=float)
    recall = correct_count / len(intended_edges)

    # A macro-association ignores short molecular gaps no longer than the
    # explicitly recorded refractory interval. This is distinct from changing
    # the KMC state or silently treating a detached pair as bonded.
    edge_occupied = np.zeros((len(observed_edges), len(intended_edges)), dtype=bool)
    for edge_index, edge in enumerate(sorted(intended_edges)):
        edge_occupied[:, edge_index] = [edge in edges for edges in observed_edges]
    saved_step = float(np.median(np.diff(steps))) if len(steps) > 1 else 1.0
    gap_frames = int(np.floor(
        refractory_time / (dt * saved_step))) if dt > 0 and saved_step > 0 else 0
    macro_occupied = edge_occupied.copy()
    if gap_frames > 0:
        for edge_index in range(macro_occupied.shape[1]):
            occupied_indices = np.flatnonzero(edge_occupied[:, edge_index])
            for left, right in zip(occupied_indices[:-1], occupied_indices[1:]):
                if right - left - 1 <= gap_frames:
                    macro_occupied[left:right + 1, edge_index] = True

    target_basin = (psi6 >= 0.9) & (gap_cv <= 0.15)
    sustained_window = min(5, len(target_basin))
    first_sustained = -1
    for frame in range(0, len(target_basin) - sustained_window + 1):
        if np.all(target_basin[frame:frame + sustained_window]):
            first_sustained = frame
            break
    tail = slice(max(1, int(0.8 * len(psi6))), None)
    row: dict[str, float | str] = {
        "path": str(path), **match.groupdict(),
        "rho": rho,
        "d_am_preferred": d_am,
        "d_aa_preferred": d_aa,
        "psi6_initial": float(psi6[0]),
        "psi6_tail": float(np.mean(psi6[tail])),
        "psi6_change": float(np.mean(psi6[tail]) - psi6[0]),
        "gap_cv_tail": float(np.mean(gap_cv[tail])),
        "radius_mean_tail": float(np.mean(radius[tail])),
        "radius_cv_tail": float(np.mean(radius_cv[tail])),
        "nematic_tail": float(np.mean(nematic[tail])),
        "distinct_bound_actins_tail": float(np.mean(distinct[tail])),
        "aa_bonds_tail": float(np.mean(aa_count[tail])),
        "aa_status1_tail": float(np.mean(aa_status1[tail])),
        "aa_status2_tail": float(np.mean(aa_status2[tail])),
        "correct_neighbor_precision_tail": float(np.mean(precision[tail])),
        "correct_neighbor_recall_tail": float(np.mean(recall[tail])),
        "correct_edge_molecular_occupancy_tail": float(np.mean(edge_occupied[tail])),
        "correct_edge_macro_occupancy_tail": float(np.mean(macro_occupied[tail])),
        "target_basin_fraction": float(np.mean(target_basin)),
        "first_sustained_target_time": (
            float(steps[first_sustained] * dt)
            if first_sustained >= 0 and dt > 0 else np.nan),
        "actin_force_cap_fraction": (
            float(np.sum(actin_force_caps[-1]) /
                  max(1.0, steps[-1] * actin_center.shape[1]))
            if actin_force_caps is not None else np.nan),
        "actin_torque_cap_fraction": (
            float(np.sum(actin_torque_caps[-1]) /
                  max(1.0, steps[-1] * actin_center.shape[1]))
            if actin_torque_caps is not None else np.nan),
    }
    return row, np.stack([y[-1], z[-1]], axis=1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    args = parser.parse_args()
    results = [one_run(path) for path in args.inputs]
    rows = sorted((result[0] for result in results),
                  key=lambda row: (row["condition"], int(row["rep"])))
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    with args.output_prefix.with_suffix(".csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    metrics = [key for key in rows[0] if key not in {
        "path", "condition", "rep", "rho", "d_am_preferred",
        "d_aa_preferred"}]
    aggregate = []
    for condition in sorted({str(row["condition"]) for row in rows}):
        subset = [row for row in rows if row["condition"] == condition]
        out: dict[str, float | str] = {
            "condition": condition,
            "rho": float(subset[0]["rho"]),
            "d_am_preferred": float(subset[0]["d_am_preferred"]),
            "d_aa_preferred": float(subset[0]["d_aa_preferred"]),
        }
        for metric in metrics:
            values = np.array([row[metric] for row in subset], dtype=float)
            out[f"{metric}_mean"] = float(values.mean())
            out[f"{metric}_std"] = float(values.std(ddof=1))
        aggregate.append(out)
    aggregate_path = args.output_prefix.with_name(
        args.output_prefix.name + "_aggregate").with_suffix(".csv")
    with aggregate_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(aggregate[0]))
        writer.writeheader(); writer.writerows(aggregate)

    rep0 = {row["condition"]: coords for row, coords in results if row["rep"] == "0"}
    conditions = sorted(rep0)
    ncols = min(5, len(conditions))
    nrows = int(np.ceil(len(conditions) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.5 * ncols, 2.5 * nrows),
                             squeeze=False)
    for axis, condition in zip(axes.flat, conditions):
        coords = rep0[condition]
        axis.scatter(coords[:, 0], coords[:, 1], c=np.arange(6), cmap="hsv", s=30)
        axis.scatter([0], [0], marker="x", color="black")
        axis.set_title(condition.replace("_", " "))
        axis.set_aspect("equal")
        axis.set_xlabel("material-frame y")
        axis.set_ylabel("material-frame z")
    for axis in axes.flat[len(conditions):]:
        axis.set_visible(False)
    fig.tight_layout()
    fig.savefig(args.output_prefix.with_suffix(".png"), dpi=200)

    labels = [str(row["condition"]) for row in aggregate]
    x = np.arange(len(labels))
    fig, axes = plt.subplots(2, 2, figsize=(8, 6), sharex=True)
    panels = [
        ("psi6_tail", r"$|\psi_6|$ (higher is better)"),
        ("gap_cv_tail", "angular-gap CV (lower is better)"),
        ("distinct_bound_actins_tail", "actins bound per myosin"),
    ]
    for axis, (metric, ylabel) in zip(axes.flat[:3], panels):
        means = [float(row[f"{metric}_mean"]) for row in aggregate]
        errors = [float(row[f"{metric}_std"]) for row in aggregate]
        axis.errorbar(x, means, yerr=errors, fmt="o", capsize=3)
        axis.set_ylabel(ylabel)
    status1 = np.array([float(row["aa_status1_tail_mean"]) for row in aggregate])
    status2 = np.array([float(row["aa_status2_tail_mean"]) for row in aggregate])
    axes.flat[3].bar(x, status1, label="status 1")
    axes.flat[3].bar(x, status2, bottom=status1, label="status 2")
    axes.flat[3].set_ylabel("AA bonds")
    axes.flat[3].legend(frameon=False)
    for axis in axes.flat:
        axis.set_xticks(x, labels, rotation=30, ha="right")
    fig.tight_layout()
    fig.savefig(args.output_prefix.with_name(
        args.output_prefix.name + "_metrics").with_suffix(".png"), dpi=200)


if __name__ == "__main__":
    main()
