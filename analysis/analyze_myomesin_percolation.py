#!/usr/bin/env python3
"""Measure sparse-myomesin graph percolation and component phase coherence."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import h5py
import numpy as np


def graph_components(n_nodes: int, pairs: set[tuple[int, int]]) -> dict:
    parent = list(range(n_nodes))
    sizes = [1] * n_nodes
    degree = np.zeros(n_nodes, dtype=int)

    def find(node: int) -> int:
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def union(first: int, second: int) -> None:
        first, second = find(first), find(second)
        if first == second:
            return
        if sizes[first] < sizes[second]:
            first, second = second, first
        parent[second] = first
        sizes[first] += sizes[second]

    for first, second in pairs:
        if first < 0 or second >= n_nodes or first >= second:
            raise ValueError(f"invalid canonical pair {(first, second)}")
        degree[first] += 1
        degree[second] += 1
        union(first, second)

    members: dict[int, list[int]] = {}
    for node in range(n_nodes):
        members.setdefault(find(node), []).append(node)
    component_sizes = sorted(
        (len(nodes) for nodes in members.values()), reverse=True
    )
    cycle_rank = len(pairs) - n_nodes + len(component_sizes)
    return {
        "members": list(members.values()),
        "degree": degree,
        "component_sizes": component_sizes,
        "component_count": len(component_sizes),
        "largest_component": component_sizes[0] if component_sizes else 0,
        "second_largest_component": (
            component_sizes[1] if len(component_sizes) > 1 else 0
        ),
        "cycle_rank": cycle_rank,
    }


def phase_metrics(
    x: np.ndarray, lx: float, components: list[list[int]]
) -> dict:
    harmonic_orders = np.asarray([
        abs(np.exp(2j * np.pi * harmonic * (x % lx) / lx).mean())
        for harmonic in (1, 2, 3, 4)
    ])
    harmonic = int(np.argmax(harmonic_orders) + 1)
    phase = np.exp(2j * np.pi * harmonic * (x % lx) / lx)
    local_order = []
    component_phase = []
    weights = []
    for nodes in components:
        mean = phase[nodes].mean()
        magnitude = abs(mean)
        local_order.append(magnitude)
        component_phase.append(mean / magnitude if magnitude else 0j)
        weights.append(len(nodes))
    weights_array = np.asarray(weights, dtype=float)
    return {
        "m_band_order": float(harmonic_orders[harmonic - 1]),
        "best_harmonic": harmonic,
        "within_component_order": float(np.average(
            np.asarray(local_order), weights=weights_array
        )),
        "between_component_phase_order": float(abs(np.average(
            np.asarray(component_phase), weights=weights_array
        ))),
    }


def analyze(path: Path, frame_stride: int = 1) -> dict:
    if frame_stride < 1:
        raise ValueError("frame_stride must be positive")
    with h5py.File(path, "r", locking=False) as handle:
        required = (
            "state/current_step",
            "myosin/center",
            "myomesin/frame_offsets",
            "myomesin/bridge_state",
        )
        missing = [name for name in required if name not in handle]
        if missing:
            raise KeyError("missing datasets: " + ", ".join(missing))
        n_frames = min(handle[name].shape[0] for name in required[:3])
        steps = handle["state/current_step"][:n_frames, 0].astype(int)
        offsets = handle["myomesin/frame_offsets"][:n_frames]
        centers = handle["myosin/center"]
        bridge_state = handle["myomesin/bridge_state"]
        n_nodes = int(centers.shape[1])
        meta = handle["meta"].attrs
        lx = float(meta.get("box_x", meta.get("Lx", 4.825)))
        dt = float(meta["dt"])
        sites = int(meta.get("myomesin_sites_per_filament", 0))
        records = []
        for frame in range(0, n_frames, frame_stride):
            step, start, count = map(int, offsets[frame])
            if step != int(steps[frame]):
                raise ValueError(f"offset step mismatch at frame {frame}")
            if start < 0 or count < 0 or start + count > bridge_state.shape[0]:
                # A live writer may have extended offsets just before the
                # bridge-state rows become visible. Exclude that partial tail.
                break
            rows = bridge_state[start:start + count]
            pairs = {
                tuple(sorted((int(row[2]), int(row[3])))) for row in rows
            }
            if len(pairs) != len(rows):
                raise ValueError(f"duplicate bridge pair at frame {frame}")
            graph = graph_components(n_nodes, pairs)
            phase = phase_metrics(
                centers[frame, :, 0], lx, graph["members"]
            )
            degree = graph["degree"]
            largest = graph["largest_component"]
            records.append({
                "frame": frame,
                "step": int(steps[frame]),
                "time_s": float(steps[frame] * dt),
                "bridge_count": len(pairs),
                "mean_degree": float(2 * len(pairs) / n_nodes),
                "covered_fraction": float(np.mean(degree > 0)),
                "degree_two_fraction": float(np.mean(degree == 2)),
                "free_site_fraction": float(
                    np.sum(np.maximum(0, sites - degree)) / (sites * n_nodes)
                ) if sites else None,
                "component_count": graph["component_count"],
                "largest_component": largest,
                "largest_component_fraction": float(largest / n_nodes),
                "second_largest_component_fraction": float(
                    graph["second_largest_component"] / n_nodes
                ),
                "cycle_rank": graph["cycle_rank"],
                "degree_one_nodes": int(np.sum(degree == 1)),
                "percolated_45pct": bool(largest >= 0.45 * n_nodes),
                **phase,
            })
    return {
        "trajectory": str(path),
        "n_myosins": n_nodes,
        "box_x_um": lx,
        "dt_s": dt,
        "myomesin_sites_per_filament": sites,
        "frames_available": n_frames,
        "frame_stride": frame_stride,
        "records": records,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trajectory", type=Path)
    parser.add_argument("--frame-stride", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--csv", type=Path)
    args = parser.parse_args()
    report = analyze(args.trajectory, args.frame_stride)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        records = report["records"]
        with args.csv.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
    final = report["records"][-1]
    print(json.dumps({
        "trajectory": report["trajectory"],
        "frames_analyzed": len(report["records"]),
        "final": final,
        "maximum_largest_component_fraction": max(
            row["largest_component_fraction"] for row in report["records"]
        ),
        "ever_percolated_45pct": any(
            row["percolated_45pct"] for row in report["records"]
        ),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
