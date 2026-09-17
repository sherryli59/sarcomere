#!/usr/bin/env python3
"""Validate persistent Z-band replication against historical Cell A."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import h5py
import numpy as np


def band_order(x: np.ndarray, lx: float) -> float:
    wrapped = np.asarray(x) % lx
    return max(
        float(abs(np.exp(2j * np.pi * harmonic * wrapped / lx).mean()))
        for harmonic in (1, 2, 3, 4)
    )


def z_band_series(path: Path) -> dict:
    with h5py.File(path, "r", locking=False) as handle:
        frame_counts = {
            name: int(handle[name].shape[0]) for name in (
                "state/current_step",
                "actin/center",
                "actin/direction",
                "myosin/center",
            )
        }
        common_frames = min(frame_counts.values())
        if common_frames == 0:
            raise ValueError("trajectory has no common committed frame")
        # A writer can expose a newly extended dataset before all peer datasets
        # have been flushed.  Restrict calculations to the common prefix and
        # retain the raw counts so the completion gate fails explicitly.
        steps = handle["state/current_step"][:common_frames, 0].astype(np.int64)
        centers = handle["actin/center"][:common_frames]
        directions = handle["actin/direction"][:common_frames]
        meta = handle["meta"].attrs
        lx = float(meta.get("box_x", meta.get("Lx", 4.825)))
        dt = float(meta["dt"])
        actin_length = float(meta.get("actin_length", 1.0))
        metadata = {
            "dt": dt,
            "k_mm": float(meta.get("k_mm", 0.0)),
            "k_aa": float(meta.get("k_aa", np.nan)),
            "myosin_steric_k": float(meta.get("myosin_steric_k", np.nan)),
            "myosin_exclusion_distance": float(
                meta.get("myosin_exclusion_distance", np.nan)
            ),
            "aa_status_mode": int(meta.get("aa_status_mode", -1)),
            "n_actins": int(centers.shape[1]),
            "n_myosins": int(handle["myosin/center"].shape[1]),
        }
    norms = np.linalg.norm(directions, axis=2, keepdims=True)
    axes = directions / np.maximum(norms, np.finfo(float).tiny)
    plus = centers[:, :, 0] + 0.5 * actin_length * axes[:, :, 0]
    minus = centers[:, :, 0] - 0.5 * actin_length * axes[:, :, 0]
    z_order = np.asarray([
        max(band_order(plus[i], lx), band_order(minus[i], lx))
        for i in range(len(steps))
    ])
    return {
        "steps": steps,
        "time_ms": steps.astype(float) * dt * 1e3,
        "z_order": z_order,
        "metadata": metadata,
        "frame_counts": frame_counts,
    }


def summarize(series: dict, tail_fraction: float = 0.2) -> dict:
    z_order = np.asarray(series["z_order"])
    if len(z_order) == 0:
        raise ValueError("trajectory has no saved frames")
    tail_start = max(0, int(np.floor((1.0 - tail_fraction) * len(z_order))))
    tail = z_order[tail_start:]
    return {
        "frames": int(len(z_order)),
        "final_step": int(series["steps"][-1]),
        "final_time_ms": float(series["time_ms"][-1]),
        "z_band_order_initial": float(z_order[0]),
        "z_band_order_final": float(z_order[-1]),
        "z_band_order_max": float(np.max(z_order)),
        "z_band_order_tail_median": float(np.median(tail)),
        "z_band_order_tail_min": float(np.min(tail)),
        "z_band_order_gain": float(z_order[-1] - z_order[0]),
        "tail_fraction": tail_fraction,
        "dataset_frame_counts": series.get("frame_counts"),
    }


def evaluate(
    reference: dict,
    candidate: dict,
    expected_step: int = 500_000,
    endpoint_min: float = 0.60,
    tail_median_min: float = 0.60,
) -> dict:
    reference_summary = summarize(reference)
    candidate_summary = summarize(candidate)
    meta = candidate["metadata"]
    expected_metadata = {
        "dt": 1e-5,
        "k_mm": 0.0,
        "k_aa": 1000.0,
        "myosin_steric_k": 1000.0,
        "myosin_exclusion_distance": 0.0,
        "aa_status_mode": 7,
        "n_actins": 800,
        "n_myosins": 320,
    }
    metadata_gates = {}
    for name, expected in expected_metadata.items():
        observed = meta.get(name)
        if isinstance(expected, float):
            passed = observed is not None and np.isclose(
                observed, expected, rtol=0.0, atol=1e-12
            )
        else:
            passed = observed == expected
        metadata_gates[name] = {
            "observed": observed, "expected": expected, "passed": bool(passed)
        }
    frame_counts = candidate.get("frame_counts")
    frame_counts_complete = (
        frame_counts is None
        or all(count == expected_step // 500 + 1
               for count in frame_counts.values())
    )
    gates = {
        "complete_5s_trajectory": (
            candidate_summary["final_step"] == expected_step
            and candidate_summary["frames"] == expected_step // 500 + 1
            and frame_counts_complete
        ),
        "cell_a_metadata": all(item["passed"] for item in metadata_gates.values()),
        "strong_endpoint_z_order": (
            candidate_summary["z_band_order_final"] >= endpoint_min
        ),
        "sustained_tail_z_order": (
            candidate_summary["z_band_order_tail_median"] >= tail_median_min
        ),
    }
    return {
        "reference": reference_summary,
        "candidate": candidate_summary,
        "candidate_metadata": meta,
        "metadata_gates": metadata_gates,
        "thresholds": {
            "expected_step": expected_step,
            "endpoint_z_order_min": endpoint_min,
            "tail_median_z_order_min": tail_median_min,
        },
        "gates": gates,
        "replicates_cell_a": all(gates.values()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--csv", type=Path)
    parser.add_argument("--expected-step", type=int, default=500_000)
    parser.add_argument("--endpoint-min", type=float, default=0.60)
    parser.add_argument("--tail-median-min", type=float, default=0.60)
    args = parser.parse_args()

    reference = z_band_series(args.reference)
    candidate = z_band_series(args.candidate)
    report = evaluate(
        reference, candidate, args.expected_step,
        args.endpoint_min, args.tail_median_min,
    )
    report.update({
        "reference_path": str(args.reference),
        "candidate_path": str(args.candidate),
    })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    rendered = json.dumps(report, indent=2)
    args.output.write_text(rendered + "\n")
    print(rendered)

    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        with args.csv.open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(("trajectory", "time_ms", "z_band_order"))
            for label, series in (("historical_cell_a", reference),
                                  ("corrected_control", candidate)):
                writer.writerows(zip(
                    [label] * len(series["time_ms"]),
                    series["time_ms"], series["z_order"],
                ))
    return 0 if report["replicates_cell_a"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
