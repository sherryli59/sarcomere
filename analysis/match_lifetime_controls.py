#!/usr/bin/env python3
"""Rank constant-lifetime controls by matched CB2 occupancy and turnover.

The engine stores cumulative status-2 binding and unbinding event counts per
actin.  Their system sums count each pair transition once because the engine
assigns half of an event to each participating actin.  This script measures a
tail-window CB2 occupancy and steady turnover rate for every run, aggregates
replicas by protocol tag, and ranks constant-lifetime protocols against a
catch-bond reference.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import h5py
import numpy as np


TAG_RE = re.compile(r"_M\d+_(?P<tag>.+)_(?P<rep>r\d+)_L")


def committed_frames(handle: h5py.File) -> int:
    available = handle["/actin/cb_status"].shape[0]
    if "/meta/committed_frames" in handle:
        committed = int(np.asarray(handle["/meta/committed_frames"])[0])
        if committed > 0:
            return min(committed, available)
    return available


def infer_saved_dt(handle: h5py.File, simulation_dt: float) -> float:
    if "/state/current_step" in handle and handle["/state/current_step"].shape[0] >= 2:
        steps = np.asarray(handle["/state/current_step"][:2]).reshape(-1)
        delta = int(steps[1] - steps[0])
        if delta > 0:
            return delta * simulation_dt
    raise ValueError("cannot infer saved-frame spacing from /state/current_step")


def protocol_identity(path: Path) -> tuple[str, str]:
    match = TAG_RE.search(path.name)
    if not match:
        return path.stem, "unknown"
    return match.group("tag"), match.group("rep")


def cumulative_total(dataset: h5py.Dataset, frame: int) -> float:
    return float(np.asarray(dataset[frame], dtype=float).sum())


def completed_lifetime_stats(handle: h5py.File) -> tuple[int, float, float, float]:
    path = "/catch_bond/completed_lifetimes"
    if path not in handle or handle[path].shape[0] == 0:
        return 0, np.nan, np.nan, np.nan
    values = np.asarray(handle[path], dtype=float)
    values = values[np.isfinite(values) & (values >= 0)]
    if values.size == 0:
        return 0, np.nan, np.nan, np.nan
    return (
        int(values.size),
        float(np.mean(values)),
        float(np.median(values)),
        float(np.percentile(values, 90)),
    )


def measure_run(
    path: Path,
    tail_fraction: float,
    stride: int,
    simulation_dt: float,
) -> dict[str, object]:
    tag, replica = protocol_identity(path)
    with h5py.File(path, "r") as handle:
        required = (
            "/actin/cb_status",
            "/actin/cumulative_chem_binding_events",
            "/actin/cumulative_chem_unbinding_events",
        )
        missing = [name for name in required if name not in handle]
        if missing:
            raise KeyError(f"{path.name}: missing {', '.join(missing)}")

        n_frames = committed_frames(handle)
        if n_frames < 2:
            raise ValueError(f"{path.name}: need at least two committed frames")
        n_actins = int(handle["/actin/cb_status"].shape[1])
        first = max(0, int(np.floor((1.0 - tail_fraction) * n_frames)))
        last = n_frames - 1
        frame_idx = np.arange(first, n_frames, max(1, stride), dtype=int)
        if frame_idx[-1] != last:
            frame_idx = np.append(frame_idx, last)

        cb2_counts = np.empty(frame_idx.size, dtype=float)
        cb_status = handle["/actin/cb_status"]
        for idx, frame in enumerate(frame_idx):
            cb2_counts[idx] = np.count_nonzero(np.asarray(cb_status[frame]).reshape(-1) == 2)

        saved_dt = infer_saved_dt(handle, simulation_dt)
        duration = (last - first) * saved_dt
        binding = handle["/actin/cumulative_chem_binding_events"]
        unbinding = handle["/actin/cumulative_chem_unbinding_events"]
        binding_events = cumulative_total(binding, last) - cumulative_total(binding, first)
        unbinding_events = cumulative_total(unbinding, last) - cumulative_total(unbinding, first)
        # At steady state these rates agree. Their mean is a lower-noise count of
        # complete turnover cycles without double-counting bind+unbind as two cycles.
        turnover_events = 0.5 * (binding_events + unbinding_events)
        turnover_rate = turnover_events / max(duration, np.finfo(float).eps)
        turnover_rate_per_actin = turnover_rate / max(n_actins, 1)
        occupancy = float(np.mean(cb2_counts) / max(n_actins, 1))
        turnover_per_cb2 = turnover_rate / max(float(np.mean(cb2_counts)), 1.0)
        n_lifetimes, lifetime_mean, lifetime_median, lifetime_p90 = completed_lifetime_stats(handle)

    return {
        "filename": path.name,
        "path": str(path.resolve()),
        "tag": tag,
        "replica": replica,
        "n_frames": n_frames,
        "tail_first_frame": first,
        "tail_last_frame": last,
        "tail_duration": duration,
        "n_actins": n_actins,
        "cb2_occupancy": occupancy,
        "mean_cb2_count": float(np.mean(cb2_counts)),
        "binding_events": binding_events,
        "unbinding_events": unbinding_events,
        "turnover_rate": turnover_rate,
        "turnover_rate_per_actin": turnover_rate_per_actin,
        "turnover_rate_per_cb2": turnover_per_cb2,
        "n_completed_lifetimes": n_lifetimes,
        "completed_lifetime_mean": lifetime_mean,
        "completed_lifetime_median": lifetime_median,
        "completed_lifetime_p90": lifetime_p90,
    }


def mean_std(values: list[float]) -> tuple[float, float]:
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if array.size == 0:
        return np.nan, np.nan
    return float(np.mean(array)), float(np.std(array, ddof=1)) if array.size > 1 else 0.0


def aggregate(records: list[dict[str, object]]) -> list[dict[str, object]]:
    by_tag: dict[str, list[dict[str, object]]] = {}
    for record in records:
        by_tag.setdefault(str(record["tag"]), []).append(record)

    metrics = (
        "cb2_occupancy",
        "turnover_rate_per_actin",
        "turnover_rate_per_cb2",
        "completed_lifetime_mean",
        "completed_lifetime_median",
        "completed_lifetime_p90",
    )
    summaries: list[dict[str, object]] = []
    for tag, rows in sorted(by_tag.items()):
        summary: dict[str, object] = {"tag": tag, "n_replicates": len(rows)}
        for metric in metrics:
            mean, std = mean_std([float(row[metric]) for row in rows])
            summary[f"{metric}_mean"] = mean
            summary[f"{metric}_std"] = std
        summaries.append(summary)
    return summaries


def rank_matches(
    summaries: list[dict[str, object]], reference_pattern: str
) -> list[dict[str, object]]:
    references = [row for row in summaries if reference_pattern in str(row["tag"])]
    if len(references) != 1:
        names = ", ".join(str(row["tag"]) for row in references) or "none"
        raise ValueError(
            f"reference pattern {reference_pattern!r} must select one tag; selected {names}"
        )
    reference = references[0]
    occ_ref = float(reference["cb2_occupancy_mean"])
    turnover_ref = float(reference["turnover_rate_per_actin_mean"])
    occ_scale = max(abs(occ_ref), 1e-12)
    turnover_scale = max(abs(turnover_ref), 1e-12)

    rankings: list[dict[str, object]] = []
    for row in summaries:
        if row is reference:
            continue
        occ_delta = (float(row["cb2_occupancy_mean"]) - occ_ref) / occ_scale
        turnover_delta = (
            float(row["turnover_rate_per_actin_mean"]) - turnover_ref
        ) / turnover_scale
        rankings.append(
            {
                "reference_tag": reference["tag"],
                "candidate_tag": row["tag"],
                "match_score": float(np.hypot(occ_delta, turnover_delta)),
                "relative_occupancy_error": occ_delta,
                "relative_turnover_error": turnover_delta,
                "candidate_cb2_occupancy": row["cb2_occupancy_mean"],
                "reference_cb2_occupancy": occ_ref,
                "candidate_turnover_rate_per_actin": row["turnover_rate_per_actin_mean"],
                "reference_turnover_rate_per_actin": turnover_ref,
                "candidate_completed_lifetime_median": row["completed_lifetime_median_mean"],
                "reference_completed_lifetime_median": reference[
                    "completed_lifetime_median_mean"
                ],
            }
        )
    return sorted(rankings, key=lambda row: float(row["match_score"]))


def write_csv(path: Path, records: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not records:
        raise ValueError(f"no records to write to {path}")
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def collect_files(inputs: list[Path]) -> list[Path]:
    files: list[Path] = []
    for item in inputs:
        if item.is_dir():
            files.extend(sorted(item.glob("*.h5")))
        else:
            files.append(item)
    return files


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--reference-pattern", default="catch")
    parser.add_argument("--tail-fraction", type=float, default=0.25)
    parser.add_argument("--stride", type=int, default=10)
    parser.add_argument("--simulation-dt", type=float, default=1e-5)
    args = parser.parse_args()

    if not 0.0 < args.tail_fraction <= 1.0:
        parser.error("--tail-fraction must lie in (0, 1]")
    if args.stride < 1:
        parser.error("--stride must be positive")
    if args.simulation_dt <= 0.0:
        parser.error("--simulation-dt must be positive")

    records = [
        measure_run(path, args.tail_fraction, args.stride, args.simulation_dt)
        for path in collect_files(args.inputs)
    ]
    summaries = aggregate(records)
    rankings = rank_matches(summaries, args.reference_pattern)
    write_csv(args.output_dir / "run_metrics.csv", records)
    write_csv(args.output_dir / "tag_summary.csv", summaries)
    write_csv(args.output_dir / "match_ranking.csv", rankings)

    print(f"wrote {len(records)} run metrics to {args.output_dir}")
    print("best occupancy/turnover matches:")
    for row in rankings[:5]:
        print(
            f"  {row['candidate_tag']}: score={float(row['match_score']):.4f}, "
            f"occupancy error={float(row['relative_occupancy_error']):+.1%}, "
            f"turnover error={float(row['relative_turnover_error']):+.1%}"
        )


if __name__ == "__main__":
    main()
