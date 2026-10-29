#!/usr/bin/env python3
"""Measure how explicit AM heads are distributed over actin partners.

The script compares matched trajectories frame-by-frame over their common
saved steps.  It distinguishes an occupied pair (at least one weak or strong
head) from a strong pair (at least one strong head), and reports both whole-
myosin and half-filament partner counts.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import h5py
import numpy as np


WEAK = 1
STRONG = 2


def committed_frames(handle: h5py.File) -> int:
    lengths = [
        handle["state/current_step"].shape[0],
        handle["state/am_head_state_current"].shape[0],
        handle["state/am_head_bound_actin_current"].shape[0],
    ]
    if "meta/committed_frames" in handle:
        lengths.append(int(np.asarray(handle["meta/committed_frames"])[0]))
    return min(lengths)


def saved_steps(handle: h5py.File) -> np.ndarray:
    count = committed_frames(handle)
    return np.asarray(handle["state/current_step"][:count], dtype=np.int64).reshape(-1)


def frame_indices_for_steps(all_steps: np.ndarray, chosen: np.ndarray) -> np.ndarray:
    indices = np.searchsorted(all_steps, chosen)
    if np.any(indices >= all_steps.size) or not np.array_equal(all_steps[indices], chosen):
        raise ValueError("chosen steps are not an exact subset of saved steps")
    return indices


def unique_count(values: np.ndarray) -> int:
    values = values[values >= 0]
    return int(np.unique(values).size) if values.size else 0


def summarize_frames(
    handle: h5py.File,
    frame_indices: np.ndarray,
    heads_per_half: int,
    n_myosins: int,
) -> dict[str, object]:
    sample_count = 0
    any_partners_sum = 0
    strong_partners_sum = 0
    any_half_partners_sum = 0
    strong_half_partners_sum = 0
    attached_heads_sum = 0
    strong_heads_sum = 0
    occupied_pairs_sum = 0
    strong_pairs_sum = 0
    occupied_half_pairs_sum = 0
    strong_half_pairs_sum = 0
    max_share_sum = 0.0
    hhi_sum = 0.0
    strong_active_myosin_count = 0
    partner_hist: Counter[int] = Counter()
    strong_partner_hist: Counter[int] = Counter()
    pair_strong_occupancy_hist: Counter[int] = Counter()
    pair_any_occupancy_hist: Counter[int] = Counter()

    states_ds = handle["state/am_head_state_current"]
    bound_ds = handle["state/am_head_bound_actin_current"]
    width = 2 * heads_per_half
    expected = n_myosins * width

    for frame in frame_indices:
        states = np.asarray(states_ds[int(frame)], dtype=np.int16).reshape(-1)
        bound = np.asarray(bound_ds[int(frame)], dtype=np.int32).reshape(-1)
        if states.size != expected or bound.size != expected:
            raise ValueError(
                f"head-array width mismatch: expected {expected}, "
                f"got state={states.size}, bound={bound.size}"
            )
        states = states.reshape(n_myosins, 2, heads_per_half)
        bound = bound.reshape(n_myosins, 2, heads_per_half)
        attached = ((states == WEAK) | (states == STRONG)) & (bound >= 0)
        strong = (states == STRONG) & (bound >= 0)

        for myosin in range(n_myosins):
            any_values = bound[myosin][attached[myosin]]
            strong_values = bound[myosin][strong[myosin]]
            any_partners = unique_count(any_values)
            strong_partners = unique_count(strong_values)
            partner_hist[any_partners] += 1
            strong_partner_hist[strong_partners] += 1
            any_partners_sum += any_partners
            strong_partners_sum += strong_partners
            attached_heads_sum += int(attached[myosin].sum())
            strong_count = int(strong[myosin].sum())
            strong_heads_sum += strong_count

            if any_values.size:
                counts = np.unique(any_values, return_counts=True)[1]
                occupied_pairs_sum += counts.size
                pair_any_occupancy_hist.update(int(value) for value in counts)
            if strong_values.size:
                counts = np.unique(strong_values, return_counts=True)[1]
                strong_pairs_sum += counts.size
                pair_strong_occupancy_hist.update(int(value) for value in counts)
                total = float(counts.sum())
                max_share_sum += float(counts.max()) / total
                hhi_sum += float(np.square(counts / total).sum())
                strong_active_myosin_count += 1

            for half in range(2):
                any_half_values = bound[myosin, half][attached[myosin, half]]
                strong_half_values = bound[myosin, half][strong[myosin, half]]
                any_half_partners_sum += unique_count(any_half_values)
                strong_half_partners_sum += unique_count(strong_half_values)
                if any_half_values.size:
                    occupied_half_pairs_sum += np.unique(any_half_values).size
                if strong_half_values.size:
                    strong_half_pairs_sum += np.unique(strong_half_values).size

        sample_count += n_myosins

    half_sample_count = 2 * sample_count
    strong_pair_heads = sum(k * v for k, v in pair_strong_occupancy_hist.items())
    any_pair_heads = sum(k * v for k, v in pair_any_occupancy_hist.items())
    strong_pair_count = sum(pair_strong_occupancy_hist.values())
    any_pair_count = sum(pair_any_occupancy_hist.values())
    return {
        "frames": int(frame_indices.size),
        "myosin_frame_samples": sample_count,
        "attached_heads_per_myosin": attached_heads_sum / sample_count,
        "strong_heads_per_myosin": strong_heads_sum / sample_count,
        "attached_actin_partners_per_myosin": any_partners_sum / sample_count,
        "strong_actin_partners_per_myosin": strong_partners_sum / sample_count,
        "attached_actin_partners_per_half": any_half_partners_sum / half_sample_count,
        "strong_actin_partners_per_half": strong_half_partners_sum / half_sample_count,
        "heads_per_occupied_pair": any_pair_heads / max(1, any_pair_count),
        "strong_heads_per_strong_pair": strong_pair_heads / max(1, strong_pair_count),
        "strong_heads_per_occupied_pair": strong_heads_sum / max(1, occupied_pairs_sum),
        "strong_pair_max_share_per_active_myosin": max_share_sum / max(1, strong_active_myosin_count),
        "strong_pair_hhi_per_active_myosin": hhi_sum / max(1, strong_active_myosin_count),
        "fraction_myosins_with_any_attachment": 1.0 - partner_hist[0] / sample_count,
        "fraction_myosins_with_strong_attachment": 1.0 - strong_partner_hist[0] / sample_count,
        "fraction_strong_pairs_with_ge2_heads": sum(
            value for key, value in pair_strong_occupancy_hist.items() if key >= 2
        ) / max(1, strong_pair_count),
        "fraction_strong_pairs_with_ge3_heads": sum(
            value for key, value in pair_strong_occupancy_hist.items() if key >= 3
        ) / max(1, strong_pair_count),
        "partner_count_hist": dict(sorted(partner_hist.items())),
        "strong_partner_count_hist": dict(sorted(strong_partner_hist.items())),
        "pair_strong_occupancy_hist": dict(sorted(pair_strong_occupancy_hist.items())),
        "pair_any_occupancy_hist": dict(sorted(pair_any_occupancy_hist.items())),
    }


def metadata(handle: h5py.File) -> tuple[float, int, int]:
    attrs = handle["meta"].attrs
    dt = float(attrs["dt"])
    heads_per_half = int(attrs["am_heads_per_half"])
    n_myosins = int(attrs.get("n_myosins", handle["myosin/center"].shape[1]))
    return dt, heads_per_half, n_myosins


def parse_pair(value: str) -> tuple[str, Path, Path]:
    parts = value.split("=", 1)
    if len(parts) != 2 or "," not in parts[1]:
        raise argparse.ArgumentTypeError("expected SEED=BASELINE,KSOFF27")
    baseline, ksoff = parts[1].split(",", 1)
    return parts[0], Path(baseline), Path(ksoff)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pair", action="append", required=True, type=parse_pair,
        help="matched trajectory pair as SEED=BASELINE,KSOFF27",
    )
    parser.add_argument("--burn-in", type=float, default=0.5)
    parser.add_argument("--tail", type=float, default=0.5)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    args = parser.parse_args()

    rows: list[dict[str, object]] = []
    nested: dict[str, object] = {"seeds": {}}
    for seed, baseline_path, ksoff_path in args.pair:
        with h5py.File(baseline_path, "r", locking=False) as baseline, h5py.File(
            ksoff_path, "r", locking=False
        ) as ksoff:
            dt0, heads0, myosins0 = metadata(baseline)
            dt1, heads1, myosins1 = metadata(ksoff)
            if (dt0, heads0, myosins0) != (dt1, heads1, myosins1):
                raise ValueError(f"metadata mismatch for seed {seed}")
            steps0 = saved_steps(baseline)
            steps1 = saved_steps(ksoff)
            common = np.intersect1d(steps0, steps1, assume_unique=True)
            times = common.astype(float) * dt0
            if common.size == 0:
                raise ValueError(f"no common saved steps for seed {seed}")
            windows = {
                "post_burn": common[times >= args.burn_in],
                "tail": common[times >= max(args.burn_in, times[-1] - args.tail)],
            }
            seed_result: dict[str, object] = {
                "baseline": str(baseline_path),
                "ksoff27": str(ksoff_path),
                "common_start": float(times[0]),
                "common_end": float(times[-1]),
                "windows": {},
            }
            for window, chosen in windows.items():
                if chosen.size == 0:
                    raise ValueError(f"empty {window} window for seed {seed}")
                index0 = frame_indices_for_steps(steps0, chosen)
                index1 = frame_indices_for_steps(steps1, chosen)
                summaries = {
                    "baseline": summarize_frames(
                        baseline, index0, heads0, myosins0
                    ),
                    "ksoff27": summarize_frames(
                        ksoff, index1, heads0, myosins0
                    ),
                }
                seed_result["windows"][window] = {
                    "start": float(chosen[0] * dt0),
                    "end": float(chosen[-1] * dt0),
                    **summaries,
                }
                scalar_keys = [
                    key for key, value in summaries["baseline"].items()
                    if isinstance(value, (int, float)) and key not in {
                        "frames", "myosin_frame_samples"
                    }
                ]
                for arm, summary in summaries.items():
                    row: dict[str, object] = {
                        "seed": seed,
                        "window": window,
                        "arm": arm,
                        "start": float(chosen[0] * dt0),
                        "end": float(chosen[-1] * dt0),
                    }
                    row.update({key: summary[key] for key in scalar_keys})
                    rows.append(row)
            nested["seeds"][seed] = seed_result

    scalar_columns = [key for key in rows[0] if key not in {"seed", "window", "arm", "start", "end"}]
    aggregate: dict[str, object] = {}
    for window in ("post_burn", "tail"):
        aggregate[window] = {}
        for arm in ("baseline", "ksoff27"):
            chosen_rows = [row for row in rows if row["window"] == window and row["arm"] == arm]
            aggregate[window][arm] = {
                key: {
                    "mean": float(np.mean([float(row[key]) for row in chosen_rows])),
                    "sd_across_seeds": float(np.std(
                        [float(row[key]) for row in chosen_rows], ddof=1
                    )),
                }
                for key in scalar_columns
            }
        aggregate[window]["paired_change"] = {
            key: {
                "mean_delta": float(np.mean([
                    float(next(row[key] for row in rows if row["seed"] == seed and row["window"] == window and row["arm"] == "ksoff27"))
                    - float(next(row[key] for row in rows if row["seed"] == seed and row["window"] == window and row["arm"] == "baseline"))
                    for seed, _, _ in args.pair
                ])),
                "mean_ratio": float(np.mean([
                    float(next(row[key] for row in rows if row["seed"] == seed and row["window"] == window and row["arm"] == "ksoff27"))
                    / max(1e-15, float(next(row[key] for row in rows if row["seed"] == seed and row["window"] == window and row["arm"] == "baseline")))
                    for seed, _, _ in args.pair
                ])),
            }
            for key in scalar_columns
        }
    nested["aggregate"] = aggregate

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(nested, indent=2, sort_keys=True) + "\n")
    with args.output_csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(aggregate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
