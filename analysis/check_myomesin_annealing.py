#!/usr/bin/env python3
"""Check M-band annealing relative to the instant each module first binds."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np


def parse_pairs(value: str) -> list[tuple[int, int]]:
    pairs = []
    for token in value.split(","):
        first, second = (int(part) for part in token.split("-", 1))
        pairs.append(tuple(sorted((first, second))))
    return pairs


def evaluate(
    trajectory: Path,
    requested_pairs: list[tuple[int, int]] | None,
    final_window_ms: float,
    minimum_reduction_fraction: float,
) -> dict:
    with h5py.File(trajectory, "r") as handle:
        if "myomesin/frame_offsets" not in handle or \
                "myomesin/bridge_state" not in handle:
            raise ValueError("trajectory has no sparse myomesin module state")
        offsets = handle["myomesin/frame_offsets"][:]
        state = handle["myomesin/bridge_state"][:]
        events = handle["myomesin/events"][:] \
            if "myomesin/events" in handle else np.empty((0, 9))
        dt = float(handle["meta"].attrs["dt"])
        force_max = float(handle["meta"].attrs["myomesin_force_max"])
        axial_k = float(handle["meta"].attrs["k_mm"])
        radial_k = float(handle["meta"].attrs["myomesin_radial_k"])

    observations: dict[tuple[int, int], list[tuple[int, float, float, float]]] = {}
    for step, start, count in offsets:
        rows = state[int(start):int(start + count)]
        for row in rows:
            pair = tuple(sorted((int(row[2]), int(row[3]))))
            observations.setdefault(pair, []).append(
                (int(step), 1e3 * abs(float(row[5])),
                 1e3 * float(row[6]), float(row[7]))
            )

    formation_events: dict[tuple[int, int], tuple[int, float]] = {}
    for row in events:
        if int(row[1]) != 1:
            continue
        pair = tuple(sorted((int(row[3]), int(row[4]))))
        formation_events.setdefault(
            pair, (int(row[0]), 1e3 * abs(float(row[5])))
        )

    pairs = requested_pairs or sorted(observations)
    final_step = int(offsets[-1, 0])
    final_start = final_step - int(round(final_window_ms * 1e-3 / dt))
    results = []
    for pair in pairs:
        values = observations.get(pair, [])
        if not values:
            results.append({
                "pair": list(pair),
                "present": False,
                "pass": False,
                "reason": "module never formed",
            })
            continue
        if pair in formation_events:
            formation_step, formation_nm = formation_events[pair]
            formation_source = "association_event"
        else:
            formation_step, formation_nm, _, _ = values[0]
            formation_source = "first_saved_frame_fallback"
        final_values = [row for row in values if row[0] >= final_start]
        if not final_values:
            results.append({
                "pair": list(pair),
                "present": True,
                "pass": False,
                "reason": "module absent throughout final window",
                "formation_step": formation_step,
                "formation_mismatch_nm": formation_nm,
            })
            continue
        final_median_nm = float(np.median([row[1] for row in final_values]))
        reduction_nm = formation_nm - final_median_nm
        reduction_fraction = reduction_nm / formation_nm if formation_nm else 0.0
        passed = reduction_nm > 0.0 and \
            reduction_fraction >= minimum_reduction_fraction
        results.append({
            "pair": list(pair),
            "present": True,
            "pass": passed,
            "formation_step": formation_step,
            "formation_time_ms": formation_step * dt * 1e3,
            "formation_mismatch_nm": formation_nm,
            "formation_source": formation_source,
            "final_window_median_nm": final_median_nm,
            "reduction_nm": reduction_nm,
            "reduction_fraction": reduction_fraction,
            "final_radial_median_nm": float(np.median(
                [row[2] for row in final_values]
            )),
            "final_force_median_pn": float(np.median(
                [row[3] for row in final_values]
            )),
            "final_samples": len(final_values),
        })

    return {
        "trajectory": str(trajectory),
        "axial_k_pn_per_um": axial_k,
        "radial_k_pn_per_um": radial_k,
        "module_force_max_pn": force_max,
        "final_window_ms": final_window_ms,
        "minimum_reduction_fraction": minimum_reduction_fraction,
        "pass": bool(results) and all(item["pass"] for item in results),
        "pairs": results,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Require each sparse M-M module to reduce its pair mismatch "
            "relative to the frame in which that module first appeared."
        )
    )
    parser.add_argument("trajectory", type=Path)
    parser.add_argument(
        "--pairs", type=parse_pairs,
        help="comma-separated expected pairs, for example 0-2,1-3",
    )
    parser.add_argument("--final-window-ms", type=float, default=10.0)
    parser.add_argument(
        "--minimum-reduction-fraction", type=float, default=0.10,
        help="required fractional reduction from formation mismatch (default 0.10)",
    )
    parser.add_argument("--json", type=Path, help="also write the result as JSON")
    args = parser.parse_args()
    if args.final_window_ms <= 0.0:
        parser.error("--final-window-ms must be positive")
    if not 0.0 <= args.minimum_reduction_fraction < 1.0:
        parser.error("--minimum-reduction-fraction must be in [0, 1)")

    result = evaluate(
        args.trajectory, args.pairs, args.final_window_ms,
        args.minimum_reduction_fraction,
    )
    if args.json:
        args.json.write_text(json.dumps(result, indent=2) + "\n")

    print(
        f"module k_axial={result['axial_k_pn_per_um']:g} pN/um, "
        f"k_radial={result['radial_k_pn_per_um']:g} pN/um, "
        f"force_max={result['module_force_max_pn']:g} pN"
    )
    for item in result["pairs"]:
        pair = "-".join(str(index) for index in item["pair"])
        if not item["present"] or "final_window_median_nm" not in item:
            print(f"M{pair}: FAIL ({item['reason']})")
            continue
        print(
            f"M{pair}: {item['formation_mismatch_nm']:.2f} -> "
            f"{item['final_window_median_nm']:.2f} nm; "
            f"change {item['reduction_nm']:+.2f} nm "
            f"({100.0 * item['reduction_fraction']:+.1f}%): "
            f"{'PASS' if item['pass'] else 'FAIL'}"
        )
    print("PASS" if result["pass"] else "FAIL")
    return 0 if result["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
