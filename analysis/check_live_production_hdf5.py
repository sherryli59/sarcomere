#!/usr/bin/env python3
"""Check the latest fully visible frame of a running production trajectory."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import h5py
import numpy as np


CORE_DATASETS = (
    "state/current_step",
    "actin/center",
    "actin/direction",
    "myosin/center",
    "myosin/direction",
)


def band_order(x: np.ndarray, lx: float) -> float:
    wrapped = np.asarray(x) % lx
    return max(
        float(abs(np.exp(2j * np.pi * harmonic * wrapped / lx).mean()))
        for harmonic in (1, 2, 3, 4)
    )


def inspect(
    path: Path,
    actin_force_cap_fraction_max: float = 0.05,
    myosin_force_cap_fraction_max: float = 0.05,
) -> dict:
    with h5py.File(path, "r", locking=False) as handle:
        missing = [name for name in CORE_DATASETS if name not in handle]
        if missing:
            return {"trajectory": str(path), "healthy": False, "issues": [
                "missing datasets: " + ", ".join(missing)
            ]}

        frame_counts = {name: int(handle[name].shape[0]) for name in CORE_DATASETS}
        common_frames = min(frame_counts.values())
        issues = []
        if common_frames == 0:
            return {
                "trajectory": str(path), "healthy": False,
                "issues": ["no committed frames are visible"],
                "frame_counts": frame_counts,
            }
        if len(set(frame_counts.values())) != 1:
            issues.append("core dataset frame counts differ during read")

        if "meta/committed_frames" in handle:
            committed_frames = int(handle["meta/committed_frames"][0])
            if committed_frames != common_frames:
                issues.append(
                    f"committed_frames is {committed_frames}, while core datasets "
                    f"expose {common_frames} frames"
                )

        index = common_frames - 1
        steps = handle["state/current_step"][:common_frames, 0].astype(np.int64)
        actin_center = handle["actin/center"][index]
        actin_direction = handle["actin/direction"][index]
        myosin_center = handle["myosin/center"][index]
        myosin_direction = handle["myosin/direction"][index]
        meta = handle["meta"].attrs
        lx = float(meta.get("box_x", meta.get("Lx", 4.825)))
        dt = float(meta["dt"])

        arrays = (actin_center, actin_direction, myosin_center, myosin_direction)
        if not all(np.isfinite(values).all() for values in arrays):
            issues.append("NaN or infinite coordinate/direction detected")
        if len(steps) > 1 and np.any(np.diff(steps) <= 0):
            issues.append("saved steps are not strictly increasing")

        mechanics_targets = [
            ("catch_bond/local_link_identity", 1),
            ("catch_bond/local_link_mechanics_v2", 2),
        ]
        # The source-resolved motor-support table is emitted only by the
        # motor_transmitted_register classifier (enum value 8). Historical
        # productive_register trajectories legitimately omit it.
        if int(meta.get("aa_status_mode", -1)) == 8:
            mechanics_targets.append(("catch_bond/local_link_motor_support", 2))
        offset_specs = (
            (
                "catch_bond/local_link_frame_offsets",
                (("catch_bond/local_link_identity", 1),),
            ),
            (
                "catch_bond/local_link_mechanics_v2_frame_offsets",
                tuple(mechanics_targets),
            ),
        )
        for offset_name, targets in offset_specs:
            if offset_name not in handle:
                continue
            offsets = handle[offset_name]
            if offsets.shape[0] != common_frames:
                issues.append(
                    f"{offset_name} has {offsets.shape[0]} rows for "
                    f"{common_frames} committed frames"
                )
                continue
            offset_rows = offsets[:common_frames]
            if not np.array_equal(offset_rows[:, 0].astype(np.int64), steps):
                issues.append(f"{offset_name} steps do not match committed steps")
            latest = offset_rows[-1]
            count = int(latest[-1])
            for target_name, start_column in targets:
                if target_name not in handle:
                    issues.append(f"missing offset target dataset: {target_name}")
                    continue
                start = int(latest[start_column])
                if start < 0 or count < 0 or start + count > handle[target_name].shape[0]:
                    issues.append(
                        f"{offset_name} latest row points beyond {target_name}"
                    )

        actin_norm = np.linalg.norm(actin_direction, axis=1)
        myosin_norm = np.linalg.norm(myosin_direction, axis=1)
        direction_error = max(
            float(np.max(np.abs(actin_norm - 1.0))),
            float(np.max(np.abs(myosin_norm - 1.0))),
        )
        if direction_error > 1e-6:
            issues.append(f"filament direction norm error is {direction_error:.3g}")
        unit_actin_direction = actin_direction / np.maximum(
            actin_norm[:, None], np.finfo(float).tiny
        )
        actin_length = float(meta.get("actin_length", 1.0))
        barbed_plus = (
            actin_center[:, 0] + 0.5 * actin_length * unit_actin_direction[:, 0]
        )
        barbed_minus = (
            actin_center[:, 0] - 0.5 * actin_length * unit_actin_direction[:, 0]
        )
        z_band_order = max(
            band_order(barbed_plus, lx), band_order(barbed_minus, lx)
        )

        caps = {}
        cap_fraction_per_step = {}
        cap_fraction_max_post_ramp = {}
        cap_fraction_max_post_ramp_step = {}
        ramp_end_step = int(meta.get("interaction_ramp_start_step", 0)) + int(
            meta.get("interaction_ramp_steps", 0)
        )
        force_cap_limits = {
            "actin": actin_force_cap_fraction_max,
            "myosin": myosin_force_cap_fraction_max,
        }
        for species in ("actin", "myosin"):
            name = f"{species}/cumulative_force_cap_events"
            if name in handle and handle[name].shape[0] > index:
                values = handle[name][index, :, 0]
                caps[species] = float(np.sum(values))
                if index > 0 and handle[name].shape[0] > index - 1:
                    previous = handle[name][index - 1, :, 0]
                    step_delta = int(steps[index] - steps[index - 1])
                    if step_delta > 0 and len(values):
                        cap_fraction_per_step[species] = float(
                            np.sum(values - previous) / (step_delta * len(values))
                        )
                        if (
                            steps[index] >= ramp_end_step
                            and cap_fraction_per_step[species]
                            > force_cap_limits[species]
                        ):
                            issues.append(
                                f"{species} force caps affect "
                                f"{100.0 * cap_fraction_per_step[species]:.2f}% "
                                "of filaments per step in the latest interval"
                            )
                if index > 0:
                    history = handle[name][:index + 1, :, 0]
                    totals = np.sum(history, axis=1)
                    step_deltas = np.diff(steps)
                    fractions = np.diff(totals) / (step_deltas * history.shape[1])
                    post_indices = np.flatnonzero(steps[1:] >= ramp_end_step)
                    if len(post_indices):
                        local_max = int(np.argmax(fractions[post_indices]))
                        interval_index = int(post_indices[local_max])
                        maximum = float(fractions[interval_index])
                        maximum_step = int(steps[interval_index + 1])
                        cap_fraction_max_post_ramp[species] = maximum
                        cap_fraction_max_post_ramp_step[species] = maximum_step
                        if maximum > force_cap_limits[species]:
                            issues.append(
                                f"{species} force caps reached "
                                f"{100.0 * maximum:.2f}% of filaments per step "
                                f"in the post-ramp interval ending at step "
                                f"{maximum_step}"
                            )

        status2 = None
        name = "actin/aa_local_status2_count"
        if name in handle and handle[name].shape[0] > index:
            status2 = float(np.sum(handle[name][index, :, 0]) / 2.0)

        module = None
        model = meta.get("myomesin_model", b"")
        if isinstance(model, bytes):
            model = model.decode()
        module_expected = (
            float(meta.get("k_mm", 0.0)) > 0.0 and model != "disabled"
        )
        if module_expected and "myomesin/frame_offsets" not in handle:
            issues.append("active myomesin model is missing frame_offsets")
        if "myomesin/frame_offsets" in handle:
            offsets = handle["myomesin/frame_offsets"]
            if offsets.shape[0] != common_frames:
                issues.append(
                    f"myomesin/frame_offsets has {offsets.shape[0]} rows for "
                    f"{common_frames} committed frames"
                )
            elif not np.array_equal(
                offsets[:common_frames, 0].astype(np.int64), steps
            ):
                issues.append("myomesin/frame_offsets steps do not match committed steps")
            stable_module_frames = min(common_frames, int(offsets.shape[0]))
            if stable_module_frames:
                _, start, count = offsets[stable_module_frames - 1]
                start, count = int(start), int(count)
                available = int(handle["myomesin/bridge_state"].shape[0])
                if start < 0 or count < 0 or start + count > available:
                    issues.append("latest myomesin frame points beyond bridge-state rows")
                else:
                    rows = handle["myomesin/bridge_state"][start:start + count]
                    finite_rows = np.all(np.isfinite(rows), axis=1)
                    if not np.all(finite_rows):
                        issues.append("NaN or infinite myomesin bridge state detected")
                    topology_rows = rows[finite_rows]
                    pairs = [
                        tuple(sorted((int(row[2]), int(row[3]))))
                        for row in topology_rows
                    ]
                    degree = np.zeros(myosin_center.shape[0], dtype=int)
                    for first, second in pairs:
                        if first < 0 or second >= len(degree) or first == second:
                            issues.append("invalid myomesin endpoint index")
                            continue
                        degree[first] += 1
                        degree[second] += 1
                    duplicates = len(pairs) - len(set(pairs))
                    max_degree = int(degree.max(initial=0))
                    if duplicates:
                        issues.append(f"{duplicates} duplicate myomesin pair rows")
                    sites = int(meta.get("myomesin_sites_per_filament", 0))
                    if max_degree > sites:
                        issues.append(
                            f"myomesin degree {max_degree} exceeds configured valency {sites}"
                        )
                    force_max = float(meta.get("myomesin_force_max", np.inf))
                    observed_force_max = (
                        float(np.max(topology_rows[:, 7]))
                        if len(topology_rows) else None
                    )
                    if (
                        observed_force_max is not None
                        and observed_force_max > force_max + 1e-9
                    ):
                        issues.append(
                            f"myomesin force {observed_force_max:.6g} pN exceeds "
                            f"configured capacity {force_max:.6g} pN"
                        )
                    module = {
                        "bridge_count": len(rows),
                        "max_degree": max_degree,
                        "duplicate_pairs": duplicates,
                        "axial_mismatch_nm_median": (
                            float(np.median(np.abs(topology_rows[:, 5])) * 1e3)
                            if len(topology_rows) else None
                        ),
                        "radial_spacing_nm_median": (
                            float(np.median(topology_rows[:, 6]) * 1e3)
                            if len(topology_rows) else None
                        ),
                        "force_pn_max": observed_force_max,
                    }

        result = {
            "trajectory": str(path),
            "healthy": not issues,
            "issues": issues,
            "visible_frames": common_frames,
            "frame_counts": frame_counts,
            "latest_step": int(steps[-1]),
            "time_ms": float(steps[-1] * dt * 1e3),
            "finite_coordinates": all(np.isfinite(values).all() for values in arrays),
            "direction_norm_error_max": direction_error,
            "m_band_order": band_order(myosin_center[:, 0], lx),
            "z_band_order": z_band_order,
            "a_band_order": band_order(actin_center[:, 0], lx),
            "status2_links": status2,
            "force_cap_events": caps,
            "force_cap_fraction_per_step_latest": cap_fraction_per_step,
            "force_cap_fraction_per_step_max_post_ramp": (
                cap_fraction_max_post_ramp
            ),
            "force_cap_fraction_max_post_ramp_step": (
                cap_fraction_max_post_ramp_step
            ),
            "myomesin": module,
        }
        return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trajectory", type=Path)
    parser.add_argument("--require-healthy", action="store_true")
    parser.add_argument("--compact", action="store_true")
    parser.add_argument(
        "--actin-force-cap-fraction-max", type=float, default=0.05
    )
    parser.add_argument(
        "--myosin-force-cap-fraction-max", type=float, default=0.05
    )
    args = parser.parse_args()
    try:
        result = inspect(
            args.trajectory,
            args.actin_force_cap_fraction_max,
            args.myosin_force_cap_fraction_max,
        )
    except (OSError, RuntimeError, KeyError, ValueError) as error:
        result = {
            "trajectory": str(args.trajectory),
            "healthy": False,
            "transient_read_error": str(error),
        }
    result["checked_utc"] = datetime.now(timezone.utc).isoformat()
    print(json.dumps(
        result,
        indent=None if args.compact else 2,
        separators=(",", ":") if args.compact else None,
    ))
    return 1 if args.require_healthy and not result.get("healthy", False) else 0


if __name__ == "__main__":
    raise SystemExit(main())
