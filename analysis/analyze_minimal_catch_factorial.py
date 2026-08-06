#!/usr/bin/env python3
"""Analyze production-like minimal AM/ACTN2 persistence and annealing tests."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/sarcomere_matplotlib")

import h5py
import matplotlib.pyplot as plt
import numpy as np


def scalar_series(handle: h5py.File, path: str) -> np.ndarray:
    values = np.asarray(handle[path], dtype=float)
    return values.reshape((len(values), -1)).sum(axis=1)


def authoritative_head_series(handle: h5py.File):
    """Count conserved weak/strong heads, including the seeded frame zero."""
    if ("/state/am_head_state_current" in handle and
            "/state/am_head_bound_actin_current" in handle):
        state = np.asarray(handle["/state/am_head_state_current"], dtype=int)
        bound = np.asarray(handle["/state/am_head_bound_actin_current"], dtype=int)
        occupied = bound >= 0
        axes = tuple(range(1, state.ndim))
        return (np.count_nonzero(occupied & (state == 1), axis=axes),
                np.count_nonzero(occupied & (state == 2), axis=axes))
    return (scalar_series(handle, "/actin/am_weak_count"),
            scalar_series(handle, "/actin/am_strong_count"))


def neighbor_head_occupancy_maxima(handle: h5py.File) -> tuple[int, int]:
    """Maximum saved weak+strong and strong occupancy per actin/myosin half."""
    if ("/state/am_head_state_current" not in handle or
            "/state/am_head_bound_actin_current" not in handle):
        return 0, 0
    state = np.asarray(handle["/state/am_head_state_current"], dtype=int)
    bound = np.asarray(
        handle["/state/am_head_bound_actin_current"], dtype=int)
    heads_per_half = int(handle["/meta"].attrs.get("am_heads_per_half", 1))
    max_occupied = 0
    max_strong = 0
    for frame in range(state.shape[0]):
        for myosin in range(state.shape[1]):
            for start in (0, heads_per_half):
                stop = start + heads_per_half
                local_bound = bound[frame, myosin, start:stop]
                local_state = state[frame, myosin, start:stop]
                for actin_index in np.unique(local_bound[local_bound >= 0]):
                    same_actin = local_bound == actin_index
                    max_occupied = max(max_occupied, int(np.count_nonzero(
                        same_actin & np.isin(local_state, (1, 2)))))
                    max_strong = max(max_strong, int(np.count_nonzero(
                        same_actin & (local_state == 2))))
    return max_occupied, max_strong


def local_series(handle: h5py.File):
    offsets = np.asarray(handle["/catch_bond/local_link_frame_offsets"], dtype=int)
    identity = np.asarray(handle["/catch_bond/local_link_identity"], dtype=int)
    mechanics = np.asarray(handle["/catch_bond/local_link_mechanics"], dtype=float)
    aa1 = np.zeros(len(offsets)); aa2 = np.zeros(len(offsets))
    bilateral = np.zeros(len(offsets)); support = np.zeros(len(offsets))
    lifetime = np.zeros(len(offsets)); occupancy = np.zeros(len(offsets))
    saved_detachment_hazard = np.zeros(len(offsets))
    for frame, (_, start, count) in enumerate(offsets):
        status = identity[start:start + count, 4]
        rows = mechanics[start:start + count]
        occupancy[frame] = count
        aa1[frame] = np.count_nonzero(status == 1)
        aa2[frame] = np.count_nonzero(status == 2)
        if count:
            bilateral[frame] = np.max(rows[:, 7])
            support[frame] = np.max(rows[:, 8])
            lifetime[frame] = np.max(rows[:, 10])
            # Column 10 stores each individual link's instantaneous lifetime.
            # KMC detachments are independent, so the frame hazard is the sum
            # of 1/tau over links, not occupancy divided by the longest tau.
            valid_lifetime = rows[:, 10][rows[:, 10] > 0]
            saved_detachment_hazard[frame] = np.sum(1.0 / valid_lifetime)
    return (offsets[:, 0], aa1, aa2, occupancy, bilateral, support, lifetime,
            saved_detachment_hazard)


def weak_dwell_statistics(events: np.ndarray, dt: float, end_step: int):
    starts: dict[tuple[int, int], int] = {}
    exposure, outcomes = [], []
    for row in events[np.argsort(events[:, 2])]:
        key = (int(row[1]), int(row[10]))
        step, old, new, reason = (int(row[2]), int(row[3]),
                                  int(row[4]), int(row[5]))
        if new == 1:
            starts[key] = step
        elif old == 1 and key in starts:
            start = starts.pop(key)
            exposure.append((step - start) * dt)
            # Geometric/neighbor ruptures are right-censored mechanical exits,
            # not realizations of the configured weak chemical hazard.
            if reason in (1, 2):
                outcomes.append(new)
    exposure.extend((end_step - start) * dt for start in starts.values())
    exposure = np.asarray(exposure, dtype=float)
    outcomes = np.asarray(outcomes, dtype=int)
    return (
        float(np.sum(exposure) / len(outcomes)) if len(outcomes) else np.nan,
        float(np.mean(outcomes == 2)) if len(outcomes) else np.nan,
        len(outcomes),
    )


def balanced_geometry(handle: h5py.File, actin_center: np.ndarray,
                      actin_direction: np.ndarray):
    """Whole-cell registration metrics for the periodic counterloaded unit."""
    myosin_center = np.asarray(handle["/myosin/center"], dtype=float)
    myosin_direction = np.asarray(handle["/myosin/direction"], dtype=float)
    attrs = handle["/meta"].attrs
    box = np.array([float(attrs.get("box_x", 5.32)),
                    float(attrs.get("box_y", 3.0)),
                    float(attrs.get("box_z", 3.0))])
    periodic = np.array([bool(attrs.get("periodic_x", 1)),
                         bool(attrs.get("periodic_y", 0)),
                         bool(attrs.get("periodic_z", 0))])
    actin_length = float(attrs.get("actin_length", 1.0))
    myosin_length = float(attrs.get("myosin_length", 1.5))
    bare_half = (0.5 * float(attrs.get("am_bare_zone_fraction", 0.10)) *
                 myosin_length)
    outer = 0.5 * myosin_length
    intended = ((0, 0, -1), (2, 0, 1), (1, 1, 1), (3, 1, -1))
    overlap = np.zeros((len(actin_center), len(intended)))
    for frame in range(len(actin_center)):
        for column, (actin_index, myosin_index, half_sign) in enumerate(intended):
            axis = myosin_direction[frame, myosin_index]
            axis /= max(np.linalg.norm(axis), 1e-30)
            displacement = (actin_center[frame, actin_index] -
                            myosin_center[frame, myosin_index])
            for dim in range(3):
                if periodic[dim] and box[dim] > 0:
                    displacement[dim] -= box[dim] * np.rint(
                        displacement[dim] / box[dim])
            axial = float(displacement @ axis)
            half_extent = 0.5 * actin_length * abs(
                float(actin_direction[frame, actin_index] @ axis))
            actin_min, actin_max = axial - half_extent, axial + half_extent
            head_min, head_max = ((-outer, -bare_half) if half_sign < 0
                                  else (bare_half, outer))
            overlap[frame, column] = max(
                0.0, min(actin_max, head_max) - max(actin_min, head_min))
    separation_vector = myosin_center[:, 0] - myosin_center[:, 1]
    for dim in range(3):
        if periodic[dim] and box[dim] > 0:
            separation_vector[:, dim] -= box[dim] * np.rint(
                separation_vector[:, dim] / box[dim])
    separation = np.linalg.norm(separation_vector, axis=1)
    status = np.asarray(handle["/actin/cb_status"], dtype=int)[:, :, 0]
    return overlap, separation, status


def analyze(path: Path) -> tuple[dict, dict]:
    variant = path.stem.split("_rep")[0]
    replica = int(path.stem.rsplit("rep", 1)[1])
    with h5py.File(path, "r", locking=False) as handle:
        dt = float(handle["/meta"].attrs["dt"])
        (steps, aa1, aa2, occupancy, bilateral, support, lifetime,
         saved_detachment_hazard) = local_series(handle)
        time = steps * dt
        weak, strong = authoritative_head_series(handle)
        max_neighbor_occupancy, max_neighbor_strong = (
            neighbor_head_occupancy_maxima(handle))
        calcium = np.asarray(handle["/actin/am_calcium"], dtype=float).reshape(len(time), -1).mean(1)
        q = np.asarray(handle["/actin/am_q"], dtype=float).reshape(len(time), -1).mean(1)
        directions = np.asarray(handle["/actin/direction"], dtype=float)
        alignment = np.abs(np.sum(directions[:, 0] * directions[:, 1], axis=1))
        angle = np.degrees(np.arccos(np.clip(alignment, 0, 1)))
        centers = np.asarray(handle["/actin/center"], dtype=float)
        com = centers.mean(axis=1)
        com_drift = np.linalg.norm(com - com[0], axis=1)
        events = np.asarray(handle.get("/actin_myo/head_transitions", np.empty((0, 11))))
        weak_dwell, weak_to_strong, n_weak_dwells = weak_dwell_statistics(
            events, dt, int(steps[-1]))
        local_events = np.asarray(
            handle.get("/catch_bond/local_link_events_v2", np.empty((0, 21))))
        aa_detach = int(np.count_nonzero(local_events[:, 4] == 2)) if len(local_events) else 0
        frame_dt = np.diff(time, prepend=time[0])
        am_detach_observed = int(np.count_nonzero(
            np.isin(events[:, 5].astype(int), (1, 3)))) if len(events) else 0
        capacity_events = (events[events[:, 5].astype(int) == 6]
                           if len(events) else np.empty((0, 11)))
        am_detach_expected = np.nan
        if "/state/am_head_detachment_rate_current" in handle:
            rates = np.asarray(
                handle["/state/am_head_detachment_rate_current"], dtype=float)
            am_detach_expected = float(np.sum(
                rates * frame_dt.reshape((-1,) + (1,) * (rates.ndim - 1))))
        expected_detach = float(np.sum(saved_detachment_hazard * frame_dt))
        force_caps = scalar_series(handle, "/actin/cumulative_force_cap_events")
        torque_caps = scalar_series(handle, "/actin/cumulative_torque_cap_events")
        rotation_caps = (scalar_series(
            handle, "/actin/cumulative_rotation_cap_events")
            if "/actin/cumulative_rotation_cap_events" in handle
            else np.zeros_like(force_caps))
        if "/myosin/cumulative_force_cap_events" in handle:
            force_caps += scalar_series(handle, "/myosin/cumulative_force_cap_events")
            torque_caps += scalar_series(handle, "/myosin/cumulative_torque_cap_events")
            if "/myosin/cumulative_rotation_cap_events" in handle:
                rotation_caps += scalar_series(
                    handle, "/myosin/cumulative_rotation_cap_events")
        cap_denominator = max(1.0, float(steps[-1]) *
                              (centers.shape[1] + handle["/myosin/center"].shape[1]))
        forced = np.asarray(handle.get("/actin_myo/forced_ruptures_v2", np.empty((0, 19))))
        row = {
            "path": str(path), "variant": variant, "replica": replica,
            "duration_s": float(time[-1]),
            "link_occupancy_fraction": float(np.mean(occupancy > 0)),
            "minimum_aa_links": int(np.min(occupancy)),
            "zero_link_frame_fraction": float(np.mean(occupancy == 0)),
            "status2_frame_fraction": float(np.mean(aa2 > 0)),
            "mean_aa1": float(np.mean(aa1)), "mean_aa2": float(np.mean(aa2)),
            "mean_weak_am": float(np.mean(weak)),
            "mean_strong_am": float(np.mean(strong)),
            "max_neighbor_head_occupancy": max_neighbor_occupancy,
            "max_neighbor_strong_occupancy": max_neighbor_strong,
            "am_capacity_detachments": len(capacity_events),
            "am_capacity_strong_detachments": int(np.count_nonzero(
                capacity_events[:, 3].astype(int) == 2))
                if len(capacity_events) else 0,
            "peak_bilateral_arm_load": float(np.max(bilateral)),
            "peak_productive_support": float(np.max(support)),
            "initial_misalignment_deg": float(angle[0]),
            "final_misalignment_deg": float(angle[-1]),
            "maximum_misalignment_deg": float(np.max(angle)),
            "misalignment_change_deg": float(angle[-1] - angle[0]),
            "bundle_com_drift_um": float(com_drift[-1]),
            "weak_chemical_timescale_s": weak_dwell,
            "weak_to_strong_fraction": weak_to_strong,
            "weak_chemical_event_count": n_weak_dwells,
            "aa_detach_observed": aa_detach,
            "aa_detach_expected_saved_hazard": expected_detach,
            "forced_am_ruptures": len(forced),
            "am_chemical_detach_observed": am_detach_observed,
            "am_chemical_detach_expected_saved_hazard": am_detach_expected,
            "am_chemical_detach_observed_over_expected": (
                float(am_detach_observed / am_detach_expected)
                if am_detach_expected > 0 else np.nan),
            "force_cap_events": float(force_caps[-1]),
            "torque_cap_events": float(torque_caps[-1]),
            "rotation_cap_events": float(rotation_caps[-1]),
            "force_cap_fraction": float(force_caps[-1] / cap_denominator),
            "torque_cap_fraction": float(torque_caps[-1] / cap_denominator),
            "rotation_cap_fraction": float(rotation_caps[-1] / cap_denominator),
            "finite": bool(all(np.isfinite(x).all() for x in
                               (centers, directions, weak, strong, bilateral))),
        }
        if variant == "balanced" and centers.shape[1] >= 4:
            intended_overlap, myosin_separation, actin_status = balanced_geometry(
                handle, centers, directions)
            row.update({
                "minimum_intended_am_overlap": float(np.min(intended_overlap)),
                "final_minimum_intended_am_overlap": float(
                    np.min(intended_overlap[-1])),
                "initial_myosin_separation": float(myosin_separation[0]),
                "final_myosin_separation": float(myosin_separation[-1]),
                "maximum_myosin_separation_change": float(np.max(
                    np.abs(myosin_separation - myosin_separation[0]))),
                "central_status2_frame_fraction": float(np.mean(
                    np.all(actin_status[:, :2] == 2, axis=1))),
                "counter_status2_frame_fraction": float(np.mean(
                    np.all(actin_status[:, 2:4] == 2, axis=1))),
                "final_strong_heads": int(strong[-1]),
                "zero_strong_frame_fraction": float(np.mean(strong == 0)),
                "final_q": float(q[-1]),
            })
        trace = dict(time=time, calcium=calcium, q=q, weak=weak, strong=strong,
                     aa1=aa1, aa2=aa2, occupancy=occupancy,
                     bilateral=bilateral, support=support, angle=angle)
    return row, trace


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--plot", type=Path, required=True)
    args = parser.parse_args()
    rows, traces = [], {}
    for path in sorted(args.directory.glob("*_rep*.h5")):
        row, trace = analyze(path)
        rows.append(row)
        traces.setdefault(row["variant"], []).append(trace)
    if not rows:
        raise SystemExit(f"no trajectories found in {args.directory}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    # Some geometries expose additional diagnostics (for example, the
    # force-balanced cell has counter-actins). Preserve the union of columns
    # instead of assuming every variant has the first row's schema.
    fieldnames = list(rows[0])
    seen = set(fieldnames)
    for row in rows[1:]:
        for field in row:
            if field not in seen:
                fieldnames.append(field)
                seen.add(field)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader(); writer.writerows(rows)

    colors = dict(perfect="#222222", balanced="#0072b2", anneal="#56b4e9", one_motor="#e69f00",
                  no_motor="#cc79a7", parallel="#009e73", wrong_polarity="#d55e00")
    fig, axes = plt.subplots(5, 1, figsize=(11, 13), sharex=True)
    for variant, members in traces.items():
        length = min(len(member["time"]) for member in members)
        time = members[0]["time"][:length] * 1e3
        mean = lambda key: np.mean([member[key][:length] for member in members], axis=0)
        color = colors.get(variant)
        if not axes[0].lines:
            axes[0].plot(time, mean("calcium"), color="#777777", label="calcium")
            axes[0].plot(time, mean("q"), color="#000000", label="mean q")
        axes[1].plot(time, mean("weak"), color=color, ls="--", label=f"{variant} weak")
        axes[1].plot(time, mean("strong"), color=color, label=f"{variant} strong")
        axes[2].plot(time, mean("aa1"), color=color, ls="--", label=f"{variant} AA1")
        axes[2].plot(time, mean("aa2"), color=color, label=f"{variant} AA2")
        axes[3].plot(time, mean("bilateral"), color=color, label=variant)
        axes[4].plot(time, mean("angle"), color=color, label=variant)
    axes[0].set_ylabel("activation"); axes[0].legend(ncol=2)
    axes[1].set_ylabel("AM heads"); axes[1].legend(fontsize=7, ncol=3)
    axes[2].set_ylabel("ACTN2 links"); axes[2].legend(fontsize=7, ncol=3)
    axes[3].set_ylabel("bilateral arm load (pN)"); axes[3].legend(ncol=3, fontsize=8)
    axes[4].set_ylabel("AA misalignment (deg)"); axes[4].set_xlabel("time (ms)")
    axes[4].legend(ncol=3, fontsize=8)
    for axis in axes: axis.grid(alpha=.2)
    fig.tight_layout(); fig.savefig(args.plot, dpi=180); plt.close(fig)

    # Machine-checkable acceptance checks. The finite/safety failures are hard
    # errors; biological contrasts remain reported for interpretation.
    if not all(row["finite"] for row in rows):
        raise SystemExit("non-finite trajectory values detected")
    if any(row["force_cap_fraction"] >= 1e-4 or
           row["torque_cap_fraction"] >= 1e-4 or
           row["rotation_cap_fraction"] >= 1e-4 for row in rows):
        raise SystemExit("force/torque/rotation cap fraction exceeded 1e-4")
    print(f"wrote {args.output}")
    print(f"wrote {args.plot}")


if __name__ == "__main__":
    main()
