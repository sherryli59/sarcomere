#!/usr/bin/env python3
"""Read-only structural and numerical health check for the explicit-titin pilot."""

from __future__ import annotations

import argparse
import json
import math

import h5py
import numpy as np


def dataset_sum(handle: h5py.File, path: str, frame: int) -> float:
    if path not in handle:
        return math.nan
    return float(np.asarray(handle[path][frame], dtype=float).sum())


def cap_fraction(
    handle: h5py.File, path: str, frame: int, steps: np.ndarray, count: int
) -> float:
    if path not in handle or frame == 0:
        return math.nan
    delta_steps = int(steps[frame] - steps[frame - 1])
    if delta_steps <= 0 or count <= 0:
        return math.nan
    current = dataset_sum(handle, path, frame)
    previous = dataset_sum(handle, path, frame - 1)
    return (current - previous) / (count * delta_steps)


def cumulative_fraction_series(
    handle: h5py.File, path: str, frame: int, steps: np.ndarray, count: int
) -> np.ndarray:
    if path not in handle or frame == 0 or count <= 0:
        return np.empty(0, dtype=float)
    values = np.asarray(handle[path][: frame + 1], dtype=float)
    totals = values.reshape(values.shape[0], -1).sum(axis=1)
    delta_steps = np.diff(steps[: frame + 1]).astype(float)
    deltas = np.diff(totals)
    valid = delta_steps > 0.0
    result = np.full_like(delta_steps, np.nan, dtype=float)
    result[valid] = deltas[valid] / (count * delta_steps[valid])
    return result


def actn2_support(
    handle: h5py.File, frame: int, n_actin: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Boolean masks over actins: any live ACTN2 link, status 1, status 2.

    Uses exactly the predicate the C++ capture gate uses -- membership in at
    least one live A-A link, of any status -- so the diagnostic and the gate
    cannot disagree about what "supported" means.
    """
    any_link = np.zeros(n_actin, dtype=bool)
    status1 = np.zeros(n_actin, dtype=bool)
    status2 = np.zeros(n_actin, dtype=bool)
    if "/catch_bond/local_link_frame_offsets" not in handle:
        return any_link, status1, status2
    if "/catch_bond/local_link_identity" not in handle:
        return any_link, status1, status2
    offsets = np.asarray(
        handle["/catch_bond/local_link_frame_offsets"], dtype=np.int64
    )
    if frame >= offsets.shape[0]:
        return any_link, status1, status2
    _, start, count = offsets[frame]
    if count <= 0:
        return any_link, status1, status2
    identity = np.asarray(
        handle["/catch_bond/local_link_identity"][start : start + count],
        dtype=np.int64,
    )
    # Columns: id_hi, id_lo, actin_i, actin_j, status, attach_step, s2_memory.
    partners = identity[:, 2:4]
    statuses = identity[:, 4]
    for column in range(2):
        ids = partners[:, column]
        valid = (ids >= 0) & (ids < n_actin)
        any_link[ids[valid]] = True
        status1[ids[valid & (statuses == 1)]] = True
        status2[ids[valid & (statuses == 2)]] = True
    return any_link, status1, status2


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("trajectory")
    parser.add_argument(
        "--saturation-occupancy",
        type=float,
        default=0.98,
        help="titin-bound fraction at or above which occupancy counts as "
        "saturated for the immediate-saturation check",
    )
    parser.add_argument(
        "--saturation-support",
        type=float,
        default=0.50,
        help="ACTN2-linked actin fraction below which saturated occupancy is "
        "judged nonspecific rather than structurally earned",
    )
    args = parser.parse_args()

    failures: list[str] = []
    warnings: list[str] = []
    try:
        with h5py.File(args.trajectory, "r", locking=False) as handle:
            steps = np.asarray(handle["/state/current_step"][:, 0], dtype=np.int64)
            frame = len(steps) - 1
            if frame < 0:
                raise RuntimeError("trajectory has no committed frames")

            actin_center = np.asarray(handle["/actin/center"][frame], dtype=float)
            myosin_center = np.asarray(handle["/myosin/center"][frame], dtype=float)
            actin_direction = np.asarray(
                handle["/actin/direction"][frame], dtype=float
            )
            myosin_direction = np.asarray(
                handle["/myosin/direction"][frame], dtype=float
            )
            finite = bool(
                np.isfinite(actin_center).all()
                and np.isfinite(myosin_center).all()
                and np.isfinite(actin_direction).all()
                and np.isfinite(myosin_direction).all()
            )
            if not finite:
                failures.append("non-finite filament coordinates or directions")

            titin_mode = int(handle["/meta"].attrs.get("titin_mode", -1))
            if titin_mode != 2:
                failures.append(f"titin_mode={titin_mode}, expected explicit_dynamic (2)")
            titin_schema = int(
                handle["/meta"].attrs.get("titin_schema_version", 0)
            )
            gate_enabled = bool(
                int(handle["/meta"].attrs.get("titin_require_live_actn2", 0))
            )
            required_titin_state = (
                "/state/titin_bound_actin_current",
                "/state/titin_distance_current",
                "/state/titin_extension_current",
                "/state/titin_force_current",
            )
            missing = [path for path in required_titin_state if path not in handle]
            if missing:
                failures.append("missing explicit-titin state: " + ", ".join(missing))

            bound_ids = np.empty(0, dtype=np.int64)
            if not missing:
                owners = np.asarray(
                    handle["/state/titin_bound_actin_current"][frame],
                    dtype=np.int64,
                )
                bound_ids = owners[owners >= 0]
                if len(bound_ids) != len(np.unique(bound_ids)):
                    failures.append("an actin is owned by more than one titin molecule")

            local_links = 0
            max_link_length = math.nan
            if "/catch_bond/local_link_frame_offsets" in handle:
                offsets = np.asarray(
                    handle["/catch_bond/local_link_frame_offsets"], dtype=np.int64
                )
                _, start, local_links = offsets[frame]
                if local_links and "/catch_bond/local_link_mechanics" in handle:
                    mechanics = np.asarray(
                        handle["/catch_bond/local_link_mechanics"][
                            start : start + local_links
                        ],
                        dtype=float,
                    )
                    max_link_length = float(np.max(mechanics[:, 3]))
                    contour = float(
                        handle["/meta"].attrs.get("aa_link_rupture_length", 0.0)
                    )
                    if contour > 0.0 and max_link_length > contour + 1.0e-9:
                        failures.append(
                            f"ACTN2 length {max_link_length:.6g} exceeds "
                            f"contour {contour:.6g} um"
                        )

            n_actin = int(actin_center.shape[0])
            n_myosin = int(myosin_center.shape[0])

            # ------------------------------------------------------------------
            # Structural occupancy gate.
            #
            # The v4 pilot was numerically healthy and scientifically dead: all
            # 1,280 actins were titin-bound by the first analyzed checkpoint
            # while only 115 ACTN2 links existed. A monitor that reports only
            # "bound_titins" scores that as success, so occupancy is judged here
            # against the ACTN2 structure that is supposed to license it, never
            # against an arbitrary target percentage.
            # ------------------------------------------------------------------
            any_link, status1_mask, status2_mask = actn2_support(
                handle, frame, n_actin
            )
            actn2_linked_fraction = float(any_link.mean()) if n_actin else math.nan
            actn2_status1_fraction = (
                float(status1_mask.mean()) if n_actin else math.nan
            )
            actn2_status2_fraction = (
                float(status2_mask.mean()) if n_actin else math.nan
            )
            titin_bound_fraction = (
                float(len(np.unique(bound_ids))) / n_actin
                if n_actin and not missing else math.nan
            )
            supported_bound_fraction = math.nan
            if not missing and len(bound_ids):
                owned = np.unique(bound_ids)
                owned = owned[(owned >= 0) & (owned < n_actin)]
                if owned.size:
                    supported_bound_fraction = float(any_link[owned].mean())

            supported_assoc = unsupported_assoc = rebind_assoc = -1
            if "/state/titin_association_counts" in handle:
                counts = np.asarray(
                    handle["/state/titin_association_counts"][frame],
                    dtype=np.int64,
                ).ravel()
                if counts.size >= 3:
                    supported_assoc = int(counts[0])
                    unsupported_assoc = int(counts[1])
                    rebind_assoc = int(counts[2])
            elif titin_mode == 2:
                warnings.append(
                    "trajectory predates titin schema 3: association "
                    "provenance is unavailable, so the unsupported-association "
                    "gate cannot be evaluated"
                )

            # An ordinary association that formed without live ACTN2 support is
            # gate leakage, not a tunable. Same-partner recaptures inside the
            # rebinding window are exempt by design and counted separately.
            if gate_enabled and unsupported_assoc > 0:
                failures.append(
                    f"{unsupported_assoc} titin associations formed without "
                    "live ACTN2 support while the capture gate was enabled"
                )

            # Immediate saturation: occupancy pinned at the monovalent cap while
            # ACTN2 structure is still sparse. This is the v4 failure signature
            # and must be reported as a scientific failure, not a healthy run.
            if (
                math.isfinite(titin_bound_fraction)
                and math.isfinite(actn2_linked_fraction)
                and titin_bound_fraction >= args.saturation_occupancy
                and actn2_linked_fraction < args.saturation_support
            ):
                failures.append(
                    f"titin occupancy saturated ({titin_bound_fraction:.3f}) "
                    f"while ACTN2 support remains sparse "
                    f"({actn2_linked_fraction:.3f}): nonspecific capture"
                )

            # Capture must be able to nucleate at all. Status 1 is the soft
            # encounter complex that licenses titin; if it is absent the gate
            # has replaced saturation with starvation.
            if (
                math.isfinite(actn2_linked_fraction)
                and actn2_linked_fraction <= 0.0
                and frame > 0
            ):
                failures.append(
                    "no actin carries a live ACTN2 link: titin capture cannot "
                    "nucleate"
                )
            elif (
                gate_enabled
                and math.isfinite(titin_bound_fraction)
                and titin_bound_fraction <= 0.0
                and math.isfinite(actn2_linked_fraction)
                and actn2_linked_fraction > 0.0
                and frame > 0
            ):
                warnings.append(
                    "ACTN2 links exist but no titin is bound: capture may be "
                    "starved rather than specific"
                )

            # Stochastic titin turnover must remain visible; a frozen bond set
            # is not the annealing mechanism the model claims.
            titin_release_events = 0
            if "/titin/events" in handle:
                events = np.asarray(handle["/titin/events"], dtype=float)
                if events.ndim == 2 and events.shape[1] >= 2:
                    titin_release_events = int((events[:, 1] == 2).sum())

            actin_force_cap_fraction = cap_fraction(
                handle, "/actin/cumulative_force_cap_events",
                frame, steps, n_actin
            )
            myosin_force_cap_fraction = cap_fraction(
                handle, "/myosin/cumulative_force_cap_events",
                frame, steps, n_myosin
            )
            actin_torque_cap_fraction = cap_fraction(
                handle, "/actin/cumulative_torque_cap_events",
                frame, steps, n_actin
            )
            myosin_torque_cap_fraction = cap_fraction(
                handle, "/myosin/cumulative_torque_cap_events",
                frame, steps, n_myosin
            )
            force_cap_fraction = (
                n_actin * actin_force_cap_fraction
                + n_myosin * myosin_force_cap_fraction
            ) / (n_actin + n_myosin)
            torque_cap_fraction = (
                n_actin * actin_torque_cap_fraction
                + n_myosin * myosin_torque_cap_fraction
            ) / (n_actin + n_myosin)
            actin_force_series = cumulative_fraction_series(
                handle, "/actin/cumulative_force_cap_events",
                frame, steps, n_actin
            )
            myosin_force_series = cumulative_fraction_series(
                handle, "/myosin/cumulative_force_cap_events",
                frame, steps, n_myosin
            )
            actin_torque_series = cumulative_fraction_series(
                handle, "/actin/cumulative_torque_cap_events",
                frame, steps, n_actin
            )
            myosin_torque_series = cumulative_fraction_series(
                handle, "/myosin/cumulative_torque_cap_events",
                frame, steps, n_myosin
            )
            force_series = (
                n_actin * actin_force_series + n_myosin * myosin_force_series
            ) / (n_actin + n_myosin)
            torque_series = (
                n_actin * actin_torque_series + n_myosin * myosin_torque_series
            ) / (n_actin + n_myosin)
            max_force_cap_fraction = (
                float(np.nanmax(force_series)) if force_series.size else math.nan
            )
            max_torque_cap_fraction = (
                float(np.nanmax(torque_series)) if torque_series.size else math.nan
            )
            for name, value in (
                ("maximum force-cap fraction", max_force_cap_fraction),
                ("maximum torque-cap fraction", max_torque_cap_fraction),
            ):
                if math.isfinite(value) and value >= 1.0e-3:
                    failures.append(
                        f"{name} is clamp-dominated ({value:.3g})"
                    )
                elif math.isfinite(value) and value >= 1.0e-4:
                    warnings.append(f"{name} is elevated ({value:.3g})")

            hard_contact_event_series = cumulative_fraction_series(
                handle, "/myosin/cumulative_hard_contact_events",
                frame, steps, n_myosin
            )
            hard_contact_correction_series = cumulative_fraction_series(
                handle, "/myosin/cumulative_hard_contact_correction",
                frame, steps, n_myosin
            )

            result = {
                "healthy": not failures,
                "warnings": warnings,
                "failures": failures,
                "frames": int(len(steps)),
                "step": int(steps[frame]),
                "simulated_ms": float(steps[frame]) * float(
                    handle["/meta"].attrs.get("dt", 1.0e-6)
                ) * 1.0e3,
                "n_actins": n_actin,
                "n_myosins": n_myosin,
                "finite": finite,
                "explicit_titin_molecules": int(
                    handle["/state/titin_bound_actin_current"].shape[1]
                ) if not missing else 0,
                "bound_titins": int(len(bound_ids)),
                "unique_titin_actin_owners": int(len(np.unique(bound_ids))),
                "titin_schema_version": titin_schema,
                "titin_require_live_actn2": gate_enabled,
                "titin_bound_fraction": titin_bound_fraction,
                "actn2_linked_actin_fraction": actn2_linked_fraction,
                "actn2_status1_actin_fraction": actn2_status1_fraction,
                "actn2_status2_actin_fraction": actn2_status2_fraction,
                "titin_bound_actins_with_actn2_support_fraction":
                    supported_bound_fraction,
                "titin_associations_actn2_supported": supported_assoc,
                "titin_associations_unsupported": unsupported_assoc,
                "titin_associations_same_partner_rebind": rebind_assoc,
                "titin_release_events": titin_release_events,
                "local_actn2_links": int(local_links),
                "max_actn2_length_um": max_link_length,
                "weak_heads": dataset_sum(handle, "/actin/am_weak_count", frame),
                "strong_heads": dataset_sum(handle, "/actin/am_strong_count", frame),
                "mean_am_overlap_um": float(
                    np.asarray(
                        handle["/actin/am_overlap_available"][frame], dtype=float
                    ).mean()
                ),
                "force_cap_fraction_last_interval": force_cap_fraction,
                "torque_cap_fraction_last_interval": torque_cap_fraction,
                "force_cap_fraction_max_interval": max_force_cap_fraction,
                "torque_cap_fraction_max_interval": max_torque_cap_fraction,
                "actin_force_cap_fraction_last_interval":
                    actin_force_cap_fraction,
                "myosin_force_cap_fraction_last_interval":
                    myosin_force_cap_fraction,
                "hard_contact_events_per_myosin_step_last_interval": float(
                    hard_contact_event_series[-1]
                ) if hard_contact_event_series.size else math.nan,
                "hard_contact_events_per_myosin_step_max_interval": float(
                    np.nanmax(hard_contact_event_series)
                ) if hard_contact_event_series.size else math.nan,
                "hard_contact_correction_um_per_myosin_step_last_interval": float(
                    hard_contact_correction_series[-1]
                ) if hard_contact_correction_series.size else math.nan,
                "hard_contact_correction_um_per_myosin_step_max_interval": float(
                    np.nanmax(hard_contact_correction_series)
                ) if hard_contact_correction_series.size else math.nan,
                "titin_events": int(handle["/titin/events"].shape[0])
                if "/titin/events" in handle else 0,
                "actn2_contour_ruptures": int(
                    handle["/catch_bond/limit_removal"].shape[0]
                ),
            }
    except (KeyError, OSError, RuntimeError, ValueError) as error:
        result = {"healthy": False, "failures": [str(error)], "warnings": []}

    print(json.dumps(result, allow_nan=True, indent=2))
    return 0 if result["healthy"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
