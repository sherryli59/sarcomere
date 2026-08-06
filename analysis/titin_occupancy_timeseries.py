#!/usr/bin/env python3
"""Per-frame titin occupancy vs ACTN2 support, and the structural gate.

The v4 pilot failed because occupancy was scored on its own: 1,280/1,280 actins
were titin-bound while only 115 ACTN2 links existed, and every numerical health
check passed. The question this script answers is not "how many titins are
bound" but "did occupancy follow the ACTN2 structure that is supposed to license
it". It is read-only.

The gate is structural, not a target percentage:

  A. no ordinary titin association ever formed without live ACTN2 support;
  B. occupancy did not jump to the monovalent cap while ACTN2 support was still
     sparse;
  C. stochastic titin turnover remained visible;
  D. status-1 ACTN2 was common enough to nucleate capture;
  E. geometry stayed finite, caps stayed below 1e-4, and no ACTN2 contour
     rupture occurred.
"""

from __future__ import annotations

import argparse
import json
import math

import h5py
import numpy as np


def actn2_masks(
    handle: h5py.File, frame: int, n_actin: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Per-actin booleans: any live ACTN2 link, status 1, status 2; plus count.

    Mirrors the C++ capture-gate predicate -- membership in at least one live
    A-A link of any status -- so the diagnostic cannot disagree with the gate.
    """
    any_link = np.zeros(n_actin, dtype=bool)
    status1 = np.zeros(n_actin, dtype=bool)
    status2 = np.zeros(n_actin, dtype=bool)
    if (
        "/catch_bond/local_link_frame_offsets" not in handle
        or "/catch_bond/local_link_identity" not in handle
    ):
        return any_link, status1, status2, 0
    offsets = np.asarray(
        handle["/catch_bond/local_link_frame_offsets"], dtype=np.int64
    )
    if frame >= offsets.shape[0]:
        return any_link, status1, status2, 0
    _, start, count = offsets[frame]
    if count <= 0:
        return any_link, status1, status2, 0
    identity = np.asarray(
        handle["/catch_bond/local_link_identity"][start : start + count],
        dtype=np.int64,
    )
    statuses = identity[:, 4]
    for column in (2, 3):
        ids = identity[:, column]
        valid = (ids >= 0) & (ids < n_actin)
        any_link[ids[valid]] = True
        status1[ids[valid & (statuses == 1)]] = True
        status2[ids[valid & (statuses == 2)]] = True
    return any_link, status1, status2, int(count)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("trajectory")
    parser.add_argument("--csv", help="optional path for the per-frame table")
    parser.add_argument("--saturation-occupancy", type=float, default=0.98)
    parser.add_argument("--saturation-support", type=float, default=0.50)
    args = parser.parse_args()

    failures: list[str] = []
    warnings: list[str] = []
    rows: list[dict] = []

    with h5py.File(args.trajectory, "r", locking=False) as handle:
        meta = handle["/meta"].attrs
        dt = float(meta.get("dt", 5.0e-7))
        gate_enabled = bool(int(meta.get("titin_require_live_actn2", 0)))
        schema = int(meta.get("titin_schema_version", 0))
        steps = np.asarray(handle["/state/current_step"][:, 0], dtype=np.int64)
        n_frames = len(steps)
        if n_frames == 0:
            raise RuntimeError("trajectory has no committed frames")
        n_actin = int(handle["/actin/center"].shape[1])

        has_counts = "/state/titin_association_counts" in handle
        has_owner = "/state/titin_bound_actin_current" in handle
        if not has_owner:
            raise RuntimeError("trajectory has no explicit-titin ownership state")
        if gate_enabled and not has_counts:
            failures.append(
                "gate is enabled but the trajectory has no association "
                "provenance; rebuild against titin schema 3"
            )

        # Titin release events are cumulative over the whole run, so turnover is
        # attributed to frames by binning event steps.
        release_steps = np.empty(0, dtype=np.int64)
        if "/titin/events" in handle:
            events = np.asarray(handle["/titin/events"], dtype=float)
            if events.ndim == 2 and events.shape[1] >= 2:
                release_steps = events[events[:, 1] == 2, 0].astype(np.int64)

        ruptures = int(handle["/catch_bond/limit_removal"].shape[0]) \
            if "/catch_bond/limit_removal" in handle else 0

        for frame in range(n_frames):
            owners = np.asarray(
                handle["/state/titin_bound_actin_current"][frame], dtype=np.int64
            )
            bound = np.unique(owners[owners >= 0])
            bound = bound[(bound >= 0) & (bound < n_actin)]
            any_link, s1, s2, link_count = actn2_masks(handle, frame, n_actin)

            supported_assoc = unsupported_assoc = rebind_assoc = -1
            if has_counts:
                counts = np.asarray(
                    handle["/state/titin_association_counts"][frame],
                    dtype=np.int64,
                ).ravel()
                if counts.size >= 3:
                    supported_assoc = int(counts[0])
                    unsupported_assoc = int(counts[1])
                    rebind_assoc = int(counts[2])

            lower = steps[frame - 1] if frame else -1
            releases = int(
                ((release_steps > lower) & (release_steps <= steps[frame])).sum()
            )

            rows.append(
                {
                    "frame": frame,
                    "step": int(steps[frame]),
                    "ms": float(steps[frame]) * dt * 1.0e3,
                    "titin_bound_actins": int(bound.size),
                    "titin_bound_fraction": bound.size / n_actin,
                    "actn2_links": link_count,
                    "actn2_linked_fraction": float(any_link.mean()),
                    "actn2_status1_fraction": float(s1.mean()),
                    "actn2_status2_fraction": float(s2.mean()),
                    "bound_with_support_fraction": (
                        float(any_link[bound].mean()) if bound.size else math.nan
                    ),
                    "assoc_supported": supported_assoc,
                    "assoc_unsupported": unsupported_assoc,
                    "assoc_rebind": rebind_assoc,
                    "titin_releases_in_interval": releases,
                }
            )

        final = rows[-1]

        # A. Gate leakage. A single unsupported ordinary association is a code
        #    failure, not a tunable.
        if gate_enabled and final["assoc_unsupported"] > 0:
            failures.append(
                f"{final['assoc_unsupported']} titin associations formed "
                "without live ACTN2 support while the gate was enabled"
            )

        # B. Immediate saturation, checked on every frame rather than only the
        #    last one: the v4 failure was already complete by the first
        #    analyzed checkpoint and a last-frame test could miss a transient.
        saturated = [
            row
            for row in rows
            if row["titin_bound_fraction"] >= args.saturation_occupancy
            and row["actn2_linked_fraction"] < args.saturation_support
        ]
        if saturated:
            first = saturated[0]
            failures.append(
                f"titin occupancy saturated at {first['ms']:.2f} ms "
                f"({first['titin_bound_fraction']:.3f}) while ACTN2 support was "
                f"sparse ({first['actn2_linked_fraction']:.3f})"
            )

        # B2. Maturation. Status 1 is a soft encounter complex; only status 2
        #     carries the axial register that produces the target A-A overlap.
        #     A run can therefore hold a large, perfectly stable ACTN2
        #     population and still contain no Z-disc at all -- two 1,280-actin
        #     production runs passed every other criterion here with 244 and 282
        #     links of which 0 and 2 were mature. Maturation is slow (status 2
        #     first appears near 25-35 ms and plateaus near 90 ms at the reduced
        #     motor drive), so this is only scored once the run is long enough
        #     to have had the chance.
        final_s2 = final["actn2_status2_fraction"]
        final_links = final["actn2_links"]
        if final["ms"] >= 40.0 and final_links > 0:
            if final_s2 <= 0.0:
                failures.append(
                    f"no status-2 ACTN2 link at {final['ms']:.0f} ms despite "
                    f"{final_links} live links: the network never matured and "
                    "carries no axial register"
                )
            elif final_s2 < 0.2 * final["actn2_linked_fraction"]:
                warnings.append(
                    f"only {final_s2:.3f} of actins hold a mature status-2 link "
                    f"against {final['actn2_linked_fraction']:.3f} linked "
                    "overall; the population is dominated by unregistered "
                    "encounter complexes"
                )
        elif final_links > 0 and final_s2 <= 0.0:
            warnings.append(
                f"no status-2 ACTN2 link yet, but the run is only "
                f"{final['ms']:.0f} ms; maturation is not expected to be "
                "resolvable before ~40 ms"
            )

        # C. Turnover must remain visible; a frozen bond set is not annealing.
        #
        # This has to be judged against the statistical power of the run, not
        # against zero. With k_off = 0.5 /s, a 1.5 ms trajectory holding 80
        # bound titins expects 0.06 releases, so "no releases" there is correct
        # physics rather than a frozen bond set. Integrate the bound population
        # over time to get the expected count and only fail when the run was
        # actually powered to see turnover.
        total_releases = sum(row["titin_releases_in_interval"] for row in rows)
        k_off = float(meta.get("titin_k_off", 0.0))
        bound_seconds = 0.0
        for index in range(1, n_frames):
            interval = float(steps[index] - steps[index - 1]) * dt
            mean_bound = 0.5 * (
                rows[index]["titin_bound_actins"]
                + rows[index - 1]["titin_bound_actins"]
            )
            bound_seconds += mean_bound * interval
        expected_releases = k_off * bound_seconds
        if total_releases == 0 and expected_releases >= 3.0:
            failures.append(
                f"no titin dissociation events although {expected_releases:.1f} "
                "were expected: turnover is not occurring"
            )
        elif total_releases == 0 and n_frames > 1:
            warnings.append(
                f"no titin dissociation events, but only "
                f"{expected_releases:.2f} were expected at k_off="
                f"{k_off:g} /s; this run is underpowered to test turnover"
            )

        # D. Capture must be able to nucleate at all. This is the failure mode
        #    the gate itself could introduce, and it is the reason status 2 is
        #    not required.
        if final["actn2_linked_fraction"] <= 0.0:
            failures.append(
                "no actin carries a live ACTN2 link: capture cannot nucleate"
            )
        elif final["actn2_status1_fraction"] <= 0.0:
            warnings.append(
                "no status-1 ACTN2 at the final frame: nucleation rests "
                "entirely on mature status-2 links"
            )
        if gate_enabled and final["titin_bound_fraction"] <= 0.0:
            failures.append(
                "the gate is enabled and no titin ever bound: capture is "
                "starved rather than specific"
            )

        # Occupancy should track support rather than lead it by a wide margin.
        if (
            final["titin_bound_fraction"] > 0.0
            and math.isfinite(final["bound_with_support_fraction"])
            and final["bound_with_support_fraction"] < 0.5
        ):
            warnings.append(
                "fewer than half of titin-bound actins currently hold an ACTN2 "
                f"link ({final['bound_with_support_fraction']:.3f}); expected "
                "under turnover, but check it is not drifting toward zero"
            )

        # E. Numerical health, so the structural verdict is not read off an
        #    unstable trajectory.
        for group, count in (("actin", n_actin),
                             ("myosin", int(handle["/myosin/center"].shape[1]))):
            path = f"/{group}/cumulative_force_cap_events"
            if path in handle and n_frames > 1:
                totals = np.asarray(handle[path], dtype=float)
                totals = totals.reshape(totals.shape[0], -1).sum(axis=1)
                dsteps = np.diff(steps).astype(float)
                frac = np.divide(
                    np.diff(totals), count * dsteps,
                    out=np.full(dsteps.shape, np.nan), where=dsteps > 0,
                )
                worst = float(np.nanmax(frac)) if frac.size else math.nan
                if math.isfinite(worst) and worst >= 1.0e-3:
                    failures.append(
                        f"{group} force-cap fraction is clamp-dominated "
                        f"({worst:.3g})"
                    )
                elif math.isfinite(worst) and worst >= 1.0e-4:
                    warnings.append(
                        f"{group} force-cap fraction elevated ({worst:.3g})"
                    )
        if ruptures > 0:
            failures.append(f"{ruptures} ACTN2 contour ruptures")

        finite = all(
            bool(np.isfinite(np.asarray(handle[path][-1], dtype=float)).all())
            for path in (
                "/actin/center", "/myosin/center",
                "/actin/direction", "/myosin/direction",
            )
        )
        if not finite:
            failures.append("non-finite filament coordinates or directions")

    if args.csv:
        import csv

        with open(args.csv, "w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    summary = {
        "trajectory": args.trajectory,
        "passed": not failures,
        "failures": failures,
        "warnings": warnings,
        "titin_schema_version": schema,
        "titin_require_live_actn2": gate_enabled,
        "frames": len(rows),
        "final_ms": rows[-1]["ms"],
        "n_actins": n_actin,
        "titin_releases_total": total_releases,
        "titin_releases_expected": expected_releases,
        "actn2_status2_fraction_final": final["actn2_status2_fraction"],
        "actn2_mature_link_share": (
            final["actn2_status2_fraction"] / final["actn2_linked_fraction"]
            if final["actn2_linked_fraction"] > 0 else float("nan")
        ),
        "actn2_contour_ruptures": ruptures,
        "first_frame": rows[0],
        "final_frame": rows[-1],
        "occupancy_trace": [
            {
                "ms": round(row["ms"], 4),
                "titin": round(row["titin_bound_fraction"], 4),
                "actn2": round(row["actn2_linked_fraction"], 4),
                "s1": round(row["actn2_status1_fraction"], 4),
                "s2": round(row["actn2_status2_fraction"], 4),
                "unsupported": row["assoc_unsupported"],
            }
            for row in rows
        ],
    }
    print(json.dumps(summary, allow_nan=True, indent=2))
    return 0 if summary["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
