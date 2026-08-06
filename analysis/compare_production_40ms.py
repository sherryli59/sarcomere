#!/usr/bin/env python3
"""Apples-to-apples comparison of the current and model-19 production runs.

All trajectories are truncated at 40 ms.  The current run has an exact 40 ms
frame; the older runs save on odd milliseconds, so their last retained frame is
39 ms.  A-M pair angles are counted once per unique head-bound actin/myosin
pair, rather than once per head.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

import h5py
import numpy as np

os.environ.setdefault("MPLCONFIGDIR", "/tmp/sarcomere_matplotlib")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


RUNS = {
    "current model 21": Path(
        "/scratch/users/xsli/sarcomere/prod_theta01_weaker_titin_20260805_r2/data/"
        "A1280_M640_thetaPerp0p1_titinK40_L0p75_seed20260805_40ms.h5"
    ),
    "older model 19, seed 1": Path(
        "/scratch/users/xsli/sarcomere/prod_titin_serial_20260802/data/"
        "PRODALIGN120_A1280_M640_h35_100ms_seed2080201.h5"
    ),
    "older model 19, seed 2": Path(
        "/scratch/users/xsli/sarcomere/prod_titin_serial_20260802/data/"
        "PRODALIGN120_A1280_M640_h35_100ms_seed2080202.h5"
    ),
}


def nematic(directions: np.ndarray) -> tuple[float, np.ndarray]:
    unit = np.asarray(directions, dtype=float)
    unit /= np.linalg.norm(unit, axis=1, keepdims=True)
    q = 1.5 * unit.T @ unit / len(unit) - 0.5 * np.eye(3)
    eigenvalues, eigenvectors = np.linalg.eigh(q)
    return float(eigenvalues[-1]), eigenvectors[:, -1]


def acute_angle_degrees(first: np.ndarray, second: np.ndarray) -> float:
    cosine = abs(float(np.dot(first, second)))
    return float(np.degrees(np.arccos(np.clip(cosine, 0.0, 1.0))))


def bound_pair_angles(
    actin_direction: np.ndarray,
    myosin_direction: np.ndarray,
    bound_actin: np.ndarray,
    head_state: np.ndarray,
) -> np.ndarray:
    active_heads = np.argwhere(
        (bound_actin >= 0) & np.isin(head_state, (1, 2))
    )
    pairs: set[tuple[int, int]] = set()
    values: list[float] = []
    for myosin_index, head_index in active_heads:
        actin_index = int(bound_actin[myosin_index, head_index])
        pair = (actin_index, int(myosin_index))
        if pair in pairs:
            continue
        pairs.add(pair)
        values.append(
            acute_angle_degrees(
                actin_direction[actin_index], myosin_direction[myosin_index]
            )
        )
    return np.asarray(values, dtype=float)


def frame_metrics(handle: h5py.File, frame: int) -> dict[str, float]:
    actin_direction = np.asarray(handle["actin/direction"][frame], dtype=float)
    myosin_direction = np.asarray(handle["myosin/direction"][frame], dtype=float)
    actin_nematic, actin_director = nematic(actin_direction)
    myosin_nematic, myosin_director = nematic(myosin_direction)
    director_angle = acute_angle_degrees(actin_director, myosin_director)

    weak = np.asarray(handle["actin/am_weak_count"][frame, :, 0], dtype=float)
    strong = np.asarray(handle["actin/am_strong_count"][frame, :, 0], dtype=float)
    am_bound = weak + strong > 0
    cb_status = np.asarray(handle["actin/cb_status"][frame, :, 0], dtype=int)
    crosslinkable = np.asarray(
        handle["actin/crosslink_ratio"][frame, :, 0], dtype=float
    )
    binding_ratio = np.asarray(
        handle["actin/myosin_binding_ratio"][frame, :, 0], dtype=float
    )
    pair_angles = bound_pair_angles(
        actin_direction,
        myosin_direction,
        np.asarray(handle["state/am_head_bound_actin_current"][frame], dtype=int),
        np.asarray(handle["state/am_head_state_current"][frame], dtype=int),
    )
    pair_median = float(np.median(pair_angles)) if len(pair_angles) else float("nan")
    pair_p90 = float(np.percentile(pair_angles, 90)) if len(pair_angles) else float("nan")
    crosslinkable_bound = (
        float(np.mean(crosslinkable[am_bound])) if np.any(am_bound) else float("nan")
    )
    binding_ratio_bound = (
        float(np.mean(binding_ratio[am_bound])) if np.any(am_bound) else float("nan")
    )
    status1_links = 0.5 * float(
        np.asarray(handle["actin/aa_local_status1_count"][frame, :, 0]).sum()
    )
    status2_links = 0.5 * float(
        np.asarray(handle["actin/aa_local_status2_count"][frame, :, 0]).sum()
    )
    linked_actins = int(np.count_nonzero(cb_status > 0))
    total_links = status1_links + status2_links
    link_endpoints_per_linked_actin = (
        2.0 * total_links / linked_actins if linked_actins else float("nan")
    )

    return {
        "actin_nematic": actin_nematic,
        "myosin_nematic": myosin_nematic,
        "am_director_angle_deg": director_angle,
        "am_pair_angle_median_deg": pair_median,
        "am_pair_angle_p90_deg": pair_p90,
        "crosslinkable_fraction_am_bound": crosslinkable_bound,
        "max_am_overlap_fraction_am_bound": binding_ratio_bound,
        "am_bound_actin_fraction": float(np.mean(am_bound)),
        "strong_actin_fraction": float(np.mean(strong > 0)),
        "weak_heads": float(weak.sum()),
        "strong_heads": float(strong.sum()),
        "status1_actins": float(np.count_nonzero(cb_status == 1)),
        "status2_actins": float(np.count_nonzero(cb_status == 2)),
        "actn2_linked_actin_fraction": float(np.mean(cb_status > 0)),
        "status2_actin_fraction": float(np.mean(cb_status == 2)),
        "status1_links": status1_links,
        "status2_links": status2_links,
        "active_actn2_links": total_links,
        "link_endpoints_per_linked_actin": link_endpoints_per_linked_actin,
        "titin_bound_fraction": float(
            np.mean(np.asarray(handle["actin/titin_bound"][frame, :, 0]) > 0)
        ),
        "mean_q": float(np.mean(handle["actin/am_q"][frame, :, 0])),
    }


def load_run(path: Path, cutoff_ms: float) -> dict[str, object]:
    with h5py.File(path, "r", locking=False) as handle:
        dt = float(handle["meta"].attrs["dt"])
        all_steps = np.asarray(handle["state/current_step"][:, 0], dtype=np.int64)
        all_time_ms = all_steps * dt * 1.0e3
        frames = np.flatnonzero(all_time_ms <= cutoff_ms + 1.0e-9)
        series = [frame_metrics(handle, int(frame)) for frame in frames]
        rupture_events = np.asarray(
            handle["actin_myo/forced_ruptures_v2"], dtype=float
        )
        rupture_events = rupture_events[
            rupture_events[:, 2] <= all_steps[frames[-1]]
        ]
    keys = list(series[0])
    arrays = {key: np.asarray([row[key] for row in series]) for key in keys}
    endpoint = {key: float(arrays[key][-1]) for key in keys}
    endpoint["endpoint_ms"] = float(all_time_ms[frames[-1]])
    endpoint["endpoint_step"] = int(all_steps[frames[-1]])
    endpoint["am_geometric_ruptures"] = int(len(rupture_events))
    endpoint["am_geometric_ruptures_weak"] = int(
        np.count_nonzero(rupture_events[:, 4] == 1)
    )
    endpoint["am_geometric_ruptures_strong"] = int(
        np.count_nonzero(rupture_events[:, 4] == 2)
    )
    tail = all_time_ms[frames] >= 30.0
    for key in (
        "actin_nematic",
        "myosin_nematic",
        "am_director_angle_deg",
        "am_pair_angle_median_deg",
        "crosslinkable_fraction_am_bound",
        "max_am_overlap_fraction_am_bound",
        "actn2_linked_actin_fraction",
        "status2_actin_fraction",
        "strong_actin_fraction",
    ):
        endpoint[f"tail30_40_mean_{key}"] = float(np.nanmean(arrays[key][tail]))
    return {
        "time_ms": all_time_ms[frames],
        "series": arrays,
        "endpoint": endpoint,
    }


def write_csv(path: Path, results: dict[str, dict[str, object]]) -> None:
    rows = []
    for label, result in results.items():
        rows.append({"run": label, **result["endpoint"]})
    fields = list(rows[0])
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def make_plot(path: Path, results: dict[str, dict[str, object]]) -> None:
    panels = (
        ("actin_nematic", "Actin nematic", "S"),
        ("myosin_nematic", "Myosin nematic", "S"),
        ("am_director_angle_deg", "A-M director angle", "degrees"),
        ("am_pair_angle_median_deg", "Median bound-pair A-M angle", "degrees"),
        (
            "crosslinkable_fraction_am_bound",
            "Free crosslinkable actin fraction",
            "fraction of actin length",
        ),
        (
            "max_am_overlap_fraction_am_bound",
            "Maximum A-M overlap fraction",
            "fraction of actin length",
        ),
        ("actn2_linked_actin_fraction", "ACTN2-linked actins", "fraction"),
        ("status2_actin_fraction", "Status-2 actins", "fraction"),
    )
    colors = ("#3366cc", "#e67e22", "#2ca25f")
    fig, axes = plt.subplots(4, 2, figsize=(11, 11), sharex=True)
    for axis, (field, title, ylabel) in zip(axes.flat, panels):
        for color, (label, result) in zip(colors, results.items()):
            axis.plot(
                result["time_ms"], result["series"][field],
                color=color, lw=2, label=label,
            )
        axis.set_title(title, loc="left")
        axis.set_ylabel(ylabel)
        axis.grid(alpha=0.2)
    for axis in axes[-1]:
        axis.set_xlabel("simulation time (ms)")
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.suptitle("Production runs compared at the same 40 ms horizon", x=0.08, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("analysis/production_40ms_comparison"))
    parser.add_argument("--cutoff-ms", type=float, default=40.0)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    results = {label: load_run(path, args.cutoff_ms) for label, path in RUNS.items()}
    write_csv(args.output / "summary.csv", results)
    make_plot(args.output / "timeseries.png", results)
    serializable = {
        label: {"path": str(RUNS[label]), **result["endpoint"]}
        for label, result in results.items()
    }
    (args.output / "summary.json").write_text(json.dumps(serializable, indent=2) + "\n")
    print(json.dumps(serializable, indent=2))


if __name__ == "__main__":
    main()
