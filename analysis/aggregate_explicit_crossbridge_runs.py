#!/usr/bin/env python3
"""Aggregate structural trends and mechanical/KMC audits for explicit AM runs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np

from summarize_am_transitions import build_dwells, read_events, read_last_step


def load_csv(path: Path) -> dict[str, np.ndarray]:
    data = np.genfromtxt(path, delimiter=",", names=True)
    return {name: np.atleast_1d(data[name]) for name in data.dtype.names or ()}


def write_rows(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("analysis_dir", type=Path)
    parser.add_argument("--dt", type=float, default=1e-5)
    parser.add_argument("--k-weak-off", type=float, default=200.0)
    parser.add_argument("--k-powerstroke", type=float, default=200.0)
    parser.add_argument("--k-strong-off", type=float, default=80.0)
    parser.add_argument("--force-scale", type=float, default=1.0)
    parser.add_argument(
        "--saved-only", action="store_true",
        help="Ignore mutable /resume event buffers for active trajectories",
    )
    parser.add_argument(
        "--structural-only", action="store_true",
        help="Skip memory-intensive production transition-event aggregation",
    )
    args = parser.parse_args()

    trajectories = sorted(args.input_dir.glob("A2000_M1000_explicit_cb_*.h5"))
    if not trajectories:
        parser.error("no production trajectories found")
    series = []
    for path in trajectories:
        csv_path = args.analysis_dir / f"{path.stem}_timeseries.csv"
        if not csv_path.exists():
            parser.error(f"missing per-run structural analysis: {csv_path}")
        series.append(load_csv(csv_path))

    fields = (
        "phi_3d_myosin", "actin_nematic_order", "myosin_nematic_order",
        "n_status1_bonds", "n_status2_bonds", "n_aa_bonds",
    )
    n = min(len(run["time"]) for run in series)
    aggregate_rows: list[dict] = []
    for frame in range(n):
        row = {
            "frame": frame,
            "step": int(series[0]["step"][frame]),
            "time": float(series[0]["time"][frame]),
        }
        for field in fields:
            values = np.asarray([run[field][frame] for run in series], dtype=float)
            row[f"mean_{field}"] = float(np.mean(values))
            row[f"sem_{field}"] = (
                float(np.std(values, ddof=1) / np.sqrt(values.size))
                if values.size > 1 else 0.0
            )
        aggregate_rows.append(row)
    write_rows(args.analysis_dir / "explicit_crossbridge_aggregate_timeseries.csv",
               aggregate_rows)

    time = np.asarray([row["time"] for row in aggregate_rows])
    fig, axes = plt.subplots(3, 1, figsize=(8.2, 8.6), sharex=True)
    plot_specs = (
        (axes[0], (("phi_3d_myosin", "#3b528b", r"Myosin $\Phi_{3D}$"),)),
        (axes[1], (("actin_nematic_order", "#21918c", "Actin"),
                   ("myosin_nematic_order", "#7e57c2", "Myosin"))),
        (axes[2], (("n_status1_bonds", "#e69f00", "Status 1"),
                   ("n_status2_bonds", "#0072b2", "Status 2"))),
    )
    for axis, specs in plot_specs:
        for field, color, label in specs:
            mean = np.asarray([row[f"mean_{field}"] for row in aggregate_rows])
            sem = np.asarray([row[f"sem_{field}"] for row in aggregate_rows])
            axis.plot(time, mean, color=color, lw=2, label=label)
            axis.fill_between(time, mean - sem, mean + sem, color=color, alpha=0.2)
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(alpha=0.2)
        axis.legend(frameon=False, ncol=2)
    axes[0].set_ylabel(r"Myosin $\Phi_{3D}$")
    axes[0].set_ylim(0.0, 1.0)
    axes[1].set_ylabel("Nematic order")
    axes[1].set_ylim(0.0, 1.0)
    axes[2].set_ylabel("ACTN2 bonds")
    axes[2].set_xlabel("Simulation time")
    fig.tight_layout()
    fig.savefig(args.analysis_dir / "explicit_crossbridge_aggregate_trends.png",
                dpi=220, bbox_inches="tight")
    plt.close(fig)

    if args.structural_only:
        print(
            f"runs={len(trajectories)} frames={n} structural_only=1 "
            f"analysis={args.analysis_dir}"
        )
        return

    kmc_rows: list[dict] = []
    load_rows: list[dict] = []
    for path in trajectories:
        events = read_events(path, include_pending=not args.saved_only)
        dwells = build_dwells(
            events, args.dt,
            read_last_step(path, include_resume=not args.saved_only),
        )
        weak_exposure = sum(float(row["dwell_time"]) for row in dwells
                            if row["state"] == "weak")
        weak_kinetic = [row for row in dwells if row["state"] == "weak" and
                        row["exit_reason"] in ("weak_detachment", "powerstroke")]
        power = sum(row["exit_reason"] == "powerstroke" for row in weak_kinetic)
        weak_rate = len(weak_kinetic) / weak_exposure if weak_exposure else np.nan
        branch = power / len(weak_kinetic) if weak_kinetic else np.nan
        strong_exposure = sum(float(row["dwell_time"]) for row in dwells
                              if row["state"] == "strong")
        strong_detach = [row for row in dwells if row["state"] == "strong" and
                         row["exit_reason"] == "strong_detachment"]
        strong_rate = len(strong_detach) / strong_exposure if strong_exposure else np.nan
        force_samples = np.asarray([
            np.nanmean([row["entry_force"], row["exit_force"]])
            for row in strong_detach
        ], dtype=float)
        expected_strong = (
            float(np.mean(args.k_strong_off * np.exp(
                -np.maximum(force_samples, 0.0) / args.force_scale)))
            if force_samples.size else np.nan
        )
        kmc_rows.append({
            "run": path.stem,
            "weak_kinetic_exits": len(weak_kinetic),
            "weak_exposure": weak_exposure,
            "weak_exit_rate": weak_rate,
            "weak_expected_rate": args.k_weak_off + args.k_powerstroke,
            "weak_rate_ratio": weak_rate / (args.k_weak_off + args.k_powerstroke),
            "powerstroke_branch": branch,
            "powerstroke_expected_branch": args.k_powerstroke /
                (args.k_weak_off + args.k_powerstroke),
            "strong_detachments": len(strong_detach),
            "strong_exposure": strong_exposure,
            "strong_detachment_rate": strong_rate,
            "strong_endpoint_expected_rate": expected_strong,
            "strong_rate_ratio_approx": strong_rate / expected_strong,
        })

        with h5py.File(path, "r", locking=False) as handle:
            loads = np.asarray(handle["/catch_bond/load_transitions"], dtype=float)
            schema = int(handle["/meta"].attrs.get("aa_load_transition_schema", 1))
        if schema != 2:
            raise RuntimeError(f"{path}: expected AA load schema 2, got {schema}")
        promotions = loads[(loads[:, 3] == 1) & (loads[:, 4] == 2)]
        demotions = loads[(loads[:, 3] == 2) & (loads[:, 4] == 1)]
        load_rows.append({
            "run": path.stem,
            "n_load_events": len(loads),
            "n_promotions": len(promotions),
            "n_demotions": len(demotions),
            "mean_motor_only_at_promotion": float(np.mean(promotions[:, 5]))
                if len(promotions) else np.nan,
            "mean_opposed_drive_at_promotion": float(np.mean(promotions[:, 6]))
                if len(promotions) else np.nan,
            "mean_transmitted_load_at_promotion": float(np.mean(promotions[:, 7]))
                if len(promotions) else np.nan,
            "mean_combined_load_at_promotion": float(np.mean(promotions[:, 8]))
                if len(promotions) else np.nan,
            "min_selected_load_at_promotion": float(np.min(promotions[:, 9]))
                if len(promotions) else np.nan,
            "max_selected_load_at_demotion": float(np.max(demotions[:, 9]))
                if len(demotions) else np.nan,
        })
    write_rows(args.analysis_dir / "explicit_crossbridge_production_kmc_audit.csv",
               kmc_rows)
    write_rows(args.analysis_dir / "explicit_crossbridge_actn2_load_audit.csv",
               load_rows)
    print(f"runs={len(trajectories)} frames={n} analysis={args.analysis_dir}")


if __name__ == "__main__":
    main()
