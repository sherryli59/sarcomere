#!/usr/bin/env python3
"""Run the autonomous overlap-mechanism diagnostic campaign."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, replace
from pathlib import Path
from typing import Dict, Iterable, List

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

from overlap_model import (
    Bond,
    MECHANISMS,
    Parameters,
    aa_overlap,
    am_overlaps,
    default_state,
    force_components,
    mechanism_parameters,
    simulate,
)


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)

INITIAL_CONDITIONS = {
    "reference": (0.10, 1.125),
    "aa_low": (0.04, 1.125),
    "aa_high": (0.16, 1.125),
    "am_deep": (0.10, 1.00),
    "am_shallow": (0.10, 1.25),
}
SEEDS = tuple(range(6))


def tail_slope(time: np.ndarray, value: np.ndarray, fraction: float = 0.35) -> float:
    start = int((1.0 - fraction) * len(time))
    if len(time) - start < 2:
        return float("nan")
    return float(np.polyfit(time[start:], value[start:], 1)[0])


def tail_mean(value: np.ndarray, fraction: float = 0.35) -> float:
    start = int((1.0 - fraction) * len(value))
    return float(np.mean(value[start:]))


def campaign_parameters(name: str) -> Parameters:
    return replace(
        mechanism_parameters(name),
        dt=5.0e-4,
        duration=3.0,
        save_every=5,
    )


def run_convergence() -> List[Dict[str, float | str | int]]:
    rows: List[Dict[str, float | str | int]] = []
    representative: Dict[str, Dict[str, np.ndarray]] = {}
    for name in MECHANISMS:
        p = campaign_parameters(name)
        for ic_name, (y0, d0) in INITIAL_CONDITIONS.items():
            for seed in SEEDS:
                initial = default_state(y0, d0)
                result = simulate(p, seed=seed, initial=initial)
                y = result["overlaps"][:, 0]
                x = np.mean(result["overlaps"][:, 1:], axis=1)
                d = 0.5 * (
                    result["coords"][:, 0] - result["coords"][:, 2]
                    + result["coords"][:, 3] - result["coords"][:, 1]
                )
                rows.append(
                    {
                        "mechanism": name,
                        "initial_condition": ic_name,
                        "seed": seed,
                        "aa_initial_um": y0,
                        "zm_initial_um": d0,
                        "aa_tail_um": tail_mean(y),
                        "am_tail_um": tail_mean(x),
                        "zm_tail_um": tail_mean(d),
                        "aa_tail_velocity_um_s": tail_slope(result["time"], y),
                        "am_tail_velocity_um_s": tail_slope(result["time"], x),
                        "mean_bonds_tail": tail_mean(result["bond_count"]),
                    }
                )
                if ic_name == "reference" and seed == 0:
                    representative[name] = result
    write_csv(OUT / "convergence_runs.csv", rows)
    summarize_convergence(rows)
    make_timeseries_plot(representative)
    make_force_decomposition_plot(representative)
    for name in ("baseline_turnover", "johann_adhesion", "register_plus_titin"):
        make_animation(name, representative[name], campaign_parameters(name))
    return rows


def summarize_convergence(rows: List[Dict[str, float | str | int]]) -> None:
    summaries: List[Dict[str, float | str | bool]] = []
    for name in MECHANISMS:
        subset = [r for r in rows if r["mechanism"] == name]
        ic_aa = []
        ic_am = []
        ic_v_aa = []
        ic_v_am = []
        for ic_name in INITIAL_CONDITIONS:
            arm = [r for r in subset if r["initial_condition"] == ic_name]
            ic_aa.append(np.mean([float(r["aa_tail_um"]) for r in arm]))
            ic_am.append(np.mean([float(r["am_tail_um"]) for r in arm]))
            ic_v_aa.append(np.mean([float(r["aa_tail_velocity_um_s"]) for r in arm]))
            ic_v_am.append(np.mean([float(r["am_tail_velocity_um_s"]) for r in arm]))
        aa_spread = float(max(ic_aa) - min(ic_aa))
        am_spread = float(max(ic_am) - min(ic_am))
        max_speed = float(max(max(abs(np.array(ic_v_aa))), max(abs(np.array(ic_v_am)))))
        mean_aa = float(np.mean(ic_aa))
        mean_am = float(np.mean(ic_am))
        # Strict diagnostic definition: convergence from both sides, negligible
        # tail drift, intact finite A-A overlap, and A-M overlap below the full
        # 0.665 um head-bearing half.
        stable_both = (
            aa_spread < 0.025
            and am_spread < 0.025
            and max_speed < 0.01
            and 0.06 < mean_aa < 0.14
            and 0.02 < mean_am < 0.655
        )
        summaries.append(
            {
                "mechanism": name,
                "mean_aa_um": mean_aa,
                "mean_am_um": mean_am,
                "aa_initial_condition_spread_um": aa_spread,
                "am_initial_condition_spread_um": am_spread,
                "max_abs_tail_velocity_um_s": max_speed,
                "stable_both_by_strict_test": stable_both,
            }
        )
    write_csv(OUT / "convergence_summary.csv", summaries)
    with (OUT / "convergence_summary.json").open("w") as handle:
        json.dump(summaries, handle, indent=2)


def run_passive_force_scans() -> None:
    rows = []
    y_values = np.linspace(0.01, 0.18, 70)
    d_values = np.linspace(0.96, 1.28, 70)
    for name in MECHANISMS:
        p = campaign_parameters(name)
        for y in y_values:
            state = default_state(float(y), 1.125)
            state.bonds = [Bond(0.10) for _ in range(p.initial_bonds)]
            comps = force_components(state, p, motor_on=False)
            for component, force in comps.items():
                rows.append(
                    {
                        "mechanism": name,
                        "scan": "aa",
                        "coordinate_um": y,
                        "component": component,
                        "generalized_force_pN": 0.5 * (force[0] - force[1]),
                    }
                )
        for d in d_values:
            state = default_state(0.10, float(d))
            state.bonds = [Bond(0.10) for _ in range(p.initial_bonds)]
            comps = force_components(state, p, motor_on=False)
            for component, force in comps.items():
                # Generalized force that increases the mean Z-to-M distance.
                generalized = 0.25 * (
                    force[0] - force[2] + force[3] - force[1]
                )
                rows.append(
                    {
                        "mechanism": name,
                        "scan": "zm",
                        "coordinate_um": d,
                        "component": component,
                        "generalized_force_pN": generalized,
                    }
                )
    write_csv(OUT / "passive_force_scans.csv", rows)
    plot_passive_force_scans(rows)


def run_creep() -> None:
    names = (
        "baseline_turnover",
        "permanent_actn2",
        "johann_adhesion",
        "zdisc_register",
        "titin",
        "register_plus_titin",
    )
    rows = []
    traces = {}
    for name in names:
        p = replace(campaign_parameters(name), duration=6.0)
        for seed in SEEDS:
            result = simulate(p, seed=100 + seed, initial=default_state())
            y = result["overlaps"][:, 0]
            x = np.mean(result["overlaps"][:, 1:], axis=1)
            rows.append(
                {
                    "mechanism": name,
                    "seed": seed,
                    "aa_tail_um": tail_mean(y, 0.25),
                    "am_tail_um": tail_mean(x, 0.25),
                    "aa_creep_velocity_um_s": tail_slope(result["time"], y, 0.5),
                    "am_creep_velocity_um_s": tail_slope(result["time"], x, 0.5),
                }
            )
            if seed == 0:
                traces[name] = result
    write_csv(OUT / "creep_runs.csv", rows)
    summaries = []
    for name in names:
        subset = [r for r in rows if r["mechanism"] == name]
        summaries.append(
            {
                "mechanism": name,
                "aa_tail_um": np.mean([r["aa_tail_um"] for r in subset]),
                "am_tail_um": np.mean([r["am_tail_um"] for r in subset]),
                "aa_creep_velocity_um_s": np.mean(
                    [r["aa_creep_velocity_um_s"] for r in subset]
                ),
                "am_creep_velocity_um_s": np.mean(
                    [r["am_creep_velocity_um_s"] for r in subset]
                ),
            }
        )
    write_csv(OUT / "creep_summary.csv", summaries)
    plot_creep(traces)


def write_csv(path: Path, rows: List[Dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def make_timeseries_plot(representative: Dict[str, Dict[str, np.ndarray]]) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    for name, result in representative.items():
        axes[0].plot(result["time"], result["overlaps"][:, 0], label=name)
        axes[1].plot(
            result["time"],
            np.mean(result["overlaps"][:, 1:], axis=1),
            label=name,
        )
    axes[0].axhline(0.10, color="black", ls=":", lw=1)
    axes[0].set_ylabel("A–A overlap (µm)")
    axes[1].set_ylabel("mean A–M overlap (µm)")
    axes[1].set_xlabel("time (s)")
    axes[0].legend(ncol=3, fontsize=8)
    fig.suptitle("Reference-start stochastic trajectories")
    fig.tight_layout()
    fig.savefig(OUT / "reference_timeseries.png", dpi=180)
    plt.close(fig)


def make_force_decomposition_plot(representative: Dict[str, Dict[str, np.ndarray]]) -> None:
    chosen = ("baseline_turnover", "johann_adhesion", "register_plus_titin")
    fig, axes = plt.subplots(len(chosen), 1, figsize=(10, 9), sharex=True)
    for ax, name in zip(axes, chosen):
        result = representative[name]
        for key, values in result.items():
            if not key.startswith("force_"):
                continue
            generalized = 0.5 * (values[:, 0] - values[:, 1])
            if np.max(np.abs(generalized)) < 1.0e-8:
                continue
            ax.plot(result["time"], generalized, label=key.removeprefix("force_"))
        ax.set_ylabel(f"{name}\nA–A force (pN)")
        ax.legend(ncol=4, fontsize=8)
    axes[-1].set_xlabel("time (s)")
    fig.suptitle("Force decomposition along the A–A coordinate")
    fig.tight_layout()
    fig.savefig(OUT / "force_decomposition.png", dpi=180)
    plt.close(fig)


def plot_passive_force_scans(rows: List[Dict]) -> None:
    chosen = ("baseline_turnover", "johann_adhesion", "zdisc_register", "titin")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for name in chosen:
        for scan, ax in (("aa", axes[0]), ("zm", axes[1])):
            subset = [r for r in rows if r["mechanism"] == name and r["scan"] == scan]
            coords = sorted({float(r["coordinate_um"]) for r in subset})
            total = []
            for coordinate in coords:
                total.append(
                    sum(
                        float(r["generalized_force_pN"])
                        for r in subset
                        if float(r["coordinate_um"]) == coordinate
                    )
                )
            ax.plot(coords, total, label=name)
    axes[0].axhline(0, color="black", lw=1)
    axes[0].set_xlabel("A–A overlap (µm)")
    axes[0].set_ylabel("passive restoring force (pN)")
    axes[1].axhline(0, color="black", lw=1)
    axes[1].set_xlabel("Z-disc–M-line distance (µm)")
    axes[1].set_ylabel("passive distance-increasing force (pN)")
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "passive_force_scans.png", dpi=180)
    plt.close(fig)


def plot_creep(traces: Dict[str, Dict[str, np.ndarray]]) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    for name, result in traces.items():
        axes[0].plot(result["time"], result["overlaps"][:, 0], label=name)
        axes[1].plot(
            result["time"],
            np.mean(result["overlaps"][:, 1:], axis=1),
            label=name,
        )
    axes[0].set_ylabel("A–A overlap (µm)")
    axes[1].set_ylabel("mean A–M overlap (µm)")
    axes[1].set_xlabel("time (s)")
    axes[0].legend(ncol=3, fontsize=8)
    fig.suptitle("Six-second turnover/creep diagnostic")
    fig.tight_layout()
    fig.savefig(OUT / "creep_timeseries.png", dpi=180)
    plt.close(fig)


def make_animation(name: str, result: Dict[str, np.ndarray], p: Parameters) -> None:
    stride = max(1, len(result["time"]) // 180)
    indices = np.arange(0, len(result["time"]), stride)
    fig, ax = plt.subplots(figsize=(10, 3.5))

    def draw(frame: int):
        ax.clear()
        idx = int(indices[frame])
        a_left, a_right, m_left, m_right = result["coords"][idx]
        y = result["overlaps"][idx, 0]
        x_left, x_right = result["overlaps"][idx, 1:]
        ax.plot([a_left - p.actin_length, a_left], [0.15, 0.15], lw=7, color="#d62728")
        ax.plot([a_right, a_right + p.actin_length], [-0.15, -0.15], lw=7, color="#ff7f0e")
        ax.plot(
            [m_left - p.myosin_half_length, m_left + p.myosin_half_length],
            [0.0, 0.0],
            lw=10,
            color="#1f77b4",
        )
        ax.plot(
            [m_right - p.myosin_half_length, m_right + p.myosin_half_length],
            [0.0, 0.0],
            lw=10,
            color="#1f77b4",
        )
        if y > 0:
            ax.axvspan(a_right, a_left, color="#2ca02c", alpha=0.25)
        ax.scatter([a_left, a_right], [0.15, -0.15], color="black", s=35, zorder=5)
        ax.set_xlim(-2.0, 2.0)
        ax.set_ylim(-0.45, 0.45)
        ax.set_yticks([])
        ax.set_xlabel("axial position (µm)")
        ax.set_title(
            f"{name}   t={result['time'][idx]:.2f} s   "
            f"A–A={y:.3f} µm   A–M=({x_left:.3f},{x_right:.3f}) µm   "
            f"ACTN2={result['bond_count'][idx]}"
        )
        return ()

    animation = FuncAnimation(fig, draw, frames=len(indices), interval=40, blit=False)
    animation.save(OUT / f"{name}.gif", writer=PillowWriter(fps=20))
    plt.close(fig)


def main() -> None:
    run_passive_force_scans()
    run_convergence()
    run_creep()
    manifest = {
        "mechanisms": list(MECHANISMS),
        "initial_conditions": INITIAL_CONDITIONS,
        "seeds": list(SEEDS),
        "base_parameters": asdict(Parameters()),
    }
    with (OUT / "campaign_manifest.json").open("w") as handle:
        json.dump(manifest, handle, indent=2)
    print(f"campaign complete: {OUT}")


if __name__ == "__main__":
    main()

