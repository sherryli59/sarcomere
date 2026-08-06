#!/usr/bin/env python3
"""Search for a basin-independent Johann-entropic + titin regime."""

from __future__ import annotations

import csv
from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from overlap_model import Parameters, default_state, simulate
from run_campaign import INITIAL_CONDITIONS, tail_mean, tail_slope


OUT = Path(__file__).resolve().parent / "results"
ENTROPIC_COUNTS = (0.0, 5.0, 10.0, 20.0, 40.0, 60.0, 80.0)
TITIN_STIFFNESSES = (100.0, 220.0, 400.0)
SEEDS = (30, 31, 32, 33)


def main() -> None:
    rows = []
    for count in ENTROPIC_COUNTS:
        for titin_k in TITIN_STIFFNESSES:
            ic_rows = []
            for ic_name, (y0, d0) in INITIAL_CONDITIONS.items():
                for seed in SEEDS:
                    p = Parameters(
                        dt=5.0e-4,
                        duration=4.0,
                        save_every=5,
                        entropic_nkbt=count,
                        entropic_force_cap=8.0,
                        titin_k=titin_k,
                    )
                    result = simulate(p, seed, default_state(y0, d0))
                    y = result["overlaps"][:, 0]
                    x = np.mean(result["overlaps"][:, 1:], axis=1)
                    ic_rows.append(
                        {
                            "initial": ic_name,
                            "aa": tail_mean(y),
                            "am": tail_mean(x),
                            "aa_v": tail_slope(result["time"], y),
                            "am_v": tail_slope(result["time"], x),
                        }
                    )
            aa_by_ic = [
                np.mean([r["aa"] for r in ic_rows if r["initial"] == ic])
                for ic in INITIAL_CONDITIONS
            ]
            am_by_ic = [
                np.mean([r["am"] for r in ic_rows if r["initial"] == ic])
                for ic in INITIAL_CONDITIONS
            ]
            velocities = [abs(r["aa_v"]) for r in ic_rows] + [
                abs(r["am_v"]) for r in ic_rows
            ]
            row = {
                "entropic_nkbt": count,
                "titin_k_pN_um": titin_k,
                "mean_aa_um": np.mean(aa_by_ic),
                "mean_am_um": np.mean(am_by_ic),
                "aa_ic_spread_um": max(aa_by_ic) - min(aa_by_ic),
                "am_ic_spread_um": max(am_by_ic) - min(am_by_ic),
                "max_abs_seed_tail_velocity_um_s": max(velocities),
            }
            row["strict_pass"] = (
                row["aa_ic_spread_um"] < 0.025
                and row["am_ic_spread_um"] < 0.025
                and row["max_abs_seed_tail_velocity_um_s"] < 0.02
                and 0.06 < row["mean_aa_um"] < 0.14
                and 0.02 < row["mean_am_um"] < 0.655
            )
            rows.append(row)
            print(row, flush=True)

    path = OUT / "entropic_titin_sweep.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for titin_k in TITIN_STIFFNESSES:
        subset = [r for r in rows if r["titin_k_pN_um"] == titin_k]
        x = [r["entropic_nkbt"] for r in subset]
        axes[0].plot(x, [r["mean_aa_um"] for r in subset], "o-", label=f"k={titin_k:g}")
        axes[1].plot(x, [r["aa_ic_spread_um"] for r in subset], "o-")
        axes[2].plot(x, [r["am_ic_spread_um"] for r in subset], "o-")
    axes[0].axhline(0.10, color="black", ls=":")
    axes[0].set_ylabel("mean A–A overlap (µm)")
    axes[1].axhline(0.025, color="black", ls=":")
    axes[1].set_ylabel("A–A initial-condition spread (µm)")
    axes[2].axhline(0.025, color="black", ls=":")
    axes[2].set_ylabel("A–M initial-condition spread (µm)")
    for ax in axes:
        ax.set_xlabel("entropic population N")
    axes[0].legend()
    fig.tight_layout()
    fig.savefig(OUT / "entropic_titin_sweep.png", dpi=180)


if __name__ == "__main__":
    main()

