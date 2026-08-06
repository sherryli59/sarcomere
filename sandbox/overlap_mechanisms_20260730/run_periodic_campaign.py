#!/usr/bin/env python3
"""Campaign on the exact periodic 2-myosin/4-actin topology."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from periodic_model import (
    aa_overlaps,
    am_overlaps,
    default_periodic_state,
    force_components,
    periodic_parameters,
    simulate_periodic,
)
from run_campaign import tail_mean, tail_slope


OUT = Path(__file__).resolve().parent / "results"
MECHANISMS = (
    "baseline_turnover",
    "permanent_actn2",
    "johann_adhesion",
    "johann_entropic",
    "titin",
    "entropic_plus_titin",
    "register_plus_titin",
)
INITIALS = {
    "reference": (0.10, 4.30),
    "aa_low": (0.04, 4.30),
    "aa_high": (0.16, 4.30),
    "am_deep": (0.10, 3.90),
    "am_shallow": (0.10, 4.70),
}
SEEDS = range(6)


def main() -> None:
    # Exact-topology mechanics checks.
    for name in MECHANISMS:
        p = periodic_parameters(name)
        state = default_periodic_state(0.10, p.box_length)
        for component, force in force_components(state, p).items():
            assert abs(float(np.sum(force))) < 1e-10, (name, component, force)
        assert np.allclose(aa_overlaps(state, p), 0.10)
        assert np.ptp(am_overlaps(state, p)) < 1e-12

    rows = []
    representative = {}
    for name in MECHANISMS:
        for initial_name, (y0, length) in INITIALS.items():
            p = periodic_parameters(name, length)
            for seed in SEEDS:
                result = simulate_periodic(
                    p, seed, default_periodic_state(y0, length)
                )
                y = np.mean(result["overlaps"][:, :2], axis=1)
                x = np.mean(result["overlaps"][:, 2:], axis=1)
                rows.append(
                    {
                        "mechanism": name,
                        "initial": initial_name,
                        "seed": seed,
                        "aa_tail_um": tail_mean(y),
                        "am_tail_um": tail_mean(x),
                        "aa_tail_velocity_um_s": tail_slope(result["time"], y),
                        "am_tail_velocity_um_s": tail_slope(result["time"], x),
                    }
                )
                if initial_name == "reference" and seed == 0:
                    representative[name] = result
                print(name, initial_name, seed, flush=True)
    with (OUT / "periodic_runs.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    summary = []
    for name in MECHANISMS:
        subset = [r for r in rows if r["mechanism"] == name]
        aa_ic = []
        am_ic = []
        aa_v = []
        am_v = []
        for initial_name in INITIALS:
            group = [r for r in subset if r["initial"] == initial_name]
            aa_ic.append(np.mean([r["aa_tail_um"] for r in group]))
            am_ic.append(np.mean([r["am_tail_um"] for r in group]))
            aa_v.append(np.mean([r["aa_tail_velocity_um_s"] for r in group]))
            am_v.append(np.mean([r["am_tail_velocity_um_s"] for r in group]))
        row = {
            "mechanism": name,
            "mean_aa_um": np.mean(aa_ic),
            "mean_am_um": np.mean(am_ic),
            "aa_ic_spread_um": max(aa_ic) - min(aa_ic),
            "am_ic_spread_um": max(am_ic) - min(am_ic),
            "max_abs_group_tail_velocity_um_s": max(
                max(abs(np.asarray(aa_v))), max(abs(np.asarray(am_v)))
            ),
        }
        row["strict_pass"] = (
            row["aa_ic_spread_um"] < 0.025
            and row["am_ic_spread_um"] < 0.025
            and row["max_abs_group_tail_velocity_um_s"] < 0.01
            and 0.06 < row["mean_aa_um"] < 0.14
            and 0.02 < row["mean_am_um"] < 0.655
        )
        summary.append(row)
    with (OUT / "periodic_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    for name, result in representative.items():
        axes[0].plot(result["time"], np.mean(result["overlaps"][:, :2], axis=1), label=name)
        axes[1].plot(result["time"], np.mean(result["overlaps"][:, 2:], axis=1), label=name)
    axes[0].axhline(0.10, color="black", ls=":")
    axes[0].set_ylabel("mean A–A overlap (µm)")
    axes[1].set_ylabel("mean A–M overlap (µm)")
    axes[1].set_xlabel("time (s)")
    axes[0].legend(ncol=3, fontsize=8)
    fig.suptitle("Exact 2-myosin/4-actin periodic topology")
    fig.tight_layout()
    fig.savefig(OUT / "periodic_timeseries.png", dpi=180)


if __name__ == "__main__":
    main()

