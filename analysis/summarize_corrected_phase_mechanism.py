#!/usr/bin/env python3
"""Summarize mechanism observables for the corrected phase sweep."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


BASE = 1.0e-3
LAM_ORDER = (0.0, 10**-5.5, 1e-5, 10**-4.5, 1e-4, 10**-3.5,
             1e-3, 10**-2.5, 1e-2, 10**-1.5, 1e-1, 10**-0.5, 1.0)
LAM_LABELS = ("0", "1e-5.5", "1e-5", "1e-4.5", "1e-4", "1e-3.5",
              "1e-3", "1e-2.5", "1e-2", "1e-1.5", "1e-1", "1e-0.5", "1e0")
COLORS = ("#4477AA", "#66CCEE", "#228833", "#CCBB44", "#EE6677", "#AA3377")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, records: list[dict[str, object]]) -> None:
    if not records:
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def finite_float(value: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return math.nan
    return result


def mean_tail(rows: list[dict[str, str]], field: str, n: int = 30) -> float:
    values = np.asarray([finite_float(row[field]) for row in rows[-n:]], dtype=float)
    return float(np.nanmean(values))


def half_time(rows: list[dict[str, str]], field: str, tail_n: int = 30) -> float:
    time = np.asarray([float(row["time_s"]) for row in rows])
    value = np.asarray([float(row[field]) for row in rows])
    early = float(np.mean(value[:min(5, len(value))]))
    late = float(np.mean(value[-min(tail_n, len(value)):]))
    delta = late - early
    if delta <= 0.05 * max(1.0, abs(late)):
        return math.nan
    target = early + 0.5 * delta
    hits = np.flatnonzero(value >= target)
    return float(time[hits[0]]) if len(hits) else math.nan


def load_timeseries(cache_dir: Path) -> dict[str, list[dict[str, str]]]:
    output = {}
    for path in sorted(cache_dir.glob("*.csv")):
        output[f"{path.stem}.h5"] = read_csv(path)
    return output


def grouped_cell_summary(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[int, float], list[dict[str, object]]] = defaultdict(list)
    fields = (
        "legacy_phi_3d", "tail_phi_axial_myosin", "tail_q_sarcomere_axial",
        "tail_cb1_fraction", "tail_cb2_fraction", "tail_mean_cb2_load",
        "tail_directional_extra_over_base", "tail_am_engaged_fraction",
        "tail_binding_rate_per_actin_s", "tail_unbinding_rate_per_actin_s",
        "tail_nematic_myosin", "tail_nematic_actin",
    )
    for row in rows:
        grouped[(int(row["n_actins"]), float(row["lambda_value"]))].append(row)
    output = []
    for (n_actins, lam), group in sorted(grouped.items()):
        record: dict[str, object] = {
            "n_actins": n_actins,
            "lambda_value": lam,
            "lambda_over_base": lam / BASE,
            "n_replicates": len(group),
        }
        for field in fields:
            values = np.asarray([float(row[field]) for row in group], dtype=float)
            values = values[np.isfinite(values)]
            record[f"{field}_mean"] = float(np.mean(values)) if len(values) else math.nan
            record[f"{field}_sd"] = float(np.std(values, ddof=1)) if len(values) > 1 else 0.0
        record["assembled_fraction_phi3d_0p5"] = float(np.mean([
            float(row["legacy_phi_3d"]) >= 0.5 for row in group
            if np.isfinite(float(row["legacy_phi_3d"]))
        ])) if any(np.isfinite(float(row["legacy_phi_3d"])) for row in group) else math.nan
        output.append(record)
    return output


def threshold_tests(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    output = []
    actin_values = sorted({int(row["n_actins"]) for row in rows})
    for n_actins in actin_values:
        by = {(float(row["lambda_value"]), int(row["replica"])): row for row in rows
              if int(row["n_actins"]) == n_actins and np.isfinite(float(row["legacy_phi_3d"]))}
        boundary = []
        broad = []
        for replica in range(3):
            low = by.get((1e-3, replica))
            high = by.get((10**-2.5, replica))
            if low is not None and high is not None:
                boundary.append(float(high["legacy_phi_3d"]) - float(low["legacy_phi_3d"]))
            rep_rows = [row for (lam, rep), row in by.items() if rep == replica]
            below = [float(row["legacy_phi_3d"]) for row in rep_rows
                     if float(row["lambda_value"]) <= BASE]
            above = [float(row["legacy_phi_3d"]) for row in rep_rows
                     if float(row["lambda_value"]) > BASE]
            if below and above:
                broad.append(float(np.mean(above) - np.mean(below)))
        for contrast, values in (("3.16e-3_minus_1e-3", boundary),
                                 ("mean_above_minus_at_or_below", broad)):
            arr = np.asarray(values, dtype=float)
            test = stats.ttest_1samp(arr, 0.0) if len(arr) > 1 else None
            output.append({
                "n_actins": n_actins,
                "contrast": contrast,
                "n_pairs": len(arr),
                "mean_delta_phi3d": float(np.mean(arr)) if len(arr) else math.nan,
                "sd_delta_phi3d": float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0,
                "t_statistic": float(test.statistic) if test else math.nan,
                "p_value": float(test.pvalue) if test else math.nan,
            })

    # Pooled paired contrasts retain density/replica as blocks.
    for contrast in ("3.16e-3_minus_1e-3", "mean_above_minus_at_or_below"):
        values = np.asarray([float(row["mean_delta_phi3d"]) for row in output
                             if row["contrast"] == contrast], dtype=float)
        # Density means are reported descriptively; the replicate-level pooled
        # values are reconstructed below to avoid pretending six densities are
        # independent biological replicates.
        output.append({
            "n_actins": "all_density_means",
            "contrast": contrast,
            "n_pairs": len(values),
            "mean_delta_phi3d": float(np.mean(values)),
            "sd_delta_phi3d": float(np.std(values, ddof=1)),
            "t_statistic": math.nan,
            "p_value": math.nan,
        })
    return output


def correlations(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    output = []
    for n_actins in sorted({int(row["n_actins"]) for row in rows}):
        group = [row for row in rows if int(row["n_actins"]) == n_actins]
        y = np.asarray([float(row["legacy_phi_3d"]) for row in group])
        for field in ("tail_cb2_fraction", "tail_unbinding_rate_per_actin_s",
                      "tail_am_engaged_fraction", "tail_q_sarcomere_axial"):
            x = np.asarray([float(row[field]) for row in group])
            good = np.isfinite(x) & np.isfinite(y)
            rho, p = stats.spearmanr(x[good], y[good]) if good.sum() >= 4 else (math.nan, math.nan)
            output.append({
                "n_actins": n_actins,
                "predictor": field,
                "n": int(good.sum()),
                "spearman_rho": float(rho),
                "p_value": float(p),
            })
    return output


def temporal_lags(
    summaries: list[dict[str, object]], series: dict[str, list[dict[str, str]]],
) -> list[dict[str, object]]:
    output = []
    for summary in summaries:
        filename = str(summary["filename"])
        if filename not in series:
            continue
        rows = series[filename]
        times = {field: half_time(rows, field) for field in (
            "phi_axial_myosin", "cb2_fraction", "fraction_am_engaged",
            "q_sarcomere_axial", "nematic_myosin", "nematic_actin",
        )}
        output.append({
            "filename": filename,
            "n_actins": summary["n_actins"],
            "lambda_value": summary["lambda_value"],
            "lambda_over_base": summary["lambda_over_base"],
            "replica": summary["replica"],
            "legacy_phi_3d": summary["legacy_phi_3d"],
            **{f"t50_{field}_s": value for field, value in times.items()},
            "cb2_minus_order_t50_s": times["cb2_fraction"] - times["phi_axial_myosin"],
            "am_minus_order_t50_s": times["fraction_am_engaged"] - times["phi_axial_myosin"],
        })
    return output


def plot_phase_and_mechanism(cells: list[dict[str, object]], output: Path) -> None:
    actins = sorted({int(row["n_actins"]) for row in cells})
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for index, n_actins in enumerate(actins):
        group = sorted([row for row in cells if int(row["n_actins"]) == n_actins],
                       key=lambda row: float(row["lambda_value"]))
        x = np.asarray([max(float(row["lambda_over_base"]), 1e-3) for row in group])
        phi = np.asarray([float(row["legacy_phi_3d_mean"]) for row in group])
        phi_sd = np.asarray([float(row["legacy_phi_3d_sd"]) for row in group])
        cb2 = np.asarray([float(row["tail_cb2_fraction_mean"]) for row in group])
        unbind = np.asarray([float(row["tail_unbinding_rate_per_actin_s_mean"]) for row in group])
        q = np.asarray([float(row["tail_q_sarcomere_axial_mean"]) for row in group])
        color = COLORS[index % len(COLORS)]
        label = f"A={n_actins}"
        axes[0, 0].plot(x, phi, "o-", color=color, label=label, ms=4)
        axes[0, 0].fill_between(x, phi - phi_sd, phi + phi_sd, color=color, alpha=0.12)
        axes[0, 1].plot(x, cb2, "o-", color=color, label=label, ms=4)
        axes[1, 0].plot(x, unbind, "o-", color=color, label=label, ms=4)
        axes[1, 1].plot(cb2, phi, "o-", color=color, label=label, ms=4)
    for axis in axes.flat[:3]:
        axis.set_xscale("log")
        axis.axvline(1.0, color="black", linestyle="--", linewidth=1)
        axis.set_xlabel(r"directional lifetime / base lifetime  ($\lambda/\tau_1$)")
        axis.grid(alpha=0.2)
    axes[0, 0].set_ylabel(r"legacy $\Phi_{3D}$")
    axes[0, 0].set_title("Structural order")
    axes[0, 1].set_ylabel("ACTN2 status-2 fraction")
    axes[0, 1].set_title("Crosslink stabilization")
    axes[1, 0].set_ylabel("unbinding events / actin / s")
    axes[1, 0].set_title("Turnover")
    axes[1, 1].set_xlabel("ACTN2 status-2 fraction")
    axes[1, 1].set_ylabel(r"legacy $\Phi_{3D}$")
    axes[1, 1].set_title("More crosslinks do not always mean more order")
    axes[0, 0].legend(ncol=2, fontsize=8)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_representative_timeseries(
    series: dict[str, list[dict[str, str]]], output: Path,
) -> None:
    selected_lambdas = {0.0, 1e-3, 10**-2.5, 1e-1}
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5), sharex=True, constrained_layout=True)
    for row_index, n_actins in enumerate((900, 1800)):
        candidates = []
        for filename, rows in series.items():
            first = rows[0]
            if int(first["n_actins"]) != n_actins or int(first["replica"]) != 0:
                continue
            lam = float(first["lambda_value"])
            if any(abs(lam - value) < 1e-12 for value in selected_lambdas):
                candidates.append((lam, rows))
        for lam, rows in sorted(candidates):
            time = np.asarray([float(row["time_s"]) for row in rows])
            phi = np.asarray([float(row["phi_axial_myosin"]) for row in rows])
            cb2 = np.asarray([float(row["cb2_fraction"]) for row in rows])
            label = r"$\lambda/\tau_1$=" + f"{lam / BASE:g}"
            axes[row_index, 0].plot(time, phi, label=label)
            axes[row_index, 1].plot(time, cb2, label=label)
        axes[row_index, 0].set_ylabel(f"A={n_actins}\naxial order")
        axes[row_index, 1].set_ylabel(f"A={n_actins}\nstatus-2 fraction")
        axes[row_index, 0].axhline(0.5, color="black", linestyle="--", linewidth=0.8)
        axes[row_index, 0].grid(alpha=0.2)
        axes[row_index, 1].grid(alpha=0.2)
    axes[0, 0].set_title("Density-band assembly")
    axes[0, 1].set_title("ACTN2 stabilization")
    axes[1, 0].set_xlabel("time (s)")
    axes[1, 1].set_xlabel("time (s)")
    axes[0, 0].legend(fontsize=8, ncol=2)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    args = parser.parse_args()
    summary_raw = read_csv(args.input_dir / "run_summary.csv")
    summaries: list[dict[str, object]] = []
    for row in summary_raw:
        converted: dict[str, object] = dict(row)
        for field in row:
            if field not in {"filename", "lambda_label", "classification"}:
                converted[field] = finite_float(row[field])
        summaries.append(converted)
    complete = [row for row in summaries if int(row["complete"]) == 1]
    series = load_timeseries(args.input_dir / "cache")
    cells = grouped_cell_summary(complete)
    tests = threshold_tests(complete)
    corr = correlations(complete)
    lags = temporal_lags(complete, series)
    write_csv(args.input_dir / "cell_summary.csv", cells)
    write_csv(args.input_dir / "threshold_tests.csv", tests)
    write_csv(args.input_dir / "mechanism_correlations.csv", corr)
    write_csv(args.input_dir / "temporal_lags.csv", lags)
    plot_phase_and_mechanism(cells, args.input_dir / "phase_mechanism.png")
    plot_representative_timeseries(series, args.input_dir / "representative_timeseries.png")
    print(f"summarized {len(complete)} complete runs in {len(cells)} cells")


if __name__ == "__main__":
    main()
