#!/usr/bin/env python3
"""Compare recorded A--A lifetimes with the load used by the kinetic law.

The corrected legacy sweep did not write load into each completed-lifetime
event.  We therefore combine two complementary observables:

1. exact completed status-2 lifetimes from ``catch_bond/completed_lifetimes``;
2. instantaneous bond load and current age reconstructed at saved checkpoints
   from active pair identities, attachment steps, directions, and f_load_cb.

The second sample is interval/length biased and is used only to diagnose the
available load range, not as an uncensored lifetime estimator.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


PATTERN = re.compile(r"A(?P<A>\d+)_M1000_s2fstab_lt(?P<label>[^_]+)_r(?P<rep>\d+)_L12x2x2\.h5")
BASE = 1.0e-3
DT = 1.0e-5


def lambda_value(label: str) -> float:
    if label == "0":
        return 0.0
    mantissa, exponent = label.split("e", 1)
    return float(mantissa) * 10.0 ** float(exponent)


def exact_geometric_median(mean_lifetime: float) -> float:
    probability = min(1.0, DT / mean_lifetime)
    return math.ceil(math.log(0.5) / math.log1p(-probability)) * DT


def sampled_load_and_age(handle: h5py.File, frame_stride: int) -> tuple[np.ndarray, np.ndarray]:
    counts = np.asarray(handle["/state/aa_pairs_current_count"]).reshape(-1)
    steps = np.asarray(handle["/state/current_step"]).reshape(-1)
    loads_out: list[np.ndarray] = []
    ages_out: list[np.ndarray] = []
    for frame in range(1, len(counts), frame_stride):
        count = int(counts[frame])
        if count <= 0:
            continue
        pairs = np.asarray(handle["/state/aa_pairs_current"][frame, :count], dtype=int)
        status = np.asarray(handle["/state/aa_status_current"][frame, :count]).reshape(-1)
        attach = np.asarray(handle["/state/aa_attach_step_current"][frame, :count]).reshape(-1)
        keep = (status == 2) & (attach >= 0)
        if not np.any(keep):
            continue
        pairs = pairs[keep]
        attach = attach[keep]
        actin_load = np.asarray(handle["/actin/f_load_cb"][frame]).reshape(-1)
        direction = np.asarray(handle["/actin/direction"][frame])
        cosine = np.abs(np.einsum(
            "ij,ij->i", direction[pairs[:, 0]], direction[pairs[:, 1]], optimize=True,
        ))
        pair_load = cosine * np.minimum(actin_load[pairs[:, 0]], actin_load[pairs[:, 1]])
        current_age = (int(steps[frame]) - attach) * DT
        valid = (pair_load >= 0.0) & np.isfinite(pair_load) & (current_age >= 0.0)
        loads_out.append(pair_load[valid])
        ages_out.append(current_age[valid])
    return np.concatenate(loads_out), np.concatenate(ages_out)


def summarize(path: Path, frame_stride: int) -> tuple[dict[str, object], list[dict[str, object]]]:
    match = PATTERN.match(path.name)
    if match is None:
        raise ValueError(path.name)
    lam = lambda_value(match["label"])
    with h5py.File(path, "r", locking=False) as handle:
        completed = np.asarray(handle["/catch_bond/completed_lifetimes"], dtype=float)
        load, age = sampled_load_and_age(handle, frame_stride)
    load_quantiles = np.quantile(load, [0.01, 0.1, 0.5, 0.9, 0.99])
    lifetime_quantiles = np.quantile(completed, [0.1, 0.5, 0.9, 0.99])
    prescribed_mean = BASE + lam * float(np.mean(load))
    rho, rho_p = stats.spearmanr(load, age)
    record: dict[str, object] = {
        "filename": path.name,
        "n_actins": int(match["A"]),
        "replica": int(match["rep"]),
        "lambda_label": match["label"],
        "lambda_value": lam,
        "lambda_over_base": lam / BASE,
        "n_completed": len(completed),
        "completed_mean_s": float(np.mean(completed)),
        "completed_q10_s": float(lifetime_quantiles[0]),
        "completed_median_s": float(lifetime_quantiles[1]),
        "completed_q90_s": float(lifetime_quantiles[2]),
        "completed_q99_s": float(lifetime_quantiles[3]),
        "n_status2_checkpoint_observations": len(load),
        "bond_load_mean": float(np.mean(load)),
        "bond_load_p01": float(load_quantiles[0]),
        "bond_load_p10": float(load_quantiles[1]),
        "bond_load_median": float(load_quantiles[2]),
        "bond_load_p90": float(load_quantiles[3]),
        "bond_load_p99": float(load_quantiles[4]),
        "fraction_bond_load_gt_0p998": float(np.mean(load > 0.998)),
        "current_age_median_s": float(np.median(age)),
        "spearman_load_vs_current_age": float(rho),
        "spearman_p_value": float(rho_p),
        "prescribed_mean_lifetime_s": prescribed_mean,
        "prescribed_geometric_median_s": exact_geometric_median(prescribed_mean),
        "observed_over_prescribed_median": float(np.median(completed)) / exact_geometric_median(prescribed_mean),
    }

    # Descriptive, checkpoint-sampled age by load bin. These are repeated,
    # length-biased observations and are not treated as completed lifetimes.
    edges = np.asarray([0.0, 0.9, 0.98, 0.99, 0.995, 0.998, 0.999, 0.9995, 1.000001])
    binned = []
    for left, right in zip(edges[:-1], edges[1:]):
        selected = (load >= left) & (load < right)
        if not np.any(selected):
            continue
        binned.append({
            "filename": path.name,
            "lambda_value": lam,
            "lambda_over_base": lam / BASE,
            "load_left": left,
            "load_right": right,
            "n_observations": int(np.sum(selected)),
            "load_mean": float(np.mean(load[selected])),
            "current_age_median_s": float(np.median(age[selected])),
            "current_age_mean_s": float(np.mean(age[selected])),
        })
    return record, binned


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot(rows: list[dict[str, object]], output: Path) -> None:
    rows = sorted(rows, key=lambda row: float(row["lambda_value"]))
    prescribed = np.asarray([float(row["prescribed_mean_lifetime_s"]) for row in rows])
    observed = np.asarray([float(row["completed_median_s"]) for row in rows])
    q10 = np.asarray([float(row["completed_q10_s"]) for row in rows])
    q90 = np.asarray([float(row["completed_q90_s"]) for row in rows])
    theoretical = np.asarray([float(row["prescribed_geometric_median_s"]) for row in rows])
    ratio = np.asarray([float(row["lambda_over_base"]) for row in rows])
    load_mean = np.asarray([float(row["bond_load_mean"]) for row in rows])
    load_p10 = np.asarray([float(row["bond_load_p10"]) for row in rows])
    load_p90 = np.asarray([float(row["bond_load_p90"]) for row in rows])

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5), constrained_layout=True)
    axes[0].fill_between(prescribed, q10, q90, alpha=0.18, label="completed 10--90%")
    axes[0].plot(prescribed, observed, "o-", label="observed completed median")
    axes[0].plot(prescribed, theoretical, "--", label="kinetic-law median")
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_xlabel("prescribed mean lifetime: base + lambda <load> (s)")
    axes[0].set_ylabel("completed A--A lifetime (s)")
    axes[0].set_title("Observed lifetime follows the law, then falls below it")
    axes[0].grid(alpha=0.2, which="both")
    axes[0].legend(fontsize=8)

    x = np.arange(len(rows))
    axes[1].fill_between(x, load_p10, load_p90, alpha=0.18, label="10--90%")
    axes[1].plot(x, load_mean, "o-", label="mean active status-2 load")
    axes[1].set_xticks(x[::2], [f"{value:g}" for value in ratio[::2]], rotation=35)
    axes[1].set_ylim(min(0.97, float(np.min(load_p10)) - 0.005), 1.0005)
    axes[1].set_xlabel("lambda / base")
    axes[1].set_ylabel("instantaneous pair load")
    axes[1].set_title("The load classifier is nearly saturated")
    axes[1].grid(alpha=0.2)
    axes[1].legend(fontsize=8)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path(
        "/scratch/users/xsli/sarcomere/data_nalam_s2fstab_base1e-3"))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--n-actins", type=int, default=2100)
    parser.add_argument("--replica", type=int, default=0)
    parser.add_argument("--frame-stride", type=int, default=10)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    files = []
    for path in args.input_dir.glob("*.h5"):
        match = PATTERN.match(path.name)
        if match and int(match["A"]) == args.n_actins and int(match["rep"]) == args.replica:
            files.append(path)
    files.sort(key=lambda path: lambda_value(PATTERN.match(path.name)["label"]))
    summaries = []
    binned = []
    for index, path in enumerate(files, 1):
        summary, bins = summarize(path, args.frame_stride)
        summaries.append(summary)
        binned.extend(bins)
        print(f"[{index}/{len(files)}] {path.name}", flush=True)
    write_csv(args.output_dir / "aa_lifetime_vs_load_summary.csv", summaries)
    write_csv(args.output_dir / "aa_checkpoint_age_by_load_bin.csv", binned)
    plot(summaries, args.output_dir / "aa_lifetime_vs_load.png")


if __name__ == "__main__":
    main()
