#!/usr/bin/env python3
"""Classify saved actin-actin bonds by status and relative polarity."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import h5py
import matplotlib.pyplot as plt
import numpy as np


CATEGORIES = (
    "status1_parallel",
    "status1_antiparallel",
    "status2_parallel",
    "status2_antiparallel",
)
COLORS = {
    "status1_parallel": "#e69f00",
    "status1_antiparallel": "#cc79a7",
    "status2_parallel": "#56b4e9",
    "status2_antiparallel": "#0072b2",
}
LABELS = {
    "status1_parallel": "status 1, parallel",
    "status1_antiparallel": "status 1, antiparallel",
    "status2_parallel": "status 2, parallel",
    "status2_antiparallel": "status 2, antiparallel",
}


def safe_fraction(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    result = np.zeros_like(numerator, dtype=float)
    np.divide(numerator, denominator, out=result, where=denominator > 0)
    return result


def classify_file(path: Path, chunk_size: int = 50) -> tuple[list[dict], dict[str, np.ndarray]]:
    rows: list[dict] = []
    late_cos = {"status1": [], "status2": []}
    with h5py.File(path, "r", locking=False) as handle:
        directions_ds = handle["/actin/direction"]
        pairs_ds = handle["/state/aa_pairs_current"]
        status_ds = handle["/state/aa_status_current"]
        steps = np.asarray(handle["/state/current_step"][:, 0], dtype=int)
        n_frames = min(len(directions_ds), len(pairs_ds), len(status_ds), len(steps))
        late_start = int(0.8 * n_frames)

        for start in range(0, n_frames, chunk_size):
            stop = min(start + chunk_size, n_frames)
            directions = np.asarray(directions_ds[start:stop], dtype=float)
            pairs = np.asarray(pairs_ds[start:stop], dtype=int)
            statuses = np.asarray(status_ds[start:stop, :, 0], dtype=np.int8)
            for local, frame in enumerate(range(start, stop)):
                valid = statuses[local] > 0
                frame_status = statuses[local, valid]
                frame_pairs = pairs[local, valid]
                if len(frame_pairs):
                    first = directions[local, frame_pairs[:, 0]]
                    second = directions[local, frame_pairs[:, 1]]
                    cos_angle = np.einsum("ij,ij->i", first, second)
                    parallel = cos_angle >= 0.0
                else:
                    cos_angle = np.empty(0, dtype=float)
                    parallel = np.empty(0, dtype=bool)
                counts = {
                    "status1_parallel": int(np.sum((frame_status == 1) & parallel)),
                    "status1_antiparallel": int(np.sum((frame_status == 1) & ~parallel)),
                    "status2_parallel": int(np.sum((frame_status == 2) & parallel)),
                    "status2_antiparallel": int(np.sum((frame_status == 2) & ~parallel)),
                }
                n_status1 = counts["status1_parallel"] + counts["status1_antiparallel"]
                n_status2 = counts["status2_parallel"] + counts["status2_antiparallel"]
                n_bonds = n_status1 + n_status2
                row = {
                    "replicate": path.stem.split("_r")[-1].split("_")[0],
                    "file": str(path),
                    "frame": frame,
                    "step": int(steps[frame]),
                    "n_bonds": n_bonds,
                    **counts,
                    "status1": n_status1,
                    "status2": n_status2,
                    "parallel": counts["status1_parallel"] + counts["status2_parallel"],
                    "antiparallel": counts["status1_antiparallel"] + counts["status2_antiparallel"],
                }
                rows.append(row)
                if frame >= late_start and len(cos_angle):
                    late_cos["status1"].append(cos_angle[frame_status == 1])
                    late_cos["status2"].append(cos_angle[frame_status == 2])

    return rows, {
        key: np.concatenate(values) if values else np.empty(0, dtype=float)
        for key, values in late_cos.items()
    }


def summarize(rep_rows: list[dict], window: str) -> dict:
    if window == "final":
        selected = rep_rows[-1:]
    elif window == "late20":
        selected = rep_rows[int(0.8 * len(rep_rows)):]
    else:
        raise ValueError(window)

    totals = {key: sum(row[key] for row in selected) for key in CATEGORIES}
    n_bonds = sum(totals.values())
    status1 = totals["status1_parallel"] + totals["status1_antiparallel"]
    status2 = totals["status2_parallel"] + totals["status2_antiparallel"]
    parallel = totals["status1_parallel"] + totals["status2_parallel"]
    antiparallel = totals["status1_antiparallel"] + totals["status2_antiparallel"]
    return {
        "replicate": rep_rows[0]["replicate"],
        "window": window,
        "n_frames": len(selected),
        "first_step": selected[0]["step"],
        "last_step": selected[-1]["step"],
        "mean_bonds_per_frame": n_bonds / len(selected),
        **{f"mean_{key}_per_frame": totals[key] / len(selected) for key in CATEGORIES},
        "fraction_status1": status1 / n_bonds if n_bonds else np.nan,
        "fraction_status2": status2 / n_bonds if n_bonds else np.nan,
        "fraction_parallel": parallel / n_bonds if n_bonds else np.nan,
        "fraction_antiparallel": antiparallel / n_bonds if n_bonds else np.nan,
        "status1_fraction_parallel": totals["status1_parallel"] / status1 if status1 else np.nan,
        "status1_fraction_antiparallel": totals["status1_antiparallel"] / status1 if status1 else np.nan,
        "status2_fraction_parallel": totals["status2_parallel"] / status2 if status2 else np.nan,
        "status2_fraction_antiparallel": totals["status2_antiparallel"] / status2 if status2 else np.nan,
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def make_plot(rows_by_rep: dict[str, list[dict]], late_cos_by_rep: dict[str, dict[str, np.ndarray]], output: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.4), constrained_layout=True)

    for replicate, rows in rows_by_rep.items():
        steps = np.array([row["step"] for row in rows], dtype=float)
        total = np.array([row["n_bonds"] for row in rows], dtype=float)
        for category in CATEGORIES:
            values = np.array([row[category] for row in rows], dtype=float)
            axes[0].plot(steps, values, color=COLORS[category], alpha=0.42, lw=0.8)
            axes[1].plot(steps, safe_fraction(values, total), color=COLORS[category], alpha=0.42, lw=0.8)

    # Overlay the replicate mean on their common saved steps.
    common_steps = sorted(set.intersection(*(
        {row["step"] for row in rows} for rows in rows_by_rep.values()
    )))
    for category in CATEGORIES:
        values = []
        fractions = []
        for rows in rows_by_rep.values():
            lookup = {row["step"]: row for row in rows}
            values.append([lookup[step][category] for step in common_steps])
            fractions.append([
                lookup[step][category] / lookup[step]["n_bonds"] if lookup[step]["n_bonds"] else 0.0
                for step in common_steps
            ])
        axes[0].plot(common_steps, np.mean(values, axis=0), color=COLORS[category], lw=2.0, label=LABELS[category])
        axes[1].plot(common_steps, np.mean(fractions, axis=0), color=COLORS[category], lw=2.0)

    bins = np.linspace(-1.0, 1.0, 51)
    combined_cos = {
        status: np.concatenate([late_cos[status] for late_cos in late_cos_by_rep.values()])
        for status in ("status1", "status2")
    }
    axes[2].hist(combined_cos["status1"], bins=bins, density=True, histtype="step", lw=2.2,
                 color="#b24c91", label="status 1")
    axes[2].hist(combined_cos["status2"], bins=bins, density=True, histtype="step", lw=2.2,
                 color="#0072b2", label="status 2")

    axes[0].set(title="Bond counts", xlabel="step", ylabel="number of bonds")
    axes[1].set(title="Composition of current bonds", xlabel="step", ylabel="fraction of bonds", ylim=(0, 1))
    axes[2].set(title="Late 20%: relative orientation", xlabel=r"$\hat u_i\cdot\hat u_j$", ylabel="probability density", xlim=(-1, 1))
    axes[0].legend(frameon=False, fontsize=8)
    axes[2].legend(frameon=False)
    fig.suptitle(r"$N_A=2100$, $\lambda=0$: actin--actin bond types (two completed replicas)")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    all_rows: list[dict] = []
    rows_by_rep: dict[str, list[dict]] = {}
    late_cos_by_rep: dict[str, dict[str, np.ndarray]] = {}
    for path in args.inputs:
        rows, late_cos = classify_file(path)
        replicate = rows[0]["replicate"]
        rows_by_rep[replicate] = rows
        late_cos_by_rep[replicate] = late_cos
        all_rows.extend(rows)

    summary_rows = []
    for rows in rows_by_rep.values():
        summary_rows.extend([summarize(rows, "final"), summarize(rows, "late20")])

    write_csv(args.output_dir / "aa_bond_types_by_frame.csv", all_rows)
    write_csv(args.output_dir / "aa_bond_types_summary.csv", summary_rows)
    make_plot(rows_by_rep, late_cos_by_rep, args.output_dir / "aa_bond_types.png")


if __name__ == "__main__":
    main()
