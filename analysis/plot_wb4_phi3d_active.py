#!/usr/bin/env python3
"""Plot stride-subsampled phi_3D trajectories for active A8000 12x4x4 runs."""

from __future__ import annotations

import argparse
import csv
import os
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/xsli-matplotlib")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def lambda_value(label: str) -> float:
    if "e" in label:
        mantissa, exponent = label.split("e", 1)
        return float(mantissa) * 10.0 ** float(exponent)
    return float(label)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    groups: dict[tuple[float, str], list[dict[str, str]]] = defaultdict(list)
    with args.csv.open(newline="") as stream:
        for row in csv.DictReader(stream):
            filename = row["filename"]
            if filename.startswith("A8000_M4000_") and filename.endswith("_L12x4x4.h5"):
                label = row["lambda"]
                groups[(lambda_value(label), label)].append(row)
    if not groups:
        raise SystemExit("No active A8000/M4000 12x4x4 rows found")
    metrics = (
        ("phi_3d_myosin", r"Myosin $\Phi_{3D}$"),
        ("phi_3d_actin_cb2", r"CB2 actin $\Phi_{3D}$"),
        ("shared_phi_3d_myosin_corrected", r"Shared-mode myosin $\Phi_{3D}$ (corrected)"),
        ("shared_phi_3d_actin_cb2_corrected", r"Shared-mode CB2 actin $\Phi_{3D}$ (corrected)"),
    )

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.5), sharex=True, constrained_layout=True)
    colors = plt.cm.viridis_r(
        [i / max(1, len(groups) - 1) for i in range(len(groups))]
    )
    for color, ((_, label), rows) in zip(colors, sorted(groups.items(), reverse=True)):
        rows.sort(key=lambda row: int(row["step"]))
        steps = [int(row["step"]) for row in rows]
        legend_label = rf"$\lambda={label}$"
        for ax, (column, title) in zip(axes.flat, metrics):
            values = [float(row[column]) for row in rows]
            ax.plot(steps, values, lw=1.7, color=color, label=legend_label)
            ax.set_title(title)
            ax.grid(alpha=0.22)
            ax.set_ylabel(r"$\Phi_{3D}$")

    for ax in axes[-1, :]:
        ax.set_xlabel("Simulation step")
    axes[0, 0].legend(ncol=2, fontsize=8, frameon=False)
    fig.suptitle(r"Active A8000/M4000, $12\times4\times4$: subsampled $\Phi_{3D}$ trajectories")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=240)
    plt.close(fig)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
