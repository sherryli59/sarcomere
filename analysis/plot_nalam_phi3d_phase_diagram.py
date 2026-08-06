#!/usr/bin/env python3
"""Render the unannotated n_actins--lambda legacy Phi_3D phase diagram."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm


LAMBDA_LABELS = (
    "0", "1e-5.5", "1e-5", "1e-4.5", "1e-4", "1e-3.5", "1e-3",
    "1e-2.5", "1e-2", "1e-1.5", "1e-1", "1e-0.5", "1e0",
)
ACTIN_COUNTS = (600, 900, 1200, 1500, 1800, 2100)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    values = {(int(row["n_actins"]), row["lambda_label"]): float(row["Phi_3D_mean"])
              for row in csv.DictReader(args.input.open())}
    data = np.full((len(ACTIN_COUNTS), len(LAMBDA_LABELS)), np.nan)
    for yi, n_actins in enumerate(ACTIN_COUNTS):
        for xi, label in enumerate(LAMBDA_LABELS):
            data[yi, xi] = values.get((n_actins, label), np.nan)

    fig, axis = plt.subplots(figsize=(11.5, 5.0), constrained_layout=True)
    image = axis.imshow(
        data,
        origin="lower",
        aspect="auto",
        cmap="coolwarm",
        norm=TwoSlopeNorm(vmin=0.0, vcenter=0.5, vmax=1.0),
        interpolation="nearest",
    )
    axis.set_xticks(range(len(LAMBDA_LABELS)), LAMBDA_LABELS, rotation=45, ha="right")
    axis.set_yticks(range(len(ACTIN_COUNTS)), ACTIN_COUNTS)
    axis.set_xlabel(r"$\lambda$")
    axis.set_ylabel(r"$n_{\mathrm{actin}}$")
    colorbar = fig.colorbar(image, ax=axis, pad=0.02)
    colorbar.set_label(r"mean legacy $\Phi_{3D}$")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    main()
