#!/usr/bin/env python3
"""Plot final-frame x-y positions, coloring each population by z coordinate."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--frame", type=int, default=-1)
    parser.add_argument("--phi", type=float, default=None)
    args = parser.parse_args()

    with h5py.File(args.input, "r", locking=False) as handle:
        n_frames = min(
            handle["/myosin/center"].shape[0],
            handle["/actin/center"].shape[0],
            handle["/actin/direction"].shape[0],
        )
        frame = args.frame if args.frame >= 0 else n_frames - 1
        myosin = np.asarray(handle["/myosin/center"][frame], dtype=float)
        actin = np.asarray(handle["/actin/center"][frame], dtype=float)
        direction = np.asarray(handle["/actin/direction"][frame], dtype=float)
        step = int(np.asarray(handle["/state/current_step"][frame]).reshape(-1)[-1])

    plus = direction[:, 0] >= 0.0
    minus = ~plus
    populations = (
        (myosin, np.ones(len(myosin), dtype=bool), "Myosin", "viridis", "o", 8),
        (actin, plus, r"Actin $+x$", "Reds", ">", 8),
        (actin, minus, r"Actin $-x$", "Blues", "<", 8),
    )
    norm = Normalize(vmin=-1.0, vmax=1.0)
    fig, axes = plt.subplots(1, 3, figsize=(15.2, 4.7), sharex=True, sharey=True,
                             constrained_layout=True)
    for axis, (coords, mask, title, cmap, marker, size) in zip(axes, populations):
        selected = coords[mask]
        points = axis.scatter(
            selected[:, 0], selected[:, 1], c=selected[:, 2], cmap=cmap, norm=norm,
            marker=marker, s=size, alpha=0.78, linewidths=0, rasterized=True,
        )
        axis.set_title(f"{title}  (n={len(selected)})")
        axis.set_xlabel("x")
        axis.set_xlim(-6.0, 6.0)
        axis.set_ylim(-1.0, 1.0)
        axis.set_aspect("equal", adjustable="box")
        colorbar = fig.colorbar(points, ax=axis, fraction=0.046, pad=0.025)
        colorbar.set_label("z")
    axes[0].set_ylabel("y")

    phi_text = "" if args.phi is None else rf", $\Phi_{{3D}}={args.phi:.3f}$"
    fig.suptitle(f"{args.input.stem}, frame {frame}, step {step}{phi_text}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    main()
