#!/usr/bin/env python3
"""Plot axial density kymographs for representative corrected-sweep runs."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d


DATA = Path("/scratch/users/xsli/sarcomere/data_nalam_s2fstab_base1e-3")
LAMBDAS = (("0", 0.0), ("1e-3", 1.0), ("1e-2.5", 10**0.5), ("1e-1", 100.0))
ACTINS = (900, 1800)
N_MYOSINS = 1000
BOX_X = 12.0
DT = 1.0e-5


def densities(path: Path, stride: int, bins: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(path, "r", locking=False) as handle:
        n = min(handle["/myosin/center"].shape[0], handle["/actin/center"].shape[0])
        idx = np.arange(0, n, stride)
        steps = np.asarray(handle["/state/current_step"][idx]).reshape(-1)
        mx = np.asarray(handle["/myosin/center"][idx, :, 0]) % BOX_X
        ax = np.asarray(handle["/actin/center"][idx, :, 0]) % BOX_X
        status = np.asarray(handle["/actin/cb_status"][idx]).reshape(len(idx), ax.shape[1])
    edges = np.linspace(0.0, BOX_X, bins + 1)
    myosin = np.asarray([np.histogram(row, edges)[0] for row in mx], dtype=float)
    cb2 = np.asarray([np.histogram(row[mask], edges)[0]
                      for row, mask in zip(ax, status == 2)], dtype=float)
    # Divide by a uniform full-occupancy expectation. Myosin values therefore
    # have mean one; CB2 values have a row mean equal to the status-2 fraction.
    myosin /= mx.shape[1] / bins
    cb2 /= ax.shape[1] / bins
    myosin = gaussian_filter1d(myosin, 1.0, axis=1, mode="wrap")
    cb2 = gaussian_filter1d(cb2, 1.0, axis=1, mode="wrap")
    return steps * DT, myosin, cb2


def plot(dataset: dict[tuple[int, str], tuple[np.ndarray, np.ndarray, np.ndarray]],
         output: Path, which: str) -> None:
    fig, axes = plt.subplots(2, 4, figsize=(14, 6.3), sharex=True, sharey=True,
                             constrained_layout=True)
    vmax = 3.0 if which == "myosin" else 1.8
    for row, n_actins in enumerate(ACTINS):
        for col, (label, ratio) in enumerate(LAMBDAS):
            time, myosin, cb2 = dataset[(n_actins, label)]
            values = myosin if which == "myosin" else cb2
            image = axes[row, col].imshow(
                values, origin="lower", aspect="auto", interpolation="nearest",
                extent=(0.0, BOX_X, time[0], time[-1]), vmin=0.0, vmax=vmax,
                cmap="magma",
            )
            axes[row, col].set_title(rf"$\lambda/\tau_1={ratio:g}$")
            if col == 0:
                axes[row, col].set_ylabel(f"A={n_actins}\ntime (s)")
            if row == 1:
                axes[row, col].set_xlabel("axial position x")
    label = "myosin density / uniform density" if which == "myosin" else "status-2 ACTN2 density / full-occupancy density"
    fig.colorbar(image, ax=axes, label=label, shrink=0.88)
    fig.suptitle("Axial assembly pathway: bands emerge, merge, or become trapped", fontsize=14)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--stride", type=int, default=5)
    parser.add_argument("--bins", type=int, default=60)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    dataset = {}
    for n_actins in ACTINS:
        for label, _ in LAMBDAS:
            filename = f"A{n_actins}_M{N_MYOSINS}_s2fstab_lt{label}_r00_L12x2x2.h5"
            dataset[(n_actins, label)] = densities(DATA / filename, args.stride, args.bins)
    plot(dataset, args.output_dir / "myosin_density_kymographs.png", "myosin")
    plot(dataset, args.output_dir / "cb2_density_kymographs.png", "cb2")


if __name__ == "__main__":
    main()
