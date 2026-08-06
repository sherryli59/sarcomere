#!/usr/bin/env python3
"""Plot the polarity structure factor S_p(q) vs mode n for directional vs nondir.

The money figure for "why nondir doesn't assemble": directional puts its polarity
power at a FINITE mode (n~5-6, the sarcomere period) = an alternating antiparallel
register; nondir piles all its power at n=1 (period = whole box) = phase
separation into two big same-polarity domains = a parallel bundle.
"""
from __future__ import annotations

import glob
import os
import re

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SCRATCH = "/scratch/users/xsli/sarcomere"
CAMPAIGNS = {"directional": "data_a2000_status2_lambda",
             "nondir": "data_a2000_status2_nondir_lambda"}
LX, NMAX = 12.0, 15
LAMBDAS = ["1e-3", "1e-2", "1e-1"]   # panels


def stepof(h):
    for k in ("/resume/current_step", "/resume_prev/current_step"):
        if k in h:
            return int(np.asarray(h[k]).reshape(-1)[-1])
    return int(np.asarray(h["/state/current_step"])[-1, 0])


def lt_of(name):
    m = re.search(r"_lt([0-9e.\-]+?)_(?:nodir_)?r\d", name)
    return m.group(1) if m else "?"


def sp_curve(path, nframes=20, minstep=290000):
    import h5py
    qs = 2 * np.pi * np.arange(0, NMAX + 1) / LX
    with h5py.File(path, "r", locking=False) as h:
        if stepof(h) < minstep:
            return None
        C, D = h["/actin/center"], h["/actin/direction"]
        n = min(C.shape[0], D.shape[0])
        acc, cnt = np.zeros(NMAX + 1), 0
        for fr in range(max(0, n - nframes), n):
            x = np.asarray(C[fr])[:, 0]
            p = np.sign(np.asarray(D[fr])[:, 0])
            acc += np.abs(np.exp(1j * np.outer(qs, x)) @ p) ** 2 / len(p)
            cnt += 1
    return acc / cnt if cnt else None


def mean_curve(dirn, lt):
    curves = []
    for p in glob.glob(f"{SCRATCH}/{dirn}/*.h5"):
        if lt_of(os.path.basename(p)) != lt:
            continue
        c = sp_curve(p)
        if c is not None:
            curves.append(c)
    return np.mean(curves, axis=0) if curves else None


def main():
    ns = np.arange(0, NMAX + 1)
    fig, axes = plt.subplots(1, len(LAMBDAS), figsize=(5 * len(LAMBDAS), 4.4),
                             squeeze=False)
    colors = {"directional": "#2a7", "nondir": "#c44"}
    for ax, lt in zip(axes[0], LAMBDAS):
        for tag, dirn in CAMPAIGNS.items():
            c = mean_curve(dirn, lt)
            if c is None:
                continue
            ax.plot(ns, c, "-o", ms=4, color=colors[tag], label=tag)
        ax.axvspan(5, 6, color="#2a7", alpha=0.07)
        ax.annotate("sarcomere\nregister", xy=(5.5, ax.get_ylim()[1]),
                    ha="center", va="top", color="#2a7", fontsize=8)
        ax.annotate("box-spanning\ndomains (bundle)", xy=(1, 0),
                    xytext=(2.5, ax.get_ylim()[1] * 0.55), color="#c44",
                    fontsize=8, arrowprops=dict(arrowstyle="->", color="#c44"))
        ax.set_xlabel("polarity mode  n   (period = Lx/n)")
        ax.set_ylabel(r"$S_p(q_n)$")
        ax.set_title(rf"$\lambda = {float(lt.replace('1e','1e')):g}$"
                     if False else rf"$\lambda$ = {lt}")
        ax.set_xticks(range(0, NMAX + 1, 2))
        ax.grid(alpha=0.25)
        ax.legend(fontsize=9)
    fig.suptitle("Polarity structure factor: directional registers, "
                 "nondirectional bundles", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig("polarity_structure.png", dpi=150)
    print("wrote polarity_structure.png")


if __name__ == "__main__":
    main()
