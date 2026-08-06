#!/usr/bin/env python3
"""Plot legacy myosin phi_3d vs constant tau1, one panel per lambda.

Reads phi3d_tau1_answer.csv (written by phi3d_tau1_answer.py). Each panel shows
the tau1 response curve against two horizontal references:
  floor   = s2only  (tau1 -> 0, no status-1 stabilization)
  ceiling = s1s2    (old load-DEPENDENT stabilization, lambda*f_load)

tau1 = 0 is drawn as a square at the far left of the log axis (a true zero has
no place on a log scale); the vertical guide marks tau1 = lambda, the
"matched-to-old" point where load-independent tau1 should meet the ceiling.
"""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# where tau1=0 is drawn on the log axis (a decade below the smallest real rung)
ZERO_X = 10 ** -5.5


def load(path):
    by_lt = defaultdict(list)
    floor, ceil = {}, {}
    with open(path) as f:
        for r in csv.DictReader(f):
            lt = r["lt"]
            by_lt[lt].append((float(r["tau1"]), float(r["phi_tau1"]),
                              int(r["n_rep"])))
            if r["phi_s2only_floor"]:
                floor[lt] = float(r["phi_s2only_floor"])
            if r["phi_s1s2_ceiling"]:
                ceil[lt] = float(r["phi_s1s2_ceiling"])
    return by_lt, floor, ceil


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--csv", type=Path, default=Path("phi3d_tau1_answer.csv"))
    p.add_argument("--out", type=Path, default=Path("phi3d_tau1.png"))
    args = p.parse_args()

    by_lt, floor, ceil = load(args.csv)
    lam = {"1e-3": 1e-3, "1e-2": 1e-2}
    lts = sorted(by_lt, key=lambda s: lam.get(s, 0.0))

    fig, axes = plt.subplots(1, len(lts), figsize=(5.4 * len(lts), 4.6),
                             squeeze=False)
    for ax, lt in zip(axes[0], lts):
        pts = sorted(by_lt[lt])
        xs = [ZERO_X if t == 0 else t for t, _, _ in pts]
        ys = [y for _, y, _ in pts]
        n1 = [n == 1 for _, _, n in pts]
        lo, hi = floor.get(lt), ceil.get(lt)

        if lo is not None:
            ax.axhline(lo, ls="--", lw=1.4, color="#c44",
                       label=f"s2only floor ({lo:.2f})")
        if hi is not None:
            ax.axhline(hi, ls="--", lw=1.4, color="#2a7",
                       label=f"s1s2 ceiling ({hi:.2f})")
            ax.axhspan(lo, hi, color="#2a7", alpha=0.05)
        # tau1 = lambda guide
        ax.axvline(lam[lt], ls=":", lw=1.2, color="#666")
        ax.annotate(r"$\tau_1=\lambda$", xy=(lam[lt], ax.get_ylim()[0]),
                    xytext=(lam[lt] * 1.15, lo - 0.06 if lo else 0.35),
                    color="#666", fontsize=9)

        ax.plot(xs, ys, "-o", color="#345", zorder=3)
        # ring the single-replica points so noisy cells are visible as such
        for x, y, one in zip(xs, ys, n1):
            if one:
                ax.plot([x], [y], "o", mfc="none", mec="#c60", mew=1.8,
                        ms=12, zorder=4)

        ax.set_xscale("log")
        ax.set_xlabel(r"status-1 lifetime  $\tau_1$")
        ax.set_ylabel(r"legacy myosin  $\phi_{3D}$")
        ax.set_title(rf"$\lambda = {lam[lt]:g}$")
        # relabel the leftmost synthetic tick as 0
        ticks = [ZERO_X, 1e-4, 1e-3, 1e-2, 1e-1]
        ticks = [t for t in ticks if t <= max(xs) * 3]
        ax.set_xticks(ticks)
        ax.set_xticklabels(["0"] + [f"{t:g}" for t in ticks[1:]])
        ax.grid(alpha=0.25)
        ax.legend(loc="lower right", fontsize=8, framealpha=0.9)

    orange = plt.Line2D([], [], marker="o", mfc="none", mec="#c60", mew=1.8,
                        ls="none", ms=10, label="single replica (noisy)")
    fig.legend(handles=[orange], loc="upper right", fontsize=8,
               framealpha=0.9)
    fig.suptitle(r"Constant load-independent $\tau_1$ recovers the s1s2 "
                 r"$\phi_{3D}$ boost", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(args.out, dpi=150)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
