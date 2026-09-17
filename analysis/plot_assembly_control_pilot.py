#!/usr/bin/env python3
"""Compact response surface for the turnover/activation assembly pilots.

One panel per outcome family, each on its own scale.  The two pilot axes
(actin turnover at full motor availability, and motor availability with
turnover off) share the `turn0` arm as their common control, so they are drawn
as two labelled bands separated by a rule rather than as one continuous sweep.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Categorical slots 1-2 of the validated default palette.  Checked with
# scripts/validate_palette.js --mode light --pairs all: all six checks PASS
# (CVD dE 24.7, normal-vision dE 33.6, both >= 3:1 on the light surface).
SERIES = ("#2a78d6", "#eb6834")
INK = "#0b0b0b"
INK_2 = "#52514e"
SURFACE = "#fcfcfb"
# Reserved status colour, used only for the placement-failure annotation.
CRITICAL = "#e34948"

TURNOVER_ARMS = ["turn0", "turn0p1", "turn0p35", "turn1"]
MOTOR_ARMS = ["healthy", "intermediate", "hcm"]
ORDER = TURNOVER_ARMS + MOTOR_ARMS
LABEL = {
    "turn0": "0\n(control)", "turn0p1": "0.1", "turn0p35": "0.35", "turn1": "1.0",
    "healthy": "healthy\n0.10", "intermediate": "intermed.\n0.30", "hcm": "HCM\n0.50",
}

PANELS = [
    (("tail_zdisc_capture_mean", "mean over Z-discs"),
     ("tail_zdisc_capture_max", "best Z-disc"),
     "Per-Z-disc capture\nlargest component / actins in that Z-disc", (0.0, 1.0)),
    (("tail_lcc_span_y", "y extent"),
     ("tail_lcc_span_z", "z extent"),
     "Transverse reach of largest component\nfraction of box spanned", (0.0, 1.0)),
    (("tail_m_order", "M order"),
     ("tail_z_order", "Z order"),
     "Axial order\nbest harmonic 1-4, tail mean", (0.0, None)),
    (("tail_on_fraction", "ON heads"),
     (None, None),
     "Motor availability\ndetached-ON + bound head fraction", (0.0, 1.0)),
    (("tail_duplicate_fraction", "duplicate ACTN2"),
     (None, None),
     "Link redundancy\nmolecules on an already-linked pair", (0.0, 1.0)),
    (("successful_turnovers", "renewals placed"),
     (None, None),
     "Actin renewal\nevents placed over the measured window", (0.0, None)),
]


def dictated_on_fraction(rows, arm):
    """Zero-tension equilibrium ON fraction actually configured for this arm."""
    for row in rows:
        if row["arm"] != arm:
            continue
        feedback = str(row.get("thick_feedback", "")).strip().lower()
        if feedback in ("0", "false", ""):
            return 1.0
        try:
            recruit = float(row["thick_recruit_rate"])
            fold = float(row["thick_fold_rate"])
        except (KeyError, ValueError, TypeError):
            return None
        return recruit / (recruit + fold) if recruit + fold > 0 else None
    return None


def collect(rows, arm, metric):
    out = []
    for row in rows:
        if row["arm"] != arm:
            continue
        value = row.get(metric, "")
        if value in ("", None):
            continue
        try:
            out.append(float(value))
        except ValueError:
            continue
    return np.array(out, dtype=float)


def draw_panel(axis, rows, arms, series, title, limits):
    metric_first = series[0][0]
    (metric_a, label_a), (metric_b, label_b) = series
    entries = [(metric_a, label_a, SERIES[0])]
    if metric_b is not None:
        entries.append((metric_b, label_b, SERIES[1]))
    span = 0.13 if len(entries) > 1 else 0.0
    seen = []
    for slot, (metric, label, color) in enumerate(entries):
        shift = (slot - (len(entries) - 1) / 2) * span
        for index, arm in enumerate(arms):
            values = collect(rows, arm, metric)
            if not len(values):
                continue
            seen.append(values)
            jitter = np.linspace(-0.04, 0.04, len(values)) if len(values) > 1 else np.zeros(1)
            axis.scatter(index + shift + jitter, values, s=18, color=color,
                         alpha=0.55, linewidths=0, zorder=2)
            axis.plot([index + shift], [values.mean()], "o", ms=9, color=color,
                      markeredgecolor=SURFACE, markeredgewidth=2.0, zorder=3,
                      label=label if index == 0 else None)
    axis.set_title(title, fontsize=9.5, color=INK, loc="left", pad=7)
    axis.set_xticks(np.arange(len(arms)))
    axis.set_xticklabels([LABEL.get(a, a) for a in arms], fontsize=8, color=INK_2)
    axis.set_xlim(-0.6, len(arms) - 0.4)
    low, high = limits
    if seen:
        flat = np.concatenate(seen)
        data_low = float(np.min(flat))
        data_high = float(np.max(flat))
        pad = 0.08 * max(data_high - data_low, 1e-9)
        if low is None:
            low = data_low - pad
        if high is None:
            high = data_high + pad if data_high > 0 else 1.0
        if high <= low:
            high = low + 1.0
    axis.set_ylim(low, high)
    axis.tick_params(axis="y", labelsize=8, colors=INK_2, length=0)
    axis.grid(axis="y", color=INK_2, alpha=0.14, linewidth=0.7)
    axis.set_axisbelow(True)
    for side in ("top", "right", "left"):
        axis.spines[side].set_visible(False)
    axis.spines["bottom"].set_color(INK_2)
    axis.spines["bottom"].set_alpha(0.35)
    if len(entries) > 1:
        axis.legend(frameon=False, fontsize=8, loc="upper left",
                    labelcolor=INK_2, handletextpad=0.35, borderaxespad=0.2)
    if metric_first == "tail_on_fraction":
        # Dictated zero-tension equilibrium, so the panel reads as a
        # calibration check rather than a bare measurement.
        for index, arm in enumerate(arms):
            target = dictated_on_fraction(rows, arm)
            if target is None:
                continue
            axis.plot([index - 0.22, index + 0.22], [target, target],
                      color=INK_2, alpha=0.55, linewidth=1.4,
                      linestyle=(0, (3, 2)), zorder=1,
                      label="dictated p_ON,0" if index == 0 else None)
        axis.legend(frameon=False, fontsize=8, loc="center left",
                    labelcolor=INK_2, handletextpad=0.35, borderaxespad=0.2)
    if metric_first == "successful_turnovers":
        # Placement failure is the silent-null mode for this axis: a renewal
        # that exhausts its attempts does nothing.  Call it out when nonzero.
        for index, arm in enumerate(arms):
            failed = collect(rows, arm, "failed_turnovers")
            if not len(failed) or failed.max() <= 0:
                continue
            placed = collect(rows, arm, "successful_turnovers")
            denom = failed.mean() + (placed.mean() if len(placed) else 0.0)
            share = failed.mean() / denom if denom > 0 else 0.0
            top = placed.mean() if len(placed) else 0.0
            axis.annotate("%.1f failed (%.1f%%)" % (failed.mean(), 100 * share),
                          xy=(index, top), xytext=(0, 11),
                          textcoords="offset points", ha="center",
                          fontsize=7.5, color=CRITICAL)

    # Separate the two pilot axes; turn0 is the shared control.
    n_turn = sum(1 for a in arms if a in TURNOVER_ARMS)
    if 0 < n_turn < len(arms):
        axis.axvline(n_turn - 0.5, color=INK_2, alpha=0.28, linewidth=0.9,
                     linestyle=(0, (4, 3)), zorder=1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("summary", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--title", default="A800/M160 assembly-control pilot")
    args = parser.parse_args()

    with args.summary.open() as stream:
        rows = list(csv.DictReader(stream))
    arms = [arm for arm in ORDER if any(row["arm"] == arm for row in rows)]
    if not arms:
        raise SystemExit("no recognised arms in %s" % args.summary)

    fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.0), constrained_layout=True)
    fig.patch.set_facecolor(SURFACE)
    for axis in axes.ravel():
        axis.set_facecolor(SURFACE)
    for axis, (a, b, title, limits) in zip(axes.ravel(), PANELS):
        draw_panel(axis, rows, arms, (a, b), title, limits)

    seeds = len({row.get("label", row["path"]) for row in rows}) // max(len(arms), 1)
    fig.suptitle(args.title, fontsize=12.5, color=INK, x=0.006, ha="left")
    n_turn = sum(1 for a in arms if a in TURNOVER_ARMS)
    caption = (
        "left of the rule: actin turnover k_turn (s^-1) at full motor availability   |   "
        "right: basal ON fraction p_ON,0 with turnover off   |   "
        "turn0 is the shared control.  "
        "Small marks are individual seeds (n=%d); ringed marks are arm means."
        % max(seeds, 1)
    ) if 0 < n_turn < len(arms) else (
        "Small marks are individual seeds; ringed marks are arm means.  "
        "Global lcc is capped near 1/n_Zdiscs, so per-Z-disc capture is used instead."
    )
    fig.text(0.006, -0.012, caption, fontsize=8, color=INK_2, ha="left")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=220, facecolor=SURFACE, bbox_inches="tight")
    svg = args.output.with_suffix(".svg")
    fig.savefig(svg, facecolor=SURFACE, bbox_inches="tight")
    print(args.output)
    print(svg)


if __name__ == "__main__":
    main()
