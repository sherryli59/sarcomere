#!/usr/bin/env python3
"""Plot and quantify heartbeat locking of AM heads and ACTN2 link states."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/sarcomere_matplotlib")

import h5py
import matplotlib.pyplot as plt
import numpy as np


def normalized(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    scale = np.std(values)
    return (values - np.mean(values)) / scale if scale > 0 else np.zeros_like(values)


def best_lag(reference: np.ndarray, response: np.ndarray, dt: float,
             period: float) -> tuple[float, float]:
    """Return lag maximizing corr(reference(t), response(t+lag))."""
    ref, rsp = normalized(reference), normalized(response)
    maximum = max(1, int(round(0.5 * period / dt)))
    best = (float("nan"), -np.inf)
    for shift in range(-maximum, maximum + 1):
        if shift < 0:
            x, y = ref[-shift:], rsp[:shift]
        elif shift > 0:
            x, y = ref[:-shift], rsp[shift:]
        else:
            x, y = ref, rsp
        if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
            continue
        correlation = float(np.corrcoef(x, y)[0, 1])
        if correlation > best[1]:
            best = (shift * dt, correlation)
    return best


def fundamental_amplitude(values: np.ndarray, time: np.ndarray,
                          period: float) -> float:
    centered = np.asarray(values, dtype=float) - np.mean(values)
    if np.std(centered) == 0:
        return 0.0
    coefficient = np.mean(centered * np.exp(-2j * np.pi * time / period))
    return float(2.0 * abs(coefficient) / np.std(centered))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trajectory", type=Path)
    parser.add_argument("--period", type=float, required=True)
    parser.add_argument("--plot", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--discard-cycles", type=float, default=1.0)
    args = parser.parse_args()

    with h5py.File(args.trajectory, "r", locking=False) as handle:
        dt_sim = float(handle["meta"].attrs["dt"])
        offsets = np.asarray(handle["catch_bond/local_link_frame_offsets"], dtype=int)
        identity = np.asarray(handle["catch_bond/local_link_identity"], dtype=int)
        steps = offsets[:, 0]
        time = steps * dt_sim
        calcium = np.asarray(handle["actin/am_calcium"], dtype=float).reshape(len(time), -1).mean(1)
        q = np.asarray(handle["actin/am_q"], dtype=float).reshape(len(time), -1).mean(1)
        state = np.asarray(handle["state/am_head_state_current"], dtype=int)
        bound = np.asarray(handle["state/am_head_bound_actin_current"], dtype=int)
        occupied = bound >= 0
        axes = tuple(range(1, state.ndim))
        weak = np.count_nonzero(occupied & (state == 1), axis=axes)
        strong = np.count_nonzero(occupied & (state == 2), axis=axes)
        aa1 = np.zeros(len(time), dtype=int)
        aa2 = np.zeros(len(time), dtype=int)
        for frame, (_, start, count) in enumerate(offsets):
            status = identity[start:start + count, 4]
            aa1[frame] = np.count_nonzero(status == 1)
            aa2[frame] = np.count_nonzero(status == 2)

    keep = time >= args.discard_cycles * args.period
    sample_dt = float(np.median(np.diff(time)))
    metrics = []
    for name, values in (("q", q), ("strong_heads", strong),
                         ("status2_links", aa2)):
        lag, correlation = best_lag(
            calcium[keep], values[keep], sample_dt, args.period)
        metrics.append({
            "signal": name,
            "heartbeat_period_s": args.period,
            "best_lag_s": lag,
            "best_lag_fraction_cycle": lag / args.period,
            "maximum_cross_correlation": correlation,
            "fundamental_amplitude_over_sd": fundamental_amplitude(
                values[keep], time[keep], args.period),
            "mean": float(np.mean(values[keep])),
            "minimum": float(np.min(values[keep])),
            "maximum": float(np.max(values[keep])),
        })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(metrics[0]))
        writer.writeheader()
        writer.writerows(metrics)

    time_ms = 1e3 * time
    fig, axes = plt.subplots(4, 1, figsize=(11, 10), sharex=True)
    for pulse in np.arange(0.0, time[-1] + args.period, args.period):
        for axis in axes:
            axis.axvline(1e3 * pulse, color="#dddddd", lw=0.8, zorder=0)
    axes[0].plot(time_ms, calcium, color="#777777", lw=2, label="calcium")
    axes[0].plot(time_ms, q, color="#111111", lw=2, label="mean thin activation q")
    axes[0].set_ylabel("activation")
    axes[0].legend(frameon=False, ncol=2)
    axes[1].plot(time_ms, weak, "--", color="#6a51a3", lw=1.8, label="weak AM")
    axes[1].plot(time_ms, strong, color="#54278f", lw=2.2, label="strong AM")
    axes[1].set_ylabel("heads")
    axes[1].legend(frameon=False, ncol=2)
    axes[2].plot(time_ms, aa1, "--", color="#e69f00", lw=1.8, label="status 1")
    axes[2].plot(time_ms, aa2, color="#009e73", lw=2.2, label="status 2")
    axes[2].set_ylabel("ACTN2 links")
    axes[2].legend(frameon=False, ncol=2)
    axes[3].plot(time_ms, normalized(calcium), color="#888888", lw=1.6,
                 label="calcium (normalized)")
    axes[3].plot(time_ms, normalized(strong), color="#54278f", lw=1.8,
                 label="strong AM (normalized)")
    axes[3].plot(time_ms, normalized(aa2), color="#009e73", lw=1.8,
                 label="status 2 (normalized)")
    axes[3].set(xlabel="time (ms)", ylabel="centered / SD")
    axes[3].legend(frameon=False, ncol=3)
    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", color="#eeeeee", lw=0.7)
    fig.suptitle(
        f"Heartbeat locking: period = {1e3 * args.period:g} ms; vertical lines mark pulse onset")
    fig.tight_layout()
    args.plot.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.plot, dpi=220, bbox_inches="tight")
    plt.close(fig)

    print(f"wrote {args.output}")
    print(f"wrote {args.plot}")
    for row in metrics:
        print(row)


if __name__ == "__main__":
    main()
