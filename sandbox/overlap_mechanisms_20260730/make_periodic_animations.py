#!/usr/bin/env python3
"""Render exact-topology mechanism animations."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

from periodic_model import DIRECTION, default_periodic_state, periodic_parameters, simulate_periodic


OUT = Path(__file__).resolve().parent / "results"


def wrapped_segments(points: np.ndarray, length: float):
    wrapped = (points + 0.5 * length) % length - 0.5 * length
    cuts = np.where(np.abs(np.diff(wrapped)) > 0.5 * length)[0] + 1
    return np.split(wrapped, cuts)


def render(name: str) -> None:
    p = periodic_parameters(name)
    result = simulate_periodic(p, 0, default_periodic_state(0.10, p.box_length))
    indices = np.linspace(0, len(result["time"]) - 1, 180, dtype=int)
    fig, ax = plt.subplots(figsize=(10, 3.8))

    def draw(frame):
        ax.clear()
        idx = indices[frame]
        q = result["coords"][idx]
        for m, y in zip(q[4:], (-0.03, 0.03)):
            points = m + np.linspace(-p.myosin_half_length, p.myosin_half_length, 80)
            for segment in wrapped_segments(points, p.box_length):
                ax.plot(segment, np.full_like(segment, y), lw=10, color="#1f77b4")
        colors = ("#d62728", "#ff7f0e", "#d62728", "#ff7f0e")
        levels = (0.18, -0.18, 0.30, -0.30)
        for i in range(4):
            points = q[i] + DIRECTION[i] * np.linspace(0, p.actin_length, 80)
            for segment in wrapped_segments(points, p.box_length):
                ax.plot(segment, np.full_like(segment, levels[i]), lw=6, color=colors[i])
            ax.scatter([q[i]], [levels[i]], color="black", s=25, zorder=5)
        y_aa = np.mean(result["overlaps"][idx, :2])
        y_am = np.mean(result["overlaps"][idx, 2:])
        ax.set_xlim(-0.5 * p.box_length, 0.5 * p.box_length)
        ax.set_ylim(-0.48, 0.48)
        ax.set_yticks([])
        ax.set_xlabel("periodic axial coordinate (µm)")
        ax.set_title(
            f"{name}   t={result['time'][idx]:.2f} s   "
            f"mean A–A={y_aa:.3f} µm   mean A–M={y_am:.3f} µm"
        )
        return ()

    animation = FuncAnimation(fig, draw, frames=len(indices), interval=40, blit=False)
    animation.save(OUT / f"periodic_{name}.gif", writer=PillowWriter(fps=20))
    plt.close(fig)


if __name__ == "__main__":
    for mechanism in ("baseline_turnover", "entropic_plus_titin", "register_plus_titin"):
        render(mechanism)

