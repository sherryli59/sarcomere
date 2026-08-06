#!/usr/bin/env python3
"""Animate the diagnostic dynamic-repeat trajectories."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

from periodic_model import DIRECTION
from variable_cell_model import (
    default_variable_state,
    simulate_variable_cell,
    variable_cell_parameters,
)


OUT = Path(__file__).resolve().parent / "results"


def wrapped_segments(points: np.ndarray, length: float):
    wrapped = (points + 0.5 * length) % length - 0.5 * length
    cuts = np.where(np.abs(np.diff(wrapped)) > 0.5 * length)[0] + 1
    return np.split(wrapped, cuts)


def render(name: str, initial_length: float) -> None:
    cp = variable_cell_parameters(name)
    p = cp.mechanics
    result = simulate_variable_cell(
        cp,
        0,
        default_variable_state(0.10, initial_length),
    )
    indices = np.linspace(0, len(result["time"]) - 1, 180, dtype=int)
    fig, ax = plt.subplots(figsize=(10, 4.2))

    def draw(frame):
        ax.clear()
        idx = indices[frame]
        q = result["coords"][idx]
        length = result["length"][idx]
        for m, y in zip(q[4:], (-0.03, 0.03)):
            points = m + np.linspace(
                -p.myosin_half_length, p.myosin_half_length, 80
            )
            for segment in wrapped_segments(points, length):
                ax.plot(
                    segment,
                    np.full_like(segment, y),
                    lw=10,
                    color="#1f77b4",
                )
        colors = ("#d62728", "#ff7f0e", "#d62728", "#ff7f0e")
        levels = (0.18, -0.18, 0.30, -0.30)
        for i in range(4):
            points = q[i] + DIRECTION[i] * np.linspace(
                0, p.actin_length, 80
            )
            for segment in wrapped_segments(points, length):
                ax.plot(
                    segment,
                    np.full_like(segment, levels[i]),
                    lw=6,
                    color=colors[i],
                )
            ax.scatter(
                [q[i]], [levels[i]], color="black", s=25, zorder=5
            )
        ax.axvline(-0.5 * length, color="0.35", ls="--", lw=1.5)
        ax.axvline(0.5 * length, color="0.35", ls="--", lw=1.5)
        y_aa = np.mean(result["overlaps"][idx, :2])
        y_am = np.mean(result["overlaps"][idx, 2:])
        ax.set_xlim(-2.65, 2.65)
        ax.set_ylim(-0.48, 0.48)
        ax.set_yticks([])
        ax.set_xlabel("axial coordinate (µm); dashed lines are cell edges")
        ax.set_title(
            f"{name}   t={result['time'][idx]:.1f} s   L={length:.3f} µm\n"
            f"mean A–A={y_aa:.3f} µm   mean A–M={y_am:.3f} µm"
        )
        return ()

    animation = FuncAnimation(
        fig, draw, frames=len(indices), interval=40, blit=False
    )
    animation.save(
        OUT / f"variable_{name}_from_{initial_length:.1f}.gif",
        writer=PillowWriter(fps=20),
    )
    plt.close(fig)


if __name__ == "__main__":
    render("register_plus_titin", 4.70)
    render("entropic_plus_series", 4.70)

