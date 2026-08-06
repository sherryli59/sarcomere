#!/usr/bin/env python3
"""Plot instantaneous restoring fields in (A--A overlap, repeat length)."""

from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from periodic_model import force_components, initialize_bonds
from variable_cell_model import (
    default_variable_state,
    external_cell_force,
    internal_cell_force_components,
    variable_cell_parameters,
)


OUT = Path(__file__).resolve().parent / "results"
MECHANISMS = (
    "johann_entropic",
    "entropic_plus_titin",
    "register_plus_titin",
    "entropic_plus_series",
)


def rates(name: str, y: float, length: float):
    cp = variable_cell_parameters(name)
    state = default_variable_state(y, length)
    p = replace(cp.mechanics, box_length=length)
    initialize_bonds(state.particles, p)
    particle_force = sum(
        force_components(state.particles, p).values(), np.zeros(6)
    )
    dy_particles = (
        particle_force[0] - particle_force[1]
    ) / p.zeta_actin
    q_cell = (
        sum(internal_cell_force_components(state, cp).values())
        + external_cell_force(state, cp)
    )
    dlength = cp.cell_mobility * q_cell
    # Affine deformation also changes the physical A--A separation.
    dy = dy_particles + y * dlength / length
    return dy, dlength


def main() -> None:
    y_values = np.linspace(0.035, 0.21, 32)
    length_values = np.linspace(3.75, 4.75, 34)
    Y, L = np.meshgrid(y_values, length_values)
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True, sharey=True)
    for ax, name in zip(axes.flat, MECHANISMS):
        DY = np.empty_like(Y)
        DL = np.empty_like(L)
        for row in range(Y.shape[0]):
            for col in range(Y.shape[1]):
                DY[row, col], DL[row, col] = rates(
                    name, Y[row, col], L[row, col]
                )
        # Normalize arrow length while retaining speed as color.
        scaled_dy = DY / 0.20
        scaled_dl = DL / 1.0
        speed = np.hypot(scaled_dy, scaled_dl)
        denom = np.maximum(speed, 1e-12)
        ax.quiver(
            Y,
            L,
            scaled_dy / denom,
            scaled_dl / denom,
            color="#24788f",
            pivot="mid",
            scale=38,
        )
        ax.contour(Y, L, DY, levels=[0], colors="#d62728", linewidths=2)
        ax.contour(
            Y, L, DL, levels=[0], colors="#1f77b4", linewidths=2
        )
        ax.scatter([0.10], [4.30], marker="+", s=100, color="black")
        ax.set_title(name.replace("_", " "))
        ax.set_xlabel("A–A overlap (µm)")
        ax.set_ylabel("repeat length (µm)")
        ax.text(
            0.038,
            4.70,
            "red: d(A–A)/dt = 0\nblue: dL/dt = 0",
            fontsize=8,
            va="top",
        )
    fig.suptitle(
        "Two-coordinate force field (stress-free ACTN2 rebinding state)"
    )
    fig.subplots_adjust(
        left=0.08, right=0.97, bottom=0.08, top=0.92, wspace=0.18, hspace=0.22
    )
    fig.savefig(OUT / "variable_cell_phase_portraits.png", dpi=190)


if __name__ == "__main__":
    main()
