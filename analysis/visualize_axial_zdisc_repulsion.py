#!/usr/bin/env python3
"""Standalone visualization of the proposed axial Z-disc/myosin-tip contact.

The model is deliberately independent of the sarcomere simulator.  It plots
the force/potential law and integrates a one-dimensional overdamped approach
under several constant motor loads.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import numpy as np


def radial_gate(rho: np.ndarray | float, r_flat: float, r_cut: float):
    """Smooth flat-top transverse gate S(rho). Distances are in micrometres."""
    rho_array = np.asarray(rho, dtype=float)
    if not 0.0 <= r_flat < r_cut:
        raise ValueError("Require 0 <= r_flat < r_cut")
    q = np.clip((rho_array - r_flat) / (r_cut - r_flat), 0.0, 1.0)
    gate = 1.0 - 3.0 * q**2 + 2.0 * q**3
    gate = np.where(rho_array <= r_flat, 1.0, gate)
    gate = np.where(rho_array >= r_cut, 0.0, gate)
    return gate


def _log_cosh(value: np.ndarray | float):
    value = np.asarray(value, dtype=float)
    return np.logaddexp(value, -value) - np.log(2.0)


def axial_contact(
    gap: np.ndarray | float,
    rho: np.ndarray | float,
    onset: float,
    stiffness: float,
    force_max: float,
    r_flat: float,
    r_cut: float,
):
    """Return (outward force [pN], potential [pN um])."""
    if stiffness <= 0.0 or force_max <= 0.0:
        raise ValueError("stiffness and force_max must be positive")
    compression = np.maximum(0.0, onset - np.asarray(gap, dtype=float))
    gate = radial_gate(rho, r_flat, r_cut)
    scaled = stiffness * compression / force_max
    force = gate * force_max * np.tanh(scaled)
    potential = gate * force_max**2 / stiffness * _log_cosh(scaled)
    return force, potential


def simulate_approach(
    motor_force: float,
    rho: float,
    onset: float,
    stiffness: float,
    force_max: float,
    r_flat: float,
    r_cut: float,
    mobility: float,
    duration: float,
    dt: float,
    initial_gap: float,
):
    """Overdamped gap dynamics: dg/dt = mobility*(F_repulsion-F_motor)."""
    steps = int(np.ceil(duration / dt)) + 1
    time = np.linspace(0.0, duration, steps)
    gap = np.empty(steps)
    force = np.empty(steps)
    gap[0] = initial_gap
    for index in range(steps - 1):
        force[index] = axial_contact(
            gap[index], rho, onset, stiffness, force_max, r_flat, r_cut
        )[0]
        gap[index + 1] = gap[index] + dt * mobility * (
            force[index] - motor_force
        )
    force[-1] = axial_contact(
        gap[-1], rho, onset, stiffness, force_max, r_flat, r_cut
    )[0]
    return time, gap, force


def build_figure(args: argparse.Namespace):
    import matplotlib.pyplot as plt

    gap = np.linspace(-0.04, 0.12, 600)
    transverse_distances = [0.080, 0.110, 0.125, 0.140]
    colors = plt.cm.viridis(np.linspace(0.12, 0.88, len(transverse_distances)))

    figure, axes = plt.subplots(1, 3, figsize=(14.5, 4.4), constrained_layout=True)
    force_axis, energy_axis, dynamics_axis = axes

    for rho, color in zip(transverse_distances, colors):
        force, potential = axial_contact(
            gap,
            rho,
            args.onset,
            args.stiffness,
            args.force_max,
            args.r_flat,
            args.r_cut,
        )
        label = f"transverse distance = {rho * 1000:.0f} nm"
        force_axis.plot(gap * 1000, force, color=color, label=label)
        energy_axis.plot(gap * 1000, potential / args.kbt, color=color)

    force_axis.axvline(args.onset * 1000, color="0.45", linestyle="--")
    force_axis.axvline(0.0, color="0.45", linestyle=":")
    force_axis.set(
        title="Axial repulsion",
        xlabel="myosin-tip clearance from Z patch (nm)",
        ylabel="outward force (pN)",
    )
    force_axis.legend(frameon=False, fontsize=8)

    energy_axis.axvline(args.onset * 1000, color="0.45", linestyle="--")
    energy_axis.axvline(0.0, color="0.45", linestyle=":")
    energy_axis.set(
        title="Finite-compliance potential",
        xlabel="myosin-tip clearance from Z patch (nm)",
        ylabel=r"potential ($k_\mathrm{B}T$)",
    )

    drive_colors = plt.cm.plasma(np.linspace(0.12, 0.88, len(args.motor_forces)))
    for drive, color in zip(args.motor_forces, drive_colors):
        time, trajectory, _ = simulate_approach(
            drive,
            args.rho,
            args.onset,
            args.stiffness,
            args.force_max,
            args.r_flat,
            args.r_cut,
            args.mobility,
            args.duration,
            args.dt,
            args.initial_gap,
        )
        dynamics_axis.plot(time, trajectory * 1000, color=color, label=f"{drive:g} pN")
    dynamics_axis.axhline(args.onset * 1000, color="0.45", linestyle="--")
    dynamics_axis.axhline(0.0, color="0.45", linestyle=":")
    dynamics_axis.set(
        title=f"Driven approach at {args.rho * 1000:.0f}-nm transverse offset",
        xlabel="time (s)",
        ylabel="myosin-tip clearance (nm)",
    )
    dynamics_axis.legend(title="motor load", frameon=False, fontsize=8)

    for axis in axes:
        axis.grid(alpha=0.18)
    figure.suptitle("Proposed local myosin-tip / ACTN2–Z-disc axial contact")
    return figure


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("axial_zdisc_repulsion.png"))
    parser.add_argument("--show", action="store_true")
    parser.add_argument("--onset", type=float, default=0.050, help="um")
    parser.add_argument("--stiffness", type=float, default=1000.0, help="pN/um")
    parser.add_argument("--force-max", type=float, default=20.0, help="pN")
    parser.add_argument("--r-flat", type=float, default=0.100, help="um")
    parser.add_argument("--r-cut", type=float, default=0.130, help="um")
    parser.add_argument("--rho", type=float, default=0.080, help="um")
    parser.add_argument("--mobility", type=float, default=0.005, help="um/(pN s)")
    parser.add_argument("--duration", type=float, default=3.0, help="s")
    parser.add_argument("--dt", type=float, default=0.001, help="s")
    parser.add_argument("--initial-gap", type=float, default=0.120, help="um")
    parser.add_argument(
        "--motor-forces",
        type=float,
        nargs="+",
        default=[5.0, 10.0, 15.0, 25.0],
        help="pN",
    )
    parser.add_argument("--kbt", type=float, default=0.00428, help="pN um")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.show:
        matplotlib.use("Agg")
    figure = build_figure(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=180, bbox_inches="tight")
    print(args.output.resolve())
    if args.show:
        import matplotlib.pyplot as plt

        plt.show()


if __name__ == "__main__":
    main()
