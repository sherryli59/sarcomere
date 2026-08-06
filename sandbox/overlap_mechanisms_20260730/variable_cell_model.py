"""Variable-length extension of the exact periodic axial sandbox.

The fixed-cell model can rearrange particles inside a repeat but cannot select
the mean Z-disc--M-line spacing: that spacing is constrained by ``Lx``.  This
module adds one slow, overdamped repeat-length coordinate.  Particle motion is
integrated as before; the cell then deforms affinely according to the internal
axial virial, an optional applied tensile load, and an optional series spring.

This is a diagnostic barostat-like model, not production simulator code.  Its
purpose is to ask whether a proposed interaction supplies a restoring force in
the otherwise missing global-length direction.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Dict

import numpy as np

from overlap_model import Parameters
from periodic_model import (
    AA_PAIRS,
    DIRECTION,
    MYOSIN_OF_ACTIN,
    PeriodicState,
    aa_overlaps,
    am_overlaps,
    default_periodic_state,
    force_components,
    initialize_bonds,
    update_bonds,
    wrap_delta,
    wrap_positions,
    zm_distances,
)


@dataclass
class CellParameters:
    mechanics: Parameters
    cell_mobility: float = 0.02  # um / (pN s)
    applied_tension: float = 0.0  # pN; positive expands the repeat
    series_k: float = 0.0  # pN / um
    series_rest_length: float = 4.30  # um
    cell_thermal_scale: float = 0.0
    min_length: float = 3.20
    max_length: float = 5.20


@dataclass
class VariableCellState:
    particles: PeriodicState
    length: float

    def copy(self) -> "VariableCellState":
        return VariableCellState(self.particles.copy(), float(self.length))


def default_variable_state(
    aa_overlap: float = 0.10,
    length: float = 4.30,
) -> VariableCellState:
    return VariableCellState(default_periodic_state(aa_overlap, length), length)


def _pair_cell_force(
    force_on_first: float,
    first: float,
    second: float,
    length: float,
) -> float:
    """Generalized force conjugate to length for one reciprocal pair."""

    return force_on_first * wrap_delta(first - second, length) / length


def internal_cell_force_components(
    state: VariableCellState,
    cp: CellParameters,
    motor_on: bool = True,
) -> Dict[str, float]:
    """Return internal generalized forces conjugate to repeat length.

    Pairwise minimum-image displacements are used so the result is independent
    of where a cross-boundary pair happens to be wrapped.
    """

    q = state.particles.q
    length = state.length
    p = replace(cp.mechanics, box_length=length)
    result: Dict[str, float] = {
        "motor": 0.0,
        "actn2_elastic": 0.0,
        "actn2_adhesion": 0.0,
        "actn2_entropic": 0.0,
        "zdisc_register": 0.0,
        "titin": 0.0,
    }

    if motor_on:
        for i, overlap in enumerate(am_overlaps(state.particles, p)):
            magnitude = min(p.motor_force_max, p.motor_density * overlap)
            force = DIRECTION[i] * magnitude
            m = 4 + MYOSIN_OF_ACTIN[i]
            result["motor"] += _pair_cell_force(force, q[i], q[m], length)

    overlaps = aa_overlaps(state.particles, p)
    for pair_index, (i, j) in enumerate(AA_PAIRS):
        y = float(overlaps[pair_index])
        elastic = sum(
            p.actn2_k * (bond.rest_overlap - y)
            for bond in state.particles.bonds[pair_index]
        )
        result["actn2_elastic"] += _pair_cell_force(
            elastic, q[i], q[j], length
        )
        if p.adhesion_gamma > 0 and 0 < y < p.aa_zone_width:
            result["actn2_adhesion"] += _pair_cell_force(
                p.adhesion_gamma, q[i], q[j], length
            )
        if p.entropic_nkbt > 0 and y > 0:
            entropic = min(
                p.entropic_force_cap,
                p.entropic_nkbt * p.kbt / max(y, 1e-4),
            )
            result["actn2_entropic"] += _pair_cell_force(
                entropic, q[i], q[j], length
            )
        if p.register_k > 0:
            register = p.register_k * (p.aa_target - y)
            result["zdisc_register"] += _pair_cell_force(
                register, q[i], q[j], length
            )

    if p.titin_k > 0:
        for i, distance in enumerate(zm_distances(state.particles, p)):
            force = p.titin_k * (distance - p.titin_rest) * DIRECTION[i]
            m = 4 + MYOSIN_OF_ACTIN[i]
            result["titin"] += _pair_cell_force(force, q[i], q[m], length)
    return result


def external_cell_force(state: VariableCellState, cp: CellParameters) -> float:
    return (
        cp.applied_tension
        - cp.series_k * (state.length - cp.series_rest_length)
    )


def change_cell_length(state: VariableCellState, new_length: float) -> None:
    """Affinely change the repeat length while preserving fractional positions."""

    old_length = state.length
    if old_length <= 0 or new_length <= 0:
        raise ValueError("cell lengths must be positive")
    state.particles.q *= new_length / old_length
    state.length = float(new_length)
    wrap_positions(state.particles.q, state.length)


def simulate_variable_cell(
    cp: CellParameters,
    seed: int,
    initial: VariableCellState | None = None,
    motor_on: bool = True,
) -> Dict[str, np.ndarray]:
    p0 = cp.mechanics
    state = (
        default_variable_state(p0.aa_target, p0.box_length)
        if initial is None
        else initial.copy()
    )
    if not any(state.particles.bonds):
        initialize_bonds(state.particles, replace(p0, box_length=state.length))

    rng = np.random.default_rng(seed)
    nsteps = int(round(p0.duration / p0.dt))
    nsave = nsteps // p0.save_every + 1
    time = np.empty(nsave)
    coords = np.empty((nsave, 6))
    overlaps = np.empty((nsave, 6))
    length = np.empty(nsave)
    zm = np.empty((nsave, 4))
    bond_count = np.empty((nsave, 2), dtype=int)
    cell_force_names = list(internal_cell_force_components(state, cp, motor_on))
    cell_forces = {name: np.empty(nsave) for name in cell_force_names}
    cell_forces["external"] = np.empty(nsave)

    mobility = np.array([1 / p0.zeta_actin] * 4 + [1 / p0.zeta_myosin] * 2)
    particle_noise = (
        p0.thermal_scale
        * np.sqrt(2 * p0.kbt * p0.dt * mobility)
    )
    cell_noise = (
        cp.cell_thermal_scale
        * np.sqrt(2 * p0.kbt * cp.cell_mobility * p0.dt)
    )

    def current_parameters() -> Parameters:
        return replace(p0, box_length=state.length)

    def save(slot: int, step: int) -> None:
        p = current_parameters()
        time[slot] = step * p0.dt
        coords[slot] = state.particles.q
        overlaps[slot, :2] = aa_overlaps(state.particles, p)
        overlaps[slot, 2:] = am_overlaps(state.particles, p)
        length[slot] = state.length
        zm[slot] = zm_distances(state.particles, p)
        bond_count[slot] = [len(group) for group in state.particles.bonds]
        internal = internal_cell_force_components(state, cp, motor_on)
        for name, value in internal.items():
            cell_forces[name][slot] = value
        cell_forces["external"][slot] = external_cell_force(state, cp)

    save(0, 0)
    slot = 1
    for step in range(1, nsteps + 1):
        p = current_parameters()
        update_bonds(state.particles, p, rng)
        particle_force = sum(
            force_components(state.particles, p, motor_on).values(),
            np.zeros(6),
        )
        state.particles.q += mobility * particle_force * p0.dt
        state.particles.q += particle_noise * rng.normal(size=6)
        wrap_positions(state.particles.q, state.length)

        internal = internal_cell_force_components(state, cp, motor_on)
        total_cell_force = sum(internal.values()) + external_cell_force(state, cp)
        new_length = (
            state.length
            + cp.cell_mobility * total_cell_force * p0.dt
            + cell_noise * rng.normal()
        )
        new_length = float(np.clip(new_length, cp.min_length, cp.max_length))
        change_cell_length(state, new_length)

        if step % p0.save_every == 0:
            save(slot, step)
            slot += 1

    result: Dict[str, np.ndarray] = {
        "time": time,
        "coords": coords,
        "overlaps": overlaps,
        "length": length,
        "zm_distances": zm,
        "bond_count": bond_count,
    }
    for name, values in cell_forces.items():
        result[f"cell_force_{name}"] = values
    return result


def variable_cell_parameters(name: str) -> CellParameters:
    base = Parameters(
        dt=1e-3,
        duration=12.0,
        save_every=5,
        box_length=4.30,
    )
    if name == "baseline_turnover":
        return CellParameters(base)
    if name == "permanent_actn2":
        return CellParameters(replace(base, permanent_bonds=True))
    if name == "johann_adhesion":
        return CellParameters(replace(base, adhesion_gamma=5.0))
    if name == "johann_entropic":
        return CellParameters(
            replace(base, entropic_nkbt=100.0, entropic_force_cap=8.0)
        )
    if name == "titin":
        return CellParameters(replace(base, titin_k=220.0))
    if name == "entropic_plus_titin":
        return CellParameters(
            replace(
                base,
                entropic_nkbt=100.0,
                entropic_force_cap=8.0,
                titin_k=220.0,
            )
        )
    if name == "register_plus_titin":
        return CellParameters(
            replace(base, register_k=180.0, titin_k=220.0)
        )
    if name == "entropic_plus_series":
        return CellParameters(
            replace(
                base,
                duration=15.0,
                entropic_nkbt=100.0,
                entropic_force_cap=8.0,
            ),
            cell_mobility=0.05,
            series_k=100.0,
            series_rest_length=4.30,
        )
    raise ValueError(name)
