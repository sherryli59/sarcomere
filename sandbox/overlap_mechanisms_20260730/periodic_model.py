"""Force-balanced periodic two-myosin/four-actin axial sandbox."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Dict, List, Tuple

import numpy as np

from overlap_model import Bond, Parameters


# Four barbed-end coordinates followed by two myosin centers.  Actin polarity
# and owning myosin reproduce C++ two_sarcomere_unit().
DIRECTION = np.array([-1.0, 1.0, -1.0, 1.0])
MYOSIN_OF_ACTIN = np.array([0, 1, 1, 0])
AA_PAIRS = ((0, 1), (2, 3))


@dataclass
class PeriodicState:
    q: np.ndarray
    bonds: List[List[Bond]] = field(default_factory=lambda: [[], []])

    def copy(self) -> "PeriodicState":
        return PeriodicState(
            self.q.copy(),
            [[Bond(b.rest_overlap) for b in group] for group in self.bonds],
        )


def wrap_delta(value: float | np.ndarray, length: float):
    return value - length * np.round(value / length)


def wrap_positions(q: np.ndarray, length: float) -> None:
    q[:] = (q + 0.5 * length) % length - 0.5 * length


def default_periodic_state(
    aa_overlap: float = 0.10,
    box_length: float = 4.30,
) -> PeriodicState:
    ls = 0.5 * box_length
    half = 0.5 * aa_overlap
    # Z-A at 0, Z-B at L/2; values are wrapped into the primary box.
    q = np.array(
        [
            half,
            -half,
            ls + half,
            ls - half,
            -0.5 * ls,
            0.5 * ls,
        ],
        dtype=float,
    )
    wrap_positions(q, box_length)
    return PeriodicState(q)


def aa_overlaps(state: PeriodicState, p: Parameters) -> np.ndarray:
    return np.array(
        [wrap_delta(state.q[i] - state.q[j], p.box_length) for i, j in AA_PAIRS]
    )


def am_overlaps(state: PeriodicState, p: Parameters) -> np.ndarray:
    values = np.empty(4)
    for i in range(4):
        m = state.q[4 + MYOSIN_OF_ACTIN[i]]
        axial_distance = DIRECTION[i] * wrap_delta(m - state.q[i], p.box_length)
        lo = axial_distance - p.myosin_half_length
        hi = axial_distance - p.bare_half_length
        values[i] = max(0.0, min(p.actin_length, hi) - max(0.0, lo))
    return values


def zm_distances(state: PeriodicState, p: Parameters) -> np.ndarray:
    values = np.empty(4)
    for i in range(4):
        m = state.q[4 + MYOSIN_OF_ACTIN[i]]
        values[i] = DIRECTION[i] * wrap_delta(m - state.q[i], p.box_length)
    return values


def initialize_bonds(state: PeriodicState, p: Parameters) -> None:
    overlaps = aa_overlaps(state, p)
    for pair in range(2):
        state.bonds[pair] = [
            Bond(float(overlaps[pair]))
            for _ in range(min(p.initial_bonds, p.n_actn2_sites))
        ]


def update_bonds(state: PeriodicState, p: Parameters, rng: np.random.Generator) -> None:
    overlaps = aa_overlaps(state, p)
    for pair in range(2):
        y = float(overlaps[pair])
        retained = []
        for bond in state.bonds[pair]:
            if p.permanent_bonds:
                retained.append(bond)
                continue
            load = abs(p.actn2_k * (bond.rest_overlap - y))
            off_fraction = 1.0
            if p.actn2_catch_force > 0:
                off_fraction = max(
                    p.actn2_min_off_fraction,
                    np.exp(-load / p.actn2_catch_force),
                )
            p_break = 1.0 - np.exp(-p.actn2_koff * off_fraction * p.dt)
            if rng.random() >= p_break:
                retained.append(bond)
        state.bonds[pair] = retained
        available = p.n_actn2_sites - len(retained)
        fraction = np.clip(y / p.aa_zone_width, 0.0, 1.0)
        p_bind = 1.0 - np.exp(-p.actn2_kon * fraction * p.dt)
        for _ in range(available):
            if rng.random() < p_bind:
                state.bonds[pair].append(Bond(y))


def force_components(
    state: PeriodicState,
    p: Parameters,
    motor_on: bool = True,
) -> Dict[str, np.ndarray]:
    components: Dict[str, np.ndarray] = {}
    motor = np.zeros(6)
    if motor_on:
        overlaps = am_overlaps(state, p)
        for i, overlap in enumerate(overlaps):
            magnitude = min(p.motor_force_max, p.motor_density * overlap)
            force = DIRECTION[i] * magnitude
            motor[i] += force
            motor[4 + MYOSIN_OF_ACTIN[i]] -= force
    components["motor"] = motor

    elastic = np.zeros(6)
    adhesion = np.zeros(6)
    entropic = np.zeros(6)
    register = np.zeros(6)
    overlaps = aa_overlaps(state, p)
    for pair_index, (i, j) in enumerate(AA_PAIRS):
        y = overlaps[pair_index]
        f_elastic = sum(p.actn2_k * (bond.rest_overlap - y) for bond in state.bonds[pair_index])
        elastic[i] += f_elastic
        elastic[j] -= f_elastic
        if p.adhesion_gamma > 0 and 0 < y < p.aa_zone_width:
            adhesion[i] += p.adhesion_gamma
            adhesion[j] -= p.adhesion_gamma
        if p.entropic_nkbt > 0 and y > 0:
            f = min(p.entropic_force_cap, p.entropic_nkbt * p.kbt / max(y, 1e-4))
            entropic[i] += f
            entropic[j] -= f
        if p.register_k > 0:
            f = p.register_k * (p.aa_target - y)
            register[i] += f
            register[j] -= f
    components["actn2_elastic"] = elastic
    components["actn2_adhesion"] = adhesion
    components["actn2_entropic"] = entropic
    components["zdisc_register"] = register

    titin = np.zeros(6)
    if p.titin_k > 0:
        for i, distance in enumerate(zm_distances(state, p)):
            force = p.titin_k * (distance - p.titin_rest) * DIRECTION[i]
            titin[i] += force
            titin[4 + MYOSIN_OF_ACTIN[i]] -= force
    components["titin"] = titin
    return components


def simulate_periodic(
    p: Parameters,
    seed: int,
    initial: PeriodicState | None = None,
) -> Dict[str, np.ndarray]:
    rng = np.random.default_rng(seed)
    state = default_periodic_state(p.aa_target, p.box_length) if initial is None else initial.copy()
    if not any(state.bonds):
        initialize_bonds(state, p)
    nsteps = int(round(p.duration / p.dt))
    nsave = nsteps // p.save_every + 1
    time = np.empty(nsave)
    coords = np.empty((nsave, 6))
    overlaps = np.empty((nsave, 6))
    bonds = np.empty((nsave, 2), dtype=int)
    mobility = np.array([1 / p.zeta_actin] * 4 + [1 / p.zeta_myosin] * 2)
    noise = p.thermal_scale * np.sqrt(2 * p.kbt * p.dt * mobility)

    def save(slot: int, step: int) -> None:
        time[slot] = step * p.dt
        coords[slot] = state.q
        overlaps[slot, :2] = aa_overlaps(state, p)
        overlaps[slot, 2:] = am_overlaps(state, p)
        bonds[slot] = [len(group) for group in state.bonds]

    save(0, 0)
    slot = 1
    for step in range(1, nsteps + 1):
        update_bonds(state, p, rng)
        total = sum(force_components(state, p).values(), np.zeros(6))
        state.q += mobility * total * p.dt + noise * rng.normal(size=6)
        wrap_positions(state.q, p.box_length)
        if step % p.save_every == 0:
            save(slot, step)
            slot += 1
    return {"time": time, "coords": coords, "overlaps": overlaps, "bond_count": bonds}


def periodic_parameters(name: str, box_length: float = 4.30) -> Parameters:
    base = Parameters()
    values = {"box_length": box_length, "dt": 5e-4, "duration": 4.0, "save_every": 5}
    if name == "baseline_turnover":
        return replace(base, **values)
    if name == "permanent_actn2":
        return replace(base, permanent_bonds=True, **values)
    if name == "johann_adhesion":
        return replace(base, adhesion_gamma=5.0, **values)
    if name == "johann_entropic":
        return replace(base, entropic_nkbt=100.0, entropic_force_cap=8.0, **values)
    if name == "titin":
        return replace(base, titin_k=220.0, **values)
    if name == "entropic_plus_titin":
        return replace(base, entropic_nkbt=100.0, entropic_force_cap=8.0, titin_k=220.0, **values)
    if name == "register_plus_titin":
        return replace(base, register_k=180.0, titin_k=220.0, **values)
    raise ValueError(name)

