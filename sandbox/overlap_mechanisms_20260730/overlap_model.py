"""Reduced stochastic sandbox for sarcomere overlap mechanisms.

This is deliberately separate from the production C++ simulator.  It retains
the four axial degrees of freedom needed to distinguish A--A overlap from A--M
overlap:

    a_left, a_right, m_left, m_right

The two ``a`` coordinates are the barbed ends meeting at a central Z-disc.
The two ``m`` coordinates are neighboring bipolar myosin centers.  Motor,
ACTN2, titin-like, and M-line forces are applied reciprocally, and all four
bodies undergo overdamped Brownian motion.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Dict, Iterable, List, Tuple

import numpy as np


@dataclass
class Parameters:
    box_length: float = 4.30
    dt: float = 2.0e-4
    duration: float = 4.0
    save_every: int = 10
    kbt: float = 0.00428  # pN um at 310 K
    zeta_actin: float = 20.0  # pN s / um
    zeta_myosin: float = 40.0
    thermal_scale: float = 1.0

    actin_length: float = 1.0
    myosin_half_length: float = 0.75
    bare_half_length: float = 0.085

    motor_density: float = 8.0  # pN / um of A-M overlap
    motor_force_max: float = 6.0  # pN per half-sarcomere

    aa_zone_width: float = 0.10
    aa_target: float = 0.10
    n_actn2_sites: int = 6
    actn2_k: float = 120.0  # pN / um per link
    actn2_kon: float = 12.0  # 1/s at full zone overlap
    actn2_koff: float = 8.0  # 1/s
    actn2_catch_force: float = 2.0  # pN; <=0 disables catch stabilization
    actn2_min_off_fraction: float = 0.05
    initial_bonds: int = 4
    permanent_bonds: bool = False

    # Population-level overlap free energies.  gamma has units pN because
    # U_adhesion = -gamma * overlap.  entropic_nkbt is the dimensionless N in
    # F = N kBT / overlap.
    adhesion_gamma: float = 0.0
    entropic_nkbt: float = 0.0
    entropic_force_cap: float = 10.0

    # A persistent Z-disc lattice/register interaction.
    register_k: float = 0.0

    # Titin-like springs connect each barbed-end/Z-disc coordinate to its
    # neighboring myosin center.
    titin_k: float = 0.0
    titin_rest: float = 1.125

    # Optional external M-line positional stiffness.  This is a boundary
    # control, not part of the proposed minimal myosin+ACTN2 system.
    m_anchor_k: float = 0.0
    m_left_anchor: float = -1.075
    m_right_anchor: float = 1.075


@dataclass
class Bond:
    rest_overlap: float


@dataclass
class State:
    q: np.ndarray
    bonds: List[Bond] = field(default_factory=list)

    def copy(self) -> "State":
        return State(self.q.copy(), [Bond(b.rest_overlap) for b in self.bonds])


def default_state(aa_overlap: float = 0.10, zm_distance: float = 1.125) -> State:
    """Symmetric two-half-sarcomere state.

    a_left - a_right is the A--A overlap.  ``zm_distance`` is the distance
    from each myosin center to the neighboring barbed-end coordinate.
    """

    a_left = 0.5 * aa_overlap
    a_right = -0.5 * aa_overlap
    m_left = a_left - zm_distance
    m_right = a_right + zm_distance
    return State(np.array([a_left, a_right, m_left, m_right], dtype=float))


def aa_overlap(q: np.ndarray) -> float:
    return float(q[0] - q[1])


def zm_distances(q: np.ndarray) -> Tuple[float, float]:
    return float(q[0] - q[2]), float(q[3] - q[1])


def interval_overlap(lo_a: float, hi_a: float, lo_b: float, hi_b: float) -> float:
    return max(0.0, min(hi_a, hi_b) - max(lo_a, lo_b))


def am_overlaps(q: np.ndarray, p: Parameters) -> Tuple[float, float]:
    a_left, a_right, m_left, m_right = q
    left = interval_overlap(
        a_left - p.actin_length,
        a_left,
        m_left + p.bare_half_length,
        m_left + p.myosin_half_length,
    )
    right = interval_overlap(
        a_right,
        a_right + p.actin_length,
        m_right - p.myosin_half_length,
        m_right - p.bare_half_length,
    )
    return left, right


def _motor_magnitudes(q: np.ndarray, p: Parameters) -> Tuple[float, float]:
    x_left, x_right = am_overlaps(q, p)
    return (
        min(p.motor_force_max, p.motor_density * x_left),
        min(p.motor_force_max, p.motor_density * x_right),
    )


def force_components(state: State, p: Parameters, motor_on: bool = True) -> Dict[str, np.ndarray]:
    """Return named forces on [a_left, a_right, m_left, m_right]."""

    q = state.q
    y = aa_overlap(q)
    components: Dict[str, np.ndarray] = {}

    motor = np.zeros(4)
    if motor_on:
        f_left, f_right = _motor_magnitudes(q, p)
        motor += np.array([-f_left, f_right, f_left, -f_right])
    components["motor"] = motor

    elastic = np.zeros(4)
    for bond in state.bonds:
        f = p.actn2_k * (bond.rest_overlap - y)
        elastic += np.array([f, -f, 0.0, 0.0])
    components["actn2_elastic"] = elastic

    adhesion = np.zeros(4)
    if p.adhesion_gamma > 0.0 and 0.0 < y < p.aa_zone_width:
        adhesion += np.array([p.adhesion_gamma, -p.adhesion_gamma, 0.0, 0.0])
    components["actn2_adhesion"] = adhesion

    entropic = np.zeros(4)
    if p.entropic_nkbt > 0.0 and y > 0.0:
        f = min(p.entropic_force_cap, p.entropic_nkbt * p.kbt / max(y, 1.0e-4))
        entropic += np.array([f, -f, 0.0, 0.0])
    components["actn2_entropic"] = entropic

    register = np.zeros(4)
    if p.register_k > 0.0:
        f = p.register_k * (p.aa_target - y)
        register += np.array([f, -f, 0.0, 0.0])
    components["zdisc_register"] = register

    titin = np.zeros(4)
    if p.titin_k > 0.0:
        d_left, d_right = zm_distances(q)
        # U = 1/2 k (d-d0)^2 for each half-sarcomere.
        f_a_left = -p.titin_k * (d_left - p.titin_rest)
        f_a_right = p.titin_k * (d_right - p.titin_rest)
        titin += np.array([f_a_left, f_a_right, -f_a_left, -f_a_right])
    components["titin"] = titin

    anchor = np.zeros(4)
    if p.m_anchor_k > 0.0:
        anchor[2] = -p.m_anchor_k * (q[2] - p.m_left_anchor)
        anchor[3] = -p.m_anchor_k * (q[3] - p.m_right_anchor)
    components["mline_anchor"] = anchor
    return components


def total_force(state: State, p: Parameters, motor_on: bool = True) -> np.ndarray:
    return sum(force_components(state, p, motor_on).values(), np.zeros(4))


def update_bonds(state: State, p: Parameters, rng: np.random.Generator) -> None:
    y = aa_overlap(state.q)
    retained: List[Bond] = []
    for bond in state.bonds:
        if p.permanent_bonds:
            retained.append(bond)
            continue
        load = abs(p.actn2_k * (bond.rest_overlap - y))
        if p.actn2_catch_force > 0.0:
            off_fraction = max(
                p.actn2_min_off_fraction,
                np.exp(-load / p.actn2_catch_force),
            )
        else:
            off_fraction = 1.0
        if rng.random() >= 1.0 - np.exp(-p.actn2_koff * off_fraction * p.dt):
            retained.append(bond)
    state.bonds = retained

    available = max(0, p.n_actn2_sites - len(state.bonds))
    overlap_fraction = np.clip(y / p.aa_zone_width, 0.0, 1.0)
    bind_probability = 1.0 - np.exp(-p.actn2_kon * overlap_fraction * p.dt)
    for _ in range(available):
        if rng.random() < bind_probability:
            # Turnover-capable ACTN2 anneals by adopting the current registry.
            state.bonds.append(Bond(rest_overlap=y))


def initialize_bonds(state: State, p: Parameters) -> None:
    state.bonds = [
        Bond(rest_overlap=aa_overlap(state.q))
        for _ in range(min(p.initial_bonds, p.n_actn2_sites))
    ]


def simulate(
    p: Parameters,
    seed: int,
    initial: State | None = None,
    motor_on: bool = True,
) -> Dict[str, np.ndarray]:
    rng = np.random.default_rng(seed)
    state = default_state() if initial is None else initial.copy()
    if not state.bonds:
        initialize_bonds(state, p)

    nsteps = int(round(p.duration / p.dt))
    nsave = nsteps // p.save_every + 1
    times = np.empty(nsave)
    coords = np.empty((nsave, 4))
    overlaps = np.empty((nsave, 3))
    bond_counts = np.empty(nsave, dtype=int)
    force_names = list(force_components(state, p, motor_on))
    recorded_forces = {name: np.empty((nsave, 4)) for name in force_names}

    mobility = np.array(
        [1.0 / p.zeta_actin, 1.0 / p.zeta_actin,
         1.0 / p.zeta_myosin, 1.0 / p.zeta_myosin]
    )
    noise_std = p.thermal_scale * np.sqrt(2.0 * p.kbt * p.dt * mobility)

    def save(slot: int, step: int) -> None:
        times[slot] = step * p.dt
        coords[slot] = state.q
        x_left, x_right = am_overlaps(state.q, p)
        overlaps[slot] = [aa_overlap(state.q), x_left, x_right]
        bond_counts[slot] = len(state.bonds)
        comps = force_components(state, p, motor_on)
        for name in force_names:
            recorded_forces[name][slot] = comps[name]

    save(0, 0)
    slot = 1
    for step in range(1, nsteps + 1):
        update_bonds(state, p, rng)
        force = total_force(state, p, motor_on)
        state.q += mobility * force * p.dt
        state.q += noise_std * rng.normal(size=4)
        if step % p.save_every == 0:
            save(slot, step)
            slot += 1

    result: Dict[str, np.ndarray] = {
        "time": times,
        "coords": coords,
        "overlaps": overlaps,
        "bond_count": bond_counts,
    }
    for name, values in recorded_forces.items():
        result[f"force_{name}"] = values
    return result


def mechanism_parameters(name: str) -> Parameters:
    base = Parameters()
    if name == "baseline_turnover":
        return base
    if name == "permanent_actn2":
        return replace(base, permanent_bonds=True)
    if name == "johann_adhesion":
        return replace(base, adhesion_gamma=5.0)
    if name == "johann_entropic":
        return replace(base, entropic_nkbt=100.0, entropic_force_cap=8.0)
    if name == "zdisc_register":
        return replace(base, register_k=180.0)
    if name == "titin":
        return replace(base, titin_k=220.0)
    if name == "adhesion_plus_titin":
        return replace(base, adhesion_gamma=5.0, titin_k=220.0)
    if name == "entropic_plus_titin":
        return replace(
            base,
            entropic_nkbt=100.0,
            entropic_force_cap=8.0,
            titin_k=220.0,
        )
    if name == "register_plus_titin":
        return replace(base, register_k=180.0, titin_k=220.0)
    if name == "fixed_mline_adhesion":
        return replace(base, adhesion_gamma=5.0, m_anchor_k=1000.0)
    raise ValueError(f"unknown mechanism: {name}")


MECHANISMS = (
    "baseline_turnover",
    "permanent_actn2",
    "johann_adhesion",
    "johann_entropic",
    "zdisc_register",
    "titin",
    "adhesion_plus_titin",
    "entropic_plus_titin",
    "register_plus_titin",
    "fixed_mline_adhesion",
)
