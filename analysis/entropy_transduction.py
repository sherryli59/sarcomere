"""Entropy production and active-work estimates for sarcomere HDF5 trajectories.

The trajectory stores saved-frame forces, active velocities, catch-bond status,
and load. This module intentionally opens required datasets by explicit path so
partially damaged or append-interrupted optional HDF5 groups do not prevent the
main analysis from running.
"""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np

try:
    import torch
except ImportError:
    torch = None


DEFAULT_BASE_DIRS = (
    Path("/scratch/users/xsli/sarcomere/data0416"),
    Path("/scratch/users/xsli/sarcomere/data0430"),
)

DEFAULT_CHEM_ENTROPY_DATASET = "/actin/cumulative_chem_entropy"
DEFAULT_CHEM_INPUT_DATASET = "/actin/cumulative_chem_input"
DEFAULT_BINDING_EVENTS_DATASET = "/actin/cumulative_chem_binding_events"
DEFAULT_UNBINDING_EVENTS_DATASET = "/actin/cumulative_chem_unbinding_events"

ENTROPY_CSV_FIELDS = [
    "Filename",
    "H5_Path",
    "A_val",
    "lt_str",
    "N_Frames",
    "Saved_dt",
    "Heat_Dissipated",
    "Active_Work",
    "Elastic_Storage",
    "Actin_Mech_Entropy",
    "Myosin_Mech_Entropy",
    "Mech_Entropy",
    "Chem_Entropy",
    "Actin_Total_Entropy",
    "Myosin_Total_Entropy",
    "Total_Entropy",
    "Actin_Entropy_Rate",
    "Myosin_Entropy_Rate",
    "Entropy_Rate",
    "Active_Work_Rate",
    "Mech_Rate",
    "Chem_Rate",
    "Total_Rate",
    "Actin_Total_Rate",
    "Myosin_Total_Rate",
    "Binding_Events",
    "Unbinding_Events",
    "IFT_Segment_Mean",
    "IFT_Segment_SEM",
]


@dataclass(frozen=True)
class EntropyTransductionResult:
    path: Path
    n_frames: int
    saved_dt: float
    beta: float
    mechanical_entropy: float
    actin_mechanical_entropy: float
    myosin_mechanical_entropy: float
    heat_dissipated: float
    active_work: float
    elastic_storage: float
    actin_active_work: float
    myosin_active_work: float
    actin_heat_dissipated: float
    myosin_heat_dissipated: float
    actin_elastic_storage: float
    myosin_elastic_storage: float
    first_law_residual: float
    first_law_relative_residual: float
    chemical_entropy: float
    chemical_input_kbt: float
    chemical_input_energy: float
    actin_total_entropy: float
    myosin_total_entropy: float
    total_entropy: float
    actin_entropy_rate: float
    myosin_entropy_rate: float
    entropy_rate: float
    active_work_rate: float
    efficiency: float
    binding_events: int
    unbinding_events: int
    visible_transition_note: str
    ift_segment_mean: float
    ift_segment_sem: float

    def as_dict(self) -> dict[str, float | int | str | Path]:
        return self.__dict__.copy()


@dataclass(frozen=True)
class EntropyRateTimeseries:
    time: np.ndarray
    mechanical_rate: np.ndarray
    actin_mechanical_rate: np.ndarray
    myosin_mechanical_rate: np.ndarray
    chemical_rate: np.ndarray
    actin_chemical_rate: np.ndarray
    total_rate: np.ndarray
    actin_total_rate: np.ndarray
    myosin_total_rate: np.ndarray
    heat_rate: np.ndarray
    actin_heat_rate: np.ndarray
    myosin_heat_rate: np.ndarray
    active_work_rate: np.ndarray
    actin_active_work_rate: np.ndarray
    myosin_active_work_rate: np.ndarray
    elastic_storage_rate: np.ndarray
    actin_elastic_storage_rate: np.ndarray
    myosin_elastic_storage_rate: np.ndarray


def _require(fh: h5py.File, path: str) -> h5py.Dataset:
    if path not in fh:
        raise KeyError(f"Required dataset {path!r} is missing")
    return fh[path]


def _has_dataset(fh: h5py.File, path: str) -> bool:
    try:
        obj = fh[path]
    except Exception:
        return False
    return isinstance(obj, h5py.Dataset)


def _read_cumulative_frame_sums(
    fh: h5py.File,
    path: str,
    n_frames: int,
) -> np.ndarray | None:
    if not _has_dataset(fh, path):
        return None
    data = np.asarray(fh[path][:n_frames], dtype=float)
    if data.shape[0] < n_frames:
        raise ValueError(f"Cumulative dataset {path!r} has fewer frames than required")
    if data.ndim == 1:
        return data
    return np.sum(data, axis=tuple(range(1, data.ndim)))


def _cumulative_deltas(
    cumulative: np.ndarray,
    frame_slice: slice,
    previous_value: float,
) -> tuple[np.ndarray, float]:
    current = cumulative[frame_slice]
    if current.size == 0:
        return np.empty(0, dtype=float), previous_value
    with_previous = np.concatenate(([previous_value], current))
    return np.diff(with_previous), float(current[-1])


def _infer_saved_dt(
    fh: h5py.File,
    simulation_dt: float,
    save_every: int | None,
) -> tuple[float, int | None]:
    if _has_dataset(fh, "/state/current_step"):
        steps = np.asarray(fh["/state/current_step"][:, 0], dtype=float)
        diffs = np.diff(steps)
        diffs = diffs[np.isfinite(diffs) & (diffs > 0)]
        if diffs.size:
            inferred_save_every = int(round(float(np.median(diffs))))
            return simulation_dt * inferred_save_every, inferred_save_every
    if save_every is None:
        raise ValueError("save_every is required when /state/current_step is unavailable")
    return simulation_dt * save_every, save_every


def _iter_frame_chunks(n_frames: int, chunk_frames: int) -> Iterable[slice]:
    for start in range(1, n_frames, chunk_frames):
        yield slice(start, min(start + chunk_frames, n_frames))


def resolve_compute_device(device: str = "auto"):
    if device == "numpy":
        return None
    if torch is None:
        if device == "auto":
            return None
        raise RuntimeError("PyTorch is not installed; use --device numpy or install torch.")
    if device == "auto":
        return torch.device("cuda") if torch.cuda.is_available() else None
    if device == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA was requested, but torch.cuda.is_available() is False.")
    return torch.device(device)


def print_compute_device_info(torch_device, requested_device: str = "auto") -> None:
    print(f"Requested compute backend: {requested_device}")
    if torch is None:
        print("PyTorch: not installed")
    else:
        print(f"PyTorch: {torch.__version__}; CUDA available: {torch.cuda.is_available()}")

    if torch_device is None:
        print("Using numpy for entropy computation.")
        return

    print(f"Using torch/{torch_device} for entropy computation.")
    if torch_device.type == "cuda":
        idx = torch_device.index if torch_device.index is not None else torch.cuda.current_device()
        props = torch.cuda.get_device_properties(idx)
        total_gb = props.total_memory / 1024**3
        capability = f"{props.major}.{props.minor}"
        print(
            f"GPU {idx}: {props.name}; "
            f"compute capability {capability}; total memory {total_gb:.2f} GiB"
        )


def _torch_float_dtype(torch_device):
    return torch.float32 if torch_device is not None and torch_device.type == "cuda" else torch.float64


def _drift_thermo_chunk(
    active_velocity: np.ndarray,
    *,
    beta: float,
    diffusion_coeff: float,
    saved_dt: float,
    torch_device=None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return per-frame heat (multiply by beta at call site to get entropy).

    Mechanical entropy from motor sliding (Sekimoto, velocity-only):

        ΔS/kB = β × γ|V|² × saved_dt = |V|²/D × saved_dt

    where V = active_velocity = the motor-prescribed velocity (actin.velocity
    or myosin.velocity in the HDF5), γ = 1/(βD) is friction, D is the
    translational diffusion coefficient.

    Why velocity-only (not including spring force F):
      Springs are conservative. Their contribution β*(F+γV)·dx = β*F·dx + β*γV·dx
      splits into (i) β*F·dx ≈ -β dU, which is an exact differential averaging
      to zero in steady state, and (ii) β*γV·dx, the motor friction term captured
      here. Including F via the trajectory estimator V·Δr/D gives negative results
      in practice because the spring force RESISTS motor sliding (V anti-correlated
      with F), and the V·F cross-term dominates before it can converge. The motor
      velocity V already encodes spring loading through the f_load factor,
      V = (1 - f_load)*v_am, so heavily loaded motors automatically contribute less.

    Returns γ|V|² × saved_dt so that entropy = beta × return_value.
    """
    if diffusion_coeff <= 0:
        raise ValueError("diffusion coefficients must be positive")

    mobility = beta * diffusion_coeff   # μ = β D
    gamma = 1.0 / mobility             # γ = 1/(β D)

    if torch_device is None:
        active_power = gamma * np.einsum(
            "fij,fij->f", active_velocity, active_velocity, optimize=True
        )
    else:
        v_tensor = torch.as_tensor(
            active_velocity,
            dtype=_torch_float_dtype(torch_device),
            device=torch_device,
        )
        with torch.no_grad():
            active_power = (
                gamma * torch.sum(v_tensor.square(), dim=(1, 2))
            ).detach().cpu().numpy()

    heat_rate = active_power
    storage_rate = np.zeros_like(active_power)
    return active_power * saved_dt, heat_rate * saved_dt, storage_rate * saved_dt


def _transition_entropy_chunk_numpy(
    status: np.ndarray,
    load: np.ndarray,
    prev_status: np.ndarray,
    *,
    k_on: float,
    base_lifetime: float,
    lifetime_coeff: float,
    directional_base_lifetime: float,
    simulation_dt: float,
    eps: float,
) -> tuple[np.ndarray, np.ndarray, int, int]:
    prev = np.empty_like(status)
    prev[0] = prev_status
    if status.shape[0] > 1:
        prev[1:] = status[:-1]

    current_load = np.maximum(load, 0.0)
    formed = (prev == 0) & (status == 2)
    broken = (prev == 2) & (status == 0)
    lifetime = base_lifetime + directional_base_lifetime + lifetime_coeff * current_load
    zero_load = current_load == 0.0
    k_off = np.where(zero_load, 1.0 / simulation_dt, 1.0 / np.maximum(lifetime, eps))
    affinity = np.log(np.maximum(k_on, eps) / np.maximum(k_off, eps))

    formed_affinity = np.where(formed, affinity, 0.0)
    broken_affinity = np.where(broken, affinity, 0.0)
    transition_entropy = np.sum(formed_affinity, axis=1) - np.sum(broken_affinity, axis=1)
    chemical_input = np.sum(np.abs(formed_affinity), axis=1) + np.sum(
        np.abs(broken_affinity), axis=1
    )
    return (
        transition_entropy,
        chemical_input,
        int(np.count_nonzero(formed)),
        int(np.count_nonzero(broken)),
    )


def _transition_entropy_chunk_torch(
    status: np.ndarray,
    load: np.ndarray,
    prev_status: np.ndarray,
    *,
    k_on: float,
    base_lifetime: float,
    lifetime_coeff: float,
    directional_base_lifetime: float,
    simulation_dt: float,
    eps: float,
    torch_device,
) -> tuple[np.ndarray, np.ndarray, int, int]:
    dtype = _torch_float_dtype(torch_device)
    status_t = torch.as_tensor(status, dtype=torch.int8, device=torch_device)
    load_t = torch.as_tensor(load, dtype=dtype, device=torch_device)
    prev_status_t = torch.as_tensor(prev_status, dtype=torch.int8, device=torch_device)

    with torch.no_grad():
        prev = torch.empty_like(status_t)
        prev[0] = prev_status_t
        if status_t.shape[0] > 1:
            prev[1:] = status_t[:-1]

        current_load = torch.clamp(load_t, min=0.0)
        formed = (prev == 0) & (status_t == 2)
        broken = (prev == 2) & (status_t == 0)
        lifetime = base_lifetime + directional_base_lifetime + lifetime_coeff * current_load
        zero_load = current_load == 0.0
        k_off = torch.where(
            zero_load,
            torch.full_like(current_load, 1.0 / simulation_dt),
            1.0 / torch.clamp(lifetime, min=eps),
        )
        affinity = torch.log(
            torch.as_tensor(max(k_on, eps), dtype=dtype, device=torch_device)
            / torch.clamp(k_off, min=eps)
        )

        formed_affinity = torch.where(formed, affinity, torch.zeros((), dtype=dtype, device=torch_device))
        broken_affinity = torch.where(broken, affinity, torch.zeros((), dtype=dtype, device=torch_device))
        transition_entropy = torch.sum(formed_affinity, dim=1) - torch.sum(broken_affinity, dim=1)
        chemical_input = torch.sum(torch.abs(formed_affinity), dim=1) + torch.sum(
            torch.abs(broken_affinity), dim=1
        )
        binding_events = int(torch.count_nonzero(formed).detach().cpu().item())
        unbinding_events = int(torch.count_nonzero(broken).detach().cpu().item())
        transition_entropy_np = transition_entropy.detach().cpu().numpy()
        chemical_input_np = chemical_input.detach().cpu().numpy()

    return transition_entropy_np, chemical_input_np, binding_events, unbinding_events


def _transition_entropy_chunk(
    status: np.ndarray,
    load: np.ndarray,
    prev_status: np.ndarray,
    *,
    k_on: float,
    base_lifetime: float,
    lifetime_coeff: float,
    directional_base_lifetime: float,
    simulation_dt: float,
    eps: float,
    torch_device=None,
) -> tuple[np.ndarray, np.ndarray, int, int]:
    if torch_device is None:
        return _transition_entropy_chunk_numpy(
            status,
            load,
            prev_status,
            k_on=k_on,
            base_lifetime=base_lifetime,
            lifetime_coeff=lifetime_coeff,
            directional_base_lifetime=directional_base_lifetime,
            simulation_dt=simulation_dt,
            eps=eps,
        )
    return _transition_entropy_chunk_torch(
        status,
        load,
        prev_status,
        k_on=k_on,
        base_lifetime=base_lifetime,
        lifetime_coeff=lifetime_coeff,
        directional_base_lifetime=directional_base_lifetime,
        simulation_dt=simulation_dt,
        eps=eps,
        torch_device=torch_device,
    )


def compute_entropy_transduction(
    h5_path: str | Path,
    *,
    beta: float = 241.0,
    simulation_dt: float = 1e-5,
    actin_diff_coeff_trans: float = 0.2,
    myosin_diff_coeff_trans: float = 0.05,
    save_every: int | None = None,
    k_on: float = 1000.0,
    base_lifetime: float = 0.0,
    lifetime_coeff: float = 10**-2.5,
    directional_base_lifetime: float = 0.0,
    load_dataset: str = "/actin/f_load_cb",
    chem_entropy_dataset: str = DEFAULT_CHEM_ENTROPY_DATASET,
    chem_input_dataset: str = DEFAULT_CHEM_INPUT_DATASET,
    binding_events_dataset: str = DEFAULT_BINDING_EVENTS_DATASET,
    unbinding_events_dataset: str = DEFAULT_UNBINDING_EVENTS_DATASET,
    include_actin_velocity: bool = True,
    chunk_frames: int = 128,
    segment_frames: int = 100,
    normalization_count: float = 1.0,
    device: str = "auto",
) -> EntropyTransductionResult:
    """Compute saved-frame entropy-production and active-work estimates.

    Heavy frame/particle reductions are vectorized. HDF5 datasets are still
    read in chunks, so CUDA memory use is controlled by ``chunk_frames``.
    With ``device="auto"``, torch/CUDA is used when available; otherwise this
    function falls back to vectorized NumPy.

    Mechanical entropy uses the motor-friction estimator:
      ΔS/kB = β × γ|V|² × saved_dt = |V|²/D × saved_dt
    where V = active motor velocity (load-adjusted via f_load). Spring forces
    are conservative and average to zero at steady state; the motor velocity
    encodes spring loading through V = (1 - f_load)*v_am. Pass
    ``include_actin_velocity=False`` to zero out actin mechanical entropy,
    leaving actin contribution as chemical-only. Chemical entropy is read from
    engine-side cumulative KMC
    event entropy when available.
    Older trajectories fall back to visible actin cb_status transitions between
    saved frames: 0 -> 2 contributes +log(k_on / k_off(load)), and 2 -> 0
    contributes the negative of that affinity. Thus:

      actin total entropy  = actin mechanical entropy + chemical entropy
      myosin total entropy = myosin mechanical entropy
    """
    path = Path(h5_path)
    if chunk_frames < 1:
        raise ValueError("chunk_frames must be positive")
    if segment_frames < 1:
        raise ValueError("segment_frames must be positive")
    if normalization_count <= 0:
        raise ValueError("normalization_count must be positive")
    torch_device = resolve_compute_device(device)

    with h5py.File(path, "r") as fh:
        actin_force = _require(fh, "/actin/force")
        actin_velocity = _require(fh, "/actin/velocity")
        myosin_force = _require(fh, "/myosin/force")
        myosin_velocity = _require(fh, "/myosin/velocity")
        cb_status = _require(fh, "/actin/cb_status")
        load_path = load_dataset if _has_dataset(fh, load_dataset) else "/actin/f_load"
        loads = _require(fh, load_path)

        n_frames = int(cb_status.shape[0])
        saved_dt, inferred_save_every = _infer_saved_dt(fh, simulation_dt, save_every)
        cumulative_chem_entropy = _read_cumulative_frame_sums(
            fh, chem_entropy_dataset, n_frames
        )
        cumulative_chem_input = _read_cumulative_frame_sums(fh, chem_input_dataset, n_frames)
        cumulative_binding_events = _read_cumulative_frame_sums(
            fh, binding_events_dataset, n_frames
        )
        cumulative_unbinding_events = _read_cumulative_frame_sums(
            fh, unbinding_events_dataset, n_frames
        )
        use_cumulative_chem = cumulative_chem_entropy is not None

        if min(
            actin_force.shape[0],
            actin_velocity.shape[0],
            myosin_force.shape[0],
            myosin_velocity.shape[0],
        ) < n_frames:
            raise ValueError("Force/velocity datasets have fewer frames than /actin/cb_status")

        actin_active_work = 0.0
        myosin_active_work = 0.0
        actin_heat_dissipated = 0.0
        myosin_heat_dissipated = 0.0
        actin_elastic_storage = 0.0
        myosin_elastic_storage = 0.0
        chemical_entropy = 0.0
        chemical_input_kbt = 0.0
        binding_events = 0
        unbinding_events = 0
        frame_entropy_chunks: list[np.ndarray] = []

        prev_status = np.asarray(cb_status[0, :, 0], dtype=np.int8)
        prev_cumulative_chem = (
            float(cumulative_chem_entropy[0]) if cumulative_chem_entropy is not None else 0.0
        )
        prev_cumulative_input = (
            float(cumulative_chem_input[0]) if cumulative_chem_input is not None else 0.0
        )
        eps = np.finfo(float).tiny

        for frame_slice in _iter_frame_chunks(n_frames, chunk_frames):
            mv = np.asarray(myosin_velocity[frame_slice], dtype=float)

            if include_actin_velocity:
                av = np.asarray(actin_velocity[frame_slice], dtype=float)
                actin_active, actin_heat, actin_storage = _drift_thermo_chunk(
                    av,
                    beta=beta,
                    diffusion_coeff=actin_diff_coeff_trans,
                    saved_dt=saved_dt,
                    torch_device=torch_device,
                )
            else:
                frame_count = frame_slice.stop - frame_slice.start
                actin_active = np.zeros(frame_count, dtype=float)
                actin_heat = np.zeros(frame_count, dtype=float)
                actin_storage = np.zeros(frame_count, dtype=float)
            myosin_active, myosin_heat, myosin_storage = _drift_thermo_chunk(
                mv,
                beta=beta,
                diffusion_coeff=myosin_diff_coeff_trans,
                saved_dt=saved_dt,
                torch_device=torch_device,
            )
            if use_cumulative_chem:
                transition_entropy_frame, prev_cumulative_chem = _cumulative_deltas(
                    cumulative_chem_entropy,
                    frame_slice,
                    prev_cumulative_chem,
                )
                if cumulative_chem_input is not None:
                    chemical_input_frame, prev_cumulative_input = _cumulative_deltas(
                        cumulative_chem_input,
                        frame_slice,
                        prev_cumulative_input,
                    )
                else:
                    chemical_input_frame = np.abs(transition_entropy_frame)
                chunk_binding = 0
                chunk_unbinding = 0
            else:
                status = np.asarray(cb_status[frame_slice, :, 0], dtype=np.int8)
                load = np.asarray(loads[frame_slice, :, 0], dtype=float)
                transition_entropy_frame, chemical_input_frame, chunk_binding, chunk_unbinding = (
                    _transition_entropy_chunk(
                        status,
                        load,
                        prev_status,
                        k_on=k_on,
                        base_lifetime=base_lifetime,
                        lifetime_coeff=lifetime_coeff,
                        directional_base_lifetime=directional_base_lifetime,
                        simulation_dt=simulation_dt,
                        eps=eps,
                        torch_device=torch_device,
                    )
                )
                prev_status = status[-1]

            actin_active_work += float(np.sum(actin_active))
            myosin_active_work += float(np.sum(myosin_active))
            actin_heat_dissipated += float(np.sum(actin_heat))
            myosin_heat_dissipated += float(np.sum(myosin_heat))
            actin_elastic_storage += float(np.sum(actin_storage))
            myosin_elastic_storage += float(np.sum(myosin_storage))
            chemical_entropy += float(np.sum(transition_entropy_frame))
            chemical_input_kbt += float(np.sum(chemical_input_frame))
            binding_events += chunk_binding
            unbinding_events += chunk_unbinding

            frame_entropy_chunks.append(
                beta * actin_heat + beta * myosin_heat + transition_entropy_frame
            )

        if cumulative_binding_events is not None:
            binding_events = int(round(float(cumulative_binding_events[-1] - cumulative_binding_events[0])))
        if cumulative_unbinding_events is not None:
            unbinding_events = int(
                round(float(cumulative_unbinding_events[-1] - cumulative_unbinding_events[0]))
            )

    # One vectorized segmented reduction replaces the Python frame loop.
    all_frame_entropy = (
        np.concatenate(frame_entropy_chunks)
        if frame_entropy_chunks
        else np.empty(0, dtype=float)
    )
    if all_frame_entropy.size:
        starts = np.arange(0, all_frame_entropy.size, segment_frames)
        segment_entropy_arr = (
            np.add.reduceat(all_frame_entropy, starts) / normalization_count
        )
    else:
        segment_entropy_arr = np.empty(0, dtype=float)

    actin_active_work /= normalization_count
    myosin_active_work /= normalization_count
    actin_heat_dissipated /= normalization_count
    myosin_heat_dissipated /= normalization_count
    actin_elastic_storage /= normalization_count
    myosin_elastic_storage /= normalization_count
    chemical_entropy /= normalization_count
    chemical_input_kbt /= normalization_count

    active_work = actin_active_work + myosin_active_work
    heat_dissipated = actin_heat_dissipated + myosin_heat_dissipated
    elastic_storage = actin_elastic_storage + myosin_elastic_storage
    actin_mechanical_entropy = beta * actin_heat_dissipated
    myosin_mechanical_entropy = beta * myosin_heat_dissipated
    mechanical_entropy = actin_mechanical_entropy + myosin_mechanical_entropy
    actin_total_entropy = actin_mechanical_entropy + chemical_entropy
    myosin_total_entropy = myosin_mechanical_entropy
    total_entropy = actin_total_entropy + myosin_total_entropy
    total_time = max((n_frames - 1) * saved_dt, np.finfo(float).eps)
    chemical_input_energy = chemical_input_kbt / beta if beta > 0 else np.nan
    efficiency = elastic_storage / active_work if active_work != 0 else np.nan
    first_law_residual = active_work - heat_dissipated - elastic_storage
    residual_scale = max(
        abs(active_work), abs(heat_dissipated), abs(elastic_storage), np.finfo(float).eps
    )
    first_law_relative_residual = first_law_residual / residual_scale

    if segment_entropy_arr.size:
        ift_terms = np.exp(-np.clip(segment_entropy_arr, -700, 700))
        ift_segment_mean = float(np.mean(ift_terms))
        ift_segment_sem = (
            float(np.std(ift_terms, ddof=1) / np.sqrt(ift_terms.size))
            if ift_terms.size > 1
            else 0.0
        )
    else:
        ift_segment_mean = np.nan
        ift_segment_sem = np.nan

    if use_cumulative_chem:
        note = (
            f"Chemical entropy is read from {chem_entropy_dataset}; the engine "
            "accumulates KMC transition affinities at event time. Per-actin values "
            "split each actin-actin pair event across both actins, so summed actin "
            "values count each transition once. Chemical entropy is assigned to "
            "actin; myosin entropy is mechanical only."
        )
    else:
        note = (
            "Chemical transitions are inferred from saved-frame /actin/cb_status changes; "
            "events that bind and unbind between saved frames are not visible. "
            "Chemical entropy is assigned to actin; myosin entropy is mechanical only."
        )
    if inferred_save_every is not None:
        note += f" Inferred save_every={inferred_save_every}."
    if not include_actin_velocity:
        note += " Actin active-velocity and force cross-term were excluded (actin entropy is chemical only)."

    return EntropyTransductionResult(
        path=path,
        n_frames=n_frames,
        saved_dt=saved_dt,
        beta=beta,
        mechanical_entropy=mechanical_entropy,
        actin_mechanical_entropy=actin_mechanical_entropy,
        myosin_mechanical_entropy=myosin_mechanical_entropy,
        heat_dissipated=heat_dissipated,
        active_work=active_work,
        elastic_storage=elastic_storage,
        actin_active_work=actin_active_work,
        myosin_active_work=myosin_active_work,
        actin_heat_dissipated=actin_heat_dissipated,
        myosin_heat_dissipated=myosin_heat_dissipated,
        actin_elastic_storage=actin_elastic_storage,
        myosin_elastic_storage=myosin_elastic_storage,
        first_law_residual=first_law_residual,
        first_law_relative_residual=first_law_relative_residual,
        chemical_entropy=chemical_entropy,
        chemical_input_kbt=chemical_input_kbt,
        chemical_input_energy=chemical_input_energy,
        actin_total_entropy=actin_total_entropy,
        myosin_total_entropy=myosin_total_entropy,
        total_entropy=total_entropy,
        actin_entropy_rate=actin_total_entropy / total_time,
        myosin_entropy_rate=myosin_total_entropy / total_time,
        entropy_rate=total_entropy / total_time,
        active_work_rate=active_work / total_time,
        efficiency=efficiency,
        binding_events=binding_events,
        unbinding_events=unbinding_events,
        visible_transition_note=note,
        ift_segment_mean=ift_segment_mean,
        ift_segment_sem=ift_segment_sem,
    )


def compute_entropy_rate_timeseries(
    h5_path: str | Path,
    *,
    beta: float = 241.0,
    simulation_dt: float = 1e-5,
    actin_diff_coeff_trans: float = 0.2,
    myosin_diff_coeff_trans: float = 0.05,
    save_every: int | None = None,
    k_on: float = 1000.0,
    base_lifetime: float = 0.0,
    lifetime_coeff: float = 10**-2.5,
    directional_base_lifetime: float = 0.0,
    load_dataset: str = "/actin/f_load_cb",
    chem_entropy_dataset: str = DEFAULT_CHEM_ENTROPY_DATASET,
    include_actin_velocity: bool = True,
    chunk_frames: int = 128,
    smooth_window: int = 25,
    normalization_count: float = 1.0,
    device: str = "auto",
) -> EntropyRateTimeseries:
    """Return time-resolved entropy production and energy rates.

    Heavy particle and transition reductions are vectorized for each HDF5
    chunk. ``device="auto"`` uses torch/CUDA if available, otherwise NumPy.
    Rates are assigned to saved-frame intervals. ``smooth_window`` applies a
    centered moving average over saved intervals; use 1 for unsmoothed rates.
    Mechanical entropy uses the motor-friction estimator β*γ|V|²*saved_dt
    for both species. Chemical entropy is actin-only; myosin total entropy is
    mechanical only. Engine-side cumulative chemical entropy is used when
    available; otherwise visible saved-frame transitions are the fallback.
    Pass ``include_actin_velocity=False`` for chemical-only actin.
    """
    path = Path(h5_path)
    if chunk_frames < 1:
        raise ValueError("chunk_frames must be positive")
    if smooth_window < 1:
        raise ValueError("smooth_window must be positive")
    if normalization_count <= 0:
        raise ValueError("normalization_count must be positive")
    torch_device = resolve_compute_device(device)

    with h5py.File(path, "r") as fh:
        actin_force = _require(fh, "/actin/force")
        actin_velocity = _require(fh, "/actin/velocity")
        myosin_force = _require(fh, "/myosin/force")
        myosin_velocity = _require(fh, "/myosin/velocity")
        cb_status = _require(fh, "/actin/cb_status")
        load_path = load_dataset if _has_dataset(fh, load_dataset) else "/actin/f_load"
        loads = _require(fh, load_path)

        n_frames = int(cb_status.shape[0])
        saved_dt, _ = _infer_saved_dt(fh, simulation_dt, save_every)
        cumulative_chem_entropy = _read_cumulative_frame_sums(
            fh, chem_entropy_dataset, n_frames
        )
        if min(
            actin_force.shape[0],
            actin_velocity.shape[0],
            myosin_force.shape[0],
            myosin_velocity.shape[0],
        ) < n_frames:
            raise ValueError("Force/velocity datasets have fewer frames than /actin/cb_status")

        n_intervals = max(n_frames - 1, 0)
        time = np.arange(1, n_frames, dtype=float) * saved_dt
        actin_mechanical_entropy = np.zeros(n_intervals, dtype=float)
        myosin_mechanical_entropy = np.zeros(n_intervals, dtype=float)
        chemical_entropy = np.zeros(n_intervals, dtype=float)
        actin_heat = np.zeros(n_intervals, dtype=float)
        myosin_heat = np.zeros(n_intervals, dtype=float)
        actin_active_work = np.zeros(n_intervals, dtype=float)
        myosin_active_work = np.zeros(n_intervals, dtype=float)
        actin_elastic_storage = np.zeros(n_intervals, dtype=float)
        myosin_elastic_storage = np.zeros(n_intervals, dtype=float)
        if cumulative_chem_entropy is not None and n_intervals:
            chemical_entropy[:] = np.diff(cumulative_chem_entropy)

        prev_status = np.asarray(cb_status[0, :, 0], dtype=np.int8)
        eps = np.finfo(float).tiny

        for frame_slice in _iter_frame_chunks(n_frames, chunk_frames):
            offset = frame_slice.start - 1
            mv = np.asarray(myosin_velocity[frame_slice], dtype=float)
            if include_actin_velocity:
                av = np.asarray(actin_velocity[frame_slice], dtype=float)
                actin_active, actin_chunk_heat, actin_storage = _drift_thermo_chunk(
                    av,
                    beta=beta,
                    diffusion_coeff=actin_diff_coeff_trans,
                    saved_dt=saved_dt,
                    torch_device=torch_device,
                )
            else:
                frame_count = frame_slice.stop - frame_slice.start
                actin_active = np.zeros(frame_count, dtype=float)
                actin_chunk_heat = np.zeros(frame_count, dtype=float)
                actin_storage = np.zeros(frame_count, dtype=float)
            myosin_active, myosin_chunk_heat, myosin_storage = _drift_thermo_chunk(
                mv,
                beta=beta,
                diffusion_coeff=myosin_diff_coeff_trans,
                saved_dt=saved_dt,
                torch_device=torch_device,
            )
            if cumulative_chem_entropy is None:
                status = np.asarray(cb_status[frame_slice, :, 0], dtype=np.int8)
                load = np.asarray(loads[frame_slice, :, 0], dtype=float)
                transition_entropy, _, _, _ = _transition_entropy_chunk(
                    status,
                    load,
                    prev_status,
                    k_on=k_on,
                    base_lifetime=base_lifetime,
                    lifetime_coeff=lifetime_coeff,
                    directional_base_lifetime=directional_base_lifetime,
                    simulation_dt=simulation_dt,
                    eps=eps,
                    torch_device=torch_device,
                )
            else:
                transition_entropy = chemical_entropy[offset:offset + actin_chunk_heat.size]

            stop = offset + actin_chunk_heat.size
            actin_heat[offset:stop] = actin_chunk_heat
            myosin_heat[offset:stop] = myosin_chunk_heat
            actin_active_work[offset:stop] = actin_active
            myosin_active_work[offset:stop] = myosin_active
            actin_elastic_storage[offset:stop] = actin_storage
            myosin_elastic_storage[offset:stop] = myosin_storage
            actin_mechanical_entropy[offset:stop] = beta * actin_chunk_heat
            myosin_mechanical_entropy[offset:stop] = beta * myosin_chunk_heat
            chemical_entropy[offset:stop] = transition_entropy
            if cumulative_chem_entropy is None:
                prev_status = status[-1]

    def rate(values: np.ndarray) -> np.ndarray:
        out = values / (saved_dt * normalization_count)
        if smooth_window <= 1 or out.size == 0:
            return out
        kernel = np.ones(smooth_window, dtype=float) / float(smooth_window)
        return np.convolve(out, kernel, mode="same")

    mechanical_entropy = actin_mechanical_entropy + myosin_mechanical_entropy
    heat = actin_heat + myosin_heat
    active_work = actin_active_work + myosin_active_work
    elastic_storage = actin_elastic_storage + myosin_elastic_storage
    actin_total_entropy = actin_mechanical_entropy + chemical_entropy
    myosin_total_entropy = myosin_mechanical_entropy

    return EntropyRateTimeseries(
        time=time,
        mechanical_rate=rate(mechanical_entropy),
        actin_mechanical_rate=rate(actin_mechanical_entropy),
        myosin_mechanical_rate=rate(myosin_mechanical_entropy),
        chemical_rate=rate(chemical_entropy),
        actin_chemical_rate=rate(chemical_entropy),
        total_rate=rate(actin_total_entropy + myosin_total_entropy),
        actin_total_rate=rate(actin_total_entropy),
        myosin_total_rate=rate(myosin_total_entropy),
        heat_rate=rate(heat),
        actin_heat_rate=rate(actin_heat),
        myosin_heat_rate=rate(myosin_heat),
        active_work_rate=rate(active_work),
        actin_active_work_rate=rate(actin_active_work),
        myosin_active_work_rate=rate(myosin_active_work),
        elastic_storage_rate=rate(elastic_storage),
        actin_elastic_storage_rate=rate(actin_elastic_storage),
        myosin_elastic_storage_rate=rate(myosin_elastic_storage),
    )


def _result_key(record: dict[str, object]) -> tuple[str, str, str]:
    return (str(record["Filename"]), str(record["A_val"]), str(record["lt_str"]))


def _validate_csv_header(path: str | Path, required_fields: list[str]) -> None:
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return

    with open(path, mode="r", newline="") as csvfile:
        fieldnames = csv.DictReader(csvfile).fieldnames or []

    missing = [field for field in required_fields if field not in fieldnames]
    if missing:
        raise ValueError(
            f"{path} is missing required output columns {missing}. "
            "Write to a new CSV or regenerate this file."
        )


def _load_existing_results(output_csv: str | Path) -> set[tuple[str, str, str]]:
    existing = set()
    output_path = Path(output_csv)
    if not output_path.exists() or output_path.stat().st_size == 0:
        return existing

    with open(output_path, mode="r", newline="") as csvfile:
        for row in csv.DictReader(csvfile):
            if not all(field in row and row[field] for field in ("Filename", "A_val", "lt_str")):
                continue
            existing.add(_result_key(row))
    return existing


def _parse_condition_from_name(name: str) -> tuple[str, int, str, str] | None:
    match = re.search(r"(.*A(\d+)_M\d+_lt([^_]+?)(?:_[^_]+)*)_r\d+(.*)", name)
    if not match:
        return None
    return match.group(1), int(match.group(2)), match.group(3), match.group(4)


def _iter_h5_inputs(paths: Iterable[str | Path]) -> Iterable[Path]:
    for path in paths:
        path = Path(path)
        if not path.exists():
            continue
        if path.is_file() and path.suffix == ".h5":
            yield path
        elif path.is_dir():
            yield from sorted(path.glob("*.h5"))


def _h5_frame_count(h5_path: Path) -> int | None:
    with h5py.File(h5_path, "r") as fh:
        candidates = []
        for dataset in (
            "/actin/cb_status",
            "/actin/force",
            "/actin/velocity",
            "/myosin/force",
            "/myosin/velocity",
        ):
            if dataset not in fh:
                return None
            candidates.append(fh[dataset].shape[0])
        return min(candidates)


def _mean_tail(values: np.ndarray, eval_window: int, stride: int = 1) -> float:
    if values.size == 0:
        return np.nan
    return float(np.mean(values[-eval_window::stride]))


def entropy_csv_record(
    h5_file: str | Path,
    *,
    filename: str | None = None,
    a_val: int | None = None,
    lt_str: str | None = None,
    eval_window: int = 300,
    eval_stride: int = 1,
    smooth_window: int = 1,
    **kwargs,
) -> dict[str, object]:
    """Return one CSV-ready entropy record for an HDF5 trajectory."""
    if eval_window < 1:
        raise ValueError("eval_window must be positive")
    if eval_stride < 1:
        raise ValueError("eval_stride must be positive")

    h5_file = Path(h5_file)
    result = compute_entropy_transduction(h5_file, **kwargs)
    ts_kwargs = dict(kwargs)
    ts_kwargs.pop("segment_frames", None)
    ts_kwargs.pop("chem_input_dataset", None)
    ts_kwargs.pop("binding_events_dataset", None)
    ts_kwargs.pop("unbinding_events_dataset", None)
    ts = compute_entropy_rate_timeseries(
        h5_file,
        smooth_window=smooth_window,
        **ts_kwargs,
    )

    return {
        "Filename": filename or h5_file.name,
        "H5_Path": str(h5_file.resolve()),
        "A_val": "" if a_val is None else a_val,
        "lt_str": "" if lt_str is None else lt_str,
        "N_Frames": result.n_frames,
        "Saved_dt": result.saved_dt,
        "Heat_Dissipated": result.heat_dissipated,
        "Active_Work": result.active_work,
        "Elastic_Storage": result.elastic_storage,
        "Actin_Mech_Entropy": result.actin_mechanical_entropy,
        "Myosin_Mech_Entropy": result.myosin_mechanical_entropy,
        "Mech_Entropy": result.mechanical_entropy,
        "Chem_Entropy": result.chemical_entropy,
        "Actin_Total_Entropy": result.actin_total_entropy,
        "Myosin_Total_Entropy": result.myosin_total_entropy,
        "Total_Entropy": result.total_entropy,
        "Actin_Entropy_Rate": result.actin_entropy_rate,
        "Myosin_Entropy_Rate": result.myosin_entropy_rate,
        "Entropy_Rate": result.entropy_rate,
        "Active_Work_Rate": result.active_work_rate,
        "Mech_Rate": _mean_tail(ts.mechanical_rate, eval_window, eval_stride),
        "Chem_Rate": _mean_tail(ts.chemical_rate, eval_window, eval_stride),
        "Total_Rate": _mean_tail(ts.total_rate, eval_window, eval_stride),
        "Actin_Total_Rate": _mean_tail(ts.actin_total_rate, eval_window, eval_stride),
        "Myosin_Total_Rate": _mean_tail(ts.myosin_total_rate, eval_window, eval_stride),
        "Binding_Events": result.binding_events,
        "Unbinding_Events": result.unbinding_events,
        "IFT_Segment_Mean": result.ift_segment_mean,
        "IFT_Segment_SEM": result.ift_segment_sem,
    }


def compute_energy_transduction_results(
    h5_files: Iterable[str | Path],
    output_csv: str | Path | None = None,
    *,
    threshold: int | None = None,
    eval_window: int = 300,
    eval_stride: int = 1,
    smooth_window: int = 1,
    infer_conditions: bool = True,
    skip_existing: bool = True,
    raise_errors: bool = False,
    progress: bool = True,
    **kwargs,
) -> list[dict[str, object]]:
    """Compute energy-transduction records for an explicit list of HDF5 files.

    This is the list-oriented counterpart to ``process_entropy_directories``:
    it does not renumber runs or scan default data directories. Each successful
    input returns one CSV-ready record. If ``output_csv`` is provided, records
    are appended as they are computed.
    """
    if eval_window < 1:
        raise ValueError("eval_window must be positive")
    if threshold is not None and threshold < 1:
        raise ValueError("threshold must be positive when provided")

    output_path = Path(output_csv) if output_csv is not None else None
    existing_results = (
        _load_existing_results(output_path) if output_path and skip_existing else set()
    )
    existing_filenames: set[str] = set()
    if (
        output_path
        and skip_existing
        and output_path.exists()
        and output_path.stat().st_size > 0
    ):
        with open(output_path, mode="r", newline="") as csvfile:
            existing_filenames = {
                row["Filename"]
                for row in csv.DictReader(csvfile)
                if "Filename" in row and row["Filename"]
            }
    records: list[dict[str, object]] = []
    processed_count = 0
    skipped_count = 0
    discarded_count = 0
    error_count = 0

    if progress:
        requested_device = str(kwargs.get("device", "auto"))
        print_compute_device_info(resolve_compute_device(requested_device), requested_device)
        print(f"{'Filename':<55} | {'Frames':<8} | {'Status'}")
        print("-" * 85)

    if output_path is not None:
        _validate_csv_header(output_path, ENTROPY_CSV_FIELDS)
        write_header = not output_path.exists() or output_path.stat().st_size == 0
        csv_handle = open(output_path, mode="a", newline="")
    else:
        write_header = False
        csv_handle = None

    try:
        writer = (
            csv.DictWriter(csv_handle, fieldnames=ENTROPY_CSV_FIELDS)
            if csv_handle
            else None
        )
        if writer is not None and write_header:
            writer.writeheader()
            csv_handle.flush()

        for h5_file in h5_files:
            path = Path(h5_file)
            name = path.name
            try:
                if not path.exists() or not path.is_file() or path.suffix != ".h5":
                    skipped_count += 1
                    if progress:
                        print(f"{name:<55} | {'-':<8} | Skipped (not an HDF5 file)")
                    continue

                n_frames = _h5_frame_count(path)
                if n_frames is None:
                    skipped_count += 1
                    if progress:
                        print(f"{name:<55} | {'-':<8} | Skipped (missing required datasets)")
                    continue
                if threshold is not None and n_frames < threshold:
                    discarded_count += 1
                    if progress:
                        print(f"{name:<55} | {n_frames:<8} | Discarded (short run)")
                    continue

                a_val: int | None = None
                lt_str: str | None = None
                if infer_conditions:
                    parsed = _parse_condition_from_name(name)
                    if parsed is not None:
                        _, a_val, lt_str, _ = parsed

                record_key = (
                    name,
                    "" if a_val is None else str(a_val),
                    "" if lt_str is None else lt_str,
                )
                if record_key in existing_results or name in existing_filenames:
                    skipped_count += 1
                    if progress:
                        print(f"{name:<55} | {n_frames:<8} | Skipped (in CSV)")
                    continue

                record = entropy_csv_record(
                    path,
                    filename=name,
                    a_val=a_val,
                    lt_str=lt_str,
                    eval_window=eval_window,
                    eval_stride=eval_stride,
                    smooth_window=smooth_window,
                    **kwargs,
                )
                records.append(record)
                processed_count += 1

                if writer is not None:
                    writer.writerow(record)
                    csv_handle.flush()
                    existing_results.add(_result_key(record))
                    existing_filenames.add(str(record["Filename"]))

                if progress:
                    print(f"{name:<55} | {n_frames:<8} | Processed")

            except Exception as exc:
                error_count += 1
                if raise_errors:
                    raise
                if progress:
                    print(f"{name:<55} | {'ERROR':<8} | {exc}")

    finally:
        if csv_handle is not None:
            csv_handle.close()

    if progress:
        print("-" * 85)
        message = f"Energy transduction extraction complete. Computed {processed_count} runs"
        if output_path is not None:
            message += f"; appended results to {output_path}"
        message += (
            f"; skipped {skipped_count}; discarded {discarded_count} short runs; "
            f"errors {error_count}."
        )
        print(message)

    return records


def process_entropy_directories(
    base_dirs: Iterable[str | Path],
    output_csv: str | Path = "entropy_transduction_results.csv",
    *,
    threshold: int = 2880,
    eval_window: int = 300,
    eval_stride: int = 1,
    smooth_window: int = 1,
    **kwargs,
) -> None:
    """Process all ``*.h5`` files in one or more directories into a CSV."""
    output_path = Path(output_csv)
    existing_results = _load_existing_results(output_path)
    condition_counters: defaultdict[str, int] = defaultdict(int)
    discarded_count = 0
    processed_count = 0
    skipped_count = 0

    requested_device = str(kwargs.get("device", "auto"))
    print_compute_device_info(resolve_compute_device(requested_device), requested_device)

    _validate_csv_header(output_path, ENTROPY_CSV_FIELDS)
    write_header = not output_path.exists() or output_path.stat().st_size == 0

    print(f"{'Original/Pseudo-Filename':<55} | {'Frames':<8} | {'Status'}")
    print("-" * 85)

    with open(output_path, mode="a", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=ENTROPY_CSV_FIELDS)
        if write_header:
            writer.writeheader()
            csvfile.flush()

        for h5_file in _iter_h5_inputs(base_dirs):
            name = h5_file.name
            try:
                n_frames = _h5_frame_count(h5_file)
                if n_frames is None:
                    continue
                if n_frames < threshold:
                    discarded_count += 1
                    continue

                parsed = _parse_condition_from_name(name)
                if parsed is None:
                    continue
                condition_base, a_val, lt_str, extension = parsed
                condition_counters[condition_base] += 1
                pseudo_name = f"{condition_base}_r{condition_counters[condition_base]:02d}{extension}"
                key = (pseudo_name, str(a_val), lt_str)

                if key in existing_results:
                    skipped_count += 1
                    print(f"{pseudo_name:<55} | {n_frames:<8} | Skipped (in CSV)")
                    continue

                record = entropy_csv_record(
                    h5_file,
                    filename=pseudo_name,
                    a_val=a_val,
                    lt_str=lt_str,
                    eval_window=eval_window,
                    eval_stride=eval_stride,
                    smooth_window=smooth_window,
                    **kwargs,
                )
                writer.writerow(record)
                csvfile.flush()
                existing_results.add(_result_key(record))
                processed_count += 1

                print(f"{pseudo_name:<55} | {n_frames:<8} | Processed")

            except Exception as exc:
                print(f"{name:<55} | {'ERROR':<8} | {exc}")

    print("-" * 85)
    print(
        f"Entropy extraction complete. Added {processed_count} runs to {output_path}; "
        f"skipped {skipped_count} already-computed runs; discarded {discarded_count} short runs."
    )


def print_summary(result: EntropyTransductionResult) -> None:
    print(f"path: {result.path}")
    print(f"frames: {result.n_frames}, saved_dt: {result.saved_dt:.6g}, beta: {result.beta:.6g}")
    print(f"heat dissipated: {result.heat_dissipated:.6g}")
    print(
        "mechanical entropy beta*Q: "
        f"{result.mechanical_entropy:.6g} "
        f"(actin {result.actin_mechanical_entropy:.6g}, "
        f"myosin {result.myosin_mechanical_entropy:.6g})"
    )
    print(
        "species total entropy: "
        f"actin {result.actin_total_entropy:.6g} "
        f"(mechanical + chemical), myosin {result.myosin_total_entropy:.6g} "
        "(mechanical only)"
    )
    print(
        "active work input: "
        f"{result.active_work:.6g} "
        f"(actin {result.actin_active_work:.6g}, myosin {result.myosin_active_work:.6g})"
    )
    print(
        "elastic storage/transduction: "
        f"{result.elastic_storage:.6g} "
        f"(actin {result.actin_elastic_storage:.6g}, myosin {result.myosin_elastic_storage:.6g})"
    )
    print(
        "first-law residual W_active - Q - dU: "
        f"{result.first_law_residual:.6g} "
        f"(relative {result.first_law_relative_residual:.3g})"
    )
    print(f"chemical entropy from visible actin transitions: {result.chemical_entropy:.6g}")
    print(f"total entropy: {result.total_entropy:.6g}")
    print(
        "entropy rate: "
        f"{result.entropy_rate:.6g} "
        f"(actin {result.actin_entropy_rate:.6g}, myosin {result.myosin_entropy_rate:.6g})"
    )
    print(
        "chemical input: "
        f"{result.chemical_input_kbt:.6g} kBT = "
        f"{result.chemical_input_energy:.6g} energy units"
    )
    print(f"transduction efficiency dU / W_active: {result.efficiency:.6g}")
    print(f"visible binding/unbinding events: {result.binding_events}/{result.unbinding_events}")
    print(
        "segment IFT proxy <exp(-Delta S)>: "
        f"{result.ift_segment_mean:.6g} +/- {result.ift_segment_sem:.3g}"
    )
    print(result.visible_transition_note)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "paths",
        nargs="*",
        type=Path,
        help=(
            "Optional HDF5 file(s) or directory/directories containing *.h5 files. "
            "Use --base_dirs for an explicit extract_structure_3d-style directory list."
        ),
    )
    parser.add_argument(
        "--base_dirs",
        "--data-root",
        dest="base_dirs",
        type=Path,
        nargs="+",
        default=None,
        help=(
            "Base directories to process. If omitted, positional paths are used; "
            "if those are also omitted, defaults to data0416 and data0430."
        ),
    )
    parser.add_argument("--out_csv", type=Path, default="entropy_transduction_results.csv")
    parser.add_argument("--threshold", type=int, default=1000)
    parser.add_argument("--eval_window", type=int, default=300)
    parser.add_argument("--eval_stride", type=int, default=10,
                        help="Stride when sampling the tail window for steady-state rate estimate"
                             " (default 10 → 30 quasi-independent samples from last 300 frames)")
    parser.add_argument("--beta", type=float, default=241.0)
    parser.add_argument("--simulation_dt", type=float, default=1e-5)
    parser.add_argument("--actin_diff_coeff_trans", type=float, default=0.2)
    parser.add_argument("--myosin_diff_coeff_trans", type=float, default=0.05)
    parser.add_argument("--save_every", type=int, default=None)
    parser.add_argument("--k_on", type=float, default=1000.0)
    parser.add_argument("--base_lifetime", type=float, default=0.0)
    parser.add_argument("--lifetime_coeff", type=float, default=10**-2.5)
    parser.add_argument("--directional_base_lifetime", type=float, default=0.0)
    parser.add_argument("--load_dataset", default="/actin/f_load_cb")
    parser.add_argument("--chem_entropy_dataset", default=DEFAULT_CHEM_ENTROPY_DATASET)
    parser.add_argument("--chem_input_dataset", default=DEFAULT_CHEM_INPUT_DATASET)
    parser.add_argument("--binding_events_dataset", default=DEFAULT_BINDING_EVENTS_DATASET)
    parser.add_argument("--unbinding_events_dataset", default=DEFAULT_UNBINDING_EVENTS_DATASET)
    parser.add_argument(
        "--no_actin_velocity",
        dest="include_actin_velocity",
        action="store_false",
        help=(
            "Exclude actin active-velocity and force cross-term; actin entropy becomes "
            "chemical only. By default both actin and myosin use the full Sekimoto formula."
        ),
    )
    parser.set_defaults(include_actin_velocity=True)
    parser.add_argument("--chunk_frames", type=int, default=128)
    parser.add_argument("--segment_frames", type=int, default=100)
    parser.add_argument("--normalization_count", type=float, default=1.0)
    parser.add_argument(
        "--device",
        choices=("auto", "numpy", "cpu", "cuda"),
        default="auto",
        help="Computation backend: auto uses CUDA when torch can access it, else NumPy.",
    )
    parser.add_argument("--smooth_window", type=int, default=1)
    args = parser.parse_args()

    paths = args.base_dirs if args.base_dirs is not None else args.paths or list(DEFAULT_BASE_DIRS)
    compute_kwargs = {
        "beta": args.beta,
        "simulation_dt": args.simulation_dt,
        "actin_diff_coeff_trans": args.actin_diff_coeff_trans,
        "myosin_diff_coeff_trans": args.myosin_diff_coeff_trans,
        "save_every": args.save_every,
        "k_on": args.k_on,
        "base_lifetime": args.base_lifetime,
        "lifetime_coeff": args.lifetime_coeff,
        "directional_base_lifetime": args.directional_base_lifetime,
        "load_dataset": args.load_dataset,
        "chem_entropy_dataset": args.chem_entropy_dataset,
        "chem_input_dataset": args.chem_input_dataset,
        "binding_events_dataset": args.binding_events_dataset,
        "unbinding_events_dataset": args.unbinding_events_dataset,
        "include_actin_velocity": args.include_actin_velocity,
        "chunk_frames": args.chunk_frames,
        "segment_frames": args.segment_frames,
        "normalization_count": args.normalization_count,
        "device": args.device,
    }

    if len(paths) == 1 and paths[0].is_file():
        print_compute_device_info(resolve_compute_device(args.device), args.device)
        print_summary(compute_entropy_transduction(paths[0], **compute_kwargs))
    else:
        process_entropy_directories(
            paths,
            args.out_csv,
            threshold=args.threshold,
            eval_window=args.eval_window,
            eval_stride=args.eval_stride,
            smooth_window=args.smooth_window,
            **compute_kwargs,
        )
