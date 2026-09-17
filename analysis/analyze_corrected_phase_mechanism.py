#!/usr/bin/env python3
"""Mechanism analysis for the corrected n_actins--lambda phase sweep.

The campaign uses a constant 1-ms lifetime for status-1 ACTN2 contacts and
adds ``lambda * f_load`` only to status-2 links.  This script reads each large
HDF5 trajectory once, caches compact per-frame observables, and then builds
run-level summaries and threshold tests.

The temporal structural observable is an inexpensive axial M/Z order parameter
over all PBC-compatible harmonics nx=1..30.  The already-computed, fully 3-D
legacy Phi_3D tail metric can be merged with ``--legacy-phi-csv``.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import h5py
import numpy as np


PATTERN = re.compile(
    r"^A(?P<n_actins>\d+)_M(?P<n_myosins>\d+)_s2fstab_lt"
    r"(?P<lambda_label>[^_]+)_r(?P<replica>\d+)_L12x2x2\.h5$"
)
BASE_LIFETIME = 1.0e-3
DT = 1.0e-5
BOX = np.array([12.0, 2.0, 2.0])

TIME_FIELDS = (
    "filename", "n_actins", "n_myosins", "lambda_label", "lambda_value",
    "replica", "frame", "step", "time_s", "phi_axial_myosin",
    "nx_axial_myosin", "phi_axial_cb2", "nx_axial_cb2",
    "phi_shared_myosin", "phi_shared_cb2", "nx_shared",
    "registration_pi", "q_sarcomere_axial", "nematic_myosin",
    "nematic_actin", "nematic_cb2", "cb1_fraction", "cb2_fraction",
    "cb_any_fraction", "mean_cb2_load", "median_cb2_load",
    "directional_extra_over_base", "fraction_cb2_directional_gt_base",
    "mean_myosin_binding_ratio", "fraction_am_engaged",
    "mean_partial_binding_ratio", "mean_crosslink_ratio",
    "binding_events_per_actin", "unbinding_events_per_actin",
    "binding_rate_per_actin_s", "unbinding_rate_per_actin_s",
    "myosin_rms_speed", "actin_rms_speed",
)

SUMMARY_FIELDS = (
    "filename", "n_actins", "n_myosins", "lambda_label", "lambda_value",
    "lambda_over_base", "replica", "checkpoint_step", "complete",
    "n_sampled_frames", "legacy_phi_3d", "legacy_theta_deg",
    "tail_phi_axial_myosin", "tail_phi_axial_cb2",
    "tail_q_sarcomere_axial", "tail_nematic_myosin", "tail_nematic_actin",
    "tail_nematic_cb2", "tail_cb1_fraction", "tail_cb2_fraction",
    "tail_mean_cb2_load", "tail_directional_extra_over_base",
    "tail_fraction_cb2_directional_gt_base", "tail_am_engaged_fraction",
    "tail_binding_rate_per_actin_s", "tail_unbinding_rate_per_actin_s",
    "onset_phi_0p4_s", "onset_phi_0p5_s", "onset_phi_0p6_s",
    "classification",
)


def lambda_value(label: str) -> float:
    if label == "0":
        return 0.0
    mantissa, exponent = label.split("e", 1)
    return float(mantissa) * 10.0 ** float(exponent)


def checkpoint_step(handle: h5py.File) -> int:
    for key in ("/resume/current_step", "/resume_prev/current_step"):
        if key in handle:
            return int(np.asarray(handle[key]).reshape(-1)[-1])
    return int(np.asarray(handle["/state/current_step"])[-1, 0])


def null_corrected(phi: np.ndarray | float, n_particles: int) -> np.ndarray:
    values = np.asarray(phi, dtype=float)
    if n_particles <= 0:
        return np.zeros_like(values)
    null = math.sqrt(math.pi) / (2.0 * math.sqrt(n_particles))
    return np.maximum(0.0, (values - null) / (1.0 - null))


def nematic_order(directions: np.ndarray, mask: np.ndarray | None = None) -> np.ndarray:
    """Largest eigenvalue of Q=(3<uu>-I)/2 for every frame."""
    output = np.zeros(directions.shape[0], dtype=float)
    if mask is None:
        second = np.einsum("fni,fnj->fij", directions, directions, optimize=True)
        second /= max(1, directions.shape[1])
        q = 1.5 * second - 0.5 * np.eye(3)[None, :, :]
        return np.linalg.eigvalsh(q)[:, -1]
    for frame in range(directions.shape[0]):
        selected = directions[frame, mask[frame]]
        if len(selected) == 0:
            continue
        second = selected.T @ selected / len(selected)
        output[frame] = np.linalg.eigvalsh(1.5 * second - 0.5 * np.eye(3))[-1]
    return output


def axial_orders(
    myosin_x: np.ndarray, actin_x: np.ndarray, cb2: np.ndarray,
    n_modes: int = 30,
) -> dict[str, np.ndarray]:
    modes = np.arange(1, n_modes + 1, dtype=int)
    kx = 2.0 * np.pi * modes / BOX[0]
    phase_m = np.exp(1j * myosin_x[:, :, None] * kx[None, None, :]).mean(axis=1)
    phi_m_all = np.abs(phase_m)
    best_m = np.argmax(phi_m_all, axis=1)
    frame_index = np.arange(len(myosin_x))

    phi_a_all = np.zeros((len(myosin_x), n_modes), dtype=float)
    phase_a_all = np.zeros((len(myosin_x), n_modes), dtype=np.complex128)
    n_cb2 = cb2.sum(axis=1)
    for frame in range(len(myosin_x)):
        selected = actin_x[frame, cb2[frame]]
        if len(selected) == 0:
            continue
        phase = np.exp(1j * selected[:, None] * kx[None, :]).mean(axis=0)
        phase_a_all[frame] = phase
        phi_a_all[frame] = np.abs(phase)
    best_a = np.argmax(phi_a_all, axis=1)

    corrected_m = null_corrected(phi_m_all, myosin_x.shape[1])
    corrected_a = np.zeros_like(phi_a_all)
    for frame, count in enumerate(n_cb2):
        corrected_a[frame] = null_corrected(phi_a_all[frame], int(count))
    score = np.sqrt(corrected_m * corrected_a)
    best_shared = np.argmax(score, axis=1)
    psi_m = phase_m[frame_index, best_shared]
    psi_a = phase_a_all[frame_index, best_shared]
    delta = np.angle(psi_m * np.conj(psi_a))
    registration_pi = -np.cos(delta)
    registration_factor = 0.5 * (1.0 + registration_pi)

    return {
        "phi_m": phi_m_all[frame_index, best_m],
        "nx_m": modes[best_m],
        "phi_a": phi_a_all[frame_index, best_a],
        "nx_a": modes[best_a],
        "phi_shared_m": phi_m_all[frame_index, best_shared],
        "phi_shared_a": phi_a_all[frame_index, best_shared],
        "nx_shared": modes[best_shared],
        "registration_pi": registration_pi,
        "q_sarc": score[frame_index, best_shared] * registration_factor,
    }


def rms_speed(centers: np.ndarray, times: np.ndarray) -> np.ndarray:
    result = np.zeros(len(times), dtype=float)
    if len(times) < 2:
        return result
    delta = np.diff(centers, axis=0)
    delta -= BOX[None, None, :] * np.round(delta / BOX[None, None, :])
    dt = np.diff(times)
    valid = dt > 0
    speed = np.zeros(len(dt), dtype=float)
    if np.any(valid):
        squared = np.sum(delta[valid] ** 2, axis=2)
        speed[valid] = np.sqrt(np.mean(squared, axis=1)) / dt[valid]
    result[1:] = speed
    result[0] = result[1] if len(result) > 1 else 0.0
    return result


def read_optional(handle: h5py.File, key: str, indices: np.ndarray, shape: tuple[int, ...]) -> np.ndarray:
    if key not in handle:
        return np.zeros(shape, dtype=float)
    return np.asarray(handle[key][indices], dtype=float).reshape(shape)


def extract_file(path: Path, stride: int) -> list[dict[str, object]]:
    match = PATTERN.match(path.name)
    if match is None:
        return []
    n_actins = int(match["n_actins"])
    n_myosins = int(match["n_myosins"])
    label = match["lambda_label"]
    lam = lambda_value(label)
    replica = int(match["replica"])

    with h5py.File(path, "r", locking=False) as handle:
        required = (
            "/actin/cb_status", "/actin/center", "/actin/direction",
            "/myosin/center", "/myosin/direction", "/state/current_step",
        )
        n_frames = min(handle[key].shape[0] for key in required)
        indices = np.arange(0, n_frames, stride, dtype=int)
        steps = np.asarray(handle["/state/current_step"][indices]).reshape(-1).astype(int)
        times = steps * DT
        status = np.asarray(handle["/actin/cb_status"][indices]).reshape(len(indices), n_actins)
        cb1 = status == 1
        cb2 = status == 2
        actin_center = np.asarray(handle["/actin/center"][indices], dtype=float)
        actin_direction = np.asarray(handle["/actin/direction"][indices], dtype=float)
        myosin_center = np.asarray(handle["/myosin/center"][indices], dtype=float)
        myosin_direction = np.asarray(handle["/myosin/direction"][indices], dtype=float)

        shape = (len(indices), n_actins)
        load = read_optional(handle, "/actin/f_load_cb", indices, shape)
        mbr = read_optional(handle, "/actin/myosin_binding_ratio", indices, shape)
        partial = read_optional(handle, "/actin/partial_binding_ratio", indices, shape)
        crosslink = read_optional(handle, "/actin/crosslink_ratio", indices, shape)
        binding = read_optional(handle, "/actin/cumulative_chem_binding_events", indices, shape)
        unbinding = read_optional(handle, "/actin/cumulative_chem_unbinding_events", indices, shape)

    order = axial_orders(myosin_center[:, :, 0], actin_center[:, :, 0], cb2)
    nem_m = nematic_order(myosin_direction)
    nem_a = nematic_order(actin_direction)
    nem_cb2 = nematic_order(actin_direction, cb2)
    n_cb2 = cb2.sum(axis=1)
    mean_load = np.divide((load * cb2).sum(axis=1), n_cb2, out=np.zeros(len(indices)), where=n_cb2 > 0)
    median_load = np.zeros(len(indices), dtype=float)
    fraction_directional = np.zeros(len(indices), dtype=float)
    for frame in range(len(indices)):
        values = load[frame, cb2[frame]]
        if len(values):
            median_load[frame] = float(np.median(values))
            fraction_directional[frame] = float(np.mean(lam * values > BASE_LIFETIME))

    bind_total = binding.sum(axis=1) / n_actins
    unbind_total = unbinding.sum(axis=1) / n_actins
    bind_rate = np.zeros(len(indices), dtype=float)
    unbind_rate = np.zeros(len(indices), dtype=float)
    if len(indices) > 1:
        delta_t = np.diff(times)
        valid = delta_t > 0
        bind_rate[1:][valid] = np.diff(bind_total)[valid] / delta_t[valid]
        unbind_rate[1:][valid] = np.diff(unbind_total)[valid] / delta_t[valid]
        bind_rate[0] = bind_rate[1]
        unbind_rate[0] = unbind_rate[1]

    speed_m = rms_speed(myosin_center, times)
    speed_a = rms_speed(actin_center, times)
    records = []
    for out_index, frame in enumerate(indices):
        records.append({
            "filename": path.name,
            "n_actins": n_actins,
            "n_myosins": n_myosins,
            "lambda_label": label,
            "lambda_value": lam,
            "replica": replica,
            "frame": int(frame),
            "step": int(steps[out_index]),
            "time_s": float(times[out_index]),
            "phi_axial_myosin": float(order["phi_m"][out_index]),
            "nx_axial_myosin": int(order["nx_m"][out_index]),
            "phi_axial_cb2": float(order["phi_a"][out_index]),
            "nx_axial_cb2": int(order["nx_a"][out_index]),
            "phi_shared_myosin": float(order["phi_shared_m"][out_index]),
            "phi_shared_cb2": float(order["phi_shared_a"][out_index]),
            "nx_shared": int(order["nx_shared"][out_index]),
            "registration_pi": float(order["registration_pi"][out_index]),
            "q_sarcomere_axial": float(order["q_sarc"][out_index]),
            "nematic_myosin": float(nem_m[out_index]),
            "nematic_actin": float(nem_a[out_index]),
            "nematic_cb2": float(nem_cb2[out_index]),
            "cb1_fraction": float(cb1[out_index].mean()),
            "cb2_fraction": float(cb2[out_index].mean()),
            "cb_any_fraction": float((cb1[out_index] | cb2[out_index]).mean()),
            "mean_cb2_load": float(mean_load[out_index]),
            "median_cb2_load": float(median_load[out_index]),
            "directional_extra_over_base": float(lam * mean_load[out_index] / BASE_LIFETIME),
            "fraction_cb2_directional_gt_base": float(fraction_directional[out_index]),
            "mean_myosin_binding_ratio": float(mbr[out_index].mean()),
            "fraction_am_engaged": float(np.mean(mbr[out_index] > 0.0)),
            "mean_partial_binding_ratio": float(partial[out_index].mean()),
            "mean_crosslink_ratio": float(crosslink[out_index].mean()),
            "binding_events_per_actin": float(bind_total[out_index]),
            "unbinding_events_per_actin": float(unbind_total[out_index]),
            "binding_rate_per_actin_s": float(bind_rate[out_index]),
            "unbinding_rate_per_actin_s": float(unbind_rate[out_index]),
            "myosin_rms_speed": float(speed_m[out_index]),
            "actin_rms_speed": float(speed_a[out_index]),
        })
    return records


def write_csv(path: Path, fields: tuple[str, ...], records: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as stream:
        return list(csv.DictReader(stream))


def sustained_onset(times: np.ndarray, values: np.ndarray, threshold: float, window: int = 5) -> float:
    if len(values) < window:
        return math.nan
    good = values >= threshold
    sustained = np.convolve(good.astype(int), np.ones(window, dtype=int), mode="valid") == window
    hits = np.flatnonzero(sustained)
    return float(times[hits[0]]) if len(hits) else math.nan


def load_legacy_phi(paths: list[Path] | None) -> dict[str, tuple[float, float]]:
    output = {}
    for path in paths or []:
        if not path.exists():
            continue
        for row in read_csv(path):
            output[row["filename"]] = (float(row["Phi_3D"]), float(row["Theta_deg"]))
    return output


def summarize_file(
    path: Path, records: list[dict[str, str]], legacy: dict[str, tuple[float, float]],
    tail_frames: int,
) -> dict[str, object]:
    first = records[0]
    numeric = {field: np.asarray([float(row[field]) for row in records])
               for field in TIME_FIELDS if field not in {
                   "filename", "lambda_label",
               }}
    tail = slice(max(0, len(records) - tail_frames), len(records))
    phi = numeric["phi_axial_myosin"]
    tail_phi = float(phi[tail].mean())
    early_phi = float(phi[:min(tail_frames, len(phi))].mean())
    if tail_phi >= 0.5:
        classification = "assembled"
    elif float(np.max(phi)) >= 0.5 and tail_phi < 0.4:
        classification = "transient"
    elif float(numeric["cb2_fraction"][tail].mean()) >= 0.9 and tail_phi < 0.5:
        classification = "bonded_disordered"
    elif tail_phi - early_phi >= 0.15:
        classification = "partial_growth"
    else:
        classification = "disordered"
    legacy_phi, legacy_theta = legacy.get(path.name, (math.nan, math.nan))
    with h5py.File(path, "r", locking=False) as handle:
        checkpoint = checkpoint_step(handle)
    lam = float(first["lambda_value"])
    return {
        "filename": path.name,
        "n_actins": int(first["n_actins"]),
        "n_myosins": int(first["n_myosins"]),
        "lambda_label": first["lambda_label"],
        "lambda_value": lam,
        "lambda_over_base": lam / BASE_LIFETIME,
        "replica": int(first["replica"]),
        "checkpoint_step": checkpoint,
        "complete": int(checkpoint >= 300000),
        "n_sampled_frames": len(records),
        "legacy_phi_3d": legacy_phi,
        "legacy_theta_deg": legacy_theta,
        "tail_phi_axial_myosin": tail_phi,
        "tail_phi_axial_cb2": float(numeric["phi_axial_cb2"][tail].mean()),
        "tail_q_sarcomere_axial": float(numeric["q_sarcomere_axial"][tail].mean()),
        "tail_nematic_myosin": float(numeric["nematic_myosin"][tail].mean()),
        "tail_nematic_actin": float(numeric["nematic_actin"][tail].mean()),
        "tail_nematic_cb2": float(numeric["nematic_cb2"][tail].mean()),
        "tail_cb1_fraction": float(numeric["cb1_fraction"][tail].mean()),
        "tail_cb2_fraction": float(numeric["cb2_fraction"][tail].mean()),
        "tail_mean_cb2_load": float(numeric["mean_cb2_load"][tail].mean()),
        "tail_directional_extra_over_base": float(numeric["directional_extra_over_base"][tail].mean()),
        "tail_fraction_cb2_directional_gt_base": float(numeric["fraction_cb2_directional_gt_base"][tail].mean()),
        "tail_am_engaged_fraction": float(numeric["fraction_am_engaged"][tail].mean()),
        "tail_binding_rate_per_actin_s": float(numeric["binding_rate_per_actin_s"][tail].mean()),
        "tail_unbinding_rate_per_actin_s": float(numeric["unbinding_rate_per_actin_s"][tail].mean()),
        "onset_phi_0p4_s": sustained_onset(numeric["time_s"], phi, 0.4),
        "onset_phi_0p5_s": sustained_onset(numeric["time_s"], phi, 0.5),
        "onset_phi_0p6_s": sustained_onset(numeric["time_s"], phi, 0.6),
        "classification": classification,
    }


def process_file(
    path_text: str, cache_dir_text: str, stride: int, tail_frames: int,
    legacy: dict[str, tuple[float, float]], include_incomplete: bool,
) -> tuple[dict[str, object], dict[str, object] | None, str]:
    path = Path(path_text)
    cache_dir = Path(cache_dir_text)
    with h5py.File(path, "r", locking=False) as handle:
        checkpoint = checkpoint_step(handle)
        committed = int(np.asarray(handle["/meta/committed_frames"]).reshape(-1)[-1])
        schema = tuple(sorted(handle.keys()))
    match = PATTERN.match(path.name)
    inventory = {
        "filename": path.name,
        "n_actins": int(match["n_actins"]),
        "lambda_label": match["lambda_label"],
        "lambda_value": lambda_value(match["lambda_label"]),
        "replica": int(match["replica"]),
        "checkpoint_step": checkpoint,
        "committed_frames": committed,
        "complete": int(checkpoint >= 300000),
        "size_bytes": path.stat().st_size,
        "groups": ";".join(schema),
    }
    if checkpoint < 300000 and not include_incomplete:
        return inventory, None, f"skip incomplete {path.name}: {checkpoint}"
    cache_path = cache_dir / f"{path.stem}.csv"
    if cache_path.exists():
        records = read_csv(cache_path)
        message = f"cached {path.name}"
    else:
        records = extract_file(path, stride)
        write_csv(cache_path, TIME_FIELDS, records)
        message = f"extracted {path.name}: {len(records)} frames"
    return inventory, summarize_file(path, records, legacy, tail_frames), message


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-dir", type=Path,
        default=Path("/scratch/users/xsli/sarcomere/data_nalam_s2fstab_base1e-3"),
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--stride", type=int, default=10)
    parser.add_argument("--tail-frames", type=int, default=30)
    parser.add_argument("--legacy-phi-csv", type=Path, action="append")
    parser.add_argument("--include-incomplete", action="store_true")
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--jobs", type=int, default=1)
    args = parser.parse_args()

    cache_dir = args.output_dir / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    legacy = load_legacy_phi(args.legacy_phi_csv)
    inventory = []
    summaries = []
    files = [path for path in sorted(args.input_dir.glob("*.h5")) if PATTERN.match(path.name)]
    if args.limit > 0:
        files = files[:args.limit]
    worker_args = [
        (str(path), str(cache_dir), args.stride, args.tail_frames, legacy,
         args.include_incomplete)
        for path in files
    ]
    if args.jobs <= 1:
        results = [process_file(*item) for item in worker_args]
        for file_index, (item_inventory, summary, message) in enumerate(results, 1):
            inventory.append(item_inventory)
            if summary is not None:
                summaries.append(summary)
            print(f"[{file_index}/{len(files)}] {message}", flush=True)
    else:
        results = []
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            futures = {executor.submit(process_file, *item): item[0] for item in worker_args}
            for file_index, future in enumerate(as_completed(futures), 1):
                result = future.result()
                results.append(result)
                print(f"[{file_index}/{len(files)}] {result[2]}", flush=True)
        for item_inventory, summary, _message in results:
            inventory.append(item_inventory)
            if summary is not None:
                summaries.append(summary)

    inventory.sort(key=lambda row: row["filename"])
    summaries.sort(key=lambda row: row["filename"])

    write_csv(args.output_dir / "inventory.csv", tuple(inventory[0].keys()), inventory)
    write_csv(args.output_dir / "run_summary.csv", SUMMARY_FIELDS, summaries)
    print(f"wrote {len(inventory)} inventory rows and {len(summaries)} summaries", flush=True)


if __name__ == "__main__":
    main()
