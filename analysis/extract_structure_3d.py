"""
Batch processing for 3D structural order.
Extracts 3D Fourier amplitude (Phi_3D), smectic tilt angle (Theta), and the
winning Fourier mode indices for myosins and for cb_status == 2 actins.
"""

import re
import csv
import h5py
import argparse
import numpy as np
from pathlib import Path
from collections import Counter, defaultdict

try:
    import torch
except ImportError:
    torch = None

# =============================================================================
# 1. 3D STRUCTURAL ANALYSIS
# =============================================================================
def _harmonic_values(values, default_start, default_stop, name):
    if values is None:
        return np.arange(default_start, default_stop + 1, dtype=int)
    arr = np.asarray(values, dtype=int)
    if arr.size == 0:
        raise ValueError(f"{name} list must not be empty.")
    if name == "nx" and np.any(arr == 0):
        raise ValueError("nx values must be nonzero; use positive stripe counts.")
    return arr


def build_k_matrix(
    Lx=12.0,
    Ly=2.0,
    Lz=2.0,
    n_max_x=30,
    n_max_y=5,
    n_max_z=5,
    nx_values=None,
    ny_values=None,
    nz_values=None,
):
    # Generate integer harmonics allowed by the PBC box. Shape is (K, 3).
    nx = _harmonic_values(nx_values, 1, n_max_x, "nx")
    ny = _harmonic_values(ny_values, -n_max_y, n_max_y, "ny")
    nz = _harmonic_values(nz_values, -n_max_z, n_max_z, "nz")

    KX, KY, KZ = np.meshgrid(
        2 * np.pi * nx / Lx,
        2 * np.pi * ny / Ly,
        2 * np.pi * nz / Lz,
        indexing="ij"
    )
    return np.stack([KX.ravel(), KY.ravel(), KZ.ravel()], axis=1)


def build_mode_matrix(
    n_max_x=30,
    n_max_y=5,
    n_max_z=5,
    nx_values=None,
    ny_values=None,
    nz_values=None,
):
    # Integer harmonic labels matching build_k_matrix row order.
    nx = _harmonic_values(nx_values, 1, n_max_x, "nx")
    ny = _harmonic_values(ny_values, -n_max_y, n_max_y, "ny")
    nz = _harmonic_values(nz_values, -n_max_z, n_max_z, "nz")

    NX, NY, NZ = np.meshgrid(nx, ny, nz, indexing="ij")
    return np.stack([NX.ravel(), NY.ravel(), NZ.ravel()], axis=1)


def theta_from_k(k_star):
    kx_star, ky_star, kz_star = k_star
    k_perp_mag = np.sqrt(ky_star**2 + kz_star**2)
    return float(np.degrees(np.arctan2(k_perp_mag, np.abs(kx_star))))


def fourier_order_3d(
    x,
    y,
    z,
    Lx=12.0,
    Ly=2.0,
    Lz=2.0,
    n_max_x=30,
    n_max_y=5,
    n_max_z=5,
    nx_values=None,
    ny_values=None,
    nz_values=None,
):
    """
    Calculates the true 3D tilt-invariant structural order and smectic tilt angle
    for one frame. Kept for compatibility; batch extraction uses the vectorized
    implementations below.
    """
    coords = np.column_stack([x, y, z])
    results = fourier_order_3d_batch_numpy(
        coords[None, :, :], Lx=Lx, Ly=Ly, Lz=Lz,
        n_max_x=n_max_x, n_max_y=n_max_y, n_max_z=n_max_z,
        nx_values=nx_values, ny_values=ny_values, nz_values=nz_values,
    )
    return results[0]


def fourier_order_3d_batch_numpy(
    coords_batch,
    Lx=12.0,
    Ly=2.0,
    Lz=2.0,
    n_max_x=30,
    n_max_y=5,
    n_max_z=5,
    nx_values=None,
    ny_values=None,
    nz_values=None,
):
    """Compute Phi_3D and theta for a batch of frames with NumPy."""
    coords_batch = np.asarray(coords_batch, dtype=np.float64)
    if coords_batch.size == 0 or coords_batch.shape[1] == 0:
        return [
            {"phi_3d_max": 0.0, "theta_deg": 0.0, "nx": 0, "ny": 0, "nz": 0}
            for _ in range(coords_batch.shape[0])
        ]

    k_matrix = build_k_matrix(
        Lx, Ly, Lz, n_max_x, n_max_y, n_max_z,
        nx_values=nx_values, ny_values=ny_values, nz_values=nz_values,
    )
    mode_matrix = build_mode_matrix(
        n_max_x, n_max_y, n_max_z,
        nx_values=nx_values, ny_values=ny_values, nz_values=nz_values,
    )
    phases = np.einsum("kc,fnc->fkn", k_matrix, coords_batch, optimize=True)
    real = np.cos(phases).mean(axis=2)
    imag = np.sin(phases).mean(axis=2)
    phi_amp_3d = np.sqrt(real**2 + imag**2)
    idx_star = np.argmax(phi_amp_3d, axis=1)

    return [
        {
            "phi_3d_max": float(phi_amp_3d[i, idx]),
            "theta_deg": theta_from_k(k_matrix[idx]),
            "nx": int(mode_matrix[idx, 0]),
            "ny": int(mode_matrix[idx, 1]),
            "nz": int(mode_matrix[idx, 2]),
        }
        for i, idx in enumerate(idx_star)
    ]


def fourier_order_3d_batch_torch(
    coords_batch,
    device,
    dtype=None,
    Lx=12.0,
    Ly=2.0,
    Lz=2.0,
    n_max_x=30,
    n_max_y=5,
    n_max_z=5,
    nx_values=None,
    ny_values=None,
    nz_values=None,
):
    """Compute Phi_3D and theta for a batch of frames with PyTorch on CPU/GPU."""
    if torch is None:
        raise RuntimeError("PyTorch is not installed; use --device numpy or install torch.")

    if dtype is None:
        dtype = torch.float32 if device.type == "cuda" else torch.float64

    k_matrix_np = build_k_matrix(
        Lx, Ly, Lz, n_max_x, n_max_y, n_max_z,
        nx_values=nx_values, ny_values=ny_values, nz_values=nz_values,
    )
    mode_matrix = build_mode_matrix(
        n_max_x, n_max_y, n_max_z,
        nx_values=nx_values, ny_values=ny_values, nz_values=nz_values,
    )
    if len(coords_batch) == 0 or coords_batch.shape[1] == 0:
        return [
            {"phi_3d_max": 0.0, "theta_deg": 0.0, "nx": 0, "ny": 0, "nz": 0}
            for _ in range(len(coords_batch))
        ]

    coords = torch.as_tensor(coords_batch, dtype=dtype, device=device)
    k_matrix = torch.as_tensor(k_matrix_np, dtype=dtype, device=device)

    with torch.no_grad():
        phases = torch.einsum("kc,fnc->fkn", k_matrix, coords)
        real = torch.cos(phases).mean(dim=2)
        imag = torch.sin(phases).mean(dim=2)
        phi_amp_3d = torch.sqrt(real.square() + imag.square())
        max_phi, idx_star = torch.max(phi_amp_3d, dim=1)

    max_phi_np = max_phi.detach().cpu().numpy()
    idx_star_np = idx_star.detach().cpu().numpy()
    if device.type == "cuda":
        torch.cuda.empty_cache()

    return [
        {
            "phi_3d_max": float(max_phi_np[i]),
            "theta_deg": theta_from_k(k_matrix_np[int(idx)]),
            "nx": int(mode_matrix[int(idx), 0]),
            "ny": int(mode_matrix[int(idx), 1]),
            "nz": int(mode_matrix[int(idx), 2]),
        }
        for i, idx in enumerate(idx_star_np)
    ]


def resolve_compute_device(device):
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


def print_compute_device_info(torch_device):
    if torch_device is None:
        return
    print(f"Using torch/{torch_device} for Fourier computation.")
    if torch_device.type == "cuda":
        idx = torch_device.index if torch_device.index is not None else torch.cuda.current_device()
        props = torch.cuda.get_device_properties(idx)
        total_gb = props.total_memory / 1024**3
        capability = f"{props.major}.{props.minor}"
        print(
            f"GPU {idx}: {props.name}; "
            f"compute capability {capability}; total memory {total_gb:.2f} GiB"
        )


def fourier_order_3d_frames(coords_batch, torch_device=None, nx_values=None, ny_values=None, nz_values=None):
    if torch_device is None:
        return fourier_order_3d_batch_numpy(
            coords_batch,
            nx_values=nx_values,
            ny_values=ny_values,
            nz_values=nz_values,
        )
    return fourier_order_3d_batch_torch(
        coords_batch,
        torch_device,
        nx_values=nx_values,
        ny_values=ny_values,
        nz_values=nz_values,
    )

# =============================================================================
# 2. BATCH PROCESSING PIPELINE
# =============================================================================
CSV_FIELDS = ["Filename", "H5_Path", "A_val", "lt_str", "Phi_3D", "Theta_deg", "nx", "ny", "nz"]
PER_FRAME_CSV_FIELDS = [
    "Filename",
    "H5_Path",
    "A_val",
    "lt_str",
    "Frame",
    "Phi_3D",
    "Theta_deg",
    "nx",
    "ny",
    "nz",
]
ACTIN_CB2_EXTRA_FIELDS = [
    "Phi_3D_plus",
    "Theta_deg_plus",
    "nx_plus",
    "ny_plus",
    "nz_plus",
    "Phi_3D_minus",
    "Theta_deg_minus",
    "nx_minus",
    "ny_minus",
    "nz_minus",
    "n_actin_cb2_plus",
    "n_actin_cb2_minus",
]
ACTIN_CB2_CSV_FIELDS = CSV_FIELDS + ACTIN_CB2_EXTRA_FIELDS
ACTIN_CB2_PER_FRAME_CSV_FIELDS = PER_FRAME_CSV_FIELDS + ACTIN_CB2_EXTRA_FIELDS


def result_key(record):
    return (record["Filename"], str(record["A_val"]), record["lt_str"])


def validate_csv_header(path, required_fields):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return

    with open(path, mode='r', newline='') as csvfile:
        fieldnames = csv.DictReader(csvfile).fieldnames or []

    missing = [field for field in required_fields if field not in fieldnames]
    if missing:
        raise ValueError(
            f"{path} is missing required output columns {missing}. "
            "Write to a new CSV or regenerate this file so H5_Path can be stored."
        )


def load_existing_results(output_csv):
    existing = set()
    output_path = Path(output_csv)
    if not output_path.exists() or output_path.stat().st_size == 0:
        return existing

    with open(output_path, mode='r', newline='') as csvfile:
        for row in csv.DictReader(csvfile):
            if not all(field in row and row[field] for field in ("Filename", "A_val", "lt_str")):
                continue
            existing.add(result_key(row))
    return existing


def default_per_frame_csv(output_csv):
    output_path = Path(output_csv)
    return str(output_path.with_name(f"{output_path.stem}_per_frame{output_path.suffix}"))


def default_actin_cb2_csv(output_csv):
    output_path = Path(output_csv)
    return str(output_path.with_name(f"{output_path.stem}_actin_cb2{output_path.suffix}"))


def default_actin_cb2_per_frame_csv(actin_cb2_csv):
    output_path = Path(actin_cb2_csv)
    return str(output_path.with_name(f"{output_path.stem}_per_frame{output_path.suffix}"))


def most_common_mode(modes):
    if not modes:
        return 0, 0, 0
    return Counter(modes).most_common(1)[0][0]


def zero_order_result():
    return {"phi_3d_max": 0.0, "theta_deg": 0.0, "nx": 0, "ny": 0, "nz": 0}


def compute_single_frame_order(coords, torch_device, nx_values=None, ny_values=None, nz_values=None):
    coords = np.asarray(coords, dtype=np.float64)
    if coords.size == 0:
        return zero_order_result()
    return fourier_order_3d_frames(
        coords[None, :, :],
        torch_device,
        nx_values=nx_values,
        ny_values=ny_values,
        nz_values=nz_values,
    )[0]


def average_population_results(result_plus, result_minus):
    return {
        "phi_3d_max": 0.5 * (result_plus["phi_3d_max"] + result_minus["phi_3d_max"]),
        "theta_deg": 0.5 * (result_plus["theta_deg"] + result_minus["theta_deg"]),
    }


def summarize_numeric(values):
    return float(np.mean(values)) if values else 0.0


def summarize_extra(records, field):
    return summarize_numeric([record[field] for record in records])


def result_mode(result):
    return int(result["nx"]), int(result["ny"]), int(result["nz"])


def nonzero_modes(modes):
    return [mode for mode in modes if mode != (0, 0, 0)]


def extract_structure_data(
    base_dirs,
    output_csv,
    per_frame_csv=None,
    actin_cb2_csv=None,
    actin_cb2_per_frame_csv=None,
    threshold=2880,
    eval_window=300,
    stride=10,
    device="auto",
    batch_size=8,
    nx_values=None,
    ny_values=None,
    nz_values=None,
):
    existing_results = load_existing_results(output_csv)
    actin_cb2_path = Path(actin_cb2_csv or default_actin_cb2_csv(output_csv))
    existing_actin_cb2_results = load_existing_results(actin_cb2_path)
    condition_counters = defaultdict(int)
    discarded_count = 0
    processed_count = 0
    actin_cb2_processed_count = 0
    skipped_count = 0
    actin_cb2_skipped_count = 0
    torch_device = resolve_compute_device(device)
    if torch_device is None:
        print("Using numpy for Fourier computation.")
    else:
        print_compute_device_info(torch_device)

    if batch_size <= 0:
        raise ValueError("batch_size must be positive.")

    output_path = Path(output_csv)
    per_frame_path = Path(per_frame_csv or default_per_frame_csv(output_csv))
    actin_cb2_per_frame_path = Path(
        actin_cb2_per_frame_csv or default_actin_cb2_per_frame_csv(actin_cb2_path)
    )
    validate_csv_header(output_path, CSV_FIELDS)
    validate_csv_header(per_frame_path, PER_FRAME_CSV_FIELDS)
    validate_csv_header(actin_cb2_path, ACTIN_CB2_CSV_FIELDS)
    validate_csv_header(actin_cb2_per_frame_path, ACTIN_CB2_PER_FRAME_CSV_FIELDS)
    write_header = not output_path.exists() or output_path.stat().st_size == 0
    write_per_frame_header = not per_frame_path.exists() or per_frame_path.stat().st_size == 0
    write_actin_cb2_header = not actin_cb2_path.exists() or actin_cb2_path.stat().st_size == 0
    write_actin_cb2_per_frame_header = (
        not actin_cb2_per_frame_path.exists() or actin_cb2_per_frame_path.stat().st_size == 0
    )

    print(f"{'Original/Pseudo-Filename':<55} | {'Frames':<8} | {'Status'}")
    print("-" * 85)

    with open(output_path, mode='a', newline='') as csvfile, \
            open(per_frame_path, mode='a', newline='') as per_frame_csvfile, \
            open(actin_cb2_path, mode='a', newline='') as actin_cb2_csvfile, \
            open(actin_cb2_per_frame_path, mode='a', newline='') as actin_cb2_per_frame_csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=CSV_FIELDS)
        per_frame_writer = csv.DictWriter(per_frame_csvfile, fieldnames=PER_FRAME_CSV_FIELDS)
        actin_cb2_writer = csv.DictWriter(actin_cb2_csvfile, fieldnames=ACTIN_CB2_CSV_FIELDS)
        actin_cb2_per_frame_writer = csv.DictWriter(
            actin_cb2_per_frame_csvfile,
            fieldnames=ACTIN_CB2_PER_FRAME_CSV_FIELDS,
        )
        if write_header:
            writer.writeheader()
            csvfile.flush()
        if write_per_frame_header:
            per_frame_writer.writeheader()
            per_frame_csvfile.flush()
        if write_actin_cb2_header:
            actin_cb2_writer.writeheader()
            actin_cb2_csvfile.flush()
        if write_actin_cb2_per_frame_header:
            actin_cb2_per_frame_writer.writeheader()
            actin_cb2_per_frame_csvfile.flush()

        for d in base_dirs:
            if not d.exists():
                continue

            for h5_file in d.glob("*.h5"):
                name = h5_file.name
                try:
                    with h5py.File(h5_file, 'r') as fh:
                        has_myosin = "/myosin/center" in fh
                        has_actin_cb2 = all(
                            dataset in fh
                            for dataset in ("/actin/center", "/actin/direction", "/actin/cb_status")
                        )
                        if not has_myosin and not has_actin_cb2:
                            continue
                        n_frame_candidates = []
                        if has_myosin:
                            n_frame_candidates.append(fh["/myosin/center"].shape[0])
                        if has_actin_cb2:
                            n_frame_candidates.extend([
                                fh["/actin/center"].shape[0],
                                fh["/actin/direction"].shape[0],
                                fh["/actin/cb_status"].shape[0],
                            ])
                        n_frames = min(n_frame_candidates)

                    if n_frames < threshold:
                        discarded_count += 1
                        continue

                    # Extract A_val and lt_str from filename. Some sweeps include
                    # tags such as "_nodir" between lt and the replicate marker.
                    match = re.search(r'(.*A(\d+)_M\d+_lt([^_]+?)(?:_[^_]+)*)_r\d+(.*)', name)
                    if not match:
                        continue

                    condition_base, A_val, lt_str, extension = match.group(1), int(match.group(2)), match.group(3), match.group(4)
                    condition_counters[condition_base] += 1
                    pseudo_name = f"{condition_base}_r{condition_counters[condition_base]:02d}{extension}"
                    key = (pseudo_name, str(A_val), lt_str)

                    do_myosin = has_myosin and key not in existing_results
                    do_actin_cb2 = has_actin_cb2 and key not in existing_actin_cb2_results

                    if has_myosin and key in existing_results:
                        skipped_count += 1
                    if has_actin_cb2 and key in existing_actin_cb2_results:
                        actin_cb2_skipped_count += 1
                    if not do_myosin and not do_actin_cb2:
                        print(f"{pseudo_name:<55} | {n_frames:<8} | Skipped (in CSV)")
                        continue

                    # Average over the last eval_window frames, sampled every stride frames.
                    start_frame = max(0, n_frames - eval_window)
                    frame_indices = list(range(start_frame, n_frames, stride))

                    if do_myosin:
                        temp_phi3d, temp_theta, temp_modes = [], [], []
                        per_frame_records = []

                        with h5py.File(h5_file, 'r') as fh:
                            centers_ds = fh["/myosin/center"]
                            for batch_start in range(0, len(frame_indices), batch_size):
                                batch_indices = frame_indices[batch_start:batch_start + batch_size]
                                coords_batch = np.asarray(centers_ds[batch_indices])
                                batch_results = fourier_order_3d_frames(
                                    coords_batch,
                                    torch_device,
                                    nx_values=nx_values,
                                    ny_values=ny_values,
                                    nz_values=nz_values,
                                )

                                for idx, res in zip(batch_indices, batch_results):
                                    temp_phi3d.append(res["phi_3d_max"])
                                    temp_theta.append(res["theta_deg"])
                                    mode = result_mode(res)
                                    temp_modes.append(mode)
                                    per_frame_records.append({
                                        "Filename": pseudo_name,
                                        "H5_Path": str(h5_file.resolve()),
                                        "A_val": A_val,
                                        "lt_str": lt_str,
                                        "Frame": idx,
                                        "Phi_3D": res["phi_3d_max"],
                                        "Theta_deg": res["theta_deg"],
                                        "nx": mode[0],
                                        "ny": mode[1],
                                        "nz": mode[2],
                                    })

                        nx_mode, ny_mode, nz_mode = most_common_mode(temp_modes)
                        record = {
                            "Filename": pseudo_name,
                            "H5_Path": str(h5_file.resolve()),
                            "A_val": A_val,
                            "lt_str": lt_str,
                            "Phi_3D": np.mean(temp_phi3d),
                            "Theta_deg": np.mean(temp_theta),
                            "nx": nx_mode,
                            "ny": ny_mode,
                            "nz": nz_mode,
                        }
                        writer.writerow(record)
                        per_frame_writer.writerows(per_frame_records)
                        csvfile.flush()
                        per_frame_csvfile.flush()
                        existing_results.add(result_key(record))
                        processed_count += 1

                    if do_actin_cb2:
                        temp_phi3d, temp_theta, temp_modes = [], [], []
                        plus_modes, minus_modes = [], []
                        actin_cb2_per_frame_records = []

                        with h5py.File(h5_file, 'r') as fh:
                            actin_centers_ds = fh["/actin/center"]
                            actin_dirs_ds = fh["/actin/direction"]
                            cb_status_ds = fh["/actin/cb_status"]
                            dir_0 = np.asarray(actin_dirs_ds[0], dtype=float)
                            plus_end_mask = dir_0[:, 0] > 0
                            minus_end_mask = dir_0[:, 0] < 0

                            for idx in frame_indices:
                                centers = np.asarray(actin_centers_ds[idx], dtype=float)
                                cb2_mask = np.asarray(cb_status_ds[idx]).reshape(-1) == 2
                                coords_plus = centers[plus_end_mask & cb2_mask]
                                coords_minus = centers[minus_end_mask & cb2_mask]

                                res_plus = compute_single_frame_order(
                                    coords_plus,
                                    torch_device,
                                    nx_values=nx_values,
                                    ny_values=ny_values,
                                    nz_values=nz_values,
                                )
                                res_minus = compute_single_frame_order(
                                    coords_minus,
                                    torch_device,
                                    nx_values=nx_values,
                                    ny_values=ny_values,
                                    nz_values=nz_values,
                                )
                                res_avg = average_population_results(res_plus, res_minus)
                                mode_plus = result_mode(res_plus)
                                mode_minus = result_mode(res_minus)
                                frame_mode = most_common_mode(nonzero_modes([mode_plus, mode_minus]))

                                temp_phi3d.append(res_avg["phi_3d_max"])
                                temp_theta.append(res_avg["theta_deg"])
                                temp_modes.extend(nonzero_modes([mode_plus, mode_minus]))
                                plus_modes.append(mode_plus)
                                minus_modes.append(mode_minus)

                                actin_cb2_per_frame_records.append({
                                    "Filename": pseudo_name,
                                    "H5_Path": str(h5_file.resolve()),
                                    "A_val": A_val,
                                    "lt_str": lt_str,
                                    "Frame": idx,
                                    "Phi_3D": res_avg["phi_3d_max"],
                                    "Theta_deg": res_avg["theta_deg"],
                                    "nx": frame_mode[0],
                                    "ny": frame_mode[1],
                                    "nz": frame_mode[2],
                                    "Phi_3D_plus": res_plus["phi_3d_max"],
                                    "Theta_deg_plus": res_plus["theta_deg"],
                                    "nx_plus": mode_plus[0],
                                    "ny_plus": mode_plus[1],
                                    "nz_plus": mode_plus[2],
                                    "Phi_3D_minus": res_minus["phi_3d_max"],
                                    "Theta_deg_minus": res_minus["theta_deg"],
                                    "nx_minus": mode_minus[0],
                                    "ny_minus": mode_minus[1],
                                    "nz_minus": mode_minus[2],
                                    "n_actin_cb2_plus": int(coords_plus.shape[0]),
                                    "n_actin_cb2_minus": int(coords_minus.shape[0]),
                                })

                        nx_mode, ny_mode, nz_mode = most_common_mode(temp_modes)
                        nx_plus, ny_plus, nz_plus = most_common_mode(nonzero_modes(plus_modes))
                        nx_minus, ny_minus, nz_minus = most_common_mode(nonzero_modes(minus_modes))
                        actin_cb2_record = {
                            "Filename": pseudo_name,
                            "H5_Path": str(h5_file.resolve()),
                            "A_val": A_val,
                            "lt_str": lt_str,
                            "Phi_3D": np.mean(temp_phi3d),
                            "Theta_deg": np.mean(temp_theta),
                            "nx": nx_mode,
                            "ny": ny_mode,
                            "nz": nz_mode,
                            "Phi_3D_plus": summarize_extra(actin_cb2_per_frame_records, "Phi_3D_plus"),
                            "Theta_deg_plus": summarize_extra(actin_cb2_per_frame_records, "Theta_deg_plus"),
                            "nx_plus": nx_plus,
                            "ny_plus": ny_plus,
                            "nz_plus": nz_plus,
                            "Phi_3D_minus": summarize_extra(actin_cb2_per_frame_records, "Phi_3D_minus"),
                            "Theta_deg_minus": summarize_extra(actin_cb2_per_frame_records, "Theta_deg_minus"),
                            "nx_minus": nx_minus,
                            "ny_minus": ny_minus,
                            "nz_minus": nz_minus,
                            "n_actin_cb2_plus": summarize_extra(actin_cb2_per_frame_records, "n_actin_cb2_plus"),
                            "n_actin_cb2_minus": summarize_extra(actin_cb2_per_frame_records, "n_actin_cb2_minus"),
                        }
                        actin_cb2_writer.writerow(actin_cb2_record)
                        actin_cb2_per_frame_writer.writerows(actin_cb2_per_frame_records)
                        actin_cb2_csvfile.flush()
                        actin_cb2_per_frame_csvfile.flush()
                        existing_actin_cb2_results.add(result_key(actin_cb2_record))
                        actin_cb2_processed_count += 1

                    status_parts = []
                    if do_myosin:
                        status_parts.append("myosin")
                    if do_actin_cb2:
                        status_parts.append("actin_cb2")
                    print(f"{pseudo_name:<55} | {n_frames:<8} | Processed {', '.join(status_parts)}")

                except Exception as e:
                    print(f"{name:<55} | {'ERROR':<8} | {e}")

    print("-" * 85)
    print(
        f"Extraction Complete. Added {processed_count} runs to {output_csv}; "
        f"skipped {skipped_count} already-computed runs; discarded {discarded_count} short runs."
    )
    print(f"Per-frame results for processed runs saved to {per_frame_path}.")
    print(
        f"Actin cb_status==2 extraction added {actin_cb2_processed_count} runs to {actin_cb2_path}; "
        f"skipped {actin_cb2_skipped_count} already-computed runs."
    )
    print(f"Actin cb_status==2 per-frame results saved to {actin_cb2_per_frame_path}.")

# =============================================================================
# 3. MAIN EXECUTION
# =============================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract 3D structural metrics from HDF5.")
    parser.add_argument("--out_csv", type=str, default="structural_phase_data.csv")
    parser.add_argument("--per_frame_csv", type=str, default=None,
                        help="CSV for sampled per-frame metrics; defaults to <out_csv stem>_per_frame.csv")
    parser.add_argument("--actin_cb2_csv", type=str, default=None,
                        help="CSV for cb_status==2 actin metrics; defaults to <out_csv stem>_actin_cb2.csv")
    parser.add_argument("--actin_cb2_per_frame_csv", type=str, default=None,
                        help="CSV for sampled per-frame cb_status==2 actin metrics; defaults to <actin_cb2_csv stem>_per_frame.csv")
    parser.add_argument("--eval_window", type=int, default=300,
                        help="Number of trailing frames to average over")
    parser.add_argument("--stride", type=int, default=10,
                        help="Frame stride within the trailing evaluation window")
    parser.add_argument("--device", choices=("auto", "cuda", "cpu", "numpy"), default="auto",
                        help="Compute backend. auto uses CUDA if PyTorch can see a GPU, otherwise NumPy.")
    parser.add_argument("--batch_size", type=int, default=128,
                        help="Number of sampled frames to process per vectorized batch")
    parser.add_argument("--nx-list", type=int, nargs="+", default=None,
                        help="Explicit nx harmonics to scan, e.g. --nx-list 1 2 3 4. Defaults to 1..30.")
    parser.add_argument("--ny-list", type=int, nargs="+", default=None,
                        help="Explicit ny harmonics to scan, e.g. --ny-list -1 0 1. Defaults to -5..5.")
    parser.add_argument("--nz-list", type=int, nargs="+", default=None,
                        help="Explicit nz harmonics to scan, e.g. --nz-list -1 0 1. Defaults to -5..5.")
    args = parser.parse_args()

    base_dirs = [
        #Path("/scratch/users/xsli/sarcomere/2d_sweep_finished"),
        #Path("/scratch/users/xsli/sarcomere/data0416"),
        #Path("/scratch/users/xsli/sarcomere/data0430")
        Path("/scratch/users/xsli/sarcomere/data0423")

    ]

    extract_structure_data(
        base_dirs,
        args.out_csv,
        per_frame_csv=args.per_frame_csv,
        actin_cb2_csv=args.actin_cb2_csv,
        actin_cb2_per_frame_csv=args.actin_cb2_per_frame_csv,
        eval_window=args.eval_window,
        stride=args.stride,
        device=args.device,
        batch_size=args.batch_size,
        nx_values=args.nx_list,
        ny_values=args.ny_list,
        nz_values=args.nz_list,
    )
