"""Extract axial/3D band order, transverse coherence, and M/Z registration."""

import argparse
import csv
import re
from pathlib import Path

import h5py
import numpy as np


FIELDS = [
    "filename", "h5_path", "lambda", "frame", "step", "nx",
    "n_myosin", "n_actin_cb2",
    "phi_x_myosin", "phi_x_actin_cb2",
    "phi_x_myosin_corrected", "phi_x_actin_cb2_corrected",
    "phi_3d_myosin", "theta_3d_myosin_deg",
    "nx_3d_myosin", "ny_3d_myosin", "nz_3d_myosin",
    "phi_3d_actin_cb2", "theta_3d_actin_cb2_deg",
    "nx_3d_actin_cb2", "ny_3d_actin_cb2", "nz_3d_actin_cb2",
    "shared_nx_3d", "shared_ny_3d", "shared_nz_3d", "shared_theta_3d_deg",
    "shared_phi_3d_myosin", "shared_phi_3d_actin_cb2",
    "shared_phi_3d_myosin_corrected", "shared_phi_3d_actin_cb2_corrected",
    "shared_local_amp_myosin", "shared_local_amp_actin_cb2",
    "shared_coherence_yz_myosin", "shared_coherence_yz_actin_cb2",
    "shared_phase_delta_m_minus_z", "shared_registration_pi",
    "shared_joint_band_amplitude",
    "local_amp_myosin", "local_amp_actin_cb2",
    "coherence_yz_myosin", "coherence_yz_actin_cb2",
    "phase_myosin", "phase_actin_cb2", "phase_delta_m_minus_z",
    "registration_pi", "joint_band_amplitude",
    "q_sarcomere_axial", "q_sarcomere",
]


def build_3d_modes(lx, ly, lz, n_max_x, n_max_perp):
    """Return PBC-compatible wavevectors and their integer mode labels."""
    nx = np.arange(1, n_max_x + 1, dtype=int)
    n_perp = np.arange(-n_max_perp, n_max_perp + 1, dtype=int)
    nx_grid, ny_grid, nz_grid = np.meshgrid(nx, n_perp, n_perp, indexing="ij")
    mode_matrix = np.column_stack(
        (nx_grid.ravel(), ny_grid.ravel(), nz_grid.ravel())
    )
    k_matrix = 2.0 * np.pi * mode_matrix / np.array([lx, ly, lz])
    return k_matrix, mode_matrix


def fourier_order_3d(coords, k_matrix, mode_matrix):
    """Find the strongest 3D Fourier mode for one particle population."""
    if len(coords) == 0:
        return {
            "phi": 0.0, "theta_deg": 0.0,
            "nx": 0, "ny": 0, "nz": 0,
        }

    phases = np.asarray(coords, dtype=float) @ k_matrix.T
    order = np.exp(1j * phases).mean(axis=0)
    index = int(np.argmax(np.abs(order)))
    k_star = k_matrix[index]
    mode_star = mode_matrix[index]
    theta = np.degrees(
        np.arctan2(np.hypot(k_star[1], k_star[2]), abs(k_star[0]))
    )
    return {
        "phi": float(abs(order[index])),
        "theta_deg": float(theta),
        "nx": int(mode_star[0]),
        "ny": int(mode_star[1]),
        "nz": int(mode_star[2]),
    }


def complex_orders_3d(coords, k_matrix):
    if len(coords) == 0:
        return np.zeros(len(k_matrix), dtype=np.complex128)
    phases = np.asarray(coords, dtype=float) @ k_matrix.T
    return np.exp(1j * phases).mean(axis=0)


def complex_order(coords, kx):
    if len(coords) == 0:
        return 0j
    return np.exp(1j * kx * coords[:, 0]).mean()


def null_corrected(phi, n_particles):
    if n_particles <= 0:
        return 0.0
    null = np.sqrt(np.pi) / (2.0 * np.sqrt(n_particles))
    return max(0.0, (phi - null) / (1.0 - null))


def null_corrected_array(phi, n_particles):
    if n_particles <= 0:
        return np.zeros_like(phi, dtype=float)
    null = np.sqrt(np.pi) / (2.0 * np.sqrt(n_particles))
    return np.maximum(0.0, (np.asarray(phi) - null) / (1.0 - null))


def transverse_metrics(coords, kx, ly, lz, ny_bins, nz_bins):
    if len(coords) == 0:
        return 0.0, 0.0
    iy = np.floor((coords[:, 1] + ly / 2) * ny_bins / ly).astype(int) % ny_bins
    iz = np.floor((coords[:, 2] + lz / 2) * nz_bins / lz).astype(int) % nz_bins
    cell = iy * nz_bins + iz
    n_cells = ny_bins * nz_bins
    phase = np.exp(1j * kx * coords[:, 0])
    counts = np.bincount(cell, minlength=n_cells)
    sums_r = np.bincount(cell, weights=phase.real, minlength=n_cells)
    sums_i = np.bincount(cell, weights=phase.imag, minlength=n_cells)
    valid = counts > 0
    local = np.zeros(n_cells, dtype=np.complex128)
    local[valid] = (sums_r[valid] + 1j * sums_i[valid]) / counts[valid]
    local_amp = float(np.sum(counts[valid] * np.abs(local[valid])) / len(coords))
    global_amp = abs(phase.mean())
    coherence = float(global_amp / local_amp) if local_amp > 0 else 0.0
    return local_amp, min(coherence, 1.0)


def transverse_metrics_3d(coords, k, ly, lz, ny_bins, nz_bins):
    if len(coords) == 0:
        return 0.0, 0.0
    iy = np.floor((coords[:, 1] + ly / 2) * ny_bins / ly).astype(int) % ny_bins
    iz = np.floor((coords[:, 2] + lz / 2) * nz_bins / lz).astype(int) % nz_bins
    cell = iy * nz_bins + iz
    n_cells = ny_bins * nz_bins
    phase = np.exp(1j * (np.asarray(coords, dtype=float) @ k))
    counts = np.bincount(cell, minlength=n_cells)
    sums_r = np.bincount(cell, weights=phase.real, minlength=n_cells)
    sums_i = np.bincount(cell, weights=phase.imag, minlength=n_cells)
    occupied = counts > 0
    local = np.zeros(n_cells, dtype=np.complex128)
    local[occupied] = (
        sums_r[occupied] + 1j * sums_i[occupied]
    ) / counts[occupied]
    local_amp = float(
        np.sum(counts[occupied] * np.abs(local[occupied])) / len(coords)
    )
    global_amp = abs(phase.mean())
    coherence = float(global_amp / local_amp) if local_amp > 0 else 0.0
    return local_amp, min(coherence, 1.0)


def shared_sarcomeric_order(
    myosin, actin_cb2, k_matrix, mode_matrix, ly, lz, ny_bins, nz_bins,
):
    """Evaluate both populations at one jointly selected 3D wavevector."""
    if len(myosin) == 0 or len(actin_cb2) == 0:
        return {
            "nx": 0, "ny": 0, "nz": 0, "theta_deg": 0.0,
            "phi_m": 0.0, "phi_z": 0.0,
            "phi_m_corrected": 0.0, "phi_z_corrected": 0.0,
            "local_amp_m": 0.0, "local_amp_z": 0.0,
            "coherence_m": 0.0, "coherence_z": 0.0,
            "delta": 0.0, "registration_pi": 0.0,
            "joint": 0.0, "q_sarcomere": 0.0,
        }

    psi_m = complex_orders_3d(myosin, k_matrix)
    psi_z = complex_orders_3d(actin_cb2, k_matrix)
    phi_m = np.abs(psi_m)
    phi_z = np.abs(psi_z)
    corrected_m = null_corrected_array(phi_m, len(myosin))
    corrected_z = null_corrected_array(phi_z, len(actin_cb2))
    selection_score = np.sqrt(corrected_m * corrected_z)

    if np.max(selection_score) > 0:
        index = int(np.argmax(selection_score))
    else:
        index = int(np.argmax(np.sqrt(phi_m * phi_z)))

    k_star = k_matrix[index]
    mode_star = mode_matrix[index]
    amp_m, coherence_m = transverse_metrics_3d(
        myosin, k_star, ly, lz, ny_bins, nz_bins
    )
    amp_z, coherence_z = transverse_metrics_3d(
        actin_cb2, k_star, ly, lz, ny_bins, nz_bins
    )
    delta = float(np.angle(psi_m[index] * np.conj(psi_z[index])))
    registration_pi = float(-np.cos(delta))
    registration_factor = 0.5 * (1.0 + registration_pi)
    joint = float(selection_score[index])
    coherence = float(np.sqrt(coherence_m * coherence_z))
    theta = np.degrees(
        np.arctan2(np.hypot(k_star[1], k_star[2]), abs(k_star[0]))
    )
    return {
        "nx": int(mode_star[0]),
        "ny": int(mode_star[1]),
        "nz": int(mode_star[2]),
        "theta_deg": float(theta),
        "phi_m": float(phi_m[index]),
        "phi_z": float(phi_z[index]),
        "phi_m_corrected": float(corrected_m[index]),
        "phi_z_corrected": float(corrected_z[index]),
        "local_amp_m": amp_m,
        "local_amp_z": amp_z,
        "coherence_m": coherence_m,
        "coherence_z": coherence_z,
        "delta": delta,
        "registration_pi": registration_pi,
        "joint": joint,
        "q_sarcomere": joint * coherence * registration_factor,
    }


def analyze_file(
    path, stride, nx, lx, ly, lz, ny_bins, nz_bins, save_every,
    phi3d_n_max_x, phi3d_n_max_perp, shared_all_nx=False,
):
    match = re.search(r"_lt([^_]+)_", path.name)
    lambda_label = match.group(1) if match else "unknown"
    kx = 2 * np.pi * nx / lx
    k_matrix, mode_matrix = build_3d_modes(
        lx, ly, lz, phi3d_n_max_x, phi3d_n_max_perp
    )
    shared_n_max_x = phi3d_n_max_x if shared_all_nx else nx
    shared_k_matrix, shared_mode_matrix = build_3d_modes(
        lx, ly, lz, shared_n_max_x, phi3d_n_max_perp
    )
    if not shared_all_nx:
        shared_mask = shared_mode_matrix[:, 0] == nx
        shared_k_matrix = shared_k_matrix[shared_mask]
        shared_mode_matrix = shared_mode_matrix[shared_mask]
    records = []
    with h5py.File(path, "r") as handle:
        required = ("/myosin/center", "/actin/center", "/actin/cb_status")
        if not all(key in handle for key in required):
            return records
        n_frames = min(handle[key].shape[0] for key in required)
        for frame in range(0, n_frames, stride):
            myosin = np.asarray(handle["/myosin/center"][frame], dtype=float)
            actin = np.asarray(handle["/actin/center"][frame], dtype=float)
            cb2 = np.asarray(handle["/actin/cb_status"][frame]).reshape(-1) == 2
            actin_cb2 = actin[cb2]

            psi_m = complex_order(myosin, kx)
            psi_z = complex_order(actin_cb2, kx)
            phi_m, phi_z = abs(psi_m), abs(psi_z)
            phi3d_m = fourier_order_3d(myosin, k_matrix, mode_matrix)
            phi3d_z = fourier_order_3d(actin_cb2, k_matrix, mode_matrix)
            shared = shared_sarcomeric_order(
                myosin, actin_cb2, shared_k_matrix, shared_mode_matrix,
                ly, lz, ny_bins, nz_bins,
            )
            amp_m, coh_m = transverse_metrics(myosin, kx, ly, lz, ny_bins, nz_bins)
            amp_z, coh_z = transverse_metrics(actin_cb2, kx, ly, lz, ny_bins, nz_bins)
            delta = float(np.angle(psi_m * np.conj(psi_z))) if len(actin_cb2) else 0.0
            registration_pi = float(-np.cos(delta)) if len(actin_cb2) else 0.0
            joint = float(np.sqrt(phi_m * phi_z))
            registration_factor = max(0.0, 0.5 * (1.0 + registration_pi))
            q_sarc_axial = joint * np.sqrt(coh_m * coh_z) * registration_factor

            records.append({
                "filename": path.name, "h5_path": str(path.resolve()),
                "lambda": lambda_label, "frame": frame, "step": frame * save_every,
                "nx": nx, "n_myosin": len(myosin), "n_actin_cb2": len(actin_cb2),
                "phi_x_myosin": phi_m, "phi_x_actin_cb2": phi_z,
                "phi_x_myosin_corrected": null_corrected(phi_m, len(myosin)),
                "phi_x_actin_cb2_corrected": null_corrected(phi_z, len(actin_cb2)),
                "phi_3d_myosin": phi3d_m["phi"],
                "theta_3d_myosin_deg": phi3d_m["theta_deg"],
                "nx_3d_myosin": phi3d_m["nx"],
                "ny_3d_myosin": phi3d_m["ny"],
                "nz_3d_myosin": phi3d_m["nz"],
                "phi_3d_actin_cb2": phi3d_z["phi"],
                "theta_3d_actin_cb2_deg": phi3d_z["theta_deg"],
                "nx_3d_actin_cb2": phi3d_z["nx"],
                "ny_3d_actin_cb2": phi3d_z["ny"],
                "nz_3d_actin_cb2": phi3d_z["nz"],
                "shared_nx_3d": shared["nx"],
                "shared_ny_3d": shared["ny"],
                "shared_nz_3d": shared["nz"],
                "shared_theta_3d_deg": shared["theta_deg"],
                "shared_phi_3d_myosin": shared["phi_m"],
                "shared_phi_3d_actin_cb2": shared["phi_z"],
                "shared_phi_3d_myosin_corrected": shared["phi_m_corrected"],
                "shared_phi_3d_actin_cb2_corrected": shared["phi_z_corrected"],
                "shared_local_amp_myosin": shared["local_amp_m"],
                "shared_local_amp_actin_cb2": shared["local_amp_z"],
                "shared_coherence_yz_myosin": shared["coherence_m"],
                "shared_coherence_yz_actin_cb2": shared["coherence_z"],
                "shared_phase_delta_m_minus_z": shared["delta"],
                "shared_registration_pi": shared["registration_pi"],
                "shared_joint_band_amplitude": shared["joint"],
                "local_amp_myosin": amp_m, "local_amp_actin_cb2": amp_z,
                "coherence_yz_myosin": coh_m, "coherence_yz_actin_cb2": coh_z,
                "phase_myosin": float(np.angle(psi_m)),
                "phase_actin_cb2": float(np.angle(psi_z)) if len(actin_cb2) else 0.0,
                "phase_delta_m_minus_z": delta, "registration_pi": registration_pi,
                "joint_band_amplitude": joint,
                "q_sarcomere_axial": q_sarc_axial,
                "q_sarcomere": shared["q_sarcomere"],
            })
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--stride", type=int, default=10)
    parser.add_argument(
        "--nx", type=int, default=5,
        help="Expected axial harmonic; fixed while selecting the shared tilted mode.",
    )
    parser.add_argument("--box", type=float, nargs=3, default=(12.0, 4.0, 4.0))
    parser.add_argument("--yz-bins", type=int, nargs=2, default=(4, 4))
    parser.add_argument(
        "--phi3d-n-max-x", type=int, default=30,
        help="Largest positive axial harmonic included in the 3D mode search.",
    )
    parser.add_argument(
        "--phi3d-n-max-perp", type=int, default=5,
        help="Largest absolute y/z harmonic included in the 3D mode search.",
    )
    parser.add_argument(
        "--shared-all-nx", action="store_true",
        help="Jointly select the M/Z shared mode over all axial harmonics instead "
             "of restricting it to --nx (recommended for history protocols).",
    )
    parser.add_argument("--save-every", type=int, default=200)
    args = parser.parse_args()

    if args.nx < 1:
        parser.error("--nx must be at least 1")
    if args.phi3d_n_max_x < 1:
        parser.error("--phi3d-n-max-x must be at least 1")
    if args.phi3d_n_max_perp < 0:
        parser.error("--phi3d-n-max-perp must be nonnegative")

    files = []
    for item in args.inputs:
        files.extend(sorted(item.glob("*.h5")) if item.is_dir() else [item])
    all_records = []
    for path in files:
        records = analyze_file(
            path, args.stride, args.nx, *args.box, *args.yz_bins, args.save_every,
            args.phi3d_n_max_x, args.phi3d_n_max_perp, args.shared_all_nx,
        )
        all_records.extend(records)
        print(f"{path.name}: {len(records)} sampled frames", flush=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(all_records)
    print(f"wrote {len(all_records)} rows to {args.output}", flush=True)


if __name__ == "__main__":
    main()
