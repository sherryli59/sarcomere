#!/usr/bin/env python3
"""Plot structural and ACTN2-state trends for one force-projection trajectory."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np

from extract_structure_3d import fourier_order_3d_batch_numpy


def nematic_order(directions: np.ndarray) -> float:
    norms = np.linalg.norm(directions, axis=1, keepdims=True)
    unit = directions / np.where(norms > 0.0, norms, 1.0)
    q_tensor = 1.5 * (unit.T @ unit) / len(unit) - 0.5 * np.eye(3)
    return float(np.linalg.eigvalsh(q_tensor)[-1])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("trajectory", type=Path)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--dt", type=float, default=1e-5)
    parser.add_argument("--box", type=float, nargs=3, default=(12.0, 2.0, 2.0))
    parser.add_argument(
        "--frame-stride", type=int, default=1,
        help="Analyze every Nth committed frame, always including the latest frame",
    )
    args = parser.parse_args()
    if args.frame_stride < 1:
        parser.error("--frame-stride must be at least 1")

    with h5py.File(args.trajectory, "r", locking=False) as handle:
        committed = int(np.asarray(handle["/meta/committed_frames"]).reshape(-1)[0])
        frame_indices = np.arange(0, committed, args.frame_stride, dtype=int)
        if committed and (frame_indices.size == 0 or frame_indices[-1] != committed - 1):
            frame_indices = np.append(frame_indices, committed - 1)
        steps = np.asarray(handle["/state/current_step"][frame_indices]).reshape(-1)
        myosin_center = np.asarray(handle["/myosin/center"][frame_indices])
        myosin_direction = np.asarray(handle["/myosin/direction"][frame_indices])
        actin_direction = np.asarray(handle["/actin/direction"][frame_indices])
        pair_counts = np.asarray(
            handle["/state/aa_pairs_current_count"][frame_indices]
        ).reshape(-1)
        pair_status = handle["/state/aa_status_current"]

        n_selected = len(frame_indices)
        status1 = np.empty(n_selected, dtype=int)
        status2 = np.empty(n_selected, dtype=int)
        for frame, (source_frame, count) in enumerate(zip(frame_indices, pair_counts)):
            values = np.asarray(pair_status[source_frame, : int(count), 0])
            status1[frame] = int(np.count_nonzero(values == 1))
            status2[frame] = int(np.count_nonzero(values == 2))

    lx, ly, lz = args.box
    phi_3d = np.empty(n_selected, dtype=float)
    for frame in range(n_selected):
        # This is the established/legacy all-myosin Phi_3D definition: maximum
        # Fourier amplitude over nx=1..30 and ny,nz=-5..5.
        result = fourier_order_3d_batch_numpy(
            myosin_center[frame : frame + 1], Lx=lx, Ly=ly, Lz=lz
        )[0]
        phi_3d[frame] = result["phi_3d_max"]

    actin_nematic = np.asarray([nematic_order(frame) for frame in actin_direction])
    myosin_nematic = np.asarray([nematic_order(frame) for frame in myosin_direction])
    time = steps * args.dt

    prefix = args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    csv_path = Path(f"{prefix}_timeseries.csv")
    with csv_path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ["frame", "step", "time", "phi_3d_myosin", "actin_nematic_order",
             "myosin_nematic_order",
             "n_status1_bonds", "n_status2_bonds", "n_aa_bonds"]
        )
        for frame in range(n_selected):
            writer.writerow(
                [frame, int(steps[frame]), time[frame], phi_3d[frame],
                 actin_nematic[frame], myosin_nematic[frame],
                 status1[frame], status2[frame], int(pair_counts[frame])]
            )

    fig, axes = plt.subplots(3, 1, figsize=(8.0, 8.4), sharex=True)
    axes[0].plot(time, phi_3d, color="#3b528b", lw=2)
    axes[0].set_ylabel(r"Myosin $\Phi_{3D}$")
    axes[0].set_ylim(0.0, 1.0)
    axes[1].plot(time, actin_nematic, color="#21918c", lw=2, label="actin")
    axes[1].plot(time, myosin_nematic, color="#7e57c2", lw=2, label="myosin")
    axes[1].set_ylabel("Nematic order")
    axes[1].set_ylim(0.0, 1.0)
    axes[1].legend(frameon=False, ncol=2)
    axes[2].plot(time, status1, color="#e69f00", lw=2, label="status 1")
    axes[2].plot(time, status2, color="#0072b2", lw=2, label="status 2")
    axes[2].set_ylabel("ACTN2 bonds")
    axes[2].set_xlabel("Simulation time")
    axes[2].legend(frameon=False, ncol=2)
    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(alpha=0.2)
    fig.tight_layout()
    png_path = Path(f"{prefix}_trends.png")
    fig.savefig(png_path, dpi=220, bbox_inches="tight")
    plt.close(fig)

    print(
        f"frames={n_selected}/{committed} stride={args.frame_stride} "
        f"last_step={int(steps[-1])}"
    )
    print(f"csv={csv_path}")
    print(f"plot={png_path}")
    print(
        "final: "
        f"phi_3d={phi_3d[-1]:.6f} actin_nematic={actin_nematic[-1]:.6f} "
        f"myosin_nematic={myosin_nematic[-1]:.6f} "
        f"status1={status1[-1]} status2={status2[-1]}"
    )


if __name__ == "__main__":
    main()
