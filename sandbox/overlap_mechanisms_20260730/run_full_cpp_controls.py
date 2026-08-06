#!/usr/bin/env python3
"""Short controls in the inherited full 3-D C++ simulator.

These are diagnostic-sized runs only.  They use the existing two-sarcomere
force-balanced initializer and existing default-off register/titin controls;
the production source tree is not modified.
"""

from __future__ import annotations

import argparse
import csv
import os
import shlex
import subprocess
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "results" / "full_cpp"
BASE_CONFIG = Path("/scratch/users/xsli/sarcomere/prodv8_base.txt")

ARMS = {
    "baseline": ["--aa_axial_register_ratio=0", "--titin_k=0"],
    "register": ["--aa_axial_register_ratio=1", "--titin_k=0"],
    "titin": ["--aa_axial_register_ratio=0", "--titin_k=20"],
    "register_plus_titin": ["--aa_axial_register_ratio=1", "--titin_k=20"],
}
INITIALS = {
    "deep": ("3.9", "0.10"),
    "reference": ("4.3", "0.10"),
    "shallow": ("4.7", "0.10"),
}
SEEDS = (7300, 7301, 7302)


def pbc_delta(value: np.ndarray, length: float) -> np.ndarray:
    return value - length * np.round(value / length)


def analyze(path: Path, arm: str, initial: str, seed: int, dt: float) -> dict:
    with h5py.File(path, "r") as handle:
        centers = np.asarray(handle["actin/center"])
        directions = np.asarray(handle["actin/direction"])
        barbed = centers - 0.5 * directions
        box_length = float(path.name.split("__Lx-")[1].split("__")[0])
        y0 = np.abs(pbc_delta(barbed[:, 0, 0] - barbed[:, 1, 0], box_length))
        y1 = np.abs(pbc_delta(barbed[:, 2, 0] - barbed[:, 3, 0], box_length))
        aa = 0.5 * (y0 + y1)
        am = np.mean(np.asarray(handle["actin/am_overlap_available"])[..., 0], axis=1)
        steps = np.asarray(handle["state/current_step"])[:, 0]
        time = steps * dt
        start = int(0.65 * len(time))
        aa_v = float(np.polyfit(time[start:], aa[start:], 1)[0])
        am_v = float(np.polyfit(time[start:], am[start:], 1)[0])
        actin_caps = np.asarray(handle["actin/cumulative_force_cap_events"])[-1, :, 0]
        status2 = np.asarray(handle["actin/aa_local_status2_count"])[start:, :, 0]
    return {
        "arm": arm,
        "initial": initial,
        "seed": seed,
        "aa_tail_um": float(np.mean(aa[start:])),
        "am_tail_um": float(np.mean(am[start:])),
        "aa_tail_velocity_um_s": aa_v,
        "am_tail_velocity_um_s": am_v,
        "status2_links_tail": float(0.5 * np.mean(np.sum(status2, axis=1))),
        "max_actin_force_cap_fraction": float(np.max(actin_caps) / steps[-1]),
        "file": str(path),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--nsteps", type=int, default=60000)
    args = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    base = shlex.split(BASE_CONFIG.read_text())
    dt = 1.0e-6
    rows = []
    for arm, arm_args in ARMS.items():
        for initial, (lx, overlap) in INITIALS.items():
            for seed in SEEDS:
                path = OUTPUT / f"{arm}__{initial}__seed-{seed}__Lx-{lx}__diag.h5"
                if not path.exists():
                    command = [
                        str(args.binary),
                        *base,
                        "--initial_structure=two_sarcomere_unit",
                        "--n_actins=4",
                        "--n_myosins=2",
                        f"--Lx={lx}",
                        "--Ly=0.4",
                        "--Lz=0.4",
                        "--pbc_mask=111",
                        "--aa_spacing_ratio=1.0",
                        f"--init_aa_overlap={overlap}",
                        "--aa_barbed_zone=0.10",
                        "--actin_steric_k=2000",
                        f"--dt={dt}",
                        f"--nsteps={args.nsteps}",
                        "--save_every=300",
                        "--pre_phase_steps=0",
                        "--interaction_ramp_steps=2000",
                        f"--seed={seed}",
                        f"--filename={path}",
                        *arm_args,
                    ]
                    log = path.with_suffix(".log")
                    with log.open("w") as handle:
                        subprocess.run(
                            command,
                            check=True,
                            stdout=handle,
                            stderr=subprocess.STDOUT,
                            env={
                                **os.environ,
                                "OMP_NUM_THREADS": "1",
                                "HDF5_USE_FILE_LOCKING": "FALSE",
                            },
                        )
                rows.append(analyze(path, arm, initial, seed, dt))
                print(f"finished {arm} {initial} seed={seed}", flush=True)

    summary_path = OUTPUT / "summary.csv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    positions = np.arange(len(ARMS))
    for offset, initial in zip((-0.22, 0.0, 0.22), INITIALS):
        aa_means = []
        am_means = []
        aa_std = []
        am_std = []
        for arm in ARMS:
            subset = [r for r in rows if r["arm"] == arm and r["initial"] == initial]
            aa_values = [r["aa_tail_um"] for r in subset]
            am_values = [r["am_tail_um"] for r in subset]
            aa_means.append(np.mean(aa_values))
            am_means.append(np.mean(am_values))
            aa_std.append(np.std(aa_values))
            am_std.append(np.std(am_values))
        axes[0].errorbar(positions + offset, aa_means, yerr=aa_std, fmt="o", label=initial)
        axes[1].errorbar(positions + offset, am_means, yerr=am_std, fmt="o", label=initial)
    for ax in axes:
        ax.set_xticks(positions, ARMS, rotation=20, ha="right")
        ax.legend()
    axes[0].axhline(0.10, color="black", ls=":")
    axes[0].set_ylabel("tail A–A overlap (µm)")
    axes[1].set_ylabel("tail available A–M overlap (µm)")
    fig.suptitle("Full C++ force-balanced diagnostic controls")
    fig.tight_layout()
    fig.savefig(OUTPUT / "full_cpp_controls.png", dpi=180)
    print(f"wrote {summary_path}")


if __name__ == "__main__":
    main()
