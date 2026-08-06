#!/usr/bin/env python3
"""Structural, numerical, and KMC audit for sparse-ACTN2 assembly pilots."""

from __future__ import annotations

import argparse
import csv
import os
import re
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import h5py
import matplotlib.pyplot as plt
import numpy as np

from analyze_local_aa_coupling_controls import local_state_exposure
from extract_structure_3d import fourier_order_3d_frames


PATTERN = re.compile(
    r"A(?P<na>\d+)_M(?P<nm>\d+)_(?P<condition>aa_off|one_nocatch|two_nocatch|two_status2|two_status1)"
    r"_r(?P<rep>\d+)_seed(?P<seed>\d+)_L.*\.h5$"
)


def nematic(directions: np.ndarray) -> np.ndarray:
    eye = np.eye(3)
    values = []
    for frame in directions:
        q = 1.5 * np.einsum("ni,nj->ij", frame, frame) / len(frame) - 0.5 * eye
        values.append(np.linalg.eigvalsh(q)[-1])
    return np.asarray(values)


def checkpoint(handle: h5py.File) -> int:
    values = []
    for key in ("/resume/current_step", "/resume_prev/current_step", "/state/current_step"):
        if key in handle:
            values.append(int(np.asarray(handle[key]).reshape(-1)[-1]))
    return max(values, default=0)


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--dt", type=float, default=1e-6)
    parser.add_argument("--pre-steps", type=int, default=50000)
    parser.add_argument("--phi-stride", type=int, default=10)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    per_run = []
    series: defaultdict[str, list[dict]] = defaultdict(list)
    for path in sorted(args.input_dir.glob("*.h5")):
        match = PATTERN.fullmatch(path.name)
        if not match:
            continue
        with h5py.File(path, "r", locking=False) as h:
            step = checkpoint(h)
            ac = np.asarray(h["/actin/center"], dtype=float)
            ad = np.asarray(h["/actin/direction"], dtype=float)
            mc = np.asarray(h["/myosin/center"], dtype=float)
            md = np.asarray(h["/myosin/direction"], dtype=float)
            finite = all(np.isfinite(x).all() for x in (ac, ad, mc, md))
            an, mn = nematic(ad), nematic(md)
            selected = np.arange(0, len(mc), max(1, args.phi_stride))
            phi = np.asarray([
                row["phi_3d_max"]
                for row in fourier_order_3d_frames(mc[selected], None)
            ])
            status1 = np.asarray(h["/actin/aa_local_status1_count"], dtype=float).reshape(len(ac), -1).sum(1) / 2
            status2 = np.asarray(h["/actin/aa_local_status2_count"], dtype=float).reshape(len(ac), -1).sum(1) / 2
            weak = np.asarray(h["/actin/am_weak_count"], dtype=float).reshape(len(ac), -1).sum(1)
            strong = np.asarray(h["/actin/am_strong_count"], dtype=float).reshape(len(ac), -1).sum(1)
            caps = {}
            for species, count in (("actin", int(match.group("na"))), ("myosin", int(match.group("nm")))):
                for kind in ("force", "torque"):
                    key = f"/{species}/cumulative_{kind}_cap_events"
                    total = float(np.asarray(h[key])[-1].sum())
                    caps[f"{species}_{kind}_cap_fraction"] = total / max(1, count * (args.pre_steps + step))
            exposure1, exposure2, promotions, demotions, detach1, detach2 = local_state_exposure(h, args.dt)
            associations = 0
            if "/catch_bond/local_link_events" in h:
                events = np.asarray(h["/catch_bond/local_link_events"], dtype=float)
                associations = int(np.count_nonzero(events[:, 4] == 1)) if len(events) else 0
            times = np.linspace(0.0, step * args.dt, len(ac))
        condition = match.group("condition")
        for index, time in enumerate(times):
            series[condition].append({
                "replicate": int(match.group("rep")), "frame": index, "time": time,
                "actin_nematic": an[index], "myosin_nematic": mn[index],
                "status1": status1[index], "status2": status2[index],
                "weak_heads": weak[index], "strong_heads": strong[index],
            })
        link_time = exposure1 + exposure2
        per_run.append({
            "condition": condition, "replicate": int(match.group("rep")),
            "checkpoint_step": step, "finite": int(finite), "frames": len(ac),
            "final_actin_nematic": an[-1], "final_myosin_nematic": mn[-1],
            "mean_legacy_phi3d": float(phi.mean()), "final_status1": status1[-1],
            "final_status2": status2[-1], "mean_weak_heads": float(weak.mean()),
            "mean_strong_heads": float(strong.mean()), "associations": associations,
            "status2_fraction": exposure2 / link_time if link_time else np.nan,
            "status1_off_rate": detach1 / exposure1 if exposure1 else np.nan,
            "status2_off_rate": detach2 / exposure2 if exposure2 else np.nan,
            "promotions": promotions, "demotions": demotions, **caps,
        })

    write_csv(args.output_dir / "local_aa_assembly_per_run.csv", per_run)
    if not per_run:
        print("no readable pilot trajectories"); return
    summary = []
    numeric = [k for k in per_run[0] if k not in ("condition", "replicate")]
    for condition in sorted({row["condition"] for row in per_run}):
        rows = [row for row in per_run if row["condition"] == condition]
        out = {"condition": condition, "n_runs": len(rows)}
        for field in numeric:
            values = np.asarray([row[field] for row in rows], dtype=float)
            valid = values[np.isfinite(values)]
            out[f"mean_{field}"] = float(valid.mean()) if len(valid) else np.nan
            out[f"sem_{field}"] = float(valid.std(ddof=1) / np.sqrt(len(valid))) if len(valid) > 1 else np.nan
        summary.append(out)
    write_csv(args.output_dir / "local_aa_assembly_summary.csv", summary)

    order = [x for x in ("aa_off", "one_nocatch", "two_nocatch", "two_status2", "two_status1") if x in series]
    fig, axes = plt.subplots(3, 1, figsize=(8.4, 9.0), sharex=True)
    colors = dict(zip(order, plt.cm.viridis(np.linspace(0.05, 0.95, len(order)))))
    for condition in order:
        rows = series[condition]
        reps = sorted({int(r["replicate"]) for r in rows})
        for axis, field in zip(axes, ("myosin_nematic", "actin_nematic", "status2")):
            traces = []
            for rep in reps:
                traces.append(np.asarray([r[field] for r in rows if r["replicate"] == rep]))
            n = min(map(len, traces)); values = np.asarray([x[:n] for x in traces])
            time = np.asarray([r["time"] for r in rows if r["replicate"] == reps[0]])[:n]
            mean = values.mean(0); sem = values.std(0, ddof=1) / np.sqrt(len(values)) if len(values) > 1 else np.zeros(n)
            axis.plot(time, mean, label=condition, color=colors[condition])
            axis.fill_between(time, mean-sem, mean+sem, color=colors[condition], alpha=.16)
    axes[0].set_ylabel("Myosin nematic")
    axes[1].set_ylabel("Actin nematic")
    axes[2].set_ylabel("Status-2 links"); axes[2].set_xlabel("Production time (s)")
    for axis in axes: axis.grid(alpha=.2); axis.spines[["top", "right"]].set_visible(False)
    axes[0].legend(frameon=False, ncol=2)
    fig.tight_layout(); fig.savefig(args.output_dir / "local_aa_assembly_trends.png", dpi=220)
    print(f"analyzed {len(per_run)} trajectories in {args.output_dir}")


if __name__ == "__main__":
    main()
