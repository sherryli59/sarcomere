#!/usr/bin/env python3
"""Audit completeness and matched outcomes of the six-run 300k production array."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np


SEEDS = (831704, 831705, 831706)


def load_json(path: Path):
    with path.open() as stream:
        return json.load(stream)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    runs = []
    integrity_issues = []
    by_arm_seed = {}
    for arm in ("control", "module40"):
        for seed in SEEDS:
            tag = f"A800_M320_{arm}_seed{seed}_750ms"
            trajectory = args.root / "data" / f"{tag}.h5"
            summary_path = args.root / "analysis" / f"{tag}.json"
            live_path = args.root / "analysis" / f"{tag}.live_final.json"
            record = {
                "arm": arm,
                "seed": seed,
                "trajectory": str(trajectory),
                "exists": trajectory.is_file(),
            }
            if not trajectory.is_file():
                integrity_issues.append(f"missing trajectory: {trajectory.name}")
            else:
                try:
                    with h5py.File(trajectory, "r", locking=False) as handle:
                        steps = handle["state/current_step"][:, 0]
                        meta = handle["meta"].attrs
                        record["frames"] = int(len(steps))
                        record["latest_step"] = int(steps[-1]) if len(steps) else None
                        record["dt_s"] = float(meta["dt"])
                        record["n_actins"] = int(handle["actin/center"].shape[1])
                        record["n_myosins"] = int(handle["myosin/center"].shape[1])
                        record["k_mm"] = float(meta.get("k_mm", -1.0))
                        record["k_aa"] = float(meta.get("k_aa", -1.0))
                        record["myomesin_force_max"] = float(
                            meta.get("myomesin_force_max", -1.0)
                        )
                        record["myomesin_radial_k"] = float(
                            meta.get("myomesin_radial_k", -1.0)
                        )
                        record["myomesin_sites_per_filament"] = int(
                            meta.get("myomesin_sites_per_filament", -1)
                        )
                        model = meta.get("myomesin_model", b"")
                        record["myomesin_model"] = (
                            model.decode() if isinstance(model, bytes) else str(model)
                        )
                        record["step_monotonic"] = bool(
                            len(steps) and np.all(np.diff(steps) > 0)
                        )
                except (OSError, KeyError, ValueError) as error:
                    record["hdf5_error"] = str(error)
                if record.get("latest_step") != 300000:
                    integrity_issues.append(
                        f"{trajectory.name}: latest step is {record.get('latest_step')}"
                    )
                if record.get("frames") != 151:
                    integrity_issues.append(
                        f"{trajectory.name}: frame count is {record.get('frames')}, expected 151"
                    )
                if record.get("dt_s") != 2.5e-6:
                    integrity_issues.append(
                        f"{trajectory.name}: dt is {record.get('dt_s')}, expected 2.5e-6"
                    )
                expected_k_mm = 240.0 if arm == "module40" else 0.0
                expected_model = (
                    "sparse_local_dynamic_finite_valency_bare_zone_modules"
                    if arm == "module40" else "disabled"
                )
                expected_values = {
                    "n_actins": 800,
                    "n_myosins": 320,
                    "k_mm": expected_k_mm,
                    "k_aa": 3000.0,
                    "myomesin_force_max": 40.0,
                    "myomesin_radial_k": 300.0,
                    "myomesin_sites_per_filament": 2,
                    "myomesin_model": expected_model,
                }
                for name, expected in expected_values.items():
                    observed = record.get(name)
                    if isinstance(expected, float):
                        matches = observed is not None and np.isclose(
                            observed, expected, rtol=0.0, atol=1e-12
                        )
                    else:
                        matches = observed == expected
                    if not matches:
                        integrity_issues.append(
                            f"{trajectory.name}: {name} is {observed!r}, "
                            f"expected {expected!r}"
                        )
                if not record.get("step_monotonic", False):
                    integrity_issues.append(
                        f"{trajectory.name}: saved steps are not strictly increasing"
                    )

            if summary_path.is_file():
                record["summary"] = load_json(summary_path)
            else:
                integrity_issues.append(f"missing summary: {summary_path.name}")
            if live_path.is_file():
                record["live_final"] = load_json(live_path)
                if not record["live_final"].get("healthy", False):
                    integrity_issues.append(f"unhealthy final state: {trajectory.name}")
            else:
                integrity_issues.append(f"missing final health record: {live_path.name}")
            runs.append(record)
            by_arm_seed[(arm, seed)] = record

    paired = []
    for seed in SEEDS:
        control = by_arm_seed[("control", seed)].get("summary", {})
        module = by_arm_seed[("module40", seed)].get("summary", {})
        control_structure = control.get("structural", {})
        module_structure = module.get("structural", {})
        comparison = {"seed": seed}
        for metric in ("m_band_order", "z_band_order", "a_band_order"):
            first = control_structure.get(metric)
            second = module_structure.get(metric)
            comparison[f"control_{metric}"] = first
            comparison[f"module_{metric}"] = second
            comparison[f"delta_{metric}"] = (
                second - first if first is not None and second is not None else None
            )
        comparison["control_status2_links"] = control_structure.get(
            "status2_links_final"
        )
        comparison["module_status2_links"] = module_structure.get(
            "status2_links_final"
        )
        comparison["control_motor_pair_annealing"] = control.get(
            "motor_pair_annealing"
        )
        comparison["module_motor_pair_annealing"] = module.get(
            "motor_pair_annealing"
        )
        comparison["module_topology"] = module.get("topology")
        comparison["module_gates"] = module.get("gates")
        paired.append(comparison)

    m_deltas = [
        item["delta_m_band_order"] for item in paired
        if item.get("delta_m_band_order") is not None
    ]
    result = {
        "root": str(args.root),
        "integrity_pass": not integrity_issues,
        "integrity_issues": integrity_issues,
        "runs": runs,
        "paired_seed_comparisons": paired,
        "matched_summary": {
            "m_band_delta_median": (
                float(np.median(m_deltas)) if m_deltas else None
            ),
            "m_band_fraction_improved": (
                float(np.mean(np.asarray(m_deltas) > 0)) if m_deltas else None
            ),
        },
    }
    rendered = json.dumps(result, indent=2)
    print(rendered)
    output = args.output or args.root / "analysis" / "production_300k_audit.json"
    output.write_text(rendered + "\n")
    return 0 if result["integrity_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
