#!/usr/bin/env python3
"""Audit the matched 5-s Cell A control and Cell A + sparse-myomesin pair."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np


ARMS = {"control": 0.0, "module40": 240.0}
SEED = 20260805


def load_json(path: Path) -> dict:
    with path.open() as stream:
        return json.load(stream)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--dt", type=float, default=1e-5)
    parser.add_argument("--nsteps", type=int, default=500_000)
    parser.add_argument("--save-every", type=int, default=500)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    expected_frames = args.nsteps // args.save_every + 1
    duration_ms = args.nsteps * args.dt * 1e3
    duration_tag = f"{duration_ms:.6g}ms"
    integrity_issues: list[str] = []
    records: dict[str, dict] = {}

    for arm, expected_k_mm in ARMS.items():
        tag = f"A800_M320_cellA_{arm}_seed{SEED}_{duration_tag}"
        trajectory = args.root / "data" / f"{tag}.h5"
        summary_path = args.root / "analysis" / f"{tag}.json"
        live_path = args.root / "analysis" / f"{tag}.live_final.json"
        record = {
            "arm": arm,
            "trajectory": str(trajectory),
            "summary_path": str(summary_path),
            "live_final_path": str(live_path),
        }
        if not trajectory.is_file():
            integrity_issues.append(f"missing trajectory: {trajectory.name}")
        else:
            try:
                with h5py.File(trajectory, "r", locking=False) as handle:
                    steps = handle["state/current_step"][:, 0].astype(int)
                    meta = handle["meta"].attrs
                    model = meta.get("myomesin_model", b"")
                    if isinstance(model, bytes):
                        model = model.decode()
                    record.update({
                        "frames": int(len(steps)),
                        "latest_step": int(steps[-1]) if len(steps) else None,
                        "strictly_increasing_steps": bool(
                            len(steps) and np.all(np.diff(steps) > 0)
                        ),
                        "dt_s": float(meta.get("dt", -1.0)),
                        "n_actins": int(handle["actin/center"].shape[1]),
                        "n_myosins": int(handle["myosin/center"].shape[1]),
                        "k_mm": float(meta.get("k_mm", -1.0)),
                        "k_aa": float(meta.get("k_aa", -1.0)),
                        "myosin_steric_k": float(
                            meta.get("myosin_steric_k", -1.0)
                        ),
                        "myosin_exclusion_distance": float(
                            meta.get("myosin_exclusion_distance", -1.0)
                        ),
                        "aa_status_mode": int(meta.get("aa_status_mode", -1)),
                        "aa_status2_force_on": float(
                            meta.get("aa_status2_force_on", -1.0)
                        ),
                        "aa_status2_force_off": float(
                            meta.get("aa_status2_force_off", -1.0)
                        ),
                        "aa_energy_gated_association": int(
                            meta.get("aa_energy_gated_association", -1)
                        ),
                        "myomesin_model": str(model),
                        "myomesin_force_max": float(
                            meta.get("myomesin_force_max", -1.0)
                        ),
                    })
            except (OSError, KeyError, ValueError) as error:
                record["hdf5_error"] = str(error)

            expected_model = (
                "sparse_local_dynamic_finite_valency_bare_zone_modules"
                if arm == "module40" else "disabled"
            )
            expected = {
                "latest_step": args.nsteps,
                "frames": expected_frames,
                "strictly_increasing_steps": True,
                "dt_s": args.dt,
                "n_actins": 800,
                "n_myosins": 320,
                "k_mm": expected_k_mm,
                "k_aa": 1000.0,
                "myosin_steric_k": 1000.0,
                "myosin_exclusion_distance": 0.0,
                "aa_status_mode": 7,
                "aa_status2_force_on": 0.25,
                "aa_status2_force_off": 0.10,
                "aa_energy_gated_association": 0,
                "myomesin_model": expected_model,
                "myomesin_force_max": 40.0,
            }
            for name, wanted in expected.items():
                observed = record.get(name)
                if isinstance(wanted, float):
                    matches = observed is not None and np.isclose(
                        observed, wanted, rtol=0.0, atol=1e-12
                    )
                else:
                    matches = observed == wanted
                if not matches:
                    integrity_issues.append(
                        f"{trajectory.name}: {name}={observed!r}, expected {wanted!r}"
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
            integrity_issues.append(f"missing health record: {live_path.name}")
        records[arm] = record

    control = records["control"].get("summary", {})
    module = records["module40"].get("summary", {})
    control_structure = control.get("structural", {})
    module_structure = module.get("structural", {})
    paired = {}
    for metric in ("m_band_order", "z_band_order", "a_band_order"):
        first = control_structure.get(metric)
        second = module_structure.get(metric)
        paired[f"control_{metric}"] = first
        paired[f"module_{metric}"] = second
        paired[f"delta_{metric}"] = (
            second - first if first is not None and second is not None else None
        )
    paired["control_status2_links"] = control_structure.get("status2_links_final")
    paired["module_status2_links"] = module_structure.get("status2_links_final")
    paired["control_z_band_order_tail_median"] = control_structure.get(
        "z_band_order_tail_median"
    )
    paired["module_z_band_order_tail_median"] = module_structure.get(
        "z_band_order_tail_median"
    )

    topology = module.get("topology", {})
    module_gates = module.get("gates", {})
    scientific_gates = {
        # Historical Cell A ends at 0.708 and has a last-20% median of 0.699.
        # Requiring 0.60 for both makes "replicates Cell A" a strong,
        # persistent-order statement without requiring the same noisy
        # multithreaded trajectory.
        "control_reproduces_cell_a_z_band": bool(
            paired.get("control_z_band_order") is not None
            and paired["control_z_band_order"] >= 0.60
            and paired.get("control_z_band_order_tail_median") is not None
            and paired["control_z_band_order_tail_median"] >= 0.60
        ),
        "module_reduces_bridge_mismatch": bool(
            topology.get("annealing_modules_evaluated", 0) > 0
            and topology.get("annealing_reduction_fraction_median") is not None
            and topology["annealing_reduction_fraction_median"] > 0
            and topology.get("annealing_fraction_positive", 0) >= 0.5
        ),
        "m_band_order_improves": bool(
            paired.get("delta_m_band_order") is not None
            and paired["delta_m_band_order"] > 0
        ),
        "z_band_not_materially_damaged": bool(
            paired.get("delta_z_band_order") is not None
            and paired["delta_z_band_order"] >= -0.10
        ),
        "finite_local_topology": bool(
            module_gates.get("finite_valency", False)
            and module_gates.get("unique_pair_edges", False)
            and module_gates.get("no_system_spanning_component", False)
        ),
    }
    result = {
        "root": str(args.root),
        "expected_duration_ms": duration_ms,
        "integrity_pass": not integrity_issues,
        "integrity_issues": integrity_issues,
        "runs": records,
        "paired_outcome": paired,
        "module_topology": topology,
        "module_gates": module_gates,
        "scientific_gates": scientific_gates,
        "scientific_pass": all(scientific_gates.values()),
    }
    rendered = json.dumps(result, indent=2)
    print(rendered)
    output = args.output or args.root / "analysis" / "cellA_sparse_myomesin_audit.json"
    output.write_text(rendered + "\n")
    return 0 if result["integrity_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
