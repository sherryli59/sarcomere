#!/usr/bin/env python3
"""Audit every evidence item required to close the Cell A replication goal."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from compare_cell_a_backward_compat import compare


EXPECTED_TREATMENT = {
    "frames": 1001,
    "latest_step": 500_000,
    "dt_s": 1e-5,
    "n_actins": 800,
    "n_myosins": 320,
    "k_mm": 240.0,
    "k_aa": 1000.0,
    "myosin_steric_k": 1000.0,
    "myosin_exclusion_distance": 0.0,
    "aa_status_mode": 7,
    "myomesin_model": "sparse_local_dynamic_finite_valency_bare_zone_modules",
    "myomesin_force_max": 40.0,
}


def load_json(path: Path) -> dict:
    with path.open() as stream:
        return json.load(stream)


def parse_params(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if not line or line.lstrip().startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.strip()] = value.strip()
    return values


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def matches(observed: object, expected: object) -> bool:
    if isinstance(expected, float):
        return observed is not None and bool(np.isclose(
            observed, expected, rtol=0.0, atol=1e-12
        ))
    return observed == expected


def evaluate(
    replication: dict,
    paired: dict,
    corrected_health: dict,
    params: dict[str, str],
    observed_binary_sha256: str,
    expected_binary: str,
    expected_binary_sha256: str,
    treatment_comparison: list[dict],
) -> dict:
    module_record = paired.get("runs", {}).get("module40", {})
    module_health = module_record.get("live_final", {})
    module_gates = paired.get("module_gates", {})
    treatment_config = {
        key: {
            "observed": module_record.get(key),
            "expected": expected,
            "passed": matches(module_record.get(key), expected),
        }
        for key, expected in EXPECTED_TREATMENT.items()
    }
    exact_failures = [
        record for record in treatment_comparison if not record.get("equal", False)
    ]
    exact_treatment = (
        len(treatment_comparison) == 28 and not exact_failures
    )
    provenance = {
        "params_executable": params.get("executable"),
        "expected_executable": expected_binary,
        "params_sha256": params.get("executable_sha256"),
        "expected_sha256": expected_binary_sha256,
        "observed_sha256": observed_binary_sha256,
    }
    provenance["passed"] = bool(
        provenance["params_executable"] == expected_binary
        and provenance["params_sha256"] == expected_binary_sha256
        and provenance["observed_sha256"] == expected_binary_sha256
    )
    gates = {
        "corrected_control_replicates_cell_a": bool(
            replication.get("replicates_cell_a", False)
            and all(replication.get("gates", {}).values())
        ),
        "corrected_control_hdf5_healthy": bool(
            corrected_health.get("healthy", False)
            and corrected_health.get("latest_step") == 500_000
            and corrected_health.get("visible_frames") == 1001
            and corrected_health.get("finite_coordinates", False)
        ),
        "corrected_binary_provenance": bool(provenance["passed"]),
        "paired_run_integrity": bool(paired.get("integrity_pass", False)),
        "treatment_configuration_retained": all(
            record["passed"] for record in treatment_config.values()
        ),
        "treatment_hdf5_healthy": bool(module_health.get("healthy", False)),
        "treatment_sparse_topology_valid": all(
            bool(module_gates.get(name, False)) for name in (
                "finite_valency",
                "unique_pair_edges",
                "no_system_spanning_component",
                "population_annealing",
            )
        ),
        "treatment_path_exact_across_patch": exact_treatment,
    }
    return {
        "goal": (
            "new k_mm=0 control reproduces historical Cell A while retaining "
            "the sparse-myomesin treatment"
        ),
        "gates": gates,
        "goal_complete": all(gates.values()),
        "corrected_control_replication": replication,
        "corrected_control_health": corrected_health,
        "corrected_binary_provenance": provenance,
        "treatment_configuration": treatment_config,
        "treatment_module_gates": module_gates,
        "treatment_scientific_observations": paired.get("scientific_gates", {}),
        "treatment_exact_comparison": {
            "datasets_compared": len(treatment_comparison),
            "failures": exact_failures,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replication", type=Path, required=True)
    parser.add_argument("--paired-audit", type=Path, required=True)
    parser.add_argument("--corrected-health", type=Path, required=True)
    parser.add_argument("--corrected-params", type=Path, required=True)
    parser.add_argument("--expected-binary", type=Path, required=True)
    parser.add_argument("--expected-binary-sha256", required=True)
    parser.add_argument("--treatment-reference", type=Path, required=True)
    parser.add_argument("--treatment-candidate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    treatment_comparison = compare(
        args.treatment_reference, args.treatment_candidate
    )
    report = evaluate(
        load_json(args.replication),
        load_json(args.paired_audit),
        load_json(args.corrected_health),
        parse_params(args.corrected_params),
        sha256(args.expected_binary),
        str(args.expected_binary),
        args.expected_binary_sha256,
        treatment_comparison,
    )
    report["evidence_paths"] = {
        "replication": str(args.replication),
        "paired_audit": str(args.paired_audit),
        "corrected_health": str(args.corrected_health),
        "corrected_params": str(args.corrected_params),
        "treatment_reference": str(args.treatment_reference),
        "treatment_candidate": str(args.treatment_candidate),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    rendered = json.dumps(report, indent=2)
    args.output.write_text(rendered + "\n")
    print(rendered)
    return 0 if report["goal_complete"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
