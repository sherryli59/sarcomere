"""Tests for the final Cell A goal-completion auditor."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_cell_a_goal_completion import EXPECTED_TREATMENT, evaluate


SHA = "a" * 64
EXE = "/scratch/frozen/sarcomere"


def valid_inputs() -> tuple:
    replication = {
        "replicates_cell_a": True,
        "gates": {
            "complete_5s_trajectory": True,
            "cell_a_metadata": True,
            "strong_endpoint_z_order": True,
            "sustained_tail_z_order": True,
        },
    }
    module = dict(EXPECTED_TREATMENT)
    module["live_final"] = {"healthy": True}
    paired = {
        "integrity_pass": True,
        "runs": {"module40": module},
        "module_gates": {
            "finite_valency": True,
            "unique_pair_edges": True,
            "no_system_spanning_component": True,
            "population_annealing": True,
        },
        "scientific_gates": {"m_band_order_improves": False},
    }
    health = {
        "healthy": True,
        "latest_step": 500_000,
        "visible_frames": 1001,
        "finite_coordinates": True,
    }
    params = {"executable": EXE, "executable_sha256": SHA}
    comparison = [
        {"path": f"dataset/{index}", "equal": True} for index in range(28)
    ]
    return replication, paired, health, params, comparison


def run(inputs: tuple) -> dict:
    replication, paired, health, params, comparison = inputs
    return evaluate(
        replication, paired, health, params, SHA, EXE, SHA, comparison
    )


def test_all_required_evidence_passes() -> None:
    report = run(valid_inputs())
    assert report["goal_complete"]
    assert all(report["gates"].values())


def test_global_treatment_improvement_is_observed_but_not_retention_gate() -> None:
    report = run(valid_inputs())
    assert not report["treatment_scientific_observations"][
        "m_band_order_improves"
    ]
    assert report["goal_complete"]


def test_incomplete_corrected_control_fails() -> None:
    inputs = list(valid_inputs())
    inputs[2] = deepcopy(inputs[2])
    inputs[2]["latest_step"] = 499_500
    report = run(tuple(inputs))
    assert not report["gates"]["corrected_control_hdf5_healthy"]
    assert not report["goal_complete"]


def test_treatment_topology_or_exact_state_failure_blocks_completion() -> None:
    inputs = list(valid_inputs())
    inputs[1] = deepcopy(inputs[1])
    inputs[1]["module_gates"]["finite_valency"] = False
    inputs[4] = deepcopy(inputs[4])
    inputs[4][5]["equal"] = False
    report = run(tuple(inputs))
    assert not report["gates"]["treatment_sparse_topology_valid"]
    assert not report["gates"]["treatment_path_exact_across_patch"]
    assert not report["goal_complete"]


def test_wrong_binary_or_treatment_configuration_blocks_completion() -> None:
    inputs = list(valid_inputs())
    inputs[1] = deepcopy(inputs[1])
    inputs[1]["runs"]["module40"]["k_mm"] = 0.0
    inputs[3] = deepcopy(inputs[3])
    inputs[3]["executable_sha256"] = "b" * 64
    report = run(tuple(inputs))
    assert not report["gates"]["treatment_configuration_retained"]
    assert not report["gates"]["corrected_binary_provenance"]
    assert not report["goal_complete"]
