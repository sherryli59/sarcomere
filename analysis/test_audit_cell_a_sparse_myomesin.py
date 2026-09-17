"""Tests for the matched Cell A + sparse-myomesin final auditor."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

import h5py
import numpy as np


HERE = Path(__file__).resolve().parent
AUDITOR = HERE / "audit_cell_a_sparse_myomesin.py"


def write_run(root: Path, arm: str, k_mm: float, m_order: float) -> None:
    tag = f"A800_M320_cellA_{arm}_seed20260805_0.005ms"
    for directory in ("data", "analysis"):
        (root / directory).mkdir(parents=True, exist_ok=True)
    trajectory = root / "data" / f"{tag}.h5"
    with h5py.File(trajectory, "w") as handle:
        meta = handle.create_group("meta")
        meta.attrs.update({
            "dt": 2.5e-6,
            "k_mm": k_mm,
            "k_aa": 1000.0,
            "myosin_steric_k": 1000.0,
            "myosin_exclusion_distance": 0.0,
            "aa_status_mode": 7,
            "aa_status2_force_on": 0.25,
            "aa_status2_force_off": 0.10,
            "aa_energy_gated_association": 0,
            "myomesin_model": (
                "sparse_local_dynamic_finite_valency_bare_zone_modules"
                if arm == "module40" else "disabled"
            ),
            "myomesin_force_max": 40.0,
        })
        state = handle.create_group("state")
        state.create_dataset("current_step", data=np.array([[0], [2]]))
        actin = handle.create_group("actin")
        actin.create_dataset("center", shape=(2, 800, 3), dtype=float)
        myosin = handle.create_group("myosin")
        myosin.create_dataset("center", shape=(2, 320, 3), dtype=float)
    summary = {
        "structural": {
            "m_band_order": m_order,
            "z_band_order": 0.7,
            "z_band_order_tail_median": 0.68,
            "a_band_order": 0.4,
            "status2_links_final": 10,
        }
    }
    if arm == "module40":
        summary.update({
            "topology": {
                "annealing_modules_evaluated": 3,
                "annealing_reduction_fraction_median": 0.2,
                "annealing_fraction_positive": 2 / 3,
            },
            "gates": {
                "finite_valency": True,
                "unique_pair_edges": True,
                "no_system_spanning_component": True,
            },
        })
    (root / "analysis" / f"{tag}.json").write_text(json.dumps(summary))
    (root / "analysis" / f"{tag}.live_final.json").write_text(
        json.dumps({"healthy": True})
    )


def test_complete_pair_separates_integrity_and_scientific_result(tmp_path):
    write_run(tmp_path, "control", 0.0, 0.2)
    write_run(tmp_path, "module40", 240.0, 0.4)
    output = tmp_path / "audit.json"
    result = subprocess.run(
        [
            sys.executable,
            str(AUDITOR),
            str(tmp_path),
            "--dt", "2.5e-6",
            "--nsteps", "2",
            "--save-every", "2",
            "--output", str(output),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    report = json.loads(output.read_text())
    assert report["integrity_pass"] is True
    assert report["scientific_pass"] is True
    assert report["scientific_gates"]["control_reproduces_cell_a_z_band"] is True
    assert np.isclose(report["paired_outcome"]["delta_m_band_order"], 0.2)


def test_missing_arm_is_an_integrity_failure(tmp_path):
    write_run(tmp_path, "control", 0.0, 0.2)
    result = subprocess.run(
        [
            sys.executable,
            str(AUDITOR),
            str(tmp_path),
            "--dt", "2.5e-6",
            "--nsteps", "2",
            "--save-every", "2",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 1
    report = json.loads(result.stdout)
    assert report["integrity_pass"] is False
    assert any("module40" in issue for issue in report["integrity_issues"])
