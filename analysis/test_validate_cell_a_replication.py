"""Unit tests for the historical Cell A phenotype verifier."""

from __future__ import annotations

from pathlib import Path
import sys

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from validate_cell_a_replication import evaluate, summarize, z_band_series


def fixture(final: float, tail: float, final_step: int = 500_000) -> dict:
    z_order = np.concatenate((np.linspace(0.05, tail, 800), np.full(200, tail)))
    z_order = np.append(z_order, final)
    return {
        "steps": np.arange(1001) * (final_step // 1000),
        "time_ms": np.arange(1001) * 5.0,
        "z_order": z_order,
        "metadata": {
            "dt": 1e-5,
            "k_mm": 0.0,
            "k_aa": 1000.0,
            "myosin_steric_k": 1000.0,
            "myosin_exclusion_distance": 0.0,
            "aa_status_mode": 7,
            "n_actins": 800,
            "n_myosins": 320,
        },
    }


def test_strong_persistent_candidate_passes() -> None:
    report = evaluate(fixture(0.708, 0.699), fixture(0.65, 0.63))
    assert report["replicates_cell_a"] is True
    assert all(report["gates"].values())


def test_transient_endpoint_fails_tail_gate() -> None:
    report = evaluate(fixture(0.708, 0.699), fixture(0.65, 0.40))
    assert report["gates"]["strong_endpoint_z_order"] is True
    assert report["gates"]["sustained_tail_z_order"] is False
    assert report["replicates_cell_a"] is False


def test_wrong_active_myomesin_control_fails_metadata() -> None:
    candidate = fixture(0.65, 0.63)
    candidate["metadata"]["k_mm"] = 240.0
    report = evaluate(fixture(0.708, 0.699), candidate)
    assert report["gates"]["cell_a_metadata"] is False
    assert report["replicates_cell_a"] is False


def test_summary_uses_final_twenty_percent() -> None:
    series = fixture(0.65, 0.63)
    summary = summarize(series)
    assert np.isclose(summary["z_band_order_tail_median"], 0.63)
    assert summary["frames"] == 1001


def test_uneven_hdf_frames_emit_failed_completion_gate(tmp_path: Path) -> None:
    path = tmp_path / "partial.h5"
    with h5py.File(path, "w") as handle:
        handle.create_dataset(
            "state/current_step", data=np.array([[0], [500], [1000]])
        )
        handle.create_dataset("actin/center", data=np.zeros((2, 800, 3)))
        directions = np.zeros((3, 800, 3))
        directions[:, :, 0] = 1.0
        handle.create_dataset("actin/direction", data=directions)
        handle.create_dataset("myosin/center", data=np.zeros((3, 320, 3)))
        meta = handle.create_group("meta")
        meta.attrs.update({
            "dt": 1e-5,
            "box_x": 4.825,
            "actin_length": 1.0,
            "k_mm": 0.0,
            "k_aa": 1000.0,
            "myosin_steric_k": 1000.0,
            "myosin_exclusion_distance": 0.0,
            "aa_status_mode": 7,
        })
    candidate = z_band_series(path)
    assert candidate["frame_counts"]["actin/center"] == 2
    assert len(candidate["steps"]) == 2
    report = evaluate(fixture(0.708, 0.699), candidate)
    assert report["gates"]["complete_5s_trajectory"] is False
    assert report["replicates_cell_a"] is False
