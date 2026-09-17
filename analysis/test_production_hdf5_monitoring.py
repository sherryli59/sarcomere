#!/usr/bin/env python3
"""Unit tests for the 300k production HDF5 integrity gates."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_live_production_hdf5 import inspect
from check_production_myomesin_pilot import motor_pair_event_counts


SEEDS = (831704, 831705, 831706)


def write_live_fixture(path: Path, actin_cap_delta: float | list[float]) -> None:
    cap_values = (
        [0.0, float(actin_cap_delta)]
        if np.isscalar(actin_cap_delta) else list(actin_cap_delta)
    )
    n_frames = len(cap_values)
    steps = np.arange(n_frames, dtype=np.int64) * 2000
    with h5py.File(path, "w") as handle:
        meta = handle.create_group("meta")
        meta.attrs["dt"] = 2.5e-6
        meta.attrs["box_x"] = 4.825
        meta.attrs["interaction_ramp_start_step"] = 0
        meta.attrs["interaction_ramp_steps"] = 0
        meta.create_dataset(
            "committed_frames", data=np.asarray([n_frames], dtype=np.int32)
        )
        state = handle.create_group("state")
        state.create_dataset("current_step", data=steps[:, None])
        actin = handle.create_group("actin")
        myosin = handle.create_group("myosin")
        for group, count in ((actin, 2), (myosin, 1)):
            group.create_dataset("center", data=np.zeros((n_frames, count, 3)))
            directions = np.zeros((n_frames, count, 3))
            directions[:, :, 0] = 1.0
            group.create_dataset("direction", data=directions)
        actin_caps = np.zeros((n_frames, 2, 1))
        actin_caps[:, 0, 0] = cap_values
        actin.create_dataset(
            "cumulative_force_cap_events", data=actin_caps,
        )
        myosin.create_dataset(
            "cumulative_force_cap_events", data=np.zeros((n_frames, 1, 1))
        )
        actin.create_dataset(
            "aa_local_status2_count", data=np.zeros((n_frames, 2, 1))
        )


def write_audit_fixture(root: Path) -> None:
    for directory in ("data", "analysis"):
        (root / directory).mkdir(parents=True, exist_ok=True)
    steps = np.arange(0, 300001, 2000, dtype=np.int64)[:, None]
    for arm in ("control", "module40"):
        for seed in SEEDS:
            tag = f"A800_M320_{arm}_seed{seed}_750ms"
            with h5py.File(root / "data" / f"{tag}.h5", "w") as handle:
                state = handle.create_group("state")
                state.create_dataset("current_step", data=steps)
                meta = handle.create_group("meta")
                meta.attrs["dt"] = 2.5e-6
                meta.attrs["k_mm"] = 240.0 if arm == "module40" else 0.0
                meta.attrs["k_aa"] = 3000.0
                meta.attrs["myomesin_force_max"] = 40.0
                meta.attrs["myomesin_radial_k"] = 300.0
                meta.attrs["myomesin_sites_per_filament"] = 2
                meta.attrs["myomesin_model"] = (
                    "sparse_local_dynamic_finite_valency_bare_zone_modules"
                    if arm == "module40" else "disabled"
                )
                actin = handle.create_group("actin")
                actin.create_dataset("center", shape=(151, 800, 3), dtype="f4")
                myosin = handle.create_group("myosin")
                myosin.create_dataset("center", shape=(151, 320, 3), dtype="f4")
            summary = {
                "structural": {
                    "m_band_order": 0.2 if arm == "module40" else 0.1,
                    "z_band_order": 0.0,
                    "a_band_order": 0.0,
                    "status2_links_final": 1,
                },
                "topology": {},
                "gates": {},
            }
            (root / "analysis" / f"{tag}.json").write_text(json.dumps(summary))
            (root / "analysis" / f"{tag}.live_final.json").write_text(
                json.dumps({"healthy": True, "latest_step": 300000})
            )


class LiveCheckerTests(unittest.TestCase):
    @staticmethod
    def add_local_link_offsets(path: Path, status_mode: int) -> None:
        with h5py.File(path, "r+") as handle:
            handle["meta"].attrs["aa_status_mode"] = status_mode
            catch_bond = handle.create_group("catch_bond")
            catch_bond.create_dataset(
                "local_link_frame_offsets",
                data=np.asarray([[0, 0, 1], [2000, 1, 1]], dtype=np.int64),
            )
            catch_bond.create_dataset(
                "local_link_mechanics_v2_frame_offsets",
                data=np.asarray(
                    [[0, 0, 0, 1], [2000, 1, 1, 1]], dtype=np.int64
                ),
            )
            catch_bond.create_dataset(
                "local_link_identity", data=np.zeros((2, 2))
            )
            catch_bond.create_dataset(
                "local_link_mechanics_v2", data=np.zeros((2, 2))
            )

    def test_below_five_percent_cap_fraction_passes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "trajectory.h5"
            # 100 / (2000 steps * 2 actins) = 2.5%.
            write_live_fixture(path, 100.0)
            result = inspect(path)
            self.assertTrue(result["healthy"])
            self.assertAlmostEqual(result["z_band_order"], 1.0)
            self.assertAlmostEqual(
                result["force_cap_fraction_per_step_latest"]["actin"], 0.025
            )

    def test_above_five_percent_cap_fraction_fails(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "trajectory.h5"
            # 220 / (2000 steps * 2 actins) = 5.5%.
            write_live_fixture(path, 220.0)
            result = inspect(path)
            self.assertFalse(result["healthy"])
            self.assertTrue(any("5.50%" in issue for issue in result["issues"]))

    def test_reference_specific_cap_limit_can_accept_known_profile(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "trajectory.h5"
            write_live_fixture(path, 220.0)
            result = inspect(path, actin_force_cap_fraction_max=0.06)
            self.assertTrue(result["healthy"], result["issues"])

    def test_committed_frame_marker_mismatch_fails(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "trajectory.h5"
            write_live_fixture(path, 0.0)
            with h5py.File(path, "r+") as handle:
                handle["meta/committed_frames"][0] = 1
            result = inspect(path)
            self.assertFalse(result["healthy"])
            self.assertTrue(any("committed_frames" in issue for issue in result["issues"]))

    def test_historical_post_ramp_cap_spike_fails_after_recovery(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "trajectory.h5"
            # First interval is 5.5%; latest interval has recovered to zero.
            write_live_fixture(path, [0.0, 220.0, 220.0])
            result = inspect(path)
            self.assertFalse(result["healthy"])
            self.assertEqual(
                result["force_cap_fraction_per_step_latest"]["actin"], 0.0
            )
            self.assertAlmostEqual(
                result["force_cap_fraction_per_step_max_post_ramp"]["actin"],
                0.055,
            )

    def test_active_myomesin_model_requires_frame_offsets(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "trajectory.h5"
            write_live_fixture(path, 0.0)
            with h5py.File(path, "r+") as handle:
                handle["meta"].attrs["k_mm"] = 240.0
                handle["meta"].attrs["myomesin_model"] = (
                    "sparse_local_dynamic_finite_valency_bare_zone_modules"
                )
            result = inspect(path)
            self.assertFalse(result["healthy"])
            self.assertTrue(any("missing frame_offsets" in issue for issue in result["issues"]))

    def test_productive_register_does_not_require_motor_support_table(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "trajectory.h5"
            write_live_fixture(path, 0.0)
            self.add_local_link_offsets(path, status_mode=7)
            result = inspect(path)
            self.assertTrue(result["healthy"], result["issues"])

    def test_motor_transmitted_register_requires_motor_support_table(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "trajectory.h5"
            write_live_fixture(path, 0.0)
            self.add_local_link_offsets(path, status_mode=8)
            result = inspect(path)
            self.assertFalse(result["healthy"])
            self.assertTrue(
                any("local_link_motor_support" in issue for issue in result["issues"])
            )


class AnalyzerEventTests(unittest.TestCase):
    def test_identity_release_excludes_empty_detachment_notification(self) -> None:
        events = np.full((3, 20), -1.0)
        events[0, 3] = 1
        events[0, 8:10] = (4, 5)
        events[1, 3] = 5
        events[1, 4:6] = (4, 5)
        events[2, 3] = 5
        counts = motor_pair_event_counts(events)
        self.assertEqual(counts["selected"], 1)
        self.assertEqual(counts["released_by_actn2_detachment"], 1)
        self.assertEqual(counts["actn2_detachment_without_held_identity"], 1)


class FinalAuditTests(unittest.TestCase):
    def run_audit(self, root: Path) -> subprocess.CompletedProcess[str]:
        script = Path(__file__).with_name("audit_production_300k.py")
        return subprocess.run(
            [sys.executable, str(script), str(root)],
            check=False,
            capture_output=True,
            text=True,
        )

    def test_complete_matched_six_run_fixture_passes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_audit_fixture(root)
            completed = self.run_audit(root)
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            result = json.loads(
                (root / "analysis" / "production_300k_audit.json").read_text()
            )
            self.assertTrue(result["integrity_pass"])
            self.assertEqual(len(result["runs"]), 6)

    def test_wrong_frame_count_fails(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_audit_fixture(root)
            path = root / "data" / "A800_M320_control_seed831704_750ms.h5"
            with h5py.File(path, "r+") as handle:
                del handle["state/current_step"]
                handle["state"].create_dataset(
                    "current_step",
                    data=np.arange(0, 300001, 4000, dtype=np.int64)[:, None],
                )
            completed = self.run_audit(root)
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("expected 151", completed.stdout)

    def test_wrong_arm_parameter_fails(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_audit_fixture(root)
            path = root / "data" / "A800_M320_module40_seed831704_750ms.h5"
            with h5py.File(path, "r+") as handle:
                handle["meta"].attrs["k_mm"] = 0.0
            completed = self.run_audit(root)
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("k_mm is 0.0, expected 240.0", completed.stdout)


if __name__ == "__main__":
    unittest.main()
