#!/usr/bin/env python3
"""Tests for the exact historical Cell A trajectory comparator."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
import sys

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_cell_a_backward_compat import (
    CORE_PATHS,
    OPTIONAL_TREATMENT_PATHS,
    compare,
)


def shape_for(path: str) -> tuple[int, ...]:
    if path == "state/current_step":
        return (3, 1)
    if path.startswith("actin/"):
        return (3, 2, 1) if any(
            token in path
            for token in (
                "status", "count", "am_q",
            )
        ) else (3, 2, 3)
    if path.startswith("myosin/"):
        return (3, 1, 3)
    if "am_head" in path:
        return (3, 1, 2)
    return (3, 2)


def write_trajectory(path: Path) -> None:
    with h5py.File(path, "w") as handle:
        for name in CORE_PATHS:
            shape = shape_for(name)
            dtype = np.int32 if any(
                token in name
                for token in ("current_step", "status", "count", "state_current", "bound_actin")
            ) else np.float64
            handle.create_dataset(name, data=np.zeros(shape, dtype=dtype))


def add_treatment_state(path: Path) -> None:
    with h5py.File(path, "r+") as handle:
        handle.create_dataset(
            "myomesin/bridge_state", data=np.zeros((2, 8), dtype=np.float64)
        )
        handle.create_dataset(
            "myomesin/events", data=np.zeros((1, 9), dtype=np.float64)
        )
        handle.create_dataset(
            "myomesin/frame_offsets", data=np.zeros((3, 3), dtype=np.int32)
        )
        handle.create_dataset(
            "myomesin/next_id_by_frame", data=np.zeros((3, 2), dtype=np.int32)
        )


class CompatibilityComparatorTest(unittest.TestCase):
    def test_identical_trajectories_pass(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            old = Path(directory) / "old.h5"
            new = Path(directory) / "new.h5"
            write_trajectory(old)
            write_trajectory(new)
            self.assertTrue(all(record["equal"] for record in compare(old, new)))

    def test_reports_first_different_frame(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            old = Path(directory) / "old.h5"
            new = Path(directory) / "new.h5"
            write_trajectory(old)
            write_trajectory(new)
            with h5py.File(new, "r+") as handle:
                handle["actin/center"][1, 0, 2] = 0.125
            records = {record["path"]: record for record in compare(old, new)}
            mismatch = records["actin/center"]
            self.assertFalse(mismatch["equal"])
            self.assertEqual(mismatch["first_frame"], 1)
            self.assertEqual(mismatch["unequal"], 1)
            self.assertAlmostEqual(mismatch["max_abs_difference"], 0.125)

    def test_active_treatment_state_is_compared_exactly(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            old = Path(directory) / "old.h5"
            new = Path(directory) / "new.h5"
            write_trajectory(old)
            write_trajectory(new)
            add_treatment_state(old)
            add_treatment_state(new)
            with h5py.File(new, "r+") as handle:
                handle["myomesin/bridge_state"][1, 5] = 0.01
            records = {record["path"]: record for record in compare(old, new)}
            self.assertEqual(
                set(OPTIONAL_TREATMENT_PATHS),
                set(records).intersection(OPTIONAL_TREATMENT_PATHS),
            )
            self.assertFalse(records["myomesin/bridge_state"]["equal"])
            self.assertEqual(records["myomesin/bridge_state"]["first_frame"], 1)

    def test_treatment_dataset_missing_from_one_side_fails(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            old = Path(directory) / "old.h5"
            new = Path(directory) / "new.h5"
            write_trajectory(old)
            write_trajectory(new)
            add_treatment_state(old)
            records = {record["path"]: record for record in compare(old, new)}
            self.assertFalse(records["myomesin/events"]["equal"])
            self.assertIn("missing", records["myomesin/events"]["reason"])


if __name__ == "__main__":
    unittest.main()
