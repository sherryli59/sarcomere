import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np

from analysis.check_myomesin_annealing import evaluate


class MyomesinAnnealingCheckTest(unittest.TestCase):
    def write_trajectory(
        self, path: Path, formation_nm: float, saved_mismatches_nm: list[float]
    ) -> None:
        with h5py.File(path, "w") as handle:
            meta = handle.create_group("meta")
            meta.attrs["dt"] = 1e-3
            meta.attrs["k_mm"] = 120.0
            meta.attrs["myomesin_radial_k"] = 300.0
            meta.attrs["myomesin_force_max"] = 20.0
            group = handle.create_group("myomesin")
            offsets = [[0, 0, 0]]
            rows = []
            for frame, mismatch_nm in enumerate(saved_mismatches_nm, start=1):
                offsets.append([frame, len(rows), 1])
                rows.append([
                    0, 1, 0, 2, 0, mismatch_nm * 1e-3, 0.11, 5.0
                ])
            group.create_dataset("frame_offsets", data=np.asarray(offsets))
            group.create_dataset("bridge_state", data=np.asarray(rows))
            group.create_dataset("events", data=np.asarray([[
                0, 1, 1, 0, 2, formation_nm * 1e-3, 0.11, 0.0, 0.0
            ]]))

    def test_uses_exact_association_event_as_baseline(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "anneals.h5"
            self.write_trajectory(path, 150.0, [100.0, 90.0, 80.0])
            result = evaluate(path, [(0, 2)], 2.0, 0.10)
            pair = result["pairs"][0]
            self.assertEqual(pair["formation_source"], "association_event")
            self.assertAlmostEqual(pair["formation_mismatch_nm"], 150.0)
            self.assertAlmostEqual(pair["final_window_median_nm"], 90.0)
            self.assertTrue(result["pass"])

    def test_rejects_final_mismatch_above_formation_mismatch(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "worsens.h5"
            self.write_trajectory(path, 100.0, [105.0, 110.0, 115.0])
            result = evaluate(path, [(0, 2)], 2.0, 0.0)
            self.assertLess(result["pairs"][0]["reduction_nm"], 0.0)
            self.assertFalse(result["pass"])


if __name__ == "__main__":
    unittest.main()
