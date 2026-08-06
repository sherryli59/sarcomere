import unittest

import numpy as np

from overlap_model import (
    Bond,
    Parameters,
    State,
    aa_overlap,
    am_overlaps,
    default_state,
    force_components,
    mechanism_parameters,
    simulate,
    total_force,
)


class OverlapModelTests(unittest.TestCase):
    def test_default_geometry(self):
        state = default_state(aa_overlap=0.10, zm_distance=1.125)
        self.assertAlmostEqual(aa_overlap(state.q), 0.10)
        left, right = am_overlaps(state.q, Parameters())
        self.assertAlmostEqual(left, right)
        self.assertGreater(left, 0.0)
        self.assertLess(left, Parameters().myosin_half_length)

    def test_internal_forces_sum_to_zero(self):
        p = mechanism_parameters("register_plus_titin")
        state = default_state(aa_overlap=0.06, zm_distance=1.08)
        state.bonds = [Bond(0.10), Bond(0.08)]
        for name, force in force_components(state, p).items():
            if name == "mline_anchor":
                continue
            self.assertAlmostEqual(float(force.sum()), 0.0, places=12, msg=name)

    def test_motor_increases_am_overlap_direction(self):
        p = Parameters(thermal_scale=0.0, duration=0.01, actn2_kon=0.0, initial_bonds=0)
        initial = default_state()
        before = np.mean(am_overlaps(initial.q, p))
        result = simulate(p, seed=1, initial=initial)
        after = np.mean(result["overlaps"][-1, 1:])
        self.assertGreater(after, before)

    def test_turnover_rebinding_relaxes_reference(self):
        p = Parameters(
            actn2_koff=1.0e6,
            actn2_kon=1.0e6,
            actn2_catch_force=0.0,
            dt=1.0e-4,
            duration=1.0e-4,
            save_every=1,
        )
        state = default_state(aa_overlap=0.04)
        state.bonds = [Bond(0.10) for _ in range(p.n_actn2_sites)]
        result = simulate(p, seed=3, initial=state, motor_on=False)
        # The exact number is stochastic, but the large rates make replacement
        # overwhelmingly likely and therefore erase most stored elastic force.
        final_force = result["force_actn2_elastic"][-1, 0]
        self.assertLess(abs(final_force), 5.0)

    def test_adhesion_force_increases_aa_overlap(self):
        p = Parameters(adhesion_gamma=4.0)
        state = default_state(aa_overlap=0.05)
        state.bonds = []
        f = force_components(state, p, motor_on=False)["actn2_adhesion"]
        self.assertGreater(f[0], 0.0)
        self.assertLess(f[1], 0.0)

    def test_titin_resists_over_shortening(self):
        p = Parameters(titin_k=100.0, titin_rest=1.125)
        state = default_state(zm_distance=1.00)
        state.bonds = []
        f = force_components(state, p, motor_on=False)["titin"]
        self.assertGreater(f[0], 0.0)
        self.assertLess(f[2], 0.0)

    def test_deterministic_seed(self):
        p = Parameters(duration=0.02)
        first = simulate(p, seed=91)
        second = simulate(p, seed=91)
        np.testing.assert_array_equal(first["coords"], second["coords"])
        np.testing.assert_array_equal(first["bond_count"], second["bond_count"])


if __name__ == "__main__":
    unittest.main()
