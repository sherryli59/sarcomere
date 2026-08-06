import unittest
from dataclasses import replace

import numpy as np

from periodic_model import (
    aa_overlaps,
    am_overlaps,
    default_periodic_state,
    force_components,
    periodic_parameters,
    simulate_periodic,
    wrap_delta,
)
from variable_cell_model import (
    change_cell_length,
    default_variable_state,
    internal_cell_force_components,
    simulate_variable_cell,
    variable_cell_parameters,
)


class PeriodicModelTests(unittest.TestCase):
    def test_default_exact_geometry(self):
        p = periodic_parameters("baseline_turnover")
        state = default_periodic_state(0.10, p.box_length)
        np.testing.assert_allclose(aa_overlaps(state, p), [0.10, 0.10])
        self.assertLess(np.ptp(am_overlaps(state, p)), 1e-12)

    def test_every_internal_component_is_reciprocal(self):
        for name in (
            "baseline_turnover",
            "permanent_actn2",
            "johann_adhesion",
            "johann_entropic",
            "titin",
            "entropic_plus_titin",
            "register_plus_titin",
        ):
            p = periodic_parameters(name)
            state = default_periodic_state(0.07, p.box_length)
            for component, forces in force_components(state, p).items():
                self.assertAlmostEqual(
                    float(np.sum(forces)), 0.0, places=11, msg=(name, component)
                )

    def test_minimum_image_is_periodic(self):
        length = 4.3
        self.assertAlmostEqual(wrap_delta(4.2, length), -0.1)
        self.assertAlmostEqual(wrap_delta(-4.2, length), 0.1)

    def test_fixed_cell_seed_is_deterministic(self):
        p = replace(
            periodic_parameters("baseline_turnover"),
            duration=0.01,
            save_every=1,
        )
        first = simulate_periodic(p, 17)
        second = simulate_periodic(p, 17)
        np.testing.assert_array_equal(first["coords"], second["coords"])


class VariableCellTests(unittest.TestCase):
    def test_affine_cell_change_preserves_fractional_coordinates(self):
        state = default_variable_state(0.10, 4.30)
        before = state.particles.q / state.length
        change_cell_length(state, 4.70)
        np.testing.assert_allclose(state.particles.q / state.length, before)

    def test_motor_cell_force_is_contractile(self):
        cp = variable_cell_parameters("baseline_turnover")
        state = default_variable_state(0.10, 4.30)
        forces = internal_cell_force_components(state, cp)
        self.assertLess(forces["motor"], 0.0)

    def test_compressed_titin_cell_force_is_expansile(self):
        cp = variable_cell_parameters("titin")
        state = default_variable_state(0.10, 3.90)
        forces = internal_cell_force_components(state, cp, motor_on=False)
        self.assertGreater(forces["titin"], 0.0)

    def test_variable_cell_seed_is_deterministic(self):
        cp = variable_cell_parameters("entropic_plus_titin")
        cp = replace(
            cp,
            mechanics=replace(cp.mechanics, duration=0.01, save_every=1),
        )
        first = simulate_variable_cell(cp, 23)
        second = simulate_variable_cell(cp, 23)
        np.testing.assert_array_equal(first["coords"], second["coords"])
        np.testing.assert_array_equal(first["length"], second["length"])


if __name__ == "__main__":
    unittest.main()

