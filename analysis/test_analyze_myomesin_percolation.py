"""Unit tests for sparse-myomesin percolation metrics."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_myomesin_percolation import graph_components, phase_metrics


def test_two_cycles_and_isolates_are_not_percolated() -> None:
    pairs = {
        (0, 1), (1, 2), (0, 2),
        (3, 4), (4, 5), (3, 5),
    }
    graph = graph_components(8, pairs)
    assert graph["component_sizes"] == [3, 3, 1, 1]
    assert graph["largest_component"] == 3
    assert graph["component_count"] == 4
    assert graph["cycle_rank"] == 2
    assert np.count_nonzero(graph["degree"] == 2) == 6


def test_chain_has_free_ends_and_zero_cycle_rank() -> None:
    graph = graph_components(4, {(0, 1), (1, 2), (2, 3)})
    assert graph["component_sizes"] == [4]
    assert graph["cycle_rank"] == 0
    assert np.count_nonzero(graph["degree"] == 1) == 2


def test_five_phased_components_are_locally_ordered_globally_cancelled() -> None:
    lx = 4.0
    x = np.repeat(np.arange(5) * lx / 5.0, 2)
    components = [[2 * index, 2 * index + 1] for index in range(5)]
    metrics = phase_metrics(
        x, lx, components
    )
    assert np.isclose(metrics["within_component_order"], 1.0)
    assert np.isclose(metrics["m_band_order"], 0.0, atol=1e-12)
    assert np.isclose(
        metrics["between_component_phase_order"], 0.0, atol=1e-12
    )


def test_invalid_noncanonical_pair_is_rejected() -> None:
    try:
        graph_components(3, {(2, 1)})
    except ValueError as error:
        assert "invalid canonical pair" in str(error)
    else:
        raise AssertionError("invalid pair did not raise")
