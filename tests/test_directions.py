from __future__ import annotations

import numpy as np

from adaptive_pareto.directions import density_direction


def test_hybrid_axis_only_gates_strongly_dominant_cells() -> None:
    gradient = np.array([[9.0, 1.0], [3.0, 2.0], [0.0, 0.0], [8.0, 1.0]])
    preconditioner = np.ones(4)
    trust = np.ones(4)
    movable = np.array([True, True, True, False])
    direction, single_axis = density_direction(
        gradient, preconditioner, trust, movable,
        mode="hybrid_axis", axis_dominance_ratio=3.0,
    )
    assert single_axis.tolist() == [True, False, False, False]
    np.testing.assert_allclose(direction[0], [-1.0, 0.0])
    np.testing.assert_allclose(np.linalg.norm(direction[1]), 1.0)
    np.testing.assert_allclose(direction[2], [0.0, 0.0])
    np.testing.assert_allclose(direction[3], [0.0, 0.0])


def test_dominant_axis_normalises_each_cell_to_a_fair_l2_distance() -> None:
    gradient = np.array([[3.0, 4.0]])
    direction, single_axis = density_direction(
        gradient, np.ones(1), np.ones(1), np.array([True]),
        mode="dominant_axis", axis_dominance_ratio=1.0,
    )
    assert single_axis.tolist() == [True]
    np.testing.assert_allclose(direction, [[0.0, -1.0]])
