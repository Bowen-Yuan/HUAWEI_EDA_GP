"""Small deterministic tests for parser-independent mathematical invariants."""

from __future__ import annotations

import numpy as np

from gpplacer.model.config import SolverConfig
from gpplacer.model.types import PlacementDB, SiteRow
from gpplacer.oracle.grid import build_density_grid
from gpplacer.oracle.evaluation import DualObjectiveOracle
from gpplacer.oracle.objective import ObjectiveOracle
from gpplacer.oracle.hpwl import hpwl_value_gradient
from gpplacer.solver.controller import PlacementSolver
from gpplacer.solver.state import PlacementState


def _tiny_db() -> PlacementDB:
    return PlacementDB(
        source_aux="tiny.aux", node_names=("a", "b", "fixed"), net_names=("n0",),
        width=np.array([2.0, 2.0, 2.0]), height=np.array([2.0, 2.0, 2.0]),
        fixed=np.array([False, False, True]),
        initial_centres=np.array([[0.0, 0.0], [0.0, 0.0], [9.0, 9.0]]),
        orientation=("N", "N", "N"),
        pin_node=np.array([0, 1, 2], dtype=np.int64),
        pin_offset=np.zeros((3, 2)), net_start=np.array([0, 3], dtype=np.int64),
        net_weight=np.array([1.0]), node_net_start=np.array([0, 1, 2, 3], dtype=np.int64),
        node_nets=np.array([0, 0, 0], dtype=np.int64),
        rows=(SiteRow(0.0, 0.0, 10.0, 2.0, 1.0), SiteRow(0.0, 2.0, 10.0, 2.0, 1.0),
              SiteRow(0.0, 4.0, 10.0, 2.0, 1.0), SiteRow(0.0, 6.0, 10.0, 2.0, 1.0),
              SiteRow(0.0, 8.0, 10.0, 2.0, 1.0)),
    )


def test_hpwl_and_tied_subgradient_are_deterministic() -> None:
    centres = np.array([[1.0, 1.0], [4.0, 1.0], [4.0, 3.0]])
    value, gradient, extrema = hpwl_value_gradient(
        centres, np.array([0, 1, 2], dtype=np.int64), np.zeros((3, 2)),
        np.array([0, 3], dtype=np.int64), np.array([1.0]),
    )
    assert value == 5.0
    assert np.isclose(gradient[:, 0].sum(), 0.0)
    assert np.isclose(gradient[:, 1].sum(), 0.0)
    assert extrema.shape == (1, 4)


def test_exact_overflow_uses_target_point_eight() -> None:
    db = _tiny_db()
    grid = build_density_grid(db, rho_target=0.8, bin_rows=1)
    oracle = ObjectiveOracle(db, grid, density_lambda=1.0)
    centres = db.initial_centres.copy()
    centres[:2] = (1.0, 1.0)
    result = oracle.evaluate(centres)
    assert result.density_linear > 0.0
    assert result.density_surrogate_squared > 0.0
    assert result.active_bin_count > 0


def test_strict_density_charges_movable_area_in_zero_capacity_bins() -> None:
    db = _tiny_db()
    grid = build_density_grid(db, rho_target=0.8, bin_rows=1)
    centres = db.initial_centres.copy()
    # The fixed 2x2 object occupies the top-right bin.  The legacy metric
    # intentionally ignores movable occupancy there; the strict protocol does
    # not, which is the invariant needed by the adaptive search controller.
    centres[0] = (9.0, 9.0)
    result = DualObjectiveOracle(db, grid).evaluate(centres)
    assert result.strict.linear > result.legacy.linear
    assert result.strict.zero_capacity_occupancy > 0.0
    assert result.strict.active_bin_count >= result.legacy.active_bin_count


def test_solver_runs_all_initial_components_on_tiny_case() -> None:
    db = _tiny_db()
    config = SolverConfig(max_iterations=3, max_seconds=30.0, snapshot_every=1, threads=1)
    result = PlacementSolver(db, config).solve()
    assert np.isfinite(result.oracle.objective)
    assert result.iterations == 3
    np.testing.assert_allclose(result.centres[db.fixed], db.initial_centres[db.fixed])


def test_limited_memory_bundle_proposes_a_finite_candidate() -> None:
    db = _tiny_db()
    solver = PlacementSolver(db, SolverConfig(coarse_iterations=0, threads=1))
    solver.oracle = ObjectiveOracle(db, solver.grid, density_lambda=1.0)
    solver._oracle_seconds = 0.0
    solver._oracle_calls = 0
    centres = np.array([[1.0, 1.0], [3.0, 1.0], [9.0, 9.0]])
    result = solver._evaluate(centres)
    state = PlacementState(centres, result, centres.copy(), result, proximal_scale=2.0)
    solver._add_cut(state, result)
    candidate, predicted_value = solver._bundle_candidate(state, maximum_displacement=1.0)
    assert np.isfinite(candidate).all()
    assert np.isfinite(predicted_value)
    assert np.max(np.abs(candidate - centres)) <= 1.0 + 1e-12
