from __future__ import annotations

import numpy as np

from gpplacer.model.types import PlacementDB, SiteRow
from gpplacer.oracle.grid import build_density_grid

from hardplacer.config import HardSolverConfig
from hardplacer.solver import SafeAnchorSolver


def tiny_db() -> PlacementDB:
    return PlacementDB(
        source_aux="tiny.aux", node_names=("a", "b", "c", "fixed"), net_names=("n0", "n1"),
        width=np.full(4, 2.0), height=np.full(4, 2.0), fixed=np.array([False, False, False, True]),
        initial_centres=np.array([[1., 1.], [1., 1.], [3., 1.], [9., 9.]]), orientation=("N",) * 4,
        pin_node=np.array([0, 1, 3, 1, 2], dtype=np.int64), pin_offset=np.zeros((5, 2)),
        net_start=np.array([0, 3, 5], dtype=np.int64), net_weight=np.ones(2),
        node_net_start=np.array([0, 1, 3, 4, 5]), node_nets=np.array([0, 0, 1, 1, 0]),
        rows=tuple(SiteRow(0., y, 10., 2., 1.) for y in (0., 2., 4., 6., 8.)),
    )


def test_hard_budget_never_rises_after_anchor() -> None:
    db = tiny_db(); grid = build_density_grid(db, rho_target=.8, bin_rows=1)
    config = HardSolverConfig(max_iterations=10, anchor_repair_rounds=4, target_overflow_percent=10.0)
    result = SafeAnchorSolver(db, grid, config).solve(db.initial_centres)
    assert all(float(row["density_linear"]) <= float(row["budget_percent"]) / 100 * grid.available_area.sum() + 1e-8
               for row in result.records)
    # A supplied 10% cap remains a cap even when the constructed anchor is
    # better; it must not be implicitly tightened to the anchor's value.
    assert all(abs(float(row["budget_percent"]) - 10.0) < 1e-8 for row in result.records)


def test_optional_square_diagnostic_can_be_enforced() -> None:
    db = tiny_db(); grid = build_density_grid(db, rho_target=.8, bin_rows=1)
    config = HardSolverConfig(max_iterations=10, anchor_repair_rounds=4,
                              target_overflow_percent=10.0,
                              enforce_challenge_square_budget=True)
    result = SafeAnchorSolver(db, grid, config).solve(db.initial_centres)
    assert all(float(row["challenge_square_density"]) <= float(row["challenge_square_budget"]) + 1e-8
               for row in result.records)


def test_trust_region_is_componentwise_and_bounded() -> None:
    db = tiny_db(); grid = build_density_grid(db, rho_target=.8, bin_rows=1)
    config = HardSolverConfig(trust_initial_rows=1., trust_min_rows=.25, trust_max_rows=2.)
    solver = SafeAnchorSolver(db, grid, config)
    trust = np.full(db.node_count, db.row_height)
    displacement = np.zeros((db.node_count, 2)); displacement[0, 0] = 1.
    solver._update_trust(trust, displacement, False)
    assert trust[0] == .5 * db.row_height
    assert trust[1] == db.row_height
    solver._update_trust(trust, displacement, True)
    assert .25 * db.row_height <= trust[0] <= 2. * db.row_height


def test_fixed_cells_and_candidates_remain_finite() -> None:
    db = tiny_db(); grid = build_density_grid(db, rho_target=.8, bin_rows=1)
    result = SafeAnchorSolver(db, grid, HardSolverConfig(max_iterations=3)).solve(db.initial_centres)
    assert np.isfinite(result.centres).all()
    np.testing.assert_allclose(result.centres[db.fixed], db.initial_centres[db.fixed])
