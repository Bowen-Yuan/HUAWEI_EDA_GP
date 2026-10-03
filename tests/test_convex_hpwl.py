from __future__ import annotations

import numpy as np

from adaptive_pareto.convex_hpwl import solve_convex_hpwl
from tests.test_core import _tiny_db


def test_convex_hpwl_solver_records_best_so_far_monotonically() -> None:
    db = _tiny_db()
    result = solve_convex_hpwl(db, db.initial_centres, c_rows=1.0, max_steps=8, max_seconds=10.0)
    assert result.records
    best = np.asarray([float(row["best_hpwl"]) for row in result.records])
    assert np.all(np.diff(best) <= 1.0e-12)
    assert result.best_hpwl <= result.records[0]["current_hpwl"] + 1.0e-12
