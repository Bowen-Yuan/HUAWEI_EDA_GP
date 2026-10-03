from __future__ import annotations

import numpy as np

from adaptive_pareto.config import AdaptiveConfig
from adaptive_pareto.solver import AdaptiveParetoSolver
from tests.test_core import _tiny_db


def test_adaptive_solver_runs_and_records_dual_metrics() -> None:
    db = _tiny_db()
    config = AdaptiveConfig()
    config.runtime.max_seconds = 20.0
    config.joint.steps = 3
    config.rescue.hpwl_enter_norm = 10.0
    config.rescue.hpwl_exit_norm = 10.0
    config.rescue.overflow_enter_percent = 1000.0
    config.rescue.overflow_exit_percent = 1000.0
    result = AdaptiveParetoSolver(db, config).run(db.initial_centres)
    assert result.records
    assert result.status == "completed"
    assert np.isfinite(result.oracle.hpwl)
    assert "strict_overflow_percent" in result.records[0]
    assert "legacy_overflow_percent" in result.records[0]
