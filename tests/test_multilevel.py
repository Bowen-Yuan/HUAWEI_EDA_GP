from __future__ import annotations

from adaptive_pareto.config import AdaptiveConfig
from adaptive_pareto.multilevel import build_hierarchy, generate_candidates, run_portfolio
from tests.test_core import _tiny_db


def test_hierarchy_and_candidate_generation_are_deterministic_on_tiny_case() -> None:
    db = _tiny_db()
    config = AdaptiveConfig()
    config.hierarchy.stop_nodes = 1
    config.hierarchy.max_levels = 2
    config.portfolio.heavy_edge_candidates = 1
    config.portfolio.partition_candidates = 1
    config.portfolio.capacity_candidates = 1
    config.portfolio.connection_candidates = 1
    config.portfolio.spectral_candidates = 1
    levels = build_hierarchy(db, config)
    assert levels
    coarse = levels[-1].coarse
    first = generate_candidates(coarse, config)
    second = generate_candidates(coarse, config)
    assert [item.generator for item in first] == [item.generator for item in second]
    assert len(first) == 5


def test_portfolio_returns_a_finalist_on_tiny_case() -> None:
    db = _tiny_db()
    config = AdaptiveConfig()
    config.runtime.max_seconds = 30.0
    config.hierarchy.stop_nodes = 1
    config.hierarchy.max_levels = 1
    config.hierarchy.level_refine_steps = (1,)
    config.portfolio.heavy_edge_candidates = 1
    config.portfolio.partition_candidates = 1
    config.portfolio.capacity_candidates = 1
    config.portfolio.connection_candidates = 1
    config.portfolio.spectral_candidates = 1
    config.portfolio.halving_steps = (1,)
    config.portfolio.halving_keep = (1,)
    config.joint.steps = 1
    config.rescue.hpwl_enter_norm = 10.0
    config.rescue.hpwl_exit_norm = 10.0
    config.rescue.overflow_enter_percent = 1000.0
    config.rescue.overflow_exit_percent = 1000.0
    result = run_portfolio(db, config)
    assert result.best.records
    assert result.candidate_rows
    assert result.hierarchy_rows
