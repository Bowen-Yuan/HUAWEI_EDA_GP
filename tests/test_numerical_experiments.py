from __future__ import annotations

from adaptive_pareto.config import AdaptiveConfig
from adaptive_pareto.density_ablation import solve_density_mode
from adaptive_pareto.two_phase import solve_two_phase
from tests.test_core import _tiny_db


def test_two_phase_experiment_records_both_phases() -> None:
    db = _tiny_db()
    config = AdaptiveConfig()
    result = solve_two_phase(
        db, db.initial_centres, config, hpwl_target=0.0, hpwl_c_rows=1.0,
        hpwl_max_steps=2, joint_steps=2, joint_c_rows=1.0,
        lambda_max=2.0, lambda_power=1.0, max_seconds=20.0, density_sample_every=1,
    )
    assert {row["phase"] for row in result.records} == {"hpwl_rescue", "joint"}


def test_density_direction_modes_use_identical_step_budget() -> None:
    db = _tiny_db()
    config = AdaptiveConfig()
    results = [
        solve_density_mode(db, db.initial_centres, config, mode=mode, c_rows=1.0, steps=2, max_seconds=20.0)
        for mode in ("full_2d", "dominant_axis", "hybrid_axis")
    ]
    assert all(len(result.records) == 2 for result in results)
