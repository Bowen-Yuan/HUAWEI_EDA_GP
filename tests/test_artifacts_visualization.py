from __future__ import annotations

from pathlib import Path
import time

from adaptive_pareto.artifacts import load_checkpoint, make_run_dir, write_checkpoint, write_run
from adaptive_pareto.config import AdaptiveConfig
from adaptive_pareto.evaluator import AdaptiveEvaluator
from adaptive_pareto.solver import AdaptiveParetoSolver
from adaptive_pareto.visualize import render_run
from tests.test_core import _tiny_db


def test_artifacts_and_offline_visualisations_are_created() -> None:
    db = _tiny_db()
    config = AdaptiveConfig()
    config.runtime.max_seconds = 20.0
    config.joint.steps = 2
    config.rescue.hpwl_enter_norm = 10.0
    config.rescue.hpwl_exit_norm = 10.0
    config.rescue.overflow_enter_percent = 1000.0
    config.rescue.overflow_exit_percent = 1000.0
    config.visualization.static_formats = ("png",)
    evaluator = AdaptiveEvaluator(db, config.evaluation)
    result = AdaptiveParetoSolver(db, config).run(db.initial_centres)
    root = make_run_dir(Path.cwd() / "tmp" / f"pytest_adaptive_visualization_{time.time_ns()}", "tiny.aux")
    write_run(root, db=db, aux="tiny.aux", initial=db.initial_centres, evaluator=evaluator, config=config, result=result)
    figures = render_run(root)
    assert figures
    assert all(path.exists() for path in figures)
    assert (root / "report.html").exists()
    assert (root / "solutions" / "solution_final.pl").exists()


def test_checkpoint_round_trip_and_resume() -> None:
    db = _tiny_db()
    config = AdaptiveConfig()
    config.runtime.max_seconds = 20.0
    config.runtime.checkpoint_every_steps = 1
    config.joint.steps = 2
    config.rescue.hpwl_enter_norm = config.rescue.hpwl_exit_norm = 10.0
    config.rescue.overflow_enter_percent = config.rescue.overflow_exit_percent = 1000.0
    root = make_run_dir(Path.cwd() / "tmp" / f"pytest_adaptive_resume_{time.time_ns()}", "tiny.aux")
    solver = AdaptiveParetoSolver(db, config)
    solver.run(db.initial_centres, on_checkpoint=lambda state: write_checkpoint(root, state))
    checkpoint = load_checkpoint(root)
    assert checkpoint.iteration == 2
    resumed = AdaptiveParetoSolver(db, config).run(checkpoint.centres, resume=checkpoint)
    assert len(resumed.records) == 2
