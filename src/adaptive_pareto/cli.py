"""Command-line entry points for adaptive Pareto experiments."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

from gpplacer.io.bookshelf import load_bookshelf
from gpplacer.io.placement import read_placement_centres
from gpplacer.multilevel.initialization import capacity_aware_seed

from adaptive_pareto.artifacts import (
    initialize_run, load_checkpoint, make_checkpoint_callback, make_run_dir, write_run, write_table,
)
from adaptive_pareto.config import AdaptiveConfig
from adaptive_pareto.compare import compare_runs
from adaptive_pareto.convex_hpwl import solve_convex_hpwl, write_convex_hpwl_run
from adaptive_pareto.density_ablation import run_density_ablation
from adaptive_pareto.evaluator import AdaptiveEvaluator
from adaptive_pareto.solver import AdaptiveParetoSolver
from adaptive_pareto.multilevel import run_portfolio
from adaptive_pareto.visualize import render_run
from adaptive_pareto.two_phase import solve_two_phase, write_two_phase_run


def _config(path: str | None, output_root: str | None) -> AdaptiveConfig:
    config = AdaptiveConfig.from_json(path) if path else AdaptiveConfig()
    if output_root is not None:
        config.output_root = output_root
    config.validate()
    return config


def _solve(args: argparse.Namespace) -> int:
    resume_root = Path(args.resume_run).resolve() if args.resume_run else None
    resume = load_checkpoint(resume_root) if resume_root is not None else None
    if resume_root is not None:
        metadata = json.loads((resume_root / "metadata.json").read_text(encoding="utf-8"))
        aux = args.aux or metadata["source_aux"]
        config = _config(args.config, args.output_root) if args.config else AdaptiveConfig.from_mapping(metadata["config"])
        if args.output_root is not None:
            config.output_root = args.output_root
        root = resume_root
    else:
        if not args.aux:
            raise ValueError("--aux is required unless --resume-run is supplied.")
        aux = args.aux
        config = _config(args.config, args.output_root)
        root = make_run_dir(config.output_root, aux)
    db = load_bookshelf(aux)
    evaluator = AdaptiveEvaluator(db, config.evaluation)
    initialize_run(
        root, db=db, aux=aux, config=config,
        input_placement=str(Path(args.placement).resolve()) if args.placement else None,
    )
    checkpoint_callback = make_checkpoint_callback(root, evaluator, config)
    portfolio_rows = None
    hierarchy_rows = None
    portfolio_elapsed = None
    if resume is not None:
        initial = resume.centres
        result = AdaptiveParetoSolver(db, config).run(
            initial, resume=resume, on_checkpoint=checkpoint_callback,
        )
    elif args.placement:
        initial = read_placement_centres(args.placement, db)
        result = AdaptiveParetoSolver(db, config).run(
            initial, on_checkpoint=checkpoint_callback,
        )
    elif config.portfolio.enabled:
        portfolio = run_portfolio(db, config)
        result = portfolio.best
        input_points = [point for point in result.archive.points if point.phase == "input"]
        initial = input_points[0].centres if input_points else result.centres
        portfolio_rows = portfolio.candidate_rows
        hierarchy_rows = portfolio.hierarchy_rows
        portfolio_elapsed = portfolio.elapsed_seconds
    else:
        initial = capacity_aware_seed(db, evaluator.grid, config.runtime.seed)
        result = AdaptiveParetoSolver(db, config).run(
            initial, on_checkpoint=checkpoint_callback,
        )
    write_run(root, db=db, aux=aux, initial=initial, evaluator=evaluator, config=config, result=result,
              extra_metadata={
                  "input_placement": str(Path(args.placement).resolve()) if args.placement else None,
                  "portfolio_elapsed_seconds": portfolio_elapsed,
              })
    if portfolio_rows is not None:
        write_table(root / "candidates.csv", portfolio_rows)
    if hierarchy_rows is not None:
        write_table(root / "hierarchy.csv", hierarchy_rows)
    if config.visualization.enabled:
        rendering_started = time.perf_counter()
        render_run(root)
        metadata_path = root / "metadata.json"
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        metadata["visualization_elapsed_seconds"] = time.perf_counter() - rendering_started
        metadata_path.write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({
        "run_dir": str(root), "status": result.status, "hpwl": result.oracle.hpwl,
        "legacy_overflow_percent": evaluator.overflow_percent(result.oracle, strict=False),
        "strict_overflow_percent": evaluator.overflow_percent(result.oracle, strict=True),
    }, indent=2))
    return 0


def _evaluate(args: argparse.Namespace) -> int:
    config = _config(args.config, None)
    db = load_bookshelf(args.aux)
    evaluator = AdaptiveEvaluator(db, config.evaluation)
    centres = read_placement_centres(args.placement, db)
    result = evaluator.evaluate(centres)
    print(json.dumps({
        "hpwl": result.hpwl,
        "normalized_hpwl": evaluator.normalized_hpwl(result),
        "legacy_overflow_percent": evaluator.overflow_percent(result, strict=False),
        "strict_overflow_percent": evaluator.overflow_percent(result, strict=True),
        "zero_capacity_occupancy": result.strict.zero_capacity_occupancy,
    }, indent=2))
    return 0


def _visualize(args: argparse.Namespace) -> int:
    written = render_run(args.run_dir)
    print(json.dumps({"figures": [str(path) for path in written], "report": str(Path(args.run_dir) / "report.html")}, indent=2))
    return 0


def _compare(args: argparse.Namespace) -> int:
    output = compare_runs(args.run_dir, args.output)
    print(json.dumps({"comparison_dir": str(output)}, indent=2))
    return 0


def _hpwl_baseline(args: argparse.Namespace) -> int:
    config = _config(args.config, None)
    db = load_bookshelf(args.aux)
    initial = read_placement_centres(args.placement, db)
    evaluator = AdaptiveEvaluator(db, config.evaluation)
    result = solve_convex_hpwl(
        db, initial, c_rows=args.c_rows, max_steps=args.max_steps, max_seconds=args.max_seconds,
    )
    root = write_convex_hpwl_run(
        args.output_root, db=db, initial=initial, result=result, evaluator=evaluator,
        c_rows=args.c_rows, max_steps=args.max_steps, max_seconds=args.max_seconds,
    )
    print(json.dumps({"run_dir": str(root), "status": result.status, "best_hpwl": result.best_hpwl}, indent=2))
    return 0


def _two_phase(args: argparse.Namespace) -> int:
    config = _config(args.config, None)
    db = load_bookshelf(args.aux)
    initial = read_placement_centres(args.placement, db)
    parameters = {
        "hpwl_target": args.hpwl_target, "hpwl_c_rows": args.hpwl_c_rows,
        "hpwl_max_steps": args.hpwl_max_steps, "joint_steps": args.joint_steps,
        "joint_c_rows": args.joint_c_rows, "lambda_max": args.lambda_max,
        "lambda_power": args.lambda_power, "lambda_start": args.lambda_start,
        "density_metric": args.density_metric,
        "max_seconds": args.max_seconds,
    }
    result = solve_two_phase(db, initial, config, **parameters)
    root = write_two_phase_run(
        args.output_root, db=db, initial=initial, config=config, result=result, parameters=parameters,
    )
    print(json.dumps({"run_dir": str(root), "status": result.status, "transition_hpwl": result.transition_hpwl}, indent=2))
    return 0


def _density_ablation(args: argparse.Namespace) -> int:
    config = _config(args.config, None)
    db = load_bookshelf(args.aux)
    initial = read_placement_centres(args.placement, db)
    root = run_density_ablation(
        args.output_root, db=db, initial=initial, config=config, c_rows=args.c_rows,
        steps=args.steps, max_seconds_per_mode=args.max_seconds_per_mode,
        gradient_objective=args.gradient_objective, modes=tuple(args.modes),
    )
    print(json.dumps({"run_dir": str(root)}, indent=2))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(prog="adaptive-pareto")
    commands = parser.add_subparsers(required=True)
    solve = commands.add_parser("solve")
    solve.add_argument("--aux")
    solve.add_argument("--placement")
    solve.add_argument("--config")
    solve.add_argument("--output-root")
    solve.add_argument("--resume-run", help="Existing run directory containing checkpoints/latest.*")
    solve.set_defaults(handler=_solve)
    evaluate = commands.add_parser("evaluate")
    evaluate.add_argument("--aux", required=True)
    evaluate.add_argument("--placement", required=True)
    evaluate.add_argument("--config")
    evaluate.set_defaults(handler=_evaluate)
    visualize = commands.add_parser("visualize-run")
    visualize.add_argument("--run-dir", required=True)
    visualize.set_defaults(handler=_visualize)
    compare = commands.add_parser("compare-runs")
    compare.add_argument("--run-dir", action="append", required=True, help="Repeat for every completed run")
    compare.add_argument("--output", required=True)
    compare.set_defaults(handler=_compare)
    baseline = commands.add_parser("hpwl-baseline")
    baseline.add_argument("--aux", required=True)
    baseline.add_argument("--placement", required=True)
    baseline.add_argument("--config")
    baseline.add_argument("--c-rows", type=float, default=8.0)
    baseline.add_argument("--max-steps", type=int, default=5000)
    baseline.add_argument("--max-seconds", type=float, default=600.0)
    baseline.add_argument("--output-root", default="adaptive_pareto_soft/results_hpwl_convex")
    baseline.set_defaults(handler=_hpwl_baseline)
    two_phase = commands.add_parser("two-phase")
    two_phase.add_argument("--aux", required=True)
    two_phase.add_argument("--placement", required=True)
    two_phase.add_argument("--config")
    two_phase.add_argument("--hpwl-target", type=float, default=40.0e6)
    two_phase.add_argument("--hpwl-c-rows", type=float, default=128.0)
    two_phase.add_argument("--hpwl-max-steps", type=int, default=15000)
    two_phase.add_argument("--joint-steps", type=int, default=3000)
    two_phase.add_argument("--joint-c-rows", type=float, default=32.0)
    two_phase.add_argument("--lambda-max", type=float, default=16.0)
    two_phase.add_argument("--lambda-start", type=float, default=0.0)
    two_phase.add_argument("--density-metric", choices=("strict", "legacy"), default="strict")
    two_phase.add_argument("--lambda-power", type=float, default=2.0)
    two_phase.add_argument("--max-seconds", type=float, default=2400.0)
    two_phase.add_argument("--output-root", default="adaptive_pareto_soft/results_two_phase")
    two_phase.set_defaults(handler=_two_phase)
    density = commands.add_parser("density-ablation")
    density.add_argument("--aux", required=True)
    density.add_argument("--placement", required=True)
    density.add_argument("--config")
    density.add_argument("--c-rows", type=float, default=16.0)
    density.add_argument("--steps", type=int, default=1000)
    density.add_argument("--max-seconds-per-mode", type=float, default=600.0)
    density.add_argument("--gradient-objective", choices=("linear", "squared"), default="squared")
    density.add_argument("--modes", nargs="+", choices=("full_2d", "dominant_axis", "hybrid_axis"),
                         default=("full_2d", "dominant_axis", "hybrid_axis"))
    density.add_argument("--output-root", default="adaptive_pareto_soft/results_density_ablation")
    density.set_defaults(handler=_density_ablation)
    args = parser.parse_args()
    return args.handler(args)


if __name__ == "__main__":
    raise SystemExit(main())
