"""Command-line entry points for solving, inspection, plotting, and benchmark review."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time
from dataclasses import asdict

from gpplacer.benchmark.manifest import evaluate_manifest
from gpplacer.io.bookshelf import load_bookshelf
from gpplacer.io.placement import read_placement_centres, write_placement
from gpplacer.experimental.anchored_search import (
    AnchoredSearch, AnchoredSearchConfig, load_anchor, select_lowest_overflow_solution,
)
from gpplacer.experimental.density_repair import DensityRepair, DensityRepairConfig
from gpplacer.metrics.logging import RunLogger
from gpplacer.model.config import SolverConfig
from gpplacer.oracle.grid import build_density_grid
from gpplacer.oracle.objective import ObjectiveOracle
from gpplacer.official import run_official_evaluation
from gpplacer.solver.controller import PlacementSolver
from gpplacer.solver.state import PlacementState
from gpplacer.visualization.animation import animate_run
from gpplacer.visualization.plots import plot_placement, plot_run


def _config(args: argparse.Namespace) -> SolverConfig:
    config = SolverConfig.from_json(args.config) if args.config else SolverConfig()
    if getattr(args, "threads", None) is not None:
        config.threads = args.threads
    if getattr(args, "max_iterations", None) is not None:
        config.max_iterations = args.max_iterations
    return config


def _solve(args: argparse.Namespace) -> int:
    config = _config(args)
    db = load_bookshelf(args.aux)
    logger = RunLogger(config.output_root, db, config)
    solver = PlacementSolver(db, config)
    final_state: PlacementState | None = None

    def observe(record: dict[str, object], state: PlacementState) -> None:
        nonlocal final_state
        logger.observe(record, state)
        final_state = state

    try:
        result = solver.solve(observe)
        output = Path(args.output) if args.output else logger.run_dir / "solution.pl"
        write_placement(output, db, result.centres)
        if final_state is not None:
            final_state.centres, final_state.oracle = result.centres, result.oracle
        logger.finish(
            status="completed",
            summary={
                "elapsed_seconds": result.elapsed_seconds, "iterations": result.iterations,
                "recoveries": result.recoveries, "density_lambda": result.density_lambda,
                "best_hpwl": result.oracle.hpwl,
                "best_density_linear": result.oracle.density_linear,
                "best_density_surrogate_squared": result.oracle.density_surrogate_squared,
                "best_objective": result.oracle.objective,
                "solution": str(output.resolve()),
            }, state=final_state,
        )
    except Exception as exc:
        logger.finish(status="failed", summary={"error": repr(exc)}, state=final_state)
        raise
    print(json.dumps({"run_dir": str(logger.run_dir), "solution": str(output),
                      "hpwl": result.oracle.hpwl, "overflow": result.oracle.density_linear}, indent=2))
    if args.plot:
        plot_run(logger.run_dir)
        plot_placement(args.aux, output, logger.run_dir / "figures",
                       rho_target=config.rho_target, bin_rows=config.bin_rows)
    return 0


def _inspect(args: argparse.Namespace) -> int:
    db = load_bookshelf(args.aux)
    print(json.dumps({
        "nodes": db.node_count, "fixed": int(db.fixed.sum()), "nets": db.net_count,
        "pins": db.pin_count, "rows": len(db.rows), "core_bounds": db.core_bounds,
        "row_height": db.row_height,
    }, indent=2))
    return 0


def _evaluate(args: argparse.Namespace) -> int:
    """Evaluate an existing placement without running an optimisation step."""
    db = load_bookshelf(args.aux)
    centres = read_placement_centres(args.placement, db)
    if not (centres[db.fixed] == db.initial_centres[db.fixed]).all():
        raise ValueError("The supplied placement moves one or more fixed objects.")
    grid = build_density_grid(db, rho_target=args.rho_target, bin_rows=args.bin_rows)
    result = ObjectiveOracle(db, grid, args.density_lambda).evaluate(centres)
    print(json.dumps({
        "hpwl": result.hpwl, "density_linear": result.density_linear,
        "density_surrogate_squared": result.density_surrogate_squared, "objective": result.objective,
        "active_bins": result.active_bin_count,
    }, indent=2))
    return 0


def _official_evaluate(args: argparse.Namespace) -> int:
    """Run the original ISPD scripts; no project-local metric is substituted."""
    try:
        result = run_official_evaluation(
            args.aux, args.placement, density_target=args.density_target, perl=args.perl,
            hpwl_script=args.hpwl_script, density_script=args.density_script,
        )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(f"official-evaluate: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result.to_dict(), indent=2))
    return 0


def _visualize_run(args: argparse.Namespace) -> int:
    for path in plot_run(args.run_dir):
        print(path.with_suffix(".png"))
    return 0


def _visualize_placement(args: argparse.Namespace) -> int:
    for path in plot_placement(args.aux, args.placement, args.output_dir,
                               rho_target=args.rho_target, bin_rows=args.bin_rows):
        print(path.with_suffix(".png"))
    return 0


def _animate_run(args: argparse.Namespace) -> int:
    output = animate_run(
        args.run_dir, args.output, fps=args.fps, max_points=args.max_points, dpi=args.dpi,
    )
    print(output)
    return 0


def _benchmark(args: argparse.Namespace) -> int:
    print(json.dumps(evaluate_manifest(args.manifest), indent=2))
    return 0


def _anchor_search(args: argparse.Namespace) -> int:
    """Run the isolated overflow-safe anchor-and-pressure experiment."""
    db = load_bookshelf(args.aux)
    anchor_path = Path(args.anchor_placement) if args.anchor_placement else select_lowest_overflow_solution(
        args.anchor_summary, args.overflow_percent_limit,
    )
    anchor_centres = load_anchor(anchor_path, db)
    if not (anchor_centres[db.fixed] == db.initial_centres[db.fixed]).all():
        raise ValueError("The anchor placement moves one or more fixed objects.")
    base_config = SolverConfig(
        rho_target=args.rho_target, bin_rows=args.bin_rows, max_iterations=args.max_iterations,
        snapshot_every=args.snapshot_every, output_root=args.output_root,
    )
    search_config = AnchoredSearchConfig(
        max_iterations=args.max_iterations, initial_step=args.initial_step,
        max_displacement_rows=args.max_displacement_rows, stall_limit=args.stall_limit,
        pressure_levels=args.pressure_levels, pressure_steps=args.pressure_steps,
        pressure_increment_area=args.pressure_increment_area,
        exchange_pair_limit=args.exchange_pair_limit,
        overflow_percent_limit=args.overflow_percent_limit,
    )
    grid = build_density_grid(db, rho_target=base_config.rho_target, bin_rows=base_config.bin_rows)
    search = AnchoredSearch(db, grid, search_config)
    logger = RunLogger(base_config.output_root, db, base_config)
    initial_oracle = search.oracle.evaluate(anchor_centres)
    final_state = PlacementState(
        centres=anchor_centres.copy(), oracle=initial_oracle,
        best_centres=anchor_centres.copy(), best_oracle=initial_oracle,
    )
    logger.save_snapshot("anchor", final_state)
    start = time.perf_counter()

    def observe(record: dict[str, object], centres, oracle) -> None:
        final_state.centres = centres.copy()
        final_state.oracle = oracle
        final_state.best_centres = centres.copy()
        final_state.best_oracle = oracle
        logger.observe(record, final_state)

    try:
        result = search.run(anchor_centres, observe)
        output = logger.run_dir / "solution.pl"
        write_placement(output, db, result.centres)
        elapsed = time.perf_counter() - start
        logger.finish(
            status="completed",
            summary={
                "elapsed_seconds": elapsed,
                "iterations": args.max_iterations,
                "anchor_placement": str(anchor_path.resolve()),
                "anchor_density_linear": initial_oracle.density_linear,
                "best_hpwl": result.oracle.hpwl,
                "best_density_linear": result.oracle.density_linear,
                "best_density_surrogate_squared": result.oracle.density_surrogate_squared,
                "pressure_probes": result.pressure_probes,
                "accepted_exchanges": result.accepted_exchanges,
                "anchored_search_config": search.config_dict(search_config),
                "solution": str(output.resolve()),
            },
            state=final_state,
        )
    except Exception as exc:
        logger.finish(status="failed", summary={"error": repr(exc)}, state=final_state)
        raise
    print(json.dumps({
        "run_dir": str(logger.run_dir), "anchor": str(anchor_path),
        "hpwl": result.oracle.hpwl, "overflow": result.oracle.density_linear,
        "overflow_percent": 100.0 * result.oracle.density_linear / float(grid.available_area.sum()),
        "accepted_exchanges": result.accepted_exchanges,
    }, indent=2))
    if args.plot:
        plot_run(logger.run_dir)
        plot_placement(args.aux, output, logger.run_dir / "figures",
                       rho_target=args.rho_target, bin_rows=args.bin_rows)
    return 0


def _repair_density(args: argparse.Namespace) -> int:
    """Run the isolated large-step density repair before anchored compression."""
    db = load_bookshelf(args.aux)
    input_centres = read_placement_centres(args.placement, db)
    if not (input_centres[db.fixed] == db.initial_centres[db.fixed]).all():
        raise ValueError("The supplied placement moves one or more fixed objects.")
    base_config = SolverConfig(
        rho_target=args.rho_target, bin_rows=args.bin_rows, max_iterations=args.max_iterations,
        snapshot_every=args.snapshot_every, output_root=args.output_root,
    )
    repair_config = DensityRepairConfig(
        max_iterations=args.max_iterations, target_overflow_percent=args.target_overflow_percent,
        initial_step=args.initial_step, max_displacement_rows=args.max_displacement_rows,
        evacuation_rounds=args.evacuation_rounds,
        evacuation_max_moves_per_round=args.evacuation_max_moves_per_round,
        evacuation_max_radius=args.evacuation_max_radius,
    )
    grid = build_density_grid(db, rho_target=base_config.rho_target, bin_rows=base_config.bin_rows)
    repair = DensityRepair(db, grid, repair_config)
    logger = RunLogger(base_config.output_root, db, base_config)
    initial_oracle = repair.oracle.evaluate(input_centres)
    final_state = PlacementState(
        centres=input_centres.copy(), oracle=initial_oracle,
        best_centres=input_centres.copy(), best_oracle=initial_oracle,
    )
    logger.save_snapshot("input", final_state)
    start = time.perf_counter()

    def observe(record: dict[str, object], centres, oracle) -> None:
        final_state.centres, final_state.oracle = centres.copy(), oracle
        logger.observe(record, final_state)

    try:
        result = repair.run(input_centres, observe)
        final_state.centres, final_state.oracle = result.centres.copy(), result.oracle
        final_state.best_centres, final_state.best_oracle = result.centres.copy(), result.oracle
        output = logger.run_dir / "solution.pl"
        write_placement(output, db, result.centres)
        logger.finish(
            status="completed",
            summary={
                "elapsed_seconds": time.perf_counter() - start,
                "iterations": result.iterations,
                "input_placement": str(Path(args.placement).resolve()),
                "target_overflow_percent": args.target_overflow_percent,
                "reached_target": result.reached_target,
                "best_hpwl": result.oracle.hpwl,
                "best_density_linear": result.oracle.density_linear,
                "best_density_surrogate_squared": result.oracle.density_surrogate_squared,
                "density_repair_config": asdict(repair_config),
                "solution": str(output.resolve()),
            },
            state=final_state,
        )
    except Exception as exc:
        logger.finish(status="failed", summary={"error": repr(exc)}, state=final_state)
        raise
    print(json.dumps({
        "run_dir": str(logger.run_dir), "hpwl": result.oracle.hpwl,
        "overflow": result.oracle.density_linear,
        "overflow_percent": 100.0 * result.oracle.density_linear / float(grid.available_area.sum()),
        "reached_target": result.reached_target,
    }, indent=2))
    if args.plot:
        plot_run(logger.run_dir)
        plot_placement(args.aux, output, logger.run_dir / "figures",
                       rho_target=args.rho_target, bin_rows=args.bin_rows)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(prog="gpplacer")
    commands = parser.add_subparsers(required=True)
    solve = commands.add_parser("solve", help="run four-stage nonsmooth placement")
    solve.add_argument("--aux", required=True)
    solve.add_argument("--output")
    solve.add_argument("--config")
    solve.add_argument("--threads", type=int)
    solve.add_argument("--max-iterations", type=int)
    solve.add_argument("--plot", action="store_true", help="create offline figures after solve")
    solve.set_defaults(handler=_solve)
    inspect = commands.add_parser("inspect", help="print verified Bookshelf case metadata")
    inspect.add_argument("--aux", required=True)
    inspect.set_defaults(handler=_inspect)
    evaluate = commands.add_parser("evaluate", help="calculate exact metrics for an existing placement")
    evaluate.add_argument("--aux", required=True)
    evaluate.add_argument("--placement", required=True)
    evaluate.add_argument("--rho-target", type=float, default=0.8)
    evaluate.add_argument("--bin-rows", type=int, default=8)
    evaluate.add_argument("--density-lambda", type=float, default=0.0)
    evaluate.set_defaults(handler=_evaluate)
    official = commands.add_parser(
        "official-evaluate", help="run the unmodified ISPD HPWL and density scripts",
    )
    official.add_argument("--aux", required=True)
    official.add_argument("--placement", required=True)
    official.add_argument("--density-target", required=True, type=float)
    official.add_argument("--perl", default="perl", help="Perl executable used to launch the official scripts")
    official.add_argument("--hpwl-script", help="override dataset/hpwl.pl/hpwl.pl")
    official.add_argument("--density-script", help="override the vendored ISPD 2006 density script")
    official.set_defaults(handler=_official_evaluate)
    run_plot = commands.add_parser("visualize-run", help="plot saved scalar run logs")
    run_plot.add_argument("--run-dir", required=True)
    run_plot.set_defaults(handler=_visualize_run)
    animate = commands.add_parser("animate-run", help="render a placement/density/state GIF or MP4 from snapshots")
    animate.add_argument("--run-dir", required=True)
    animate.add_argument("--output", help="animation path; defaults to figures/optimization_animation.gif")
    animate.add_argument("--fps", type=int, default=5)
    animate.add_argument("--max-points", type=int, default=100_000,
                         help="maximum deterministically sampled movable cells per frame")
    animate.add_argument("--dpi", type=int, default=120)
    animate.set_defaults(handler=_animate_run)
    placement_plot = commands.add_parser("visualize-placement", help="plot placement and density")
    placement_plot.add_argument("--aux", required=True)
    placement_plot.add_argument("--placement", required=True)
    placement_plot.add_argument("--output-dir", required=True)
    placement_plot.add_argument("--rho-target", type=float, default=0.8)
    placement_plot.add_argument("--bin-rows", type=int, default=8)
    placement_plot.set_defaults(handler=_visualize_placement)
    benchmark = commands.add_parser("benchmark", help="check a completed eight-case baseline manifest")
    benchmark.add_argument("--manifest", required=True)
    benchmark.set_defaults(handler=_benchmark)
    anchor_search = commands.add_parser(
        "anchor-search", help="experimental overflow-safe HPWL compression from a placement anchor",
    )
    anchor_search.add_argument("--aux", required=True)
    anchor_search.add_argument("--anchor-placement", help="complete Bookshelf .pl used as the safe anchor")
    anchor_search.add_argument(
        "--anchor-summary", default="runs/adaptec1/adaptec1_valid_results_summary.csv",
        help="CSV used to choose the lowest-overflow anchor when --anchor-placement is omitted",
    )
    anchor_search.add_argument("--rho-target", type=float, default=0.8)
    anchor_search.add_argument("--bin-rows", type=int, default=8)
    anchor_search.add_argument("--output-root", default="runs_anchor")
    anchor_search.add_argument("--max-iterations", type=int, default=200)
    anchor_search.add_argument("--snapshot-every", type=int, default=20)
    anchor_search.add_argument("--initial-step", type=float, default=1.0)
    anchor_search.add_argument("--max-displacement-rows", type=float, default=4.0)
    anchor_search.add_argument("--stall-limit", type=int, default=8)
    anchor_search.add_argument("--pressure-levels", type=int, default=3)
    anchor_search.add_argument("--pressure-steps", type=int, default=6)
    anchor_search.add_argument("--pressure-increment-area", type=float, default=10000.0)
    anchor_search.add_argument("--exchange-pair-limit", type=int, default=24)
    anchor_search.add_argument(
        "--overflow-percent-limit", type=float,
        help="hard cap for final safe anchors; auto-selects the lowest-HPWL anchor within this cap",
    )
    anchor_search.add_argument("--plot", action="store_true")
    anchor_search.set_defaults(handler=_anchor_search)
    repair = commands.add_parser(
        "repair-density", help="experimental large-step density-only repair from an existing placement",
    )
    repair.add_argument("--aux", required=True)
    repair.add_argument("--placement", required=True)
    repair.add_argument("--rho-target", type=float, default=0.8)
    repair.add_argument("--bin-rows", type=int, default=8)
    repair.add_argument("--target-overflow-percent", type=float, default=6.0)
    repair.add_argument("--output-root", default="runs_density_repair")
    repair.add_argument("--max-iterations", type=int, default=400)
    repair.add_argument("--snapshot-every", type=int, default=20)
    repair.add_argument("--initial-step", type=float, default=8.0)
    repair.add_argument("--max-displacement-rows", type=float, default=64.0)
    repair.add_argument("--evacuation-rounds", type=int, default=20)
    repair.add_argument("--evacuation-max-moves-per-round", type=int, default=10_000)
    repair.add_argument("--evacuation-max-radius", type=int, default=112)
    repair.add_argument("--plot", action="store_true")
    repair.set_defaults(handler=_repair_density)
    args = parser.parse_args()
    return args.handler(args)


if __name__ == "__main__":
    raise SystemExit(main())
