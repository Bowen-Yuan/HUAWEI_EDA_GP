"""CLI and reproducible run artifacts for :mod:`hardplacer`."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np

from gpplacer.benchmark.manifest import evaluate_manifest
from gpplacer.io.bookshelf import load_bookshelf
from gpplacer.io.placement import read_placement_centres, write_placement
from gpplacer.oracle.grid import build_density_grid
from gpplacer.solver.steps import project_centres

from .config import HardSolverConfig
from .solver import SafeAnchorSolver
from .visualize import plot_layout, plot_run


_FIELDS = ("iteration", "event", "hpwl", "density_linear", "density_surrogate_squared", "overflow_percent",
           "challenge_square_density",
           "budget_percent", "active_bins", "step", "trust_p50", "trust_p90", "serious_step",
           "challenge_square_budget",
           "bundle_size", "pressure")


def _config(args: argparse.Namespace) -> HardSolverConfig:
    cfg = HardSolverConfig.from_json(args.config) if args.config else HardSolverConfig()
    for name in ("max_iterations", "max_seconds", "target_overflow_percent", "output_root"):
        value = getattr(args, name, None)
        if value is not None: setattr(cfg, name, value)
    return cfg


def _solve(args: argparse.Namespace) -> int:
    cfg = _config(args); db = load_bookshelf(args.aux)
    grid = build_density_grid(db, rho_target=cfg.rho_target, bin_rows=cfg.bin_rows)
    bootstrap = SafeAnchorSolver(db, grid, cfg)
    initial = read_placement_centres(args.initial_placement, db) if args.initial_placement else bootstrap.capacity_seed()
    if not np.array_equal(initial[db.fixed], db.initial_centres[db.fixed]):
        raise ValueError("Initial placement moves fixed objects.")
    initial = project_centres(db, initial)
    run_dir = Path(cfg.output_root) / Path(args.aux).stem / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_dir.mkdir(parents=True, exist_ok=False)
    snapshots: dict[str, np.ndarray] = {}
    records: list[dict[str, object]] = []
    solver = bootstrap
    initial_result = solver.oracle.evaluate(initial)
    snapshots["initial_centres"], snapshots["initial_occupancy"] = initial, initial_result.occupancy

    def observe(record: dict[str, object], centres: np.ndarray, result, trust: np.ndarray) -> None:
        records.append(record)

    result = solver.solve(initial, observe)
    snapshots["anchor_centres"], snapshots["anchor_occupancy"] = result.anchor_centres, result.anchor_oracle.occupancy
    snapshots["final_centres"], snapshots["final_occupancy"] = result.centres, result.oracle.occupancy
    np.savez_compressed(run_dir / "snapshots.npz", **snapshots)
    with (run_dir / "iterations.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_FIELDS); writer.writeheader(); writer.writerows(records)
    solution = run_dir / "solution.pl"; write_placement(solution, db, result.centres)
    payload = {"status": "completed", "source_aux": str(Path(args.aux).resolve()), "config": cfg.to_dict(),
               "elapsed_seconds": result.elapsed_seconds, "iterations": result.iterations,
               "anchor_repairs": result.anchor_repairs, "accepted_pressure_steps": result.accepted_pressure_steps,
               "initial_hpwl": initial_result.hpwl, "initial_overflow_percent": solver._overflow_percent(initial_result),
               "initial_challenge_square_density": solver.challenge_square_density(initial_result),
               "anchor_hpwl": result.anchor_oracle.hpwl, "anchor_overflow_percent": solver._overflow_percent(result.anchor_oracle),
               "anchor_challenge_square_density": solver.challenge_square_density(result.anchor_oracle),
               "final_hpwl": result.oracle.hpwl, "final_overflow_percent": solver._overflow_percent(result.oracle),
               "final_challenge_square_density": solver.challenge_square_density(result.oracle),
               "solution": str(solution.resolve())}
    (run_dir / "metadata.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    if args.plot:
        plot_run(run_dir); plot_layout(run_dir)
    print(json.dumps({"run_dir": str(run_dir), "solution": str(solution), "hpwl": result.oracle.hpwl,
                      "overflow_percent": solver._overflow_percent(result.oracle)}, indent=2))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(prog="hardplacer")
    commands = parser.add_subparsers(required=True)
    solve = commands.add_parser("solve", help="run safe-anchor nonsmooth placement")
    solve.add_argument("--aux", required=True); solve.add_argument("--initial-placement")
    solve.add_argument("--config"); solve.add_argument("--max-iterations", type=int)
    solve.add_argument("--max-seconds", type=float); solve.add_argument("--target-overflow-percent", type=float)
    solve.add_argument("--output-root"); solve.add_argument("--plot", action="store_true")
    solve.set_defaults(handler=_solve)
    visual = commands.add_parser("visualize-run", help="redraw run figures")
    visual.add_argument("--run-dir", required=True)
    visual.set_defaults(handler=lambda a: (plot_run(a.run_dir), plot_layout(a.run_dir), 0)[-1])
    bench = commands.add_parser("benchmark", help="evaluate an eight-case manifest")
    bench.add_argument("--manifest", required=True)
    bench.set_defaults(handler=lambda a: (print(json.dumps(evaluate_manifest(a.manifest), indent=2)), 0)[-1])
    args = parser.parse_args(); return args.handler(args)


if __name__ == "__main__":
    raise SystemExit(main())
