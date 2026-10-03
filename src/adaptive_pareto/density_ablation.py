"""Controlled density-direction ablation for full-2D, single-axis, and hybrid descent."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
import time

import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

from gpplacer.io.placement import write_placement
from gpplacer.model.types import PlacementDB
from gpplacer.oracle.density import strict_density_squared_gradient
from gpplacer.solver.steps import build_preconditioner, project_centres

from adaptive_pareto.config import AdaptiveConfig
from adaptive_pareto.directions import density_direction
from adaptive_pareto.evaluator import AdaptiveEvaluator


@dataclass(slots=True)
class DensityModeResult:
    mode: str
    centres: np.ndarray
    records: list[dict[str, float | int | str]]
    elapsed_seconds: float


def solve_density_mode(
    db: PlacementDB, initial: np.ndarray, config: AdaptiveConfig, *, mode: str,
    c_rows: float, steps: int, max_seconds: float, gradient_objective: str = "linear",
) -> DensityModeResult:
    evaluator = AdaptiveEvaluator(db, config.evaluation)
    preconditioner = build_preconditioner(db)
    current = project_centres(db, initial)
    result = evaluator.evaluate(current)
    trust = np.full(db.node_count, c_rows * db.row_height)
    trust[db.fixed] = 0.0
    records: list[dict[str, float | int | str]] = []
    started = time.perf_counter()
    for iteration in range(1, steps + 1):
        if time.perf_counter() - started >= max_seconds:
            break
        if gradient_objective == "squared":
            grid = evaluator.grid
            gradient = strict_density_squared_gradient(
                current, db.width, db.height, db.movable, grid.x_edges, grid.y_edges,
                grid.available_area, grid.rho_target, result.strict.occupancy,
            )
            current_objective = result.strict.squared + result.strict.zero_capacity_occupancy
        else:
            gradient = result.strict.gradient
            current_objective = result.strict.linear
        direction, single_axis = density_direction(
            gradient, preconditioner, trust, db.movable,
            mode=mode, axis_dominance_ratio=config.direction.axis_dominance_ratio,
        )
        step = c_rows * db.row_height / np.sqrt(iteration)
        candidate = project_centres(db, current + step * direction)
        trial = evaluator.evaluate(candidate)
        trial_objective = (trial.strict.squared + trial.strict.zero_capacity_occupancy
                           if gradient_objective == "squared" else trial.strict.linear)
        accepted = trial_objective < current_objective - 1.0e-8
        if accepted:
            current, result = candidate, trial
        records.append({
            "mode": mode, "iteration": iteration, "elapsed_seconds": time.perf_counter() - started,
            "hpwl": result.hpwl, "strict_density_linear": result.strict.linear,
            "density_descent_objective": (result.strict.squared + result.strict.zero_capacity_occupancy
                                           if gradient_objective == "squared" else result.strict.linear),
            "strict_overflow_percent": evaluator.overflow_percent(result, strict=True),
            "legacy_overflow_percent": evaluator.overflow_percent(result, strict=False),
            "zero_capacity_occupancy": result.strict.zero_capacity_occupancy,
            "step": step, "accepted": int(accepted),
            "single_axis_fraction": float(np.mean(single_axis[db.movable])),
        })
    return DensityModeResult(mode, current, records, time.perf_counter() - started)


def run_density_ablation(
    output_root: str | Path, *, db: PlacementDB, initial: np.ndarray, config: AdaptiveConfig,
    c_rows: float, steps: int, max_seconds_per_mode: float,
    gradient_objective: str = "linear",
    modes: tuple[str, ...] = ("full_2d", "dominant_axis", "hybrid_axis"),
) -> Path:
    root = Path(output_root) / Path(db.source_aux).stem / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    (root / "solutions").mkdir(parents=True)
    (root / "figures").mkdir()
    write_placement(root / "solutions" / "solution_input.pl", db, initial)
    results = [
        solve_density_mode(
            db, initial, config, mode=mode, c_rows=c_rows,
            steps=steps, max_seconds=max_seconds_per_mode, gradient_objective=gradient_objective,
        ) for mode in modes
    ]
    rows = [row for result in results for row in result.records]
    with (root / "iterations.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    evaluator = AdaptiveEvaluator(db, config.evaluation)
    initial_eval = evaluator.evaluate(initial)
    summary: dict[str, dict[str, float | int]] = {}
    for result in results:
        write_placement(root / "solutions" / f"solution_{result.mode}.pl", db, result.centres)
        final = evaluator.evaluate(result.centres)
        summary[result.mode] = {
            "hpwl": final.hpwl, "strict_overflow_percent": evaluator.overflow_percent(final, strict=True),
            "legacy_overflow_percent": evaluator.overflow_percent(final, strict=False),
            "strict_density_linear": final.strict.linear,
            "zero_capacity_occupancy": final.strict.zero_capacity_occupancy,
            "accepted_steps": sum(int(row["accepted"]) for row in result.records),
            "steps": len(result.records), "elapsed_seconds": result.elapsed_seconds,
        }
    meta = {
        "method": "exact-oracle density coordinate-direction ablation",
        "source_aux": str(Path(db.source_aux).resolve()), "created_utc": datetime.now(timezone.utc).isoformat(),
        "c_rows": c_rows, "steps": steps, "max_seconds_per_mode": max_seconds_per_mode,
        "gradient_objective": gradient_objective,
        "axis_dominance_ratio": config.direction.axis_dominance_ratio,
        "input": {
            "hpwl": initial_eval.hpwl,
            "strict_overflow_percent": evaluator.overflow_percent(initial_eval, strict=True),
            "strict_density_linear": initial_eval.strict.linear,
        },
        "results": summary,
    }
    (root / "metadata.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
    render_density_ablation(root)
    return root


def render_density_ablation(root: str | Path) -> list[Path]:
    root = Path(root)
    with (root / "iterations.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    colors = {"full_2d": "#0072B2", "dominant_axis": "#D55E00", "hybrid_axis": "#009E73"}
    fig, axes = plt.subplots(2, 2, figsize=(10.2, 6.2), constrained_layout=True)
    interactive = go.Figure()
    for mode, color in colors.items():
        local = [row for row in rows if row["mode"] == mode]
        if not local:
            continue
        value = lambda key: np.asarray([float(row[key]) for row in local])
        x, overflow = value("iteration"), value("strict_overflow_percent")
        axes[0, 0].plot(x, overflow, color=color, label=mode)
        axes[0, 1].plot(x, value("density_descent_objective"), color=color, label=mode)
        axes[1, 0].plot(x, (value("hpwl") - value("hpwl")[0]) / value("hpwl")[0] * 100.0, color=color, label=mode)
        window = min(30, len(local))
        rate = np.convolve(value("accepted"), np.ones(window) / window, mode="same") * 100.0
        axes[1, 1].plot(x, rate, color=color, label=mode)
        interactive.add_trace(go.Scatter(x=x, y=overflow, name=mode))
    axes[0, 0].set(xlabel="Iteration", ylabel="Strict overflow (%)", title="Overflow reduction")
    axes[0, 1].set(xlabel="Iteration", ylabel="Descent objective", title="Exact density objective")
    axes[1, 0].set(xlabel="Iteration", ylabel="HPWL change (%)", title="Wirelength side effect")
    axes[1, 1].set(xlabel="Iteration", ylabel="Rolling acceptance (%)", title="Exact-oracle acceptance")
    for axis in axes.ravel(): axis.grid(alpha=.18); axis.legend()
    written: list[Path] = []
    for suffix in ("png", "pdf"):
        path = root / "figures" / f"01_density_direction_ablation.{suffix}"
        fig.savefig(path, dpi=220 if suffix == "png" else None, bbox_inches="tight"); written.append(path)
    plt.close(fig)
    interactive.update_layout(title="Density direction ablation", xaxis_title="Iteration", yaxis_title="Strict overflow (%)")
    meta = json.loads((root / "metadata.json").read_text(encoding="utf-8"))
    chart = pio.to_html(interactive, include_plotlyjs="inline", full_html=False)
    (root / "report.html").write_text(
        "<!doctype html><meta charset='utf-8'><title>Density direction ablation</title>"
        "<style>body{font-family:Arial;margin:24px;max-width:1400px}pre{background:#f5f5f5;padding:12px}</style>"
        f"<h1>Full-2D vs short-axis density descent</h1>{chart}<pre>{json.dumps(meta, ensure_ascii=False, indent=2)}</pre>",
        encoding="utf-8",
    )
    return written
