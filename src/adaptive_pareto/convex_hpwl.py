"""Pure convex HPWL baseline with diminishing projected subgradient steps."""

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
from gpplacer.oracle.hpwl import hpwl_value_gradient
from gpplacer.solver.steps import project_centres

from adaptive_pareto.evaluator import AdaptiveEvaluator


@dataclass(slots=True)
class ConvexHPWLResult:
    best_centres: np.ndarray
    best_hpwl: float
    current_centres: np.ndarray
    current_hpwl: float
    records: list[dict[str, float | int]]
    elapsed_seconds: float
    status: str


def solve_convex_hpwl(
    db: PlacementDB, initial: np.ndarray, *, c_rows: float, max_steps: int, max_seconds: float,
) -> ConvexHPWLResult:
    """Run a normalized projected subgradient method and retain best-so-far.

    The scalar normalization preserves the exact HPWL subgradient direction.
    Its maximum per-cell displacement is ``row_height*c_rows/sqrt(k)``.
    No density value, state machine, adaptive penalty, or monotonic acceptance
    is consulted by this diagnostic solver.
    """
    current = project_centres(db, initial)
    current_hpwl, gradient, _ = hpwl_value_gradient(
        current, db.pin_node, db.pin_offset, db.net_start, db.net_weight,
    )
    best = current.copy()
    best_hpwl = current_hpwl
    initial_hpwl = current_hpwl
    started = time.perf_counter()
    records: list[dict[str, float | int]] = []
    status = "step_limit"
    for iteration in range(1, max_steps + 1):
        if time.perf_counter() - started >= max_seconds:
            status = "timeout"
            break
        gradient = gradient.copy()
        gradient[db.fixed] = 0.0
        node_norm = np.linalg.norm(gradient, axis=1)
        maximum = float(node_norm[db.movable].max(initial=0.0))
        l2_norm = float(np.linalg.norm(gradient[db.movable]))
        if maximum <= 1.0e-15:
            status = "stationary"
            break
        max_displacement = db.row_height * c_rows / np.sqrt(iteration)
        alpha = max_displacement / maximum
        current = project_centres(db, current - alpha * gradient)
        current_hpwl, gradient, _ = hpwl_value_gradient(
            current, db.pin_node, db.pin_offset, db.net_start, db.net_weight,
        )
        improved = current_hpwl < best_hpwl
        if improved:
            best_hpwl = current_hpwl
            best = current.copy()
        records.append({
            "iteration": iteration, "elapsed_seconds": time.perf_counter() - started,
            "current_hpwl": current_hpwl, "best_hpwl": best_hpwl,
            "best_improvement_percent": (initial_hpwl - best_hpwl) / max(initial_hpwl, 1.0) * 100.0,
            "max_displacement": max_displacement, "alpha": alpha,
            "gradient_l2": l2_norm, "gradient_max_node": maximum,
            "improved_best": int(improved),
        })
    elapsed = time.perf_counter() - started
    return ConvexHPWLResult(best, best_hpwl, current, current_hpwl, records, elapsed, status)


def write_convex_hpwl_run(
    output_root: str | Path, *, db: PlacementDB, initial: np.ndarray, result: ConvexHPWLResult,
    evaluator: AdaptiveEvaluator, c_rows: float, max_steps: int, max_seconds: float,
) -> Path:
    root = Path(output_root) / Path(db.source_aux).stem / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    (root / "solutions").mkdir(parents=True)
    (root / "snapshots").mkdir()
    (root / "figures").mkdir()
    write_placement(root / "solutions" / "solution_input.pl", db, initial)
    write_placement(root / "solutions" / "solution_best_hpwl.pl", db, result.best_centres)
    np.savez_compressed(root / "snapshots" / "input.npz", centres=initial)
    np.savez_compressed(root / "snapshots" / "best.npz", centres=result.best_centres)
    with (root / "iterations.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result.records[0]) if result.records else ["iteration"])
        writer.writeheader()
        writer.writerows(result.records)
    input_eval = evaluator.evaluate(initial)
    best_eval = evaluator.evaluate(result.best_centres)
    metadata = {
        "method": "normalized projected HPWL subgradient",
        "objective": "HPWL only; density and state machine disabled",
        "source_aux": str(Path(db.source_aux).resolve()), "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": result.status, "c_rows": c_rows, "max_steps": max_steps, "max_seconds": max_seconds,
        "elapsed_seconds": result.elapsed_seconds, "records": len(result.records),
        "initial_hpwl": input_eval.hpwl, "best_hpwl": best_eval.hpwl,
        "initial_strict_overflow_percent": evaluator.overflow_percent(input_eval, strict=True),
        "best_strict_overflow_percent": evaluator.overflow_percent(best_eval, strict=True),
    }
    (root / "metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    render_convex_hpwl_run(root)
    return root


def render_convex_hpwl_run(run_dir: str | Path) -> list[Path]:
    root = Path(run_dir)
    with (root / "iterations.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return []
    value = lambda key: np.asarray([float(row[key]) for row in rows])
    iteration, current, best = value("iteration"), value("current_hpwl") / 1e6, value("best_hpwl") / 1e6
    fig, axes = plt.subplots(2, 2, figsize=(10.0, 6.0), constrained_layout=True)
    axes[0, 0].plot(iteration, current, color="#9ecae1", linewidth=.7, label="Current")
    axes[0, 0].plot(iteration, best, color="#0072B2", linewidth=1.2, label="Best-so-far")
    axes[0, 0].set(xlabel="Iteration", ylabel="HPWL (M)", title="Pure convex HPWL convergence")
    axes[0, 0].legend()
    axes[0, 1].plot(iteration, value("best_improvement_percent"), color="#009E73")
    axes[0, 1].set(xlabel="Iteration", ylabel="Improvement (%)", title="Best improvement from input")
    axes[1, 0].plot(iteration, value("max_displacement"), color="#E69F00")
    axes[1, 0].set(xlabel="Iteration", ylabel="Maximum displacement", title=r"Diminishing $c/\sqrt{k}$ schedule")
    axes[1, 1].plot(iteration, value("gradient_l2"), color="#D55E00")
    axes[1, 1].set_yscale("log")
    axes[1, 1].set(xlabel="Iteration", ylabel="Gradient L2", title="HPWL subgradient norm")
    for axis in axes.ravel():
        axis.grid(alpha=.18)
    written: list[Path] = []
    for suffix in ("png", "pdf"):
        path = root / "figures" / f"01_convex_hpwl_convergence.{suffix}"
        fig.savefig(path, dpi=220 if suffix == "png" else None, bbox_inches="tight")
        written.append(path)
    plt.close(fig)
    interactive = go.Figure()
    interactive.add_trace(go.Scatter(x=iteration, y=current, name="Current HPWL (M)", line={"color": "#9ecae1"}))
    interactive.add_trace(go.Scatter(x=iteration, y=best, name="Best-so-far HPWL (M)", line={"color": "#0072B2"}))
    interactive.update_layout(title="Pure convex HPWL convergence", xaxis_title="Iteration", yaxis_title="HPWL (M)")
    meta = json.loads((root / "metadata.json").read_text(encoding="utf-8"))
    chart = pio.to_html(interactive, include_plotlyjs="inline", full_html=False)
    (root / "report.html").write_text(
        "<!doctype html><meta charset='utf-8'><title>Convex HPWL baseline</title>"
        "<style>body{font-family:Arial;margin:24px;max-width:1400px}pre{background:#f5f5f5;padding:12px}</style>"
        f"<h1>Convex HPWL baseline</h1>{chart}<h2>Configuration and result</h2>"
        f"<pre>{json.dumps(meta, ensure_ascii=False, indent=2)}</pre>", encoding="utf-8",
    )
    return written
