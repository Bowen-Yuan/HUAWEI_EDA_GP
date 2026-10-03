"""Numerical experiment: convex HPWL rescue followed by lambda continuation."""

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

from adaptive_pareto.archive import StratifiedParetoArchive
from adaptive_pareto.config import AdaptiveConfig
from adaptive_pareto.evaluator import AdaptiveEvaluator


@dataclass(slots=True)
class TwoPhaseResult:
    transition_centres: np.ndarray
    final_centres: np.ndarray
    target_centres: np.ndarray
    records: list[dict[str, float | int | str]]
    archive: StratifiedParetoArchive
    elapsed_seconds: float
    transition_hpwl: float
    status: str


def _max_node_norm(gradient: np.ndarray, movable: np.ndarray) -> float:
    return float(np.linalg.norm(gradient[movable], axis=1).max(initial=0.0))


def solve_two_phase(
    db: PlacementDB, initial: np.ndarray, config: AdaptiveConfig, *, hpwl_target: float,
    hpwl_c_rows: float, hpwl_max_steps: int, joint_steps: int, joint_c_rows: float,
    lambda_max: float, lambda_power: float, max_seconds: float, density_sample_every: int = 50,
    lambda_start: float = 0.0, density_metric: str = "strict",
) -> TwoPhaseResult:
    """Run nonmonotone convex HPWL rescue, then raw-gradient lambda continuation."""
    if density_metric not in {"strict", "legacy"}:
        raise ValueError(f"unsupported density metric: {density_metric}")
    evaluator = AdaptiveEvaluator(db, config.evaluation)
    archive = StratifiedParetoArchive(
        per_band=config.joint.archive_per_band, maximum=config.joint.archive_maximum,
    )
    current = project_centres(db, initial)
    input_result = evaluator.evaluate(current)
    archive.consider(
        current, input_result,
        strict_overflow_percent=evaluator.overflow_percent(input_result, strict=True),
        legacy_overflow_percent=evaluator.overflow_percent(input_result, strict=False),
        phase="input", iteration=0,
    )
    hpwl, gradient, _ = hpwl_value_gradient(current, db.pin_node, db.pin_offset, db.net_start, db.net_weight)
    best_hpwl, best_centres = hpwl, current.copy()
    records: list[dict[str, float | int | str]] = []
    started = time.perf_counter()
    status = "completed"
    for iteration in range(1, hpwl_max_steps + 1):
        if time.perf_counter() - started >= max_seconds:
            status = "timeout_hpwl"
            break
        gradient = gradient.copy()
        gradient[db.fixed] = 0.0
        maximum = _max_node_norm(gradient, db.movable)
        if maximum <= 1.0e-15:
            status = "stationary_hpwl"
            break
        displacement = db.row_height * hpwl_c_rows / np.sqrt(iteration)
        current = project_centres(db, current - displacement / maximum * gradient)
        hpwl, gradient, _ = hpwl_value_gradient(current, db.pin_node, db.pin_offset, db.net_start, db.net_weight)
        if hpwl < best_hpwl:
            best_hpwl, best_centres = hpwl, current.copy()
        strict = legacy = density = float("nan")
        if iteration == 1 or iteration % density_sample_every == 0 or best_hpwl <= hpwl_target:
            sampled = evaluator.evaluate(current)
            strict = evaluator.overflow_percent(sampled, strict=True)
            legacy = evaluator.overflow_percent(sampled, strict=False)
            density = sampled.strict.linear
        records.append({
            "iteration": iteration, "phase": "hpwl_rescue", "hpwl": hpwl,
            "best_hpwl": best_hpwl, "strict_overflow_percent": strict,
            "legacy_overflow_percent": legacy, "strict_density_linear": density,
            "lambda_dimensionless": 0.0, "lambda_effective": 0.0,
            "max_displacement": displacement, "elapsed_seconds": time.perf_counter() - started,
        })
        if best_hpwl <= hpwl_target:
            break

    transition_centres = best_centres.copy()
    current_result = evaluator.evaluate(transition_centres)
    transition_hpwl = current_result.hpwl
    strict = evaluator.overflow_percent(current_result, strict=True)
    legacy = evaluator.overflow_percent(current_result, strict=False)
    archive.consider(
        transition_centres, current_result, strict_overflow_percent=strict,
        legacy_overflow_percent=legacy, phase="transition", iteration=len(records),
    )
    movable = db.movable
    h_norm = float(np.linalg.norm(current_result.hpwl_gradient[movable]))
    density_result = current_result.strict if density_metric == "strict" else current_result.legacy
    d_norm = float(np.linalg.norm(density_result.gradient[movable]))
    density_scale = h_norm / max(d_norm, 1.0e-12)
    for local_step in range(1, joint_steps + 1):
        if time.perf_counter() - started >= max_seconds:
            status = "timeout_joint"
            break
        progress = local_step / max(joint_steps, 1)
        lambda_dimensionless = lambda_start + (lambda_max - lambda_start) * progress ** lambda_power
        lambda_effective = density_scale * lambda_dimensionless
        density_result = current_result.strict if density_metric == "strict" else current_result.legacy
        gradient = current_result.hpwl_gradient + lambda_effective * density_result.gradient
        gradient = gradient.copy()
        gradient[db.fixed] = 0.0
        maximum = _max_node_norm(gradient, movable)
        if maximum <= 1.0e-15:
            status = "stationary_joint"
            break
        displacement = db.row_height * joint_c_rows / np.sqrt(local_step)
        candidate = project_centres(db, transition_centres - displacement / maximum * gradient)
        current_result = evaluator.evaluate(candidate)
        transition_centres = candidate
        strict = evaluator.overflow_percent(current_result, strict=True)
        legacy = evaluator.overflow_percent(current_result, strict=False)
        archive.consider(
            candidate, current_result, strict_overflow_percent=strict,
            legacy_overflow_percent=legacy, phase="joint", iteration=len(records) + 1,
        )
        records.append({
            "iteration": len(records) + 1, "phase": "joint", "hpwl": current_result.hpwl,
            "best_hpwl": min(best_hpwl, current_result.hpwl),
            "strict_overflow_percent": strict, "legacy_overflow_percent": legacy,
            "strict_density_linear": current_result.strict.linear,
            "lambda_dimensionless": lambda_dimensionless, "lambda_effective": lambda_effective,
            "max_displacement": displacement, "elapsed_seconds": time.perf_counter() - started,
        })
    target = archive.target_representative(
        config.evaluation.target_overflow_low, config.evaluation.target_overflow_high,
    )
    return TwoPhaseResult(
        best_centres, transition_centres, target.centres, records, archive,
        time.perf_counter() - started, transition_hpwl, status,
    )


def write_two_phase_run(
    output_root: str | Path, *, db: PlacementDB, initial: np.ndarray, config: AdaptiveConfig,
    result: TwoPhaseResult, parameters: dict[str, float | int | str],
) -> Path:
    root = Path(output_root) / Path(db.source_aux).stem / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    (root / "solutions").mkdir(parents=True)
    (root / "figures").mkdir()
    write_placement(root / "solutions" / "solution_input.pl", db, initial)
    write_placement(root / "solutions" / "solution_transition_hpwl.pl", db, result.transition_centres)
    write_placement(root / "solutions" / "solution_joint_final.pl", db, result.final_centres)
    write_placement(root / "solutions" / "solution_best_target.pl", db, result.target_centres)
    with (root / "iterations.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result.records[0]))
        writer.writeheader(); writer.writerows(result.records)
    evaluator = AdaptiveEvaluator(db, config.evaluation)
    initial_eval, transition_eval = evaluator.evaluate(initial), evaluator.evaluate(result.transition_centres)
    final_eval, target_eval = evaluator.evaluate(result.final_centres), evaluator.evaluate(result.target_centres)
    summary = lambda value: {
        "hpwl": value.hpwl, "strict_overflow_percent": evaluator.overflow_percent(value, strict=True),
        "legacy_overflow_percent": evaluator.overflow_percent(value, strict=False),
        "zero_capacity_occupancy": value.strict.zero_capacity_occupancy,
    }
    meta = {
        "status": result.status, "method": "nonmonotone convex HPWL rescue then raw-gradient lambda continuation",
        "source_aux": str(Path(db.source_aux).resolve()), "created_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": result.elapsed_seconds, "parameters": parameters,
        "density_gradient_scale": next((float(row["lambda_effective"]) / float(row["lambda_dimensionless"])
                                        for row in result.records if float(row["lambda_dimensionless"]) > 0), 0.0),
        "input": summary(initial_eval), "transition": summary(transition_eval),
        "joint_final": summary(final_eval), "best_target": summary(target_eval),
    }
    (root / "metadata.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
    render_two_phase_run(root, config)
    return root


def render_two_phase_run(root: str | Path, config: AdaptiveConfig) -> list[Path]:
    root = Path(root)
    with (root / "iterations.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    value = lambda key: np.asarray([float(row[key]) for row in rows])
    x, hpwl = value("iteration"), value("hpwl") / 1e6
    strict, legacy = value("strict_overflow_percent"), value("legacy_overflow_percent")
    phases = np.asarray([row["phase"] for row in rows])
    joint = phases == "joint"
    transition = int(x[np.flatnonzero(joint)[0]]) if np.any(joint) else int(x[-1])
    fig, axes = plt.subplots(2, 2, figsize=(10.2, 6.2), constrained_layout=True)
    axes[0, 0].plot(x, hpwl, color="#0072B2"); axes[0, 0].axvline(transition, color="#555", linestyle="--")
    axes[0, 0].set(xlabel="Iteration", ylabel="HPWL (M)", title="HPWL rescue and joint convergence")
    sampled = np.isfinite(legacy)
    axes[0, 1].plot(x[sampled], legacy[sampled], color="#D55E00", label="Legacy")
    axes[0, 1].plot(x[sampled], strict[sampled], color="#999999", linestyle="--", label="Strict")
    axes[0, 1].axhspan(config.evaluation.target_overflow_low, config.evaluation.target_overflow_high, color="#009E73", alpha=.13)
    axes[0, 1].axvline(transition, color="#555", linestyle="--")
    axes[0, 1].set(xlabel="Iteration", ylabel="Overflow (%)", title="Overflow trajectory (Legacy primary)")
    axes[0, 1].legend()
    lam = value("lambda_dimensionless")
    axes[1, 0].plot(x, lam, color="#7A5195"); axes[1, 0].set(xlabel="Iteration", ylabel="Dimensionless lambda", title="Lambda continuation")
    axes[1, 1].plot(legacy[joint], hpwl[joint], color="#009E73")
    axes[1, 1].axvspan(config.evaluation.target_overflow_low, config.evaluation.target_overflow_high, color="#009E73", alpha=.13)
    axes[1, 1].set(xlabel="Legacy overflow (%)", ylabel="HPWL (M)", title="Joint Pareto trajectory")
    for axis in axes.ravel(): axis.grid(alpha=.18)
    written: list[Path] = []
    for suffix in ("png", "pdf"):
        path = root / "figures" / f"01_two_phase_convergence.{suffix}"
        fig.savefig(path, dpi=220 if suffix == "png" else None, bbox_inches="tight"); written.append(path)
    plt.close(fig)
    interactive = go.Figure()
    interactive.add_trace(go.Scatter(x=x, y=hpwl, name="HPWL (M)"))
    interactive.add_trace(go.Scatter(x=x[sampled], y=legacy[sampled], name="Legacy overflow (%)", yaxis="y2"))
    interactive.add_trace(go.Scatter(x=x[sampled], y=strict[sampled], name="Strict overflow (%)", yaxis="y2", line={"dash": "dash"}))
    interactive.update_layout(
        title="Two-phase convergence", xaxis_title="Iteration", yaxis_title="HPWL (M)",
        yaxis2={"title": "Overflow (%)", "overlaying": "y", "side": "right"},
    )
    meta = json.loads((root / "metadata.json").read_text(encoding="utf-8"))
    chart = pio.to_html(interactive, include_plotlyjs="inline", full_html=False)
    (root / "report.html").write_text(
        "<!doctype html><meta charset='utf-8'><title>Two-phase experiment</title>"
        "<style>body{font-family:Arial;margin:24px;max-width:1400px}pre{background:#f5f5f5;padding:12px}</style>"
        f"<h1>HPWL rescue → lambda continuation</h1>{chart}<pre>{json.dumps(meta, ensure_ascii=False, indent=2)}</pre>",
        encoding="utf-8",
    )
    return written
