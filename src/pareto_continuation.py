"""Adaptive-penalty Pareto experiment, independent of hard-constraint solvers.

The controller first samples the HPWL/density trade-off through short pure
phases, then follows an adaptive scalarisation ``HPWL + lambda * D``.  Overflow
is never used as a hard acceptance constraint; every reported point is instead
kept or discarded by its current scalar merit and a bounded nondominated archive.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
import time

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np

from gpplacer.io.bookshelf import load_bookshelf
from gpplacer.io.placement import read_placement_centres, write_placement
from gpplacer.oracle.density import density_value_gradient
from gpplacer.oracle.grid import build_density_grid
from gpplacer.oracle.objective import ObjectiveOracle, OracleResult
from gpplacer.solver.steps import build_preconditioner, project_centres


@dataclass(slots=True)
class ParetoConfig:
    rho_target: float = .8
    bin_rows: int = 8
    warm_cycles: int = 3
    warm_hpwl_steps: int = 110
    warm_density_steps: int = 35
    joint_steps: int = 520
    target_overflow_percent: float = 25.0
    overflow_band_percent: float = 5.0
    lambda_gain: float = .70
    initial_step: float = 4.0
    min_step: float = 1.0 / 256.0
    trust_rows: float = 4.0
    max_trust_rows: float = 16.0
    archive_size: int = 12
    max_seconds: float = 900.0
    output_root: str = "adaptive_pareto_soft/results"


@dataclass(slots=True)
class ParetoPoint:
    centres: np.ndarray
    hpwl: float
    density: float
    overflow_percent: float
    phase: str
    iteration: int


class BoundedParetoArchive:
    """Keep a small set of exact nondominated placements without huge memory use."""

    def __init__(self, maximum: int) -> None:
        self.maximum = maximum
        self.points: list[ParetoPoint] = []

    def consider(self, centres: np.ndarray, result: OracleResult, overflow: float,
                 phase: str, iteration: int) -> bool:
        # H and D are both minimised. Equal points are ignored deterministically.
        if any(point.hpwl <= result.hpwl and point.density <= result.density_linear for point in self.points):
            return False
        self.points = [point for point in self.points
                       if not (result.hpwl <= point.hpwl and result.density_linear <= point.density)]
        self.points.append(ParetoPoint(centres.copy(), result.hpwl, result.density_linear, overflow, phase, iteration))
        if len(self.points) > self.maximum:
            ordered = sorted(self.points, key=lambda point: point.density)
            pick = np.linspace(0, len(ordered) - 1, self.maximum).round().astype(int)
            self.points = [ordered[index] for index in np.unique(pick)]
        return True

    def target_representative(self, target: float, band: float) -> ParetoPoint:
        """Return the shortest archived point in the requested soft band.

        If the archive has no point in the band, choose the closest one. This
        is a reporting choice only; it never alters acceptance or adds a hard
        constraint to the optimisation.
        """
        def key(point: ParetoPoint) -> tuple[float, float]:
            distance = max(abs(point.overflow_percent - target) - band, 0.0)
            return distance, point.hpwl
        return min(self.points, key=key)


class AdaptiveParetoSolver:
    def __init__(self, db, grid, cfg: ParetoConfig) -> None:
        self.db, self.grid, self.cfg = db, grid, cfg
        self.oracle = ObjectiveOracle(db, grid, density_lambda=0.0)
        self.preconditioner = build_preconditioner(db)
        self.available = float(grid.available_area.sum())
        self.node_degree = np.diff(db.node_net_start).astype(np.float64)
        self.archive = BoundedParetoArchive(cfg.archive_size)
        self.records: list[dict[str, object]] = []
        self.iteration = 0

    def overflow(self, result: OracleResult) -> float:
        return 100.0 * result.density_linear / max(self.available, 1e-12)

    def density_gradient(self, centres: np.ndarray) -> np.ndarray:
        _, _, gradient, _, _ = density_value_gradient(
            centres, self.db.width, self.db.height, self.db.movable,
            self.grid.x_edges, self.grid.y_edges, self.grid.available_area, self.grid.rho_target)
        gradient[self.db.fixed] = 0.0
        return gradient

    def candidate(self, centres: np.ndarray, direction: np.ndarray, step: float, trust: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        displacement = step * direction
        np.clip(displacement, -trust[:, None], trust[:, None], out=displacement)
        displacement[self.db.fixed] = 0.0
        return project_centres(self.db, centres + displacement), displacement

    def record(self, phase: str, event: str, result: OracleResult, step: float, trust: np.ndarray,
               lambda_value: float, accepted: bool, centres: np.ndarray) -> None:
        self.iteration += 1
        overflow = self.overflow(result)
        archived = self.archive.consider(centres, result, overflow, phase, self.iteration) if accepted else False
        movable_trust = trust[self.db.movable]
        self.records.append({"iteration": self.iteration, "phase": phase, "event": event,
                             "hpwl": result.hpwl, "density_linear": result.density_linear,
                             "overflow_percent": overflow, "lambda": lambda_value, "step": step,
                             "trust_p50": float(np.quantile(movable_trust, .5)),
                             "trust_p90": float(np.quantile(movable_trust, .9)),
                             "accepted": int(accepted), "archived": int(archived)})

    def centre_bins(self, centres: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        ix = np.searchsorted(self.grid.x_edges, centres[:, 0], side="right") - 1
        iy = np.searchsorted(self.grid.y_edges, centres[:, 1], side="right") - 1
        return np.clip(ix, 0, self.grid.nx - 1), np.clip(iy, 0, self.grid.ny - 1)

    def density_evacuation(self, centres: np.ndarray, result: OracleResult, limit: int = 6_000) -> np.ndarray | None:
        """Density-only plateau escape. Exact D reduction remains mandatory."""
        capacity = self.grid.rho_target * self.grid.available_area
        excess = np.maximum(result.occupancy - capacity, 0.0)
        slack = capacity - result.occupancy
        sources = np.flatnonzero(excess.ravel() > 0.0)
        targets = np.flatnonzero(slack.ravel() > 0.0)
        if not len(sources) or not len(targets):
            return None
        ix, iy = self.centre_bins(centres); flat = iy * self.grid.nx + ix
        area = self.db.width * self.db.height
        candidate, approximate_slack = centres.copy(), slack.copy()
        targets = targets[np.argsort(approximate_slack.ravel()[targets])[::-1]]
        moved = 0
        for source in sources[np.argsort(excess.ravel()[sources])[::-1]]:
            if moved >= limit:
                break
            sy, sx = divmod(int(source), self.grid.nx)
            source_credit = 0.0
            members = np.flatnonzero(self.db.movable & (flat == source))
            members = members[np.argsort(self.node_degree[members] * np.maximum(area[members], 1.0), kind="stable")]
            for node in members:
                if moved >= limit or source_credit >= excess[sy, sx]:
                    break
                feasible = targets[approximate_slack.ravel()[targets] >= area[node]]
                if not len(feasible):
                    break
                target = int(feasible[0]); ty, tx = divmod(target, self.grid.nx)
                candidate[node] = (.5 * (self.grid.x_edges[tx] + self.grid.x_edges[tx + 1]),
                                   .5 * (self.grid.y_edges[ty] + self.grid.y_edges[ty + 1]))
                approximate_slack[ty, tx] -= area[node]
                source_credit += area[node]
                moved += 1
        return project_centres(self.db, candidate) if moved else None

    def pure_phase(self, centres: np.ndarray, result: OracleResult, trust: np.ndarray, *, phase: str,
                   mode: str, steps: int, started: float) -> tuple[np.ndarray, OracleResult]:
        step = self.cfg.initial_step
        for _ in range(steps):
            if time.perf_counter() - started >= self.cfg.max_seconds:
                break
            gradient = result.gradient if mode == "hpwl" else self.density_gradient(centres)
            direction = -gradient / self.preconditioner[:, None]
            candidate, displacement = self.candidate(centres, direction, step, trust)
            trial = self.oracle.evaluate(candidate)
            accepted = trial.hpwl < result.hpwl - 1e-8 if mode == "hpwl" else trial.density_linear < result.density_linear - 1e-8
            event = f"{mode}_accept" if accepted else f"{mode}_reject"
            if not accepted and mode == "density":
                evacuated = self.density_evacuation(centres, result)
                if evacuated is not None:
                    relocated = self.oracle.evaluate(evacuated)
                    if relocated.density_linear < result.density_linear - 1e-8:
                        candidate, trial, accepted, event = evacuated, relocated, True, "density_evacuation"
                        displacement = candidate - centres
            moved = np.any(np.abs(displacement) > 1e-12, axis=1) & self.db.movable
            if accepted:
                centres, result = candidate, trial
                trust[moved] = np.minimum(trust[moved] * 1.07, self.cfg.max_trust_rows * self.db.row_height)
                step = min(step * 1.08, self.cfg.initial_step * 16.0)
            else:
                trust[moved] = np.maximum(trust[moved] * .65, .25 * self.db.row_height)
                step *= .5
                if step < self.cfg.min_step:
                    step = self.cfg.initial_step
                    event = f"{mode}_restart"
            self.record(phase, event, result, step, trust, 0.0, accepted, centres)
        return centres, result

    def joint_phase(self, centres: np.ndarray, result: OracleResult, trust: np.ndarray, started: float) -> tuple[np.ndarray, OracleResult]:
        step = self.cfg.initial_step
        for _ in range(self.cfg.joint_steps):
            if time.perf_counter() - started >= self.cfg.max_seconds:
                break
            g_hpwl = result.gradient
            g_density = self.density_gradient(centres)
            h_norm = float(np.linalg.norm(g_hpwl[self.db.movable]))
            d_norm = float(np.linalg.norm(g_density[self.db.movable]))
            base_lambda = h_norm / max(d_norm, 1e-12)
            pressure = (self.overflow(result) - self.cfg.target_overflow_percent) / max(self.cfg.overflow_band_percent, 1e-9)
            lambda_value = base_lambda * float(np.exp(self.cfg.lambda_gain * np.clip(pressure, -4.0, 4.0)))
            gradient = g_hpwl + lambda_value * g_density
            direction = -gradient / self.preconditioner[:, None]
            candidate, displacement = self.candidate(centres, direction, step, trust)
            trial = self.oracle.evaluate(candidate)
            current_merit = result.hpwl + lambda_value * result.density_linear
            trial_merit = trial.hpwl + lambda_value * trial.density_linear
            accepted = trial_merit < current_merit - 1e-8
            moved = np.any(np.abs(displacement) > 1e-12, axis=1) & self.db.movable
            if accepted:
                centres, result = candidate, trial
                trust[moved] = np.minimum(trust[moved] * 1.06, self.cfg.max_trust_rows * self.db.row_height)
                step = min(step * 1.06, self.cfg.initial_step * 16.0)
                event = "joint_accept"
            else:
                trust[moved] = np.maximum(trust[moved] * .65, .25 * self.db.row_height)
                step *= .5
                event = "joint_reject"
                if step < self.cfg.min_step:
                    step = self.cfg.initial_step
                    event = "joint_restart"
            self.record("joint", event, result, step, trust, lambda_value, accepted, centres)
        return centres, result

    def run(self, initial: np.ndarray) -> tuple[np.ndarray, OracleResult, ParetoPoint, list[dict[str, object]], float]:
        centres = project_centres(self.db, initial)
        result = self.oracle.evaluate(centres)
        trust = np.full(self.db.node_count, self.cfg.trust_rows * self.db.row_height)
        trust[self.db.fixed] = 0.0
        self.archive.consider(centres, result, self.overflow(result), "input", 0)
        started = time.perf_counter()
        for cycle in range(1, self.cfg.warm_cycles + 1):
            centres, result = self.pure_phase(centres, result, trust, phase=f"warm_hpwl_{cycle}",
                                               mode="hpwl", steps=self.cfg.warm_hpwl_steps, started=started)
            centres, result = self.pure_phase(centres, result, trust, phase=f"warm_density_{cycle}",
                                               mode="density", steps=self.cfg.warm_density_steps, started=started)
        centres, result = self.joint_phase(centres, result, trust, started)
        representative = self.archive.target_representative(
            self.cfg.target_overflow_percent, self.cfg.overflow_band_percent)
        return centres, result, representative, self.records, time.perf_counter() - started


def plot(root: Path, rows: list[dict[str, object]], archive: BoundedParetoArchive, final: OracleResult, final_overflow: float) -> None:
    figures = root / "figures"; figures.mkdir(exist_ok=True)
    x = np.array([float(row["iteration"]) for row in rows]); hpwl = np.array([float(row["hpwl"]) for row in rows])
    overflow = np.array([float(row["overflow_percent"]) for row in rows]); lam = np.array([float(row["lambda"]) for row in rows])
    phase = [str(row["phase"]) for row in rows]
    color = ["#0072B2" if "hpwl" in item else "#D55E00" if "density" in item else "#009E73" for item in phase]
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.4), constrained_layout=True, dpi=180)
    axes[0].plot(x, hpwl, color="#555555", linewidth=.9); axes[0].scatter(x, hpwl, c=color, s=7)
    axes[0].yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value / 1e9:.3f}"))
    axes[0].set(xlabel="record", ylabel="HPWL (×10⁹)", title="Warm-up and adaptive joint descent")
    axes[1].plot(x, overflow, color="#D55E00", linewidth=1.1); axes[1].axhspan(20, 30, color="#009E73", alpha=.10)
    axes[1].set(xlabel="record", ylabel="overflow (%)", title="Overflow trajectory; no hard cap")
    axes[2].plot(x, np.maximum(lam, 1e-12), color="#7A5195", linewidth=1.1)
    axes[2].set_yscale("log"); axes[2].set(xlabel="record", ylabel="adaptive λ (log)", title="Penalty adaptation")
    fig.savefig(figures / "adaptive_convergence.png", dpi=220); plt.close(fig)
    fig, axis = plt.subplots(figsize=(5.1, 3.8), constrained_layout=True, dpi=180)
    axis.plot(overflow, hpwl / 1e9, color="#BBBBBB", linewidth=.8, zorder=1)
    axis.scatter(overflow, hpwl / 1e9, c=color, s=7, alpha=.65, zorder=2)
    points = sorted(archive.points, key=lambda point: point.overflow_percent)
    axis.plot([point.overflow_percent for point in points], [point.hpwl / 1e9 for point in points],
              "o-", color="#000000", linewidth=1.4, markersize=4, label="bounded nondominated archive")
    axis.scatter(final_overflow, final.hpwl / 1e9, marker="*", s=90, color="#E45756", label="final")
    axis.set(xlabel="linear overflow (%)", ylabel="HPWL (×10⁹)", title="Observed HPWL-overflow frontier")
    axis.legend(fontsize=7); fig.savefig(figures / "pareto_frontier.png", dpi=220); plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Adaptive-penalty Pareto placement experiment")
    parser.add_argument("--aux", required=True); parser.add_argument("--placement", required=True)
    parser.add_argument("--warm-cycles", type=int); parser.add_argument("--joint-steps", type=int)
    parser.add_argument("--target-overflow-percent", type=float); parser.add_argument("--output-root")
    args = parser.parse_args(); cfg = ParetoConfig()
    if args.warm_cycles is not None: cfg.warm_cycles = args.warm_cycles
    if args.joint_steps is not None: cfg.joint_steps = args.joint_steps
    if args.target_overflow_percent is not None: cfg.target_overflow_percent = args.target_overflow_percent
    if args.output_root is not None: cfg.output_root = args.output_root
    db = load_bookshelf(args.aux); grid = build_density_grid(db, rho_target=cfg.rho_target, bin_rows=cfg.bin_rows)
    initial = read_placement_centres(args.placement, db)
    solver = AdaptiveParetoSolver(db, grid, cfg)
    initial = project_centres(db, initial); initial_result = solver.oracle.evaluate(initial)
    final_centres, final, representative, rows, elapsed = solver.run(initial)
    root = Path(cfg.output_root) / Path(args.aux).stem / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    root.mkdir(parents=True, exist_ok=False)
    write_placement(root / "solution_final.pl", db, final_centres)
    write_placement(root / "solution_pareto_target.pl", db, representative.centres)
    with (root / "iterations.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    archive_rows = [{"hpwl": point.hpwl, "density_linear": point.density, "overflow_percent": point.overflow_percent,
                     "phase": point.phase, "iteration": point.iteration} for point in solver.archive.points]
    with (root / "pareto_archive.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(archive_rows[0])); writer.writeheader(); writer.writerows(archive_rows)
    metadata = {"status": "completed", "source_aux": str(Path(args.aux).resolve()),
                "input_placement": str(Path(args.placement).resolve()), "config": asdict(cfg),
                "initial_hpwl": initial_result.hpwl, "initial_overflow_percent": solver.overflow(initial_result),
                "final_hpwl": final.hpwl, "final_overflow_percent": solver.overflow(final),
                "target_representative_hpwl": representative.hpwl,
                "target_representative_overflow_percent": representative.overflow_percent,
                "records": len(rows), "archive_points": len(solver.archive.points), "elapsed_seconds": elapsed,
                "overflow_hard_constraint": False}
    (root / "metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    plot(root, rows, solver.archive, final, solver.overflow(final))
    print(json.dumps({"run_dir": str(root), "final_hpwl": final.hpwl, "final_overflow_percent": solver.overflow(final),
                      "target_representative_hpwl": representative.hpwl,
                      "target_representative_overflow_percent": representative.overflow_percent}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
