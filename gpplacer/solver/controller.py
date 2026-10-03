"""Controller implementing the report's Explore -> Bundle -> Refine flow."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
import time

import numpy as np
from numba import set_num_threads
import psutil

from gpplacer.model.config import SolverConfig
from gpplacer.model.types import PlacementDB
from gpplacer.multilevel.initialization import cluster_capacity_seed, connectivity_aware_seed
from gpplacer.multilevel.partition_seed import recursive_partition_seed
from gpplacer.multilevel.coarsen import coarsen_once, uncoarsen_centres
from gpplacer.oracle.density import density_value_gradient
from gpplacer.oracle.grid import DensityGrid, build_density_grid
from gpplacer.oracle.hpwl import hpwl_value_gradient
from gpplacer.oracle.objective import ObjectiveOracle, OracleResult
from gpplacer.solver.state import BundleCut, Phase, PlacementState
from gpplacer.solver.steps import active_change, build_preconditioner, limited_descent, project_centres
from gpplacer.solver.evacuation import evacuate_overfull_bins


Observer = Callable[[dict[str, object], PlacementState], None]


@dataclass(slots=True)
class SolveResult:
    """Final state and summary metrics for one completed solver invocation."""

    centres: np.ndarray
    oracle: OracleResult
    elapsed_seconds: float
    iterations: int
    recoveries: int
    density_lambda: float
    grid: DensityGrid


class PlacementSolver:
    """Safe CPU implementation of the planned four-stage nonsmooth algorithm."""

    def __init__(self, db: PlacementDB, config: SolverConfig) -> None:
        self.db = db
        self.config = config
        self.grid = build_density_grid(
            db, rho_target=config.rho_target, bin_rows=config.bin_rows,
        )
        self.preconditioner = build_preconditioner(db)

    def _calibrate_lambda(self, centres: np.ndarray) -> float:
        """Balance HPWL and density subgradient scales once at run start."""
        hpwl, hpwl_gradient, _ = hpwl_value_gradient(
            centres, self.db.pin_node, self.db.pin_offset,
            self.db.net_start, self.db.net_weight,
        )
        _, _, density_gradient, _, _ = density_value_gradient(
            centres, self.db.width, self.db.height, self.db.movable,
            self.grid.x_edges, self.grid.y_edges, self.grid.available_area,
            self.grid.rho_target,
        )
        numerator = float(np.linalg.norm(hpwl_gradient[self.db.movable], ord=1))
        denominator = float(np.linalg.norm(density_gradient[self.db.movable], ord=1))
        if denominator <= 1e-12:
            # A density-free seed is unusual; retain a finite, documented
            # scale instead of introducing a case-specific constant.
            return max(1.0, hpwl / max(float(self.db.movable.sum()), 1.0))
        return float(np.clip(numerator / denominator, 1e-4, 1e4))

    def _multilevel_initial_centres(self) -> np.ndarray:
        """Solve a compact coarse hypergraph before restoring fine coordinates."""
        coarsening = coarsen_once(self.db, self.config.coarsen_degree_limit)
        coarse = coarsening.coarse
        if coarse.node_count >= self.db.node_count or self.config.coarse_iterations == 0:
            return self._finish_initialisation(
                connectivity_aware_seed(self.db, self.config.initialization_rounds),
            )
        coarse_grid = build_density_grid(
            coarse, rho_target=self.config.rho_target, bin_rows=self.config.bin_rows,
        )
        centres = project_centres(
            coarse, connectivity_aware_seed(coarse, self.config.initialization_rounds),
        )
        coarse_lambda = self.config.density_lambda if self.config.density_lambda is not None else 1.0
        oracle = ObjectiveOracle(coarse, coarse_grid, coarse_lambda)
        preconditioner = build_preconditioner(coarse)
        result = oracle.evaluate(centres)
        step = self.config.initial_step
        maximum_displacement = self.config.max_displacement_rows * coarse.row_height
        for _ in range(self.config.coarse_iterations):
            candidate = limited_descent(
                coarse, centres, result.gradient, preconditioner, step, maximum_displacement,
            )
            candidate_result = oracle.evaluate(candidate)
            if self._safe(candidate_result, candidate) and candidate_result.objective < result.objective:
                centres, result = candidate, candidate_result
                step *= 1.1
            else:
                step *= 0.5
        fine_centres = uncoarsen_centres(coarsening, centres, self.config.seed)
        return self._finish_initialisation(fine_centres)

    def _finish_initialisation(self, centres: np.ndarray) -> np.ndarray:
        """Select an opt-in capacity seed while retaining the old default."""
        if self.config.initialization_strategy == "connection":
            return project_centres(self.db, centres)
        if self.config.initialization_strategy == "cluster_capacity":
            return project_centres(self.db, cluster_capacity_seed(
                self.db, self.grid, centres, self.config.cluster_macro_bins,
            ))
        if self.config.initialization_strategy == "recursive_partition":
            return project_centres(self.db, recursive_partition_seed(
                self.db, self.grid, centres,
                leaf_bins=self.config.partition_leaf_bins,
                refine_passes=self.config.partition_refine_passes,
                degree_limit=self.config.coarsen_degree_limit,
            ))
        raise ValueError(f"Unknown initialization strategy: {self.config.initialization_strategy!r}")

    def _safe(self, result: OracleResult, centres: np.ndarray) -> bool:
        numerically_safe = bool(
            np.isfinite(result.objective)
            and np.isfinite(result.gradient).all()
            and np.isfinite(centres).all()
        )
        if not numerically_safe:
            return False
        if self.config.max_memory_mb is None:
            return True
        return psutil.Process().memory_info().rss <= self.config.max_memory_mb * 1024**2

    def _evaluate(self, centres: np.ndarray) -> OracleResult:
        """Time the exact fine-level Oracle for logs and runtime attribution."""
        started = time.perf_counter()
        result = self.oracle.evaluate(centres)
        self._oracle_seconds += time.perf_counter() - started
        self._oracle_calls += 1
        return result

    def _add_cut(self, state: PlacementState, result: OracleResult,
                 value_at_center: float | None = None) -> None:
        value = result.objective if value_at_center is None else value_at_center
        state.bundle.append(
            BundleCut(float(value), result.gradient.astype(np.float32, copy=True), state.iteration)
        )
        while len(state.bundle) > self.config.bundle_size:
            first, second = state.bundle.pop(0), state.bundle.pop(0)
            aggregate = BundleCut(
                value_at_center=(first.value_at_center + second.value_at_center) / 2.0,
                gradient=((first.gradient.astype(np.float64) + second.gradient.astype(np.float64)) / 2.0)
                .astype(np.float32),
                source_iteration=min(first.source_iteration, second.source_iteration),
            )
            state.bundle.insert(0, aggregate)

    def _move_bundle_center(self, state: PlacementState, displacement: np.ndarray) -> None:
        """Re-express every stored affine cut at a newly accepted centre."""
        for cut in state.bundle:
            cut.value_at_center += float(np.sum(cut.gradient * displacement, dtype=np.float64))

    def _bundle_candidate(self, state: PlacementState, maximum_displacement: float) -> tuple[np.ndarray, float]:
        """Approximately solve the small-memory proximal bundle subproblem.

        The primal variable is high-dimensional, but the maximum of at most
        eight affine cuts is handled with a short deterministic active-cut loop.
        Original objective acceptance, rather than this local model, decides
        every serious step.
        """
        displacement = np.zeros_like(state.centres)
        trust = max(state.proximal_scale, 1e-8)
        for inner in range(12):
            scores = np.asarray([
                cut.value_at_center + float(np.sum(cut.gradient * displacement, dtype=np.float64))
                for cut in state.bundle
            ])
            active_cut = state.bundle[int(np.argmax(scores))]
            model_gradient = active_cut.gradient + (self.preconditioner[:, None] / trust) * displacement
            learning_rate = trust / (2.0 + inner)
            displacement -= learning_rate * model_gradient / self.preconditioner[:, None]
            displacement[self.db.fixed] = 0.0
            np.clip(displacement, -maximum_displacement, maximum_displacement, out=displacement)
        candidate = project_centres(self.db, state.centres + displacement)
        actual_displacement = candidate - state.centres
        predicted_value = max(
            cut.value_at_center + float(np.sum(cut.gradient * actual_displacement, dtype=np.float64))
            for cut in state.bundle
        )
        return candidate, predicted_value

    def _explore_step(self, state: PlacementState, maximum_displacement: float) -> tuple[np.ndarray, OracleResult, int]:
        """Try a preconditioned subgradient step with original-objective backtracking."""
        step = state.step
        for backtrack in range(7):
            candidate = limited_descent(
                self.db, state.centres, state.oracle.gradient, self.preconditioner,
                step, maximum_displacement,
            )
            candidate_result = self._evaluate(candidate)
            if self._safe(candidate_result, candidate) and candidate_result.objective < state.oracle.objective:
                state.step = min(step * 1.1, self.config.initial_step * 4.0)
                return candidate, candidate_result, backtrack
            step *= 0.5
        state.step = max(step, 1e-8)
        return state.centres.copy(), state.oracle, 6

    def _refine_step(self, state: PlacementState, maximum_displacement: float) -> tuple[np.ndarray, OracleResult]:
        """Apply a guarded local correction to the highest-sensitivity movable cells."""
        sensitivity = np.linalg.norm(state.oracle.gradient, axis=1)
        movable = np.flatnonzero(self.db.movable)
        count = max(1, int(len(movable) * self.config.refine_fraction))
        chosen = movable[np.argpartition(sensitivity[movable], -count)[-count:]]
        mask = np.zeros(self.db.node_count, dtype=bool)
        mask[chosen] = True
        candidate = limited_descent(
            self.db, state.centres, state.oracle.gradient, self.preconditioner,
            state.step * 0.5, maximum_displacement * 0.5, mask,
        )
        candidate_result = self._evaluate(candidate)
        if self._safe(candidate_result, candidate) and candidate_result.objective < state.oracle.objective:
            return candidate, candidate_result
        return state.centres.copy(), state.oracle

    def solve(self, observer: Observer | None = None) -> SolveResult:
        """Run all enabled stages until an iteration or wall-time budget is exhausted."""
        set_num_threads(self.config.threads)
        centres = self._multilevel_initial_centres()
        if self.config.enable_bin_evacuator:
            centres, self.evacuation_stats = evacuate_overfull_bins(
                self.db, self.grid, centres,
                rounds=self.config.evacuation_rounds,
                max_moves_per_round=self.config.evacuation_max_moves_per_round,
                max_radius=self.config.evacuation_max_radius,
            )
        density_lambda = self.config.density_lambda
        if density_lambda is None:
            density_lambda = self._calibrate_lambda(centres)
        self.oracle = ObjectiveOracle(self.db, self.grid, density_lambda)
        self._oracle_seconds = 0.0
        self._oracle_calls = 0
        total_available_area = float(self.grid.available_area.sum())
        initial = self._evaluate(centres)
        if not self._safe(initial, centres):
            raise RuntimeError("The initial placement produced an unsafe objective state.")
        state = PlacementState(
            centres=centres, oracle=initial, best_centres=centres.copy(), best_oracle=initial,
            step=self.config.initial_step,
            proximal_scale=self.config.proximal_scale * self.db.row_height,
        )
        objective_history = [initial.objective]
        previous_extrema, previous_bins = initial.net_extrema, initial.active_bins
        previous_gradient = initial.gradient.copy()
        started = time.perf_counter()
        if observer is not None:
            observer({
                "iteration": 0, "elapsed_seconds": 0.0, "phase": state.phase.value,
                "hpwl": state.oracle.hpwl, "density_linear": state.oracle.density_linear,
                "overflow_ratio": state.oracle.density_linear / total_available_area,
                "overflow_percent": 100.0 * state.oracle.density_linear / total_available_area,
                "density_surrogate_squared": state.oracle.density_surrogate_squared, "objective": state.oracle.objective,
                "step": state.step, "proximal_scale": state.proximal_scale,
                "active_nets_change": 0.0, "active_bins_change": 0.0,
                "active_bins": state.oracle.active_bin_count, "bundle_size": 0,
                "backtracks": 0, "serious_step": 0, "predicted_decrease": None,
                "recoveries": 0,
            }, state)
        while state.iteration < self.config.max_iterations:
            elapsed = time.perf_counter() - started
            if elapsed >= self.config.max_seconds:
                break
            state.iteration += 1
            iteration_started = time.perf_counter()
            oracle_seconds_before, oracle_calls_before = self._oracle_seconds, self._oracle_calls
            maximum_displacement = max(
                self.config.min_displacement_rows * self.db.row_height,
                self.config.max_displacement_rows * self.db.row_height
                * (1.0 - state.iteration / max(self.config.max_iterations * 1.2, 1.0)),
            )
            backtracks = 0
            predicted = None
            serious = None
            if state.phase is Phase.EXPLORE:
                candidate, candidate_result, backtracks = self._explore_step(state, maximum_displacement)
                serious = candidate_result.objective < state.oracle.objective
            elif state.phase is Phase.STABILIZE:
                if not state.bundle:
                    self._add_cut(state, state.oracle)
                candidate, predicted_value = self._bundle_candidate(state, maximum_displacement)
                candidate_result = self._evaluate(candidate)
                predicted = state.oracle.objective - predicted_value
                actual = state.oracle.objective - candidate_result.objective
                ratio = actual / max(predicted, 1e-12)
                serious = bool(actual > 0.0 and ratio >= self.config.serious_ratio)
                if serious:
                    displacement = candidate - state.centres
                    self._move_bundle_center(state, displacement)
                    state.proximal_scale *= 1.2
                    state.null_steps = 0
                else:
                    value_here = candidate_result.objective + float(
                        np.sum(candidate_result.gradient * (state.centres - candidate), dtype=np.float64)
                    )
                    self._add_cut(state, candidate_result, value_here)
                    state.proximal_scale *= 0.5
                    state.null_steps += 1
                    candidate, candidate_result = state.centres.copy(), state.oracle
                    minimum_proximal_scale = self.config.bundle_min_proximal_scale_rows * self.db.row_height
                    if (state.null_steps >= self.config.bundle_null_limit
                            or state.proximal_scale <= minimum_proximal_scale):
                        state.phase = Phase.EXPLORE
                        state.bundle.clear()
            else:
                candidate, candidate_result = self._refine_step(state, maximum_displacement)
                serious = candidate_result.objective < state.oracle.objective
                if not serious:
                    state.phase = Phase.STABILIZE
                    state.stable_rounds = 0
            if not self._safe(candidate_result, candidate):
                state.recoveries += 1
                state.step *= 0.5
                state.proximal_scale *= 0.5
                candidate, candidate_result = state.checkpoint()
            if serious:
                state.centres, state.oracle = candidate, candidate_result
                if state.oracle.objective < state.best_oracle.objective:
                    state.best_centres, state.best_oracle = state.centres.copy(), state.oracle
            net_change, bin_change = active_change(
                previous_extrema, previous_bins, state.oracle.net_extrema, state.oracle.active_bins,
            )
            previous_extrema, previous_bins = state.oracle.net_extrema, state.oracle.active_bins
            previous_norm = float(np.linalg.norm(previous_gradient))
            current_norm = float(np.linalg.norm(state.oracle.gradient))
            direction_cosine = 1.0 if previous_norm == 0.0 or current_norm == 0.0 else float(
                np.sum(previous_gradient * state.oracle.gradient, dtype=np.float64)
                / (previous_norm * current_norm)
            )
            previous_gradient = state.oracle.gradient.copy()
            objective_history.append(state.oracle.objective)
            if state.phase is Phase.EXPLORE and len(objective_history) > self.config.explore_window:
                old = objective_history[-self.config.explore_window - 1]
                decrease = (old - state.oracle.objective) / max(abs(old), 1.0)
                # Active-set instability is the report-aligned Bundle trigger.
                # Direction reversal/backtracking is useful additional evidence,
                # but making it mandatory can indefinitely suppress Bundle on a
                # smoothly descending yet nonsmooth trajectory.
                directional_instability = backtracks > 0 or direction_cosine < 0.0
                should_stabilize = (
                    decrease <= self.config.switch_rel_improvement
                    and max(net_change, bin_change) >= self.config.switch_active_change
                    and 100.0 * state.oracle.density_linear / total_available_area
                    <= self.config.bundle_overflow_percent_limit
                    and (
                        directional_instability
                        or not self.config.require_directional_instability
                    )
                )
                if should_stabilize:
                    state.phase = Phase.STABILIZE
                    self._add_cut(state, state.oracle)
            if state.phase is Phase.STABILIZE and max(net_change, bin_change) <= self.config.active_stable_change:
                state.stable_rounds += 1
                if self.config.enable_refine and state.stable_rounds >= self.config.active_stable_rounds:
                    state.phase = Phase.REFINE
            else:
                state.stable_rounds = 0
            if observer is not None:
                observer({
                    "iteration": state.iteration, "elapsed_seconds": time.perf_counter() - started,
                    "phase": state.phase.value, "hpwl": state.oracle.hpwl,
                    "density_linear": state.oracle.density_linear,
                    "overflow_ratio": state.oracle.density_linear / total_available_area,
                    "overflow_percent": 100.0 * state.oracle.density_linear / total_available_area,
                    "density_surrogate_squared": state.oracle.density_surrogate_squared,
                    "objective": state.oracle.objective, "step": state.step,
                    "proximal_scale": state.proximal_scale, "active_nets_change": net_change,
                    "active_bins_change": bin_change, "active_bins": state.oracle.active_bin_count,
                    "bundle_size": len(state.bundle), "backtracks": backtracks,
                    "serious_step": int(bool(serious)), "predicted_decrease": predicted,
                    "recoveries": state.recoveries,
                    "direction_cosine": direction_cosine,
                    "iteration_seconds": time.perf_counter() - iteration_started,
                    "oracle_seconds": self._oracle_seconds - oracle_seconds_before,
                    "oracle_calls": self._oracle_calls - oracle_calls_before,
                }, state)
        elapsed = time.perf_counter() - started
        return SolveResult(
            centres=state.best_centres, oracle=state.best_oracle, elapsed_seconds=elapsed,
            iterations=state.iteration, recoveries=state.recoveries,
            density_lambda=density_lambda, grid=self.grid,
        )
