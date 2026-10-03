"""Large-step, density-only repair used before anchored HPWL compression.

The routine deliberately ignores HPWL during this first phase.  It starts
from a low-HPWL/high-overflow placement and accepts a candidate solely when
the exact linear overflow decreases.  The output can then be used as the
safe anchor for :mod:`gpplacer.experimental.anchored_search`.
"""

from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Callable

import numpy as np

from gpplacer.model.types import PlacementDB
from gpplacer.oracle.density import density_value_gradient
from gpplacer.oracle.grid import DensityGrid
from gpplacer.oracle.objective import ObjectiveOracle, OracleResult
from gpplacer.solver.evacuation import evacuate_overfull_bins
from gpplacer.solver.steps import project_centres


@dataclass(frozen=True, slots=True)
class DensityRepairConfig:
    """Controls for a deliberately aggressive density-only phase."""

    max_iterations: int = 400
    target_overflow_percent: float = 6.0
    initial_step: float = 8.0
    max_displacement_rows: float = 64.0
    min_step: float = 1.0 / 128.0
    max_backtracks: int = 10
    evacuation_rounds: int = 20
    evacuation_max_moves_per_round: int = 10_000
    evacuation_max_radius: int = 112


@dataclass(slots=True)
class DensityRepairResult:
    """Best density state encountered by the repair trajectory."""

    centres: np.ndarray
    oracle: OracleResult
    iterations: int
    reached_target: bool


class DensityRepair:
    """Exact-overflow descent with large but projected cell displacements."""

    def __init__(self, db: PlacementDB, grid: DensityGrid, config: DensityRepairConfig) -> None:
        if not 0.0 < config.target_overflow_percent <= 100.0:
            raise ValueError("target_overflow_percent must lie in (0, 100].")
        self.db = db
        self.grid = grid
        self.config = config
        self.oracle = ObjectiveOracle(db, grid, density_lambda=0.0)
        self.total_available_area = float(grid.available_area.sum())
        self.target_density = self.total_available_area * config.target_overflow_percent / 100.0

    def run(
        self,
        initial_centres: np.ndarray,
        observer: Callable[[dict[str, object], np.ndarray, OracleResult], None] | None = None,
    ) -> DensityRepairResult:
        centres = initial_centres.copy()
        result = self.oracle.evaluate(centres)
        best_centres, best_result = centres.copy(), result
        step = self.config.initial_step
        used_evacuation = False
        for iteration in range(1, self.config.max_iterations + 1):
            if best_result.density_linear <= self.target_density + 1e-9:
                return DensityRepairResult(best_centres, best_result, iteration - 1, True)
            candidate, candidate_result, accepted, backtracks = self._step(centres, result, step)
            if accepted:
                centres, result = candidate, candidate_result
                step = min(step * 1.2, self.config.initial_step * 16.0)
                event = "density_repair_accept"
                if result.density_linear < best_result.density_linear:
                    best_centres, best_result = centres.copy(), result
            else:
                # A zero density subgradient is common when whole cells sit
                # inside active bins.  One deterministic capacity evacuation
                # crosses those bin boundaries; it is accepted only if exact D
                # improves, and never changes the default solver.
                if not used_evacuation and backtracks == 0:
                    used_evacuation = True
                    evacuated, _ = evacuate_overfull_bins(
                        self.db, self.grid, centres,
                        rounds=self.config.evacuation_rounds,
                        max_moves_per_round=self.config.evacuation_max_moves_per_round,
                        max_radius=self.config.evacuation_max_radius,
                    )
                    evacuated_result = self.oracle.evaluate(evacuated)
                    if evacuated_result.density_linear < result.density_linear - 1e-6:
                        centres, result = evacuated, evacuated_result
                        event = "density_evacuation_accept"
                        if result.density_linear < best_result.density_linear:
                            best_centres, best_result = centres.copy(), result
                    else:
                        step = max(step * 0.5, self.config.min_step)
                        event = "density_evacuation_reject"
                else:
                    step = max(step * 0.5, self.config.min_step)
                    event = "density_repair_reject"
            if observer is not None:
                observer(self._record(iteration, event, result, step, backtracks), centres, result)
        return DensityRepairResult(
            best_centres, best_result, self.config.max_iterations,
            best_result.density_linear <= self.target_density + 1e-9,
        )

    def _step(
        self, centres: np.ndarray, result: OracleResult, step: float,
    ) -> tuple[np.ndarray, OracleResult, bool, int]:
        """Backtrack a pure density-subgradient move using exact D acceptance."""
        _, _, density_gradient, _, _ = density_value_gradient(
            centres, self.db.width, self.db.height, self.db.movable,
            self.grid.x_edges, self.grid.y_edges, self.grid.available_area,
            self.grid.rho_target,
        )
        if not np.any(density_gradient[self.db.movable]):
            return centres, result, False, 0
        trial_step = step
        for backtracks in range(self.config.max_backtracks + 1):
            displacement = -trial_step * density_gradient
            displacement[self.db.fixed] = 0.0
            maximum = float(np.max(np.abs(displacement)))
            limit = self.config.max_displacement_rows * self.db.row_height
            if maximum > limit:
                displacement *= limit / maximum
            candidate = project_centres(self.db, centres + displacement)
            candidate_result = self.oracle.evaluate(candidate)
            if candidate_result.density_linear < result.density_linear - 1e-6:
                return candidate, candidate_result, True, backtracks
            trial_step *= 0.5
        return centres, result, False, self.config.max_backtracks

    def _record(
        self, iteration: int, event: str, result: OracleResult, step: float, backtracks: int,
    ) -> dict[str, object]:
        return {
            "iteration": iteration,
            "phase": event,
            "hpwl": result.hpwl,
            "density_linear": result.density_linear,
            "overflow_ratio": result.density_linear / self.total_available_area,
            "overflow_percent": 100.0 * result.density_linear / self.total_available_area,
            "density_surrogate_squared": result.density_surrogate_squared,
            "objective": result.density_linear,
            "step": step,
            "proximal_scale": "",
            "active_nets_change": "",
            "active_bins_change": "",
            "active_bins": result.active_bin_count,
            "bundle_size": "",
            "backtracks": backtracks,
            "serious_step": "",
            "predicted_decrease": "",
            "recoveries": "",
            "iteration_seconds": "",
            "oracle_seconds": "",
            "oracle_calls": "",
            "direction_cosine": "",
        }
