"""Adaptive-controller view of the shared dual density oracle."""

from __future__ import annotations

import numpy as np

from gpplacer.model.types import PlacementDB
from gpplacer.oracle.evaluation import DualObjectiveOracle, DualOracleResult
from gpplacer.oracle.grid import DensityGrid, build_density_grid

from adaptive_pareto.config import EvaluationConfig


class AdaptiveEvaluator:
    """Own the fixed protocol and its normalised metrics for one benchmark."""

    def __init__(self, db: PlacementDB, config: EvaluationConfig) -> None:
        self.db = db
        self.config = config
        self.grid: DensityGrid = build_density_grid(
            db, rho_target=config.rho_target, bin_rows=config.bin_rows,
        )
        self.oracle = DualObjectiveOracle(db, self.grid)
        self.available_area = float(self.grid.available_area.sum())
        left, bottom, right, top = db.core_bounds
        self.hpwl_scale = float(db.net_weight.sum()) * ((right - left) + (top - bottom))

    def evaluate(self, centres: np.ndarray) -> DualOracleResult:
        return self.oracle.evaluate(centres)

    def overflow_percent(self, result: DualOracleResult, *, strict: bool = True) -> float:
        density = result.strict if strict else result.legacy
        return 100.0 * density.linear / max(self.available_area, 1.0e-12)

    def normalized_hpwl(self, result: DualOracleResult) -> float:
        return result.hpwl / max(self.hpwl_scale, 1.0e-12)
