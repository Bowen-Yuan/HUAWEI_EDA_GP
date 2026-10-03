"""Fixed, explicit density-evaluation protocols for placement experiments."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

import numpy as np

from gpplacer.model.types import PlacementDB
from gpplacer.oracle.density import density_value_gradient, strict_density_value_gradient
from gpplacer.oracle.grid import DensityGrid
from gpplacer.oracle.hpwl import hpwl_value_gradient


class DensityProtocol(str, Enum):
    """Historical optimization proxies, not ISPD contest evaluation protocols."""

    LEGACY = "legacy"
    STRICT = "strict"


@dataclass(slots=True)
class DensityEvaluation:
    """One density metric plus all state needed for exact search decisions."""

    linear: float
    squared: float
    gradient: np.ndarray
    occupancy: np.ndarray
    active_bins: np.ndarray
    zero_capacity_occupancy: float

    @property
    def active_bin_count(self) -> int:
        return int(np.count_nonzero(self.active_bins))


@dataclass(slots=True)
class DualOracleResult:
    """Exact HPWL with legacy and strict density views of one placement."""

    hpwl: float
    hpwl_gradient: np.ndarray
    net_extrema: np.ndarray
    legacy: DensityEvaluation
    strict: DensityEvaluation


class DualObjectiveOracle:
    """Evaluate both legacy optimization proxies without changing old callers.

    Formal reports must call :func:`gpplacer.official.run_official_evaluation`
    after legalization instead of exposing either proxy as an overflow result.
    """

    def __init__(self, db: PlacementDB, grid: DensityGrid) -> None:
        self.db = db
        self.grid = grid

    def evaluate(self, centres: np.ndarray) -> DualOracleResult:
        hpwl, hpwl_gradient, extrema = hpwl_value_gradient(
            centres, self.db.pin_node, self.db.pin_offset,
            self.db.net_start, self.db.net_weight,
        )
        legacy_value, legacy_squared, legacy_gradient, occupancy, legacy_active = density_value_gradient(
            centres, self.db.width, self.db.height, self.db.movable,
            self.grid.x_edges, self.grid.y_edges, self.grid.available_area,
            self.grid.rho_target,
        )
        strict_value, strict_squared, strict_gradient, strict_occupancy, strict_active, zero_occupancy = (
            strict_density_value_gradient(
                centres, self.db.width, self.db.height, self.db.movable,
                self.grid.x_edges, self.grid.y_edges, self.grid.available_area,
                self.grid.rho_target,
            )
        )
        hpwl_gradient[self.db.fixed] = 0.0
        legacy_gradient[self.db.fixed] = 0.0
        strict_gradient[self.db.fixed] = 0.0
        return DualOracleResult(
            hpwl=hpwl,
            hpwl_gradient=hpwl_gradient,
            net_extrema=extrema,
            legacy=DensityEvaluation(
                legacy_value, legacy_squared, legacy_gradient, occupancy,
                legacy_active, 0.0,
            ),
            strict=DensityEvaluation(
                strict_value, strict_squared, strict_gradient, strict_occupancy,
                strict_active, zero_occupancy,
            ),
        )
