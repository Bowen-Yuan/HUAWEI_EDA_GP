"""Composition of exact HPWL and exact linear-overflow density oracles."""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from gpplacer.model.types import PlacementDB
from gpplacer.oracle.density import density_value_gradient
from gpplacer.oracle.grid import DensityGrid
from gpplacer.oracle.hpwl import hpwl_value_gradient


@dataclass(slots=True)
class OracleResult:
    """Objective decomposition and all state required by phase decisions."""

    hpwl: float
    density_linear: float
    density_surrogate_squared: float
    objective: float
    gradient: np.ndarray
    net_extrema: np.ndarray
    active_bins: np.ndarray
    occupancy: np.ndarray

    @property
    def active_bin_count(self) -> int:
        return int(np.count_nonzero(self.active_bins))


class ObjectiveOracle:
    """Evaluate the un-smoothed placement objective and a Clarke subgradient."""

    def __init__(self, db: PlacementDB, grid: DensityGrid, density_lambda: float) -> None:
        if density_lambda < 0:
            raise ValueError("density_lambda must be non-negative")
        self.db = db
        self.grid = grid
        self.density_lambda = float(density_lambda)

    def evaluate(self, centres: np.ndarray) -> OracleResult:
        """Evaluate exact HPWL plus linear overflow at one placement state."""
        hpwl, hpwl_gradient, extrema = hpwl_value_gradient(
            centres, self.db.pin_node, self.db.pin_offset,
            self.db.net_start, self.db.net_weight,
        )
        density, surrogate_squared, density_gradient, occupancy, active = density_value_gradient(
            centres, self.db.width, self.db.height, self.db.movable,
            self.grid.x_edges, self.grid.y_edges, self.grid.available_area,
            self.grid.rho_target,
        )
        gradient = hpwl_gradient + self.density_lambda * density_gradient
        gradient[self.db.fixed] = 0.0
        return OracleResult(
            hpwl=hpwl, density_linear=density, density_surrogate_squared=surrogate_squared,
            objective=hpwl + self.density_lambda * density,
            gradient=gradient, net_extrema=extrema, active_bins=active,
            occupancy=occupancy,
        )
