"""Capacity-aware density-bin construction from legal rows and fixed objects."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from gpplacer.model.types import PlacementDB


@dataclass(frozen=True, slots=True)
class DensityGrid:
    """Regular bins with exact usable area derived from irregular site rows."""

    x_edges: np.ndarray
    y_edges: np.ndarray
    available_area: np.ndarray
    rho_target: float

    @property
    def nx(self) -> int:
        return len(self.x_edges) - 1

    @property
    def ny(self) -> int:
        return len(self.y_edges) - 1

    @property
    def bin_width(self) -> float:
        return float(self.x_edges[1] - self.x_edges[0])

    @property
    def bin_height(self) -> float:
        return float(self.y_edges[1] - self.y_edges[0])

    @property
    def bounds(self) -> tuple[float, float, float, float]:
        return (float(self.x_edges[0]), float(self.y_edges[0]),
                float(self.x_edges[-1]), float(self.y_edges[-1]))


def _overlap_1d(low_a: float, high_a: float, low_b: float, high_b: float) -> float:
    return max(0.0, min(high_a, high_b) - max(low_a, low_b))


def _add_rectangle(area: np.ndarray, x_edges: np.ndarray, y_edges: np.ndarray,
                   left: float, bottom: float, right: float, top: float, sign: float) -> None:
    """Accumulate exact rectangle overlap into a regular bin array."""
    nx, ny = len(x_edges) - 1, len(y_edges) - 1
    dx, dy = x_edges[1] - x_edges[0], y_edges[1] - y_edges[0]
    ix0 = max(0, int(math.floor((left - x_edges[0]) / dx)))
    ix1 = min(nx - 1, int(math.floor(np.nextafter(right, -np.inf) - x_edges[0]) / dx))
    iy0 = max(0, int(math.floor((bottom - y_edges[0]) / dy)))
    iy1 = min(ny - 1, int(math.floor(np.nextafter(top, -np.inf) - y_edges[0]) / dy))
    if ix0 > ix1 or iy0 > iy1:
        return
    for iy in range(iy0, iy1 + 1):
        height = _overlap_1d(bottom, top, y_edges[iy], y_edges[iy + 1])
        if height == 0:
            continue
        for ix in range(ix0, ix1 + 1):
            width = _overlap_1d(left, right, x_edges[ix], x_edges[ix + 1])
            area[iy, ix] += sign * width * height


def build_density_grid(db: PlacementDB, *, rho_target: float, bin_rows: int) -> DensityGrid:
    """Build bins whose capacity respects legal rows and fixed-object blockage.

    ``bin_rows`` controls a square bin height measured in site rows.  Bins remain
    regular even when row segments are staggered; their capacity then naturally
    records partial and unavailable regions.
    """
    if not 0.0 < rho_target <= 1.0:
        raise ValueError("rho_target must lie in (0, 1].")
    if bin_rows < 1:
        raise ValueError("bin_rows must be positive.")
    left, bottom, right, top = db.core_bounds
    edge = db.row_height * bin_rows
    nx = math.ceil((right - left) / edge)
    ny = math.ceil((top - bottom) / edge)
    x_edges = left + np.arange(nx + 1, dtype=np.float64) * edge
    y_edges = bottom + np.arange(ny + 1, dtype=np.float64) * edge
    x_edges[-1], y_edges[-1] = right, top
    available = np.zeros((ny, nx), dtype=np.float64)
    for row in db.rows:
        _add_rectangle(available, x_edges, y_edges, row.x, row.y,
                       row.x + row.width, row.y + row.height, 1.0)
    for node in np.flatnonzero(db.fixed):
        centre = db.initial_centres[node]
        _add_rectangle(
            available, x_edges, y_edges,
            centre[0] - db.width[node] / 2.0, centre[1] - db.height[node] / 2.0,
            centre[0] + db.width[node] / 2.0, centre[1] + db.height[node] / 2.0,
            -1.0,
        )
    # Fixed rectangles can overlap outside the usable row region.  Numerical
    # subtraction must never produce a negative physical capacity.
    np.maximum(available, 0.0, out=available)
    return DensityGrid(x_edges, y_edges, available, rho_target)
