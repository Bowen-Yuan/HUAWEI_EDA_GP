"""Exact rectangular cell-bin overlap density oracle without smoothing."""

from __future__ import annotations

import numpy as np
from numba import njit


@njit(cache=True)
def _cell_range(low: float, high: float, origin: float, step: float, count: int) -> tuple[int, int]:
    first = int(np.floor((low - origin) / step))
    last = int(np.floor((np.nextafter(high, -np.inf) - origin) / step))
    if first < 0:
        first = 0
    if last >= count:
        last = count - 1
    return first, last


@njit(cache=True)
def density_value_gradient(
    centres: np.ndarray, width: np.ndarray, height: np.ndarray, movable: np.ndarray,
    x_edges: np.ndarray, y_edges: np.ndarray, available: np.ndarray, rho_target: float,
) -> tuple[float, float, np.ndarray, np.ndarray, np.ndarray]:
    """Return linear overflow, squared diagnostic, subgradient, occupancy, active bins."""
    ny, nx = available.shape
    occupancy = np.zeros((ny, nx), dtype=np.float64)
    x0, y0 = x_edges[0], y_edges[0]
    dx, dy = x_edges[1] - x_edges[0], y_edges[1] - y_edges[0]
    for node in range(centres.shape[0]):
        if not movable[node]:
            continue
        left = centres[node, 0] - width[node] / 2.0
        right = centres[node, 0] + width[node] / 2.0
        bottom = centres[node, 1] - height[node] / 2.0
        top = centres[node, 1] + height[node] / 2.0
        ix0, ix1 = _cell_range(left, right, x0, dx, nx)
        iy0, iy1 = _cell_range(bottom, top, y0, dy, ny)
        if ix0 > ix1 or iy0 > iy1:
            continue
        for iy in range(iy0, iy1 + 1):
            overlap_y = min(top, y_edges[iy + 1]) - max(bottom, y_edges[iy])
            if overlap_y <= 0.0:
                continue
            for ix in range(ix0, ix1 + 1):
                overlap_x = min(right, x_edges[ix + 1]) - max(left, x_edges[ix])
                if overlap_x > 0.0:
                    occupancy[iy, ix] += overlap_x * overlap_y
    active = np.zeros((ny, nx), dtype=np.bool_)
    linear = 0.0
    squared = 0.0
    for iy in range(ny):
        for ix in range(nx):
            capacity = available[iy, ix]
            if capacity <= 0.0:
                continue
            excess_area = occupancy[iy, ix] - rho_target * capacity
            if excess_area > 0.0:
                active[iy, ix] = True
                linear += excess_area
                density_excess = excess_area / capacity
                squared += capacity * density_excess * density_excess
    gradient = np.zeros_like(centres)
    for node in range(centres.shape[0]):
        if not movable[node]:
            continue
        left = centres[node, 0] - width[node] / 2.0
        right = centres[node, 0] + width[node] / 2.0
        bottom = centres[node, 1] - height[node] / 2.0
        top = centres[node, 1] + height[node] / 2.0
        ix0, ix1 = _cell_range(left, right, x0, dx, nx)
        iy0, iy1 = _cell_range(bottom, top, y0, dy, ny)
        for iy in range(iy0, iy1 + 1):
            bin_bottom, bin_top = y_edges[iy], y_edges[iy + 1]
            overlap_y = min(top, bin_top) - max(bottom, bin_bottom)
            if overlap_y <= 0.0:
                continue
            d_overlap_y = 0.0
            if top < bin_top:
                d_overlap_y += 1.0
            if bottom > bin_bottom:
                d_overlap_y -= 1.0
            for ix in range(ix0, ix1 + 1):
                if not active[iy, ix]:
                    continue
                bin_left, bin_right = x_edges[ix], x_edges[ix + 1]
                overlap_x = min(right, bin_right) - max(left, bin_left)
                if overlap_x <= 0.0:
                    continue
                d_overlap_x = 0.0
                if right < bin_right:
                    d_overlap_x += 1.0
                if left > bin_left:
                    d_overlap_x -= 1.0
                gradient[node, 0] += d_overlap_x * overlap_y
                gradient[node, 1] += overlap_x * d_overlap_y
    return linear, squared, gradient, occupancy, active


@njit(cache=True)
def strict_density_value_gradient(
    centres: np.ndarray, width: np.ndarray, height: np.ndarray, movable: np.ndarray,
    x_edges: np.ndarray, y_edges: np.ndarray, available: np.ndarray, rho_target: float,
) -> tuple[float, float, np.ndarray, np.ndarray, np.ndarray, float]:
    """Return the strict capacity-overflow value and a Clarke subgradient.

    This matches :func:`density_value_gradient` on bins with usable capacity,
    while additionally charging every movable-cell area that falls in a bin
    whose usable area is zero.  The legacy oracle deliberately skips those
    bins for backwards compatibility; callers that make feasibility or search
    decisions should use this strict variant instead.
    """
    ny, nx = available.shape
    occupancy = np.zeros((ny, nx), dtype=np.float64)
    x0, y0 = x_edges[0], y_edges[0]
    dx, dy = x_edges[1] - x_edges[0], y_edges[1] - y_edges[0]
    for node in range(centres.shape[0]):
        if not movable[node]:
            continue
        left = centres[node, 0] - width[node] / 2.0
        right = centres[node, 0] + width[node] / 2.0
        bottom = centres[node, 1] - height[node] / 2.0
        top = centres[node, 1] + height[node] / 2.0
        ix0, ix1 = _cell_range(left, right, x0, dx, nx)
        iy0, iy1 = _cell_range(bottom, top, y0, dy, ny)
        if ix0 > ix1 or iy0 > iy1:
            continue
        for iy in range(iy0, iy1 + 1):
            overlap_y = min(top, y_edges[iy + 1]) - max(bottom, y_edges[iy])
            if overlap_y <= 0.0:
                continue
            for ix in range(ix0, ix1 + 1):
                overlap_x = min(right, x_edges[ix + 1]) - max(left, x_edges[ix])
                if overlap_x > 0.0:
                    occupancy[iy, ix] += overlap_x * overlap_y

    active = np.zeros((ny, nx), dtype=np.bool_)
    linear = 0.0
    squared = 0.0
    zero_capacity_occupancy = 0.0
    for iy in range(ny):
        for ix in range(nx):
            capacity = available[iy, ix]
            if capacity <= 0.0:
                if occupancy[iy, ix] > 0.0:
                    active[iy, ix] = True
                    linear += occupancy[iy, ix]
                    zero_capacity_occupancy += occupancy[iy, ix]
                continue
            excess_area = occupancy[iy, ix] - rho_target * capacity
            if excess_area > 0.0:
                active[iy, ix] = True
                linear += excess_area
                density_excess = excess_area / capacity
                squared += capacity * density_excess * density_excess

    gradient = np.zeros_like(centres)
    for node in range(centres.shape[0]):
        if not movable[node]:
            continue
        left = centres[node, 0] - width[node] / 2.0
        right = centres[node, 0] + width[node] / 2.0
        bottom = centres[node, 1] - height[node] / 2.0
        top = centres[node, 1] + height[node] / 2.0
        ix0, ix1 = _cell_range(left, right, x0, dx, nx)
        iy0, iy1 = _cell_range(bottom, top, y0, dy, ny)
        for iy in range(iy0, iy1 + 1):
            bin_bottom, bin_top = y_edges[iy], y_edges[iy + 1]
            overlap_y = min(top, bin_top) - max(bottom, bin_bottom)
            if overlap_y <= 0.0:
                continue
            d_overlap_y = 0.0
            if top < bin_top:
                d_overlap_y += 1.0
            if bottom > bin_bottom:
                d_overlap_y -= 1.0
            for ix in range(ix0, ix1 + 1):
                if not active[iy, ix]:
                    continue
                bin_left, bin_right = x_edges[ix], x_edges[ix + 1]
                overlap_x = min(right, bin_right) - max(left, bin_left)
                if overlap_x <= 0.0:
                    continue
                d_overlap_x = 0.0
                if right < bin_right:
                    d_overlap_x += 1.0
                if left > bin_left:
                    d_overlap_x -= 1.0
                gradient[node, 0] += d_overlap_x * overlap_y
                gradient[node, 1] += overlap_x * d_overlap_y
    return linear, squared, gradient, occupancy, active, zero_capacity_occupancy


@njit(cache=True)
def strict_density_squared_gradient(
    centres: np.ndarray, width: np.ndarray, height: np.ndarray, movable: np.ndarray,
    x_edges: np.ndarray, y_edges: np.ndarray, available: np.ndarray, rho_target: float,
    occupancy: np.ndarray,
) -> np.ndarray:
    """Gradient of squared positive-capacity excess plus linear blocked occupancy.

    Positive-capacity bins are weighted by ``2*excess/capacity`` so gradients
    do not cancel merely because every neighbouring bin is active. Blocked bins
    retain the strict linear-area term because normalized squared capacity is
    undefined when available area is zero.
    """
    ny, nx = available.shape
    x0, y0 = x_edges[0], y_edges[0]
    dx, dy = x_edges[1] - x_edges[0], y_edges[1] - y_edges[0]
    gradient = np.zeros_like(centres)
    for node in range(centres.shape[0]):
        if not movable[node]:
            continue
        left = centres[node, 0] - width[node] / 2.0
        right = centres[node, 0] + width[node] / 2.0
        bottom = centres[node, 1] - height[node] / 2.0
        top = centres[node, 1] + height[node] / 2.0
        ix0, ix1 = _cell_range(left, right, x0, dx, nx)
        iy0, iy1 = _cell_range(bottom, top, y0, dy, ny)
        for iy in range(iy0, iy1 + 1):
            bin_bottom, bin_top = y_edges[iy], y_edges[iy + 1]
            overlap_y = min(top, bin_top) - max(bottom, bin_bottom)
            if overlap_y <= 0.0:
                continue
            d_overlap_y = 0.0
            if top < bin_top:
                d_overlap_y += 1.0
            if bottom > bin_bottom:
                d_overlap_y -= 1.0
            for ix in range(ix0, ix1 + 1):
                capacity = available[iy, ix]
                if capacity <= 0.0:
                    if occupancy[iy, ix] <= 0.0:
                        continue
                    weight = 1.0
                else:
                    excess = occupancy[iy, ix] - rho_target * capacity
                    if excess <= 0.0:
                        continue
                    weight = 2.0 * excess / capacity
                bin_left, bin_right = x_edges[ix], x_edges[ix + 1]
                overlap_x = min(right, bin_right) - max(left, bin_left)
                if overlap_x <= 0.0:
                    continue
                d_overlap_x = 0.0
                if right < bin_right:
                    d_overlap_x += 1.0
                if left > bin_left:
                    d_overlap_x -= 1.0
                gradient[node, 0] += weight * d_overlap_x * overlap_y
                gradient[node, 1] += weight * overlap_x * d_overlap_y
    return gradient
