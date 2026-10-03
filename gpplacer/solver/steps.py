"""Projection, preconditioning, and active-set comparisons used by all stages."""

from __future__ import annotations

import numpy as np

from gpplacer.model.types import PlacementDB


def build_preconditioner(db: PlacementDB) -> np.ndarray:
    """Build a bounded diagonal scale from weighted incidence and cell area."""
    degree = np.zeros(db.node_count, dtype=np.float64)
    for net in range(db.net_count):
        nodes = np.unique(db.pin_node[db.net_start[net]:db.net_start[net + 1]])
        degree[nodes] += db.net_weight[net]
    area_scale = db.width * db.height / (db.row_height * db.row_height)
    diagonal = 1.0 + degree + 0.1 * area_scale
    return np.clip(diagonal, 1.0, 1.0e4)


def project_centres(db: PlacementDB, centres: np.ndarray) -> np.ndarray:
    """Project movable centres into the global legal-row bounding rectangle."""
    projected = centres.copy()
    left, bottom, right, top = db.core_bounds
    movable = db.movable
    projected[movable, 0] = np.clip(
        projected[movable, 0], left + db.width[movable] / 2.0,
        right - db.width[movable] / 2.0,
    )
    projected[movable, 1] = np.clip(
        projected[movable, 1], bottom + db.height[movable] / 2.0,
        top - db.height[movable] / 2.0,
    )
    projected[db.fixed] = db.initial_centres[db.fixed]
    return projected


def limited_descent(
    db: PlacementDB, centres: np.ndarray, gradient: np.ndarray, preconditioner: np.ndarray,
    step: float, max_displacement: float, mask: np.ndarray | None = None,
) -> np.ndarray:
    """Apply scaled subgradient descent with an infinity-norm displacement cap."""
    direction = -gradient / preconditioner[:, None]
    if mask is not None:
        direction = direction.copy()
        direction[~mask] = 0.0
    direction[db.fixed] = 0.0
    displacement = step * direction
    maximum = float(np.max(np.abs(displacement)))
    if maximum > max_displacement:
        displacement *= max_displacement / maximum
    return project_centres(db, centres + displacement)


def active_change(previous_extrema: np.ndarray, previous_bins: np.ndarray,
                  current_extrema: np.ndarray, current_bins: np.ndarray) -> tuple[float, float]:
    """Return report-style boundary-net and overflow-bin change rates."""
    net_change = float(np.mean(np.any(previous_extrema != current_extrema, axis=1)))
    old = previous_bins.ravel()
    new = current_bins.ravel()
    denominator = max(1, int(np.count_nonzero(old)))
    bin_change = float(np.count_nonzero(old != new) / denominator)
    return net_change, bin_change
