"""Density-direction policies with fair per-cell displacement normalisation."""

from __future__ import annotations

import numpy as np


def density_direction(
    gradient: np.ndarray,
    preconditioner: np.ndarray,
    trust: np.ndarray,
    movable: np.ndarray,
    *, mode: str, axis_dominance_ratio: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return a descent direction and a per-cell single-axis mask.

    Axis decisions use the same preconditioned, trust-limited first-order
    quantity that later controls displacement.  Returned directions are scaled
    to unit L2 norm for cells with a non-zero density direction, so comparing a
    coordinate step with a two-dimensional step cannot gain distance merely by
    having two active coordinates.
    """
    if mode not in {"full_2d", "dominant_axis", "hybrid_axis"}:
        raise ValueError(f"Unknown density direction mode {mode!r}")
    direction = -gradient / preconditioner[:, None]
    direction = direction.copy()
    direction[~movable] = 0.0
    scores = np.abs(direction) * trust[:, None]
    x_score, y_score = scores[:, 0], scores[:, 1]
    minimum = np.minimum(x_score, y_score)
    maximum = np.maximum(x_score, y_score)
    ratio = maximum / np.maximum(minimum, 1.0e-18)
    choose_axis = (mode == "dominant_axis") | ((mode == "hybrid_axis") & (ratio > axis_dominance_ratio))
    single_axis = choose_axis & movable & (maximum > 0.0)
    choose_x = single_axis & (x_score >= y_score)
    choose_y = single_axis & ~choose_x
    direction[choose_x, 1] = 0.0
    direction[choose_y, 0] = 0.0
    length = np.linalg.norm(direction, axis=1)
    normalise = movable & (length > 0.0)
    direction[normalise] /= length[normalise, None]
    return direction, single_axis
