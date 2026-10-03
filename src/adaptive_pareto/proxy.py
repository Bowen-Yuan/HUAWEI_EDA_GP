"""Tiny dependency-free ridge proxy used only to rank short candidate trials."""

from __future__ import annotations

import numpy as np


class RidgeProxy:
    """Standardised closed-form ridge regression with deterministic fallback."""

    def __init__(self, ridge: float = 1.0e-3) -> None:
        self.ridge = ridge
        self.mean: np.ndarray | None = None
        self.scale: np.ndarray | None = None
        self.weights: np.ndarray | None = None
        self.bias: float = 0.0

    def fit(self, features: np.ndarray, targets: np.ndarray) -> "RidgeProxy":
        self.mean = features.mean(axis=0)
        self.scale = np.maximum(features.std(axis=0), 1.0e-9)
        design = (features - self.mean) / self.scale
        augmented = np.column_stack((np.ones(len(design)), design))
        regularizer = np.eye(augmented.shape[1]) * self.ridge
        regularizer[0, 0] = 0.0
        coefficients = np.linalg.pinv(augmented.T @ augmented + regularizer) @ augmented.T @ targets
        self.bias = float(coefficients[0])
        self.weights = coefficients[1:]
        return self

    def predict(self, features: np.ndarray) -> np.ndarray:
        if self.mean is None or self.scale is None or self.weights is None:
            raise RuntimeError("RidgeProxy must be fitted before prediction.")
        return self.bias + ((features - self.mean) / self.scale) @ self.weights
