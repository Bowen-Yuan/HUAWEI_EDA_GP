"""Mutable solver state and compact Bundle-cut records."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
import numpy as np

from gpplacer.oracle.objective import OracleResult


class Phase(str, Enum):
    EXPLORE = "explore"
    STABILIZE = "stabilize"
    REFINE = "refine"


@dataclass(slots=True)
class BundleCut:
    """One local affine model stored relative to the current serious centre."""

    value_at_center: float
    gradient: np.ndarray
    source_iteration: int


@dataclass(slots=True)
class PlacementState:
    """All mutable numerical state needed for rollback and stage transitions."""

    centres: np.ndarray
    oracle: OracleResult
    best_centres: np.ndarray
    best_oracle: OracleResult
    phase: Phase = Phase.EXPLORE
    iteration: int = 0
    step: float = 1.0
    proximal_scale: float = 1.0
    null_steps: int = 0
    recoveries: int = 0
    stable_rounds: int = 0
    bundle: list[BundleCut] = field(default_factory=list)

    def checkpoint(self) -> tuple[np.ndarray, OracleResult]:
        """Return independent data for rollback after an unsafe candidate."""
        return self.centres.copy(), self.oracle
