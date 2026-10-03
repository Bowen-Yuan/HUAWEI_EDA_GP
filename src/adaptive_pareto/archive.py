"""Stratified strict-overflow Pareto archive with recoverable placements."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from gpplacer.oracle.evaluation import DualOracleResult


@dataclass(slots=True)
class ArchivePoint:
    centres: np.ndarray
    hpwl: float
    strict_overflow_percent: float
    legacy_overflow_percent: float
    phase: str
    iteration: int


class StratifiedParetoArchive:
    """Keep strict-nondominated points while preserving coverage by overflow band."""

    def __init__(self, *, band_width_percent: float = 5.0, per_band: int = 4, maximum: int = 32) -> None:
        self.band_width_percent = band_width_percent
        self.per_band = per_band
        self.maximum = maximum
        self.points: list[ArchivePoint] = []

    @staticmethod
    def _dominates(left: ArchivePoint, right: ArchivePoint) -> bool:
        return (
            left.hpwl <= right.hpwl
            and left.strict_overflow_percent <= right.strict_overflow_percent
            and (left.hpwl < right.hpwl or left.strict_overflow_percent < right.strict_overflow_percent)
        )

    def consider(
        self, centres: np.ndarray, result: DualOracleResult, *, strict_overflow_percent: float,
        legacy_overflow_percent: float, phase: str, iteration: int,
    ) -> bool:
        point = ArchivePoint(
            centres.copy(), result.hpwl, strict_overflow_percent, legacy_overflow_percent, phase, iteration,
        )
        if any(self._dominates(old, point) or (
            old.hpwl == point.hpwl and old.strict_overflow_percent == point.strict_overflow_percent
        ) for old in self.points):
            return False
        self.points = [old for old in self.points if not self._dominates(point, old)]
        self.points.append(point)
        self._trim()
        return True

    def _trim(self) -> None:
        if not self.points:
            return
        kept: list[ArchivePoint] = []
        bands: dict[int, list[ArchivePoint]] = {}
        for point in self.points:
            key = int(np.floor(point.strict_overflow_percent / self.band_width_percent))
            bands.setdefault(key, []).append(point)
        for points in bands.values():
            points.sort(key=lambda item: item.hpwl)
            kept.extend(points[:self.per_band])
        if len(kept) > self.maximum:
            ordered = sorted(kept, key=lambda item: item.strict_overflow_percent)
            selected = np.linspace(0, len(ordered) - 1, self.maximum).round().astype(int)
            kept = [ordered[index] for index in np.unique(selected)]
        self.points = kept

    def target_representative(self, low: float, high: float) -> ArchivePoint:
        if not self.points:
            raise ValueError("Cannot select a representative from an empty archive.")
        def key(point: ArchivePoint) -> tuple[float, float]:
            distance = max(low - point.strict_overflow_percent, point.strict_overflow_percent - high, 0.0)
            return distance, point.hpwl
        return min(self.points, key=key)
