"""Deterministic state selection for guarded single-objective rescue blocks."""

from __future__ import annotations

from enum import Enum

from adaptive_pareto.config import RescueConfig


class SolverPhase(str, Enum):
    HPWL_RESCUE = "hpwl_rescue"
    DENSITY_RESCUE = "density_rescue"
    DUAL_EMERGENCY = "dual_emergency"
    JOINT = "joint"
    PLATEAU_ESCAPE = "plateau_escape"
    REFINE = "refine"


def choose_phase(normalized_hpwl: float, strict_overflow: float, config: RescueConfig) -> SolverPhase:
    high_hpwl = normalized_hpwl > config.hpwl_enter_norm
    high_overflow = strict_overflow > config.overflow_enter_percent
    if high_hpwl and high_overflow:
        return SolverPhase.DUAL_EMERGENCY
    if high_hpwl:
        return SolverPhase.HPWL_RESCUE
    if high_overflow:
        return SolverPhase.DENSITY_RESCUE
    return SolverPhase.JOINT


def dual_priority(normalized_hpwl: float, strict_overflow: float, config: RescueConfig) -> SolverPhase:
    """Choose the more severe objective for one short dual-emergency block."""
    hpwl_span = max(config.hpwl_enter_norm - config.hpwl_exit_norm, 1.0e-12)
    overflow_span = max(config.overflow_enter_percent - config.overflow_exit_percent, 1.0e-12)
    hpwl_severity = (normalized_hpwl - config.hpwl_exit_norm) / hpwl_span
    overflow_severity = (strict_overflow - config.overflow_exit_percent) / overflow_span
    return SolverPhase.HPWL_RESCUE if hpwl_severity >= overflow_severity else SolverPhase.DENSITY_RESCUE
