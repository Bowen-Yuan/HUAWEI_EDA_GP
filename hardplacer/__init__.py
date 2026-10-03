"""Overflow-safe, native nonsmooth placement control layer."""

from .config import HardSolverConfig
from .solver import HardSolveResult, SafeAnchorSolver

__all__ = ["HardSolverConfig", "HardSolveResult", "SafeAnchorSolver"]
