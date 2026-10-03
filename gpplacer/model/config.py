"""Configuration objects for the placement solver.

Defaults are deliberately case-independent.  Tune them through a JSON file
rather than adding per-case constants to solver code.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
import json


@dataclass(slots=True)
class SolverConfig:
    """Numerical, budget, and observability controls for one solver run."""

    rho_target: float = 0.8
    bin_rows: int = 8
    density_lambda: float | None = None
    coarse_iterations: int = 2
    coarsen_degree_limit: int = 64
    initialization_rounds: int = 16
    initialization_strategy: str = "connection"
    cluster_macro_bins: int = 16
    partition_leaf_bins: int = 16
    partition_refine_passes: int = 2
    enable_bin_evacuator: bool = True
    evacuation_rounds: int = 2
    evacuation_max_moves_per_round: int = 10_000
    evacuation_max_radius: int = 112
    max_iterations: int = 200
    max_seconds: float = 600.0
    max_memory_mb: int | None = None
    threads: int = 1
    seed: int = 17
    snapshot_every: int = 10
    max_displacement_rows: float = 4.0
    min_displacement_rows: float = 0.25
    initial_step: float = 1.0
    explore_window: int = 8
    switch_rel_improvement: float = 1e-4
    switch_active_change: float = 0.05
    require_directional_instability: bool = False
    bundle_size: int = 8
    bundle_null_limit: int = 6
    bundle_overflow_percent_limit: float = 15.0
    bundle_min_proximal_scale_rows: float = 0.05
    serious_ratio: float = 0.10
    proximal_scale: float = 1.0
    active_stable_change: float = 0.01
    active_stable_rounds: int = 5
    enable_refine: bool = True
    refine_fraction: float = 0.15
    full_oracle_period: int = 10
    output_root: str = "runs"

    @classmethod
    def from_json(cls, path: str | Path) -> "SolverConfig":
        """Load a config file while rejecting unknown keys early."""
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        known = set(cls.__dataclass_fields__)
        unknown = set(payload) - known
        if unknown:
            raise ValueError(f"Unknown solver configuration keys: {sorted(unknown)}")
        return cls(**payload)

    def to_dict(self) -> dict[str, object]:
        """Return JSON-serialisable configuration data."""
        return asdict(self)
