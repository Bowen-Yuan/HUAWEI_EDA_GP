"""Configuration for the safe-anchor solver."""
from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path


@dataclass(slots=True)
class HardSolverConfig:
    rho_target: float = 0.8
    bin_rows: int = 8
    max_iterations: int = 120
    max_seconds: float = 600.0
    initial_step: float = 1.0
    min_step: float = 1.0 / 64.0
    trust_initial_rows: float = 2.0
    trust_min_rows: float = 0.25
    trust_max_rows: float = 4.0
    trust_expand: float = 1.2
    trust_shrink: float = 0.5
    target_overflow_percent: float | None = None
    # The PDF-style squared density is useful for diagnosis, but differs in
    # scale and normalisation from the capacity-aware linear overflow gate.
    # Keep it opt-in as a second hard constraint.
    enforce_challenge_square_budget: bool = False
    # Feasibility construction is deliberately independent of HPWL: reach the
    # requested overflow target first, then preserve it as a hard constraint.
    anchor_hpwl_ratio: float | None = None
    anchor_repair_rounds: int = 24
    anchor_evacuation_rounds: int = 20
    anchor_max_moves_per_round: int = 10_000
    anchor_max_radius: int = 112
    repair_candidates: int = 8
    repair_target_bins: int = 8
    bundle_size: int = 4
    stall_limit: int = 6
    pressure_percent_points: tuple[float, ...] = (0.10, 0.25, 0.50)
    pressure_steps: int = 2
    pressure_repair_rounds: int = 8
    snapshot_every: int = 10
    seed: int = 17
    output_root: str = "code-hard/runs"

    @classmethod
    def from_json(cls, path: str | Path) -> "HardSolverConfig":
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
        known = set(cls.__dataclass_fields__)
        unknown = set(raw) - known
        if unknown:
            raise ValueError(f"Unknown configuration keys: {sorted(unknown)}")
        if "pressure_percent_points" in raw:
            raw["pressure_percent_points"] = tuple(raw["pressure_percent_points"])
        return cls(**raw)

    def to_dict(self) -> dict[str, object]:
        return asdict(self)
