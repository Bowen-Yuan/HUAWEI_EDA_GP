"""Validated configuration for the adaptive Pareto placement controller."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
from pathlib import Path
from typing import Any


@dataclass(slots=True)
class EvaluationConfig:
    rho_target: float = 0.8
    bin_rows: int = 8
    target_overflow_low: float = 20.0
    target_overflow_high: float = 30.0


@dataclass(slots=True)
class RescueConfig:
    hpwl_enter_norm: float = 0.15
    hpwl_exit_norm: float = 0.08
    overflow_enter_percent: float = 45.0
    overflow_exit_percent: float = 35.0
    stable_exit_steps: int = 5
    block_steps: int = 20
    max_blocks: int = 10
    hpwl_rescue_overflow_increase_percent: float = 2.0
    density_rescue_hpwl_increase_ratio: float = 0.08
    density_rescue_total_hpwl_ratio: float = 2.0


@dataclass(slots=True)
class DirectionConfig:
    mode: str = "hybrid_axis"
    axis_dominance_ratio: float = 3.0


@dataclass(slots=True)
class EvacuationConfig:
    enabled: bool = True
    max_radius_bins: int = 16
    max_batch_moves: int = 256
    max_move_arrows: int = 2_000


@dataclass(slots=True)
class HierarchyConfig:
    enabled: bool = True
    stop_nodes: int = 12_000
    max_levels: int = 6
    min_reduction_ratio: float = 0.15
    degree_limit: int = 64
    level_refine_steps: tuple[int, ...] = (20, 30, 50, 80)


@dataclass(slots=True)
class PortfolioConfig:
    enabled: bool = True
    heavy_edge_candidates: int = 8
    partition_candidates: int = 8
    capacity_candidates: int = 8
    connection_candidates: int = 4
    spectral_candidates: int = 4
    halving_steps: tuple[int, ...] = (10, 30, 80)
    halving_keep: tuple[int, ...] = (16, 8, 4)
    max_per_generator_fraction: float = 0.5
    proxy_enabled: bool = True
    proxy_ridge: float = 1.0e-3


@dataclass(slots=True)
class JointConfig:
    steps: int = 520
    initial_step: float = 4.0
    min_step: float = 1.0 / 256.0
    trust_rows: float = 4.0
    max_trust_rows: float = 16.0
    lambda_gain: float = 0.7
    lambda_ema_beta: float = 0.9
    lambda_min: float = 1.0e-6
    lambda_max: float = 1.0e6
    lambda_change_min: float = 0.25
    lambda_change_max: float = 4.0
    archive_per_band: int = 4
    archive_maximum: int = 32


@dataclass(slots=True)
class RuntimeConfig:
    max_seconds: float = 1800.0
    max_steps: int | None = None
    seed: int = 17
    checkpoint_every_steps: int = 50


@dataclass(slots=True)
class VisualizationConfig:
    enabled: bool = True
    static_formats: tuple[str, ...] = ("png", "pdf")
    html: bool = True
    snapshot_every_accepted: int = 25
    snapshot_every_seconds: float = 60.0
    max_periodic_snapshots: int = 40
    max_cells_static: int = 100_000
    max_cells_html: int = 100_000
    max_trajectory_points_html: int = 5_000
    max_move_arrows: int = 2_000
    dpi: int = 300


@dataclass(slots=True)
class AdaptiveConfig:
    evaluation: EvaluationConfig = field(default_factory=EvaluationConfig)
    rescue: RescueConfig = field(default_factory=RescueConfig)
    direction: DirectionConfig = field(default_factory=DirectionConfig)
    evacuation: EvacuationConfig = field(default_factory=EvacuationConfig)
    hierarchy: HierarchyConfig = field(default_factory=HierarchyConfig)
    portfolio: PortfolioConfig = field(default_factory=PortfolioConfig)
    joint: JointConfig = field(default_factory=JointConfig)
    runtime: RuntimeConfig = field(default_factory=RuntimeConfig)
    visualization: VisualizationConfig = field(default_factory=VisualizationConfig)
    output_root: str = "adaptive_pareto_soft/results"

    @classmethod
    def from_json(cls, path: str | Path) -> "AdaptiveConfig":
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        return cls.from_mapping(payload)

    @classmethod
    def from_mapping(cls, payload: dict[str, Any]) -> "AdaptiveConfig":
        if not isinstance(payload, dict):
            raise ValueError("Adaptive configuration must be a JSON object.")
        known = set(cls.__dataclass_fields__)
        unknown = set(payload) - known
        if unknown:
            raise ValueError(f"Unknown adaptive configuration keys: {sorted(unknown)}")
        config = cls()
        for name, value in payload.items():
            current = getattr(config, name)
            if hasattr(current, "__dataclass_fields__"):
                if not isinstance(value, dict):
                    raise ValueError(f"{name} must be an object.")
                nested_unknown = set(value) - set(current.__dataclass_fields__)
                if nested_unknown:
                    raise ValueError(f"Unknown {name} keys: {sorted(nested_unknown)}")
                for key, nested_value in value.items():
                    if isinstance(getattr(current, key), tuple):
                        nested_value = tuple(nested_value)
                    setattr(current, key, nested_value)
            else:
                setattr(config, name, value)
        config.validate()
        return config

    def validate(self) -> None:
        if not 0.0 < self.evaluation.rho_target <= 1.0:
            raise ValueError("evaluation.rho_target must lie in (0, 1].")
        if self.evaluation.bin_rows < 1:
            raise ValueError("evaluation.bin_rows must be positive.")
        if not 0.0 <= self.evaluation.target_overflow_low <= self.evaluation.target_overflow_high:
            raise ValueError("Target overflow band is invalid.")
        if self.direction.mode not in {"full_2d", "dominant_axis", "hybrid_axis"}:
            raise ValueError("direction.mode must be full_2d, dominant_axis, or hybrid_axis.")
        if self.direction.axis_dominance_ratio < 1.0:
            raise ValueError("direction.axis_dominance_ratio must be at least one.")
        if len(self.portfolio.halving_steps) != len(self.portfolio.halving_keep):
            raise ValueError("portfolio.halving_steps and halving_keep must have the same length.")
        if any(value <= 0 for value in self.portfolio.halving_steps + self.portfolio.halving_keep):
            raise ValueError("Portfolio budgets and survivor counts must be positive.")
        if self.runtime.max_seconds <= 0.0:
            raise ValueError("runtime.max_seconds must be positive.")
        if self.runtime.max_steps is not None and self.runtime.max_steps <= 0:
            raise ValueError("runtime.max_steps must be positive when supplied.")

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)
