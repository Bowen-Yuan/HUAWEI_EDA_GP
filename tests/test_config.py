from __future__ import annotations

from adaptive_pareto.config import AdaptiveConfig


def test_default_adaptive_config_is_valid() -> None:
    config = AdaptiveConfig()
    config.validate()
    assert config.evaluation.rho_target == 0.8
    assert config.visualization.html
