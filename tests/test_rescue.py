from __future__ import annotations

from adaptive_pareto.config import RescueConfig
from adaptive_pareto.rescue import SolverPhase, choose_phase, dual_priority


def test_state_machine_selects_expected_rescue_modes() -> None:
    config = RescueConfig()
    assert choose_phase(.20, 20.0, config) is SolverPhase.HPWL_RESCUE
    assert choose_phase(.02, 50.0, config) is SolverPhase.DENSITY_RESCUE
    assert choose_phase(.20, 50.0, config) is SolverPhase.DUAL_EMERGENCY
    assert choose_phase(.02, 20.0, config) is SolverPhase.JOINT
    assert dual_priority(.20, 50.0, config) is SolverPhase.HPWL_RESCUE
