"""Overflow-safe HPWL compression with temporary pressure and pair exchanges.

This is deliberately separate from the default Explore/Bundle/Refine solver.
It starts from a complete placement, retains a *safe anchor* whose linear
overflow never exceeds the starting value, and discards every pressure probe
that cannot be repaired back to that bound.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import csv
from collections.abc import Callable

import numpy as np

from gpplacer.io.placement import read_placement_centres
from gpplacer.model.types import PlacementDB
from gpplacer.oracle.density import density_value_gradient
from gpplacer.oracle.grid import DensityGrid
from gpplacer.oracle.objective import ObjectiveOracle, OracleResult
from gpplacer.solver.steps import build_preconditioner, project_centres


@dataclass(frozen=True, slots=True)
class AnchoredSearchConfig:
    """Controls for the isolated anchor-and-pressure experiment."""

    max_iterations: int = 200
    initial_step: float = 1.0
    max_displacement_rows: float = 4.0
    min_step: float = 1.0 / 64.0
    stall_limit: int = 8
    pressure_levels: int = 3
    pressure_steps: int = 6
    pressure_increment_area: float = 10_000.0
    exchange_source_cells: int = 8
    exchange_target_cells: int = 32
    exchange_pair_limit: int = 24
    hpwl_improvement_epsilon: float = 1e-6
    overflow_percent_limit: float | None = None


@dataclass(slots=True)
class AnchoredSearchResult:
    """Best safe placement and diagnostics from an anchored search."""

    centres: np.ndarray
    oracle: OracleResult
    records: list[dict[str, object]]
    pressure_probes: int
    accepted_exchanges: int


def select_lowest_overflow_solution(
    summary_csv: str | Path, overflow_percent_limit: float | None = None,
) -> Path:
    """Select an anchor, optionally preferring the shortest solution below a cap."""
    summary_path = Path(summary_csv)
    run_root = summary_path.parent
    with summary_path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    candidates: list[tuple[float, float, Path]] = []
    for row in rows:
        solution = run_root / row["run_id"] / "solution.pl"
        if solution.exists():
            overflow_percent = float(row["overflow_percent"])
            hpwl = float(row["hpwl"])
            if overflow_percent_limit is None or overflow_percent <= overflow_percent_limit:
                candidates.append((overflow_percent, hpwl, solution))
    if not candidates:
        raise ValueError(f"No readable solution within the requested overflow cap in {summary_path}")
    if overflow_percent_limit is None:
        return min(candidates, key=lambda item: item[0])[2]
    return min(candidates, key=lambda item: item[1])[2]


class AnchoredSearch:
    """Greedy HPWL compression constrained by an immutable overflow budget."""

    def __init__(self, db: PlacementDB, grid: DensityGrid, config: AnchoredSearchConfig) -> None:
        self.db = db
        self.grid = grid
        self.config = config
        # ``density_lambda=0`` provides the exact HPWL subgradient while the
        # oracle still evaluates density for the hard acceptance test.
        self.oracle = ObjectiveOracle(db, grid, density_lambda=0.0)
        self.preconditioner = build_preconditioner(db)
        self.total_available_area = float(grid.available_area.sum())

    @staticmethod
    def config_dict(config: AnchoredSearchConfig) -> dict[str, object]:
        return asdict(config)

    def run(
        self,
        anchor_centres: np.ndarray,
        observer: Callable[[dict[str, object], np.ndarray, OracleResult], None] | None = None,
    ) -> AnchoredSearchResult:
        safe_centres = anchor_centres.copy()
        safe_oracle = self.oracle.evaluate(safe_centres)
        safe_budget = safe_oracle.density_linear
        if self.config.overflow_percent_limit is not None:
            if not 0.0 < self.config.overflow_percent_limit <= 100.0:
                raise ValueError("overflow_percent_limit must lie in (0, 100].")
            safe_budget = self.total_available_area * self.config.overflow_percent_limit / 100.0
            if safe_oracle.density_linear > safe_budget + 1e-9:
                raise ValueError("The supplied anchor already violates overflow_percent_limit.")
        step = self.config.initial_step
        records: list[dict[str, object]] = []
        stalled = 0
        pressure_probes = 0
        accepted_exchanges = 0

        for iteration in range(1, self.config.max_iterations + 1):
            candidate, candidate_oracle, accepted, backtracks = self._compress_once(
                safe_centres, safe_oracle, safe_budget, step,
            )
            if accepted:
                safe_centres, safe_oracle = candidate, candidate_oracle
                stalled = 0
                step = min(step * 1.1, self.config.initial_step)
                event = "compress_accept"
            else:
                stalled += 1
                step = max(step * 0.5, self.config.min_step)
                event = "compress_reject"

            if stalled >= self.config.stall_limit:
                repaired, probe_count = self._pressure_and_exchange(safe_centres, safe_oracle, safe_budget)
                pressure_probes += probe_count
                if repaired is not None:
                    safe_centres, safe_oracle = repaired
                    accepted_exchanges += 1
                    stalled = 0
                    step = self.config.initial_step
                    event = "exchange_accept"
                else:
                    # No temporary solution is ever retained after a failed repair.
                    stalled = 0
                    event = "pressure_discard"

            record = self._record(
                iteration=iteration, event=event, oracle=safe_oracle,
                step=step, backtracks=backtracks, safe_budget=safe_budget,
            )
            records.append(record)
            if observer is not None:
                observer(record, safe_centres, safe_oracle)
        return AnchoredSearchResult(
            centres=safe_centres, oracle=safe_oracle, records=records,
            pressure_probes=pressure_probes, accepted_exchanges=accepted_exchanges,
        )

    def _compress_once(
        self, centres: np.ndarray, oracle: OracleResult, density_budget: float, step: float,
        *, constrain_density_direction: bool = True,
    ) -> tuple[np.ndarray, OracleResult, bool, int]:
        """Backtrack an HPWL-only move until it respects the hard overflow budget."""
        trial_step = step
        for backtracks in range(9):
            use_tangent_projection = (
                constrain_density_direction
                and oracle.density_linear >= density_budget * 0.995
            )
            candidate = self._feasible_descent_candidate(
                centres, oracle, trial_step, use_tangent_projection,
            )
            candidate_oracle = self.oracle.evaluate(candidate)
            if (
                candidate_oracle.density_linear <= density_budget + 1e-9
                and candidate_oracle.hpwl < oracle.hpwl - self.config.hpwl_improvement_epsilon
            ):
                return candidate, candidate_oracle, True, backtracks
            trial_step *= 0.5
        return centres, oracle, False, 8

    def _feasible_descent_candidate(
        self, centres: np.ndarray, oracle: OracleResult, step: float, constrain_density_direction: bool,
    ) -> np.ndarray:
        """Project an HPWL descent direction into the linearised density cone.

        The pure HPWL direction is ``-g_hpwl / M``.  If it would increase the
        exact linear-density model to first order, its component along
        ``g_density / M`` is removed so that ``g_density.T @ direction <= 0``.
        Exact evaluation remains the final, non-linear hard-constraint check.
        """
        direction = -oracle.gradient / self.preconditioner[:, None]
        if constrain_density_direction:
            _, _, density_gradient, _, _ = density_value_gradient(
                centres, self.db.width, self.db.height, self.db.movable,
                self.grid.x_edges, self.grid.y_edges, self.grid.available_area,
                self.grid.rho_target,
            )
            density_direction = density_gradient / self.preconditioner[:, None]
            directional_density_change = float(np.sum(density_gradient * direction))
            density_norm_sq = float(np.sum(density_gradient * density_direction))
            if directional_density_change > 0.0 and density_norm_sq > 1e-12:
                direction -= (directional_density_change / density_norm_sq) * density_direction
        direction[self.db.fixed] = 0.0
        displacement = step * direction
        maximum = float(np.max(np.abs(displacement)))
        limit = self.config.max_displacement_rows * self.db.row_height
        if maximum > limit:
            displacement *= limit / maximum
        return project_centres(self.db, centres + displacement)

    def _pressure_and_exchange(
        self, safe_centres: np.ndarray, safe_oracle: OracleResult, safe_budget: float,
    ) -> tuple[tuple[np.ndarray, OracleResult] | None, int]:
        """Explore relaxed budgets, then require an exact exchange-based repair."""
        probes = 0
        for level in range(1, self.config.pressure_levels + 1):
            relaxed_budget = safe_budget + level * self.config.pressure_increment_area
            probe_centres, probe_oracle = safe_centres, safe_oracle
            for _ in range(self.config.pressure_steps):
                candidate, candidate_oracle, accepted, _ = self._compress_once(
                    probe_centres, probe_oracle, relaxed_budget, self.config.initial_step,
                    constrain_density_direction=False,
                )
                probes += 1
                if not accepted:
                    break
                probe_centres, probe_oracle = candidate, candidate_oracle
                if candidate_oracle.density_linear <= safe_budget + 1e-9:
                    # This is already safe and shorter; it is a valid new anchor.
                    return (candidate, candidate_oracle), probes
            repaired = self._best_safe_exchange(probe_centres, probe_oracle, safe_oracle, safe_budget)
            if repaired is not None:
                return repaired, probes
        return None, probes

    def _best_safe_exchange(
        self, probe_centres: np.ndarray, probe_oracle: OracleResult,
        safe_oracle: OracleResult, safe_budget: float,
    ) -> tuple[np.ndarray, OracleResult] | None:
        """Try bounded area-similar cell swaps and accept only exact safe improvements."""
        source_nodes, target_nodes = self._exchange_node_sets(probe_centres, probe_oracle)
        best: tuple[np.ndarray, OracleResult] | None = None
        evaluated = 0
        for source in source_nodes:
            ordered_targets = sorted(
                target_nodes,
                key=lambda target: abs(
                    self.db.width[source] * self.db.height[source]
                    - self.db.width[target] * self.db.height[target]
                ),
            )
            for target in ordered_targets:
                if source == target:
                    continue
                candidate = probe_centres.copy()
                candidate[source], candidate[target] = probe_centres[target], probe_centres[source]
                candidate_oracle = self.oracle.evaluate(candidate)
                evaluated += 1
                if (
                    candidate_oracle.density_linear <= safe_budget + 1e-9
                    and candidate_oracle.hpwl < safe_oracle.hpwl - self.config.hpwl_improvement_epsilon
                    and (best is None or candidate_oracle.hpwl < best[1].hpwl)
                ):
                    best = (candidate, candidate_oracle)
                if evaluated >= self.config.exchange_pair_limit:
                    return best
        return best

    def _exchange_node_sets(
        self, centres: np.ndarray, oracle: OracleResult,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Choose cells in the most overfull and least occupied bins deterministically."""
        available = self.grid.available_area
        occupancy_ratio = np.divide(
            oracle.occupancy, available, out=np.full_like(available, np.inf), where=available > 0.0,
        )
        source_bins = np.flatnonzero(oracle.active_bins.ravel())
        target_bins = np.flatnonzero((available > 0.0).ravel())
        source_bins = source_bins[np.argsort(occupancy_ratio.ravel()[source_bins])[::-1]]
        target_bins = target_bins[np.argsort(occupancy_ratio.ravel()[target_bins])]
        source = self._nodes_in_bins(centres, source_bins, self.config.exchange_source_cells)
        target = self._nodes_in_bins(centres, target_bins, self.config.exchange_target_cells)
        return source, target

    def _nodes_in_bins(self, centres: np.ndarray, flat_bins: np.ndarray, limit: int) -> np.ndarray:
        nx = self.grid.nx
        dx, dy = self.grid.bin_width, self.grid.bin_height
        x0, y0 = self.grid.x_edges[0], self.grid.y_edges[0]
        ix = np.clip(((centres[:, 0] - x0) // dx).astype(np.int64), 0, nx - 1)
        iy = np.clip(((centres[:, 1] - y0) // dy).astype(np.int64), 0, self.grid.ny - 1)
        node_bins = iy * nx + ix
        chosen: list[int] = []
        for flat_bin in flat_bins:
            nodes = np.flatnonzero(self.db.movable & (node_bins == flat_bin))
            if nodes.size:
                # Larger cells are more likely to repair an overfull bin.
                areas = self.db.width[nodes] * self.db.height[nodes]
                for node in nodes[np.argsort(areas)[::-1]]:
                    chosen.append(int(node))
                    if len(chosen) >= limit:
                        return np.asarray(chosen, dtype=np.int64)
        return np.asarray(chosen, dtype=np.int64)

    def _record(
        self, *, iteration: int, event: str, oracle: OracleResult, step: float,
        backtracks: int, safe_budget: float,
    ) -> dict[str, object]:
        return {
            "iteration": iteration,
            "phase": event,
            "hpwl": oracle.hpwl,
            "density_linear": oracle.density_linear,
            "overflow_ratio": oracle.density_linear / self.total_available_area,
            "overflow_percent": 100.0 * oracle.density_linear / self.total_available_area,
            "density_surrogate_squared": oracle.density_surrogate_squared,
            "objective": oracle.hpwl,
            "step": step,
            "proximal_scale": "",
            "active_nets_change": "",
            "active_bins_change": "",
            "active_bins": oracle.active_bin_count,
            "bundle_size": "",
            "backtracks": backtracks,
            "serious_step": "",
            "predicted_decrease": "",
            "recoveries": "",
            "iteration_seconds": "",
            "oracle_seconds": "",
            "oracle_calls": "",
            "direction_cosine": "",
            "safe_overflow_budget": safe_budget,
        }


def load_anchor(path: str | Path, db: PlacementDB) -> np.ndarray:
    """Load a complete Bookshelf placement as a centre-coordinate anchor."""
    return read_placement_centres(path, db)
