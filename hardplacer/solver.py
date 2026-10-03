"""Safe-anchor preconditioned subgradient solver without objective smoothing.

The continuous phase minimizes exact HPWL under an immutable exact-overflow
budget.  Discrete, local moves cross bin-boundary plateaus in the overlap
oracle; all proposed states are validated by the same exact oracle.
"""
from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Callable
import heapq
import time

import numpy as np

from gpplacer.model.types import PlacementDB
from gpplacer.oracle.density import density_value_gradient
from gpplacer.oracle.grid import DensityGrid
from gpplacer.oracle.objective import ObjectiveOracle, OracleResult
from gpplacer.solver.steps import build_preconditioner, project_centres

from .config import HardSolverConfig

Observer = Callable[[dict[str, object], np.ndarray, OracleResult, np.ndarray], None]


@dataclass(slots=True)
class HardSolveResult:
    centres: np.ndarray
    oracle: OracleResult
    anchor_centres: np.ndarray
    anchor_oracle: OracleResult
    iterations: int
    elapsed_seconds: float
    anchor_repairs: int
    accepted_pressure_steps: int
    records: list[dict[str, object]]


class SafeAnchorSolver:
    """A compact safe-anchor / trust-region / bundle placement controller."""

    def __init__(self, db: PlacementDB, grid: DensityGrid, config: HardSolverConfig) -> None:
        self.db, self.grid, self.config = db, grid, config
        self.oracle = ObjectiveOracle(db, grid, density_lambda=0.0)
        self.preconditioner = build_preconditioner(db)
        self.available_area = float(grid.available_area.sum())
        self._node_degree = np.diff(db.node_net_start).astype(np.float64)
        # The challenge statement writes rho_b as overlap area / geometric bin
        # area.  Keep that scalar separate from the capacity-aware operational
        # overflow used by the <=6% baseline requirement.
        _, _, _, self._fixed_occupancy, _ = density_value_gradient(
            db.initial_centres, db.width, db.height, db.fixed,
            grid.x_edges, grid.y_edges, np.ones_like(grid.available_area), 0.0,
        )
        self._geometric_bin_area = np.outer(np.diff(grid.y_edges), np.diff(grid.x_edges))

    def _overflow_percent(self, result: OracleResult) -> float:
        return 100.0 * result.density_linear / max(self.available_area, 1e-12)

    def challenge_square_density(self, result: OracleResult) -> float:
        """Return the PDF-style Σ_b(max(rho_b-rho_t, 0))² diagnostic.

        It uses geometric bin area and includes fixed instances.  It is logged
        separately because it is not dimensionally equivalent to the linear,
        usable-capacity-normalised overflow percentage used as the 6% gate.
        """
        rho = (result.occupancy + self._fixed_occupancy) / self._geometric_bin_area
        return float(np.square(np.maximum(rho - self.grid.rho_target, 0.0)).sum())

    def _record(self, iteration: int, event: str, result: OracleResult, budget: float,
                square_budget: float, step: float, trust: np.ndarray, *, serious: bool = False,
                bundle_size: int = 0, pressure: bool = False) -> dict[str, object]:
        movable = trust[self.db.movable]
        return {
            "iteration": iteration, "event": event, "hpwl": result.hpwl,
            "density_linear": result.density_linear,
            "density_surrogate_squared": result.density_surrogate_squared,
            "challenge_square_density": self.challenge_square_density(result),
            "overflow_percent": self._overflow_percent(result),
            "budget_percent": 100.0 * budget / max(self.available_area, 1e-12),
            "challenge_square_budget": square_budget,
            "active_bins": result.active_bin_count, "step": step,
            "trust_p50": float(np.quantile(movable, 0.5)),
            "trust_p90": float(np.quantile(movable, 0.9)),
            "serious_step": int(serious), "bundle_size": bundle_size,
            "pressure": int(pressure),
        }

    def _incident_targets(self, centres: np.ndarray) -> np.ndarray:
        pins = centres[self.db.pin_node] + self.db.pin_offset
        net_sum = np.add.reduceat(pins, self.db.net_start[:-1], axis=0)
        net_mean = net_sum / np.maximum(np.diff(self.db.net_start)[:, None], 1)
        target = np.zeros_like(centres)
        for node in np.flatnonzero(self.db.movable):
            nets = self.db.node_nets[self.db.node_net_start[node]:self.db.node_net_start[node + 1]]
            if len(nets):
                target[node] = net_mean[nets].mean(axis=0)
            else:
                target[node] = centres[node]
        return target

    def _centre_bins(self, centres: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        ix = np.searchsorted(self.grid.x_edges, centres[:, 0], side="right") - 1
        iy = np.searchsorted(self.grid.y_edges, centres[:, 1], side="right") - 1
        return np.clip(ix, 0, self.grid.nx - 1), np.clip(iy, 0, self.grid.ny - 1)

    def capacity_seed(self) -> np.ndarray:
        """Create an independent, capacity-balanced starting placement.

        ISPD Bookshelf ``.pl`` files can contain all movable objects outside the
        legal core.  Instead of clipping that invalid cloud onto one boundary,
        assign cells to legal bin centres by remaining target capacity.  This is
        a deterministic feasibility seed, not a reused historical placement.
        """
        centres = self.db.initial_centres.copy()
        capacity = (self.grid.rho_target * self.grid.available_area).ravel()
        heap = [(-float(value), int(index)) for index, value in enumerate(capacity) if value > 0.0]
        heapq.heapify(heap)
        area = self.db.width * self.db.height
        movable = np.flatnonzero(self.db.movable)
        movable = movable[np.argsort(area[movable])[::-1]]
        for node in movable:
            if not heap:
                break
            remaining, flat = heapq.heappop(heap)
            remaining = -remaining
            iy, ix = divmod(flat, self.grid.nx)
            centres[node] = (0.5 * (self.grid.x_edges[ix] + self.grid.x_edges[ix + 1]),
                             0.5 * (self.grid.y_edges[iy] + self.grid.y_edges[iy + 1]))
            heapq.heappush(heap, (-(remaining - float(area[node])), flat))
        return project_centres(self.db, centres)

    def _nearest_slack_bin(self, ix: int, iy: int, target: np.ndarray, area: float,
                           remaining: np.ndarray) -> tuple[int, int] | None:
        """Find a nearby centre-bin with capacity, preferring net-centroid proximity."""
        for radius in range(1, self.config.anchor_max_radius + 1):
            x0, x1 = max(0, ix - radius), min(self.grid.nx - 1, ix + radius)
            y0, y1 = max(0, iy - radius), min(self.grid.ny - 1, iy + radius)
            candidates: list[tuple[int, int]] = []
            for tx in range(x0, x1 + 1):
                for ty in (y0, y1):
                    if remaining[ty, tx] >= area:
                        candidates.append((tx, ty))
            for ty in range(y0 + 1, y1):
                for tx in (x0, x1):
                    if remaining[ty, tx] >= area:
                        candidates.append((tx, ty))
            if candidates:
                return min(candidates, key=lambda p: (
                    (0.5 * (self.grid.x_edges[p[0]] + self.grid.x_edges[p[0] + 1]) - target[0]) ** 2
                    + (0.5 * (self.grid.y_edges[p[1]] + self.grid.y_edges[p[1] + 1]) - target[1]) ** 2))
        return None

    def _bulk_density_candidate(self, centres: np.ndarray) -> tuple[np.ndarray, int]:
        """Independent density-first construction using centre-bin capacity flow.

        This phase intentionally has no HPWL acceptance condition.  It creates a
        feasible anchor first; HPWL is handled only after the anchor exists.
        """
        result = centres.copy(); area = self.db.width * self.db.height
        ix, iy = self._centre_bins(result); flat = iy * self.grid.nx + ix
        capacity = self.grid.rho_target * self.grid.available_area
        load = np.bincount(flat[self.db.movable], weights=area[self.db.movable],
                           minlength=self.grid.nx * self.grid.ny).reshape(self.grid.ny, self.grid.nx)
        remaining = capacity - load
        overflow = np.maximum(-remaining, 0.0)
        sources = np.flatnonzero(overflow.ravel() > 0.0)
        if not len(sources):
            return result, 0
        targets = self._incident_targets(result)
        moved = 0
        for source in sources[np.argsort(overflow.ravel()[sources])[::-1]]:
            if moved >= self.config.anchor_max_moves_per_round:
                break
            sy, sx = divmod(int(source), self.grid.nx)
            members = np.flatnonzero(self.db.movable & (flat == source))
            members = members[np.argsort(self._node_degree[members], kind="stable")]
            for node in members:
                if remaining[sy, sx] >= 0.0 or moved >= self.config.anchor_max_moves_per_round:
                    break
                destination = self._nearest_slack_bin(sx, sy, targets[node], float(area[node]), remaining)
                if destination is None:
                    continue
                tx, ty = destination
                result[node] = (0.5 * (self.grid.x_edges[tx] + self.grid.x_edges[tx + 1]),
                                0.5 * (self.grid.y_edges[ty] + self.grid.y_edges[ty + 1]))
                remaining[sy, sx] += area[node]; remaining[ty, tx] -= area[node]
                moved += 1
        return project_centres(self.db, result), moved

    def _repair_candidate(self, centres: np.ndarray, result: OracleResult) -> np.ndarray | None:
        """One local capacity move selected by HPWL-damage-aware geometry."""
        ix, iy = self._centre_bins(centres)
        flat = iy * self.grid.nx + ix
        excess = np.maximum(result.occupancy - self.grid.rho_target * self.grid.available_area, 0.0)
        sources = np.flatnonzero(excess.ravel() > 0.0)
        slack = self.grid.rho_target * self.grid.available_area - result.occupancy
        targets = np.flatnonzero((slack > 0.0).ravel())
        if not len(sources) or not len(targets):
            return None
        targets = targets[np.argsort(slack.ravel()[targets])[::-1]][:self.config.repair_target_bins]
        source = int(sources[np.argmax(excess.ravel()[sources])])
        nodes = np.flatnonzero(self.db.movable & (flat == source))
        if not len(nodes):
            return None
        area = self.db.width * self.db.height
        # Low-degree / small displacement cells are cheaper to detach.
        nodes = nodes[np.argsort(self._node_degree[nodes] / np.maximum(area[nodes], 1.0), kind="stable")]
        targets_xy = self._incident_targets(centres)
        best: tuple[float, int, int] | None = None
        for node in nodes[:self.config.repair_candidates]:
            for target in targets:
                ty, tx = divmod(int(target), self.grid.nx)
                if slack[ty, tx] + 1e-9 < min(area[node], self.grid.bin_width * self.grid.bin_height):
                    continue
                cx = 0.5 * (self.grid.x_edges[tx] + self.grid.x_edges[tx + 1])
                cy = 0.5 * (self.grid.y_edges[ty] + self.grid.y_edges[ty + 1])
                damage = self._node_degree[node] * float(np.hypot(cx - targets_xy[node, 0], cy - targets_xy[node, 1]))
                score = damage / max(min(area[node], excess.ravel()[source]), 1.0)
                if best is None or score < best[0]:
                    best = (score, int(node), int(target))
        if best is None:
            return None
        _, node, target = best
        ty, tx = divmod(target, self.grid.nx)
        candidate = centres.copy()
        candidate[node] = (0.5 * (self.grid.x_edges[tx] + self.grid.x_edges[tx + 1]),
                           0.5 * (self.grid.y_edges[ty] + self.grid.y_edges[ty + 1]))
        return project_centres(self.db, candidate)

    def _exchange_candidate(self, centres: np.ndarray, result: OracleResult) -> np.ndarray | None:
        """Try one bounded area-similar swap when a direct move is unavailable."""
        ix, iy = self._centre_bins(centres)
        flat = iy * self.grid.nx + ix
        excess = np.maximum(result.occupancy - self.grid.rho_target * self.grid.available_area, 0.0)
        source_bins = np.flatnonzero(excess.ravel() > 0.0)
        if not len(source_bins):
            return None
        source_bin = int(source_bins[np.argmax(excess.ravel()[source_bins])])
        sources = np.flatnonzero(self.db.movable & (flat == source_bin))
        if not len(sources):
            return None
        slack = self.grid.rho_target * self.grid.available_area - result.occupancy
        target_bins = np.flatnonzero(slack.ravel() > 0.0)
        targets: list[int] = []
        for bin_id in target_bins[np.argsort(slack.ravel()[target_bins])[::-1]][:self.config.repair_target_bins]:
            targets.extend(np.flatnonzero(self.db.movable & (flat == bin_id))[:2].tolist())
        if not targets:
            return None
        area = self.db.width * self.db.height
        source = int(sources[np.argmin(self._node_degree[sources] / np.maximum(area[sources], 1.0))])
        target = min(targets, key=lambda node: abs(area[node] - area[source]))
        candidate = centres.copy()
        candidate[source], candidate[target] = centres[target], centres[source]
        return project_centres(self.db, candidate)

    def _build_anchor(self, centres: np.ndarray, result: OracleResult, target_budget: float) -> tuple[np.ndarray, OracleResult, int]:
        """Reduce exact overflow without allowing destructive HPWL expansion."""
        repairs = 0
        best_centres, best = centres.copy(), result
        # Stage A: independently reduce overflow with large capacity-aware moves.
        for _ in range(self.config.anchor_evacuation_rounds):
            if best.density_linear <= target_budget + 1e-9:
                break
            candidate, moved = self._bulk_density_candidate(best_centres)
            if moved == 0:
                break
            evaluated = self.oracle.evaluate(candidate)
            if evaluated.density_linear < best.density_linear - 1e-9:
                best_centres, best = candidate, evaluated
                repairs += moved
            else:
                break
        # Stage B: small local repairs can improve the exact overlap metric when
        # centre-bin assignment is already near-feasible.
        cap = np.inf if self.config.anchor_hpwl_ratio is None else result.hpwl * self.config.anchor_hpwl_ratio
        for _ in range(self.config.anchor_repair_rounds):
            if best.density_linear <= target_budget + 1e-9:
                break
            candidate = self._repair_candidate(best_centres, best)
            if candidate is None:
                candidate = self._exchange_candidate(best_centres, best)
            if candidate is None:
                break
            evaluated = self.oracle.evaluate(candidate)
            if evaluated.density_linear < best.density_linear - 1e-9 and evaluated.hpwl <= cap:
                best_centres, best = candidate, evaluated
                repairs += 1
            else:
                break
        return best_centres, best, repairs

    def _direction(self, result: OracleResult, centres: np.ndarray, step: float,
                   trust: np.ndarray, *, gradient: np.ndarray | None = None,
                   constrain_density: bool = True) -> np.ndarray:
        grad = result.gradient if gradient is None else gradient
        direction = -grad / self.preconditioner[:, None]
        if constrain_density:
            _, _, density_grad, _, _ = density_value_gradient(
                centres, self.db.width, self.db.height, self.db.movable,
                self.grid.x_edges, self.grid.y_edges, self.grid.available_area, self.grid.rho_target,
            )
            ddir = density_grad / self.preconditioner[:, None]
            directional = float(np.sum(density_grad * direction))
            norm_sq = float(np.sum(density_grad * ddir))
            if directional > 0.0 and norm_sq > 1e-12:
                direction -= directional / norm_sq * ddir
        displacement = step * direction
        np.clip(displacement, -trust[:, None], trust[:, None], out=displacement)
        displacement[self.db.fixed] = 0.0
        return project_centres(self.db, centres + displacement)

    def _update_trust(self, trust: np.ndarray, displacement: np.ndarray, accepted: bool) -> None:
        moved = np.any(np.abs(displacement) > 1e-12, axis=1) & self.db.movable
        lower = self.config.trust_min_rows * self.db.row_height
        upper = self.config.trust_max_rows * self.db.row_height
        factor = self.config.trust_expand if accepted else self.config.trust_shrink
        trust[moved] = np.clip(trust[moved] * factor, lower, upper)

    def _bundle_gradient(self, cuts: list[np.ndarray]) -> np.ndarray:
        # Equal-weight aggregate is a deliberately bounded-memory bundle model.
        return np.mean(np.stack(cuts, axis=0), axis=0)

    def _pressure_escape(self, centres: np.ndarray, result: OracleResult, budget: float, square_budget: float,
                         trust: np.ndarray, step: float) -> tuple[np.ndarray, OracleResult, bool]:
        base_hpwl = result.hpwl
        for points in self.config.pressure_percent_points:
            relaxed = budget + points / 100.0 * self.available_area
            probe_c, probe_r = centres.copy(), result
            for _ in range(self.config.pressure_steps):
                candidate = self._direction(probe_r, probe_c, step, trust, constrain_density=False)
                evaluated = self.oracle.evaluate(candidate)
                if evaluated.hpwl < probe_r.hpwl and evaluated.density_linear <= relaxed:
                    probe_c, probe_r = candidate, evaluated
                else:
                    break
            repaired_c, repaired_r, _ = self._build_anchor(probe_c, probe_r, budget)
            square_ok = (not self.config.enforce_challenge_square_budget
                         or self.challenge_square_density(repaired_r) <= square_budget + 1e-9)
            if (repaired_r.density_linear <= budget + 1e-9
                    and square_ok
                    and repaired_r.hpwl < base_hpwl):
                return repaired_c, repaired_r, True
        return centres, result, False

    def solve(self, initial_centres: np.ndarray, observer: Observer | None = None) -> HardSolveResult:
        if initial_centres.shape != (self.db.node_count, 2):
            raise ValueError("initial_centres has an invalid shape")
        centres = project_centres(self.db, initial_centres)
        result = self.oracle.evaluate(centres)
        target = (self.config.target_overflow_percent / 100.0 * self.available_area
                  if self.config.target_overflow_percent is not None else result.density_linear)
        anchor_c, anchor_r, repairs = self._build_anchor(centres, result, target)
        if anchor_r.density_linear > target + 1e-9:
            raise RuntimeError(
                "Unable to construct an anchor within the requested overflow cap: "
                f"{self._overflow_percent(anchor_r):.4f}% > "
                f"{100.0 * target / max(self.available_area, 1e-12):.4f}%"
            )
        centres, result = anchor_c.copy(), anchor_r
        # When supplied, the requested percentage is the hard compression cap.
        # The safe anchor may be better (including exactly 0%), but must not
        # silently turn a <=6% requirement into a strict-zero requirement.
        budget = target if self.config.target_overflow_percent is not None else result.density_linear
        square_budget = self.challenge_square_density(result)
        trust = np.full(self.db.node_count, self.config.trust_initial_rows * self.db.row_height)
        trust[self.db.fixed] = 0.0
        cuts: list[np.ndarray] = []
        records: list[dict[str, object]] = []
        step, stalled, pressure_accepts = self.config.initial_step, 0, 0
        started = time.perf_counter()
        for iteration in range(1, self.config.max_iterations + 1):
            if time.perf_counter() - started >= self.config.max_seconds:
                break
            use_bundle = len(cuts) >= 2 and stalled >= self.config.stall_limit
            grad = self._bundle_gradient(cuts) if use_bundle else None
            candidate = self._direction(result, centres, step, trust, gradient=grad)
            evaluated = self.oracle.evaluate(candidate)
            square_ok = (not self.config.enforce_challenge_square_budget
                         or self.challenge_square_density(evaluated) <= square_budget + 1e-9)
            accepted = bool(evaluated.hpwl < result.hpwl - 1e-8
                            and evaluated.density_linear <= budget + 1e-9
                            and square_ok)
            displacement = candidate - centres
            if accepted:
                centres, result = candidate, evaluated
                if self.config.target_overflow_percent is None:
                    budget = min(budget, result.density_linear)
                if self.config.enforce_challenge_square_budget:
                    square_budget = min(square_budget, self.challenge_square_density(result))
                self._update_trust(trust, displacement, True)
                step = min(step * 1.1, self.config.initial_step * 4.0)
                stalled = 0
                event = "bundle_serious" if use_bundle else "compress_serious"
            else:
                self._update_trust(trust, displacement, False)
                cuts.append(evaluated.gradient.astype(np.float32, copy=True))
                if len(cuts) > self.config.bundle_size:
                    cuts[0] = ((cuts[0].astype(np.float64) + cuts[1].astype(np.float64)) / 2.0).astype(np.float32)
                    del cuts[1]
                step = max(step * 0.5, self.config.min_step)
                stalled += 1
                event = "bundle_null" if use_bundle else "compress_null"
                if stalled >= self.config.stall_limit * 2:
                    escaped_c, escaped_r, escaped = self._pressure_escape(centres, result, budget, square_budget, trust, step)
                    if escaped:
                        centres, result = escaped_c, escaped_r
                        if self.config.target_overflow_percent is None:
                            budget = min(budget, result.density_linear)
                        if self.config.enforce_challenge_square_budget:
                            square_budget = min(square_budget, self.challenge_square_density(result))
                        pressure_accepts += 1
                        stalled, step, event = 0, self.config.initial_step, "pressure_repaired"
                    else:
                        stalled, event = 0, "pressure_discard"
            record = self._record(iteration, event, result, budget, square_budget, step, trust,
                                  serious=accepted, bundle_size=len(cuts), pressure=event.startswith("pressure"))
            records.append(record)
            if observer:
                observer(record, centres, result, trust)
        return HardSolveResult(centres, result, anchor_c, anchor_r, len(records), time.perf_counter() - started,
                               repairs, pressure_accepts, records)
