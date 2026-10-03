"""Guarded rescue and adaptive Pareto joint descent on the strict evaluator."""

from __future__ import annotations

from dataclasses import dataclass
import time
from typing import Callable

import numpy as np

from gpplacer.model.types import PlacementDB
from gpplacer.oracle.evaluation import DualOracleResult
from gpplacer.solver.steps import build_preconditioner, project_centres

from adaptive_pareto.archive import ArchivePoint, StratifiedParetoArchive
from adaptive_pareto.config import AdaptiveConfig
from adaptive_pareto.directions import density_direction
from adaptive_pareto.evacuation import local_evacuation_candidate
from adaptive_pareto.evaluator import AdaptiveEvaluator
from adaptive_pareto.rescue import SolverPhase, choose_phase, dual_priority


@dataclass(slots=True)
class AdaptiveSolveResult:
    centres: np.ndarray
    oracle: DualOracleResult
    archive: StratifiedParetoArchive
    records: list[dict[str, object]]
    evacuation_moves: list[dict[str, float | int | str]]
    elapsed_seconds: float
    status: str


@dataclass(slots=True)
class SolverCheckpoint:
    """Complete recoverable state at an accepted/rejected controller boundary."""

    centres: np.ndarray
    trust: np.ndarray
    archive_points: list[ArchivePoint]
    records: list[dict[str, object]]
    evacuation_moves: list[dict[str, float | int | str]]
    iteration: int
    accepted: int
    lambda_value: float | None
    hpwl_norm_ema: float | None
    density_norm_ema: float | None
    phase: str
    joint_stable: int
    stalled_rescue_blocks: int
    rescue_step: int
    joint_steps: int
    rejected_density: int
    entry_hpwl: float
    block_hpwl: float
    block_overflow: float
    step: float
    elapsed_seconds: float


class AdaptiveParetoSolver:
    """Single-seed controller; hierarchy and portfolio feed it candidate seeds."""

    def __init__(self, db: PlacementDB, config: AdaptiveConfig) -> None:
        config.validate()
        self.db = db
        self.config = config
        self.evaluator = AdaptiveEvaluator(db, config.evaluation)
        self.preconditioner = build_preconditioner(db)
        self.archive = StratifiedParetoArchive(
            per_band=config.joint.archive_per_band, maximum=config.joint.archive_maximum,
        )
        self.records: list[dict[str, object]] = []
        self.evacuation_moves: list[dict[str, float | int | str]] = []
        self.iteration = 0
        self.accepted = 0
        self._lambda: float | None = None
        self._hpwl_norm_ema: float | None = None
        self._density_norm_ema: float | None = None

    def _metrics(self, result: DualOracleResult) -> tuple[float, float, float]:
        return (
            self.evaluator.normalized_hpwl(result),
            self.evaluator.overflow_percent(result, strict=True),
            self.evaluator.overflow_percent(result, strict=False),
        )

    def _candidate(self, centres: np.ndarray, direction: np.ndarray, step: float, trust: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        displacement = step * direction
        np.clip(displacement, -trust[:, None], trust[:, None], out=displacement)
        displacement[self.db.fixed] = 0.0
        return project_centres(self.db, centres + displacement), displacement

    def _adaptive_lambda(self, result: DualOracleResult) -> float:
        movable = self.db.movable
        h_norm = float(np.linalg.norm(result.hpwl_gradient[movable]))
        d_norm = float(np.linalg.norm(result.strict.gradient[movable]))
        beta = self.config.joint.lambda_ema_beta
        self._hpwl_norm_ema = h_norm if self._hpwl_norm_ema is None else beta * self._hpwl_norm_ema + (1.0 - beta) * h_norm
        self._density_norm_ema = d_norm if self._density_norm_ema is None else beta * self._density_norm_ema + (1.0 - beta) * d_norm
        base = self._hpwl_norm_ema / max(self._density_norm_ema, self._hpwl_norm_ema * 1.0e-12, 1.0e-12)
        _, strict_overflow, _ = self._metrics(result)
        target = .5 * (self.config.evaluation.target_overflow_low + self.config.evaluation.target_overflow_high)
        band = max(.5 * (self.config.evaluation.target_overflow_high - self.config.evaluation.target_overflow_low), 1.0e-9)
        pressure = (strict_overflow - target) / band
        value = base * float(np.exp(self.config.joint.lambda_gain * np.clip(pressure, -4.0, 4.0)))
        value = float(np.clip(value, self.config.joint.lambda_min, self.config.joint.lambda_max))
        if self._lambda is not None:
            value = float(np.clip(
                value,
                self._lambda * self.config.joint.lambda_change_min,
                self._lambda * self.config.joint.lambda_change_max,
            ))
        self._lambda = value
        return value

    def _record(
        self, *, phase: SolverPhase, event: str, result: DualOracleResult, step: float,
        trust: np.ndarray, lambda_value: float, accepted: bool, centres: np.ndarray,
        single_axis_fraction: float = 0.0, elapsed_seconds: float = 0.0,
    ) -> None:
        self.iteration += 1
        normalized_hpwl, strict_overflow, legacy_overflow = self._metrics(result)
        archived = False
        if accepted:
            archived = self.archive.consider(
                centres, result, strict_overflow_percent=strict_overflow,
                legacy_overflow_percent=legacy_overflow, phase=phase.value, iteration=self.iteration,
            )
        movable_trust = trust[self.db.movable]
        self.records.append({
            "iteration": self.iteration,
            "phase": phase.value,
            "event": event,
            "hpwl": result.hpwl,
            "normalized_hpwl": normalized_hpwl,
            "legacy_density_linear": result.legacy.linear,
            "legacy_overflow_percent": legacy_overflow,
            "strict_density_linear": result.strict.linear,
            "strict_overflow_percent": strict_overflow,
            "zero_capacity_occupancy": result.strict.zero_capacity_occupancy,
            "strict_active_bins": result.strict.active_bin_count,
            "legacy_active_bins": result.legacy.active_bin_count,
            "lambda": lambda_value,
            "step": step,
            "trust_p50": float(np.quantile(movable_trust, .5)),
            "trust_p90": float(np.quantile(movable_trust, .9)),
            "single_axis_fraction": single_axis_fraction,
            "accepted": int(accepted),
            "archived": int(archived),
            "elapsed_seconds": elapsed_seconds,
        })

    def _accept(self, *, phase: SolverPhase, current: DualOracleResult, trial: DualOracleResult,
                lambda_value: float, block_hpwl: float, block_overflow: float, entry_hpwl: float) -> bool:
        _, current_overflow, _ = self._metrics(current)
        _, trial_overflow, _ = self._metrics(trial)
        if phase is SolverPhase.HPWL_RESCUE:
            return (
                trial.hpwl < current.hpwl - 1.0e-8
                and trial_overflow <= self.config.rescue.overflow_enter_percent + 1.0e-9
                and trial_overflow <= block_overflow + self.config.rescue.hpwl_rescue_overflow_increase_percent + 1.0e-9
            )
        if phase is SolverPhase.DENSITY_RESCUE:
            return (
                trial.strict.linear < current.strict.linear - 1.0e-8
                and trial.hpwl <= block_hpwl * (1.0 + self.config.rescue.density_rescue_hpwl_increase_ratio) + 1.0e-8
                and trial.hpwl <= entry_hpwl * self.config.rescue.density_rescue_total_hpwl_ratio + 1.0e-8
            )
        current_merit = current.hpwl + lambda_value * current.strict.linear
        trial_merit = trial.hpwl + lambda_value * trial.strict.linear
        return trial_merit < current_merit - 1.0e-8

    def _direction(self, phase: SolverPhase, result: DualOracleResult, trust: np.ndarray) -> tuple[np.ndarray, float, float]:
        if phase is SolverPhase.HPWL_RESCUE:
            return -result.hpwl_gradient / self.preconditioner[:, None], 0.0, 0.0
        density_direction_value, single_axis = density_direction(
            result.strict.gradient, self.preconditioner, trust, self.db.movable,
            mode=self.config.direction.mode, axis_dominance_ratio=self.config.direction.axis_dominance_ratio,
        )
        fraction = float(np.mean(single_axis[self.db.movable]))
        if phase is SolverPhase.DENSITY_RESCUE:
            return density_direction_value, 0.0, fraction
        lambda_value = self._adaptive_lambda(result)
        hpwl_direction = -result.hpwl_gradient / self.preconditioner[:, None]
        return hpwl_direction + lambda_value * density_direction_value, lambda_value, fraction

    def _attempt_evacuation(
        self, centres: np.ndarray, result: DualOracleResult, *, phase: SolverPhase,
        block_hpwl: float, block_overflow: float, entry_hpwl: float,
    ) -> tuple[np.ndarray, DualOracleResult, bool]:
        candidate, moves = local_evacuation_candidate(
            self.db, self.evaluator.grid, centres, result,
            max_radius_bins=self.config.evacuation.max_radius_bins,
            max_moves=self.config.evacuation.max_batch_moves,
        )
        if candidate is None:
            return centres, result, False
        trial = self.evaluator.evaluate(candidate)
        accepted = self._accept(
            phase=SolverPhase.DENSITY_RESCUE, current=result, trial=trial, lambda_value=0.0,
            block_hpwl=block_hpwl, block_overflow=block_overflow, entry_hpwl=entry_hpwl,
        )
        if accepted:
            for move in moves:
                move["iteration"] = self.iteration + 1
                move["phase"] = phase.value
                move["accepted"] = 1
            self.evacuation_moves.extend(moves)
            return candidate, trial, True
        for move in moves:
            move["iteration"] = self.iteration + 1
            move["phase"] = phase.value
            move["accepted"] = 0
        self.evacuation_moves.extend(moves)
        return centres, result, False

    def _checkpoint(
        self, *, centres: np.ndarray, trust: np.ndarray, phase: SolverPhase, joint_stable: int,
        stalled_rescue_blocks: int, rescue_step: int, joint_steps: int, rejected_density: int,
        entry_hpwl: float, block_hpwl: float, block_overflow: float, step: float, elapsed_seconds: float,
    ) -> SolverCheckpoint:
        return SolverCheckpoint(
            centres=centres.copy(), trust=trust.copy(), archive_points=[
                ArchivePoint(point.centres.copy(), point.hpwl, point.strict_overflow_percent,
                             point.legacy_overflow_percent, point.phase, point.iteration)
                for point in self.archive.points
            ],
            records=list(self.records), evacuation_moves=list(self.evacuation_moves), iteration=self.iteration,
            accepted=self.accepted, lambda_value=self._lambda, hpwl_norm_ema=self._hpwl_norm_ema,
            density_norm_ema=self._density_norm_ema, phase=phase.value, joint_stable=joint_stable,
            stalled_rescue_blocks=stalled_rescue_blocks, rescue_step=rescue_step, joint_steps=joint_steps,
            rejected_density=rejected_density, entry_hpwl=entry_hpwl, block_hpwl=block_hpwl,
            block_overflow=block_overflow, step=step, elapsed_seconds=elapsed_seconds,
        )

    def run(
        self, initial: np.ndarray, *, resume: SolverCheckpoint | None = None,
        on_checkpoint: Callable[[SolverCheckpoint], None] | None = None,
    ) -> AdaptiveSolveResult:
        if resume is None:
            centres = project_centres(self.db, initial)
            result = self.evaluator.evaluate(centres)
            trust = np.full(self.db.node_count, self.config.joint.trust_rows * self.db.row_height)
            trust[self.db.fixed] = 0.0
            _, strict_overflow, legacy_overflow = self._metrics(result)
            self.archive.consider(
                centres, result, strict_overflow_percent=strict_overflow, legacy_overflow_percent=legacy_overflow,
                phase="input", iteration=0,
            )
            phase = SolverPhase.JOINT
            joint_stable = stalled_rescue_blocks = rescue_step = joint_steps = rejected_density = 0
            entry_hpwl = block_hpwl = result.hpwl
            block_overflow = strict_overflow
            step = self.config.joint.initial_step
            elapsed_before = 0.0
        else:
            centres = project_centres(self.db, resume.centres)
            result = self.evaluator.evaluate(centres)
            trust = resume.trust.copy()
            trust[self.db.fixed] = 0.0
            self.archive.points = resume.archive_points
            self.records = list(resume.records)
            self.evacuation_moves = list(resume.evacuation_moves)
            self.iteration = resume.iteration
            self.accepted = resume.accepted
            self._lambda = resume.lambda_value
            self._hpwl_norm_ema = resume.hpwl_norm_ema
            self._density_norm_ema = resume.density_norm_ema
            phase = SolverPhase(resume.phase)
            joint_stable, stalled_rescue_blocks = resume.joint_stable, resume.stalled_rescue_blocks
            rescue_step, joint_steps, rejected_density = resume.rescue_step, resume.joint_steps, resume.rejected_density
            entry_hpwl, block_hpwl, block_overflow, step = resume.entry_hpwl, resume.block_hpwl, resume.block_overflow, resume.step
            elapsed_before = resume.elapsed_seconds
        started = time.perf_counter()
        starting_iteration = self.iteration
        status = "completed"

        while elapsed_before + time.perf_counter() - started < self.config.runtime.max_seconds:
            if (self.config.runtime.max_steps is not None
                    and self.iteration - starting_iteration >= self.config.runtime.max_steps):
                status = "step_limit"
                break
            normalized_hpwl, strict_overflow, _ = self._metrics(result)
            emergency = choose_phase(normalized_hpwl, strict_overflow, self.config.rescue)
            if emergency is SolverPhase.DUAL_EMERGENCY:
                requested_phase = dual_priority(normalized_hpwl, strict_overflow, self.config.rescue)
            elif emergency is not SolverPhase.JOINT:
                requested_phase = emergency
            elif normalized_hpwl <= self.config.rescue.hpwl_exit_norm and strict_overflow <= self.config.rescue.overflow_exit_percent:
                joint_stable += 1
                requested_phase = SolverPhase.JOINT if joint_stable >= self.config.rescue.stable_exit_steps else phase
            else:
                joint_stable = 0
                requested_phase = SolverPhase.HPWL_RESCUE if normalized_hpwl > self.config.rescue.hpwl_exit_norm else SolverPhase.DENSITY_RESCUE

            if requested_phase is not phase:
                phase = requested_phase
                rescue_step = 0
                block_hpwl, block_overflow = result.hpwl, strict_overflow
                step = self.config.joint.initial_step
            if phase is not SolverPhase.JOINT:
                rescue_step += 1
                if rescue_step > self.config.rescue.block_steps:
                    if phase is SolverPhase.HPWL_RESCUE:
                        improved = result.hpwl < block_hpwl * (1.0 - 1.0e-4)
                    else:
                        improved = result.strict.linear < (block_overflow / 100.0 * self.evaluator.available_area) * (1.0 - 1.0e-4)
                    stalled_rescue_blocks = 0 if improved else stalled_rescue_blocks + 1
                    rescue_step = 1
                    block_hpwl, block_overflow = result.hpwl, strict_overflow
                    if stalled_rescue_blocks >= self.config.rescue.max_blocks:
                        status = "rescue_stalled"
                        break
            elif joint_steps >= self.config.joint.steps:
                break

            direction, lambda_value, single_axis_fraction = self._direction(phase, result, trust)
            candidate, displacement = self._candidate(centres, direction, step, trust)
            trial = self.evaluator.evaluate(candidate)
            accepted = self._accept(
                phase=phase, current=result, trial=trial, lambda_value=lambda_value,
                block_hpwl=block_hpwl, block_overflow=block_overflow, entry_hpwl=entry_hpwl,
            )
            event = f"{phase.value}_accept" if accepted else f"{phase.value}_reject"
            moved = np.any(np.abs(displacement) > 1.0e-12, axis=1) & self.db.movable
            if accepted:
                centres, result = candidate, trial
                trust[moved] = np.minimum(trust[moved] * 1.06, self.config.joint.max_trust_rows * self.db.row_height)
                step = min(step * 1.06, self.config.joint.initial_step * 16.0)
                self.accepted += 1
                rejected_density = 0
            else:
                trust[moved] = np.maximum(trust[moved] * .65, .25 * self.db.row_height)
                step *= .5
                if phase is SolverPhase.DENSITY_RESCUE:
                    rejected_density += 1
                if step < self.config.joint.min_step:
                    step = self.config.joint.initial_step
                    event = f"{phase.value}_restart"
            if phase is SolverPhase.DENSITY_RESCUE and rejected_density >= 5 and self.config.evacuation.enabled:
                evacuated_centres, evacuated_result, evacuated = self._attempt_evacuation(
                    centres, result, phase=phase, block_hpwl=block_hpwl,
                    block_overflow=block_overflow, entry_hpwl=entry_hpwl,
                )
                if evacuated:
                    centres, result = evacuated_centres, evacuated_result
                    event = "density_local_evacuation"
                    accepted = True
                    rejected_density = 0
            self._record(
                phase=phase, event=event, result=result, step=step, trust=trust,
                lambda_value=lambda_value, accepted=accepted, centres=centres,
                single_axis_fraction=single_axis_fraction,
                elapsed_seconds=elapsed_before + time.perf_counter() - started,
            )
            if on_checkpoint is not None and self.iteration % self.config.runtime.checkpoint_every_steps == 0:
                on_checkpoint(self._checkpoint(
                    centres=centres, trust=trust, phase=phase, joint_stable=joint_stable,
                    stalled_rescue_blocks=stalled_rescue_blocks, rescue_step=rescue_step, joint_steps=joint_steps,
                    rejected_density=rejected_density, entry_hpwl=entry_hpwl, block_hpwl=block_hpwl,
                    block_overflow=block_overflow, step=step,
                    elapsed_seconds=elapsed_before + time.perf_counter() - started,
                ))
            if phase is SolverPhase.JOINT:
                joint_steps += 1
        else:
            status = "timeout"
        elapsed = elapsed_before + time.perf_counter() - started
        if elapsed >= self.config.runtime.max_seconds and status == "completed":
            status = "timeout"
        if on_checkpoint is not None:
            on_checkpoint(self._checkpoint(
                centres=centres, trust=trust, phase=phase, joint_stable=joint_stable,
                stalled_rescue_blocks=stalled_rescue_blocks, rescue_step=rescue_step, joint_steps=joint_steps,
                rejected_density=rejected_density, entry_hpwl=entry_hpwl, block_hpwl=block_hpwl,
                block_overflow=block_overflow, step=step, elapsed_seconds=elapsed,
            ))
        return AdaptiveSolveResult(
            centres, result, self.archive, self.records, self.evacuation_moves, elapsed, status,
        )
