"""Recursive hierarchy, diverse coarse seeds, and deterministic successive halving."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
import time

import numpy as np

from gpplacer.model.types import PlacementDB
from gpplacer.multilevel.coarsen import Coarsening, coarsen_once, uncoarsen_centres
from gpplacer.multilevel.initialization import (
    capacity_aware_seed, cluster_capacity_seed, connectivity_aware_seed,
)
from gpplacer.multilevel.partition_seed import recursive_partition_seed
from gpplacer.solver.steps import project_centres

from adaptive_pareto.config import AdaptiveConfig
from adaptive_pareto.evaluator import AdaptiveEvaluator
from adaptive_pareto.proxy import RidgeProxy
from adaptive_pareto.solver import AdaptiveParetoSolver, AdaptiveSolveResult


@dataclass(slots=True)
class Candidate:
    candidate_id: int
    generator: str
    seed: int
    centres: np.ndarray


@dataclass(slots=True)
class CandidateOutcome:
    candidate: Candidate
    result: AdaptiveSolveResult
    stage: str
    elapsed_seconds: float


@dataclass(slots=True)
class PortfolioResult:
    best: AdaptiveSolveResult
    candidate_rows: list[dict[str, object]]
    hierarchy_rows: list[dict[str, object]]
    elapsed_seconds: float


def build_hierarchy(db: PlacementDB, config: AdaptiveConfig) -> list[Coarsening]:
    """Build fine→coarse levels while preserving a reversible parent mapping."""
    levels: list[Coarsening] = []
    current = db
    for _ in range(config.hierarchy.max_levels):
        if current.node_count <= config.hierarchy.stop_nodes:
            break
        coarsening = coarsen_once(current, config.hierarchy.degree_limit)
        reduction = 1.0 - coarsening.coarse.node_count / current.node_count
        if reduction < config.hierarchy.min_reduction_ratio:
            break
        levels.append(coarsening)
        current = coarsening.coarse
    return levels


def _jitter_targets(targets: np.ndarray, db: PlacementDB, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    jitter = rng.normal(0.0, .35 * db.row_height, size=targets.shape)
    result = targets.copy()
    result[db.movable] += jitter[db.movable]
    return project_centres(db, result)


def _approximate_spectral_targets(db: PlacementDB, seed: int) -> np.ndarray:
    """Map two random hypergraph-Laplacian sketches to a deterministic 2-D seed.

    This avoids materialising a node-by-node Laplacian: random net signals are
    accumulated through the existing node-to-net CSR incidence, then centred
    and rank-mapped to the core.  It is deliberately an inexpensive *approximate*
    spectral initializer, not a replacement for a costly exact eigensolve.
    """
    rng = np.random.default_rng(seed)
    signals = rng.normal(size=(db.net_count, 2))
    summed = np.empty((db.node_count, 2), dtype=np.float64)
    for node in range(db.node_count):
        nets = db.node_nets[db.node_net_start[node]:db.node_net_start[node + 1]]
        summed[node] = signals[nets].mean(axis=0) if len(nets) else 0.0
    left, bottom, right, top = db.core_bounds
    targets = db.initial_centres.copy()
    movable = np.flatnonzero(db.movable)
    for axis, (lower, upper) in enumerate(((left, right), (bottom, top))):
        order = movable[np.argsort(summed[movable, axis], kind="stable")]
        targets[order, axis] = lower + (np.arange(len(order)) + .5) / max(len(order), 1) * (upper - lower)
    targets[db.fixed] = db.initial_centres[db.fixed]
    return targets


def _static_features(db: PlacementDB, evaluator: AdaptiveEvaluator, centres: np.ndarray) -> np.ndarray:
    result = evaluator.evaluate(centres)
    occupancy = result.strict.occupancy
    positive = evaluator.grid.available_area > 0.0
    utilisation = np.divide(occupancy, evaluator.grid.available_area, out=np.zeros_like(occupancy), where=positive)
    movable = np.flatnonzero(db.movable)
    span = np.ptp(centres[movable], axis=0) if len(movable) else np.zeros(2)
    core_width = evaluator.grid.x_edges[-1] - evaluator.grid.x_edges[0]
    core_height = evaluator.grid.y_edges[-1] - evaluator.grid.y_edges[0]
    return np.asarray([
        np.log1p(result.hpwl), evaluator.overflow_percent(result, strict=True),
        evaluator.overflow_percent(result, strict=False), float(np.var(utilisation[positive])),
        float(span[0] / max(core_width, 1.0)), float(span[1] / max(core_height, 1.0)),
    ])


def generate_candidates(db: PlacementDB, config: AdaptiveConfig) -> list[Candidate]:
    """Produce the configured portfolio on one coarse placement database."""
    evaluator = AdaptiveEvaluator(db, config.evaluation)
    candidates: list[Candidate] = []
    identifier = 0
    base_targets = connectivity_aware_seed(db, rounds=16)
    rng = np.random.default_rng(config.runtime.seed)
    for _ in range(config.portfolio.heavy_edge_candidates):
        seed = int(rng.integers(0, 2**31 - 1))
        targets = _jitter_targets(base_targets, db, seed)
        macro_bins = int(rng.choice(np.asarray([8, 12, 16, 20], dtype=np.int64)))
        centres = cluster_capacity_seed(db, evaluator.grid, targets, macro_bins)
        candidates.append(Candidate(identifier, "randomized_heavy_edge", seed, project_centres(db, centres)))
        identifier += 1
    for _ in range(config.portfolio.partition_candidates):
        seed = int(rng.integers(0, 2**31 - 1))
        targets = _jitter_targets(base_targets, db, seed)
        leaf_bins = int(rng.choice(np.asarray([4, 8, 16], dtype=np.int64)))
        centres = recursive_partition_seed(
            db, evaluator.grid, targets, leaf_bins=leaf_bins,
            refine_passes=2, degree_limit=config.hierarchy.degree_limit,
        )
        candidates.append(Candidate(identifier, "randomized_partition", seed, project_centres(db, centres)))
        identifier += 1
    for _ in range(config.portfolio.capacity_candidates):
        seed = int(rng.integers(0, 2**31 - 1))
        candidates.append(Candidate(identifier, "capacity_random", seed, capacity_aware_seed(db, evaluator.grid, seed)))
        identifier += 1
    for _ in range(config.portfolio.connection_candidates):
        seed = int(rng.integers(0, 2**31 - 1))
        candidates.append(Candidate(identifier, "connection", seed, _jitter_targets(base_targets, db, seed)))
        identifier += 1
    for _ in range(config.portfolio.spectral_candidates):
        seed = int(rng.integers(0, 2**31 - 1))
        targets = _approximate_spectral_targets(db, seed)
        macro_bins = int(rng.choice(np.asarray([8, 12, 16, 20], dtype=np.int64)))
        centres = cluster_capacity_seed(db, evaluator.grid, targets, macro_bins)
        candidates.append(Candidate(identifier, "approximate_spectral", seed, project_centres(db, centres)))
        identifier += 1
    return candidates


def _score_values(hpwl: float, strict_overflow: float, config: AdaptiveConfig) -> tuple[float, float]:
    low, high = config.evaluation.target_overflow_low, config.evaluation.target_overflow_high
    distance = max(low - strict_overflow, strict_overflow - high, 0.0)
    return distance, hpwl


def _select_diverse(
    outcomes: list[CandidateOutcome], keep: int, config: AdaptiveConfig,
    proxy_scores: dict[int, float] | None = None,
) -> list[CandidateOutcome]:
    if keep >= len(outcomes):
        return outcomes
    records: list[tuple[tuple[float, float], CandidateOutcome]] = []
    for outcome in outcomes:
        result = outcome.result.oracle
        # Results carry strict density but their evaluator is candidate-local;
        # each outcome has persisted the percentage in the row below.
        strict_percent = float(outcome.result.records[-1]["strict_overflow_percent"])
        direct = _score_values(result.hpwl, strict_percent, config)
        score = ((proxy_scores[outcome.candidate.candidate_id], direct[1])
                 if proxy_scores is not None else direct)
        records.append((score, outcome))
    records.sort(key=lambda item: item[0])
    selected: list[CandidateOutcome] = []
    by_generator: dict[str, int] = {}
    maximum = max(1, int(np.ceil(keep * config.portfolio.max_per_generator_fraction)))
    for _, outcome in records:
        generator = outcome.candidate.generator
        if by_generator.get(generator, 0) >= maximum:
            continue
        selected.append(outcome)
        by_generator[generator] = by_generator.get(generator, 0) + 1
        if len(selected) == keep:
            return selected
    selected_ids = {outcome.candidate.candidate_id for outcome in selected}
    for _, outcome in records:
        if outcome.candidate.candidate_id not in selected_ids:
            selected.append(outcome)
            selected_ids.add(outcome.candidate.candidate_id)
            if len(selected) == keep:
                break
    return selected


def _child_config(config: AdaptiveConfig, *, steps: int, seconds: float, seed: int) -> AdaptiveConfig:
    child = deepcopy(config)
    child.hierarchy.enabled = False
    child.portfolio.enabled = False
    child.joint.steps = steps
    child.runtime.max_steps = steps
    child.runtime.max_seconds = max(seconds, 1.0)
    child.runtime.seed = seed
    return child


def run_portfolio(db: PlacementDB, config: AdaptiveConfig) -> PortfolioResult:
    """Run coarse exploration, prune deterministically, then uncoarsen finalists."""
    started = time.perf_counter()
    levels = build_hierarchy(db, config) if config.hierarchy.enabled else []
    coarse_db = levels[-1].coarse if levels else db
    hierarchy_rows: list[dict[str, object]] = [{
        "level": 0, "nodes": db.node_count, "nets": db.net_count, "pins": db.pin_count,
    }]
    for index, level in enumerate(levels, start=1):
        hierarchy_rows.append({
            "level": index, "nodes": level.coarse.node_count, "nets": level.coarse.net_count,
            "pins": level.coarse.pin_count, "parent_nodes": level.fine.node_count,
        })
    candidates = generate_candidates(coarse_db, config)
    coarse_evaluator = AdaptiveEvaluator(coarse_db, config.evaluation)
    candidate_features = {candidate.candidate_id: _static_features(coarse_db, coarse_evaluator, candidate.centres) for candidate in candidates}
    outcomes: list[CandidateOutcome] = []
    candidate_rows: list[dict[str, object]] = []
    budgets = list(zip(config.portfolio.halving_steps, config.portfolio.halving_keep, strict=True))
    for stage_index, (steps, keep) in enumerate(budgets, start=1):
        stage_outcomes: list[CandidateOutcome] = []
        for candidate in candidates:
            remaining = config.runtime.max_seconds - (time.perf_counter() - started)
            if remaining <= 0.0:
                break
            per_candidate = remaining / max(len(candidates), 1)
            child = _child_config(config, steps=steps, seconds=per_candidate, seed=candidate.seed)
            solver = AdaptiveParetoSolver(coarse_db, child)
            before = time.perf_counter()
            result = solver.run(candidate.centres)
            outcome = CandidateOutcome(candidate, result, f"halving_{stage_index}", time.perf_counter() - before)
            stage_outcomes.append(outcome)
            strict_percent = float(result.records[-1]["strict_overflow_percent"]) if result.records else float("inf")
            candidate_rows.append({
                "candidate_id": candidate.candidate_id, "generator": candidate.generator, "seed": candidate.seed,
                "stage": outcome.stage, "status": result.status, "hpwl": result.oracle.hpwl,
                "strict_overflow_percent": strict_percent,
                "legacy_overflow_percent": float(result.records[-1]["legacy_overflow_percent"]) if result.records else float("inf"),
                "elapsed_seconds": outcome.elapsed_seconds, "survived": 0,
            })
        proxy_scores: dict[int, float] | None = None
        if config.portfolio.proxy_enabled and len(stage_outcomes) >= 3:
            features = np.vstack([candidate_features[outcome.candidate.candidate_id] for outcome in stage_outcomes])
            raw_targets = np.asarray([
                _score_values(outcome.result.oracle.hpwl, float(outcome.result.records[-1]["strict_overflow_percent"]), config)[0]
                + np.log1p(outcome.result.oracle.hpwl) * 1.0e-3
                for outcome in stage_outcomes
            ])
            predicted = RidgeProxy(config.portfolio.proxy_ridge).fit(features, raw_targets).predict(features)
            proxy_scores = {}
            for outcome, value in zip(stage_outcomes, predicted, strict=True):
                candidate_features[outcome.candidate.candidate_id] = np.append(candidate_features[outcome.candidate.candidate_id], value)
                proxy_scores[outcome.candidate.candidate_id] = float(value)
            for row, value in zip(candidate_rows[-len(stage_outcomes):], predicted, strict=True):
                row["proxy_score"] = float(value)
        outcomes = _select_diverse(stage_outcomes, keep, config, proxy_scores)
        survivor_ids = {outcome.candidate.candidate_id for outcome in outcomes}
        for row in candidate_rows:
            if row["stage"] == f"halving_{stage_index}" and row["candidate_id"] in survivor_ids:
                row["survived"] = 1
        candidates = [
            Candidate(outcome.candidate.candidate_id, outcome.candidate.generator, outcome.candidate.seed, outcome.result.centres)
            for outcome in outcomes
        ]
        if not candidates:
            break
    finalists: list[CandidateOutcome] = []
    for candidate in candidates:
        centres = candidate.centres
        current_db = coarse_db
        final_result: AdaptiveSolveResult | None = None
        # Expand from coarsest to fine.  Each expansion is followed by the
        # configured local refinement budget for that finer level.
        for reverse_index, coarsening in enumerate(reversed(levels)):
            centres = uncoarsen_centres(coarsening, centres, candidate.seed + reverse_index)
            current_db = coarsening.fine
            steps = config.hierarchy.level_refine_steps[min(reverse_index, len(config.hierarchy.level_refine_steps) - 1)]
            remaining = config.runtime.max_seconds - (time.perf_counter() - started)
            if remaining <= 0.0:
                break
            child = _child_config(config, steps=steps, seconds=remaining / max(len(candidates), 1), seed=candidate.seed)
            final_result = AdaptiveParetoSolver(current_db, child).run(centres)
            centres = final_result.centres
        if not levels:
            remaining = config.runtime.max_seconds - (time.perf_counter() - started)
            child = _child_config(config, steps=config.joint.steps, seconds=remaining, seed=candidate.seed)
            final_result = AdaptiveParetoSolver(db, child).run(centres)
        elif current_db is db:
            remaining = config.runtime.max_seconds - (time.perf_counter() - started)
            if remaining > 0.0:
                child = _child_config(config, steps=config.joint.steps, seconds=remaining / max(len(candidates), 1), seed=candidate.seed)
                final_result = AdaptiveParetoSolver(db, child).run(centres)
        if final_result is not None:
            finalist = CandidateOutcome(candidate, final_result, "final", final_result.elapsed_seconds)
            finalists.append(finalist)
            candidate_rows.append({
                "candidate_id": candidate.candidate_id, "generator": candidate.generator, "seed": candidate.seed,
                "stage": "final", "status": final_result.status, "hpwl": final_result.oracle.hpwl,
                "strict_overflow_percent": float(final_result.records[-1]["strict_overflow_percent"]) if final_result.records else float("inf"),
                "legacy_overflow_percent": float(final_result.records[-1]["legacy_overflow_percent"]) if final_result.records else float("inf"),
                "elapsed_seconds": final_result.elapsed_seconds, "survived": 1,
            })
    if not finalists:
        raise RuntimeError("No portfolio candidate completed before the deadline.")
    best = _select_diverse(finalists, 1, config)[0].result
    return PortfolioResult(
        best=best, candidate_rows=candidate_rows, hierarchy_rows=hierarchy_rows,
        elapsed_seconds=time.perf_counter() - started,
    )
