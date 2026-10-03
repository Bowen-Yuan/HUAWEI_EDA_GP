"""Persistent run artifacts for reproducibility, resume, and offline plots."""

from __future__ import annotations

import csv
from dataclasses import fields
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any, Callable

import numpy as np

from gpplacer.io.placement import write_placement
from gpplacer.model.types import PlacementDB

from adaptive_pareto.config import AdaptiveConfig
from adaptive_pareto.evaluator import AdaptiveEvaluator
from adaptive_pareto.archive import ArchivePoint
from adaptive_pareto.solver import AdaptiveSolveResult, SolverCheckpoint


def make_run_dir(root: str | Path, aux: str | Path) -> Path:
    candidate = Path(root) / Path(aux).stem / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    candidate.mkdir(parents=True, exist_ok=False)
    (candidate / "solutions").mkdir()
    (candidate / "snapshots").mkdir()
    (candidate / "figures").mkdir()
    (candidate / "checkpoints").mkdir()
    return candidate


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(rows[0])
    for row in rows[1:]:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_table(path: str | Path, rows: list[dict[str, object]]) -> None:
    """Public deterministic CSV writer for candidate and hierarchy summaries."""
    _write_csv(Path(path), rows)


def initialize_run(root: Path, *, db: PlacementDB, aux: str | Path, config: AdaptiveConfig, input_placement: str | None) -> None:
    """Persist enough provenance before the first solver step to permit recovery."""
    path = root / "metadata.json"
    if path.exists():
        return
    path.write_text(json.dumps({
        "status": "running", "source_aux": str(Path(aux).resolve()),
        "created_utc": datetime.now(timezone.utc).isoformat(), "node_count": db.node_count,
        "fixed_count": int(db.fixed.sum()), "net_count": db.net_count, "pin_count": db.pin_count,
        "config": config.to_dict(), "input_placement": input_placement,
        "evaluation_protocol": {
            "legacy": "positive-capacity linear overflow",
            "strict": "legacy plus movable occupancy in zero-capacity bins",
        },
    }, ensure_ascii=False, indent=2), encoding="utf-8")


def write_checkpoint(root: Path, checkpoint: SolverCheckpoint) -> None:
    """Atomically replace the latest standalone-solver recovery checkpoint."""
    directory = root / "checkpoints"
    directory.mkdir(exist_ok=True)
    archive_centres = (np.stack([point.centres for point in checkpoint.archive_points])
                       if checkpoint.archive_points else np.empty((0, *checkpoint.centres.shape), dtype=np.float64))
    np.savez_compressed(directory / "latest.npz", centres=checkpoint.centres, trust=checkpoint.trust,
                        archive_centres=archive_centres)
    archive = [{
        "hpwl": point.hpwl, "strict_overflow_percent": point.strict_overflow_percent,
        "legacy_overflow_percent": point.legacy_overflow_percent, "phase": point.phase,
        "iteration": point.iteration,
    } for point in checkpoint.archive_points]
    state = {
        field.name: getattr(checkpoint, field.name) for field in fields(checkpoint)
        if field.name not in {"centres", "trust", "archive_points", "records", "evacuation_moves"}
    }
    state.update({"archive": archive, "records": checkpoint.records, "evacuation_moves": checkpoint.evacuation_moves})
    (directory / "latest.json").write_text(json.dumps(state, ensure_ascii=False), encoding="utf-8")


def load_checkpoint(root: str | Path) -> SolverCheckpoint:
    """Load a checkpoint created by :func:`write_checkpoint`."""
    directory = Path(root) / "checkpoints"
    arrays = np.load(directory / "latest.npz")
    state = json.loads((directory / "latest.json").read_text(encoding="utf-8"))
    centres = arrays["centres"]
    archive_points = [
        ArchivePoint(arrays["archive_centres"][index], item["hpwl"], item["strict_overflow_percent"],
                     item["legacy_overflow_percent"], item["phase"], item["iteration"])
        for index, item in enumerate(state.pop("archive"))
    ]
    return SolverCheckpoint(centres=centres, trust=arrays["trust"], archive_points=archive_points,
                            records=state.pop("records"), evacuation_moves=state.pop("evacuation_moves"), **state)


def make_checkpoint_callback(
    root: Path, evaluator: AdaptiveEvaluator, config: AdaptiveConfig,
) -> Callable[[SolverCheckpoint], None]:
    """Checkpoint every callback and retain bounded accepted/time-spaced layouts."""
    snapshot_count = len(list((root / "snapshots").glob("progress_*.npz")))
    last_accepted = 0
    last_time = 0.0

    def callback(state: SolverCheckpoint) -> None:
        nonlocal snapshot_count, last_accepted, last_time
        write_checkpoint(root, state)
        if snapshot_count >= config.visualization.max_periodic_snapshots or not state.records:
            return
        accepted = state.accepted
        due_accepts = accepted - last_accepted >= config.visualization.snapshot_every_accepted
        due_time = state.elapsed_seconds - last_time >= config.visualization.snapshot_every_seconds
        if not (due_accepts or due_time):
            return
        _snapshot(root / "snapshots" / f"progress_{state.iteration:06d}.npz", state.centres, evaluator)
        snapshot_count += 1
        last_accepted, last_time = accepted, state.elapsed_seconds

    return callback


def _snapshot(path: Path, centres: np.ndarray, evaluator: AdaptiveEvaluator) -> None:
    result = evaluator.evaluate(centres)
    np.savez_compressed(
        path,
        centres=centres,
        occupancy=result.strict.occupancy,
        legacy_active_bins=result.legacy.active_bins,
        strict_active_bins=result.strict.active_bins,
        hpwl=result.hpwl,
        legacy_density_linear=result.legacy.linear,
        strict_density_linear=result.strict.linear,
        zero_capacity_occupancy=result.strict.zero_capacity_occupancy,
    )


def write_run(
    root: Path, *, db: PlacementDB, aux: str | Path, initial: np.ndarray,
    evaluator: AdaptiveEvaluator, config: AdaptiveConfig, result: AdaptiveSolveResult,
    extra_metadata: dict[str, Any] | None = None,
) -> None:
    """Write all solver-owned artifacts before plotting them offline."""
    write_placement(root / "solutions" / "solution_final.pl", db, result.centres)
    if result.archive.points:
        target = result.archive.target_representative(
            config.evaluation.target_overflow_low, config.evaluation.target_overflow_high,
        )
        write_placement(root / "solutions" / "solution_best_target.pl", db, target.centres)
        best_hpwl = min(result.archive.points, key=lambda point: point.hpwl)
        best_density = min(result.archive.points, key=lambda point: point.strict_overflow_percent)
        write_placement(root / "solutions" / "solution_min_hpwl.pl", db, best_hpwl.centres)
        write_placement(root / "solutions" / "solution_min_strict_overflow.pl", db, best_density.centres)
    _write_csv(root / "iterations.csv", result.records)
    _write_csv(root / "evacuation_moves.csv", result.evacuation_moves)
    archive_rows = [
        {
            "hpwl": point.hpwl,
            "strict_overflow_percent": point.strict_overflow_percent,
            "legacy_overflow_percent": point.legacy_overflow_percent,
            "phase": point.phase,
            "iteration": point.iteration,
        }
        for point in result.archive.points
    ]
    _write_csv(root / "pareto_archive.csv", archive_rows)
    input_snapshot = root / "snapshots" / "input.npz"
    if not input_snapshot.exists():
        _snapshot(input_snapshot, initial, evaluator)
    _snapshot(root / "snapshots" / "final.npz", result.centres, evaluator)
    if result.archive.points:
        _snapshot(root / "snapshots" / "best_target.npz", target.centres, evaluator)
    metadata: dict[str, Any] = {
        "status": result.status,
        "source_aux": str(Path(aux).resolve()),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "node_count": db.node_count,
        "fixed_count": int(db.fixed.sum()),
        "net_count": db.net_count,
        "pin_count": db.pin_count,
        "config": config.to_dict(),
        "evaluation_protocol": {
            "legacy": "positive-capacity linear overflow",
            "strict": "legacy plus movable occupancy in zero-capacity bins",
        },
        "elapsed_seconds": result.elapsed_seconds,
        "records": len(result.records),
        "archive_points": len(result.archive.points),
        "final_hpwl": result.oracle.hpwl,
        "final_legacy_overflow_percent": evaluator.overflow_percent(result.oracle, strict=False),
        "final_strict_overflow_percent": evaluator.overflow_percent(result.oracle, strict=True),
        "final_zero_capacity_occupancy": result.oracle.strict.zero_capacity_occupancy,
    }
    if extra_metadata:
        metadata.update(extra_metadata)
    (root / "metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
