"""Low-overhead, reproducible run artifacts shared by solvers and plotters."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
from typing import Any

import numpy as np
import psutil

from gpplacer.model.config import SolverConfig
from gpplacer.model.types import PlacementDB
from gpplacer.solver.state import PlacementState


class RunLogger:
    """Write metadata, scalar iteration records, and selected placement snapshots."""

    _FIELDS = (
        "iteration", "elapsed_seconds", "phase", "hpwl", "density_linear",
        "overflow_ratio", "overflow_percent",
        "density_surrogate_squared", "objective", "step", "proximal_scale",
        "active_nets_change", "active_bins_change", "active_bins", "bundle_size",
        "backtracks", "serious_step", "predicted_decrease", "recoveries", "rss_bytes",
        "iteration_seconds", "oracle_seconds", "oracle_calls",
        "direction_cosine",
    )

    def __init__(self, root: str | Path, db: PlacementDB, config: SolverConfig) -> None:
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        case = Path(db.source_aux).stem
        self.run_dir = Path(root) / case / stamp
        self.snapshot_dir = self.run_dir / "snapshots"
        self.snapshot_dir.mkdir(parents=True, exist_ok=False)
        self.db = db
        self.config = config
        self._csv = (self.run_dir / "iterations.csv").open("w", newline="", encoding="utf-8")
        self._writer = csv.DictWriter(self._csv, fieldnames=self._FIELDS)
        self._writer.writeheader()
        self._write_metadata({"status": "running"})

    def _write_metadata(self, extra: dict[str, Any]) -> None:
        payload: dict[str, Any] = {
            "source_aux": self.db.source_aux,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "node_count": self.db.node_count,
            "fixed_count": int(self.db.fixed.sum()),
            "net_count": self.db.net_count,
            "pin_count": self.db.pin_count,
            "row_count": len(self.db.rows),
            "python": platform.python_version(),
            "config": self.config.to_dict(),
        }
        payload.update(extra)
        (self.run_dir / "metadata.json").write_text(
            json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8",
        )

    def save_snapshot(self, name: str, state: PlacementState) -> None:
        """Persist enough state to regenerate spatial figures without rerunning."""
        np.savez_compressed(
            self.snapshot_dir / f"{name}.npz", centres=state.centres,
            occupancy=state.oracle.occupancy, active_bins=state.oracle.active_bins,
        )

    def observe(self, record: dict[str, object], state: PlacementState) -> None:
        """Solver callback: append scalars and sample state at configured cadence."""
        row = dict(record)
        row["rss_bytes"] = psutil.Process().memory_info().rss
        self._writer.writerow({field: row.get(field, "") for field in self._FIELDS})
        self._csv.flush()
        iteration = int(record["iteration"])
        if iteration == 1 or iteration % self.config.snapshot_every == 0:
            self.save_snapshot(f"iteration_{iteration:05d}", state)

    def finish(self, *, status: str, summary: dict[str, Any], state: PlacementState | None = None) -> None:
        """Close the scalar log and atomically leave final result metadata."""
        if state is not None:
            self.save_snapshot("best", state)
        self._csv.close()
        self._write_metadata({"status": status, **summary})
