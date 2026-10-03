"""Independent long-run HPWL-only experiment for a complete placement.

Overflow is logged only as a diagnostic.  It has no role in candidate
generation, backtracking, trust-region updates, or acceptance.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
import time

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np

from gpplacer.io.bookshelf import load_bookshelf
from gpplacer.io.placement import read_placement_centres, write_placement
from gpplacer.oracle.grid import build_density_grid
from gpplacer.oracle.objective import ObjectiveOracle, OracleResult
from gpplacer.solver.steps import build_preconditioner, project_centres


@dataclass(slots=True)
class HPWLOnlyConfig:
    rho_target: float = .8
    bin_rows: int = 8
    iterations: int = 1_000
    max_seconds: float = 900.0
    initial_step: float = 4.0
    min_step: float = 1.0 / 256.0
    initial_trust_rows: float = 4.0
    max_trust_rows: float = 16.0
    output_root: str = "code-hard/runs_hpwl_only_from_25"


class HPWLOnlySolver:
    def __init__(self, db, grid, cfg: HPWLOnlyConfig) -> None:
        self.db, self.grid, self.cfg = db, grid, cfg
        self.oracle = ObjectiveOracle(db, grid, density_lambda=0.0)
        self.preconditioner = build_preconditioner(db)
        self.available = float(grid.available_area.sum())

    def overflow_percent(self, result: OracleResult) -> float:
        return 100.0 * result.density_linear / max(self.available, 1e-12)

    def run(self, initial: np.ndarray) -> tuple[np.ndarray, OracleResult, list[dict[str, object]], float]:
        centres = project_centres(self.db, initial)
        result = self.oracle.evaluate(centres)
        trust = np.full(self.db.node_count, self.cfg.initial_trust_rows * self.db.row_height)
        trust[self.db.fixed] = 0.0
        step = self.cfg.initial_step
        records: list[dict[str, object]] = []
        started = time.perf_counter()
        for iteration in range(1, self.cfg.iterations + 1):
            if time.perf_counter() - started >= self.cfg.max_seconds:
                break
            # This is deliberately only -M^{-1} grad(HPWL). No density
            # direction, tangent projection, overflow cap, or density repair.
            direction = -result.gradient / self.preconditioner[:, None]
            displacement = step * direction
            np.clip(displacement, -trust[:, None], trust[:, None], out=displacement)
            displacement[self.db.fixed] = 0.0
            candidate = project_centres(self.db, centres + displacement)
            trial = self.oracle.evaluate(candidate)
            accepted = trial.hpwl < result.hpwl - 1e-8
            moved = np.any(np.abs(displacement) > 1e-12, axis=1) & self.db.movable
            if accepted:
                centres, result = candidate, trial
                trust[moved] = np.minimum(trust[moved] * 1.08, self.cfg.max_trust_rows * self.db.row_height)
                step = min(step * 1.08, self.cfg.initial_step * 16.0)
                event = "hpwl_accept"
            else:
                trust[moved] = np.maximum(trust[moved] * .60, .25 * self.db.row_height)
                step *= .5
                event = "hpwl_reject"
                # Restart only the HPWL step scale after a prolonged local
                # plateau; no density information enters this decision.
                if step < self.cfg.min_step:
                    step = self.cfg.initial_step
                    event = "hpwl_restart"
            records.append({
                "iteration": iteration, "event": event, "hpwl": result.hpwl,
                "density_linear": result.density_linear,
                "overflow_percent": self.overflow_percent(result), "step": step,
                "trust_p50": float(np.quantile(trust[self.db.movable], .5)),
                "trust_p90": float(np.quantile(trust[self.db.movable], .9)),
                "accepted": int(accepted),
            })
        return centres, result, records, time.perf_counter() - started


def plot(root: Path, records: list[dict[str, object]]) -> None:
    figures = root / "figures"; figures.mkdir(exist_ok=True)
    x = np.array([float(row["iteration"]) for row in records])
    hpwl = np.array([float(row["hpwl"]) for row in records])
    overflow = np.array([float(row["overflow_percent"]) for row in records])
    accepted = np.array([int(row["accepted"]) for row in records], dtype=bool)
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.4), constrained_layout=True, dpi=180)
    axes[0].plot(x, hpwl, color="#0072B2", linewidth=1.3)
    axes[0].scatter(x[accepted], hpwl[accepted], color="#009E73", s=7, label="accepted HPWL step")
    axes[0].yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value / 1e9:.3f}"))
    axes[0].set(xlabel="iteration", ylabel="exact HPWL (×10⁹)", title="HPWL-only long-run convergence")
    axes[0].legend(fontsize=8)
    axes[1].plot(x, overflow, color="#D55E00", linewidth=1.3)
    axes[1].set(xlabel="iteration", ylabel="linear overflow (%)", title="Overflow diagnostic (not constrained)")
    fig.savefig(figures / "hpwl_only_convergence.png", dpi=220); plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Long-run HPWL-only experiment")
    parser.add_argument("--aux", required=True); parser.add_argument("--placement", required=True)
    parser.add_argument("--iterations", type=int); parser.add_argument("--max-seconds", type=float)
    parser.add_argument("--output-root"); args = parser.parse_args()
    cfg = HPWLOnlyConfig()
    if args.iterations is not None: cfg.iterations = args.iterations
    if args.max_seconds is not None: cfg.max_seconds = args.max_seconds
    if args.output_root is not None: cfg.output_root = args.output_root
    db = load_bookshelf(args.aux); grid = build_density_grid(db, rho_target=cfg.rho_target, bin_rows=cfg.bin_rows)
    initial = read_placement_centres(args.placement, db)
    solver = HPWLOnlySolver(db, grid, cfg)
    initial = project_centres(db, initial); initial_result = solver.oracle.evaluate(initial)
    centres, final, records, elapsed = solver.run(initial)
    root = Path(cfg.output_root) / Path(args.aux).stem / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    root.mkdir(parents=True, exist_ok=False)
    write_placement(root / "solution.pl", db, centres)
    np.savez_compressed(root / "snapshots.npz", initial_centres=initial, final_centres=centres,
                        initial_occupancy=initial_result.occupancy, final_occupancy=final.occupancy)
    with (root / "iterations.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0])); writer.writeheader(); writer.writerows(records)
    metadata = {"status": "completed", "source_aux": str(Path(args.aux).resolve()),
                "input_placement": str(Path(args.placement).resolve()), "config": asdict(cfg),
                "initial_hpwl": initial_result.hpwl, "initial_overflow_percent": solver.overflow_percent(initial_result),
                "final_hpwl": final.hpwl, "final_overflow_percent": solver.overflow_percent(final),
                "iterations": len(records), "elapsed_seconds": elapsed,
                "overflow_used_for_acceptance": False}
    (root / "metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    plot(root, records)
    print(json.dumps({"run_dir": str(root), "hpwl": final.hpwl,
                      "overflow_percent_diagnostic": solver.overflow_percent(final)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
