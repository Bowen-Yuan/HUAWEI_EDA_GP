"""Independent alternating overflow/HPWL placement prototype.

This module deliberately does not import :mod:`hardplacer`.  It starts from a
given low-HPWL/high-overflow placement and repeatedly runs three *separate*
acceptance policies:

    density-only repair -> HPWL-only compression under an overflow cap
    -> density-only repair while preserving the just-earned HPWL.

All values and gradients are exact (cell-bin overlap and HPWL); no smoothing
surrogate is used.
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
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np

from gpplacer.io.bookshelf import load_bookshelf
from gpplacer.io.placement import read_placement_centres, write_placement
from gpplacer.oracle.density import density_value_gradient
from gpplacer.oracle.grid import build_density_grid
from gpplacer.oracle.objective import ObjectiveOracle, OracleResult
from gpplacer.solver.steps import build_preconditioner, project_centres


@dataclass(slots=True)
class AlternatingConfig:
    rho_target: float = .8
    bin_rows: int = 8
    overflow_cap_percent: float = 6.0
    cycles: int = 3
    density_iterations: int = 320
    hpwl_iterations: int = 120
    post_density_iterations: int = 80
    hpwl_guard_percent: float = .25
    initial_step: float = 1.0
    trust_rows: float = 4.0
    max_seconds: float = 900.0
    output_root: str = "code-hard/runs_alternating"


class AlternatingSolver:
    def __init__(self, db, grid, config: AlternatingConfig) -> None:
        self.db, self.grid, self.cfg = db, grid, config
        self.oracle = ObjectiveOracle(db, grid, density_lambda=0.0)
        self.precond = build_preconditioner(db)
        self.available = float(grid.available_area.sum())
        self.node_degree = np.diff(db.node_net_start).astype(np.float64)
        self.trust = np.full(db.node_count, config.trust_rows * db.row_height)
        self.trust[db.fixed] = 0.0
        self.records: list[dict[str, object]] = []
        self.iteration = 0

    def overflow(self, result: OracleResult) -> float:
        return 100.0 * result.density_linear / max(self.available, 1e-12)

    def _record(self, phase: str, event: str, result: OracleResult) -> None:
        self.iteration += 1
        self.records.append({"iteration": self.iteration, "phase": phase, "event": event,
                             "hpwl": result.hpwl, "density_linear": result.density_linear,
                             "overflow_percent": self.overflow(result),
                             "overflow_cap_percent": self.cfg.overflow_cap_percent,
                             "active_bins": result.active_bin_count})

    def _candidate(self, centres: np.ndarray, direction: np.ndarray, step: float) -> np.ndarray:
        delta = step * direction
        np.clip(delta, -self.trust[:, None], self.trust[:, None], out=delta)
        delta[self.db.fixed] = 0.0
        return project_centres(self.db, centres + delta)

    def _density_gradient(self, centres: np.ndarray) -> np.ndarray:
        _, _, gradient, _, _ = density_value_gradient(
            centres, self.db.width, self.db.height, self.db.movable,
            self.grid.x_edges, self.grid.y_edges, self.grid.available_area, self.grid.rho_target)
        gradient[self.db.fixed] = 0.0
        return gradient

    def _centre_bins(self, centres: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        ix = np.searchsorted(self.grid.x_edges, centres[:, 0], side="right") - 1
        iy = np.searchsorted(self.grid.y_edges, centres[:, 1], side="right") - 1
        return np.clip(ix, 0, self.grid.nx - 1), np.clip(iy, 0, self.grid.ny - 1)

    def _evacuate_candidate(self, centres: np.ndarray, result: OracleResult, maximum_moves: int = 12_000) -> np.ndarray | None:
        """Build one density-only batch relocation across overlap plateaus.

        The assignment uses centre bins only to propose a large move, but the
        exact overlap oracle remains the sole acceptance test.  HPWL, nets and
        historical placements do not participate in this phase.
        """
        capacity = self.grid.rho_target * self.grid.available_area
        excess = np.maximum(result.occupancy - capacity, 0.0)
        slack = capacity - result.occupancy
        sources = np.flatnonzero(excess.ravel() > 0.0)
        targets = np.flatnonzero(slack.ravel() > 0.0)
        if not len(sources) or not len(targets):
            return None
        ix, iy = self._centre_bins(centres); flat = iy * self.grid.nx + ix
        area = self.db.width * self.db.height
        candidate = centres.copy(); approximate_slack = slack.copy()
        targets = targets[np.argsort(approximate_slack.ravel()[targets])[::-1]]
        moved = 0
        for source in sources[np.argsort(excess.ravel()[sources])[::-1]]:
            if moved >= maximum_moves:
                break
            members = np.flatnonzero(self.db.movable & (flat == source))
            # Move small/low-degree objects first: this reduces the chance that
            # one large cell consumes the remaining slack in a target bin.
            members = members[np.argsort(self.node_degree[members] * np.maximum(area[members], 1.0), kind="stable")]
            for node in members:
                if moved >= maximum_moves:
                    break
                valid = targets[approximate_slack.ravel()[targets] >= area[node]]
                if not len(valid):
                    break
                target = int(valid[0]); ty, tx = divmod(target, self.grid.nx)
                candidate[node] = (.5 * (self.grid.x_edges[tx] + self.grid.x_edges[tx + 1]),
                                   .5 * (self.grid.y_edges[ty] + self.grid.y_edges[ty + 1]))
                approximate_slack[ty, tx] -= area[node]
                moved += 1
        return project_centres(self.db, candidate) if moved else None

    def density_phase(self, centres: np.ndarray, result: OracleResult, iterations: int,
                      *, hpwl_limit: float | None, phase: str, started: float,
                      stop_at_cap: bool) -> tuple[np.ndarray, OracleResult]:
        """Minimise only exact linear overflow; HPWL is merely an optional guard."""
        step = self.cfg.initial_step
        for _ in range(iterations):
            if time.perf_counter() - started >= self.cfg.max_seconds:
                break
            density_grad = self._density_gradient(centres)
            direction = -density_grad / self.precond[:, None]
            candidate = self._candidate(centres, direction, step)
            trial = self.oracle.evaluate(candidate)
            accepted = (trial.density_linear < result.density_linear - 1e-8
                        and (hpwl_limit is None or trial.hpwl <= hpwl_limit + 1e-8))
            if accepted:
                centres, result = candidate, trial
                step = min(step * 1.15, self.cfg.initial_step * 8.0)
                self._record(phase, "density_accept", result)
                if stop_at_cap and self.overflow(result) <= self.cfg.overflow_cap_percent + 1e-9:
                    break
            else:
                step *= .5
                self._record(phase, "density_reject", result)
                # Exact overlap is piecewise linear.  When its local
                # subgradient is stuck, cross several bin boundaries in a
                # density-only batch and accept only a verified reduction.
                evacuation = self._evacuate_candidate(centres, result)
                if evacuation is not None:
                    relocated = self.oracle.evaluate(evacuation)
                    if (relocated.density_linear < result.density_linear - 1e-8
                            and (hpwl_limit is None or relocated.hpwl <= hpwl_limit + 1e-8)):
                        centres, result = evacuation, relocated
                        step = self.cfg.initial_step
                        self._record(phase, "density_evacuation", result)
                        if stop_at_cap and self.overflow(result) <= self.cfg.overflow_cap_percent + 1e-9:
                            break
                if step < 1.0 / 128.0:
                    break
        return centres, result

    def hpwl_phase(self, centres: np.ndarray, result: OracleResult, iterations: int,
                   started: float, phase: str) -> tuple[np.ndarray, OracleResult]:
        """Minimise only HPWL, accepting candidates solely inside the fixed cap."""
        step = self.cfg.initial_step
        cap = self.cfg.overflow_cap_percent
        for _ in range(iterations):
            if time.perf_counter() - started >= self.cfg.max_seconds:
                break
            direction = -result.gradient / self.precond[:, None]
            # On an active density boundary, remove the first-order component
            # that increases exact linear overlap before trust-region clipping.
            density_grad = self._density_gradient(centres)
            directional = float(np.sum(density_grad * direction))
            denom = float(np.sum(density_grad * density_grad / self.precond[:, None]))
            if directional > 0.0 and denom > 1e-12:
                direction -= directional / denom * density_grad / self.precond[:, None]
            candidate = self._candidate(centres, direction, step)
            trial = self.oracle.evaluate(candidate)
            accepted = trial.hpwl < result.hpwl - 1e-8 and self.overflow(trial) <= cap + 1e-9
            if accepted:
                centres, result = candidate, trial
                step = min(step * 1.1, self.cfg.initial_step * 4.0)
                self._record(phase, "hpwl_accept", result)
            else:
                step *= .5
                self._record(phase, "hpwl_reject", result)
                if step < 1.0 / 128.0:
                    break
        return centres, result

    def solve(self, centres: np.ndarray) -> tuple[np.ndarray, OracleResult, list[dict[str, object]], float]:
        centres = project_centres(self.db, centres)
        result = self.oracle.evaluate(centres)
        self._record("input", "seed", result)
        started = time.perf_counter()
        # The initial repair is density-only and establishes feasibility.  The
        # following loop is intentionally strict alternation: HPWL, density,
        # HPWL, density -- never a blended acceptance policy.
        centres, result = self.density_phase(centres, result, self.cfg.density_iterations,
                                             hpwl_limit=None, phase="density_init", started=started,
                                             stop_at_cap=True)
        if self.overflow(result) > self.cfg.overflow_cap_percent + 1e-9:
            raise RuntimeError(f"density phase stopped at {self.overflow(result):.4f}%, above {self.cfg.overflow_cap_percent:.4f}% cap")
        for cycle in range(1, self.cfg.cycles + 1):
            centres, result = self.hpwl_phase(centres, result, self.cfg.hpwl_iterations, started, f"hpwl_{cycle}")
            guard = result.hpwl * (1.0 + self.cfg.hpwl_guard_percent / 100.0)
            centres, result = self.density_phase(centres, result, self.cfg.post_density_iterations,
                                                 hpwl_limit=guard, phase=f"density_{cycle}", started=started,
                                                 stop_at_cap=False)
        return centres, result, self.records, time.perf_counter() - started


def plot_run(root: Path, records: list[dict[str, object]], initial_centres: np.ndarray,
             final_centres: np.ndarray, initial: OracleResult, final: OracleResult, grid, db) -> None:
    """Write convergence and placement/density figures for the isolated run."""
    figures = root / "figures"; figures.mkdir(exist_ok=True)
    plt.rcParams.update({"figure.dpi": 180, "savefig.dpi": 220, "axes.grid": True, "grid.alpha": .2})
    iteration = np.array([float(row["iteration"]) for row in records])
    hpwl = np.array([float(row["hpwl"]) for row in records])
    overflow = np.array([float(row["overflow_percent"]) for row in records])
    overflow_cap = float(records[0]["overflow_cap_percent"])
    phase = [str(row["phase"]) for row in records]
    colors = ["#0072B2" if name.startswith("hpwl") else "#D55E00" for name in phase]
    # Consecutive phase spans make every density/HPWL alternation visible,
    # instead of visually merging all accepted and rejected records together.
    spans: list[tuple[int, int, str]] = []
    start = 0
    for index in range(1, len(phase) + 1):
        if index == len(phase) or phase[index] != phase[start]:
            spans.append((start, index, phase[start])); start = index
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.4), constrained_layout=True)
    for begin, end, name in spans:
        colour = "#0072B2" if name.startswith("hpwl") else "#D55E00"
        axes[0].plot(iteration[begin:end], hpwl[begin:end], color=colour, linewidth=2.2)
        axes[1].plot(iteration[begin:end], overflow[begin:end], color=colour, linewidth=2.2)
        if begin:
            for axis in axes:
                axis.axvline(iteration[begin], color="#888888", linewidth=.55, linestyle=":")
    axes[0].scatter(iteration, hpwl, c=colors, s=10, zorder=3)
    axes[0].yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value / 1e9:.3f}"))
    axes[0].set(xlabel="phase record", ylabel="exact HPWL (×10⁹)", title="Full alternating HPWL trajectory")
    # The initial feasibility repair changes the scale by orders of magnitude.
    # A zoomed inset exposes the later up/down movement of alternating phases.
    feasible_start = next((begin for begin, _, name in spans if name.startswith("hpwl")), 0)
    zoom_x, zoom_y = iteration[feasible_start:], hpwl[feasible_start:]
    inset = inset_axes(axes[0], width="48%", height="45%", loc="lower right", borderpad=1.1)
    for begin, end, name in spans:
        if end <= feasible_start:
            continue
        colour = "#0072B2" if name.startswith("hpwl") else "#D55E00"
        local_begin = max(begin, feasible_start)
        inset.plot(iteration[local_begin:end], hpwl[local_begin:end] / 1e9, color=colour, linewidth=1.5)
    inset.set(xlim=(zoom_x.min(), zoom_x.max()), ylim=(zoom_y.min() / 1e9 - .00005, zoom_y.max() / 1e9 + .00005),
              title="feasible-region zoom")
    inset.tick_params(labelsize=6); inset.grid(alpha=.2)
    axes[1].scatter(iteration, overflow, c=colors, s=10, zorder=3)
    axes[1].axhline(overflow_cap, color="#222222", linestyle="--", label=f"{overflow_cap:g}% cap")
    axes[1].set(xlabel="phase record", ylabel="linear overflow (%)", title="Overflow convergence")
    handles = [Line2D([], [], marker="o", linestyle="", color="#0072B2", label="HPWL phase"),
               Line2D([], [], marker="o", linestyle="", color="#D55E00", label="density phase"),
               Line2D([], [], linestyle="--", color="#222222", label=f"{records[0]['overflow_cap_percent']}% cap")]
    axes[1].legend(handles=handles, fontsize=8)
    fig.savefig(figures / "convergence.png"); plt.close(fig)
    movable = np.flatnonzero(db.movable)
    if len(movable) > 50_000:
        movable = np.random.default_rng(0).choice(movable, 50_000, replace=False)
    fig, axes = plt.subplots(2, 2, figsize=(8, 7), constrained_layout=True)
    for col, (name, centres, oracle) in enumerate((("seed", initial_centres, initial), ("final", final_centres, final))):
        axes[0, col].scatter(centres[movable, 0], centres[movable, 1], s=.3, alpha=.35, rasterized=True)
        axes[0, col].set(title=f"{name} placement", xlabel="x", ylabel="y", aspect="equal")
        density = np.divide(oracle.occupancy, grid.available_area, out=np.full_like(oracle.occupancy, np.nan), where=grid.available_area > 0)
        image = axes[1, col].pcolormesh(grid.x_edges, grid.y_edges, density, shading="auto", cmap="YlOrRd", vmin=0, vmax=2)
        axes[1, col].contour((grid.x_edges[:-1] + grid.x_edges[1:]) / 2, (grid.y_edges[:-1] + grid.y_edges[1:]) / 2,
                             np.nan_to_num(density, nan=-1), levels=[grid.rho_target], colors=["#0072B2"], linewidths=.6)
        axes[1, col].set(title=f"{name} density", xlabel="x", ylabel="y", aspect="equal")
    fig.colorbar(image, ax=axes[1], label="occupancy density")
    fig.savefig(figures / "placement_density.png"); plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Independent alternating overflow/HPWL prototype")
    parser.add_argument("--aux", required=True); parser.add_argument("--seed-placement", required=True)
    parser.add_argument("--cycles", type=int); parser.add_argument("--overflow-cap-percent", type=float)
    parser.add_argument("--output-root"); args = parser.parse_args()
    cfg = AlternatingConfig()
    if args.cycles is not None: cfg.cycles = args.cycles
    if args.overflow_cap_percent is not None: cfg.overflow_cap_percent = args.overflow_cap_percent
    if args.output_root is not None: cfg.output_root = args.output_root
    db = load_bookshelf(args.aux); grid = build_density_grid(db, rho_target=cfg.rho_target, bin_rows=cfg.bin_rows)
    centres = read_placement_centres(args.seed_placement, db)
    solver = AlternatingSolver(db, grid, cfg)
    initial = solver.oracle.evaluate(project_centres(db, centres))
    initial_centres = project_centres(db, centres)
    final_centres, final, records, elapsed = solver.solve(initial_centres)
    root = Path(cfg.output_root) / Path(args.aux).stem / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    root.mkdir(parents=True, exist_ok=False)
    write_placement(root / "solution.pl", db, final_centres)
    with (root / "iterations.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0])); writer.writeheader(); writer.writerows(records)
    metadata = {"status": "completed", "source_aux": str(Path(args.aux).resolve()),
                "seed_placement": str(Path(args.seed_placement).resolve()), "config": asdict(cfg),
                "initial_hpwl": initial.hpwl, "initial_overflow_percent": solver.overflow(initial),
                "final_hpwl": final.hpwl, "final_overflow_percent": solver.overflow(final),
                "iterations": len(records), "elapsed_seconds": elapsed}
    (root / "metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    plot_run(root, records, initial_centres, final_centres, initial, final, grid, db)
    print(json.dumps({"run_dir": str(root), "hpwl": final.hpwl, "overflow_percent": solver.overflow(final)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
