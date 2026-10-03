"""Reproducible numerical figures for safe-anchor runs."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from gpplacer.io.bookshelf import load_bookshelf
from gpplacer.oracle.grid import build_density_grid


def _style() -> None:
    plt.rcParams.update({"font.family": "DejaVu Serif", "font.size": 9, "figure.dpi": 300,
                         "savefig.dpi": 300, "savefig.bbox": "tight", "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True, "grid.alpha": .16})


def _save(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path.with_suffix(".png"), dpi=300)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def _rows(run_dir: Path) -> list[dict[str, str]]:
    with (run_dir / "iterations.csv").open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def plot_run(run_dir: str | Path) -> list[Path]:
    """Plot convergence and Pareto trajectory from a completed run."""
    _style(); root = Path(run_dir); rows = _rows(root); out = root / "figures"; out.mkdir(exist_ok=True)
    values = lambda key: np.asarray([float(row[key]) for row in rows])
    x, hpwl, overflow = values("iteration"), values("hpwl"), values("overflow_percent")
    budget, trust50, trust90 = values("budget_percent"), values("trust_p50"), values("trust_p90")
    serious, pressure = values("serious_step"), values("pressure")
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 2.9), constrained_layout=True)
    axes[0].plot(x, hpwl, color="#0072B2", label="Exact HPWL")
    axes[0].scatter(x[serious > 0], hpwl[serious > 0], s=12, color="#009E73", label="Serious step")
    axes[0].set(xlabel="Iteration", ylabel="HPWL", title="Constrained wirelength convergence")
    axes[1].plot(x, overflow, color="#D55E00", label="Exact overflow (%)")
    axes[1].plot(x, budget, color="#333333", linestyle="--", label="Configured overflow cap")
    axes[1].scatter(x[pressure > 0], overflow[pressure > 0], s=18, marker="D", color="#CC79A7", label="Pressure event")
    axes[1].set(xlabel="Iteration", ylabel="Overflow (%)", title="Overflow remains below configured cap")
    axes[1].ticklabel_format(axis="y", style="plain", useOffset=False)
    for ax in axes: ax.legend(fontsize=7)
    convergence = out / "convergence"; _save(fig, convergence)
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 2.9), constrained_layout=True)
    axes[0].plot(x, trust50, color="#E69F00", label="Trust radius p50")
    axes[0].plot(x, trust90, color="#56B4E9", label="Trust radius p90")
    axes[0].set(xlabel="Iteration", ylabel="Distance", title="Component-wise trust regions")
    axes[0].legend(fontsize=7)
    axes[1].plot(overflow, hpwl, color="#666666", marker="o", markersize=2.5)
    axes[1].scatter(overflow[-1], hpwl[-1], color="#E76F51", s=30, label="Final")
    axes[1].set(xlabel="Overflow (%)", ylabel="HPWL", title="HPWL-overflow trajectory")
    axes[1].legend(fontsize=7)
    diagnostics = out / "trust_and_pareto"; _save(fig, diagnostics)
    return [convergence, diagnostics]


def plot_layout(run_dir: str | Path, *, max_points: int = 50_000) -> Path:
    """Render input, safe anchor, and final placement/density distributions."""
    _style(); root = Path(run_dir); meta = json.loads((root / "metadata.json").read_text(encoding="utf-8"))
    db = load_bookshelf(meta["source_aux"]); cfg = meta["config"]
    snapshots = np.load(root / "snapshots.npz")
    grid = build_density_grid(db, rho_target=float(cfg["rho_target"]), bin_rows=int(cfg["bin_rows"]))
    rng = np.random.default_rng(0); movable = np.flatnonzero(db.movable)
    if len(movable) > max_points: movable = rng.choice(movable, max_points, replace=False)
    fig, axes = plt.subplots(2, 3, figsize=(10.2, 5.6), constrained_layout=True)
    image = None
    for col, key, title in zip(range(3), ("initial", "anchor", "final"), ("Initial", "Safe anchor", "Final")):
        centres, occupancy = snapshots[f"{key}_centres"], snapshots[f"{key}_occupancy"]
        axes[0, col].scatter(centres[movable, 0], centres[movable, 1], s=.25, alpha=.32, color="#0072B2", rasterized=True)
        axes[0, col].set(title=title, xlabel="x", ylabel="y", aspect="equal")
        density = np.divide(occupancy, grid.available_area, out=np.full_like(occupancy, np.nan), where=grid.available_area > 0)
        image = axes[1, col].pcolormesh(grid.x_edges, grid.y_edges, density, shading="auto", cmap="YlOrRd", vmin=0)
        axes[1, col].contour((grid.x_edges[:-1] + grid.x_edges[1:]) / 2, (grid.y_edges[:-1] + grid.y_edges[1:]) / 2,
                             np.nan_to_num(density, nan=-1), levels=[grid.rho_target], colors=["#0072B2"], linewidths=.8)
        axes[1, col].set(title=f"{title} density", xlabel="x", ylabel="y", aspect="equal")
    fig.colorbar(image, ax=axes[1], shrink=.78, label="Occupancy density")
    target = root / "figures" / "placement_and_density"; _save(fig, target)
    return target
