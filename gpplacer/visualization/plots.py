"""Publication-ready data-driven figures generated from saved run artifacts."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np

from gpplacer.io.bookshelf import load_bookshelf
from gpplacer.io.placement import read_placement_centres
from gpplacer.model.config import SolverConfig
from gpplacer.oracle.grid import build_density_grid


_COLORS = {
    "hpwl": "#0072B2", "overflow": "#D55E00", "objective": "#009E73",
    "explore": "#56B4E9", "stabilize": "#E69F00", "refine": "#CC79A7",
}


def _style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9, "axes.labelsize": 9,
        "axes.titlesize": 10, "axes.titleweight": "bold", "legend.fontsize": 8,
        "legend.frameon": False, "figure.dpi": 150, "savefig.dpi": 300,
        "savefig.bbox": "tight", "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.alpha": 0.18, "lines.linewidth": 1.7,
    })


def _records(path: Path) -> dict[str, np.ndarray]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"No iteration records in {path}")
    numeric = [
        "iteration", "elapsed_seconds", "hpwl", "density_linear", "density_surrogate_squared",
        "objective", "step", "proximal_scale", "active_nets_change", "active_bins_change",
        "active_bins", "bundle_size", "backtracks", "serious_step", "recoveries", "rss_bytes",
        "iteration_seconds", "oracle_seconds", "oracle_calls",
        "direction_cosine",
    ]
    result = {
        name: np.asarray([float(row.get(name, "")) if row.get(name, "") else np.nan for row in rows])
        for name in numeric
    }
    result["phase"] = np.asarray([row.get("phase", "unknown") for row in rows])
    return result


def _save(fig: plt.Figure, target: Path) -> None:
    fig.savefig(target.with_suffix(".pdf"))
    fig.savefig(target.with_suffix(".png"), dpi=300)
    plt.close(fig)


def plot_run(run_dir: str | Path) -> list[Path]:
    """Create convergence, phase, activity, and resource figures from one run."""
    _style()
    root = Path(run_dir)
    data = _records(root / "iterations.csv")
    figures = root / "figures"
    figures.mkdir(exist_ok=True)
    written: list[Path] = []
    x = data["iteration"]
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 2.8))
    axes[0].plot(x, data["hpwl"], color=_COLORS["hpwl"], label="HPWL")
    axes[0].set(xlabel="Iteration", ylabel="HPWL", title="Wirelength convergence")
    axes[1].plot(x, data["density_linear"], color=_COLORS["overflow"], label="Linear overflow")
    axes[1].plot(x, data["density_surrogate_squared"], color="#666666", linestyle="--", label="Squared surrogate diagnostic")
    axes[1].set(xlabel="Iteration", ylabel="Density penalty", title="Overflow diagnostics")
    for axis in axes:
        axis.legend()
    target = figures / "convergence"
    _save(fig, target); written.append(target)

    fig, axes = plt.subplots(1, 2, figsize=(8.0, 2.8))
    axes[0].plot(x, data["active_nets_change"], label="Boundary-net change", color="#0072B2")
    axes[0].plot(x, data["active_bins_change"], label="Overflow-bin change", color="#D55E00")
    axes[0].plot(x, data["direction_cosine"], label="Direction cosine", color="#009E73", linestyle="--")
    axes[0].set(xlabel="Iteration", ylabel="Change rate", title="Active-set stability")
    axes[0].legend()
    axes[1].step(x, data["step"], where="mid", label="Explore step", color="#009E73")
    axes[1].plot(x, data["proximal_scale"], label="Proximal scale", color="#CC79A7")
    axes[1].scatter(x[data["serious_step"] > 0], data["proximal_scale"][data["serious_step"] > 0],
                    s=12, color="#222222", label="Serious step", zorder=3)
    axes[1].set(xlabel="Iteration", ylabel="Scale", title="Step and Bundle control")
    axes[1].legend()
    target = figures / "active_sets_and_steps"
    _save(fig, target); written.append(target)

    fig, axes = plt.subplots(1, 2, figsize=(8.0, 2.8))
    elapsed = data["elapsed_seconds"]
    axes[0].plot(elapsed, data["rss_bytes"] / 1024**2, color="#0072B2")
    axes[0].set(xlabel="Wall time (s)", ylabel="RSS (MiB)", title="Memory trace")
    axes[1].plot(x, data["iteration_seconds"], color="#D55E00", label="Total iteration")
    axes[1].plot(x, data["oracle_seconds"], color="#0072B2", label="Oracle")
    axes[1].set(xlabel="Iteration", ylabel="Seconds", title="Iteration time")
    axes[1].legend()
    target = figures / "memory_and_iteration_time"
    _save(fig, target); written.append(target)

    fig, axis = plt.subplots(figsize=(4.0, 2.8))
    phase_code = {"explore": 0, "stabilize": 1, "refine": 2}
    phases = np.asarray([phase_code.get(item, -1) for item in data["phase"]])
    axis.step(x, phases, where="mid", color="#444444")
    axis.set(yticks=[0, 1, 2], yticklabels=["Explore", "Bundle", "Refine"],
             xlabel="Iteration", title="Stage timeline")
    target = figures / "stage_timeline"
    _save(fig, target); written.append(target)
    written.extend(plot_snapshot_comparison(root))
    return written


def plot_snapshot_comparison(run_dir: str | Path) -> list[Path]:
    """Show initial and best density maps from exact solver snapshots."""
    _style()
    root = Path(run_dir)
    metadata = json.loads((root / "metadata.json").read_text(encoding="utf-8"))
    initial_path = root / "snapshots" / "iteration_00000.npz"
    best_path = root / "snapshots" / "best.npz"
    source_aux = Path(metadata["source_aux"])
    if not initial_path.exists() or not best_path.exists() or not source_aux.exists():
        return []
    initial = np.load(initial_path)
    best = np.load(best_path)
    db = load_bookshelf(source_aux)
    config = metadata["config"]
    grid = build_density_grid(db, rho_target=float(config["rho_target"]), bin_rows=int(config["bin_rows"]))
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.4), constrained_layout=True)
    for axis, snapshot, title in zip(axes, (initial, best), ("Initial density", "Best density")):
        density = np.divide(snapshot["occupancy"], grid.available_area,
                            out=np.full_like(snapshot["occupancy"], np.nan),
                            where=grid.available_area > 0)
        image = axis.pcolormesh(grid.x_edges, grid.y_edges, density, shading="auto", cmap="YlOrRd", vmin=0)
        axis.contour((grid.x_edges[:-1] + grid.x_edges[1:]) / 2.0,
                     (grid.y_edges[:-1] + grid.y_edges[1:]) / 2.0,
                     np.nan_to_num(density, nan=-1), levels=[grid.rho_target],
                     colors=["#0072B2"], linewidths=1.0)
        axis.set_aspect("equal", adjustable="box")
        axis.set(title=title, xlabel="x", ylabel="y")
    fig.colorbar(image, ax=axes, label="Occupancy density", shrink=0.85)
    target = root / "figures" / "initial_vs_best_density"
    _save(fig, target)
    return [target]


def plot_placement(aux_path: str | Path, placement_path: str | Path,
                   output_dir: str | Path, *, rho_target: float = 0.8, bin_rows: int = 8,
                   max_points: int = 100_000) -> list[Path]:
    """Render sampled cells and an exact-bin density map for a placement file."""
    _style()
    db = load_bookshelf(aux_path)
    centres = read_placement_centres(placement_path, db)
    grid = build_density_grid(db, rho_target=rho_target, bin_rows=bin_rows)
    from gpplacer.oracle.density import density_value_gradient
    _, _, _, occupancy, active = density_value_gradient(
        centres, db.width, db.height, db.movable, grid.x_edges, grid.y_edges,
        grid.available_area, grid.rho_target,
    )
    density = np.divide(occupancy, grid.available_area, out=np.full_like(occupancy, np.nan),
                        where=grid.available_area > 0)
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)
    movable = np.flatnonzero(db.movable)
    if len(movable) > max_points:
        movable = rng.choice(movable, size=max_points, replace=False)
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.5), constrained_layout=True)
    axes[0].scatter(centres[movable, 0], centres[movable, 1], s=0.4, alpha=0.35,
                    color="#0072B2", rasterized=True, label="Movable cells")
    fixed = np.flatnonzero(db.fixed)
    axes[0].scatter(centres[fixed, 0], centres[fixed, 1], s=4, color="#D55E00",
                    label="Fixed objects")
    axes[0].set_aspect("equal", adjustable="box")
    axes[0].set(title="Placement sample", xlabel="x", ylabel="y")
    axes[0].legend(markerscale=4)
    image = axes[1].pcolormesh(grid.x_edges, grid.y_edges, density, shading="auto", cmap="YlOrRd")
    axes[1].contour((grid.x_edges[:-1] + grid.x_edges[1:]) / 2.0,
                    (grid.y_edges[:-1] + grid.y_edges[1:]) / 2.0,
                    np.nan_to_num(density, nan=-1), levels=[rho_target], colors=["#0072B2"], linewidths=1.0)
    axes[1].set_aspect("equal", adjustable="box")
    axes[1].set(title=f"Density map (threshold = {rho_target:.2f})", xlabel="x", ylabel="y")
    fig.colorbar(image, ax=axes[1], label="Occupancy density")
    target = output / "placement_and_density"
    _save(fig, target)
    written = [target]
    pin_positions = centres[db.pin_node] + db.pin_offset
    starts = db.net_start[:-1]
    min_x = np.minimum.reduceat(pin_positions[:, 0], starts)
    max_x = np.maximum.reduceat(pin_positions[:, 0], starts)
    min_y = np.minimum.reduceat(pin_positions[:, 1], starts)
    max_y = np.maximum.reduceat(pin_positions[:, 1], starts)
    non_degenerate = (max_x > min_x) & (max_y > min_y)
    box_count = int(np.count_nonzero(non_degenerate))
    segments = np.empty((box_count * 4, 2, 2), dtype=np.float32)
    indices = np.flatnonzero(non_degenerate)
    left, right, bottom, top = min_x[indices], max_x[indices], min_y[indices], max_y[indices]
    segments[0::4, 0] = np.column_stack((left, bottom)); segments[0::4, 1] = np.column_stack((right, bottom))
    segments[1::4, 0] = np.column_stack((right, bottom)); segments[1::4, 1] = np.column_stack((right, top))
    segments[2::4, 0] = np.column_stack((right, top)); segments[2::4, 1] = np.column_stack((left, top))
    segments[3::4, 0] = np.column_stack((left, top)); segments[3::4, 1] = np.column_stack((left, bottom))
    difference = np.zeros((grid.ny + 1, grid.nx + 1), dtype=np.int64)
    ix0 = np.clip(np.searchsorted(grid.x_edges, min_x[indices], side="right") - 1, 0, grid.nx - 1)
    ix1 = np.clip(np.searchsorted(grid.x_edges, np.nextafter(max_x[indices], -np.inf), side="right") - 1, 0, grid.nx - 1)
    iy0 = np.clip(np.searchsorted(grid.y_edges, min_y[indices], side="right") - 1, 0, grid.ny - 1)
    iy1 = np.clip(np.searchsorted(grid.y_edges, np.nextafter(max_y[indices], -np.inf), side="right") - 1, 0, grid.ny - 1)
    np.add.at(difference, (iy0, ix0), 1)
    np.add.at(difference, (iy1 + 1, ix0), -1)
    np.add.at(difference, (iy0, ix1 + 1), -1)
    np.add.at(difference, (iy1 + 1, ix1 + 1), 1)
    coverage = difference.cumsum(axis=0).cumsum(axis=1)[:-1, :-1]
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.5), constrained_layout=True)
    axes[0].add_collection(LineCollection(
        segments, colors="#3B3B3B", linewidths=0.08, alpha=0.015, rasterized=True,
    ))
    axes[0].set_xlim(grid.x_edges[0], grid.x_edges[-1]); axes[0].set_ylim(grid.y_edges[0], grid.y_edges[-1])
    axes[0].set_aspect("equal", adjustable="box")
    axes[0].set(title=f"All net bounding boxes ({box_count:,})", xlabel="x", ylabel="y")
    image = axes[1].pcolormesh(grid.x_edges, grid.y_edges, np.log1p(coverage), shading="auto", cmap="magma")
    axes[1].set_aspect("equal", adjustable="box")
    axes[1].set(title="Bounding-box overlap coverage", xlabel="x", ylabel="y")
    fig.colorbar(image, ax=axes[1], label="log(1 + covering net boxes)")
    target = output / "net_bounding_boxes_and_overlap"
    _save(fig, target)
    written.append(target)

    # A full-chip rendering of hundreds of thousands of boxes necessarily
    # becomes a hairball.  Following the reference paper's geometric
    # comparison style, use length-stratified representative panels so that
    # individual rectangles are actually readable.  The preceding heatmap
    # remains the exact all-net overview.
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.4), constrained_layout=True)
    lengths = (max_x - min_x) + (max_y - min_y)
    valid_lengths = lengths[indices]
    quantiles = np.quantile(valid_lengths, [0.0, 0.25, 0.50, 0.75, 1.0])
    rng = np.random.default_rng(0)
    colors = ("#0072B2", "#009E73", "#E69F00", "#D55E00")
    for panel, low, high, color in zip(axes.ravel(), quantiles[:-1], quantiles[1:], colors):
        selected_indices = indices[(valid_lengths >= low) & (valid_lengths <= high)]
        if len(selected_indices) > 600:
            selected_indices = rng.choice(selected_indices, size=600, replace=False)
        local_segments = np.empty((len(selected_indices) * 4, 2, 2), dtype=np.float32)
        local_left, local_right = min_x[selected_indices], max_x[selected_indices]
        local_bottom, local_top = min_y[selected_indices], max_y[selected_indices]
        local_segments[0::4, 0] = np.column_stack((local_left, local_bottom)); local_segments[0::4, 1] = np.column_stack((local_right, local_bottom))
        local_segments[1::4, 0] = np.column_stack((local_right, local_bottom)); local_segments[1::4, 1] = np.column_stack((local_right, local_top))
        local_segments[2::4, 0] = np.column_stack((local_right, local_top)); local_segments[2::4, 1] = np.column_stack((local_left, local_top))
        local_segments[3::4, 0] = np.column_stack((local_left, local_top)); local_segments[3::4, 1] = np.column_stack((local_left, local_bottom))
        panel.add_collection(LineCollection(
            local_segments, colors=color, linewidths=0.35, alpha=0.38, rasterized=True,
        ))
        panel.set(xlim=(grid.x_edges[0], grid.x_edges[-1]), ylim=(grid.y_edges[0], grid.y_edges[-1]),
                  aspect="equal", title=f"HPWL span [{low:.1f}, {high:.1f}]: {len(selected_indices):,} nets",
                  xlabel="x", ylabel="y")
    target = output / "net_bounding_box_length_strata"
    _save(fig, target)
    written.append(target)
    return written
