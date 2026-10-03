"""Reproducible optimisation animations rendered from saved run snapshots."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from matplotlib import animation
from matplotlib.collections import PatchCollection
from matplotlib.patches import Patch, Rectangle
import matplotlib.pyplot as plt
import numpy as np

from gpplacer.io.bookshelf import load_bookshelf
from gpplacer.oracle.grid import build_density_grid


_PHASE_COLORS = {
    "explore": "#56B4E9",
    "stabilize": "#E69F00",
    "refine": "#CC79A7",
}


def _records(path: Path) -> dict[int, dict[str, str]]:
    # ``utf-8-sig`` also accepts ordinary UTF-8 and makes historical Windows
    # run logs with a BOM use the same ``iteration`` column name.
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return {int(row["iteration"]): row for row in csv.DictReader(handle)}


def _snapshot_paths(root: Path) -> list[tuple[int, Path]]:
    paths: list[tuple[int, Path]] = []
    for path in (root / "snapshots").glob("iteration_*.npz"):
        try:
            paths.append((int(path.stem.rsplit("_", 1)[1]), path))
        except ValueError:
            continue
    paths.sort(key=lambda item: item[0])
    if not paths:
        raise ValueError(f"No iteration snapshots found in {root / 'snapshots'}")
    return paths


def _as_float(record: dict[str, str], name: str) -> float:
    value = record.get(name, "")
    return float(value) if value else float("nan")


def _step_label(iteration: int, record: dict[str, str]) -> str:
    if iteration == 0:
        return "Initial state"
    if record.get("serious_step") == "1":
        return "Serious / accepted step"
    if record.get("phase") == "stabilize":
        return "Null / rejected Bundle step"
    if record.get("phase") == "refine":
        return "Rejected refine step"
    return "Guarded step"


def animate_run(
    run_dir: str | Path,
    output: str | Path | None = None,
    *,
    fps: int = 5,
    max_points: int = 100_000,
    dpi: int = 120,
) -> Path:
    """Create a GIF or MP4 showing placement, density, and solver state.

    The function reads only saved artifacts, so rendering does not affect the
    measured solver runtime.  Large designs use a deterministic movable-cell
    sample while retaining the complete saved density grid.
    """
    if fps <= 0:
        raise ValueError("fps must be positive")
    if max_points <= 0:
        raise ValueError("max_points must be positive")

    root = Path(run_dir)
    metadata = json.loads((root / "metadata.json").read_text(encoding="utf-8"))
    snapshots = _snapshot_paths(root)
    records = _records(root / "iterations.csv")
    db = load_bookshelf(metadata["source_aux"])
    config = metadata["config"]
    grid = build_density_grid(
        db, rho_target=float(config["rho_target"]), bin_rows=int(config["bin_rows"]),
    )
    rho_target = grid.rho_target

    if output is None:
        target = root / "figures" / "optimization_animation.gif"
    else:
        target = Path(output)
    if target.suffix.lower() not in {".gif", ".mp4"}:
        raise ValueError("Animation output must use a .gif or .mp4 suffix")
    target.parent.mkdir(parents=True, exist_ok=True)

    movable = np.flatnonzero(db.movable)
    if len(movable) > max_points:
        movable = np.random.default_rng(0).choice(movable, size=max_points, replace=False)
    fixed = np.flatnonzero(db.fixed)

    density_peak = rho_target
    for _, path in snapshots:
        with np.load(path) as snapshot:
            occupancy = snapshot["occupancy"]
            density = np.divide(
                occupancy, grid.available_area, out=np.full_like(occupancy, np.nan),
                where=grid.available_area > 0,
            )
            if np.isfinite(density).any():
                density_peak = max(density_peak, float(np.nanpercentile(density, 99.5)))
    density_peak = max(density_peak, rho_target * 1.05)

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 10,
        "axes.titleweight": "bold", "figure.dpi": dpi, "savefig.dpi": dpi,
        "axes.spines.top": False, "axes.spines.right": False,
    })
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 4.0))
    fig.subplots_adjust(left=0.055, right=0.965, bottom=0.13, top=0.88, wspace=0.34)
    placement_axis, density_axis, state_axis = axes
    x0, x1, y0, y1 = grid.x_edges[0], grid.x_edges[-1], grid.y_edges[0], grid.y_edges[-1]

    placement_axis.set(
        title=f"Cell placement (sample <= {max_points:,})", xlabel="x", ylabel="y",
        xlim=(x0, x1), ylim=(y0, y1), aspect="equal",
    )
    movable_artist = placement_axis.scatter(
        [], [], s=0.55, alpha=0.35, color="#0072B2", rasterized=True, label="Movable cells",
    )
    if len(fixed):
        fixed_snapshot = np.load(snapshots[0][1])
        try:
            fixed_centres = fixed_snapshot["centres"][fixed]
        finally:
            fixed_snapshot.close()
        fixed_rectangles = [
            Rectangle(
                (fixed_centres[idx, 0] - db.width[cell] / 2.0,
                 fixed_centres[idx, 1] - db.height[cell] / 2.0),
                db.width[cell], db.height[cell],
            )
            for idx, cell in enumerate(fixed)
        ]
        placement_axis.add_collection(PatchCollection(
            fixed_rectangles, facecolor="#D55E00", edgecolor="#8C2D04",
            linewidth=0.15, alpha=0.70, rasterized=True,
        ))
    placement_axis.legend(
        handles=[
            Patch(facecolor="#0072B2", alpha=0.55, label="Movable cells"),
            Patch(facecolor="#D55E00", alpha=0.70, label="Fixed objects"),
        ],
        loc="upper right",
    )

    density_artist = density_axis.imshow(
        np.zeros((grid.ny, grid.nx)), extent=(x0, x1, y0, y1), origin="lower",
        aspect="equal", cmap="YlOrRd", vmin=0.0, vmax=density_peak,
    )
    density_axis.set(title="Bin density and overflow threshold", xlabel="x", ylabel="y")
    fig.colorbar(density_artist, ax=density_axis, label="Occupancy density", shrink=0.82)
    contour_holder: list[object] = []

    state_axis.set_axis_off()
    state_text = state_axis.text(
        0.05, 0.96, "", transform=state_axis.transAxes, va="top", ha="left",
        family="monospace", fontsize=10,
        bbox={"boxstyle": "round,pad=0.7", "facecolor": "#F7F7F7", "edgecolor": "#D0D0D0"},
    )
    banner = state_axis.text(
        0.05, 0.12, "", transform=state_axis.transAxes, va="bottom", ha="left",
        fontsize=10, fontweight="bold",
    )
    state_axis.set_title("Optimisation state")

    def update(frame_index: int):
        iteration, path = snapshots[frame_index]
        record = records.get(iteration, {})
        with np.load(path) as snapshot:
            centres = snapshot["centres"]
            occupancy = snapshot["occupancy"]
        movable_artist.set_offsets(centres[movable])
        density = np.divide(
            occupancy, grid.available_area, out=np.full_like(occupancy, np.nan),
            where=grid.available_area > 0,
        )
        density_artist.set_data(density)
        for collection in contour_holder:
            collection.remove()
        contour_holder.clear()
        if np.nanmax(density) >= rho_target:
            contours = density_axis.contour(
                (grid.x_edges[:-1] + grid.x_edges[1:]) / 2.0,
                (grid.y_edges[:-1] + grid.y_edges[1:]) / 2.0,
                np.nan_to_num(density, nan=-1.0), levels=[rho_target],
                colors=["#0072B2"], linewidths=1.0,
            )
            contour_holder.extend(contours.collections)

        phase = record.get("phase", "initial")
        phase_color = _PHASE_COLORS.get(phase, "#555555")
        state_text.set_text(
            f"Iteration     {iteration}\n"
            f"Phase         {phase.title()}\n\n"
            f"HPWL          {_as_float(record, 'hpwl'):,.3f}\n"
            f"Overflow      {_as_float(record, 'overflow_percent'):.3f}%\n"
            f"Objective     {_as_float(record, 'objective'):,.3f}\n\n"
            f"Active bins   {_as_float(record, 'active_bins'):,.0f}\n"
            f"Step scale    {_as_float(record, 'step'):.3g}\n"
            f"Backtracks    {_as_float(record, 'backtracks'):,.0f}"
        )
        banner.set_text(_step_label(iteration, record))
        banner.set_color(phase_color)
        return movable_artist, density_artist, state_text, banner, *contour_holder

    movie = animation.FuncAnimation(
        fig, update, frames=len(snapshots), interval=1000 / fps, blit=False,
    )
    if target.suffix.lower() == ".gif":
        writer = animation.PillowWriter(fps=fps)
    else:
        if not animation.writers.is_available("ffmpeg"):
            plt.close(fig)
            raise RuntimeError("MP4 export requires ffmpeg; use a .gif output or install ffmpeg")
        writer = animation.FFMpegWriter(fps=fps, codec="libx264")
    movie.save(target, writer=writer, dpi=dpi)
    plt.close(fig)
    return target
