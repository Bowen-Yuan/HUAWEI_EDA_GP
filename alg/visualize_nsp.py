#!/usr/bin/env python3
"""Render optional NSP C++ snapshots as a placement/density/status animation."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from matplotlib import animation
from matplotlib.collections import PatchCollection
from matplotlib.patches import Patch, Rectangle
import matplotlib.pyplot as plt
import numpy as np


PHASE_COLOURS = {"initial": "#555555", "hpwl": "#56B4E9", "coarse": "#E69F00",
                 "medium": "#009E73", "fine": "#CC79A7"}


def read_nodes(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    widths, heights, terminals = [], [], []
    for raw in path.read_text(encoding="ascii").splitlines():
        fields = raw.split()
        if len(fields) < 3 or fields[0] in {"UCLA", "NumNodes", "NumTerminals"}:
            continue
        try:
            widths.append(float(fields[1])); heights.append(float(fields[2]))
        except ValueError:
            continue
        terminals.append(len(fields) >= 4 and fields[3].startswith("terminal"))
    return np.asarray(widths), np.asarray(heights), np.asarray(terminals, dtype=bool)


def read_bounds(path: Path) -> tuple[float, float, float, float]:
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if fields and fields[0] == "bounds":
            return tuple(map(float, fields[1:5]))  # type: ignore[return-value]
    raise ValueError(f"Missing bounds in {path}")


def render(snapshot_dir: Path, benchmark: Path, output: Path, fps: int, max_points: int) -> Path:
    widths, heights, fixed = read_nodes(benchmark.with_suffix(".nodes"))
    bounds = read_bounds(snapshot_dir / "metadata.txt")
    snapshots = sorted(snapshot_dir.glob("snapshot_*.bin"))
    if not snapshots:
        raise ValueError(f"No snapshots in {snapshot_dir}")
    with (snapshot_dir / "metrics.csv").open(newline="", encoding="utf-8") as handle:
        metrics = {int(row["iteration"]): row for row in csv.DictReader(handle)}
    node_count = len(widths)
    movable = np.flatnonzero(~fixed)
    if len(movable) > max_points:
        movable = np.random.default_rng(0).choice(movable, max_points, replace=False)
    first = np.fromfile(snapshots[0], dtype=np.float32).reshape(node_count, 2)
    x0, y0, x1, y1 = bounds
    fixed_ids = np.flatnonzero(fixed)
    rectangles = [Rectangle((first[i, 0], first[i, 1]), widths[i], heights[i]) for i in fixed_ids]

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "figure.dpi": 110})
    bins = 96
    density_maps: list[np.ndarray] = []
    bin_area = ((x1 - x0) / bins) * ((y1 - y0) / bins)
    for path in snapshots:
        positions = np.fromfile(path, dtype=np.float32).reshape(node_count, 2)
        centres = positions + np.column_stack((widths / 2, heights / 2))
        area, _, _ = np.histogram2d(
            centres[~fixed, 1], centres[~fixed, 0], bins=bins,
            range=[[y0, y1], [x0, x1]], weights=(widths[~fixed] * heights[~fixed]),
        )
        density_maps.append(area / bin_area)
    density_vmax = max(0.8, float(np.nanpercentile(np.stack(density_maps), 99.5)))

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    fig.subplots_adjust(left=0.055, right=0.965, bottom=0.13, top=0.88, wspace=0.34)
    placement, density, status = axes
    placement.set(title=f"Cell placement (sample <= {max_points:,})", xlim=(x0, x1), ylim=(y0, y1), aspect="equal")
    cells_artist = placement.scatter([], [], s=0.5, alpha=0.35, color="#0072B2", rasterized=True)
    placement.add_collection(PatchCollection(rectangles, facecolor="#D55E00", edgecolor="#8C2D04", linewidth=0.15, alpha=0.7, rasterized=True))
    placement.legend(handles=[Patch(facecolor="#0072B2", label="Movable cells"), Patch(facecolor="#D55E00", label="Fixed objects")], loc="upper right")
    heatmap = density.imshow(
        np.zeros((bins, bins)), extent=(x0, x1, y0, y1), origin="lower",
        aspect="equal", cmap="YlOrRd", vmin=0.0, vmax=density_vmax,
    )
    density.set(title="Movable-cell density", xlabel="x", ylabel="y")
    fig.colorbar(heatmap, ax=density, label="Area per bin", shrink=0.82)
    status.set_axis_off(); status.set_title("Optimisation state")
    text = status.text(0.05, 0.95, "", va="top", transform=status.transAxes, family="monospace", fontsize=10,
                       bbox={"boxstyle": "round,pad=0.7", "facecolor": "#F7F7F7", "edgecolor": "#D0D0D0"})
    banner = status.text(0.05, 0.12, "", transform=status.transAxes, fontweight="bold")

    def update(frame: int):
        path = snapshots[frame]
        iteration = int(path.stem.rsplit("_", 1)[1])
        positions = np.fromfile(path, dtype=np.float32).reshape(node_count, 2)
        centres = positions + np.column_stack((widths / 2, heights / 2))
        cells_artist.set_offsets(centres[movable])
        heatmap.set_data(density_maps[frame])
        row = metrics.get(iteration, {})
        phase = row.get("phase", "initial")
        text.set_text(f"Iteration     {iteration}\nPhase         {phase.title()}\n\n"
                      f"HPWL          {float(row.get('hpwl', 'nan')):,.3f}\n"
                      f"Overflow      {float(row.get('overflow', 'nan')):.5f}\n"
                      f"Density pen.  {float(row.get('density_penalty', 'nan')):,.3f}\n"
                      f"Lambda        {float(row.get('lambda', 'nan')):.4g}")
        banner.set_text("Snapshot exported by optional C++ visualizer")
        banner.set_color(PHASE_COLOURS.get(phase, "#555555"))
        return cells_artist, heatmap, text, banner

    output.parent.mkdir(parents=True, exist_ok=True)
    movie = animation.FuncAnimation(fig, update, frames=len(snapshots), interval=1000 / fps, blit=False)
    movie.save(output, writer=animation.PillowWriter(fps=fps))
    plt.close(fig)
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshots", type=Path, required=True)
    parser.add_argument("--benchmark", type=Path, required=True, help="Bookshelf base path without .nodes")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--fps", type=int, default=5)
    parser.add_argument("--max-points", type=int, default=100_000)
    args = parser.parse_args()
    output = args.output or args.snapshots / "optimization_animation.gif"
    print(render(args.snapshots, args.benchmark, output, args.fps, args.max_points))


if __name__ == "__main__":
    main()
