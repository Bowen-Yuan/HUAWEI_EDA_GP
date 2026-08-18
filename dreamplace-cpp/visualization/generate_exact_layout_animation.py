#!/usr/bin/env python3
"""Render real DREAMPlace C++ global snapshots and legalization stages."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Patch
import numpy as np


def read_nodes(path: Path):
    widths = {}
    heights = {}
    fixed = {}
    for raw in path.read_text(encoding="ascii", errors="ignore").splitlines():
        line = raw.strip()
        if not line or line.startswith("UCLA") or line.startswith("#") or line.startswith("Num"):
            continue
        fields = line.split()
        if len(fields) < 3:
            continue
        try:
            widths[fields[0]] = float(fields[1])
            heights[fields[0]] = float(fields[2])
        except ValueError:
            continue
        fixed[fields[0]] = any("terminal" in token.lower() for token in fields[3:])
    return widths, heights, fixed


def read_pl(path: Path, names, widths, heights):
    values = np.full((len(names), 2), np.nan, dtype=np.float64)
    index = {name: i for i, name in enumerate(names)}
    for raw in path.read_text(encoding="ascii", errors="ignore").splitlines():
        fields = raw.split()
        if len(fields) < 3 or fields[0] not in index:
            continue
        try:
            left = float(fields[1])
            bottom = float(fields[2])
        except ValueError:
            continue
        i = index[fields[0]]
        values[i] = (left + widths[fields[0]] * 0.5,
                     bottom + heights[fields[0]] * 0.5)
    if np.isnan(values).any():
        missing = int(np.isnan(values[:, 0]).sum())
        raise ValueError(f"{path} is missing {missing} node locations")
    return values


def read_metrics(path: Path):
    with path.open(newline="", encoding="ascii") as stream:
        rows = list(csv.DictReader(stream))
    return rows


def read_nets(path: Path):
    nets = []
    current = []
    for raw in path.read_text(encoding="ascii", errors="ignore").splitlines():
        line = raw.strip()
        if line.startswith("NetDegree"):
            if len(current) >= 2:
                nets.append(current)
            current = []
            continue
        fields = line.split()
        if len(fields) >= 5 and fields[0] in _NODE_INDEX and fields[2] == ":":
            try:
                current.append((fields[0], float(fields[3]), float(fields[4])))
            except ValueError:
                pass
    if len(current) >= 2:
        nets.append(current)
    return nets


_NODE_INDEX = {}


def exact_hpwl(positions, nets):
    total = 0.0
    for net in nets:
        xs = [positions[_NODE_INDEX[name], 0] + ox for name, ox, _ in net]
        ys = [positions[_NODE_INDEX[name], 1] + oy for name, _, oy in net]
        total += max(xs) - min(xs) + max(ys) - min(ys)
    return total


def read_summary(path: Path):
    data = {}
    for line in path.read_text(encoding="ascii").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            data[key] = value
    return data


def render(raw_base: Path, run_dir: Path, output: Path, still: Path, fps: int,
           all_cells: bool, max_display_cells: int):
    dataset_name = raw_base.name
    nodes_path = raw_base.with_suffix(".nodes")
    widths, heights, fixed_map = read_nodes(nodes_path)
    names = list(widths)
    global _NODE_INDEX
    _NODE_INDEX = {name: i for i, name in enumerate(names)}
    movable_names = [name for name in names if not fixed_map[name]]
    fixed_names = [name for name in names if fixed_map[name]]
    movable_idx = [_NODE_INDEX[name] for name in movable_names]
    display_stride = (1 if all_cells else
                      max(1, int(np.ceil(len(movable_idx) /
                                         max_display_cells))))
    display_idx = movable_idx[::display_stride]
    point_size = 0.18 if all_cells else 0.32
    point_alpha = 0.42 if all_cells else 0.56

    snapshots_dir = run_dir / "snapshots"
    global_paths = sorted(snapshots_dir.glob("global_[0-9][0-9][0-9][0-9].pl"))
    if not global_paths:
        raise FileNotFoundError(f"no global snapshots in {snapshots_dir}")
    stage_names = ["greedy", "abacus", "adjacent", "k_reorder",
                   "global_swap", "independent_set"]
    stage_paths = [snapshots_dir / f"legal_{name}.pl" for name in stage_names]
    if not all(path.exists() for path in stage_paths):
        missing = [path.name for path in stage_paths if not path.exists()]
        raise FileNotFoundError(f"missing legalization snapshots: {missing}")
    selected_path = snapshots_dir / "global_selected.pl"
    if selected_path.exists():
        global_paths.append(selected_path)

    metrics = read_metrics(run_dir / "global_metrics.csv")
    summary = read_summary(run_dir / "summary.txt")
    bounds = np.array([459.0, 459.0, 11151.0, 11139.0])
    # Use the actual placement bounds from all snapshot coordinates.
    first = read_pl(global_paths[0], names, widths, heights)
    bounds[0] = min(bounds[0], float(first[:, 0].min()))
    bounds[1] = min(bounds[1], float(first[:, 1].min()))
    bounds[2] = max(bounds[2], float(first[:, 0].max()))
    bounds[3] = max(bounds[3], float(first[:, 1].max()))

    global_frames = []
    for path in global_paths:
        match = re.search(r"global_(\d{4})", path.name)
        iteration = (int(match.group(1)) if match else
                     int(summary["selected_legal_iteration"]))
        if path.name == "global_selected.pl":
            hpwl = float(summary["gp_hpwl"])
            overflow = float(summary["gp_overflow"])
            label = "Best feasible GP"
            lambda_value = float(metrics[min(iteration, len(metrics) - 1)]["lambda_effective"])
        else:
            row = metrics[iteration]
            hpwl = float(row["exact_hpwl"])
            overflow = float(row["overflow"])
            lambda_value = float(row["lambda_effective"])
            label = f"Global placement, iter {iteration}"
        global_frames.append((path, label, hpwl, overflow, lambda_value, None,
                              iteration))

    final_iteration = int(metrics[-1]["iteration"])
    stage_frames = []
    for path, name in zip(stage_paths, stage_names):
        stage_label = f"Legalization: {name.replace('_', ' ').title()}"
        if name == "independent_set":
            stage_label += " (legal result)"
        stage_frames.append((path, stage_label, 0.0,
                             float(summary["gp_overflow"]),
                             float(metrics[-1]["lambda_effective"]),
                             stage_names.index(name), final_iteration))

    frames = global_frames + stage_frames
    positions = [read_pl(path, names, widths, heights) for path, *_ in frames]
    # Most legalization-stage HPWL values are already emitted by the solver.
    # Only the adjacent-swap stage lacks a dedicated summary field, so avoid
    # recomputing exact HPWL for all six large placements in Python.
    stage_summary_keys = {
        "greedy": "greedy_hpwl",
        "abacus": "abacus_hpwl",
        "k_reorder": "k_reorder_hpwl",
        "global_swap": "global_swap_hpwl",
        "independent_set": "independent_set_hpwl",
    }
    stage_values = []
    nets = None
    for stage_offset, name in enumerate(stage_names):
        summary_key = stage_summary_keys.get(name)
        if summary_key is not None:
            stage_values.append(float(summary[summary_key]) / 1e6)
            continue
        if nets is None:
            nets = read_nets(raw_base.with_suffix(".nets"))
        position_offset = len(global_frames) + stage_offset
        stage_values.append(exact_hpwl(positions[position_offset], nets) / 1e6)
    for i, value in enumerate(stage_values):
        path, label, _, overflow_value, lambda_value, stage_index, iteration = frames[len(global_frames) + i]
        frames[len(global_frames) + i] = (path, label, value * 1e6,
                                          overflow_value, lambda_value,
                                          stage_index, iteration)
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.18,
        "figure.facecolor": "white",
    })
    fig = plt.figure(figsize=(13.2, 7.8), dpi=100, constrained_layout=True)
    grid = fig.add_gridspec(3, 2, width_ratios=(1.3, 1.0), height_ratios=(1, 1, 1.05))
    ax_place = fig.add_subplot(grid[:, 0])
    ax_hpwl = fig.add_subplot(grid[0, 1])
    ax_overflow = fig.add_subplot(grid[1, 1], sharex=ax_hpwl)
    ax_stage = fig.add_subplot(grid[2, 1])

    iterations = np.array([float(row["iteration"]) for row in metrics])
    hpwl = np.array([float(row["exact_hpwl"]) for row in metrics]) / 1e6
    overflow = np.array([float(row["overflow"]) for row in metrics]) * 100.0
    lambda_values = np.maximum(
        np.array([float(row["lambda_effective"]) for row in metrics]), 1e-15)
    ax_hpwl.plot(iterations, hpwl, color="#0072B2", lw=1.5)
    ax_overflow.plot(iterations, overflow, color="#D55E00", lw=1.5)
    ax_overflow.axhline(7.0, color="#555555", ls="--", lw=1.0)
    ax_hpwl.set_ylabel("Exact HPWL (M)")
    ax_overflow.set_ylabel("Overflow (%)")
    ax_hpwl.set_title("Convergence")
    ax_overflow.set_title("Density constraint")
    ax_stage.set_title("Legalization HPWL (M)")
    ax_hpwl.set_xlim(0, final_iteration)
    ax_overflow.set_xlim(0, final_iteration)
    cursor_hpwl = ax_hpwl.axvline(0, color="#333333", lw=1.0)
    cursor_overflow = ax_overflow.axvline(0, color="#333333", lw=1.0)

    fixed_patches = []
    for name in fixed_names:
        i = _NODE_INDEX[name]
        x, y = positions[0][i]
        fixed_patches.append(Rectangle(
            (x - widths[name] / 2, y - heights[name] / 2),
            widths[name], heights[name]))
    fixed_collection = PatchCollection(
        fixed_patches, facecolor="#E69F00", alpha=0.78,
        edgecolor="#8C5A00", linewidth=0.2)
    ax_place.add_collection(fixed_collection)
    movable_scatter = ax_place.scatter(
        positions[0][display_idx, 0], positions[0][display_idx, 1],
        s=point_size, c="#0072B2", alpha=point_alpha, linewidths=0,
        rasterized=True)
    ax_place.set_xlim(bounds[0], bounds[2])
    ax_place.set_ylim(bounds[1], bounds[3])
    ax_place.set_aspect("equal", adjustable="box")
    ax_place.set_xlabel("x")
    ax_place.set_ylabel("y")
    placement_title = ax_place.set_title("")
    ax_place.legend(handles=[
        Patch(facecolor="#0072B2", alpha=0.7, label="Movable cells"),
        Patch(facecolor="#E69F00", alpha=0.78, label="Fixed macros"),
    ], loc="upper right", frameon=True, fontsize=8)

    stage_bars = ax_stage.barh(
        np.arange(len(stage_names)), stage_values, color="#9AA5B1")
    ax_stage.set_yticks(np.arange(len(stage_names)))
    ax_stage.set_yticklabels(
        [name.replace("_", " ").title() for name in stage_names], fontsize=8)
    ax_stage.set_ylim(len(stage_names) - 0.5, -0.5)
    ax_stage.set_xlim(0, max(max(stage_values), 1.0) * 1.12)
    ax_stage.set_xlabel("HPWL (M)")
    ax_stage.grid(axis="x", alpha=0.18)
    for i, value in enumerate(stage_values):
        ax_stage.text(value + 0.1, i, f"{value:.3f}", va="center", fontsize=7)
    figure_title = fig.suptitle("", fontsize=12, fontweight="bold")

    def draw(frame_index: int):
        path, label, frame_hpwl, frame_overflow, frame_lambda, stage_index, iteration = frames[frame_index]
        pos = positions[frame_index]
        movable_scatter.set_offsets(pos[display_idx])
        placement_title.set_text(label)

        cursor_hpwl.set_xdata([iteration, iteration])
        cursor_overflow.set_xdata([iteration, iteration])
        colors = ["#9AA5B1"] * len(stage_names)
        if stage_index is not None:
            colors[:stage_index] = ["#009E73"] * stage_index
            colors[stage_index] = "#D55E00"
        for bar, color in zip(stage_bars, colors):
            bar.set_color(color)
        status = " | LEGAL" if stage_index == len(stage_names) - 1 else ""
        figure_title.set_text(
            f"{dataset_name} | exact nonsmooth HPWL + electrostatic density | "
            f"HPWL {frame_hpwl / 1e6:.3f}M | overflow {frame_overflow * 100:.3f}% | "
            f"lambda {frame_lambda:.3e}{status}")

    output.parent.mkdir(parents=True, exist_ok=True)
    still.parent.mkdir(parents=True, exist_ok=True)
    writer = PillowWriter(fps=fps)
    with writer.saving(fig, output, dpi=100):
        for frame_index in range(len(global_frames)):
            draw(frame_index)
            writer.grab_frame()
        for frame_index in range(len(global_frames), len(frames)):
            for _ in range(2):
                draw(frame_index)
                writer.grab_frame()
        for _ in range(12):
            draw(len(frames) - 1)
            writer.grab_frame()
    draw(len(frames) - 1)
    fig.savefig(still, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("raw_base", type=Path)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--still", type=Path, required=True)
    parser.add_argument("--fps", type=int, default=6)
    parser.add_argument(
        "--all-cells", action="store_true",
        help="draw every movable cell instead of deterministic display sampling")
    parser.add_argument(
        "--max-display-cells", type=int, default=12000,
        help="maximum movable cells to draw when sampling (default: 12000)")
    args = parser.parse_args()
    if args.max_display_cells <= 0:
        parser.error("--max-display-cells must be positive")
    render(args.raw_base, args.run_dir, args.output, args.still, args.fps,
           args.all_cells, args.max_display_cells)


if __name__ == "__main__":
    main()
