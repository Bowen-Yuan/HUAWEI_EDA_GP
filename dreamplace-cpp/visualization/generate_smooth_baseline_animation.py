#!/usr/bin/env python3
"""Render all-cell smooth DREAMPlace C++ convergence and legalization."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter
from matplotlib.collections import PatchCollection
from matplotlib.patches import Patch, Rectangle
import numpy as np


def read_nodes(path: Path):
    widths, heights, fixed = {}, {}, {}
    for raw in path.read_text(encoding="ascii", errors="ignore").splitlines():
        fields = raw.split()
        if len(fields) < 3 or fields[0] in {"UCLA", "NumNodes", "NumTerminals"}:
            continue
        try:
            widths[fields[0]] = float(fields[1])
            heights[fields[0]] = float(fields[2])
        except ValueError:
            continue
        fixed[fields[0]] = any("terminal" in token.lower() for token in fields[3:])
    return widths, heights, fixed


def read_pl(path: Path, names, widths, heights, index):
    values = np.full((len(names), 2), np.nan, dtype=np.float64)
    for raw in path.read_text(encoding="ascii", errors="ignore").splitlines():
        fields = raw.split()
        if len(fields) < 3 or fields[0] not in index:
            continue
        try:
            left, bottom = float(fields[1]), float(fields[2])
        except ValueError:
            continue
        i = index[fields[0]]
        values[i] = (left + 0.5 * widths[fields[0]],
                     bottom + 0.5 * heights[fields[0]])
    if np.isnan(values).any():
        raise ValueError(f"{path} has missing node coordinates")
    return values


def read_key_values(path: Path):
    result = {}
    for line in path.read_text(encoding="ascii", errors="ignore").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            result[key] = value
    return result


def read_metrics(path: Path):
    with path.open(newline="", encoding="ascii") as stream:
        return list(csv.DictReader(stream))


def render(raw_base: Path, run_dir: Path, output: Path, still: Path, fps: int):
    widths, heights, fixed_map = read_nodes(raw_base.with_suffix(".nodes"))
    names = list(widths)
    index = {name: i for i, name in enumerate(names)}
    movable_idx = np.asarray([index[name] for name in names if not fixed_map[name]])
    fixed_names = [name for name in names if fixed_map[name]]

    snapshots = run_dir / "snapshots"
    global_paths = sorted(snapshots.glob("global_[0-9][0-9][0-9][0-9].pl"))
    selected = snapshots / "global_selected.pl"
    if not global_paths or not selected.exists():
        raise FileNotFoundError(f"missing global snapshots in {snapshots}")
    global_paths.append(selected)

    stage_names = ["greedy", "abacus", "adjacent", "k_reorder",
                   "global_swap", "independent_set"]
    stage_paths = [snapshots / f"legal_{name}.pl" for name in stage_names]
    missing = [path.name for path in stage_paths if not path.exists()]
    if missing:
        raise FileNotFoundError(f"missing legal snapshots: {missing}")

    metrics = read_metrics(run_dir / "global_metrics.csv")
    summary = read_key_values(run_dir / "summary.txt")
    benchmark = summary.get("benchmark", raw_base.name)
    all_paths = global_paths + stage_paths
    positions = [read_pl(path, names, widths, heights, index) for path in all_paths]

    all_xy = np.vstack((positions[0], positions[-1]))
    margin_x = max(1.0, 0.015 * (all_xy[:, 0].max() - all_xy[:, 0].min()))
    margin_y = max(1.0, 0.015 * (all_xy[:, 1].max() - all_xy[:, 1].min()))
    bounds = (all_xy[:, 0].min() - margin_x, all_xy[:, 0].max() + margin_x,
              all_xy[:, 1].min() - margin_y, all_xy[:, 1].max() + margin_y)

    iterations = np.asarray([float(row["iteration"]) for row in metrics])
    hpwl_curve = np.asarray([float(row["exact_hpwl"]) for row in metrics]) / 1e6
    overflow_curve = np.asarray([float(row["overflow"]) for row in metrics]) * 100
    lambda_curve = np.maximum(
        np.asarray([float(row["lambda"]) for row in metrics]), 1e-30)

    stage_keys = ["greedy_hpwl", "abacus_hpwl", "detailed_hpwl",
                  "k_reorder_hpwl", "global_swap_hpwl", "independent_set_hpwl"]
    stage_values = []
    fallback = float(summary["detailed_hpwl"])
    for key in stage_keys:
        stage_values.append(float(summary.get(key, fallback)) / 1e6)

    frames = []
    for path in global_paths[:-1]:
        match = re.search(r"global_(\d{4})", path.name)
        iteration = int(match.group(1))
        row = metrics[min(iteration, len(metrics) - 1)]
        frames.append((f"Global placement, iter {iteration}", iteration,
                       float(row["exact_hpwl"]) / 1e6,
                       100 * float(row["overflow"]), None))
    selected_iter = int(float(summary.get("selected_iteration", len(metrics) - 1)))
    if selected_iter < 0:
        selected_iter = len(metrics) - 1
    frames.append(("Selected global placement", selected_iter,
                   float(summary["gp_hpwl"]) / 1e6,
                   100 * float(summary["gp_overflow"]), None))
    for i, name in enumerate(stage_names):
        frames.append((f"Legalization: {name.replace('_', ' ').title()}",
                       selected_iter, stage_values[i],
                       100 * float(summary["gp_overflow"]), i))

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.titlesize": 10, "axes.titleweight": "bold",
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.alpha": 0.16,
    })
    fig = plt.figure(figsize=(13.2, 7.8), dpi=100, constrained_layout=True)
    grid = fig.add_gridspec(3, 2, width_ratios=(1.35, 1.0))
    ax_place = fig.add_subplot(grid[:, 0])
    ax_hpwl = fig.add_subplot(grid[0, 1])
    ax_overflow = fig.add_subplot(grid[1, 1], sharex=ax_hpwl)
    ax_stage = fig.add_subplot(grid[2, 1])

    ax_hpwl.plot(iterations, hpwl_curve, color="#0072B2", lw=1.6)
    ax_overflow.plot(iterations, overflow_curve, color="#D55E00", lw=1.6)
    ax_overflow.axhline(7.0, color="#555555", ls="--", lw=1.0)
    ax_hpwl.set_ylabel("Exact HPWL (M)")
    ax_overflow.set_ylabel("Overflow (%)")
    ax_overflow.set_xlabel("Iteration")
    ax_hpwl.set_title("Smooth GP convergence")
    ax_overflow.set_title("Electrostatic density")
    cursor_h = ax_hpwl.axvline(0, color="#333333", lw=1.0)
    cursor_o = ax_overflow.axvline(0, color="#333333", lw=1.0)

    first = positions[0]
    fixed_patches = []
    for name in fixed_names:
        i = index[name]
        fixed_patches.append(Rectangle(
            (first[i, 0] - 0.5 * widths[name], first[i, 1] - 0.5 * heights[name]),
            widths[name], heights[name]))
    ax_place.add_collection(PatchCollection(
        fixed_patches, facecolor="#E69F00", edgecolor="#8C5A00",
        alpha=0.78, linewidth=0.18))
    point_size = 0.14 if len(movable_idx) > 350000 else 0.20
    movable = ax_place.scatter(
        first[movable_idx, 0], first[movable_idx, 1], s=point_size,
        c="#0072B2", alpha=0.42, linewidths=0, rasterized=True)
    ax_place.set_xlim(bounds[0], bounds[1])
    ax_place.set_ylim(bounds[2], bounds[3])
    ax_place.set_aspect("equal", adjustable="box")
    ax_place.set_xlabel("x")
    ax_place.set_ylabel("y")
    placement_title = ax_place.set_title("")
    ax_place.legend(handles=[
        Patch(facecolor="#0072B2", label="Movable cells"),
        Patch(facecolor="#E69F00", label="Fixed macros")
    ], loc="upper right", fontsize=8)

    bars = ax_stage.barh(np.arange(len(stage_names)), stage_values,
                         color="#9AA5B1")
    ax_stage.set_yticks(np.arange(len(stage_names)))
    ax_stage.set_yticklabels([name.replace("_", " ").title()
                              for name in stage_names], fontsize=8)
    ax_stage.invert_yaxis()
    ax_stage.set_xlim(0, max(stage_values) * 1.14)
    ax_stage.set_xlabel("Exact HPWL (M)")
    ax_stage.set_title("Legalization and detailed placement")
    for i, value in enumerate(stage_values):
        ax_stage.text(value + 0.15, i, f"{value:.3f}", va="center", fontsize=7)

    title = fig.suptitle("", fontsize=12, fontweight="bold")

    def draw(frame_index):
        label, iteration, hpwl, overflow, stage = frames[frame_index]
        movable.set_offsets(positions[frame_index][movable_idx])
        placement_title.set_text(label)
        cursor_h.set_xdata([iteration, iteration])
        cursor_o.set_xdata([iteration, iteration])
        colors = ["#9AA5B1"] * len(stage_names)
        if stage is not None:
            for i in range(stage):
                colors[i] = "#009E73"
            colors[stage] = "#D55E00"
        for bar, color in zip(bars, colors):
            bar.set_color(color)
        metric_row = metrics[min(iteration, len(metrics) - 1)]
        lambda_value = float(metric_row["lambda"])
        legal = " | LEGAL" if stage == len(stage_names) - 1 else ""
        title.set_text(
            f"{benchmark} | weighted-average + BB-Nesterov + electrostatic density | "
            f"exact HPWL {hpwl:.3f}M | overflow {overflow:.3f}% | "
            f"lambda {lambda_value:.3e}{legal}")

    output.parent.mkdir(parents=True, exist_ok=True)
    still.parent.mkdir(parents=True, exist_ok=True)
    writer = PillowWriter(fps=fps)
    with writer.saving(fig, output, dpi=100):
        for frame_index in range(len(frames)):
            repeats = 2 if frame_index >= len(global_paths) else 1
            for _ in range(repeats):
                draw(frame_index)
                writer.grab_frame()
        for _ in range(max(6, 2 * fps)):
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
    args = parser.parse_args()
    render(args.raw_base, args.run_dir, args.output, args.still, args.fps)


if __name__ == "__main__":
    main()
