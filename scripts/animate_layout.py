#!/usr/bin/env python3
"""Render GP snapshots as a layout GIF without changing placement data."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


def split_fields(line: str):
    return line.split("#", 1)[0].split()


def read_nodes(path: Path):
    nodes = {}
    for line in path.read_text(encoding="ascii", errors="ignore").splitlines():
        fields = split_fields(line)
        if len(fields) < 3 or fields[0] in {"UCLA", "NumNodes", "NumTerminals"}:
            continue
        try:
            nodes[fields[0]] = (float(fields[1]), float(fields[2]),
                                len(fields) >= 4 and fields[3].startswith("terminal"),
                                len(fields) >= 4 and fields[3] == "terminal_NI")
        except ValueError:
            pass
    return nodes


def read_rows(path: Path):
    rows = []
    y = height = origin = spacing = None
    for line in path.read_text(encoding="ascii", errors="ignore").splitlines():
        if "Coordinate" in line:
            y = float(line.split(":", 1)[1])
        elif "Height" in line:
            height = float(line.split(":", 1)[1])
        elif "Sitespacing" in line:
            spacing = float(line.split(":", 1)[1])
        elif "SubrowOrigin" in line:
            fields = line.split()
            origin = float(fields[2])
            rows.append((origin, y, int(fields[5]) * spacing, height))
    return (min(x for x, _, _, _ in rows), min(y for _, y, _, _ in rows),
            max(x + w for x, _, w, _ in rows), max(y + h for _, y, _, h in rows))


def read_placement(path: Path, nodes, sampled):
    positions = {}
    for line in path.read_text(encoding="ascii", errors="ignore").splitlines():
        fields = split_fields(line)
        if len(fields) < 3 or fields[0] not in nodes:
            continue
        name = fields[0]
        width, height, fixed, fixed_ni = nodes[name]
        if fixed or name in sampled:
            positions[name] = (float(fields[1]), float(fields[2]), width, height,
                               fixed, fixed_ni)
    return positions


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark", required=True, type=Path,
                        help="Bookshelf base path without extension")
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--fps", type=int, default=12)
    parser.add_argument("--max-cells", type=int, default=50000)
    parser.add_argument("--append-placement", type=Path,
                        help="Optional additional placement used as one final frame")
    parser.add_argument("--append-iteration", type=int,
                        help="Iteration label for --append-placement")
    parser.add_argument("--append-hpwl", type=float,
                        help="Exact HPWL for the optional final frame")
    parser.add_argument("--append-overflow", type=float,
                        help="Exact overflow fraction for the optional final frame")
    args = parser.parse_args()
    if args.append_placement is not None and (args.append_iteration is None
                                               or args.append_hpwl is None
                                               or args.append_overflow is None):
        raise ValueError("append placement requires iteration, HPWL, and overflow")

    nodes = read_nodes(args.benchmark.with_suffix(".nodes"))
    movable = sorted(name for name, item in nodes.items() if not item[2])
    stride = max(1, len(movable) // args.max_cells)
    sampled = set(movable[::stride])
    bounds = read_rows(args.benchmark.with_suffix(".scl"))
    def snapshot_iteration(path):
        match = re.search(r"(?:global|iter)_(\d+)$", path.stem)
        return int(match.group(1)) if match else 10**12

    snapshots = list((args.run_dir / "snapshots").glob("global_*.pl"))
    snapshots.extend((args.run_dir / "snapshots").glob("iter_*.pl"))
    snapshots = sorted(set(snapshots), key=snapshot_iteration)
    if (args.run_dir / "global.pl").is_file():
        snapshots.append(args.run_dir / "global.pl")
    if not snapshots:
        raise ValueError("run contains no snapshots or global.pl")
    if args.append_placement is not None:
        if args.append_iteration is None:
            raise ValueError("--append-iteration is required with --append-placement")
        snapshots.append(args.append_placement)
    frames = [read_placement(path, nodes, sampled) for path in snapshots]

    metrics = {}
    metrics_path = args.run_dir / "global_metrics.csv"
    if metrics_path.is_file():
        with metrics_path.open(newline="", encoding="ascii") as stream:
            metrics = {int(row["iteration"]): row for row in csv.DictReader(stream)}
    homotopy_path = args.run_dir / "homotopy_metrics.csv"
    if homotopy_path.is_file():
        with homotopy_path.open(newline="", encoding="ascii") as stream:
            metrics = {int(row["iteration"]): row for row in csv.DictReader(stream)}

    fig, ax = plt.subplots(figsize=(9.0, 7.0), constrained_layout=True)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(bounds[0], bounds[2])
    ax.set_ylim(bounds[1], bounds[3])
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    scatter = ax.scatter([], [], s=0.35, color="#0072B2", alpha=0.55,
                         linewidths=0)
    fixed_patches = [
        Rectangle((x, y), width, height)
        for x, y, width, height, fixed, fixed_ni in frames[0].values()
        if fixed and not fixed_ni
    ]
    fixed_collection = PatchCollection(fixed_patches, facecolor="#D55E00",
                                       edgecolor="#8B1A1A", alpha=0.35,
                                       linewidth=0.35)
    ax.add_collection(fixed_collection)

    def update(frame_index):
        frame = frames[frame_index]
        centers = [(x + 0.5 * width, y + 0.5 * height)
                   for x, y, width, height, fixed, _ in frame.values() if not fixed]
        scatter.set_offsets(centers)
        path = snapshots[frame_index]
        if frame_index == len(snapshots) - 1 and args.append_placement is not None:
            ax.set_title(f"Iteration {args.append_iteration} | HPWL {args.append_hpwl / 1e6:.3f}M | "
                         f"overflow {100 * args.append_overflow:.2f}%")
        elif snapshot_iteration(path) < 10**12:
            iteration = snapshot_iteration(path)
            row = metrics.get(iteration)
            if row:
                overflow_key = "overflow" if "overflow" in row else "exact_overflow"
                ax.set_title(f"Iteration {iteration} | HPWL {float(row['exact_hpwl']) / 1e6:.3f}M | "
                             f"overflow {100 * float(row[overflow_key]):.2f}%")
            else:
                ax.set_title(f"Iteration {iteration}")
        else:
            ax.set_title("Selected global placement")
        return (scatter,)

    animation = FuncAnimation(fig, update, frames=len(frames), interval=1000 / args.fps,
                              blit=False)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    animation.save(args.output, writer=PillowWriter(fps=args.fps))


if __name__ == "__main__":
    main()
