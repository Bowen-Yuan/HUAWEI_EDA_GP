#!/usr/bin/env python3
"""Render visualize-module snapshots as a density/convergence animation."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter
from matplotlib.colors import LogNorm
import numpy as np


COLORS = {
    "hpwl": "#0072B2",
    "overflow": "#D55E00",
    "lambda": "#009E73",
    "bundle_weight": "#CC79A7",
    "bundle_norm": "#E69F00",
    "cursor": "#3A3A3A",
}


def read_metadata(path: Path) -> tuple[int, tuple[float, float, float, float]]:
    lines = path.read_text(encoding="ascii").splitlines()
    nodes = int(lines[0].split()[1])
    bounds = tuple(float(value) for value in lines[1].split()[1:])
    return nodes, bounds  # type: ignore[return-value]


def read_metrics(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="ascii") as stream:
        return list(csv.DictReader(stream))


def read_convergence(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline="", encoding="ascii") as stream:
        rows = list(csv.DictReader(stream))
    return {
        "iteration": np.arange(1, len(rows) + 1, dtype=float),
        "hpwl": np.asarray([float(row["hpwl"]) for row in rows]),
        "overflow": np.asarray([float(row["avg_overflow"]) for row in rows]),
        "lambda": np.asarray([float(row["lambda"]) for row in rows]),
    }


def read_bundle(path: Path) -> dict[str, np.ndarray] | None:
    if not path.exists():
        return None
    with path.open(newline="", encoding="ascii") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        return None
    return {
        "iteration": np.asarray([float(row["global_iteration"]) for row in rows]),
        "newest_weight": np.asarray([float(row["newest_weight"]) for row in rows]),
        "norm_ratio": np.asarray([float(row["norm_ratio"]) for row in rows]),
    }


def render(
    input_dir: Path,
    output_gif: Path,
    output_still: Path,
    fps: int,
    overflow_limit: float,
    benchmark: str,
) -> None:
    nodes, bounds = read_metadata(input_dir / "metadata.txt")
    xl, yl, xh, yh = bounds
    metrics = read_metrics(input_dir / "metrics.csv")
    convergence = read_convergence(input_dir / "nsp_convergence.csv")
    bundle = read_bundle(input_dir / "bundle_strategy.csv")
    snapshots = [input_dir / f"snapshot_{int(row['iteration']):06d}.bin" for row in metrics]
    if not all(path.exists() for path in snapshots):
        missing = [path.name for path in snapshots if not path.exists()]
        raise FileNotFoundError(f"Missing snapshots: {missing[:5]}")

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.18,
        "figure.facecolor": "white",
    })
    plot_rows = 4 if bundle is not None else 3
    fig_height = 8.4 if bundle is not None else 7.2
    fig = plt.figure(figsize=(13.2, fig_height), dpi=100, constrained_layout=True)
    grid = fig.add_gridspec(plot_rows, 2, width_ratios=(1.38, 1.0))
    ax_density = fig.add_subplot(grid[:, 0])
    ax_hpwl = fig.add_subplot(grid[0, 1])
    ax_overflow = fig.add_subplot(grid[1, 1], sharex=ax_hpwl)
    ax_lambda = fig.add_subplot(grid[2, 1], sharex=ax_hpwl)
    ax_bundle = (fig.add_subplot(grid[3, 1], sharex=ax_hpwl)
                 if bundle is not None else None)

    iterations = convergence["iteration"]
    hpwl_m = convergence["hpwl"] / 1.0e6
    overflow_pct = convergence["overflow"] * 100.0
    effective_lambda = np.maximum(convergence["lambda"], 1.0e-12)
    ax_hpwl.plot(iterations, hpwl_m, color=COLORS["hpwl"], lw=1.5)
    ax_overflow.plot(iterations, overflow_pct, color=COLORS["overflow"], lw=1.5)
    overflow_limit_pct = overflow_limit * 100.0
    ax_overflow.axhline(
        overflow_limit_pct,
        color="#555555",
        ls="--",
        lw=1.0,
        label=f"{overflow_limit_pct:g}% reference limit",
    )
    ax_lambda.plot(iterations, effective_lambda, color=COLORS["lambda"], lw=1.5)
    ax_lambda.set_yscale("log")
    ax_hpwl.set_ylabel("Exact HPWL (M)")
    ax_overflow.set_ylabel("Overflow (%)")
    ax_lambda.set_ylabel("Effective lambda")
    if ax_bundle is None:
        ax_lambda.set_xlabel("Global iteration")
    ax_hpwl.set_title("Convergence history")
    ax_overflow.legend(loc="upper right", frameon=False, fontsize=8)
    axes = [ax_hpwl, ax_overflow, ax_lambda]
    if ax_bundle is not None and bundle is not None:
        ax_bundle.plot(
            bundle["iteration"], bundle["newest_weight"],
            color=COLORS["bundle_weight"], lw=1.25, label="Newest-cut weight",
        )
        ax_bundle.plot(
            bundle["iteration"], bundle["norm_ratio"],
            color=COLORS["bundle_norm"], lw=1.25, label="Aggregate/current norm",
        )
        ax_bundle.set_ylabel("Bundle state")
        ax_bundle.set_xlabel("Global iteration")
        ax_bundle.set_ylim(-0.03, 1.08)
        ax_bundle.legend(loc="upper right", frameon=False, fontsize=7, ncol=2)
        axes.append(ax_bundle)

    for axis in axes[:-1]:
        axis.tick_params(labelbottom=False)
    for axis in axes:
        axis.set_xlim(0, len(iterations))

    cursors = [axis.axvline(0, color=COLORS["cursor"], lw=1.1, alpha=0.75)
               for axis in axes]
    bins = (180, 180)

    def draw_frame(frame_index: int) -> None:
        row = metrics[frame_index]
        iteration = int(row["iteration"])
        data = np.fromfile(snapshots[frame_index], dtype=np.float32)
        if data.size != nodes * 2:
            raise ValueError(f"Unexpected size in {snapshots[frame_index]}")
        positions = data.reshape(nodes, 2)
        hist, _, _ = np.histogram2d(
            positions[:, 0], positions[:, 1], bins=bins,
            range=((xl, xh), (yl, yh)),
        )
        ax_density.clear()
        positive = hist[hist > 0]
        vmax = max(2.0, float(np.percentile(positive, 99.5))) if positive.size else 2.0
        ax_density.imshow(
            hist.T, origin="lower", extent=(xl, xh, yl, yh),
            cmap="magma", norm=LogNorm(vmin=1.0, vmax=vmax), aspect="auto",
            interpolation="nearest",
        )
        ax_density.set_title(f"Cell density - iteration {iteration} ({row['phase']})")
        ax_density.set_xlabel("x")
        ax_density.set_ylabel("y")
        ax_density.grid(False)
        for cursor in cursors:
            cursor.set_xdata([iteration, iteration])
        hpwl = float(row["hpwl"]) / 1.0e6
        overflow = float(row["overflow"]) * 100.0
        lam = float(row["lambda"])
        method = "Bundle HPWL + electrostatic density" if bundle is not None else \
                 "Exact nonsmooth HPWL + electrostatic density"
        benchmark_prefix = f"{benchmark} | " if benchmark else ""
        fig.suptitle(
            f"{benchmark_prefix}{method} | "
            f"HPWL {hpwl:.3f}M | overflow {overflow:.3f}% | lambda {lam:.3e}",
            fontsize=12, fontweight="bold",
        )

    writer = PillowWriter(fps=fps)
    with writer.saving(fig, output_gif, dpi=100):
        for frame_index in range(len(snapshots)):
            draw_frame(frame_index)
            writer.grab_frame()
    draw_frame(len(snapshots) - 1)
    fig.savefig(output_still, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--still", type=Path, required=True)
    parser.add_argument("--fps", type=int, default=10)
    parser.add_argument("--overflow-limit", type=float, default=0.06)
    parser.add_argument("--benchmark", default="")
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.still.parent.mkdir(parents=True, exist_ok=True)
    render(
        args.input_dir,
        args.output,
        args.still,
        args.fps,
        args.overflow_limit,
        args.benchmark,
    )


if __name__ == "__main__":
    main()
