#!/usr/bin/env python3
"""Render ePlace initialization and final NSP placement/density side by side."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "ispd2005" / "adaptec1" / "adaptec1"
FINAL = ROOT / "experiments" / "final_default_eplace" / "adaptec1.nsp.pl"
OUTPUT = ROOT / "visualizations" / "final_default_eplace"


def read_nodes(path: Path) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray]:
    names: list[str] = []
    widths: list[float] = []
    heights: list[float] = []
    fixed: list[bool] = []
    for raw in path.read_text(encoding="ascii").splitlines():
        fields = raw.split()
        if len(fields) < 3 or fields[0] in {"UCLA", "NumNodes", "NumTerminals"}:
            continue
        try:
            widths.append(float(fields[1]))
            heights.append(float(fields[2]))
        except ValueError:
            continue
        names.append(fields[0])
        fixed.append(len(fields) >= 4 and fields[3].startswith("terminal"))
    return names, np.asarray(widths), np.asarray(heights), np.asarray(fixed)


def read_placement(path: Path, index: dict[str, int], count: int) -> np.ndarray:
    positions = np.full((count, 2), np.nan, dtype=float)
    for raw in path.read_text(encoding="ascii").splitlines():
        fields = raw.split()
        if len(fields) < 3 or fields[0] not in index:
            continue
        try:
            positions[index[fields[0]]] = (float(fields[1]), float(fields[2]))
        except ValueError:
            continue
    if np.isnan(positions).any():
        raise ValueError(f"Incomplete placement: {path}")
    return positions


def generate() -> tuple[Path, Path]:
    names, widths, heights, fixed = read_nodes(BENCH.with_suffix(".nodes"))
    index = {name: i for i, name in enumerate(names)}
    initial = read_placement(BENCH.with_suffix(".eplace-ip.pl"), index, len(names))
    final = read_placement(FINAL, index, len(names))
    for positions in (initial, final):
        positions[:, 0] += 0.5 * widths
        positions[:, 1] += 0.5 * heights

    bounds = (459.0, 11151.0, 459.0, 11139.0)
    sample = np.random.default_rng(0).choice(np.flatnonzero(~fixed), 60000, replace=False)
    bins = 96
    extent = (bounds[0], bounds[1], bounds[2], bounds[3])
    maps = []
    for positions in (initial, final):
        area, _, _ = np.histogram2d(
            positions[~fixed, 1], positions[~fixed, 0], bins=bins,
            range=[[bounds[2], bounds[3]], [bounds[0], bounds[1]]],
            weights=widths[~fixed] * heights[~fixed],
        )
        bin_area = ((bounds[1] - bounds[0]) / bins) * ((bounds[3] - bounds[2]) / bins)
        maps.append(area / bin_area)
    vmax = float(np.percentile(np.concatenate([item.ravel() for item in maps]), 99.0))

    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"],
        "font.size": 9.5, "axes.titlesize": 10.5, "axes.titleweight": "bold",
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.dpi": 160, "savefig.dpi": 300, "savefig.bbox": "tight",
    })
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 6.1))
    labels = ("ePlace initialization", "Final alg-electronic placement")
    for column, (positions, density, label) in enumerate(zip((initial, final), maps, labels)):
        axes[0, column].scatter(positions[sample, 0], positions[sample, 1], s=0.25,
                                alpha=0.32, color="#0072B2", rasterized=True)
        axes[0, column].set_title(label)
        axes[0, column].set_xlim(bounds[0], bounds[1])
        axes[0, column].set_ylim(bounds[2], bounds[3])
        axes[0, column].set_aspect("equal")
        axes[0, column].set_xlabel("x")
        axes[0, column].set_ylabel("y")
        image = axes[1, column].imshow(density, origin="lower", extent=extent,
                                       cmap="YlOrRd", vmin=0.0, vmax=vmax,
                                       aspect="equal", rasterized=True)
        axes[1, column].set_title(f"{label}: movable density")
        axes[1, column].set_xlabel("x")
        axes[1, column].set_ylabel("y")
    fig.colorbar(image, ax=axes[1, :], label="Movable area / bin area",
                 shrink=0.82, pad=0.03)
    fig.suptitle("Electrostatic Spreading on ISPD 2005 adaptec1",
                 fontsize=12, fontweight="bold")
    fig.subplots_adjust(left=0.08, right=0.89, bottom=0.07, top=0.92,
                        hspace=0.28, wspace=0.24)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    png = OUTPUT / "placement_comparison.png"
    pdf = OUTPUT / "placement_comparison.pdf"
    fig.savefig(png)
    fig.savefig(pdf)
    plt.close(fig)
    return png, pdf


if __name__ == "__main__":
    print("\n".join(map(str, generate())))
