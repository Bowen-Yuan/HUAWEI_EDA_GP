#!/usr/bin/env python3
"""Generate reproducible figures for the progressive legalization report."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
FIGURES = Path(__file__).resolve().parent

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "DejaVu Serif"],
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "legend.fontsize": 7.5,
    "legend.frameon": False,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.18,
    "grid.linestyle": "-",
    "lines.linewidth": 1.6,
})

BLUE = "#0072B2"
GREEN = "#009E73"
ORANGE = "#E69F00"
VERMILION = "#D55E00"
GRAY = "#8C8C8C"
LIGHT_GRAY = "#C8CDD2"


def read_csv(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    columns: dict[str, np.ndarray] = {}
    for key in rows[0]:
        try:
            columns[key] = np.asarray([float(row[key]) for row in rows])
        except ValueError:
            columns[key] = np.asarray([row[key] for row in rows], dtype=object)
    return columns


def read_summary(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            result[key] = value
    return result


def shade_stages(ax: plt.Axes) -> None:
    ax.axvspan(0, 100, color=LIGHT_GRAY, alpha=0.22, linewidth=0)
    ax.axvspan(100, 575, color="#56B4E9", alpha=0.10, linewidth=0)
    ax.axvspan(575, 1000, color="#009E73", alpha=0.08, linewidth=0)
    ax.axvline(783, color=VERMILION, linestyle="--", linewidth=1.0)


def convergence_figure(data: dict[str, np.ndarray]) -> None:
    iteration = data["iteration"]
    fig, axes = plt.subplots(3, 1, figsize=(6.75, 5.8), sharex=True)
    for ax in axes:
        shade_stages(ax)

    axes[0].plot(iteration, data["exact_hpwl"] / 1e6, color=BLUE)
    axes[0].axhspan(75, 80, color=ORANGE, alpha=0.14, linewidth=0)
    axes[0].set_ylabel("Exact HPWL (M)")
    axes[0].set_title("Progressive exact-HPWL trajectory (adaptec1, 512x512)")

    axes[1].plot(iteration, 100 * data["overflow"], color=GREEN)
    axes[1].axhspan(7, 8, color=ORANGE, alpha=0.18, linewidth=0)
    axes[1].set_ylabel("Overflow (%)")

    positive_lambda = np.maximum(data["lambda_effective"], 1e-14)
    axes[2].semilogy(iteration, positive_lambda, color=VERMILION)
    axes[2].set_ylabel("Effective lambda")
    axes[2].set_xlabel("Global-placement iteration")

    axes[0].text(35, axes[0].get_ylim()[1] - 2, "HPWL only", ha="center", fontsize=8)
    axes[0].text(330, axes[0].get_ylim()[1] - 2, "density spreading", ha="center", fontsize=8)
    axes[0].text(680, axes[0].get_ylim()[1] - 2, "progressive legality", ha="center", fontsize=8)
    axes[0].text(790, axes[0].get_ylim()[0] + 2, "AMSGrad reset", color=VERMILION, fontsize=7)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_progressive_convergence.pdf")
    fig.savefig(FIGURES / "fig_progressive_convergence.png")
    plt.close(fig)


def legality_proxy_figure(data: dict[str, np.ndarray]) -> None:
    mask = data["iteration"] >= 575
    iteration = data["iteration"][mask]
    fig, axes = plt.subplots(1, 3, figsize=(6.75, 2.25))
    axes[0].plot(iteration, data["macro_overlap_cells"][mask], color=VERMILION)
    axes[0].set_ylabel("Overlapping cells")
    axes[0].set_xlabel("Iteration")
    axes[0].set_title("Fixed-macro overlap")

    axes[1].plot(iteration, data["row_distance"][mask], color=BLUE)
    axes[1].set_ylabel("Mean row distance")
    axes[1].set_xlabel("Iteration")
    axes[1].set_title("Continuous row proxy")

    axes[2].plot(iteration, 100 * data["segment_overflow"][mask], color=GREEN)
    axes[2].set_ylabel("Segment overflow (%)")
    axes[2].set_xlabel("Iteration")
    axes[2].set_title("Obstacle-free capacity proxy")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_legality_proxies.pdf")
    fig.savefig(FIGURES / "fig_legality_proxies.png")
    plt.close(fig)


def ablation_figure() -> None:
    runs = [
        ("Old exact\n800", "exact_adam_512_800_refine_bundle_dp"),
        ("Pure BandDual\n800", "progressive_a1_512_800_full"),
        ("Hybrid\n800", "progressive_a1_512_800_hybrid"),
        ("Hybrid+refine\n1000", "progressive_a1_512_1000_refine"),
        ("+ sampling/bundle\n1000", "progressive_a1_512_1000_bundle"),
    ]
    gp_values = []
    legal_values = []
    for _, directory in runs:
        summary = read_summary(ROOT / "output" / directory / "summary.txt")
        gp_values.append(float(summary["gp_hpwl"]) / 1e6)
        legal_values.append(float(summary["detailed_hpwl"]) / 1e6)

    x = np.arange(len(runs))
    width = 0.36
    fig, ax = plt.subplots(figsize=(6.75, 2.8))
    gp = ax.bar(x - width / 2, gp_values, width, color=BLUE, label="Pre-legal GP")
    legal = ax.bar(x + width / 2, legal_values, width, color=ORANGE, label="Final legal")
    ax.axhspan(75, 80, color=GREEN, alpha=0.08, linewidth=0)
    ax.set_ylabel("Exact HPWL (M)")
    ax.set_xticks(x)
    ax.set_xticklabels([label for label, _ in runs])
    ax.set_ylim(70, max(legal_values) + 5)
    ax.legend(ncol=2, loc="upper right")
    ax.set_title("Ablation of density control and feasible-band refinement")
    for bars in (gp, legal):
        for bar in bars:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.35,
                    f"{bar.get_height():.2f}", ha="center", va="bottom", fontsize=7)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_ablation.pdf")
    fig.savefig(FIGURES / "fig_ablation.png")
    plt.close(fig)


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    data = read_csv(
        ROOT / "output" / "progressive_a1_512_1000_refine" /
        "global_metrics.csv"
    )
    convergence_figure(data)
    legality_proxy_figure(data)
    ablation_figure()


if __name__ == "__main__":
    main()
