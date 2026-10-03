#!/usr/bin/env python3
"""Generate publication figures for DREAMPlace-style legalization experiments."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "experiments" / "runs"
OUTPUT = ROOT / "to_human" / "dreamplace_legalization_800"

CONFIGS = [
    {
        "name": "adaptec1",
        "run": "final_dp_a1_800_ov050_d32",
        "baseline": 70.3,
        "phase_boundary": 150,
        "color": "#0072B2",
        "init": "center",
    },
    {
        "name": "adaptec2",
        "run": "final_dp_a2_800_ov050_d16",
        "baseline": 79.3,
        "phase_boundary": None,
        "color": "#D55E00",
        "init": "ePlace-GP warm-start",
    },
]


def style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9.2,
        "axes.titlesize": 10.2,
        "axes.titleweight": "bold",
        "axes.labelsize": 9.2,
        "legend.fontsize": 8.0,
        "legend.frameon": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.18,
        "lines.linewidth": 1.45,
        "figure.dpi": 160,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    })


def load_convergence(run: str) -> dict[str, np.ndarray]:
    with (RUNS / run / "nsp_convergence.csv").open(
        newline="", encoding="ascii"
    ) as handle:
        rows = list(csv.DictReader(handle))
    return {
        "step": np.arange(1, len(rows) + 1),
        "phase": np.asarray([int(row["phase"]) for row in rows]),
        "hpwl": np.asarray([float(row["hpwl"]) / 1.0e6 for row in rows]),
        "overflow": np.asarray([float(row["avg_overflow"]) * 100 for row in rows]),
        "lambda": np.maximum(
            np.asarray([float(row["lambda"]) for row in rows]), 1.0e-12
        ),
    }


def save(fig: plt.Figure, stem: str) -> tuple[Path, Path]:
    png = OUTPUT / f"{stem}.png"
    pdf = OUTPUT / f"{stem}.pdf"
    fig.savefig(png)
    fig.savefig(pdf)
    plt.close(fig)
    return png, pdf


def convergence_figure() -> tuple[Path, Path]:
    fig, axes = plt.subplots(3, 2, figsize=(7.2, 6.8), sharex="col")
    for col, cfg in enumerate(CONFIGS):
        data = load_convergence(cfg["run"])
        step = data["step"]
        axes[0, col].plot(step, data["hpwl"], color=cfg["color"])
        axes[0, col].axhline(
            cfg["baseline"], color="#009E73", linestyle="--",
            linewidth=1.15, label="DREAMPlace pre-legal baseline",
        )
        axes[1, col].plot(step, data["overflow"], color=cfg["color"])
        axes[1, col].axhline(
            5.0, color="#3A3A3A", linestyle="--", linewidth=1.05,
            label="5% feasible-trust boundary",
        )
        axes[2, col].plot(step, data["lambda"], color="#CC79A7")
        axes[2, col].set_yscale("log")
        if cfg["phase_boundary"] is not None:
            for axis in axes[:, col]:
                axis.axvline(
                    cfg["phase_boundary"], color="#777777", linestyle=":",
                    linewidth=1.0,
                )
        axes[0, col].set_title(f"{cfg['name']} ({cfg['init']})")
        axes[0, col].legend(loc="best")
        axes[1, col].legend(loc="best")
        axes[2, col].set_xlabel("Total optimizer step")
    axes[0, 0].set_ylabel("Exact HPWL (million)")
    axes[1, 0].set_ylabel("Density overflow (%)")
    axes[2, 0].set_ylabel("Effective density weight")
    for row in range(3):
        axes[row, 1].tick_params(labelleft=True)
    fig.suptitle(
        "800-step exact-HPWL optimization with trajectory lambda control",
        fontsize=11.5, fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return save(fig, "fig_dp_800_convergence")


def tradeoff_figure() -> tuple[Path, Path]:
    tradeoffs = {
        "adaptec1": {
            "pre": [68.183, 69.727, 70.171],
            "legal": [150.517, 146.691, 145.597],
            "overflow": [6.8495, 5.3450, 4.9267],
            "baseline": 70.3,
            "color": "#0072B2",
        },
        "adaptec2": {
            "pre": [73.610, 76.927, 78.782, 82.478],
            "legal": [116.295, 104.144, 99.111, 91.006],
            "overflow": [6.9999, 5.4999, 5.0000, 4.5000],
            "baseline": 79.3,
            "color": "#D55E00",
        },
    }
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.15))
    for ax, (name, values) in zip(axes, tradeoffs.items()):
        pre = np.asarray(values["pre"])
        legal = np.asarray(values["legal"])
        ax.plot(pre, legal, marker="o", markersize=5.2, color=values["color"])
        for x, y, overflow in zip(pre, legal, values["overflow"]):
            ax.annotate(
                f"{overflow:.2f}%", (x, y), xytext=(4, 5),
                textcoords="offset points", fontsize=7.8,
            )
        ax.axvline(
            values["baseline"], color="#009E73", linestyle="--",
            linewidth=1.1, label="DREAMPlace pre-legal HPWL",
        )
        ax.set_title(name)
        ax.set_xlabel("Pre-legal exact HPWL (million)")
        ax.set_ylabel("Greedy/Abacus HPWL (million)")
        ax.legend(loc="best")
    fig.suptitle(
        "Density target versus legalization cost (labels show achieved overflow)",
        fontsize=11.2, fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return save(fig, "fig_dp_density_legalization_tradeoff")


def stage_figure() -> tuple[Path, Path]:
    labels = ["Pre-legal", "Greedy", "Abacus\ncandidate", "Exact reorder"]
    a1 = [70.171, 146.234, 145.597, 130.724]
    a2 = [78.782, 99.111, 99.489, 96.523]
    x = np.arange(len(labels))
    width = 0.35
    fig, ax = plt.subplots(figsize=(7.0, 3.25))
    bars1 = ax.bar(x - width / 2, a1, width, color="#0072B2", label="adaptec1")
    bars2 = ax.bar(x + width / 2, a2, width, color="#D55E00", label="adaptec2")
    bars2[2].set_hatch("///")
    bars2[2].set_alpha(0.52)
    for bars in (bars1, bars2):
        for bar in bars:
            ax.text(
                bar.get_x() + bar.get_width() / 2, bar.get_height() + 1.2,
                f"{bar.get_height():.1f}", ha="center", va="bottom", fontsize=7.5,
            )
    ax.annotate(
        "rejected; Greedy retained", xy=(x[2] + width / 2, a2[2]),
        xytext=(15, 18), textcoords="offset points", fontsize=7.5,
        arrowprops={"arrowstyle": "->", "color": "#555555", "lw": 0.8},
    )
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Exact HPWL (million)")
    ax.set_ylim(0, 165)
    ax.set_title("Legalization and detailed-placement stage ablation")
    ax.legend(loc="upper left", ncol=2)
    fig.tight_layout()
    return save(fig, "fig_dp_legalization_stage_ablation")


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    style()
    for path in (*convergence_figure(), *tradeoff_figure(), *stage_figure()):
        print(path)


if __name__ == "__main__":
    main()
