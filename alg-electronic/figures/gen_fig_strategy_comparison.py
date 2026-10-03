#!/usr/bin/env python3
"""Compare the best feasible point from each lambda/optimizer strategy."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "visualizations" / "strategy_comparison"
RUNS = [
    ("2.0x update", "eplace_default_3000"),
    ("1.5x update", "eplace_up150_3000"),
    ("Target 5.99%", "eplace_up150_target0599_3000"),
    ("P0 restore, 5.98%", "eplace_p0restore_up150_target0598_3000"),
    ("P0 restore, 1.6x", "eplace_p0restore_up160_target0599_3000"),
    ("Final: P0 restore, 1.5x", "final_default_eplace"),
]


def best_feasible(directory: str) -> tuple[float, float]:
    path = ROOT / "experiments" / directory / "nsp_convergence.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = [row for row in csv.DictReader(handle)
                if int(row["phase"]) == 3 and float(row["avg_overflow"]) <= 0.06]
    if not rows:
        raise ValueError(f"No feasible point in {path}")
    row = min(rows, key=lambda item: float(item["hpwl"]))
    return float(row["hpwl"]) / 1e6, 100.0 * float(row["avg_overflow"])


def generate() -> tuple[Path, Path]:
    values = [(label, *best_feasible(directory)) for label, directory in RUNS]
    labels = [item[0] for item in values]
    hpwl = [item[1] for item in values]
    overflow = [item[2] for item in values]

    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif"],
        "font.size": 9.5,
        "axes.titlesize": 10.5,
        "axes.titleweight": "bold",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.16,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    })
    colors = ["#B0BEC5"] * len(labels)
    colors[-1] = "#D55E00"

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(7.0, 3.3))
    y = range(len(labels))
    bars = ax0.barh(y, hpwl, color=colors, height=0.62)
    ax0.set_yticks(list(y), labels)
    ax0.invert_yaxis()
    ax0.axvline(70.0, color="#0072B2", linestyle="--", linewidth=1.1)
    ax0.text(70.35, 0.15, "70M baseline", color="#0072B2", rotation=90,
             va="top", fontsize=8)
    ax0.set_xlabel("Best feasible HPWL (million)")
    ax0.set_title("Wirelength Progress")
    ax0.set_xlim(min(hpwl) - 0.15, max(hpwl) + 0.35)
    for bar, value in zip(bars, hpwl):
        ax0.text(value + 0.35, bar.get_y() + bar.get_height() / 2,
                 f"{value:.2f}", va="center", fontsize=8)

    points = ax1.scatter(hpwl, overflow, c=colors, s=55,
                         edgecolors="white", linewidths=0.7, zorder=3)
    del points
    ax1.axhline(6.0, color="#555555", linestyle="--", linewidth=1.1)
    ax1.axvline(70.0, color="#0072B2", linestyle="--", linewidth=1.1)
    for label, x, y_value in values:
        if label in {"2.0x update", "Final: P0 restore, 1.5x"}:
            ax1.annotate(label, (x, y_value), xytext=(4, 5),
                         textcoords="offset points", fontsize=8)
    ax1.set_xlabel("HPWL (million)")
    ax1.set_ylabel("Overflow (%)")
    ax1.set_title("Feasible Trade-off")
    ax1.set_ylim(min(overflow) - 0.03, 6.04)

    fig.suptitle("Lambda-Strategy Ablation on adaptec1",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    OUTPUT.mkdir(parents=True, exist_ok=True)
    png = OUTPUT / "strategy_comparison.png"
    pdf = OUTPUT / "strategy_comparison.pdf"
    fig.savefig(png)
    fig.savefig(pdf)
    plt.close(fig)
    return png, pdf


if __name__ == "__main__":
    print("\n".join(map(str, generate())))
