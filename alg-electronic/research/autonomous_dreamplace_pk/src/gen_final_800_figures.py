#!/usr/bin/env python3
"""Generate final 800-step convergence and DREAMPlace comparison figures."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "experiments" / "runs"
OUTPUT = ROOT / "to_human" / "final_800"

CONFIGS = [
    {
        "name": "adaptec1",
        "label": "adaptec1 (center initialization)",
        "run": "final_a1_800_center_ftrust_visual",
        "reference": 70.3,
        "color": "#0072B2",
    },
    {
        "name": "adaptec2",
        "label": "adaptec2 (ePlace-GP warm-start)",
        "run": "final_a2_800_eplacegp_ftrust_visual",
        "reference": 79.3,
        "color": "#D55E00",
    },
]


def load(run: str) -> dict[str, np.ndarray]:
    path = RUNS / run / "nsp_convergence.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    phase = np.asarray([int(row["phase"]) for row in rows])
    hpwl = np.asarray([float(row["hpwl"]) / 1e6 for row in rows])
    overflow = np.asarray([float(row["avg_overflow"]) for row in rows])
    overflow[phase == 0] = np.nan
    feasible = overflow <= 0.07
    best = np.full_like(hpwl, np.nan)
    incumbent = np.inf
    for i, ok in enumerate(feasible):
        if ok:
            incumbent = min(incumbent, hpwl[i])
        if np.isfinite(incumbent):
            best[i] = incumbent
    return {
        "step": np.arange(len(rows)),
        "phase": phase,
        "hpwl": hpwl,
        "overflow": overflow,
        "best": best,
    }


def style() -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif"],
        "font.size": 9.5,
        "axes.titlesize": 10.5,
        "axes.titleweight": "bold",
        "axes.labelsize": 9.5,
        "legend.fontsize": 8.2,
        "legend.frameon": False,
        "figure.dpi": 160,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.16,
        "lines.linewidth": 1.6,
    })


def convergence() -> tuple[Path, Path]:
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.2), sharex="col")
    for col, cfg in enumerate(CONFIGS):
        data = load(cfg["run"])
        step = data["step"]
        hpwl_ax = axes[0, col]
        overflow_ax = axes[1, col]
        hpwl_ax.plot(step, data["hpwl"], color="#9AA0A6", alpha=0.72,
                     linewidth=1.0, label="Current exact HPWL")
        hpwl_ax.plot(step, data["best"], color=cfg["color"], linewidth=2.1,
                     label="Best feasible HPWL")
        hpwl_ax.axhline(cfg["reference"], color="#009E73", linestyle="--",
                        linewidth=1.2, label="DREAMPlace reference")
        fine_steps = np.flatnonzero(data["phase"] == 3)
        density_start = int(fine_steps[0]) if fine_steps.size else 0
        if density_start > 0:
            hpwl_ax.axvline(density_start, color="#444444", linestyle=":",
                            linewidth=1.0)
        hpwl_ax.set_title(cfg["label"])
        hpwl_ax.set_ylabel("HPWL (million)")
        hpwl_ax.legend(loc="best")

        overflow_ax.plot(step, data["overflow"], color=cfg["color"])
        overflow_ax.axhline(0.07, color="#009E73", linestyle="--",
                            linewidth=1.2, label="7% feasibility limit")
        if density_start > 0:
            overflow_ax.axvline(density_start, color="#444444", linestyle=":",
                                linewidth=1.0, label="Density phase starts")
        overflow_ax.set_xlabel("Total iteration")
        overflow_ax.set_ylabel("Overflow")
        overflow_ax.set_ylim(bottom=0)
        overflow_ax.legend(loc="upper right")

    fig.suptitle("Exact non-smooth HPWL + electrostatic density: 800-step convergence",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    png = OUTPUT / "fig_final_800_convergence.png"
    pdf = OUTPUT / "fig_final_800_convergence.pdf"
    fig.savefig(png)
    fig.savefig(pdf)
    plt.close(fig)
    return png, pdf


def comparison() -> tuple[Path, Path]:
    measured = []
    overflows = []
    for cfg in CONFIGS:
        data = load(cfg["run"])
        feasible = np.flatnonzero(data["overflow"] <= 0.07)
        best = feasible[np.argmin(data["hpwl"][feasible])]
        measured.append(data["hpwl"][best])
        overflows.append(data["overflow"][best])

    x = np.arange(2)
    width = 0.34
    fig, ax = plt.subplots(figsize=(6.75, 3.0))
    refs = [cfg["reference"] for cfg in CONFIGS]
    ax.bar(x - width / 2, refs, width, color="#B0BEC5", label="DREAMPlace reference")
    bars = ax.bar(x + width / 2, measured, width,
                  color=[cfg["color"] for cfg in CONFIGS], label="800-step result")
    for bar, hpwl, overflow in zip(bars, measured, overflows):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.8,
                f"{hpwl:.3f}M\n{overflow * 100:.4f}% ov.",
                ha="center", va="bottom", fontsize=8.2)
    ax.set_xticks(x)
    ax.set_xticklabels(["adaptec1\ncenter", "adaptec2\nwarm-start"])
    ax.set_ylabel("Exact all-net HPWL (million)")
    ax.set_ylim(0, max(refs) * 1.18)
    ax.set_title("Best feasible result within 800 total iterations")
    ax.legend(loc="upper left")
    fig.tight_layout()
    png = OUTPUT / "fig_final_800_baseline_comparison.png"
    pdf = OUTPUT / "fig_final_800_baseline_comparison.pdf"
    fig.savefig(png)
    fig.savefig(pdf)
    plt.close(fig)
    return png, pdf


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    style()
    for path in (*convergence(), *comparison()):
        print(path)


if __name__ == "__main__":
    main()
