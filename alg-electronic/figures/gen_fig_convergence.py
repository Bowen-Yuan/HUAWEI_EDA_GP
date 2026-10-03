#!/usr/bin/env python3
"""Generate publication-quality convergence curves for alg-electronic."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "experiments" / "final_default_eplace" / "nsp_convergence.csv"
OUTPUT_DIR = ROOT / "visualizations" / "final_default_eplace"

PHASES = {
    0: ("HPWL only", "#56B4E9"),
    1: ("Coarse", "#E69F00"),
    2: ("Medium", "#009E73"),
    3: ("Fine", "#CC79A7"),
}


def load_metrics(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"No convergence records in {path}")

    data: dict[str, np.ndarray] = {
        "step": np.arange(len(rows), dtype=float),
        "phase": np.asarray([int(row["phase"]) for row in rows]),
    }
    for key in ("hpwl", "objective", "avg_overflow", "dpen", "lambda", "g_rms", "lr"):
        data[key] = np.asarray([float(row[key]) for row in rows], dtype=float)

    # Phase 0 does not evaluate density. Its logged zeros mean "not measured".
    phase_zero = data["phase"] == 0
    data["avg_overflow"][phase_zero] = np.nan
    data["dpen"][phase_zero] = np.nan
    data["lambda"][data["lambda"] <= 0.0] = np.nan
    return data


def add_phase_regions(axes: np.ndarray, phase: np.ndarray) -> None:
    changes = np.flatnonzero(np.diff(phase)) + 1
    starts = np.r_[0, changes]
    ends = np.r_[changes, len(phase)]
    for ax in axes.flat:
        for start, end in zip(starts, ends):
            phase_id = int(phase[start])
            ax.axvspan(start, end - 1, color=PHASES[phase_id][1], alpha=0.055, zorder=0)
        for boundary in changes:
            ax.axvline(boundary, color="#777777", linewidth=0.8,
                       linestyle=(0, (3, 3)), alpha=0.65, zorder=1)


def generate() -> tuple[Path, Path]:
    data = load_metrics(INPUT)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif"],
        "font.size": 9.5,
        "axes.titlesize": 10.5,
        "axes.titleweight": "bold",
        "axes.labelsize": 9.5,
        "legend.fontsize": 8.5,
        "legend.frameon": False,
        "figure.dpi": 160,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.18,
        "grid.linestyle": "-",
        "lines.linewidth": 1.7,
    })

    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.3), sharex=True)
    add_phase_regions(axes, data["phase"])
    step = data["step"]

    axes[0, 0].plot(step, data["hpwl"] / 1e6, color="#0072B2")
    axes[0, 0].axhline(70.0, color="#555555", linestyle="--", linewidth=1.0,
                       label="Constraint = 70M")
    axes[0, 0].set_ylabel("HPWL (million)")
    axes[0, 0].set_title("Wirelength")

    axes[0, 1].plot(step, data["avg_overflow"], color="#D55E00")
    axes[0, 1].axhline(0.06, color="#555555", linestyle="--", linewidth=1.0,
                       label="Constraint = 0.06")
    axes[0, 1].set_ylabel("Average overflow")
    axes[0, 1].set_title("Baseline-Compatible Overflow")
    axes[0, 1].legend(loc="upper right")

    axes[1, 0].plot(step, data["dpen"] / 1e9, color="#009E73")
    axes[1, 0].set_yscale("log")
    axes[1, 0].set_ylabel("Electric energy (billion, log scale)")
    axes[1, 0].set_xlabel("Global iteration")
    axes[1, 0].set_title("Electrostatic Density Energy")

    axes[1, 1].plot(step, data["lambda"], color="#CC79A7")
    axes[1, 1].set_yscale("log")
    axes[1, 1].set_ylabel("Effective lambda (log scale)")
    axes[1, 1].set_xlabel("Global iteration")
    axes[1, 1].set_title("Density Weight")

    phase_handles = [
        Patch(facecolor=color, alpha=0.20, edgecolor="none", label=label)
        for label, color in PHASES.values()
    ]
    fig.legend(handles=phase_handles, loc="upper center", ncol=4,
               bbox_to_anchor=(0.5, 0.955), columnspacing=1.4, handlelength=1.5)
    feasible = np.flatnonzero((data["phase"] == 3) &
                              (data["avg_overflow"] <= 0.06))
    if feasible.size:
        best = feasible[np.argmin(data["hpwl"][feasible])]
        axes[0, 0].scatter([step[best]], [data["hpwl"][best] / 1e6],
                           color="#D55E00", s=28, zorder=5,
                           label="Best feasible")
        axes[0, 0].legend(loc="upper right")

    fig.suptitle("Final alg-electronic Run on ISPD 2005 adaptec1",
                 fontsize=12, fontweight="bold", y=0.995)
    fig.subplots_adjust(left=0.11, right=0.98, bottom=0.10, top=0.86,
                        hspace=0.36, wspace=0.30)

    png = OUTPUT_DIR / "convergence_curves.png"
    pdf = OUTPUT_DIR / "convergence_curves.pdf"
    fig.savefig(png)
    fig.savefig(pdf)
    plt.close(fig)
    return png, pdf


if __name__ == "__main__":
    png_path, pdf_path = generate()
    print(png_path)
    print(pdf_path)
