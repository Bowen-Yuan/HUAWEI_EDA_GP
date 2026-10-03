#!/usr/bin/env python3
"""Generate the 3000-step optimizer/lambda ablation figure."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures"

METHODS = [
    "Continuous lambda (best)",
    "Projected HPWL (+/-0.03%)",
    "Projected HPWL (+/-0.20%)",
    "Lambda target 5.99%",
    "Trajectory PI-D",
    "One-sided projection",
    "Adam",
    "Local momentum restart",
]
HPWL = np.array([79.423, 79.484, 79.720, 79.897, 80.526, 84.329, 87.351, 91.105])
OVERFLOW = np.array([5.9985, 5.9704, 5.8010, 5.9998, 6.0000, 2.6910, 5.9999, 5.9868])


def main() -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif"],
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.16,
        "savefig.bbox": "tight",
    })
    y = np.arange(len(METHODS))
    colors = ["#D55E00"] + ["#9AA7B0"] * (len(METHODS) - 1)
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.75), gridspec_kw={"width_ratios": (1.35, 1.0)})

    bars = axes[0].barh(y, HPWL, color=colors, height=0.62, edgecolor="white")
    axes[0].set_yticks(y, METHODS)
    axes[0].invert_yaxis()
    axes[0].set_xlim(78.5, 92.5)
    axes[0].set_xlabel("Best feasible exact HPWL (million)")
    axes[0].set_title("Wirelength quality")
    for bar, value in zip(bars, HPWL):
        axes[0].text(value + 0.12, bar.get_y() + bar.get_height() / 2,
                     f"{value:.3f}", va="center", fontsize=7.5)

    bars = axes[1].barh(y, OVERFLOW, color=colors, height=0.62, edgecolor="white")
    axes[1].axvline(6.0, color="#333333", ls="--", lw=1.0, label="6% constraint")
    axes[1].set_yticks(y, [""] * len(METHODS))
    axes[1].invert_yaxis()
    axes[1].set_xlim(0, 6.45)
    axes[1].set_xlabel("Overflow (%)")
    axes[1].set_title("Constraint outcome")
    axes[1].legend(loc="lower right", frameon=False, fontsize=7.5)
    for bar, value in zip(bars, OVERFLOW):
        axes[1].text(max(0.08, value - 0.08), bar.get_y() + bar.get_height() / 2,
                     f"{value:.3f}", ha="right", va="center", fontsize=7.3,
                     color="white" if value > 1.0 else "#333333")

    fig.suptitle("3000-step exact-HPWL electrostatic placement ablations", fontweight="bold")
    fig.savefig(OUT / "fig_exact_hpwl_ablation.pdf")
    fig.savefig(OUT / "fig_exact_hpwl_ablation.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
