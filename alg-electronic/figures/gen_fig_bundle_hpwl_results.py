#!/usr/bin/env python3
"""Generate convergence and ablation plots for the exact-HPWL bundle study."""

from pathlib import Path
import csv

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BASELINE = ROOT / "experiments" / "exact_hpwl_baseline_up115_target0598_3000"
BUNDLE = ROOT / "experiments" / "bundle_size2_3000_repeat2"
OUT = ROOT / "figures"

BLUE = "#0072B2"
VERMILION = "#D55E00"
GREEN = "#009E73"
GRAY = "#9AA7B0"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="ascii") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    baseline = read_csv(BASELINE / "nsp_convergence.csv")
    bundle = read_csv(BUNDLE / "nsp_convergence.csv")
    bundle_log = read_csv(BUNDLE / "bundle_strategy.csv")
    x0 = np.arange(1, len(baseline) + 1)
    x1 = np.arange(1, len(bundle) + 1)

    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif"],
        "font.size": 8.5,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.16,
        "legend.frameon": False,
        "savefig.bbox": "tight",
    })
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.1))

    axes[0, 0].plot(x0, [float(r["hpwl"]) / 1e6 for r in baseline],
                    color=BLUE, lw=1.5, label="Heavy-ball baseline")
    axes[0, 0].plot(x1, [float(r["hpwl"]) / 1e6 for r in bundle],
                    color=VERMILION, lw=1.5, label="Two-cut bundle")
    axes[0, 0].set_ylabel("Exact HPWL (million)")
    axes[0, 0].set_title("Wirelength trajectory")
    axes[0, 0].legend(loc="upper right")

    axes[0, 1].plot(x0, [100 * float(r["avg_overflow"]) for r in baseline],
                    color=BLUE, lw=1.5)
    axes[0, 1].plot(x1, [100 * float(r["avg_overflow"]) for r in bundle],
                    color=VERMILION, lw=1.5)
    axes[0, 1].axhline(6.0, color="#444444", ls="--", lw=1.0, label="6% constraint")
    axes[0, 1].set_ylabel("Overflow (%)")
    axes[0, 1].set_title("Density feasibility")
    axes[0, 1].legend(loc="upper right")

    bx = np.array([int(r["global_iteration"]) for r in bundle_log])
    axes[1, 0].plot(bx, [float(r["newest_weight"]) for r in bundle_log],
                    color=GREEN, lw=1.2, label="Newest-cut weight")
    axes[1, 0].plot(bx, [float(r["norm_ratio"]) for r in bundle_log],
                    color="#CC79A7", lw=1.2, label="Aggregate/current norm")
    axes[1, 0].set_ylim(0, 1.05)
    axes[1, 0].set_ylabel("Bundle diagnostic")
    axes[1, 0].set_xlabel("Global iteration")
    axes[1, 0].set_title("Two-cut aggregation")
    axes[1, 0].legend(loc="lower right")

    labels = [
        "No bundle", "6 cuts", "4 cuts", "3 cuts", "2 cuts (run 1)",
        "2 cuts (run 2)", "Cosine mix", "Decayed cosine mix",
    ]
    hpwl = np.array([79.423, 83.318, 76.522, 76.164, 76.120, 75.677, 76.920, 76.146])
    feasible = np.array([True, False, True, True, True, True, True, True])
    order = np.argsort(hpwl)[::-1]
    colors = [VERMILION if labels[i] == "2 cuts (run 2)" else
              (GRAY if feasible[i] else "#D0D5D8") for i in order]
    y = np.arange(len(labels))
    bars = axes[1, 1].barh(y, hpwl[order], color=colors, height=0.62, edgecolor="white")
    axes[1, 1].set_yticks(y, [labels[i] + ("*" if not feasible[i] else "") for i in order])
    axes[1, 1].set_xlim(74.5, 84.5)
    axes[1, 1].set_xlabel("Best feasible HPWL (million)")
    axes[1, 1].set_title("Bundle structure ablation")
    for bar, value in zip(bars, hpwl[order]):
        axes[1, 1].text(value + 0.08, bar.get_y() + bar.get_height() / 2,
                        f"{value:.3f}", va="center", fontsize=7.2)
    axes[1, 1].text(0.99, 0.02, "* infeasible at 1800 steps", transform=axes[1, 1].transAxes,
                    ha="right", fontsize=7, color="#555555")

    for axis in axes.flat[:3]:
        for boundary in (150, 250, 400):
            axis.axvline(boundary, color="#777777", lw=0.6, alpha=0.35)
    axes[0, 0].set_xlabel("Global iteration")
    axes[0, 1].set_xlabel("Global iteration")
    fig.suptitle("Exact nonsmooth HPWL: proximal bundle optimization", fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT / "fig_bundle_hpwl_results.pdf")
    fig.savefig(OUT / "fig_bundle_hpwl_results.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
