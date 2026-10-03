"""Create reproducible summary figures for the end-to-end placement report."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
plt.rcParams.update({"font.family": "DejaVu Serif", "font.size": 9, "axes.titlesize": 10,
                     "axes.labelsize": 9, "legend.fontsize": 8, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.alpha": .18,
                     "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight"})

def main() -> None:
    names = ["Full init", "HPWL rescue", "Joint lambda", "Short axis", "Alternating"]
    hpwl_m = np.array([127.623, 57.490, 53.998, 77.618, 693.380])
    legacy_o = np.array([21.745, 40.694, 38.548, 38.945, 18.620])
    colors = ["#0072B2", "#E69F00", "#D55E00", "#009E73", "#CC79A7"]
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.65), constrained_layout=True)
    ax = axes[0]
    offsets = [(5, 5), (4, 12), (4, -11), (4, 5), (5, 5)]
    for name, hpwl, overflow, color, offset in zip(names, hpwl_m, legacy_o, colors, offsets):
        ax.scatter(hpwl, overflow, s=48, color=color, zorder=3, label=name)
        ax.annotate(name, (hpwl, overflow), xytext=offset, textcoords="offset points", fontsize=7)
    ax.axvspan(55, 60, color="#56B4E9", alpha=.10)
    ax.axhspan(20, 30, color="#009E73", alpha=.10)
    ax.set(xlabel="HPWL (million)", ylabel="Legacy overflow (%)", title="Observed trade-off")
    ax.set_xscale("log"); ax.set_xlim(45, 1200); ax.set_ylim(0, 55)
    ax = axes[1]
    stages = ["Init", "Rescue", "Joint"]
    stage_h = np.array([127.623, 57.490, 53.998]); stage_o = np.array([21.745, 40.694, 38.548])
    x = np.arange(len(stages)); width = .36
    left = ax.bar(x - width / 2, stage_h, width, color="#0072B2", label="HPWL (M)")
    right = ax.twinx()
    bars = right.bar(x + width / 2, stage_o, width, color="#D55E00", label="Legacy overflow (%)")
    ax.set(xticks=x, xticklabels=stages, ylabel="HPWL (million)", title="New end-to-end run")
    right.set_ylabel("Legacy overflow (%)", color="#D55E00"); right.tick_params(axis="y", colors="#D55E00")
    ax.legend([left, bars], ["HPWL (M)", "Legacy overflow (%)"], loc="upper left", fontsize=7)
    for bar in list(left) + list(bars):
        target = ax if bar in left else right
        target.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f"{bar.get_height():.1f}",
                    ha="center", va="bottom", fontsize=7)
    for suffix in ("pdf", "png"):
        fig.savefig(ROOT / f"experiment_summary.{suffix}")

if __name__ == "__main__":
    main()
