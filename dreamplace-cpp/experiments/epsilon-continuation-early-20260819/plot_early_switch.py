#!/usr/bin/env python3
"""Compare late, early, and continuous-to-zero epsilon schedules."""

import csv
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]
RUNS = {
    "adaptec1": [
        ("Late -> 15%", ROOT / "output/epsilon_continuation_20260819/a1_adam200_amsgrad100_legalaware"),
        ("Early -> 15%", ROOT / "output/epsilon_continuation_early_20260819/a1_early1100_to15"),
        ("Early -> 0", ROOT / "output/epsilon_continuation_early_20260819/a1_early1100_tozero"),
    ],
    "adaptec2": [
        ("Late -> 15%", ROOT / "output/epsilon_continuation_20260819/a2_adam200_amsgrad100_legalaware"),
        ("Early -> 15%", ROOT / "output/epsilon_continuation_early_20260819/a2_early1300_to15"),
    ],
}
COLORS = {"Late -> 15%": "#0072B2", "Early -> 15%": "#D55E00", "Early -> 0": "#009E73"}
STYLES = {"Late -> 15%": "--", "Early -> 15%": "-", "Early -> 0": "-."}
OUT = ROOT / "output/epsilon_continuation_early_20260819/figures"


def load(path):
    with (path / "global_metrics.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    keys = ("iteration", "exact_hpwl", "overflow", "active_radius")
    return {key: [float(row[key]) for row in rows] for key in keys}


def main():
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"],
        "font.size": 9, "axes.titlesize": 11, "axes.titleweight": "bold",
        "axes.labelsize": 9, "legend.fontsize": 8, "legend.frameon": False,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.alpha": 0.18, "lines.linewidth": 1.55,
        "figure.dpi": 160, "savefig.dpi": 300, "savefig.bbox": "tight",
    })
    fig, axes = plt.subplots(3, 2, figsize=(7.2, 6.6), sharex="col")
    for column, (dataset, runs) in enumerate(RUNS.items()):
        for label, path in runs:
            data = load(path)
            keep = [i for i, step in enumerate(data["iteration"]) if step >= 600]
            x = [data["iteration"][i] for i in keep]
            style = dict(color=COLORS[label], linestyle=STYLES[label], label=label)
            axes[0, column].plot(x, [data["exact_hpwl"][i] / 1e6 for i in keep], **style)
            axes[1, column].plot(x, [100.0 * data["overflow"][i] for i in keep], **style)
            axes[2, column].plot(x, [data["active_radius"][i] for i in keep], **style)
        axes[0, column].set_title(dataset)
        axes[1, column].axhline(7.0, color="#333333", linewidth=1.0, linestyle=":")
        axes[1, column].set_ylim(0, 12)
        axes[2, column].set_xlabel("Iteration")
        axes[2, column].set_ylim(bottom=0)
    axes[0, 0].set_ylabel("Exact HPWL (M)")
    axes[1, 0].set_ylabel("Overflow (%)")
    axes[2, 0].set_ylabel("Epsilon radius (DBU)")
    axes[0, 1].legend(loc="upper right")
    fig.tight_layout(h_pad=0.8, w_pad=1.0)
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "fig_early_switch_comparison.pdf")
    fig.savefig(OUT / "fig_early_switch_comparison.png")


if __name__ == "__main__":
    main()
