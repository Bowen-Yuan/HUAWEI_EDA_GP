#!/usr/bin/env python3
"""Plot exact-HPWL continuation experiments from solver CSV files."""

import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


ROOT = Path(__file__).resolve().parents[2]
RUNS = {
    "adaptec1": ROOT / "output/epsilon_continuation_20260819/a1_adam200_amsgrad100_legalaware",
    "adaptec2": ROOT / "output/epsilon_continuation_20260819/a2_adam200_amsgrad100_legalaware",
}
OUT = ROOT / "output/epsilon_continuation_20260819/figures"

COLORS = {
    "gp": "#0072B2",
    "legal": "#D55E00",
    "overflow": "#009E73",
    "epsilon": "#E69F00",
    "filter": "#CC79A7",
    "fixed_bg": "#DCEAF5",
    "relative_bg": "#FFF0CC",
    "subgradient_bg": "#ECECEC",
}


def read_metrics(path):
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    numeric = {}
    for key in rows[0]:
        if key == "optimizer":
            continue
        numeric[key] = [float(row[key]) for row in rows]
    return numeric


def read_summary(path):
    result = {}
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if "=" not in line:
                continue
            key, value = line.strip().split("=", 1)
            try:
                result[key] = float(value)
            except ValueError:
                result[key] = value
    return result


def stage_spans(iterations, stages):
    spans = []
    start = 0
    current = int(stages[0])
    for index, stage in enumerate(stages[1:], 1):
        stage = int(stage)
        if stage != current:
            spans.append((start, index, current))
            start = index
            current = stage
    spans.append((start, len(iterations) - 1, current))
    return spans


def shade_stages(ax, iterations, stages):
    colors = {1: COLORS["fixed_bg"], 2: COLORS["relative_bg"], 3: COLORS["subgradient_bg"]}
    for left_index, right_index, stage in stage_spans(iterations, stages):
        left = iterations[left_index]
        right = iterations[right_index]
        ax.axvspan(left, right, color=colors[stage], alpha=0.48, linewidth=0)


def main():
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif"],
        "font.size": 9,
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
        "axes.labelsize": 9,
        "legend.fontsize": 8,
        "legend.frameon": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.18,
        "grid.linestyle": "-",
        "lines.linewidth": 1.55,
        "figure.dpi": 160,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    })
    fig, axes = plt.subplots(3, 2, figsize=(7.2, 7.0), sharex="col")
    for column, (dataset, run_dir) in enumerate(RUNS.items()):
        data = read_metrics(run_dir / "global_metrics.csv")
        summary = read_summary(run_dir / "summary.txt")
        x = data["iteration"]
        stage = data["epsilon_stage"]
        for row in range(3):
            shade_stages(axes[row, column], x, stage)

        ax = axes[0, column]
        ax.plot(x, [value / 1e6 for value in data["exact_hpwl"]], color=COLORS["gp"])
        legal_x = []
        legal_y = []
        for iteration, value in zip(x, data["quick_legal_hpwl"]):
            if value > 0:
                legal_x.append(iteration)
                legal_y.append(value / 1e6)
        ax.plot(legal_x, legal_y, color=COLORS["legal"], marker="o", markersize=3.4,
                linestyle="--", label="Quick-legal checkpoint")
        selected_iteration = summary.get("selected_legal_iteration", -1)
        final_legal = summary.get("detailed_hpwl", 0.0) / 1e6
        if selected_iteration >= 0 and final_legal > 0:
            ax.scatter([selected_iteration], [final_legal], marker="*", s=58,
                       color="#000000", zorder=5, label="Selected final legal")
        ax.set_title(dataset)
        ax.set_ylabel("Exact HPWL (M)")

        ax = axes[1, column]
        ax.plot(x, [100.0 * value for value in data["overflow"]], color=COLORS["overflow"])
        ax.axhline(7.0, color="#333333", linestyle="--", linewidth=1.0)
        ax.set_ylabel("Overflow (%)")
        ax.set_ylim(bottom=0)

        ax = axes[2, column]
        ax.plot(x, [100.0 * value for value in data["effective_span_ratio"]],
                color=COLORS["epsilon"], label="Effective span cap")
        filter_values = [100.0 * value if int(s) == 3 else float("nan")
                         for value, s in zip(data["previous_exact_filter_fraction"], stage)]
        ax.plot(x, filter_values, color=COLORS["filter"], alpha=0.75,
                label="Accepted step fraction")
        ax.set_ylabel("Control value (%)")
        ax.set_xlabel("Iteration")
        ax.set_ylim(-1, 105)

    stage_legend = [
        Patch(facecolor=COLORS["fixed_bg"], label="Fixed epsilon"),
        Patch(facecolor=COLORS["relative_bg"], label="Per-net continuation"),
        Patch(facecolor=COLORS["subgradient_bg"], label="Exact subgradient"),
        Line2D([0], [0], color=COLORS["gp"], label="GP exact HPWL"),
        Line2D([0], [0], color=COLORS["legal"], linestyle="--", marker="o",
               markersize=3.4, label="Quick-legal HPWL"),
    ]
    fig.legend(handles=stage_legend, loc="upper center", ncol=3,
               bbox_to_anchor=(0.5, 1.01))
    fig.tight_layout(rect=(0, 0, 1, 0.955), h_pad=1.0, w_pad=1.1)
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "fig_epsilon_continuation_convergence.pdf")
    fig.savefig(OUT / "fig_epsilon_continuation_convergence.png")


if __name__ == "__main__":
    main()
