#!/usr/bin/env python3
"""Plot frozen state-triggered epsilon convergence without rerunning placement."""
from pathlib import Path
import csv
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "figures"
OUT.mkdir(parents=True, exist_ok=True)
RUNS = {
    "adaptec1": ROOT / "output/state_triggered_epsilon_20260819/ablation/state15/adaptec1/global_metrics.csv",
    "adaptec2": ROOT / "output/state_triggered_epsilon_20260819/ablation/state15/adaptec2/global_metrics.csv",
    "adaptec3": ROOT / "output/state_triggered_epsilon_20260819/frozen15/adaptec3/global_metrics.csv",
    "adaptec4": ROOT / "output/state_triggered_epsilon_20260819/frozen15/adaptec4/global_metrics.csv",
    "bigblue1": ROOT / "output/state_triggered_epsilon_20260819/frozen15/bigblue1/global_metrics.csv",
    "bigblue2": ROOT / "output/state_triggered_epsilon_20260819/frozen15/bigblue2/global_metrics.csv",
}
COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
plt.rcParams.update({
    "font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"],
    "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 10,
    "axes.titleweight": "bold", "legend.fontsize": 7.5,
    "legend.frameon": False, "figure.dpi": 150, "savefig.dpi": 300,
    "savefig.bbox": "tight", "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.alpha": 0.18,
})

def load(path):
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    it = [int(float(r["iteration"])) for r in rows]
    hpwl = [float(r["exact_hpwl"]) / 1e6 for r in rows]
    overflow = [100.0 * float(r["overflow"]) for r in rows]
    stage = [int(float(r["epsilon_stage"])) for r in rows]
    return it, hpwl, overflow, stage

def plot(metric, ylabel, filename, ylim=None):
    fig, axes = plt.subplots(2, 3, figsize=(9.2, 5.2), sharey=False)
    for ax, (name, path), color in zip(axes.flat, RUNS.items(), COLORS):
        it, hpwl, overflow, stage = load(path)
        values = hpwl if metric == "hpwl" else overflow
        ax.plot(it, values, color=color, linewidth=1.35)
        transitions = [it[i] for i in range(1, len(stage)) if stage[i] != stage[i-1]]
        for transition in transitions:
            ax.axvline(transition, color="#555555", linestyle="--", linewidth=0.7, alpha=0.7)
        ax.set_title(name)
        ax.set_xlabel("Iteration")
        ax.set_ylabel(ylabel)
        if ylim is not None:
            ax.set_ylim(*ylim)
        ax.grid(True, axis="both")
    fig.suptitle("State-triggered epsilon-active placement", y=1.01, fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT / f"{filename}.png")
    fig.savefig(OUT / f"{filename}.pdf")
    plt.close(fig)

plot("hpwl", "Exact HPWL (M)", "fig_state_triggered_hpwl")
plot("overflow", "Overflow (%)", "fig_state_triggered_overflow", (0, 105))
print(f"saved figures to {OUT}")
