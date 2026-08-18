#!/usr/bin/env python3
"""Create reproducible plots for the adaptive epsilon ISPD experiment."""
from pathlib import Path
import csv
import math

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "adaptive_epsilon_ispd2005_20260814"
FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(parents=True, exist_ok=True)

DATASETS = ["adaptec1", "adaptec2", "adaptec3", "adaptec4",
            "bigblue1", "bigblue2", "bigblue3_retry3", "bigblue4"]
LABELS = {"bigblue3_retry3": "bigblue3"}
PAPER = {"adaptec1": 73.22, "adaptec2": 82.22, "adaptec3": 193.72,
         "adaptec4": 174.08, "bigblue1": 89.38, "bigblue2": 136.54,
         "bigblue3_retry3": 303.90, "bigblue4": 743.75}

plt.rcParams.update({
    "font.family": "DejaVu Serif", "font.size": 9,
    "axes.titlesize": 10, "axes.labelsize": 9,
    "legend.fontsize": 7.5, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True,
    "grid.alpha": 0.18, "figure.dpi": 160,
})

summary = {}
for name in DATASETS:
    values = {}
    with (OUT / name / "summary.txt").open(encoding="utf-8") as stream:
        for line in stream:
            if "=" in line:
                key, value = line.rstrip().split("=", 1)
                values[key] = value
    summary[name] = values

fig, axes = plt.subplots(4, 2, figsize=(10.0, 12.0), sharex=False)
for ax, name in zip(axes.flat, DATASETS):
    curve_name = "bigblue3" if name == "bigblue3_retry3" else name
    rows = []
    with (OUT / curve_name / "global_metrics.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    x = [int(row["iteration"]) for row in rows]
    hpwl = [float(row["exact_hpwl"]) / 1e6 for row in rows]
    overflow = [100.0 * float(row["overflow"]) for row in rows]
    ax.plot(x, hpwl, color="#0072B2", lw=1.35, label="exact GP HPWL (M)")
    ax.set_title(LABELS.get(name, name))
    ax.set_xlabel("iteration")
    ax.set_ylabel("HPWL (M)")
    ax.grid(True, alpha=0.18)
    ax2 = ax.twinx()
    ax2.plot(x, overflow, color="#D55E00", lw=1.1, ls="--", label="overflow (%)")
    ax2.axhspan(6.5, 7.0, color="#009E73", alpha=0.09)
    ax2.axhline(7.0, color="#009E73", lw=0.7)
    ax2.set_ylabel("overflow (%)")
    ax2.set_ylim(bottom=0)
    if name == "bigblue4":
        ax2.text(0.98, 0.78, "not feasible", transform=ax2.transAxes,
                 ha="right", color="#D55E00", fontsize=8)
fig.suptitle("Exact nonsmooth HPWL and electrostatic overflow convergence", y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.985))
fig.savefig(FIG / "adaptive_epsilon_convergence.pdf", bbox_inches="tight")
fig.savefig(FIG / "adaptive_epsilon_convergence.png", dpi=260, bbox_inches="tight")
plt.close(fig)

labels = [LABELS.get(name, name) for name in DATASETS]
ours = [float(summary[name]["detailed_hpwl"]) / 1e6 for name in DATASETS]
paper = [PAPER[name] for name in DATASETS]
fig, ax = plt.subplots(figsize=(9.2, 4.0))
pos = list(range(len(labels)))
width = 0.36
ax.bar([p - width / 2 for p in pos], paper, width, label="DREAMPlace V100 (paper)", color="#B0BEC5")
ax.bar([p + width / 2 for p in pos], ours, width, label="this C++ flow: final legal", color="#E69F00")
ax.set_xticks(pos, labels)
ax.set_ylabel("HPWL (M), lower is better")
ax.set_title("Final legal HPWL versus DREAMPlace Table II V100 reference")
ax.legend(frameon=False, ncol=2, loc="upper left")
ax.grid(axis="y", alpha=0.18)
fig.tight_layout()
fig.savefig(FIG / "adaptive_epsilon_vs_dreamplace.pdf", bbox_inches="tight")
fig.savefig(FIG / "adaptive_epsilon_vs_dreamplace.png", dpi=260, bbox_inches="tight")
plt.close(fig)

print("wrote", FIG / "adaptive_epsilon_convergence.pdf")
print("wrote", FIG / "adaptive_epsilon_vs_dreamplace.pdf")
for name in DATASETS:
    legal = float(summary[name]["detailed_hpwl"]) / 1e6
    ref = PAPER[name]
    print(name, "legal=%.3fM reference=%.2fM degradation=%.2f%%" %
          (legal, ref, 100.0 * (legal / ref - 1.0)))
