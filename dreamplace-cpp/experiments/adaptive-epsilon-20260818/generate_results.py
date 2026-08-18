from pathlib import Path
import csv

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / "output" / "adaptive_epsilon_20260818_performance"
OUT = RUN / "visualization"
OUT.mkdir(parents=True, exist_ok=True)
DATASETS = ["adaptec1", "adaptec2", "adaptec3", "adaptec4", "bigblue1", "bigblue2"]

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.18,
})

curves = {}
for dataset in DATASETS:
    path = RUN / dataset / "global_metrics.csv"
    with path.open(newline="", encoding="ascii") as stream:
        rows = list(csv.DictReader(stream))
    curves[dataset] = {
        "iteration": [int(row["iteration"]) for row in rows],
        "hpwl": [float(row["exact_hpwl"]) / 1.0e6 for row in rows],
        "overflow": [100.0 * float(row["overflow"]) for row in rows],
    }

fig, axes = plt.subplots(3, 2, figsize=(11.5, 9.0), constrained_layout=True)
for ax, dataset in zip(axes.flat, DATASETS):
    curve = curves[dataset]
    ax2 = ax.twinx()
    ax.plot(curve["iteration"], curve["hpwl"], color="#0072B2", lw=1.35,
            label="exact HPWL")
    ax2.plot(curve["iteration"], curve["overflow"], color="#D55E00", lw=1.1,
             ls="--", label="overflow")
    ax2.axhspan(6.5, 7.0, color="#009E73", alpha=0.10)
    ax2.axhline(7.0, color="#009E73", lw=0.8)
    ax.set_title(dataset)
    ax.set_xlabel("iteration")
    ax.set_ylabel("HPWL (M)")
    ax2.set_ylabel("overflow (%)")
    ax2.set_ylim(bottom=0)
fig.suptitle("Epsilon-active exact HPWL and electrostatic overflow convergence")
fig.savefig(OUT / "convergence_curves.pdf", bbox_inches="tight")
fig.savefig(OUT / "convergence_curves.png", dpi=240, bbox_inches="tight")
plt.close(fig)

fig, axes = plt.subplots(2, 1, figsize=(10.5, 7.2), constrained_layout=True)
colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
hpwl_lines = []
overflow_lines = []
for dataset, color in zip(DATASETS, colors):
    curve = curves[dataset]
    line_hpwl, = axes[0].plot([], [], lw=1.5, color=color, label=dataset)
    line_overflow, = axes[1].plot([], [], lw=1.5, color=color, label=dataset)
    hpwl_lines.append(line_hpwl)
    overflow_lines.append(line_overflow)
axes[0].set_title("Exact HPWL convergence")
axes[0].set_ylabel("HPWL (M)")
axes[1].set_title("Electrostatic overflow convergence")
axes[1].set_xlabel("iteration")
axes[1].set_ylabel("overflow (%)")
axes[1].axhspan(6.5, 7.0, color="#009E73", alpha=0.10)
axes[1].axhline(7.0, color="#009E73", lw=0.8, ls="--")
axes[0].legend(frameon=False, ncol=3, loc="best")
axes[1].legend(frameon=False, ncol=3, loc="best")
max_iteration = max(curve["iteration"][-1] for curve in curves.values())
all_hpwl = [value for curve in curves.values() for value in curve["hpwl"]]
all_overflow = [value for curve in curves.values() for value in curve["overflow"]]
axes[0].set_xlim(0, max_iteration)
axes[0].set_ylim(min(all_hpwl) * 0.95, max(all_hpwl) * 1.03)
axes[1].set_xlim(0, max_iteration)
axes[1].set_ylim(0, max(15.0, max(all_overflow) * 1.05))

def update(frame):
    for dataset, line_hpwl, line_overflow in zip(DATASETS, hpwl_lines, overflow_lines):
        curve = curves[dataset]
        end = min(frame + 1, len(curve["iteration"]))
        line_hpwl.set_data(curve["iteration"][:end], curve["hpwl"][:end])
        line_overflow.set_data(curve["iteration"][:end], curve["overflow"][:end])
    fig.suptitle(f"Epsilon-active convergence, iteration {frame}")
    return hpwl_lines + overflow_lines

frames = max(len(curve["iteration"]) for curve in curves.values())
animation = FuncAnimation(fig, update, frames=frames, interval=35, blit=False)
animation.save(OUT / "convergence_animation.gif", writer=PillowWriter(fps=24))
plt.close(fig)
print(OUT / "convergence_curves.pdf")
print(OUT / "convergence_curves.png")
print(OUT / "convergence_animation.gif")
