#!/usr/bin/env python3
"""Generate convergence and placement-density figures for exact-HPWL runs."""

from pathlib import Path
import csv

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "figures"
BEST_3000 = ROOT / "experiments" / "exact_hpwl_baseline_up115_target0598_3000"
BEST_800 = ROOT / "experiments" / "exact_hpwl_adam_up150_direct_800_nocool"

BLUE = "#0072B2"
ORANGE = "#E69F00"
VERMILION = "#D55E00"
GRAY = "#707070"


plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "DejaVu Serif"],
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.titleweight": "bold",
    "axes.labelsize": 9,
    "legend.fontsize": 8,
    "legend.frameon": False,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.18,
    "lines.linewidth": 1.5,
})


def load_convergence(path):
    rows = []
    with path.open(newline="") as handle:
        for global_iter, row in enumerate(csv.DictReader(handle)):
            rows.append({
                "iter": global_iter,
                "phase": int(row["phase"]),
                "hpwl": float(row["hpwl"]),
                "overflow": float(row["avg_overflow"]),
                "lambda": float(row["lambda"]),
            })
    return {key: np.asarray([row[key] for row in rows]) for key in rows[0]}


def plot_convergence():
    long_run = load_convergence(BEST_3000 / "nsp_convergence.csv")
    short_run = load_convergence(BEST_800 / "nsp_convergence.csv")
    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.55))

    for data, label, color in (
        (long_run, "3000-step heavy-ball", BLUE),
        (short_run, "800-step Adam", ORANGE),
    ):
        axes[0].plot(data["iter"], data["hpwl"] / 1e6, label=label, color=color)
        density = data["phase"] > 0
        axes[1].plot(data["iter"][density], 100.0 * data["overflow"][density],
                     label=label, color=color)
        positive = density & (data["lambda"] > 0.0)
        axes[2].semilogy(data["iter"][positive], data["lambda"][positive],
                         label=label, color=color)

    axes[0].axhline(70.0, color=GRAY, linestyle="--", linewidth=1.0)
    axes[0].set(title="Exact HPWL", xlabel="Global iteration", ylabel="HPWL (million)")
    axes[0].legend(loc="best")
    axes[1].axhline(6.0, color=VERMILION, linestyle="--", linewidth=1.0,
                    label="6% constraint")
    axes[1].set(title="Smooth overflow", xlabel="Global iteration", ylabel="Overflow (%)")
    axes[1].set_ylim(bottom=0.0)
    axes[1].legend(loc="best")
    axes[2].set(title="Electric-density weight", xlabel="Global iteration",
                ylabel=r"Effective $\lambda$")
    axes[2].legend(loc="best")

    fig.tight_layout(w_pad=1.0)
    for extension in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"fig_exact_hpwl_convergence.{extension}")
    plt.close(fig)


def load_pl(path):
    points = []
    with path.open(errors="replace") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 3 or fields[0].startswith("#") or fields[0] == "UCLA":
                continue
            try:
                points.append((float(fields[1]), float(fields[2])))
            except ValueError:
                continue
    return np.asarray(points)


def plot_placement_density():
    initial = load_pl(ROOT / "ispd2005" / "adaptec1" / "adaptec1.eplace-ip.pl")
    final = load_pl(BEST_3000 / "adaptec1.nsp.pl")
    bounds = [[459.0, 11151.0], [459.0, 11139.0]]
    histograms = [np.histogram2d(points[:, 0], points[:, 1], bins=(96, 96),
                                 range=bounds)[0].T for points in (initial, final)]
    vmax = max(float(hist.max()) for hist in histograms)

    fig, axes = plt.subplots(1, 2, figsize=(6.75, 3.0), sharex=True, sharey=True)
    titles = ("ePlace initialization", "Legalized electrostatic result")
    image = None
    for ax, hist, title in zip(axes, histograms, titles):
        image = ax.imshow(hist, origin="lower", extent=(459, 11151, 459, 11139),
                          cmap="viridis", norm=LogNorm(vmin=1.0, vmax=vmax),
                          interpolation="nearest", aspect="equal")
        ax.set_title(title)
        ax.set_xlabel("x")
        ax.grid(False)
    axes[0].set_ylabel("y")
    cbar = fig.colorbar(image, ax=axes, fraction=0.035, pad=0.025)
    cbar.set_label("Cell-center count per bin (log scale)")
    fig.subplots_adjust(left=0.07, right=0.91, bottom=0.14, top=0.88, wspace=0.12)
    for extension in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"fig_placement_density.{extension}")
    plt.close(fig)


if __name__ == "__main__":
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    plot_convergence()
    plot_placement_density()
    print("generated exact-HPWL convergence and placement-density figures")
