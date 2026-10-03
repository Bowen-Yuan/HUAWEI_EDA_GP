"""Render the new end-to-end trajectory with legacy overflow as primary metric."""
from pathlib import Path
import csv

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
CSV = ROOT.parents[2] / "results_e2e_hpwl57_joint512" / "adaptec1" / "20260721T040215Z" / "iterations.csv"

plt.rcParams.update({"font.family": "DejaVu Serif", "font.size": 9, "axes.titlesize": 10,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.alpha": .18, "figure.dpi": 300, "savefig.dpi": 300,
                     "savefig.bbox": "tight"})


def values(rows, key):
    return np.asarray([float(row[key]) for row in rows], dtype=float)


def main() -> None:
    with CSV.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    iteration = values(rows, "iteration")
    hpwl = values(rows, "hpwl") / 1e6
    legacy = values(rows, "legacy_overflow_percent")
    lam = values(rows, "lambda_dimensionless")
    split = next((i for i, row in enumerate(rows) if row["phase"] == "joint"), len(rows))
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.0), constrained_layout=True)
    axes[0, 0].plot(iteration, hpwl, color="#0072B2")
    axes[0, 0].axvline(iteration[min(split, len(rows)-1)], color="#444", ls="--", lw=1)
    axes[0, 0].set(xlabel="Iteration", ylabel="HPWL (M)", title="HPWL rescue and joint optimization")
    finite = np.isfinite(legacy)
    axes[0, 1].plot(iteration[finite], legacy[finite], color="#D55E00")
    axes[0, 1].axhspan(20, 30, color="#009E73", alpha=.12)
    axes[0, 1].axvline(iteration[min(split, len(rows)-1)], color="#444", ls="--", lw=1)
    axes[0, 1].set(xlabel="Iteration", ylabel="Legacy overflow (%)", title="Legacy overflow trajectory")
    axes[1, 0].plot(iteration, lam, color="#6A3D9A")
    axes[1, 0].set(xlabel="Iteration", ylabel="Lambda", title="Lambda continuation")
    axes[1, 1].plot(legacy[finite], hpwl[finite], color="#009E73")
    axes[1, 1].axvspan(20, 30, color="#009E73", alpha=.12)
    axes[1, 1].set(xlabel="Legacy overflow (%)", ylabel="HPWL (M)", title="Joint trajectory")
    for suffix in ("pdf", "png"):
        fig.savefig(ROOT / f"legacy_convergence.{suffix}")


if __name__ == "__main__":
    main()
