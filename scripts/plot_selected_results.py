#!/usr/bin/env python3
"""Publication plots for selected exact nonsmooth placement runs.

The plotting code intentionally reads only exact audit columns.  Homotopy
terms are shown for diagnosis, but no smoothed density value is used as the
reported overflow.
"""

import csv
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RUNS = {
    "H160 a1 strict feasible": ROOT / "output/h160_a1_relaxed_cap_exchange",
    "H165 a2 strict feasible": ROOT / "output/h165_a2_shape_exchange",
    "H171 a4 strict feasible": ROOT / "output/h171_a4_breakpoint_recovery",
}
COLORS = ["#0072B2", "#D55E00", "#009E73"]


def read_csv(path):
    with path.open(newline="", encoding="ascii") as stream:
        return list(csv.DictReader(stream))


def load(label, run):
    """Load exact metrics while accepting all experiment-era filenames."""
    parts = []
    for filename, overflow_key in (("homotopy_metrics.csv", "exact_overflow"),
                                   ("global_metrics.csv", "overflow"),
                                   ("recovery_metrics.csv", "overflow"),
                                   ("swap_recovery_metrics.csv", "overflow"),
                                   ("post_recovery_metrics.csv", "overflow"),
                                   ("post_swap_metrics.csv", "overflow")):
        path = run / filename
        if not path.is_file():
            continue
        rows = read_csv(path)
        if not rows:
            continue
        hpwl_key = "exact_hpwl" if "exact_hpwl" in rows[0] else "hpwl"
        key = overflow_key if overflow_key in rows[0] else (
            "exact_overflow" if "exact_overflow" in rows[0] else "overflow")
        parts.append((filename, rows, hpwl_key, key))
    if not parts:
        raise FileNotFoundError(f"no exact metric CSV in {run}")
    hpwl, overflow, splits = [], [], []
    for filename, rows, hpwl_key, overflow_key in parts:
        splits.append(len(hpwl))
        hpwl.extend(float(row[hpwl_key]) / 1e6 for row in rows)
        overflow.extend(100.0 * float(row[overflow_key]) for row in rows)
    return list(range(len(hpwl))), hpwl, overflow, splits


def main():
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"],
        "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
        "axes.spines.top": False, "axes.grid": True, "grid.alpha": 0.16,
        "legend.frameon": False, "figure.dpi": 300, "savefig.dpi": 300,
    })
    fig, axes = plt.subplots(len(RUNS), 1, figsize=(7.4, 2.8 * len(RUNS)),
                             constrained_layout=True, squeeze=False)
    for ax, color, (label, run) in zip(axes.flat, COLORS, RUNS.items()):
        x, hpwl, overflow, splits = load(label, run)
        right = ax.twinx()
        ax.plot(x, hpwl, color=color, linewidth=1.7, label="Exact HPWL")
        right.plot(x, overflow, color="#6B7280", linestyle="--", linewidth=1.35,
                   label="Exact overlap overflow")
        right.axhline(7.0, color="#009E73", linewidth=1.0, linestyle=":")
        right.axhspan(0, 7, color="#009E73", alpha=0.06)
        for split in splits[1:]:
            ax.axvline(split - 0.5, color="#CC79A7", linewidth=0.8, linestyle="-.")
        ax.set_title(label)
        ax.set_ylabel("HPWL (M)")
        right.set_ylabel("Overflow (%)")
        right.set_ylim(bottom=0)
        handles = ax.get_lines() + right.get_lines()[:1]
        ax.legend(handles, [h.get_label() for h in handles], loc="best")
    axes.flat[-1].set_xlabel("Iteration / exact-oracle update")
    out = ROOT / "output/visualizations"
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "selected_convergence.png", bbox_inches="tight")
    fig.savefig(out / "selected_convergence.pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()
