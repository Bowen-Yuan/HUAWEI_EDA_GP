#!/usr/bin/env python3
"""Plot concatenated exact trajectories for the strongest a1 branches."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt


def rows(path: Path, hpwl_key: str, overflow_key: str):
    with path.open(newline="", encoding="ascii") as stream:
        return [(float(r[hpwl_key]) / 1e6, 100.0 * float(r[overflow_key]))
                for r in csv.DictReader(stream)]


def load_run(run_dir: Path):
    result = []
    labels = []
    gp = run_dir / "global_metrics.csv"
    if gp.is_file():
        result.extend(rows(gp, "exact_hpwl", "overflow"))
        labels.extend(["GP"] * len(result))
    rp = run_dir / "recovery_metrics.csv"
    if rp.is_file():
        rec = rows(rp, "exact_hpwl", "overflow")
        result.extend(rec)
        labels.extend(["exact recovery"] * len(rec))
    return result, labels


def load_stages(stages):
    values = []
    labels = []
    boundaries = []
    for name, path in stages:
        start = len(values)
        current, current_labels = load_run(Path(path))
        values.extend(current)
        labels.extend([name + " / " + label for label in current_labels])
        boundaries.append((start, len(values), name))
    return values, labels, boundaries


def main():
    root = Path(__file__).resolve().parents[1]
    branches = {
        "H257 strict best": [
            ("H252 relaxed", root / "experiments/h252_a1_relaxed2_net_recovery"),
            ("H253 bridge", root / "experiments/h253_a1_relaxed2_bridge80"),
            ("H254 strong", root / "experiments/h254_a1_relaxed2_strong475"),
            ("H255 polish", root / "experiments/h255_a1_relaxed2_four_node_polish"),
            ("H257 polish", root / "experiments/h257_a1_fifth_node_line4"),
        ],
        "H271 continuous switch": [
            ("H270 GP", root / "experiments/h270_a1_inprocess_switch70_reset_adam"),
            ("H271 polish", root / "experiments/h271_a1_switch70_exact_polish"),
        ],
        "H275 low-HPWL near-feasible": [
            ("H273 GP", root / "experiments/h273_a1_switch70_step04"),
            ("H275 polish", root / "experiments/h275_a1_switch70_step04_four_recovery"),
        ],
    }
    plt.rcParams.update({
        "font.family": "serif", "font.size": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.alpha": 0.18,
    })
    fig, axes = plt.subplots(3, 1, figsize=(8.2, 8.4), sharex=False,
                             constrained_layout=True)
    for ax, (title, stages) in zip(axes, branches.items()):
        values, labels, boundaries = load_stages(stages)
        x = list(range(len(values)))
        hpwl = [v[0] for v in values]
        overflow = [v[1] for v in values]
        ax2 = ax.twinx()
        ax.plot(x, hpwl, color="#0072B2", linewidth=1.55, label="Exact HPWL")
        ax2.plot(x, overflow, color="#D55E00", linewidth=1.25,
                 linestyle="--", label="Exact overlap overflow")
        ax2.axhline(7.0, color="#009E73", linewidth=0.9)
        ax2.axhspan(0.0, 7.0, color="#009E73", alpha=0.06)
        for start, end, name in boundaries[:-1]:
            ax.axvline(end - 0.5, color="#777777", linewidth=0.65,
                       linestyle=":", alpha=0.75)
        ax.set_title(title, loc="left", fontweight="bold")
        ax.set_ylabel("HPWL (M)")
        ax2.set_ylabel("Overflow (%)")
        ax2.set_ylim(bottom=0.0)
        handles = ax.get_lines() + ax2.get_lines()[:1]
        ax.legend(handles, [h.get_label() for h in handles], frameon=False,
                  loc="best")
    axes[-1].set_xlabel("Concatenated exact-oracle stage sample")
    output = root / "to_human/a1_best_staged_convergence"
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")


if __name__ == "__main__":
    main()
