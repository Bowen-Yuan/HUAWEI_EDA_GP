#!/usr/bin/env python3
"""Plot exact HPWL and exact-overlap overflow from one or more GP runs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


def load_curve(run_dir: Path):
    rows = []
    for filename, stage in (("homotopy_metrics.csv", "homotopy"),
                            ("global_metrics.csv", "GP")):
        global_path = run_dir / filename
        if not global_path.is_file():
            continue
        with global_path.open(newline="", encoding="ascii") as stream:
            for row in csv.DictReader(stream):
                overflow_key = "exact_overflow" if "exact_overflow" in row else "overflow"
                rows.append((stage, row["exact_hpwl"], row[overflow_key]))
    for stage, filename in (("recovery", "recovery_metrics.csv"),
                            ("swap", "swap_recovery_metrics.csv"),
                            ("post-recovery", "post_recovery_metrics.csv"),
                            ("post-swap", "post_swap_metrics.csv")):
        path = run_dir / filename
        if not path.is_file():
            continue
        with path.open(newline="", encoding="ascii") as stream:
            reader = csv.DictReader(stream)
            for row in reader:
                hpwl_key = "exact_hpwl" if "exact_hpwl" in row else "hpwl"
                rows.append((stage, row[hpwl_key], row["overflow"]))
    if not rows:
        raise ValueError(f"empty metrics file: {run_dir}")
    return {
        "iteration": list(range(len(rows))),
        "stage": [row[0] for row in rows],
        "hpwl": [float(row[1]) / 1.0e6 for row in rows],
        "overflow": [100.0 * float(row[2]) for row in rows],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dirs", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.spines.top": False,
        "axes.grid": True,
        "grid.alpha": 0.18,
    })
    curves = [(path.name, load_curve(path)) for path in args.run_dirs]
    fig, axes = plt.subplots(len(curves), 1, figsize=(9.0, 3.2 * len(curves)),
                             squeeze=False, constrained_layout=True)
    for ax, (label, curve) in zip(axes.flat, curves):
        ax2 = ax.twinx()
        ax.plot(curve["iteration"], curve["hpwl"], color="#0072B2",
                linewidth=1.4, label="Exact HPWL")
        ax2.plot(curve["iteration"], curve["overflow"], color="#D55E00",
                 linewidth=1.15, linestyle="--", label="Exact overlap overflow")
        ax2.axhline(7.0, color="#009E73", linewidth=0.9)
        ax2.axhspan(0.0, 7.0, color="#009E73", alpha=0.06)
        ax.set_title(label)
        ax.set_xlabel("Exact-oracle stage sample")
        ax.set_ylabel("HPWL (M)")
        ax2.set_ylabel("Overflow (%)")
        ax2.set_ylim(bottom=0.0)
        handles = ax.get_lines() + ax2.get_lines()[:1]
        ax.legend(handles, [line.get_label() for line in handles], frameon=False)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output.with_suffix(".png"), dpi=240, bbox_inches="tight")
    fig.savefig(args.output.with_suffix(".pdf"), bbox_inches="tight")


if __name__ == "__main__":
    main()
