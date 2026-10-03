#!/usr/bin/env python3
"""Summarize and plot the eight ISPD2005 exact-HPWL bundle runs."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


BENCHMARKS = [
    "adaptec1", "adaptec2", "adaptec3", "adaptec4",
    "bigblue1", "bigblue2", "bigblue3", "bigblue4",
]

# DREAMPlace TCAD Table II, V100 float64 HPWL, in millions. These values are
# end-to-end reference context, not a controlled baseline for this solver.
DREAMPLACE_HPWL_M = {
    "adaptec1": 73.22,
    "adaptec2": 82.22,
    "adaptec3": 193.72,
    "adaptec4": 174.08,
    "bigblue1": 89.38,
    "bigblue2": 136.54,
    "bigblue3": 303.90,
    "bigblue4": 743.75,
}

COLORS = [
    "#E69F00", "#56B4E9", "#009E73", "#0072B2",
    "#D55E00", "#CC79A7", "#4D4D4D", "#F0E442",
]


def configure_matplotlib() -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif"],
        "font.size": 9.5,
        "axes.titlesize": 10.5,
        "axes.titleweight": "bold",
        "axes.labelsize": 9.5,
        "legend.fontsize": 7.5,
        "legend.frameon": False,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.18,
        "grid.linestyle": "-",
        "lines.linewidth": 1.5,
    })


def read_convergence(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline="", encoding="ascii") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 3000:
        raise ValueError(f"Expected 3000 convergence rows in {path}, got {len(rows)}")
    values = {
        # The solver's `iter` column is phase-local and resets at every phase.
        # CSV row order is the unambiguous global 3000-step trajectory.
        "iteration": np.arange(1, len(rows) + 1, dtype=int),
        "hpwl": np.asarray([float(row["hpwl"]) for row in rows]),
        "overflow": np.asarray([float(row["avg_overflow"]) for row in rows]),
        "lambda": np.asarray([float(row["lambda"]) for row in rows]),
    }
    for name, array in values.items():
        if not np.all(np.isfinite(array)):
            raise ValueError(f"Non-finite {name} values in {path}")
    return values


def parse_log(path: Path) -> dict[str, float | int | str]:
    text = path.read_text(encoding="utf-16", errors="ignore")
    if "Phase A (Subgrad):" not in text:
        text = path.read_text(encoding="utf-8", errors="ignore")

    def last_float(pattern: str) -> float:
        matches = re.findall(pattern, text, flags=re.MULTILINE)
        if not matches:
            raise ValueError(f"Pattern not found in {path}: {pattern}")
        return float(matches[-1])

    selected = re.findall(
        r"Restored best feasible fine-grid state: iter=(\d+) HPWL=([0-9.eE+-]+) overflow=([0-9.eE+-]+)",
        text,
    )
    selected_iter = int(selected[-1][0]) if selected else -1
    return {
        "hpwl": last_float(r"^\s+HPWL:\s+([0-9.eE+-]+)\s*$"),
        "overflow": last_float(r"^\s+Overflow:\s+([0-9.eE+-]+)\s*$"),
        "runtime_sec": last_float(r"^\s+Runtime:\s+([0-9.eE+-]+) sec\s*$"),
        "selected_iteration": selected_iter,
    }


def save_figure(fig: plt.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".png"), dpi=300)
    plt.close(fig)


def plot_per_benchmark(
    benchmark: str,
    data: dict[str, np.ndarray],
    result: dict[str, float | int | str],
    output_dir: Path,
    overflow_limit: float,
) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(6.75, 6.5), sharex=True)
    iteration = data["iteration"]
    axes[0].plot(iteration, data["hpwl"] / 1.0e6, color="#0072B2")
    axes[0].set_ylabel("Exact HPWL (M)")
    axes[0].set_title(f"{benchmark}: exact-HPWL bundle convergence")
    axes[1].plot(iteration, data["overflow"] * 100.0, color="#D55E00")
    axes[1].axhline(overflow_limit * 100.0, color="#555555", ls="--", lw=1.0,
                    label=f"{overflow_limit * 100:g}% reference limit")
    axes[1].set_ylabel("Overflow (%)")
    axes[1].legend(loc="upper right")
    axes[2].plot(iteration, np.maximum(data["lambda"], 1.0e-14), color="#009E73")
    axes[2].set_yscale("log")
    axes[2].set_ylabel("Effective lambda")
    axes[2].set_xlabel("Global iteration")

    selected_iteration = int(result["selected_iteration"])
    if selected_iteration > 0:
        for axis in axes:
            axis.axvline(selected_iteration, color="#CC79A7", ls=":", lw=1.1)
    fig.tight_layout()
    save_figure(fig, output_dir / f"fig_{benchmark}_convergence")


def plot_aggregate(
    all_data: dict[str, dict[str, np.ndarray]],
    output_dir: Path,
    overflow_limit: float,
) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(6.75, 5.3), sharex=True)
    for color, benchmark in zip(COLORS, BENCHMARKS):
        data = all_data[benchmark]
        hpwl_normalized = data["hpwl"] / max(data["hpwl"][0], 1.0)
        axes[0].plot(data["iteration"], hpwl_normalized, color=color, label=benchmark)
        axes[1].plot(data["iteration"], data["overflow"] * 100.0, color=color)
    axes[0].set_ylabel("HPWL / initial HPWL")
    axes[0].set_title("ISPD2005: normalized convergence across eight benchmarks")
    axes[0].legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.02))
    axes[1].axhline(overflow_limit * 100.0, color="#333333", ls="--", lw=1.0,
                    label=f"{overflow_limit * 100:g}% reference limit")
    axes[1].set_ylabel("Overflow (%)")
    axes[1].set_xlabel("Global iteration")
    axes[1].legend(loc="upper right")
    fig.tight_layout()
    save_figure(fig, output_dir / "fig_ispd2005_aggregate_convergence")


def plot_summary(results: list[dict[str, float | int | str]], output_dir: Path) -> None:
    labels = [str(row["benchmark"]) for row in results]
    hpwl_m = np.asarray([float(row["hpwl"]) / 1.0e6 for row in results])
    reference_m = np.asarray([DREAMPLACE_HPWL_M[name] for name in labels])
    overflow_pct = np.asarray([float(row["overflow"]) * 100.0 for row in results])
    runtime_min = np.asarray([float(row["runtime_sec"]) / 60.0 for row in results])
    x = np.arange(len(labels))

    fig, axes = plt.subplots(3, 1, figsize=(6.75, 7.1), sharex=True)
    width = 0.38
    axes[0].bar(x - width / 2, hpwl_m, width, color="#E76F51", label="This run (pre-LG)")
    axes[0].bar(x + width / 2, reference_m, width, color="#B0BEC5",
                label="DREAMPlace Table II (different flow)")
    axes[0].set_ylabel("HPWL (M)")
    axes[0].set_title("ISPD2005 bundle-run summary")
    axes[0].legend(ncol=2, loc="upper left")
    axes[1].bar(x, overflow_pct, width=0.62, color="#0072B2")
    axes[1].axhline(10.0, color="#333333", ls="--", lw=1.0)
    axes[1].set_ylabel("Overflow (%)")
    axes[2].bar(x, runtime_min, width=0.62, color="#009E73")
    axes[2].set_ylabel("Runtime (min)")
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(labels, rotation=30, ha="right")
    fig.tight_layout()
    save_figure(fig, output_dir / "fig_ispd2005_result_summary")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("results_root", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--summary-csv", type=Path, required=True)
    parser.add_argument("--overflow-limit", type=float, default=0.10)
    args = parser.parse_args()

    configure_matplotlib()
    all_data: dict[str, dict[str, np.ndarray]] = {}
    results: list[dict[str, float | int | str]] = []
    for benchmark in BENCHMARKS:
        run_dir = args.results_root / benchmark
        data = read_convergence(run_dir / "nsp_convergence.csv")
        result = parse_log(run_dir / "stdout.log")
        result.update({
            "benchmark": benchmark,
            "initialization": "eplace-ip.pl" if benchmark == "adaptec1" else "supplied-pl",
            "dreamplace_table2_hpwl_m": DREAMPLACE_HPWL_M[benchmark],
            "hpwl_ratio_to_dreamplace_table2":
                float(result["hpwl"]) / 1.0e6 / DREAMPLACE_HPWL_M[benchmark],
            "meets_10pct_reference": float(result["overflow"]) <= args.overflow_limit,
        })
        all_data[benchmark] = data
        results.append(result)
        plot_per_benchmark(benchmark, data, result, args.output_dir, args.overflow_limit)

    args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "benchmark", "initialization", "hpwl", "overflow", "runtime_sec",
        "selected_iteration", "meets_10pct_reference",
        "dreamplace_table2_hpwl_m", "hpwl_ratio_to_dreamplace_table2",
    ]
    with args.summary_csv.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    plot_aggregate(all_data, args.output_dir, args.overflow_limit)
    plot_summary(results, args.output_dir)


if __name__ == "__main__":
    main()
