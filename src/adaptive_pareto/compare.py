"""Cross-run comparison artifacts using each run's fixed strict evaluator outputs."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt


def compare_runs(run_dirs: list[str | Path], output: str | Path) -> Path:
    """Create a compact table and plot; never recompute incompatible metrics."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for item in run_dirs:
        root = Path(item)
        meta = json.loads((root / "metadata.json").read_text(encoding="utf-8"))
        rows.append({
            "run": str(root.resolve()), "case": root.parent.name, "status": meta["status"],
            "hpwl": float(meta["final_hpwl"]),
            "legacy_overflow_percent": float(meta["final_legacy_overflow_percent"]),
            "strict_overflow_percent": float(meta["final_strict_overflow_percent"]),
            "zero_capacity_occupancy": float(meta["final_zero_capacity_occupancy"]),
            "elapsed_seconds": float(meta["elapsed_seconds"]),
        })
    with (output / "comparison.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    labels = [f"{row['case']}\n{Path(str(row['run'])).name}" for row in rows]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.5), constrained_layout=True)
    axes[0].bar(labels, [float(row["hpwl"]) / 1e6 for row in rows], color="#0072B2")
    axes[0].set(ylabel="HPWL (M)", title="Final wirelength")
    axes[1].bar(labels, [float(row["strict_overflow_percent"]) for row in rows], color="#D55E00", label="Strict")
    axes[1].bar(labels, [float(row["legacy_overflow_percent"]) for row in rows], color="#E69F00", alpha=.55, label="Legacy")
    axes[1].set(ylabel="Overflow (%)", title="Fixed-protocol final overflow")
    axes[1].legend()
    for axis in axes:
        axis.tick_params(axis="x", rotation=18, labelsize=7)
        axis.grid(axis="y", alpha=.18)
    figure = output / "comparison.png"
    fig.savefig(figure, dpi=220, bbox_inches="tight")
    plt.close(fig)
    table = "".join("<tr>" + "".join(f"<td>{value}</td>" for value in row.values()) + "</tr>" for row in rows)
    headers = "".join(f"<th>{key}</th>" for key in rows[0])
    (output / "comparison.html").write_text(
        f"<!doctype html><meta charset='utf-8'><title>Adaptive comparison</title>"
        f"<style>body{{font-family:Arial;margin:24px}}table{{border-collapse:collapse}}td,th{{border:1px solid #ccc;padding:6px}}</style>"
        f"<h1>Adaptive Pareto comparison</h1><img src='comparison.png' width='900'><table><tr>{headers}</tr>{table}</table>",
        encoding="utf-8",
    )
    return output
