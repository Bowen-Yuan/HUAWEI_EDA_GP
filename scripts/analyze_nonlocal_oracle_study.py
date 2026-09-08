"""Create the non-local oracle convergence plots from metrics-only experiments."""
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt


def rows(root, prefix):
    result = []
    for run in sorted(Path(root).iterdir()):
        if prefix not in run.name or run.name.endswith("_finite"):
            continue
        path = run / "trajectory.csv"
        if not path.exists():
            continue
        with path.open(newline="") as stream:
            for row in csv.DictReader(stream):
                if row.get("record_type") == "iteration":
                    result.append((run.name, row))
    return result


def number(row, key):
    return float(row[key]) if row.get(key) else 0.0


def main():
    if len(sys.argv) not in (3, 4):
        raise SystemExit("usage: analyze_nonlocal_oracle_study.py EXPERIMENT_ROOT OUTPUT_DIR [RUN_PREFIX]")
    grouped = {}
    for name, row in rows(sys.argv[1], sys.argv[3] if len(sys.argv) == 4 else "nonlocal_primary_"):
        grouped.setdefault(name, []).append(row)
    output = Path(sys.argv[2]); output.mkdir(parents=True, exist_ok=True)
    for title, field, ylabel, log in [
        ("hpwl", "hpwl", "exact HPWL", False),
        ("overflow", "overflow_percent_after", "exact overflow (%)", False),
    ]:
        plt.figure(figsize=(8, 4.5))
        for name, data in grouped.items():
            plt.plot([number(r, "stage_iteration") for r in data],
                     [number(r, field) for r in data], label=name)
        if title == "overflow":
            plt.axhline(7, color="black", linestyle="--", linewidth=.8)
            plt.axhline(15, color="black", linestyle=":", linewidth=.8)
        if log: plt.yscale("log")
        plt.xlabel("iteration"); plt.ylabel(ylabel); plt.legend(fontsize=7); plt.tight_layout()
        plt.savefig(output / f"{title}.png", dpi=160); plt.close()
    plt.figure(figsize=(6, 4.5))
    for name, data in grouped.items():
        plt.plot([number(r, "overflow_percent_after") for r in data],
                 [number(r, "hpwl") for r in data], label=name)
    plt.xlabel("exact overflow (%)"); plt.ylabel("exact HPWL"); plt.legend(fontsize=7); plt.tight_layout()
    plt.savefig(output / "hpwl_overflow_phase.png", dpi=160)


if __name__ == "__main__":
    main()
