"""Summarise completed adaptec1 runs into a comparable CSV table.

The CSV log format evolved during development.  In particular, older logs do
not contain ``overflow_ratio``.  This utility recomputes that ratio as
``density_linear / total_available_area`` using the same density-grid
construction as the solver, so all completed runs remain comparable.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "code"))

from gpplacer.io.bookshelf import load_bookshelf  # noqa: E402
from gpplacer.oracle.grid import build_density_grid  # noqa: E402


OUTPUT_FIELDS = [
    "run_id",
    "created_utc",
    "run_kind",
    "iterations",
    "elapsed_seconds",
    "phase",
    "hpwl",
    "density_linear",
    "overflow_ratio",
    "overflow_percent",
    "rho_target",
    "bin_rows",
    "available_area",
    "metrics_source",
]


def _read_json(path: Path) -> dict[str, object]:
    with path.open("r", encoding="utf-8-sig") as handle:
        return json.load(handle)


def _last_iteration(path: Path) -> dict[str, str] | None:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return next(reversed(list(csv.DictReader(handle))), None)


def _number(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    if not value:
        raise ValueError(f"missing {key}")
    return float(value)


def completed_rows(run_root: Path, aux_path: Path) -> tuple[list[dict[str, object]], list[str]]:
    """Collect completed runs that have both a placement and final metrics."""
    db = load_bookshelf(aux_path)
    available_area_by_bin_rows: dict[int, float] = {}
    rows: list[dict[str, object]] = []
    excluded: list[str] = []

    for run_dir in sorted(path for path in run_root.iterdir() if path.is_dir()):
        metadata_path = run_dir / "metadata.json"
        iterations_path = run_dir / "iterations.csv"
        solution_path = run_dir / "solution.pl"
        if not (metadata_path.exists() and iterations_path.exists() and solution_path.exists()):
            excluded.append(f"{run_dir.name}: missing metadata, iterations log, or solution")
            continue

        try:
            metadata = _read_json(metadata_path)
            final = _last_iteration(iterations_path)
            if metadata.get("status") != "completed" or final is None:
                raise ValueError("run is not completed or has no iteration row")
            config = metadata.get("config", {})
            if not isinstance(config, dict):
                raise ValueError("invalid config in metadata")
            bin_rows = int(config.get("bin_rows", 8))
            rho_target = float(config.get("rho_target", 0.8))
            if bin_rows not in available_area_by_bin_rows:
                grid = build_density_grid(db, rho_target=rho_target, bin_rows=bin_rows)
                available_area_by_bin_rows[bin_rows] = float(grid.available_area.sum())
            available_area = available_area_by_bin_rows[bin_rows]
            density_linear = _number(final, "density_linear")
            overflow_ratio = density_linear / available_area
            iterations = int(metadata.get("iterations", final.get("iteration", 0)))
            rows.append({
                "run_id": run_dir.name,
                "created_utc": metadata.get("created_utc", ""),
                "run_kind": "initialization_only" if iterations == 0 else "solver_run",
                "iterations": iterations,
                "elapsed_seconds": float(metadata.get("elapsed_seconds", 0.0)),
                "phase": final.get("phase", ""),
                "hpwl": _number(final, "hpwl"),
                "density_linear": density_linear,
                "overflow_ratio": overflow_ratio,
                "overflow_percent": 100.0 * overflow_ratio,
                "rho_target": rho_target,
                "bin_rows": bin_rows,
                "available_area": available_area,
                "metrics_source": "final_iteration; overflow recomputed from density_linear",
            })
        except (OSError, ValueError, TypeError, json.JSONDecodeError) as error:
            excluded.append(f"{run_dir.name}: {error}")
    return rows, excluded


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, default=PROJECT_ROOT / "runs" / "adaptec1")
    parser.add_argument("--aux", type=Path, default=PROJECT_ROOT / "dataset" / "adaptec1" / "adaptec1.aux")
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    output = args.output or args.run_root / "adaptec1_valid_results_summary.csv"
    rows, excluded = completed_rows(args.run_root, args.aux)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} completed results to {output}")
    if excluded:
        print(f"Excluded {len(excluded)} incomplete/unreadable runs:")
        print("\n".join(excluded))


if __name__ == "__main__":
    main()
