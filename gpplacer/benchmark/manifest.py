"""Baseline manifest validation and case-by-case acceptance summary."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def evaluate_manifest(path: str | Path) -> dict[str, Any]:
    """Validate a recorded benchmark comparison without inventing missing data.

    The manifest contains ``baseline`` and ``candidate`` records per case with
    HPWL, overflow, runtime seconds, iteration count, and thread count.  This
    function only reports ratios; invoking the actual solver is the CLI's job.
    """
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    cases = payload.get("cases", [])
    if len(cases) != 8:
        raise ValueError("Formal challenge acceptance requires exactly eight cases.")
    hpwl_ratios: list[float] = []
    rows: list[dict[str, Any]] = []
    for item in cases:
        baseline, candidate = item["baseline"], item["candidate"]
        if baseline["threads"] != candidate["threads"]:
            raise ValueError(f"{item['name']}: baseline and candidate thread counts differ.")
        hpwl_ratio = candidate["hpwl"] / baseline["hpwl"]
        runtime_ratio = candidate["runtime_seconds"] / baseline["runtime_seconds"]
        iteration_ratio = candidate["iterations"] / baseline["iterations"]
        overflow_ratio = candidate["overflow"] / max(baseline["overflow"], 1e-12)
        hpwl_ratios.append(hpwl_ratio)
        rows.append({
            "name": item["name"], "hpwl_ratio": hpwl_ratio,
            "overflow_ratio": overflow_ratio, "runtime_ratio": runtime_ratio,
            "iteration_ratio": iteration_ratio,
        })
    return {
        "cases": rows,
        "average_hpwl_ratio": sum(hpwl_ratios) / len(hpwl_ratios),
        "meets_hpwl": sum(hpwl_ratios) / len(hpwl_ratios) < 1.01,
        "meets_runtime": all(row["runtime_ratio"] <= 2.0 for row in rows),
        "meets_iterations": all(row["iteration_ratio"] <= 5.0 for row in rows),
        # "Same order" is intentionally reported rather than guessed; the
        # program surfaces ratios for the challenge owner to apply its policy.
        "overflow_ratios": {row["name"]: row["overflow_ratio"] for row in rows},
    }
