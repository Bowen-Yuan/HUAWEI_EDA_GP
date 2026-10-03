"""Command-line interface for independent RBSM experiments."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np

from .bookshelf import load_bookshelf
from .metrics import evaluate
from .solver import RBSMSolver, SolverConfig, config_from_json


def _read_pl(design, path: Path) -> np.ndarray:
    index = {name: i for i, name in enumerate(design.names)}; result = design.centres.copy()
    for raw in path.read_text(encoding="ascii").splitlines():
        f = raw.split()
        if len(f) >= 5 and f[0] in index and f[3] == ":":
            i = index[f[0]]; result[i] = (float(f[1]) + design.width[i] / 2, float(f[2]) + design.height[i] / 2)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(prog="rbsm.py")
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("inspect", "solve", "evaluate"):
        p = sub.add_parser(name); p.add_argument("--aux", required=True)
    solve = sub.choices["solve"]; solve.add_argument("--config", required=True); solve.add_argument("--seed", type=int); solve.add_argument("--output-dir")
    ev = sub.choices["evaluate"]; ev.add_argument("--placement", required=True); ev.add_argument("--bins", type=int, default=512)
    bench = sub.add_parser("benchmark-a1")
    bench.add_argument("--aux", required=True); bench.add_argument("--config", required=True)
    bench.add_argument("--seeds", type=int, nargs="+", default=[2026, 2027, 2028])
    bench.add_argument("--output-root", default="runs/adaptec1")
    args = parser.parse_args(); design = load_bookshelf(args.aux)
    if args.command == "inspect":
        print(json.dumps({"nodes": len(design.names), "movable": int(design.movable.sum()), "fixed": int(design.fixed.sum()), "nets": int(design.net_weight.size), "pins": int(design.pin_node.size), "core_bounds": design.core_bounds}, indent=2)); return
    if args.command == "evaluate":
        print(json.dumps(evaluate(design, _read_pl(design, Path(args.placement)), args.bins), indent=2)); return
    if args.command == "benchmark-a1":
        base = config_from_json(args.config); root = Path(args.output_root); records = []
        for seed in args.seeds:
            config = SolverConfig(**vars(base)); config.seed = seed
            result = RBSMSolver(design, config).run(root / f"seed{seed}")
            records.append(result)
        targets = {"hpwl": 5.05e7, "density_overflow_percent": 26.8, "pair_overlap_percent": 10.36}
        best = min(records, key=lambda x: x["hpwl"])
        report = {"targets": targets, "runs": records, "best_hpwl_run": best,
                  "exact_hit": any(all(r[k] <= v for k, v in targets.items()) for r in records),
                  "comparable_5pct": any(all(r[k] <= v * 1.05 for k, v in targets.items()) for r in records)}
        (root / "benchmark_summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(json.dumps(report, indent=2)); return
    config = config_from_json(args.config)
    if args.seed is not None: config.seed = args.seed
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    output = Path(args.output_dir) if args.output_dir else Path("runs") / "adaptec1" / f"{stamp}-seed{config.seed}"
    print(json.dumps(RBSMSolver(design, config).run(output), indent=2))
