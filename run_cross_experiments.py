"""Run the two adaptec1 cross-solver experiments and compare exact metrics."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys
import time


REPO_ROOT = Path(__file__).resolve().parents[1]
CODE_ROOT = REPO_ROOT / "code"
ADAPTIVE_SRC = REPO_ROOT / "adaptive_pareto_soft" / "src"
for source_root in (str(CODE_ROOT), str(ADAPTIVE_SRC)):
    if source_root not in sys.path:
        sys.path.insert(0, source_root)

from gpplacer.io.bookshelf import load_bookshelf
from gpplacer.io.placement import read_placement_centres, write_placement
from gpplacer.multilevel.initialization import capacity_aware_seed
from gpplacer.solver.steps import project_centres
from gpplacer.official import run_official_evaluation

from adaptive_pareto.config import AdaptiveConfig
from adaptive_pareto.evaluator import AdaptiveEvaluator


def metrics(evaluator: AdaptiveEvaluator, centres) -> dict[str, float]:
    """Optimization-only proxy metrics; never use them as contest results."""
    result = evaluator.evaluate(centres)
    return {
        "hpwl": float(result.hpwl),
        "normalized_hpwl": float(evaluator.normalized_hpwl(result)),
        "legacy_overflow_percent": float(evaluator.overflow_percent(result, strict=False)),
        "strict_overflow_percent": float(evaluator.overflow_percent(result, strict=True)),
        "zero_capacity_occupancy": float(result.strict.zero_capacity_occupancy),
    }


def official_metrics(aux: Path, placement: Path, *, perl: str | None, density_target: float) -> dict[str, object]:
    """Run the supplied scripts when a Perl runtime was explicitly provided."""
    if perl is None:
        return {"status": "not_run", "reason": "pass --official-perl to run the supplied ISPD scripts"}
    try:
        return {"status": "ok", **run_official_evaluation(
            aux, placement, density_target=density_target, perl=perl,
        ).to_dict()}
    except RuntimeError as exc:
        # A global placement can be rejected by the official checker.  Preserve
        # that fact instead of falling back to the internal strict surrogate.
        return {"status": "rejected", "reason": str(exc)}


def change(before: dict[str, float], after: dict[str, float]) -> dict[str, float | None]:
    result = {}
    for key in ("hpwl", "legacy_overflow_percent", "strict_overflow_percent"):
        result[f"{key}_delta"] = after[key] - before[key]
        result[f"{key}_change_percent"] = (
            100.0 * (after[key] - before[key]) / before[key] if before[key] else None
        )
    return result


def run_logged(
    command: list[str], log_path: Path, env: dict[str, str] | None = None,
    cwd: Path = REPO_ROOT,
) -> float:
    started = time.perf_counter()
    with log_path.open("w", encoding="utf-8") as log:
        completed = subprocess.run(
            command, cwd=cwd, env=env, stdout=log, stderr=subprocess.STDOUT,
            text=True, check=False,
        )
    elapsed = time.perf_counter() - started
    if completed.returncode:
        raise RuntimeError(
            f"Command failed with exit code {completed.returncode}; see {log_path}"
        )
    return elapsed


def adaptive_run_dir(output_root: Path) -> Path:
    metadata = list(output_root.rglob("metadata.json"))
    if not metadata:
        raise RuntimeError(f"Adaptive solver did not create metadata below {output_root}")
    return max(metadata, key=lambda path: path.stat().st_mtime).parent


def markdown(summary: dict) -> str:
    lines = [
        "# adaptec1 双向交叉求解实验",
        "",
        "正式 HPWL 与 overflow 仅由数据集随附的 ISPD Perl 脚本复算；内部 strict/legacy 值只是优化代理。",
        "",
        "| 实验 | 阶段 | 官方 HPWL | Scaled Overflow per bin |",
        "|---|---:|---:|---:|---:|",
    ]
    for experiment, label in (("adaptive_init_to_alg", "Adaptive init → Alg"),
                              ("a1_default_to_adaptive", "a1 default → Adaptive")):
        payload = summary[experiment]
        stage_labels = {
            "raw_input": "raw file", "input": "solver input",
            "global_output": "pre-legalize", "output": "final",
        }
        for stage in ("raw_input", "input", "global_output", "output"):
            if stage not in payload:
                continue
            official = payload.get(f"{stage}_official", {})
            if official.get("status") == "ok":
                lines.append(
                    f"| {label} | {stage_labels[stage]} | {official['hpwl']} | "
                    f"{official['scaled_overflow_per_bin']:.6f} |"
                )
            else:
                lines.append(f"| {label} | {stage_labels[stage]} | rejected / not run | rejected / not run |")
    lines.extend(("", f"原始数据：`{summary['summary_json']}`", ""))
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--solver", type=Path, help="Path to nsp_placer or nsp_placer.exe")
    parser.add_argument(
        "--config", type=Path,
        default=REPO_ROOT / "combine" / "configs" / "a1_adaptive_from_default.json",
    )
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--alg-p0-iters", type=int, default=150)
    parser.add_argument("--alg-p1-iters", type=int, default=100)
    parser.add_argument("--alg-p2-iters", type=int, default=150)
    parser.add_argument("--skip-alg", action="store_true")
    parser.add_argument("--skip-adaptive", action="store_true")
    parser.add_argument("--official-perl", help="Perl runtime used to run the supplied official scripts")
    parser.add_argument("--official-density-target", type=float, default=0.60)
    parser.add_argument("--output-root", type=Path, default=REPO_ROOT / "combine" / "results")
    parser.add_argument(
        "--run-root", type=Path,
        help="Reuse an existing run directory, preserving any completed skipped side",
    )
    args = parser.parse_args()

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_root = (args.run_root.resolve() if args.run_root else
                args.output_root.resolve() / f"adaptec1_cross_{stamp}")
    run_root.mkdir(parents=True, exist_ok=args.run_root is not None)
    aux = REPO_ROOT / "dataset" / "adaptec1" / "adaptec1.aux"
    benchmark = aux.with_suffix("")
    default_pl = aux.with_suffix(".pl")
    config = AdaptiveConfig.from_json(args.config)
    config.runtime.seed = args.seed
    db = load_bookshelf(aux)
    evaluator = AdaptiveEvaluator(db, config.evaluation)

    summary_path = run_root / "summary.json"
    summary: dict[str, object] = (
        json.loads(summary_path.read_text(encoding="utf-8")) if summary_path.exists() else {}
    )
    summary.update({
        "benchmark": str(aux),
        "seed": args.seed,
        "evaluation": config.to_dict()["evaluation"],
    })

    adaptive_seed = capacity_aware_seed(db, evaluator.grid, args.seed)
    adaptive_seed_pl = run_root / "adaptive_capacity_seed.pl"
    write_placement(adaptive_seed_pl, db, adaptive_seed)
    alg_output = run_root / "alg_from_adaptive_seed.pl"
    alg_global_output = run_root / "alg_from_adaptive_seed_global.pl"
    alg_payload: dict[str, object] = {
        "input_placement": str(adaptive_seed_pl),
        "input": metrics(evaluator, adaptive_seed),
        "input_official": official_metrics(
            aux, adaptive_seed_pl, perl=args.official_perl, density_target=args.official_density_target,
        ),
    }
    if args.skip_alg and "adaptive_init_to_alg" in summary:
        alg_payload.update(summary["adaptive_init_to_alg"])
    if not args.skip_alg:
        solver = args.solver
        if solver is None:
            candidates = (
                REPO_ROOT / "combine" / "Alg" / "Alg" / "nsp_placer.exe",
                REPO_ROOT / "combine" / "Alg" / "Alg" / "nsp_placer",
            )
            solver = next((path for path in candidates if path.exists()), None)
        if solver is None or not solver.exists():
            raise FileNotFoundError("Alg solver binary not found; build it or pass --solver")
        command = [
            str(solver.resolve()), str(benchmark), "--init-pl", str(adaptive_seed_pl),
            "--prelegal-output", str(alg_global_output), "--output", str(alg_output),
            "--p0-iters", str(args.alg_p0_iters),
            "--p1-iters", str(args.alg_p1_iters), "--p2-iters", str(args.alg_p2_iters),
        ]
        alg_payload["runtime_seconds"] = run_logged(
            command, run_root / "alg.log", cwd=run_root,
        )
        alg_payload["global_output_placement"] = str(alg_global_output)
        alg_payload["global_output"] = metrics(
            evaluator, read_placement_centres(alg_global_output, db),
        )
        alg_payload["global_output_official"] = official_metrics(
            aux, alg_global_output, perl=args.official_perl, density_target=args.official_density_target,
        )
        alg_payload["global_change"] = change(
            alg_payload["input"], alg_payload["global_output"],
        )
        alg_payload["output_placement"] = str(alg_output)
        alg_payload["output"] = metrics(evaluator, read_placement_centres(alg_output, db))
        alg_payload["output_official"] = official_metrics(
            aux, alg_output, perl=args.official_perl, density_target=args.official_density_target,
        )
        alg_payload["change"] = change(alg_payload["input"], alg_payload["output"])
    summary["adaptive_init_to_alg"] = alg_payload

    default_centres = read_placement_centres(default_pl, db)
    projected_default = project_centres(db, default_centres)
    adaptive_payload: dict[str, object] = {
        "input_placement": str(default_pl),
        "raw_input": metrics(evaluator, default_centres),
        "raw_input_official": official_metrics(
            aux, default_pl, perl=args.official_perl, density_target=args.official_density_target,
        ),
        "input": metrics(evaluator, projected_default),
        "projection_change": change(
            metrics(evaluator, default_centres), metrics(evaluator, projected_default),
        ),
    }
    if args.skip_adaptive and "a1_default_to_adaptive" in summary:
        adaptive_payload.update(summary["a1_default_to_adaptive"])
        adaptive_payload["raw_input"] = metrics(evaluator, default_centres)
        adaptive_payload["input"] = metrics(evaluator, projected_default)
        if "output" in adaptive_payload:
            adaptive_payload["change"] = change(
                adaptive_payload["input"], adaptive_payload["output"],
            )
    if not args.skip_adaptive:
        adaptive_output_root = run_root / "adaptive_from_a1_default"
        env = os.environ.copy()
        env["PYTHONPATH"] = os.pathsep.join((str(CODE_ROOT), str(ADAPTIVE_SRC)))
        command = [
            sys.executable, "-m", "adaptive_pareto.cli", "solve", "--aux", str(aux),
            "--placement", str(default_pl), "--config", str(args.config),
            "--output-root", str(adaptive_output_root),
        ]
        adaptive_payload["runtime_seconds"] = run_logged(
            command, run_root / "adaptive.log", env=env,
        )
        adaptive_root = adaptive_run_dir(adaptive_output_root)
        adaptive_output = adaptive_root / "solutions" / "solution_final.pl"
        adaptive_payload["run_dir"] = str(adaptive_root)
        adaptive_payload["output_placement"] = str(adaptive_output)
        adaptive_payload["output"] = metrics(
            evaluator, read_placement_centres(adaptive_output, db),
        )
        adaptive_payload["output_official"] = official_metrics(
            aux, adaptive_output, perl=args.official_perl, density_target=args.official_density_target,
        )
        adaptive_payload["change"] = change(adaptive_payload["input"], adaptive_payload["output"])
    summary["a1_default_to_adaptive"] = adaptive_payload

    summary["summary_json"] = str(summary_path)
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    (run_root / "REPORT.md").write_text(markdown(summary), encoding="utf-8")
    print(json.dumps({"run_dir": str(run_root), "summary": str(summary_path)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
