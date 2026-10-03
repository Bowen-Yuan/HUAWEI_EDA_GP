"""Run-artifact and offline-plot smoke test without a display server."""

from __future__ import annotations

from pathlib import Path

from gpplacer.metrics.logging import RunLogger
from gpplacer.model.config import SolverConfig
from gpplacer.solver.controller import PlacementSolver
from gpplacer.solver.state import PlacementState
from gpplacer.visualization.animation import animate_run
from gpplacer.visualization.plots import plot_run
from tests.test_core import _tiny_db


def _write_tiny_bookshelf(root: Path) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    aux = root / "tiny.aux"
    aux.write_text("RowBasedPlacement : tiny.nodes tiny.nets tiny.wts tiny.pl tiny.scl\n", encoding="ascii")
    (root / "tiny.nodes").write_text(
        "NumNodes : 3\nNumTerminals : 1\na 2 2\nb 2 2\nfixed 2 2 terminal\n", encoding="ascii",
    )
    (root / "tiny.nets").write_text(
        "NumNets : 1\nNumPins : 3\nNetDegree : 3 n0\na I : 0 0\nb I : 0 0\nfixed I : 0 0\n",
        encoding="ascii",
    )
    (root / "tiny.wts").write_text("", encoding="ascii")
    (root / "tiny.pl").write_text(
        "a 0 0 : N\nb 0 0 : N\nfixed 8 8 : N /FIXED\n", encoding="ascii",
    )
    rows = "".join(
        "CoreRow Horizontal\n"
        f" Coordinate : {y}\n Height : 2\n Sitewidth : 1\n SubrowOrigin : 0 NumSites : 10\nEnd\n"
        for y in (0, 2, 4, 6, 8)
    )
    (root / "tiny.scl").write_text(f"NumRows : 5\n{rows}", encoding="ascii")
    return aux


def test_logger_and_plotter_create_reproducible_artifacts() -> None:
    db = _tiny_db()
    db.source_aux = str(_write_tiny_bookshelf(Path.cwd() / "runs" / "tiny_bookshelf_input"))
    config = SolverConfig(max_iterations=1, max_seconds=30.0, threads=1, bin_rows=1,
                          snapshot_every=1, output_root=str(Path.cwd() / "runs"))
    logger = RunLogger(Path.cwd() / "runs", db, config)
    last: PlacementState | None = None

    def observe(record: dict[str, object], state: PlacementState) -> None:
        nonlocal last
        logger.observe(record, state)
        last = state

    result = PlacementSolver(db, config).solve(observe)
    assert last is not None
    logger.finish(status="completed", summary={"best_objective": result.oracle.objective}, state=last)
    figures = plot_run(logger.run_dir)
    assert figures
    assert all(path.with_suffix(".png").exists() for path in figures)
    movie = animate_run(logger.run_dir, fps=2, max_points=100, dpi=60)
    assert movie.suffix == ".gif"
    assert movie.exists()
