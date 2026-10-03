"""Static and self-contained HTML visualisations for adaptive Pareto runs."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

from gpplacer.io.bookshelf import load_bookshelf
from gpplacer.oracle.grid import build_density_grid


_COLORS = {
    "hpwl": "#0072B2", "legacy": "#E69F00", "strict": "#D55E00",
    "joint": "#009E73", "rescue": "#7A5195", "plateau": "#8C564B",
}


def _style(dpi: int) -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9, "figure.dpi": 150,
        "savefig.dpi": dpi, "savefig.bbox": "tight", "axes.spines.top": False,
        "axes.spines.right": False, "axes.grid": True, "grid.alpha": .18,
    })


def _rows(root: Path) -> list[dict[str, str]]:
    path = root / "iterations.csv"
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _value(rows: list[dict[str, str]], key: str) -> np.ndarray:
    return np.asarray([float(row.get(key, "nan")) for row in rows], dtype=np.float64)


def _csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _save(fig: plt.Figure, target: Path, formats: tuple[str, ...], dpi: int) -> list[Path]:
    written: list[Path] = []
    for suffix in formats:
        path = target.with_suffix(f".{suffix}")
        fig.savefig(path, dpi=dpi if suffix == "png" else None)
        written.append(path)
    plt.close(fig)
    return written


def _phase_color(phase: str) -> str:
    if "joint" in phase:
        return _COLORS["joint"]
    if "plateau" in phase or "evacuation" in phase:
        return _COLORS["plateau"]
    if "rescue" in phase or "hpwl" in phase or "density" in phase:
        return _COLORS["rescue"]
    return "#777777"


def render_run(run_dir: str | Path) -> list[Path]:
    """Render core static figures and a fully offline interactive report."""
    root = Path(run_dir)
    meta = json.loads((root / "metadata.json").read_text(encoding="utf-8"))
    cfg = meta["config"]["visualization"]
    formats = tuple(cfg["static_formats"])
    dpi = int(cfg["dpi"])
    _style(dpi)
    figures = root / "figures"
    figures.mkdir(exist_ok=True)
    rows = _rows(root)
    if not rows:
        return []
    x = _value(rows, "iteration")
    hpwl = _value(rows, "hpwl")
    hnorm = _value(rows, "normalized_hpwl")
    legacy = _value(rows, "legacy_overflow_percent")
    strict = _value(rows, "strict_overflow_percent")
    blocked = _value(rows, "zero_capacity_occupancy")
    lam = _value(rows, "lambda")
    phase = [row["phase"] for row in rows]
    colors = [_phase_color(item) for item in phase]
    target_low = float(meta["config"]["evaluation"]["target_overflow_low"])
    target_high = float(meta["config"]["evaluation"]["target_overflow_high"])
    written: list[Path] = []

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 6.2), constrained_layout=True)
    axes[0, 0].plot(x, hpwl / 1e6, color=_COLORS["hpwl"])
    axes[0, 0].set(xlabel="Iteration", ylabel="HPWL (×10⁶)", title="Wirelength convergence")
    axes[0, 1].plot(x, legacy, color=_COLORS["legacy"], label="Legacy")
    axes[0, 1].plot(x, strict, color=_COLORS["strict"], label="Strict")
    axes[0, 1].axhspan(target_low, target_high, color=_COLORS["joint"], alpha=.12, label="Target band")
    axes[0, 1].set(xlabel="Iteration", ylabel="Overflow (%)", title="Fixed-protocol overflow")
    axes[0, 1].legend()
    axes[1, 0].plot(x, blocked, color=_COLORS["strict"])
    axes[1, 0].set(xlabel="Iteration", ylabel="Area", title="Zero-capacity movable occupancy")
    axes[1, 1].scatter(x, np.maximum(lam, 1e-12), c=colors, s=8)
    axes[1, 1].set_yscale("log")
    axes[1, 1].set(xlabel="Iteration", ylabel="λ", title="Adaptive penalty and phase")
    written.extend(_save(fig, figures / "01_convergence_overview", formats, dpi))

    archive_path = root / "pareto_archive.csv"
    archive_rows = _csv_rows(archive_path)
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.8), constrained_layout=True)
    axes[0].plot(strict, hpwl / 1e6, color="#BBBBBB", linewidth=.7)
    axes[0].scatter(strict, hpwl / 1e6, c=colors, s=9, alpha=.75)
    axes[0].axvspan(target_low, target_high, color=_COLORS["joint"], alpha=.12)
    if archive_rows:
        ax = np.asarray([float(row["strict_overflow_percent"]) for row in archive_rows])
        ay = np.asarray([float(row["hpwl"]) / 1e6 for row in archive_rows])
        order = np.argsort(ax)
        axes[0].plot(ax[order], ay[order], "o-", color="black", markersize=4, label="Archive")
        axes[0].legend()
    axes[0].set(xlabel="Strict overflow (%)", ylabel="HPWL (×10⁶)", title="Strict Pareto frontier")
    axes[1].scatter(legacy, strict, c=colors, s=10)
    limits = [min(np.nanmin(legacy), np.nanmin(strict)), max(np.nanmax(legacy), np.nanmax(strict))]
    axes[1].plot(limits, limits, "--", color="#555555", linewidth=.8)
    axes[1].set(xlabel="Legacy overflow (%)", ylabel="Strict overflow (%)", title="Protocol gap")
    written.extend(_save(fig, figures / "03_pareto_and_protocol_gap", formats, dpi))

    fig, axes = plt.subplots(2, 1, figsize=(10.0, 4.8), sharex=True, constrained_layout=True)
    phase_order = {"hpwl_rescue": 0, "density_rescue": 1, "joint": 2, "plateau_escape": 3, "refine": 4}
    codes = np.asarray([phase_order.get(item, 5) for item in phase])
    axes[0].step(x, codes, where="mid", color="#444444")
    axes[0].set(yticks=list(phase_order.values()), yticklabels=[item.replace("_", " ") for item in phase_order], title="State timeline")
    axes[1].plot(x, _value(rows, "step"), color="#0072B2", label="Step")
    axes[1].plot(x, _value(rows, "trust_p50"), color="#E69F00", label="Trust P50")
    axes[1].plot(x, _value(rows, "trust_p90"), color="#56B4E9", label="Trust P90")
    axes[1].set(xlabel="Iteration", ylabel="Distance", title="Trust-region control")
    axes[1].legend(ncol=3)
    written.extend(_save(fig, figures / "02_state_and_control", formats, dpi))

    fig, axes = plt.subplots(2, 1, figsize=(10.0, 5.2), sharex=True, constrained_layout=True)
    axes[0].plot(x, (hpwl - hpwl[0]) / max(hpwl[0], 1.0) * 100.0, color=_COLORS["hpwl"], label="HPWL change")
    axes[0].axhline(0.0, color="#555555", linewidth=.7)
    axes[0].set(ylabel="Change from first record (%)", title="Convergence progress")
    axes[0].legend()
    accepted = _value(rows, "accepted")
    rolling = np.convolve(accepted, np.ones(min(20, len(accepted))) / min(20, len(accepted)), mode="same")
    axes[1].plot(x, rolling * 100.0, color=_COLORS["joint"], label="Rolling acceptance")
    axes[1].plot(x, _value(rows, "single_axis_fraction") * 100.0, color=_COLORS["plateau"], label="Single-axis density use")
    axes[1].set(xlabel="Iteration", ylabel="Percent", ylim=(-2, 102), title="Acceptance and direction diagnostics")
    axes[1].legend(ncol=2)
    written.extend(_save(fig, figures / "04_convergence_diagnostics", formats, dpi))

    moves_path = root / "evacuation_moves.csv"
    if moves_path.exists() and moves_path.stat().st_size:
        with moves_path.open(encoding="utf-8", newline="") as handle:
            moves = list(csv.DictReader(handle))
        if moves:
            distance = np.asarray([np.hypot(float(row["dx"]), float(row["dy"])) for row in moves])
            accepted = np.asarray([int(row["accepted"]) for row in moves])
            fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.2), constrained_layout=True)
            axes[0].hist(distance[accepted > 0], bins=30, color=_COLORS["joint"], alpha=.8, label="Accepted")
            axes[0].hist(distance[accepted == 0], bins=30, color=_COLORS["strict"], alpha=.5, label="Rejected")
            axes[0].set(xlabel="Move distance", ylabel="Count", title="Local evacuation distances")
            axes[0].legend()
            axes[1].scatter([int(row["radius_bins"]) for row in moves], distance, c=accepted, cmap="RdYlGn", s=10)
            axes[1].set(xlabel="Cardinal radius (bins)", ylabel="Move distance", title="Evacuation proposals")
            written.extend(_save(fig, figures / "10_evacuation_summary", formats, dpi))

    written.extend(_render_candidate_and_hierarchy(root, figures, formats, dpi))
    written.extend(_render_spatial_diagnostics(root, meta, figures, formats, dpi))
    written.extend(_render_structure_and_runtime(root, meta, rows, figures, formats, dpi))

    manifest = {"figures": [str(path.relative_to(root)) for path in written], "report": "report.html"}
    (figures / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    if cfg.get("html", True):
        _render_html(root, meta, rows, target_low, target_high)
    return written


def _render_candidate_and_hierarchy(root: Path, figures: Path, formats: tuple[str, ...], dpi: int) -> list[Path]:
    written: list[Path] = []
    candidates = _csv_rows(root / "candidates.csv")
    if candidates:
        stages = list(dict.fromkeys(row["stage"] for row in candidates))
        fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.6), constrained_layout=True)
        counts = [sum(row["stage"] == stage for row in candidates) for stage in stages]
        survived = [sum(row["stage"] == stage and row.get("survived") == "1" for row in candidates) for stage in stages]
        positions = np.arange(len(stages))
        axes[0].bar(positions, counts, color="#BBBBBB", label="Evaluated")
        axes[0].bar(positions, survived, color=_COLORS["joint"], label="Survived")
        axes[0].set(xticks=positions, xticklabels=stages, ylabel="Candidates", title="Successive Halving funnel")
        axes[0].legend()
        generators = sorted(set(row["generator"] for row in candidates))
        palette = {name: color for name, color in zip(generators, ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00"])}
        for name in generators:
            local = [row for row in candidates if row["generator"] == name]
            axes[1].scatter(
                [float(row["strict_overflow_percent"]) for row in local],
                [float(row["hpwl"]) / 1e6 for row in local],
                label=name, color=palette[name], s=20, alpha=.75,
            )
        axes[1].set(xlabel="Strict overflow (%)", ylabel="HPWL (×10⁶)", title="Candidate quality by generator")
        axes[1].legend(fontsize=7)
        written.extend(_save(fig, figures / "05_candidate_funnel", formats, dpi))
        proxied = [row for row in candidates if row.get("proxy_score", "")]
        if proxied:
            fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.2), constrained_layout=True)
            predicted = np.asarray([float(row["proxy_score"]) for row in proxied])
            actual = np.asarray([float(row["strict_overflow_percent"]) for row in proxied])
            axes[0].scatter(predicted, actual, color=_COLORS["joint"], s=22)
            axes[0].set(xlabel="Ridge proxy score", ylabel="Strict overflow (%)", title="Proxy calibration")
            by_generator = sorted(set(row["generator"] for row in proxied))
            values = [[float(row["proxy_score"]) for row in proxied if row["generator"] == name] for name in by_generator]
            axes[1].boxplot(values, labels=by_generator)
            axes[1].set(ylabel="Proxy score", title="Predicted potential by generator")
            axes[1].tick_params(axis="x", rotation=25, labelsize=7)
            written.extend(_save(fig, figures / "06_candidate_proxy", formats, dpi))
    hierarchy = _csv_rows(root / "hierarchy.csv")
    if hierarchy:
        level = [int(row["level"]) for row in hierarchy]
        nodes = [int(row["nodes"]) for row in hierarchy]
        nets = [int(row["nets"]) for row in hierarchy]
        fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.2), constrained_layout=True)
        axes[0].plot(level, nodes, "o-", color=_COLORS["hpwl"], label="Nodes")
        axes[0].plot(level, nets, "o-", color=_COLORS["strict"], label="Nets")
        axes[0].set_yscale("log")
        axes[0].set(xlabel="Hierarchy level", ylabel="Count", title="Multilevel reduction")
        axes[0].legend()
        reduction = [1.0]
        for previous, current in zip(nodes, nodes[1:]):
            reduction.append(current / max(previous, 1))
        axes[1].bar(level, reduction, color=_COLORS["joint"])
        axes[1].set(xlabel="Hierarchy level", ylabel="Node ratio to previous level", title="Coarsening ratio", ylim=(0, 1.05))
        written.extend(_save(fig, figures / "07_multilevel_hierarchy", formats, dpi))
    return written


def _render_structure_and_runtime(
    root: Path, meta: dict[str, object], rows: list[dict[str, str]], figures: Path, formats: tuple[str, ...], dpi: int,
) -> list[Path]:
    """Render topology and measured controller throughput from persisted artifacts."""
    written: list[Path] = []
    source_aux = Path(str(meta.get("source_aux", "")))
    if source_aux.exists():
        db = load_bookshelf(source_aux)
        net_degree = np.diff(db.net_start)
        node_degree = np.diff(db.node_net_start)
        fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.2), constrained_layout=True)
        axes[0].hist(net_degree, bins=np.arange(1, min(int(net_degree.max()) + 2, 80)), color=_COLORS["hpwl"])
        axes[0].set(xlabel="Pins per net", ylabel="Count", title="Hyperedge degree distribution")
        axes[1].hist(node_degree, bins=np.arange(0, min(int(node_degree.max()) + 2, 80)), color=_COLORS["strict"])
        axes[1].set(xlabel="Incident nets per node", ylabel="Count", title="Cell connectivity distribution")
        written.extend(_save(fig, figures / "11_net_structure", formats, dpi))
    elapsed = _value(rows, "elapsed_seconds")
    if len(elapsed) and np.any(np.isfinite(elapsed)):
        fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.1), constrained_layout=True)
        axes[0].plot(_value(rows, "iteration"), elapsed, color=_COLORS["hpwl"])
        axes[0].set(xlabel="Iteration", ylabel="Wall time (s)", title="Controller runtime")
        accepted = _value(rows, "accepted")
        axes[1].plot(elapsed, np.cumsum(accepted) / np.maximum(elapsed, 1.0e-9), color=_COLORS["joint"])
        axes[1].set(xlabel="Wall time (s)", ylabel="Accepted steps / s", title="Acceptance throughput")
        written.extend(_save(fig, figures / "12_runtime_resources", formats, dpi))
    return written


def _render_spatial_diagnostics(
    root: Path, meta: dict[str, object], figures: Path, formats: tuple[str, ...], dpi: int,
) -> list[Path]:
    source_aux = Path(str(meta.get("source_aux", "")))
    snapshots = root / "snapshots"
    if not source_aux.exists() or not snapshots.exists():
        return []
    progress = sorted(path.stem for path in snapshots.glob("progress_*.npz"))
    selected_progress = [progress[len(progress) // 2]] if progress else []
    names = [name for name in (["input"] + selected_progress + ["best_target", "final"])
             if (snapshots / f"{name}.npz").exists()]
    if not names:
        return []
    db = load_bookshelf(source_aux)
    evaluation = meta["config"]["evaluation"]
    grid = build_density_grid(db, rho_target=float(evaluation["rho_target"]), bin_rows=int(evaluation["bin_rows"]))
    max_cells = int(meta["config"]["visualization"]["max_cells_static"])
    movable = np.flatnonzero(db.movable)
    if len(movable) > max_cells:
        movable = np.random.default_rng(0).choice(movable, size=max_cells, replace=False)
    fig, axes = plt.subplots(2, len(names), figsize=(3.6 * len(names), 6.2), constrained_layout=True, squeeze=False)
    image = None
    for column, name in enumerate(names):
        snapshot = np.load(snapshots / f"{name}.npz")
        centres = snapshot["centres"]
        occupancy = snapshot["occupancy"]
        axes[0, column].scatter(centres[movable, 0], centres[movable, 1], s=.25, alpha=.28, rasterized=True, color=_COLORS["hpwl"])
        fixed = np.flatnonzero(db.fixed)
        axes[0, column].scatter(centres[fixed, 0], centres[fixed, 1], s=3, color=_COLORS["strict"], rasterized=True)
        axes[0, column].set(title=name.replace("_", " "), xlabel="x", ylabel="y", aspect="equal")
        density = np.divide(occupancy, grid.available_area, out=np.full_like(occupancy, np.nan), where=grid.available_area > 0.0)
        image = axes[1, column].pcolormesh(grid.x_edges, grid.y_edges, density, shading="auto", cmap="YlOrRd", vmin=0.0, vmax=2.0)
        axes[1, column].contour(
            (grid.x_edges[:-1] + grid.x_edges[1:]) / 2.0, (grid.y_edges[:-1] + grid.y_edges[1:]) / 2.0,
            np.nan_to_num(density, nan=-1.0), levels=[grid.rho_target], colors=[_COLORS["hpwl"]], linewidths=.7,
        )
        axes[1, column].set(title=f"{name.replace('_', ' ')} density", xlabel="x", ylabel="y", aspect="equal")
    if image is not None:
        fig.colorbar(image, ax=axes[1], shrink=.78, label="Occupancy / available area")
    written = _save(fig, figures / "08_layout_stage_montage", formats, dpi)
    final = np.load(snapshots / "final.npz")
    occupancy = final["occupancy"]
    strict_active = final["strict_active_bins"]
    capacity = grid.rho_target * grid.available_area
    excess = np.maximum(occupancy - capacity, 0.0)
    zero_mask = grid.available_area <= 0.0
    fig, axes = plt.subplots(1, 3, figsize=(10.0, 3.2), constrained_layout=True)
    first = axes[0].pcolormesh(grid.x_edges, grid.y_edges, excess, shading="auto", cmap="Reds")
    axes[0].set(title="Strict positive-capacity excess", xlabel="x", ylabel="y", aspect="equal")
    fig.colorbar(first, ax=axes[0], shrink=.78)
    second = axes[1].pcolormesh(grid.x_edges, grid.y_edges, np.where(zero_mask, occupancy, np.nan), shading="auto", cmap="magma")
    axes[1].set(title="Movable occupancy in zero-capacity bins", xlabel="x", ylabel="y", aspect="equal")
    fig.colorbar(second, ax=axes[1], shrink=.78)
    axes[2].pcolormesh(grid.x_edges, grid.y_edges, strict_active.astype(float), shading="auto", cmap="Greys", vmin=0, vmax=1)
    axes[2].set(title="Strict active-bin mask", xlabel="x", ylabel="y", aspect="equal")
    written.extend(_save(fig, figures / "09_density_diagnostics", formats, dpi))
    return written


def _render_html(root: Path, meta: dict[str, object], rows: list[dict[str, str]], target_low: float, target_high: float) -> None:
    max_points = int(meta["config"]["visualization"]["max_trajectory_points_html"])
    if len(rows) > max_points:
        selected = np.linspace(0, len(rows) - 1, max_points).round().astype(int)
        rows = [rows[index] for index in np.unique(selected)]
    x = _value(rows, "iteration")
    hpwl = _value(rows, "hpwl") / 1e6
    legacy = _value(rows, "legacy_overflow_percent")
    strict = _value(rows, "strict_overflow_percent")
    phase = [row["phase"] for row in rows]
    overview = go.Figure()
    overview.add_trace(go.Scatter(x=x, y=hpwl, name="HPWL (×10⁶)", line={"color": _COLORS["hpwl"]}))
    overview.add_trace(go.Scatter(x=x, y=strict, name="Strict overflow (%)", yaxis="y2", line={"color": _COLORS["strict"]}))
    overview.add_trace(go.Scatter(x=x, y=legacy, name="Legacy overflow (%)", yaxis="y2", line={"color": _COLORS["legacy"], "dash": "dash"}))
    overview.update_layout(
        title="Adaptive Pareto overview", xaxis_title="Iteration", yaxis_title="HPWL (×10⁶)",
        yaxis2={"title": "Overflow (%)", "overlaying": "y", "side": "right"},
        shapes=[{"type": "rect", "xref": "x", "yref": "y2", "x0": float(x.min()), "x1": float(x.max()),
                 "y0": target_low, "y1": target_high, "fillcolor": _COLORS["joint"], "opacity": .12, "line_width": 0}],
    )
    pareto = go.Figure()
    pareto.add_trace(go.Scatter(
        x=strict, y=hpwl, mode="markers+lines", marker={"color": [_phase_color(item) for item in phase], "size": 6},
        text=phase, hovertemplate="Strict=%{x:.3f}%<br>HPWL=%{y:.3f}M<br>Phase=%{text}<extra></extra>", name="Trajectory",
    ))
    pareto.add_vrect(x0=target_low, x1=target_high, fillcolor=_COLORS["joint"], opacity=.12, line_width=0)
    pareto.update_layout(title="Strict Pareto trajectory", xaxis_title="Strict overflow (%)", yaxis_title="HPWL (×10⁶)")
    candidates = _csv_rows(root / "candidates.csv")
    candidate_html = "<p>No portfolio was run for this placement-seed experiment.</p>"
    if candidates:
        stages = list(dict.fromkeys(row["stage"] for row in candidates))
        candidate_figure = go.Figure()
        candidate_figure.add_trace(go.Bar(
            x=stages, y=[sum(row["stage"] == stage for row in candidates) for stage in stages], name="Evaluated", marker_color="#BBBBBB",
        ))
        candidate_figure.add_trace(go.Bar(
            x=stages, y=[sum(row["stage"] == stage and row.get("survived") == "1" for row in candidates) for stage in stages],
            name="Survived", marker_color=_COLORS["joint"],
        ))
        candidate_figure.update_layout(barmode="overlay", title="Candidate funnel", yaxis_title="Candidates")
        candidate_html = pio.to_html(candidate_figure, include_plotlyjs=False, full_html=False)
    hierarchy = _csv_rows(root / "hierarchy.csv")
    hierarchy_html = "<p>No hierarchy was run for this placement-seed experiment.</p>"
    if hierarchy:
        hierarchy_figure = go.Figure()
        hierarchy_figure.add_trace(go.Scatter(
            x=[int(row["level"]) for row in hierarchy], y=[int(row["nodes"]) for row in hierarchy],
            mode="lines+markers", name="Nodes", line={"color": _COLORS["hpwl"]},
        ))
        hierarchy_figure.add_trace(go.Scatter(
            x=[int(row["level"]) for row in hierarchy], y=[int(row["nets"]) for row in hierarchy],
            mode="lines+markers", name="Nets", line={"color": _COLORS["strict"]},
        ))
        hierarchy_figure.update_layout(title="Hierarchy reduction", xaxis_title="Level", yaxis={"type": "log", "title": "Count"})
        hierarchy_html = pio.to_html(hierarchy_figure, include_plotlyjs=False, full_html=False)
    moves = _csv_rows(root / "evacuation_moves.csv")
    evacuation_html = "<p>No local evacuation proposal was recorded.</p>"
    if moves:
        distance = [float(np.hypot(float(row["dx"]), float(row["dy"]))) for row in moves]
        evacuation_figure = go.Figure(go.Histogram(x=distance, marker_color=_COLORS["plateau"]))
        evacuation_figure.update_layout(title="Local evacuation move distances", xaxis_title="Distance", yaxis_title="Count")
        evacuation_html = pio.to_html(evacuation_figure, include_plotlyjs=False, full_html=False)
    metadata_html = "<pre>" + json.dumps(meta, ensure_ascii=False, indent=2) + "</pre>"
    overview_html = pio.to_html(overview, include_plotlyjs="inline", full_html=False)
    pareto_html = pio.to_html(pareto, include_plotlyjs=False, full_html=False)
    report_data = {
        "iteration": x.tolist(), "hpwl_million": hpwl.tolist(), "legacy_overflow": legacy.tolist(),
        "strict_overflow": strict.tolist(), "phase": phase, "metadata": meta,
    }
    (root / "report_data.json").write_text(json.dumps(report_data, ensure_ascii=False), encoding="utf-8")
    html = f"""<!doctype html><html><head><meta charset='utf-8'><title>Adaptive Pareto report</title>
<style>body{{font-family:Arial,sans-serif;margin:24px;max-width:1400px}}h1{{margin-bottom:0}}.tabs{{display:flex;gap:16px;margin:16px 0}}section{{margin:28px 0}}pre{{background:#f5f5f5;padding:12px;overflow:auto}}</style></head>
<body><h1>Adaptive Pareto report</h1><p>Status: {meta['status']} · strict overflow is the primary protocol.</p>
<div class='tabs'><a href='#overview'>Overview</a><a href='#pareto'>Pareto</a><a href='#candidates'>Candidates</a><a href='#hierarchy'>Hierarchy</a><a href='#evacuation'>Evacuation</a><a href='#figures'>Spatial/Static</a><a href='#config'>Configuration</a></div>
<section id='overview'>{overview_html}</section><section id='pareto'>{pareto_html}</section>
<section id='candidates'><h2>Candidates</h2>{candidate_html}</section><section id='hierarchy'><h2>Hierarchy</h2>{hierarchy_html}</section><section id='evacuation'><h2>Evacuation</h2>{evacuation_html}</section>
<section id='figures'><h2>Spatial and static figures</h2><p>See the figures directory for publication PNG/PDF exports, density diagnostics, layout montage and manifest.</p></section>
<section id='config'><h2>Configuration and provenance</h2>{metadata_html}</section></body></html>"""
    (root / "report.html").write_text(html, encoding="utf-8")
