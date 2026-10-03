#!/usr/bin/env python3
"""Generate the annotated figures and audit tables for the H375 report.

All plotted optimization values are read from the archived experiment CSV files.
The only manually supplied numbers are the DREAMPlace paper baseline (Table II)
and stage runtimes already recorded in each experiment's summary.txt.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
EXP = ROOT / "experiments"
OUT = Path(__file__).resolve().parent
DATA_OUT = OUT.parent / "data"

BASELINE_HPWL_M = 73.22
BASELINE_RUNTIME_S = 67.0
RUNTIME_GATE_S = 2.0 * BASELINE_RUNTIME_S

COLORS = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "red": "#D55E00",
    "purple": "#CC79A7",
    "sky": "#56B4E9",
    "yellow": "#F0E442",
    "gray": "#7A7A7A",
    "lightgray": "#D9D9D9",
}

mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9.0,
        "axes.titlesize": 10.5,
        "axes.labelsize": 9.5,
        "legend.fontsize": 8.0,
        "xtick.labelsize": 8.0,
        "ytick.labelsize": 8.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.22,
        "grid.linewidth": 0.6,
        "lines.linewidth": 1.65,
        "lines.markersize": 4.2,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.06,
    }
)


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def read_kv(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        if "=" not in raw:
            continue
        key, value = raw.split("=", 1)
        result[key.strip()] = value.strip()
    return result


def save_figure(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.pdf")
    fig.savefig(OUT / f"{stem}.png", dpi=260)
    plt.close(fig)


def contiguous_spans(labels: Iterable[str]) -> list[tuple[str, int, int]]:
    labels = list(labels)
    if not labels:
        return []
    spans: list[tuple[str, int, int]] = []
    start = 0
    current = labels[0]
    for index, label in enumerate(labels[1:], start=1):
        if label != current:
            spans.append((current, start, index - 1))
            current = label
            start = index
    spans.append((current, start, len(labels) - 1))
    return spans


def shade_stages(
    axes: Iterable[plt.Axes], spans: list[tuple[str, int, int]],
    *, label_axis: plt.Axes | None = None, y: float = 1.02,
) -> None:
    palette = ["#E8F1F8", "#FFF0D8", "#E6F4EE", "#F5EAF2"]
    for number, (name, start, end) in enumerate(spans):
        for ax in axes:
            ax.axvspan(start - 0.48, end + 0.48, color=palette[number % len(palette)],
                       alpha=0.62, zorder=-10, linewidth=0)
        if label_axis is not None:
            width = end - start + 1
            inside = y < 1.0
            label_axis.text(
                0.5 * (start + end), y, name, ha="center",
                va="top" if inside else "bottom",
                # Short one- or few-point stages otherwise crowd one another
                # at the handoff boundaries (notably H372/H375).  Keeping the
                # label anchored at the stage centre makes the annotation
                # readable in both the full-chain and exact-tail panels.
                rotation=90 if width <= 8 else 0,
                transform=label_axis.get_xaxis_transform(), fontsize=6.7,
                color="#333333", clip_on=False,
            )


def full_chain_records() -> pd.DataFrame:
    seed = read_csv(EXP / "h219_a1_hpwl_seed" / "global_metrics.csv")
    coarse = read_csv(EXP / "h219_a1_coarse" / "homotopy_metrics.csv")
    fine = read_csv(EXP / "h221_a1_fine" / "homotopy_metrics.csv")
    r252 = read_csv(EXP / "h252_a1_relaxed2_net_recovery" / "recovery_metrics.csv")
    g253 = read_csv(EXP / "h253_a1_relaxed2_bridge80" / "global_metrics.csv")
    s254 = read_kv(EXP / "h254_a1_relaxed2_strong475" / "summary.txt")
    r255 = read_csv(EXP / "h255_a1_relaxed2_four_node_polish" / "recovery_metrics.csv")
    r257 = read_csv(EXP / "h257_a1_fifth_node_line4" / "recovery_metrics.csv")
    b372 = read_csv(
        EXP / "h372_a1_surplus_only_no_axis_map" / "bisect" / "bisection_metrics.csv"
    )
    r372 = read_csv(
        EXP / "h372_a1_surplus_only_no_axis_map" / "recovery10" / "recovery_metrics.csv"
    )
    s375 = read_csv(
        EXP / "h375_a1_surplus_recovery10_exchange" / "run" / "swap_recovery_metrics.csv"
    )

    records: list[dict[str, object]] = []

    def add(group: str, label: str, hpwl: float, overflow: float, source: str) -> None:
        records.append(
            {
                "group": group,
                "label": label,
                "hpwl_m": hpwl / 1.0e6,
                "overflow_pct": 100.0 * overflow,
                "source": source,
            }
        )

    add("Init / HPWL", "raw fresh", seed.iloc[0].exact_hpwl,
        seed.iloc[0].overflow, "H219 HPWL seed row 0")
    add("Init / HPWL", "H219 seed", seed.iloc[-1].exact_hpwl,
        seed.iloc[-1].overflow, "H219 HPWL seed row 199")
    c1100 = coarse.loc[coarse.iteration == 1100].iloc[0]
    add("128-grid", "H219 coarse", c1100.exact_hpwl, c1100.exact_overflow,
        "H219 coarse row 1100")
    f182 = fine.loc[fine.iteration == 182].iloc[0]
    add("512-grid", "H221 selector", f182.exact_hpwl, f182.exact_overflow,
        "H221 fine row 182; restored at row 350")
    for _, row in r252.iterrows():
        add("cap 15%", f"H252 s{int(row.sweep) + 1}", row.exact_hpwl, row.overflow,
            f"H252 recovery sweep {int(row.sweep)}")
    h253 = g253.loc[g253.iteration == 80].iloc[0]
    add("bridge", "H253", h253.exact_hpwl, h253.overflow,
        "H253 global row 80")
    add("retighten", "H254", float(s254["gp_hpwl"]), float(s254["gp_overflow"]),
        "H254 selected row 99")
    for _, row in r255.iterrows():
        add("7% polish", f"H255 s{int(row.sweep) + 1}", row.exact_hpwl, row.overflow,
            f"H255 recovery sweep {int(row.sweep)}")
    for _, row in r257.iterrows():
        add("7% polish", "H257", row.exact_hpwl, row.overflow,
            "H257 recovery sweep 0")
    b = b372.iloc[0]
    add("capacity cut", "H372 assign", b.final_hpwl, b.final_overflow,
        "H372 bisection final")
    for _, row in r372.iterrows():
        add("exact recovery", f"H372 r{int(row.sweep) + 1}", row.exact_hpwl,
            row.overflow, f"H372 recovery10 sweep {int(row.sweep)}")
    for _, row in s375.iterrows():
        add("equal-shape", f"H375 x{int(row.sweep) + 1}", row.exact_hpwl,
            row.overflow, f"H375 exchange sweep {int(row.sweep)}")
    return pd.DataFrame.from_records(records)


def plot_full_chain() -> None:
    chain = full_chain_records()
    chain.to_csv(DATA_OUT / "h375_stage_summary.csv", index=False)
    x = np.arange(len(chain))
    spans = contiguous_spans(chain.group)
    fig, (ax_w, ax_o) = plt.subplots(
        2, 1, figsize=(11.2, 6.6), sharex=True,
        gridspec_kw={"height_ratios": [1.0, 1.0], "hspace": 0.10},
    )
    shade_stages([ax_w, ax_o], spans, label_axis=ax_w, y=0.985)
    ax_w.plot(x, chain.hpwl_m, color=COLORS["blue"], marker="o", zorder=3)
    ax_o.plot(x, chain.overflow_pct, color=COLORS["orange"], marker="o", zorder=3)
    ax_w.axhline(BASELINE_HPWL_M, color=COLORS["red"], linestyle="--", linewidth=1.25,
                 label="DREAMPlace a1 GP baseline: 73.22M")
    ax_o.axhline(7.0, color=COLORS["red"], linestyle="--", linewidth=1.25,
                 label="7% overflow target")
    ax_w.set_ylabel("Exact HPWL (M)")
    ax_o.set_ylabel("Exact overflow (%)")
    ax_o.set_xlabel("Inherited checkpoint / exact operator event")
    ax_w.set_title("H375 complete inherited checkpoint chain")
    ax_w.legend(loc="upper left", frameon=False)
    ax_o.legend(loc="upper left", frameon=False)
    ax_o.set_ylim(-1.0, 104.5)
    ax_w.set_ylim(38.0, 118.0)

    notable = {
        "raw fresh": (7, 10),
        "H219 seed": (8, 10),
        "H221 selector": (8, 10),
        "H257": (-56, -24),
        "H372 assign": (8, 10),
        "H375 x5": (-72, 12),
    }
    for label, offset in notable.items():
        idx = chain.index[chain.label == label][0]
        ax_w.annotate(
            f"{label}\n{chain.loc[idx, 'hpwl_m']:.2f}M",
            (idx, chain.loc[idx, "hpwl_m"]), xytext=offset,
            textcoords="offset points", fontsize=7.2,
            arrowprops={"arrowstyle": "-", "lw": 0.6, "color": "#555555"},
        )
    for label, offset in {"raw fresh": (9, -30), "H221 selector": (8, -34),
                          "H372 assign": (8, 12), "H375 x5": (-78, 14)}.items():
        idx = chain.index[chain.label == label][0]
        ax_o.annotate(
            f"{chain.loc[idx, 'overflow_pct']:.3f}%",
            (idx, chain.loc[idx, "overflow_pct"]), xytext=offset,
            textcoords="offset points", fontsize=7.2,
            arrowprops={"arrowstyle": "-", "lw": 0.6, "color": "#555555"},
        )
    ax_o.set_xticks(x)
    ax_o.set_xticklabels(chain.label, rotation=58, ha="right")
    fig.subplots_adjust(top=0.90, bottom=0.24, left=0.075, right=0.985)
    save_figure(fig, "fig_full_chain_convergence")


def plot_homotopy_detail() -> None:
    coarse = read_csv(EXP / "h219_a1_coarse" / "homotopy_metrics.csv")
    fine = read_csv(EXP / "h221_a1_fine" / "homotopy_metrics.csv")
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 7.2))

    ax = axes[0, 0]
    ax2 = ax.twinx()
    inherited = coarse.iteration <= 1100
    ax.plot(coarse.loc[inherited, "iteration"], coarse.loc[inherited, "exact_hpwl"] / 1e6,
            color=COLORS["blue"], label="exact HPWL")
    ax2.plot(coarse.loc[inherited, "iteration"], 100 * coarse.loc[inherited, "exact_overflow"],
             color=COLORS["orange"], label="exact overflow")
    ax.plot(coarse.loc[~inherited, "iteration"], coarse.loc[~inherited, "exact_hpwl"] / 1e6,
            color=COLORS["gray"], linestyle="--", alpha=0.8)
    ax2.plot(coarse.loc[~inherited, "iteration"], 100 * coarse.loc[~inherited, "exact_overflow"],
             color=COLORS["lightgray"], linestyle="--", alpha=0.95)
    ax.axvline(1100, color=COLORS["red"], linestyle=":", linewidth=1.2)
    ax.annotate("inherited row 1100", (1100, 60.62), xytext=(-110, -30),
                textcoords="offset points", fontsize=7.5,
                arrowprops={"arrowstyle": "->", "lw": 0.7})
    ax.set_title("H219 coarse electrostatic homotopy (128×128)")
    ax.set_xlabel("iteration")
    ax.set_ylabel("Exact HPWL (M)", color=COLORS["blue"])
    ax2.set_ylabel("Exact overflow (%)", color=COLORS["orange"])
    ax.tick_params(axis="y", colors=COLORS["blue"])
    ax2.tick_params(axis="y", colors=COLORS["orange"])
    ax2.grid(False)

    ax = axes[0, 1]
    ax2 = ax.twinx()
    selected = fine.iteration <= 182
    middle = (fine.iteration > 182) & (fine.iteration < 350)
    rejected = fine.iteration >= 350
    ax.plot(fine.loc[selected, "iteration"], fine.loc[selected, "exact_hpwl"] / 1e6,
            color=COLORS["blue"], label="inherited positive-$\\mu$ branch")
    ax2.plot(fine.loc[selected, "iteration"], 100 * fine.loc[selected, "exact_overflow"],
             color=COLORS["orange"])
    ax.plot(fine.loc[middle, "iteration"], fine.loc[middle, "exact_hpwl"] / 1e6,
            color=COLORS["gray"], linestyle="--", alpha=0.75)
    ax2.plot(fine.loc[middle, "iteration"], 100 * fine.loc[middle, "exact_overflow"],
             color=COLORS["lightgray"], linestyle="--", alpha=0.95)
    ax.plot(fine.loc[rejected, "iteration"], fine.loc[rejected, "exact_hpwl"] / 1e6,
            color=COLORS["purple"], linestyle=":", label="logged $\\mu=0$ trial")
    ax2.plot(fine.loc[rejected, "iteration"], 100 * fine.loc[rejected, "exact_overflow"],
             color=COLORS["red"], linestyle=":")
    for value, text_value in [(60, "$\\lambda$ activates"), (182, "best overflow"),
                              (350, "restore; $\\mu=0$")]:
        ax.axvline(value, color="#555555", linestyle=":", linewidth=0.9)
        ax.text(value + 3, 118.5, text_value, rotation=90, va="top", fontsize=7.0)
    ax.scatter([182, 350], [109.425511, 109.425511], color=COLORS["red"], zorder=5)
    ax.set_title("H221 fine homotopy and selector semantics (512×512)")
    ax.set_xlabel("iteration")
    ax.set_ylabel("Exact HPWL (M)", color=COLORS["blue"])
    ax2.set_ylabel("Exact overflow (%)", color=COLORS["orange"])
    ax.tick_params(axis="y", colors=COLORS["blue"])
    ax2.tick_params(axis="y", colors=COLORS["orange"])
    ax2.grid(False)
    ax.legend(loc="lower left", frameon=False)

    ax = axes[1, 0]
    ax2 = ax.twinx()
    ax.plot(fine.iteration, fine.mu_control, color=COLORS["purple"], label="$\\mu$ control")
    ax2.plot(fine.iteration, fine.lambda_control, color=COLORS["green"],
             label="$\\lambda$ control")
    ax.axvspan(0, 59, color="#F4E7CD", alpha=0.65, zorder=-5)
    ax.text(29.5, 0.93, "$\\lambda=0$", transform=ax.get_xaxis_transform(),
            ha="center", va="top", fontsize=8)
    ax.axvline(182, color="#555555", linestyle=":", linewidth=0.9)
    ax.axvline(350, color="#555555", linestyle=":", linewidth=0.9)
    ax.set_title("H221 adaptive homotopy controls")
    ax.set_xlabel("iteration")
    ax.set_ylabel("$\\mu$ control", color=COLORS["purple"])
    ax2.set_ylabel("$\\lambda$ control", color=COLORS["green"])
    ax.tick_params(axis="y", colors=COLORS["purple"])
    ax2.tick_params(axis="y", colors=COLORS["green"])
    ax2.grid(False)
    handles = [mpl.lines.Line2D([], [], color=COLORS["purple"], label="$\\mu$ control"),
               mpl.lines.Line2D([], [], color=COLORS["green"], label="$\\lambda$ control")]
    ax.legend(handles=handles, loc="center right", frameon=False)

    ax = axes[1, 1]
    for column, label, color in [
        ("hpwl_term_norm", "HPWL term", COLORS["blue"]),
        ("overlap_term_norm", "exact-overlap term", COLORS["green"]),
        ("electro_term_norm", "electrostatic term", COLORS["purple"]),
    ]:
        values = fine[column].replace(0.0, np.nan)
        ax.plot(fine.iteration, values, label=label, color=color)
    ax.axvspan(0, 59, color="#F4E7CD", alpha=0.65, zorder=-5)
    ax.axvline(182, color="#555555", linestyle=":", linewidth=0.9)
    ax.axvline(350, color="#555555", linestyle=":", linewidth=0.9)
    ax.set_yscale("log")
    ax.set_title("H221 normalized gradient-term norms")
    ax.set_xlabel("iteration")
    ax.set_ylabel("$\\ell_2$ norm (log scale)")
    ax.legend(loc="best", frameon=False)
    ax.text(30, 1.6e-6, "overflow falls with\noverlap term disabled",
            ha="center", va="center", fontsize=7.3)

    fig.suptitle("Coarse-to-fine homotopy: logged trajectories, handoffs, and force attribution",
                 fontsize=12, y=1.005)
    fig.tight_layout()
    save_figure(fig, "fig_homotopy_detail")


def plot_exact_tail() -> None:
    s252 = read_kv(EXP / "h252_a1_relaxed2_net_recovery" / "summary.txt")
    r252 = read_csv(EXP / "h252_a1_relaxed2_net_recovery" / "recovery_metrics.csv")
    g253 = read_csv(EXP / "h253_a1_relaxed2_bridge80" / "global_metrics.csv")
    g254 = read_csv(EXP / "h254_a1_relaxed2_strong475" / "global_metrics.csv")
    s254 = read_kv(EXP / "h254_a1_relaxed2_strong475" / "summary.txt")
    r255 = read_csv(EXP / "h255_a1_relaxed2_four_node_polish" / "recovery_metrics.csv")
    r257 = read_csv(EXP / "h257_a1_fifth_node_line4" / "recovery_metrics.csv")
    b372 = read_csv(
        EXP / "h372_a1_surplus_only_no_axis_map" / "bisect" / "bisection_metrics.csv"
    ).iloc[0]
    r372 = read_csv(
        EXP / "h372_a1_surplus_only_no_axis_map" / "recovery10" / "recovery_metrics.csv"
    )
    x375 = read_csv(
        EXP / "h375_a1_surplus_recovery10_exchange" / "run" / "swap_recovery_metrics.csv"
    )

    points: list[dict[str, object]] = []

    def append(stage: str, hpwl: float, overflow: float, label: str) -> None:
        points.append({"stage": stage, "hpwl_m": hpwl / 1e6,
                       "overflow_pct": 100 * overflow, "label": label})

    append("H252 cap15", float(s252["recovery_initial_hpwl"]),
           float(s252["recovery_initial_overflow"]), "audit")
    for _, row in r252.iterrows():
        append("H252 cap15", row.exact_hpwl, row.overflow, f"s{int(row.sweep)+1}")
    for _, row in g253.loc[g253.iteration <= 80].iterrows():
        append("H253 bridge", row.exact_hpwl, row.overflow, f"i{int(row.iteration)}")
    for _, row in g254.loc[g254.iteration <= 99].iterrows():
        append("H254 retighten", row.exact_hpwl, row.overflow, f"i{int(row.iteration)}")
    append("H255/257", float(s254["gp_hpwl"]), float(s254["gp_overflow"]),
           "H254 out")
    for _, row in r255.iterrows():
        append("H255/257", row.exact_hpwl, row.overflow,
               f"H255 s{int(row.sweep)+1}")
    for _, row in r257.iterrows():
        append("H255/257", row.exact_hpwl, row.overflow, "H257")
    append("H372 assign", b372.final_hpwl, b372.final_overflow, "assignment")
    for _, row in r372.iterrows():
        append("H372 rec.", row.exact_hpwl, row.overflow, f"r{int(row.sweep)+1}")
    for _, row in x375.iterrows():
        append("H375 swap", row.exact_hpwl, row.overflow, f"x{int(row.sweep)+1}")

    tail = pd.DataFrame(points)
    x = np.arange(len(tail))
    spans = contiguous_spans(tail.stage)
    fig, (ax_w, ax_o) = plt.subplots(2, 1, figsize=(11.2, 6.3), sharex=True,
                                     gridspec_kw={"hspace": 0.10})
    shade_stages([ax_w, ax_o], spans, label_axis=ax_w, y=0.985)
    ax_w.plot(x, tail.hpwl_m, color=COLORS["blue"])
    ax_o.plot(x, tail.overflow_pct, color=COLORS["orange"])
    ax_w.axhline(BASELINE_HPWL_M, color=COLORS["red"], linestyle="--", linewidth=1.1)
    ax_o.axhline(7.0, color=COLORS["red"], linestyle="--", linewidth=1.1)
    ax_o.axhline(15.0, color=COLORS["purple"], linestyle=":", linewidth=1.1)
    ax_w.set_ylabel("Exact HPWL (M)")
    ax_o.set_ylabel("Exact overflow (%)")
    ax_o.set_xlabel("Inherited exact-tail update index")
    ax_w.set_title("H252–H375 inherited exact-tail trajectory")
    ax_o.set_ylim(0.0, 31.0)
    ax_w.set_ylim(70.0, 112.0)

    key_labels = ["audit", "s2", "i80", "H254 out", "H257", "assignment", "r10", "x5"]
    used: set[str] = set()
    offsets = {
        "audit": (5, 10), "s2": (-28, -28), "i80": (-18, -30),
        "H254 out": (6, 10), "H257": (-34, -28), "assignment": (7, 10),
        "r10": (-30, -28), "x5": (-30, 10),
    }
    for idx, row in tail.iterrows():
        if row.label not in key_labels or row.label in used:
            continue
        used.add(str(row.label))
        ax_w.scatter(idx, row.hpwl_m, color=COLORS["red"], zorder=4)
        ax_w.annotate(f"{row.label}: {row.hpwl_m:.2f}M", (idx, row.hpwl_m),
                      xytext=offsets[str(row.label)], textcoords="offset points",
                      fontsize=7.1, arrowprops={"arrowstyle": "-", "lw": 0.55})
    ax_o.annotate("15% topology-recovery band", (1.5, 15.0), xytext=(10, 35),
                  textcoords="offset points", fontsize=7.5,
                  arrowprops={"arrowstyle": "->", "lw": 0.65})
    ax_o.annotate("surplus-only capacity assignment", (len(tail) - 16, 2.306),
                  xytext=(-95, 38), textcoords="offset points", fontsize=7.5,
                  arrowprops={"arrowstyle": "->", "lw": 0.65})
    fig.subplots_adjust(top=0.90, bottom=0.10, left=0.075, right=0.985)
    save_figure(fig, "fig_exact_tail")


def plot_runtime_breakdown() -> None:
    stages = [
        ("H219 HPWL seed", 5.8001669, "gradient"),
        ("H219 coarse→i1100", 64.4396803, "homotopy"),
        ("H221 fine→restore", 19.8891736, "homotopy"),
        ("H252 cap15", 8.7395387, "exact"),
        ("H253 bridge", 4.9110000, "gradient"),
        ("H254 retighten", 7.6143017, "gradient"),
        ("H255 polish", 19.0895167, "exact"),
        ("H257 polish", 3.2884923, "exact"),
        ("H372 assignment", 2.0589292, "assignment"),
        ("H372 recovery×10", 62.9419250, "exact"),
        ("H375 exchange×5", 37.6254400, "exchange"),
    ]
    runtime = pd.DataFrame(stages, columns=["stage", "seconds", "kind"])
    runtime["cumulative_seconds"] = runtime.seconds.cumsum()
    runtime.to_csv(DATA_OUT / "h375_runtime_breakdown.csv", index=False)
    h257_prefix = runtime.loc[:7, "seconds"].sum()
    post_h257 = runtime.loc[8:, "seconds"].sum()
    full = runtime.seconds.sum()

    kind_color = {
        "gradient": COLORS["blue"],
        "homotopy": COLORS["purple"],
        "exact": COLORS["green"],
        "assignment": COLORS["orange"],
        "exchange": COLORS["sky"],
    }
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(11.2, 6.6),
                                  gridspec_kw={"height_ratios": [1.15, 1.0],
                                               "hspace": 0.76})
    left = 0.0
    for _, row in runtime.iterrows():
        ax.barh([0], [row.seconds], left=left, height=0.48,
                color=kind_color[row.kind], edgecolor="white", linewidth=0.7)
        if row.seconds >= 4.5:
            ax.text(left + 0.5 * row.seconds, 0, f"{row.seconds:.1f}",
                    ha="center", va="center", color="white", fontsize=7.1,
                    fontweight="bold")
        left += row.seconds
    ax.axvline(BASELINE_RUNTIME_S, color="#333333", linestyle=":", linewidth=1.1)
    ax.axvline(RUNTIME_GATE_S, color=COLORS["red"], linestyle="--", linewidth=1.2)
    ax.axvline(h257_prefix, color="#333333", linestyle="-.", linewidth=1.1)
    ax.text(BASELINE_RUNTIME_S, 0.43, "paper GP 67 s", ha="center", va="bottom", fontsize=7.5)
    ax.text(RUNTIME_GATE_S, -0.43, "2× gate 134 s", ha="center", va="top",
            fontsize=7.5, color=COLORS["red"])
    ax.text(h257_prefix, 0.43, f"H257 prefix {h257_prefix:.2f} s", ha="center",
            va="bottom", fontsize=7.5)
    ax.text(full, 0.0, f"  full {full:.2f} s", ha="left", va="center", fontsize=8.2)
    ax.set_xlim(0, 252)
    ax.set_yticks([])
    ax.set_xlabel("Checkpoint-handoff wall-clock accounting (s)")
    ax.set_title("H375 stage-runtime waterfall")
    legend = [Patch(facecolor=kind_color[k], label=l) for k, l in [
        ("gradient", "Adam / gradient"), ("homotopy", "electrostatic homotopy"),
        ("exact", "exact candidate recovery"), ("assignment", "capacity assignment"),
        ("exchange", "equal-shape exchange")]]
    fig.legend(handles=legend, ncol=5, loc="center", bbox_to_anchor=(0.5, 0.505),
               frameon=False)

    labels = ["DREAMPlace\n40T GP", "2× runtime\ngate", "H257\nprefix",
              "post-H257\nonly", "H375\nfull chain"]
    values = [BASELINE_RUNTIME_S, RUNTIME_GATE_S, h257_prefix, post_h257, full]
    colors = [COLORS["gray"], COLORS["lightgray"], COLORS["orange"],
              COLORS["green"], COLORS["red"]]
    bars = ax2.bar(labels, values, color=colors, width=0.66)
    ax2.axhline(RUNTIME_GATE_S, color=COLORS["red"], linestyle="--", linewidth=1.0)
    ax2.set_ylabel("seconds")
    ax2.set_ylim(0, 255)
    ax2.set_title("Why the historical 103 s claim is not an end-to-end runtime", pad=10)
    for bar, value in zip(bars, values):
        ax2.text(bar.get_x() + bar.get_width() / 2, value + 5, f"{value:.2f}",
                 ha="center", va="bottom", fontsize=8.0)
    ax2.text(4.48, 18, "H221 completed run: 22.843 s;\nplot uses row-350 handoff: 19.889 s",
             ha="right", va="bottom", fontsize=7.2, color="#444444")
    fig.subplots_adjust(top=0.94, bottom=0.085, left=0.075, right=0.97)
    save_figure(fig, "fig_runtime_breakdown")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    DATA_OUT.mkdir(parents=True, exist_ok=True)
    plot_full_chain()
    plot_homotopy_detail()
    plot_exact_tail()
    plot_runtime_breakdown()
    print(f"Generated H375 report figures in {OUT}")


if __name__ == "__main__":
    main()
