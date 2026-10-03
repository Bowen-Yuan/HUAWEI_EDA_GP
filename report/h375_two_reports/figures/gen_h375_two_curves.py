#!/usr/bin/env python3
"""Generate the two compact, data-driven H375 report figures.

The values are exact checkpoint audits copied from the archived H375 chain.
Lines only connect recorded events; no smoothing or interpolation is used for
the optimization model.  Both vector PDF and 300-dpi PNG outputs are written.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


OUT = Path(__file__).resolve().parent

COLORS = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "red": "#D55E00",
    "purple": "#CC79A7",
    "gray": "#6B6B6B",
    "light_blue": "#DDEBF4",
    "light_orange": "#FFF0D8",
    "light_green": "#E5F3EC",
    "light_purple": "#F4EAF2",
}

mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9.0,
        "axes.titlesize": 11.0,
        "axes.labelsize": 9.5,
        "legend.fontsize": 8.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.22,
        "grid.linewidth": 0.6,
        "lines.linewidth": 1.8,
        "lines.markersize": 5.0,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.08,
    }
)


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.pdf")
    fig.savefig(OUT / f"{stem}.png", dpi=300)
    plt.close(fig)


def stage_bands(
    ax: plt.Axes,
    starts: list[int],
    ends: list[int],
    labels: list[str],
    *,
    draw_labels: bool = True,
    rotate_narrow: bool = True,
) -> None:
    fills = [COLORS["light_blue"], COLORS["light_orange"], COLORS["light_green"], COLORS["light_purple"]]
    for i, (start, end, label) in enumerate(zip(starts, ends, labels)):
        ax.axvspan(start - 0.45, end + 0.45, color=fills[i % len(fills)], alpha=0.72, zorder=-10)
        if draw_labels:
            width = end - start + 1
            ax.text(
                (start + end) / 2.0,
                0.985,
                label,
                transform=ax.get_xaxis_transform(),
                ha="center",
                va="top",
                rotation=90 if (width <= 1 and rotate_narrow) else 0,
                fontsize=7.2,
                color="#333333",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.58, "pad": 1.0},
                clip_on=True,
            )


def overflow_figure() -> None:
    labels = ["raw fresh", "H219 seed", "H219 coarse", "H221 fine start", "H221 selector"]
    hpwl = np.array([53.214813, 43.034762, 60.617468, 60.618692, 109.425511])
    overflow = np.array([99.966213, 96.056285, 71.707837, 77.650357, 7.751547])
    x = np.arange(len(labels))

    df = pd.DataFrame({"event": labels, "hpwl_m": hpwl, "overflow_pct": overflow})
    df.to_csv(OUT.parent / "overflow_reduction_curve.csv", index=False)

    fig, (ax_w, ax_o) = plt.subplots(2, 1, figsize=(10.5, 6.2), sharex=True,
                                     gridspec_kw={"hspace": 0.12})
    stage_labels = ["fresh + seed", "H219 coarse", "H221 fine"]
    stage_bands(ax_w, [0, 2, 3], [1, 2, 4], stage_labels, rotate_narrow=False)
    stage_bands(ax_o, [0, 2, 3], [1, 2, 4], stage_labels,
                draw_labels=False)

    ax_w.plot(x, hpwl, color=COLORS["blue"], marker="o", label="exact HPWL")
    ax_w.axhline(73.22, color=COLORS["red"], linestyle="--", linewidth=1.2,
                 label="DREAMPlace a1 GP baseline: 73.22M")
    ax_w.set_ylabel("Exact HPWL (M)")
    ax_w.set_title("Overflow reduction chain: fresh initialization to H221 selector", pad=8)
    ax_w.legend(loc="upper left", frameon=False)
    ax_w.set_ylim(35, 120)

    ax_o.plot(x, overflow, color=COLORS["orange"], marker="o", label="exact overflow")
    ax_o.axhline(7.0, color=COLORS["red"], linestyle="--", linewidth=1.2,
                 label="7% target")
    ax_o.set_ylabel("Exact overflow (%)")
    ax_o.set_xlabel("Inherited checkpoint")
    ax_o.set_ylim(0, 105)
    ax_o.legend(loc="upper right", frameon=False)
    ax_o.set_xticks(x)
    ax_o.set_xticklabels(labels, rotation=28, ha="right")

    for i, (w, o) in enumerate(zip(hpwl, overflow)):
        wo = (7, 9) if i < 3 else (6, -16)
        oo = (7, -18) if i < 3 else (7, 9)
        ax_w.annotate(f"{w:.2f}M", (i, w), xytext=wo, textcoords="offset points",
                      fontsize=7.4, arrowprops={"arrowstyle": "-", "lw": 0.55, "color": "#555"})
        ax_o.annotate(f"{o:.3f}%", (i, o), xytext=oo, textcoords="offset points",
                      fontsize=7.4, arrowprops={"arrowstyle": "-", "lw": 0.55, "color": "#555"})

    ax_o.text(4.02, 10.0, "H221 selector: 7.7515%\n(about 7%; exact tail retightening is separate)",
              fontsize=7.6, color="#444", ha="right", va="bottom")
    fig.subplots_adjust(top=0.91, bottom=0.19, left=0.085, right=0.98)
    save(fig, "overflow_reduction_curve")


def hpwl_figure() -> None:
    rows = [
        ("H252 audit", 109.425511, 7.803537, "H252"),
        ("H252 s1", 97.813945, 14.999745, "H252"),
        ("H252 s2", 94.548507, 14.999997, "H252"),
        ("H253 i80", 86.241998, 13.891233, "H253"),
        ("H254 selected", 89.549340, 6.556721, "H254"),
        ("H255 s1", 88.064828, 6.999984, "H255/H257"),
        ("H255 s2", 87.567244, 6.999999, "H255/H257"),
        ("H255 s3", 87.317935, 6.999977, "H255/H257"),
        ("H255 s4", 87.178874, 7.000000, "H255/H257"),
        ("H257", 87.119465, 6.999958, "H255/H257"),
        ("H372 assign", 93.705492, 2.306362, "H372"),
        ("H372 r1", 89.143747, 6.999985, "H372"),
        ("H372 r2", 87.641415, 6.999963, "H372"),
        ("H372 r3", 86.972538, 6.999991, "H372"),
        ("H372 r4", 86.632343, 6.999971, "H372"),
        ("H372 r5", 86.450273, 6.999985, "H372"),
        ("H372 r6", 86.343968, 6.999999, "H372"),
        ("H372 r7", 86.278326, 7.000000, "H372"),
        ("H372 r8", 86.234833, 7.000000, "H372"),
        ("H372 r9", 86.207181, 6.999993, "H372"),
        ("H372 r10", 86.189626, 6.999990, "H372"),
        ("H375 x1", 86.064268, 6.999990, "H375"),
        ("H375 x2", 86.025509, 6.999990, "H375"),
        ("H375 x3", 86.009906, 6.999990, "H375"),
        ("H375 x4", 86.002560, 6.999990, "H375"),
        ("H375 x5", 85.999318, 6.999990, "H375"),
    ]
    df = pd.DataFrame(rows, columns=["event", "hpwl_m", "overflow_pct", "stage"])
    df.to_csv(OUT.parent / "hpwl_recovery_curve.csv", index=False)
    x = np.arange(len(df))

    fig, (ax_w, ax_o) = plt.subplots(2, 1, figsize=(11.5, 6.7), sharex=True,
                                     gridspec_kw={"height_ratios": [1.15, 0.85], "hspace": 0.12})
    stage_names = ["H252 relaxed", "H253", "H254", "H255/H257 polish", "H372 recovery", "H375 exchange"]
    starts = [0, 3, 4, 5, 10, 21]
    ends = [2, 3, 4, 9, 20, 25]
    stage_bands(ax_w, starts, ends, stage_names)
    stage_bands(ax_o, starts, ends, stage_names, draw_labels=False)

    ax_w.plot(x, df.hpwl_m, color=COLORS["blue"], marker="o", markersize=3.7,
              label="exact HPWL")
    ax_w.axhline(73.22, color=COLORS["red"], linestyle="--", linewidth=1.15,
                 label="DREAMPlace baseline: 73.22M")
    ax_w.set_ylabel("Exact HPWL (M)")
    ax_w.set_title("HPWL recovery in the exact non-smooth tail", pad=8)
    ax_w.set_ylim(78, 113)
    ax_w.legend(loc="upper right", frameon=False)

    ax_o.plot(x, df.overflow_pct, color=COLORS["orange"], marker="o", markersize=3.0,
              label="exact overflow")
    ax_o.axhline(7.0, color=COLORS["red"], linestyle="--", linewidth=1.1, label="7% cap")
    ax_o.axhline(15.0, color=COLORS["purple"], linestyle=":", linewidth=1.1, label="15% relaxed cap")
    ax_o.set_ylabel("Exact overflow (%)")
    ax_o.set_xlabel("Inherited exact operator event")
    ax_o.set_ylim(0, 19)
    ax_o.legend(loc="upper right", ncol=3, frameon=False)
    ax_o.set_xticks(x)
    ax_o.set_xticklabels(df.event, rotation=62, ha="right", fontsize=7.2)

    notable = {
        "H252 s2": (6, -24), "H253 i80": (-18, -26), "H254 selected": (6, 10),
        "H257": (-25, -25), "H372 assign": (6, 10), "H372 r10": (-31, -23),
        "H375 x5": (-35, 9),
    }
    for event, offset in notable.items():
        idx = int(df.index[df.event == event][0])
        ax_w.scatter([idx], [df.loc[idx, "hpwl_m"]], color=COLORS["red"], zorder=5)
        ax_w.annotate(f"{event}\n{df.loc[idx, 'hpwl_m']:.2f}M", (idx, df.loc[idx, "hpwl_m"]),
                      xytext=offset, textcoords="offset points", fontsize=7.1,
                      arrowprops={"arrowstyle": "-", "lw": 0.55, "color": "#555"})
    ax_o.annotate("15% topology-recovery band", (1.7, 15.0), xytext=(8, 18),
                  textcoords="offset points", fontsize=7.4,
                  arrowprops={"arrowstyle": "->", "lw": 0.6})
    ax_w.text(0.01, 0.03, "H254 is retightening: density improves while HPWL temporarily rises.",
              transform=ax_w.transAxes, fontsize=7.5, color="#444", va="bottom")
    fig.subplots_adjust(top=0.90, bottom=0.24, left=0.075, right=0.985)
    save(fig, "hpwl_recovery_curve")


if __name__ == "__main__":
    overflow_figure()
    hpwl_figure()
    print("Generated overflow_reduction_curve and hpwl_recovery_curve (PDF, PNG, CSV).")
