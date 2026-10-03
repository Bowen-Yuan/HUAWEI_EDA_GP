"""Publication-style exact nonsmooth ablation and active-set plots."""
from pathlib import Path
import csv
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "experiments" / "h403_a1_h375_ablation"
plt.rcParams.update({
    "font.family": "serif", "font.size": 9, "axes.titlesize": 10,
    "axes.labelsize": 9, "legend.fontsize": 7.5, "figure.dpi": 150,
    "savefig.dpi": 300, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.18, "lines.linewidth": 1.6,
})

variants = {
    "A0 no recovery": (86.1896263567, 0.0699999030043, 0.0),
    "A1 node only": (86.1785616573, 0.0699999915872, 26.2469404),
    "A2 breakpoint": (86.1764118998, 0.0699999941276, 52.292086),
    "A3 + compact": (86.1720906014, 0.0699999812065, 107.4269283),
    "A4 + net block": (86.1682593077, 0.0699999992997, 73.4479344),
    "A5 full recovery": (86.1703997074, 0.0699999998065, 129.0573095),
    "A6 + exchange": (85.9826245280, 0.0699999998076, 80.8860393),
}

with (OUT / "ablation_summary.csv").open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["variant", "hpwl_m", "overflow", "stage_seconds"])
    for k, (hp, ov, sec) in variants.items():
        w.writerow([k, hp, ov, sec])

names = list(variants)
hpwl = [variants[k][0] for k in names]
secs = [variants[k][2] for k in names]
fig, ax = plt.subplots(figsize=(6.5, 3.2))
colors = ["#B0BEC5"] * 5 + ["#0072B2", "#D55E00"]
bars = ax.bar(range(len(names)), hpwl, color=colors, edgecolor="white")
ax.set_ylim(85.8, 86.35)
ax.text(0.99, 0.04, "paper GP = 73.22M (off-scale)", transform=ax.transAxes,
        ha="right", va="bottom", fontsize=7, color="#444444")
ax.set_xticks(range(len(names)), ["A0", "A1", "A2", "A3", "A4", "A5", "A6"])
ax.set_ylabel("Final exact HPWL (M)")
ax.set_title("H375 component ablation on adaptec1")
for b, v in zip(bars, hpwl):
    ax.text(b.get_x() + b.get_width()/2, b.get_height() + 0.012, f"{v:.3f}",
            ha="center", va="bottom", fontsize=7, rotation=90)
fig.tight_layout()
fig.savefig(OUT / "ablation_hpwl.pdf", bbox_inches="tight")
fig.savefig(OUT / "ablation_hpwl.png", bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(figsize=(6.5, 3.2))
ax.bar(range(len(names)), secs, color=colors, edgecolor="white")
ax.axhline(134, color="#444444", ls="--", lw=1.2, label="2× a1 paper GP budget")
ax.set_xticks(range(len(names)), ["A0", "A1", "A2", "A3", "A4", "A5", "A6"])
ax.set_ylabel("Stage runtime (s)")
ax.set_title("Ablation runtime")
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(OUT / "ablation_runtime.pdf", bbox_inches="tight")
fig.savefig(OUT / "ablation_runtime.png", bbox_inches="tight")
plt.close(fig)

def curve(path, label, ax1, ax2, color):
    if not path.exists():
        return
    with path.open() as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return
    x = [int(r.get("sweep", r.get("iteration", 0))) + 1 for r in rows]
    y1 = [float(r.get("exact_hpwl", r.get("hpwl", 0))) / 1e6 for r in rows]
    y2 = [float(r.get("overflow", 0)) * 100 for r in rows]
    ax1.plot(x, y1, label=label, color=color, marker="o", markevery=max(1, len(x)//5))
    ax2.plot(x, y2, label=label, color=color, marker="o", markevery=max(1, len(x)//5))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.9))
cols = ["#0072B2", "#009E73", "#E69F00", "#D55E00", "#CC79A7"]
for i, d in enumerate(["A1_node_only", "A2_breakpoint", "A3_breakpoint_compact",
                       "A4_breakpoint_netblock", "A5_full_recovery"]):
    curve(OUT / d / "recovery_metrics.csv", d.replace("_", " "), ax1, ax2, cols[i])
curve(OUT / "A6_exchange_after_A5" / "swap_recovery_metrics.csv",
      "A6 exchange", ax1, ax2, "#D55E00")
ax1.set_xlabel("Recovery sweep"); ax1.set_ylabel("Exact HPWL (M)")
ax2.set_xlabel("Recovery sweep"); ax2.set_ylabel("Exact overflow (%)")
ax1.set_title("HPWL recovery trajectory"); ax2.set_title("Overflow trajectory")
ax1.legend(frameon=False, fontsize=6, loc="best")
fig.tight_layout()
fig.savefig(OUT / "recovery_convergence.pdf", bbox_inches="tight")
fig.savefig(OUT / "recovery_convergence.png", bbox_inches="tight")
plt.close(fig)

active = ROOT / "experiments" / "h404_a1_exact_active_escape" / "raw_breakpoint" / "recovery_metrics.csv"
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.9))
curve(active, "H404 raw exact breakpoint", ax1, ax2, "#D55E00")
for ax in (ax1, ax2): ax.legend(frameon=False)
ax1.set_xlabel("Recovery sweep"); ax1.set_ylabel("Exact HPWL (M)")
ax2.set_xlabel("Recovery sweep"); ax2.set_ylabel("Exact overflow (%)")
ax1.set_title("H404 dead-zone escape"); ax2.set_title("H404 exact overflow")
fig.tight_layout()
fig.savefig(ROOT / "experiments" / "h404_a1_exact_active_escape" / "h404_convergence.pdf", bbox_inches="tight")
fig.savefig(ROOT / "experiments" / "h404_a1_exact_active_escape" / "h404_convergence.png", bbox_inches="tight")
plt.close(fig)
