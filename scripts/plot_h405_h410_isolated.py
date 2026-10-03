"""Publication-quality isolated zero-gradient mechanism ablation plots."""
from pathlib import Path
import csv
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
EXP = ROOT / "experiments"
names = {
    "H405 breakpoint": "h405_a1_breakpoint_only",
    "H406 boundary strip": "h406_a1_boundary_strip_only",
    "H407 dead-zone": "h407_a1_dead_zone_only",
    "H409 finite oracle": "h409_a1_finite_difference_only",
    "H410 exact noise": "h410_a1_controlled_noise_only",
}
colors = ["#0072B2", "#E69F00", "#D55E00", "#009E73", "#CC79A7"]

def read_csv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))

fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.0), constrained_layout=True)
for (label, folder), color in zip(names.items(), colors):
    metric_file = "density_coordinate_metrics.csv" if folder.startswith("h409") else "recovery_metrics.csv"
    rows = read_csv(EXP / folder / metric_file)
    x = [int(r["sweep"]) + 1 for r in rows]
    axes[0].plot(x, [float(r["exact_hpwl"]) / 1e6 for r in rows], marker="o", lw=2,
                 label=label, color=color)
    axes[1].plot(x, [100 * float(r["overflow"]) for r in rows], marker="o", lw=2,
                 label=label, color=color)

rows = read_csv(EXP / "h408_a1_net_share_only" / "global_metrics.csv")
x = [int(r["iteration"]) for r in rows]
axes[0].plot(x, [float(r["exact_hpwl"]) / 1e6 for r in rows], color="#56B4E9",
             lw=2, label="H408 net-share (GP)")
axes[1].plot(x, [100 * float(r["overflow"]) for r in rows], color="#56B4E9",
             lw=2, label="H408 net-share (GP)")

axes[0].set(xlabel="recovery sweep / GP iteration", ylabel="exact HPWL (M)",
            title="Raw a1: exact HPWL")
axes[1].set(xlabel="recovery sweep / GP iteration", ylabel="exact overflow (%)",
            title="Raw a1: exact rectangle-bin overflow")
for ax in axes:
    ax.grid(True, alpha=.25)
    ax.legend(fontsize=8, frameon=False)
fig.savefig(EXP / "isolated_zero_gradient_ablation.png", dpi=220)
fig.savefig(EXP / "isolated_zero_gradient_ablation.pdf")

with (EXP / "isolated_zero_gradient_summary.csv").open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["experiment", "hpwl_m", "overflow", "moves", "seconds"])
    for label, folder in names.items():
        s = {}
        for line in (EXP / folder / "summary.txt").read_text().splitlines():
            if "=" in line:
                k, v = line.split("=", 1); s[k] = v
        moves_key = "density_coordinate_moves" if folder.startswith("h409") else "recovery_moves"
        w.writerow([label, float(s.get("gp_hpwl", 0)) / 1e6,
                    float(s.get("gp_overflow", 0)), int(float(s.get(moves_key, 0))),
                    float(s.get("gp_seconds", 0))])
    s = {line.split("=",1)[0]: line.split("=",1)[1]
         for line in (EXP / "h408_a1_net_share_only" / "summary.txt").read_text().splitlines() if "=" in line}
    w.writerow(["H408 net-share (GP)", float(s["gp_hpwl"]) / 1e6,
                float(s["gp_overflow"]), "", float(s["gp_seconds"])])
