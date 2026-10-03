"""Render the saved RBSM outer-iteration snapshots as a GIF."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
from matplotlib import animation
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--snapshots",
        type=Path,
        default=ROOT / "output" / "layout_snapshots",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "output" / "layout_evolution.gif",
    )
    parser.add_argument("--fps", type=int, default=5)
    parser.add_argument("--dpi", type=int, default=110)
    return parser.parse_args()


def snapshot_paths(snapshot_dir):
    paths = sorted(snapshot_dir.glob("iteration_*.npz"))
    if not paths:
        raise RuntimeError(f"No layout snapshots found in {snapshot_dir}")
    return paths


def padded_limits(values, log_scale=False):
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if log_scale:
        values = values[values > 0]
    if values.size == 0:
        return (1.0, 10.0) if log_scale else (0.0, 1.0)
    low = float(values.min())
    high = float(values.max())
    if low == high:
        pad = max(abs(low) * 0.05, 1.0)
        return max(low - pad, 1e-12) if log_scale else low - pad, high + pad
    if log_scale:
        return low / 1.15, high * 1.15
    pad = 0.08 * (high - low)
    return low - pad, high + pad


def main():
    args = parse_args()
    if args.fps <= 0:
        raise ValueError("fps must be positive")

    paths = snapshot_paths(args.snapshots)
    with np.load(args.snapshots / "metadata.npz") as metadata:
        width = float(metadata["width"])
        height = float(metadata["height"])
        fixed_positions = metadata["fixed_positions"]

    rounds = []
    hpwl = []
    overflow = []
    outside = []
    density_peak = 1.0
    for path in paths:
        with np.load(path) as snapshot:
            rounds.append(int(snapshot["outer_round"]))
            hpwl.append(float(snapshot["hpwl_full"]))
            overflow.append(float(snapshot["overflow_region_pct"]))
            outside.append(float(snapshot["outside_pct"]))
            density_peak = max(
                density_peak,
                float(np.nanpercentile(snapshot["density"], 99.5)),
            )

    rounds = np.asarray(rounds)
    hpwl = np.asarray(hpwl)
    overflow = np.asarray(overflow)
    outside = np.asarray(outside)
    total_rounds = max(int(rounds.max()), 1)

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.2), dpi=args.dpi)
    fig.subplots_adjust(left=0.055, right=0.95, bottom=0.13, top=0.84, wspace=0.34)
    placement_axis, density_axis, metric_axis = axes

    placement_axis.set(
        title="Movable-cell placement",
        xlabel="x",
        ylabel="y",
        xlim=(0, width),
        ylim=(0, height),
        aspect="equal",
    )
    movable_artist = placement_axis.scatter(
        [], [], s=0.5, alpha=0.30, color="#176B87", rasterized=True,
        label="Movable sample",
    )
    if fixed_positions.size:
        placement_axis.scatter(
            fixed_positions[:, 0], fixed_positions[:, 1], s=5.0,
            color="#C4512D", alpha=0.75, marker="x", label="Fixed terminals",
        )
    placement_axis.legend(loc="upper right", markerscale=2.0)

    density_artist = density_axis.imshow(
        np.zeros((64, 64)),
        extent=(0, width, 0, height),
        origin="lower",
        aspect="equal",
        cmap="YlOrRd",
        vmin=0.0,
        vmax=max(density_peak, 1.05),
    )
    density_axis.set(title="Approximate bin density", xlabel="x", ylabel="y")
    fig.colorbar(density_artist, ax=density_axis, label="Density", shrink=0.82)

    overflow_axis = metric_axis.twinx()
    hpwl_line, = metric_axis.plot([], [], color="#176B87", linewidth=1.5, label="Full HPWL")
    overflow_line, = overflow_axis.plot([], [], color="#C4512D", linewidth=1.5, label="Overflow")
    metric_axis.set(
        title="Convergence through current round",
        xlabel="Outer round k",
        ylabel="Full HPWL",
        xlim=(0, total_rounds),
        ylim=padded_limits(hpwl, log_scale=True),
    )
    metric_axis.set_yscale("log")
    overflow_axis.set_ylabel("Overflow (%)")
    overflow_axis.set_ylim(padded_limits(overflow))
    metric_axis.grid(True, alpha=0.25)
    metric_axis.legend(
        [hpwl_line, overflow_line], ["Full HPWL", "Overflow"], loc="best",
    )

    def update(frame_index):
        with np.load(paths[frame_index]) as snapshot:
            movable_artist.set_offsets(snapshot["positions"])
            density_artist.set_data(snapshot["density"])
        hpwl_line.set_data(rounds[: frame_index + 1], hpwl[: frame_index + 1])
        overflow_line.set_data(rounds[: frame_index + 1], overflow[: frame_index + 1])
        fig.suptitle(
            f"RBSM adaptec1 | outer round {rounds[frame_index]}/{total_rounds} | "
            f"HPWL {hpwl[frame_index]:.3e} | overflow {overflow[frame_index]:.3f}% | "
            f"outside {outside[frame_index]:.3f}%",
            fontsize=12,
            fontweight="bold",
        )
        return movable_artist, density_artist, hpwl_line, overflow_line

    movie = animation.FuncAnimation(
        fig, update, frames=len(paths), interval=1000 / args.fps, blit=False,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    movie.save(args.output, writer=animation.PillowWriter(fps=args.fps), dpi=args.dpi)
    plt.close(fig)
    print(f"Saved layout animation: {args.output.resolve()}")


if __name__ == "__main__":
    main()
