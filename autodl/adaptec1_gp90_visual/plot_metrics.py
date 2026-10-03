import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent


def parse_args():
    parser = argparse.ArgumentParser(description="Plot RBSM convergence metrics.")
    parser.add_argument(
        "--input",
        type=Path,
        default=ROOT / "output" / "convergence_metrics.csv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "output" / "convergence_curves.png",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    data = np.genfromtxt(args.input, delimiter=",", names=True)
    data = np.atleast_1d(data)
    if data.size == 0:
        raise RuntimeError(f"No metric rows found in {args.input}")

    step = data["global_step"]
    hpwl = data["hpwl_full"]
    overflow = data["overflow_region_pct"]

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    axes[0].plot(step, hpwl, color="#176B87", linewidth=1.4)
    axes[0].set_ylabel("Full HPWL")
    axes[0].set_yscale("log")
    axes[0].grid(True, alpha=0.25)

    axes[1].plot(step, overflow, color="#C4512D", linewidth=1.4)
    axes[1].set_xlabel("Inner optimization step")
    axes[1].set_ylabel("Overflow (%)")
    axes[1].grid(True, alpha=0.25)

    fig.suptitle("RBSM convergence")
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180, bbox_inches="tight")
    print(f"Saved convergence plot: {args.output.resolve()}")


if __name__ == "__main__":
    main()
