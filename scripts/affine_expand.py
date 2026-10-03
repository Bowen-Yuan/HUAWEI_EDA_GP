#!/usr/bin/env python3
"""Apply a common affine expansion to movable Bookshelf node centers."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


def fields(line: str):
    return line.split("#", 1)[0].split()


def read_nodes(path: Path):
    result = {}
    for line in path.read_text(encoding="ascii", errors="ignore").splitlines():
        f = fields(line)
        if len(f) < 3 or f[0] in {"UCLA", "NumNodes", "NumTerminals"}:
            continue
        try:
            result[f[0]] = (float(f[1]), float(f[2]),
                            len(f) >= 4 and f[3].startswith("terminal"))
        except ValueError:
            pass
    return result


def read_core(path: Path):
    xl = float("inf")
    xh = float("-inf")
    yl = float("inf")
    yh = float("-inf")
    y = height = spacing = None
    for line in path.read_text(encoding="ascii", errors="ignore").splitlines():
        if "Coordinate" in line:
            y = float(line.split(":", 1)[1])
        elif "Height" in line:
            height = float(line.split(":", 1)[1])
        elif "Sitespacing" in line:
            spacing = float(line.split(":", 1)[1])
        elif "SubrowOrigin" in line:
            f = line.split()
            origin = float(f[2])
            sites = int(f[5])
            xl = min(xl, origin)
            xh = max(xh, origin + sites * spacing)
            yl = min(yl, y)
            yh = max(yh, y + height)
    if not (xl < xh and yl < yh):
        raise ValueError("could not parse core bounds")
    return xl, yl, xh, yh


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark", required=True, type=Path)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--scale", type=float)
    parser.add_argument("--scale-x", type=float)
    parser.add_argument("--scale-y", type=float)
    args = parser.parse_args()
    if args.scale is None and (args.scale_x is None or args.scale_y is None):
        parser.error("provide --scale or both --scale-x and --scale-y")
    sx = args.scale if args.scale is not None else args.scale_x
    sy = args.scale if args.scale is not None else args.scale_y
    nodes = read_nodes(args.benchmark.with_suffix(".nodes"))
    xl, yl, xh, yh = read_core(args.benchmark.with_suffix(".scl"))
    cx, cy = 0.5 * (xl + xh), 0.5 * (yl + yh)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.input.open(encoding="ascii", errors="ignore") as source, \
            args.output.open("w", encoding="ascii") as target:
        target.write("UCLA pl 1.0\n# affine topology-preserving expansion\n\n")
        for line in source:
            f = fields(line)
            if len(f) < 3 or f[0] in {"UCLA"} or f[0].startswith("#"):
                continue
            name = f[0]
            if name not in nodes:
                continue
            width, height, fixed = nodes[name]
            x, y = float(f[1]), float(f[2])
            if not fixed:
                px = cx + sx * (x + 0.5 * width - cx)
                py = cy + sy * (y + 0.5 * height - cy)
                px = min(max(px, xl + 0.5 * width), xh - 0.5 * width)
                py = min(max(py, yl + 0.5 * height), yh - 0.5 * height)
                x, y = px - 0.5 * width, py - 0.5 * height
            suffix = " : " + " ".join(f[3:]) if len(f) > 3 else ""
            target.write(f"{name}\t{x:.12g}\t{y:.12g}{suffix}\n")


if __name__ == "__main__":
    main()
