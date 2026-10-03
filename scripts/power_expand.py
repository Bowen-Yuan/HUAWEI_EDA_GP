#!/usr/bin/env python3
"""Apply an order-preserving power-law expansion to movable Bookshelf centers."""

from __future__ import annotations

import argparse
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
    xl = float("inf"); xh = float("-inf")
    yl = float("inf"); yh = float("-inf")
    y = height = spacing = None
    for line in path.read_text(encoding="ascii", errors="ignore").splitlines():
        if "Coordinate" in line: y = float(line.split(":", 1)[1])
        elif "Height" in line: height = float(line.split(":", 1)[1])
        elif "Sitespacing" in line: spacing = float(line.split(":", 1)[1])
        elif "SubrowOrigin" in line:
            f = line.split(); origin = float(f[2]); sites = int(f[5])
            xl = min(xl, origin); xh = max(xh, origin + sites * spacing)
            yl = min(yl, y); yh = max(yh, y + height)
    if not (xl < xh and yl < yh): raise ValueError("invalid core bounds")
    return xl, yl, xh, yh


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark", required=True, type=Path)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--exponent", required=True, type=float)
    parser.add_argument("--coverage", type=float, default=0.9)
    args = parser.parse_args()
    if not (0.0 < args.exponent < 1.0 and 0.0 < args.coverage <= 1.0):
        raise ValueError("exponent must be in (0,1), coverage in (0,1]")
    nodes = read_nodes(args.benchmark.with_suffix(".nodes"))
    xl, yl, xh, yh = read_core(args.benchmark.with_suffix(".scl"))
    cx, cy = 0.5 * (xl + xh), 0.5 * (yl + yh)
    hx, hy = 0.5 * (xh - xl), 0.5 * (yh - yl)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.input.open(encoding="ascii", errors="ignore") as source, \
            args.output.open("w", encoding="ascii") as target:
        target.write("UCLA pl 1.0\n# order-preserving power-law expansion\n\n")
        for line in source:
            f = fields(line)
            if len(f) < 3 or f[0] == "UCLA" or f[0].startswith("#"): continue
            name = f[0]
            if name not in nodes: continue
            width, height, fixed = nodes[name]
            x, y = float(f[1]), float(f[2])
            if not fixed:
                dx = x + 0.5 * width - cx; dy = y + 0.5 * height - cy
                nx = cx + (1 if dx >= 0 else -1) * hx * args.coverage * (abs(dx) / hx) ** args.exponent
                ny = cy + (1 if dy >= 0 else -1) * hy * args.coverage * (abs(dy) / hy) ** args.exponent
                nx = min(max(nx, xl + 0.5 * width), xh - 0.5 * width)
                ny = min(max(ny, yl + 0.5 * height), yh - 0.5 * height)
                x, y = nx - 0.5 * width, ny - 0.5 * height
            suffix = " : " + " ".join(f[3:]) if len(f) > 3 else ""
            target.write(f"{name}\t{x:.12g}\t{y:.12g}{suffix}\n")


if __name__ == "__main__": main()
