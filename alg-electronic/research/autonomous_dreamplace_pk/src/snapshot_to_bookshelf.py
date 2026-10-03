#!/usr/bin/env python3
"""Convert a placer visualization snapshot to a Bookshelf placement file."""

from __future__ import annotations

import argparse
import struct
from pathlib import Path


def load_nodes(path: Path) -> list[tuple[str, float, float, bool]]:
    nodes: list[tuple[str, float, float, bool]] = []
    for raw in path.read_text(encoding="ascii", errors="ignore").splitlines():
        fields = raw.split()
        if not fields or not fields[0].startswith("o") or len(fields) < 3:
            continue
        nodes.append((fields[0], float(fields[1]), float(fields[2]),
                      len(fields) >= 4 and fields[3].startswith("terminal")))
    return nodes


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--snapshot", required=True, type=Path)
    parser.add_argument("--nodes", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    nodes = load_nodes(args.nodes)
    raw = args.snapshot.read_bytes()
    expected = len(nodes) * 2 * 4
    if len(raw) != expected:
        raise ValueError(f"snapshot has {len(raw)} bytes; expected {expected}")
    coordinates = struct.iter_unpack("<ff", raw)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="ascii", newline="\n") as handle:
        handle.write("UCLA pl 1.0\n")
        handle.write("# Replayed pre-legal visualization snapshot\n")
        for (name, width, height, fixed), (cx, cy) in zip(nodes, coordinates):
            suffix = " /FIXED" if fixed else ""
            handle.write(
                f"{name}\t{cx - width / 2:.6f}\t{cy - height / 2:.6f}"
                f"\t: N{suffix}\n"
            )


if __name__ == "__main__":
    main()
