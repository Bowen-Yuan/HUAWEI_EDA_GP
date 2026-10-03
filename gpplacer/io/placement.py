"""Bookshelf placement writer."""

from __future__ import annotations

from pathlib import Path
import numpy as np

from gpplacer.model.types import PlacementDB


def read_placement_centres(path: str | Path, db: PlacementDB) -> np.ndarray:
    """Read a complete Bookshelf placement file into solver centre coordinates."""
    index = {name: i for i, name in enumerate(db.node_names)}
    lower_left = np.full((db.node_count, 2), np.nan, dtype=np.float64)
    with Path(path).open("r", encoding="ascii") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#") or line.startswith("UCLA"):
                continue
            fields = line.split()
            if len(fields) < 5 or fields[3] != ":":
                raise ValueError(f"Malformed placement record in {path}: {line}")
            if fields[0] not in index:
                raise ValueError(f"Placement file {path} contains unknown node {fields[0]!r}")
            if fields[4] != "N":
                raise ValueError("Only orientation N is currently supported.")
            lower_left[index[fields[0]]] = (float(fields[1]), float(fields[2]))
    if np.isnan(lower_left).any():
        raise ValueError(f"Placement file {path} does not assign every object.")
    return lower_left + np.column_stack((db.width / 2.0, db.height / 2.0))


def write_placement(path: str | Path, db: PlacementDB, centres: np.ndarray) -> None:
    """Write all nodes as a Bookshelf ``.pl`` file using lower-left positions."""
    if centres.shape != (db.node_count, 2):
        raise ValueError(f"Expected centres shape {(db.node_count, 2)}, got {centres.shape}")
    lower_left = centres - np.column_stack((db.width / 2.0, db.height / 2.0))
    lines = ["UCLA pl 1.0", "# Written by gpplacer"]
    for node, name in enumerate(db.node_names):
        suffix = " /FIXED" if db.fixed[node] else ""
        x, y = lower_left[node]
        lines.append(f"{name}\t{x:.6f}\t{y:.6f}\t: {db.orientation[node]}{suffix}")
    Path(path).write_text("\n".join(lines) + "\n", encoding="ascii")
