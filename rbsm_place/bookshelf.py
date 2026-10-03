"""Minimal, strict reader/writer for ISPD2005 Bookshelf placements.

The package deliberately owns this reader rather than importing the older project.
Coordinates exposed by :class:`Design` are cell centres; Bookshelf .pl uses lower left.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from typing import Optional, Union

import numpy as np


_COUNT = re.compile(r"^(Num\w+)\s*:\s*(\d+)")


@dataclass
class Design:
    aux: Path
    names: tuple[str, ...]
    width: np.ndarray
    height: np.ndarray
    fixed: np.ndarray
    centres: np.ndarray
    pin_node: np.ndarray
    pin_offset: np.ndarray
    net_start: np.ndarray
    net_weight: np.ndarray
    rows: tuple[tuple[float, float, float, float], ...]

    @property
    def movable(self) -> np.ndarray:
        return ~self.fixed

    @property
    def core_bounds(self) -> tuple[float, float, float, float]:
        return (
            min(r[0] for r in self.rows), min(r[1] for r in self.rows),
            max(r[0] + r[2] for r in self.rows), max(r[1] + r[3] for r in self.rows),
        )


def _lines(path: Path):
    with path.open("r", encoding="ascii") as handle:
        for raw in handle:
            line = raw.strip()
            if line and not line.startswith("#") and not line.startswith("UCLA"):
                yield line


def _files(aux: Path) -> dict[str, Path]:
    for line in _lines(aux):
        if ":" not in line:
            continue
        declared = [aux.parent / x for x in line.split(":", 1)[1].split()]
        result = {p.suffix: p for p in declared}
        if {".nodes", ".nets", ".wts", ".pl", ".scl"} <= result.keys():
            return result
    raise ValueError(f"{aux} is not a complete RowBasedPlacement .aux")


def load_bookshelf(aux_path: Union[str, Path]) -> Design:
    aux = Path(aux_path).resolve()
    files = _files(aux)
    expected_nodes = expected_fixed = None
    names: list[str] = []
    width: list[float] = []
    height: list[float] = []
    fixed: list[bool] = []
    for line in _lines(files[".nodes"]):
        m = _COUNT.match(line)
        if m:
            if m.group(1) == "NumNodes": expected_nodes = int(m.group(2))
            if m.group(1) == "NumTerminals": expected_fixed = int(m.group(2))
            continue
        f = line.split()
        if len(f) < 3: raise ValueError(f"Malformed node: {line}")
        names.append(f[0]); width.append(float(f[1])); height.append(float(f[2]))
        fixed.append(len(f) > 3 and f[3].startswith("terminal"))
    if expected_nodes != len(names) or expected_fixed != sum(fixed):
        raise ValueError("Bookshelf node counts do not match records")
    index = {name: i for i, name in enumerate(names)}
    w, h = np.asarray(width, np.float64), np.asarray(height, np.float64)
    ll = np.full((len(names), 2), np.nan, np.float64)
    placement_fixed = np.zeros(len(names), bool)
    for line in _lines(files[".pl"]):
        f = line.split()
        if len(f) < 5 or f[3] != ":": raise ValueError(f"Malformed placement: {line}")
        i = index[f[0]]
        if f[4] != "N": raise ValueError("Only orientation N is supported")
        ll[i] = (float(f[1]), float(f[2]))
        placement_fixed[i] = any(x.startswith("/FIXED") for x in f[5:])
    if np.isnan(ll).any(): raise ValueError("Placement misses objects")
    weights = {f[0]: float(f[1]) for f in (x.split() for x in _lines(files[".wts"])) if len(f) == 2}
    net_start = [0]; pins: list[int] = []; offsets: list[tuple[float, float]] = []; net_weights: list[float] = []
    it = iter(_lines(files[".nets"]))
    for line in it:
        if _COUNT.match(line): continue
        f = line.split()
        if len(f) < 4 or f[:2] != ["NetDegree", ":"]: raise ValueError(f"Malformed net: {line}")
        degree, net = int(f[2]), f[3]
        for _ in range(degree):
            p = next(it).split()
            if len(p) < 5 or p[2] != ":": raise ValueError("Malformed pin")
            pins.append(index[p[0]]); offsets.append((float(p[3]), float(p[4])))
        net_start.append(len(pins)); net_weights.append(weights.get(net, 1.0))
    rows: list[tuple[float, float, float, float]] = []; current: Optional[dict[str, float]] = None
    for line in _lines(files[".scl"]):
        if line.startswith("CoreRow"): current = {}
        elif line == "End" and current is not None:
            rows.append((current["x"], current["y"], current["width"], current["height"])); current = None
        elif current is not None and ":" in line:
            key, value = (x.strip() for x in line.split(":", 1))
            if key == "Coordinate": current["y"] = float(value)
            elif key == "Height": current["height"] = float(value)
            elif key == "Sitewidth": current["site"] = float(value)
            elif key == "SubrowOrigin":
                values = value.replace(":", " ").split(); current["x"] = float(values[0]); current["width"] = float(values[-1]) * current["site"]
    if not rows: raise ValueError("No placement rows")
    return Design(aux, tuple(names), w, h, np.asarray(fixed) | placement_fixed,
                  ll + np.column_stack((w / 2, h / 2)), np.asarray(pins, np.int64),
                  np.asarray(offsets, np.float64), np.asarray(net_start, np.int64),
                  np.asarray(net_weights, np.float64), tuple(rows))


def write_placement(design: Design, centres: np.ndarray, path: Union[str, Path]) -> None:
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as out:
        out.write("UCLA pl 1.0\n\n")
        for i, name in enumerate(design.names):
            x, y = centres[i] - (design.width[i] / 2, design.height[i] / 2)
            suffix = " /FIXED" if design.fixed[i] else ""
            out.write(f"{name}\t{x:.6f}\t{y:.6f}\t: N{suffix}\n")


def read_placement(design: Design, path: Union[str, Path]) -> np.ndarray:
    """Read a complete .pl into centre coordinates without changing fixed objects."""
    index = {name: i for i, name in enumerate(design.names)}
    result = design.centres.copy(); seen = np.zeros(len(design.names), bool)
    for line in _lines(Path(path)):
        f = line.split()
        if len(f) >= 5 and f[0] in index and f[3] == ":":
            i = index[f[0]]; result[i] = (float(f[1]) + design.width[i] / 2, float(f[2]) + design.height[i] / 2); seen[i] = True
    if not seen.all(): raise ValueError(f"{path} does not contain every Bookshelf object")
    result[design.fixed] = design.centres[design.fixed]
    return result
