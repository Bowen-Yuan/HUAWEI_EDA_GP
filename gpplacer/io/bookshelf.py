"""Strict, readable reader for the ISPD2005 Bookshelf row-placement format."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import re

import numpy as np

from gpplacer.model.types import PlacementDB, SiteRow


_COUNT_RE = re.compile(r"^(Num\w+)\s*:\s*(\d+)")


def _content_lines(path: Path):
    """Yield non-empty, non-comment lines without changing their token order."""
    with path.open("r", encoding="ascii") as handle:
        for raw in handle:
            line = raw.strip()
            if line and not line.startswith("#") and not line.startswith("UCLA"):
                yield line


def _aux_files(aux_path: Path) -> dict[str, Path]:
    for line in _content_lines(aux_path):
        if ":" not in line:
            continue
        _, names = line.split(":", 1)
        paths = [aux_path.parent / name for name in names.split()]
        recognised = {path.suffix: path for path in paths}
        required = {".nodes", ".nets", ".wts", ".pl", ".scl"}
        if required <= set(recognised):
            return recognised
    raise ValueError(f"{aux_path} does not declare a RowBasedPlacement input set")


def _read_nodes(path: Path):
    expected_nodes = expected_terminals = None
    names: list[str] = []
    width: list[float] = []
    height: list[float] = []
    terminal: list[bool] = []
    for line in _content_lines(path):
        match = _COUNT_RE.match(line)
        if match:
            if match.group(1) == "NumNodes":
                expected_nodes = int(match.group(2))
            elif match.group(1) == "NumTerminals":
                expected_terminals = int(match.group(2))
            continue
        fields = line.split()
        if len(fields) < 3:
            raise ValueError(f"Malformed node record in {path}: {line}")
        names.append(fields[0])
        width.append(float(fields[1]))
        height.append(float(fields[2]))
        terminal.append(len(fields) >= 4 and fields[3].startswith("terminal"))
    if expected_nodes is not None and expected_nodes != len(names):
        raise ValueError(f"{path}: NumNodes={expected_nodes}, read {len(names)}")
    if expected_terminals is not None and expected_terminals != sum(terminal):
        raise ValueError(f"{path}: NumTerminals={expected_terminals}, read {sum(terminal)}")
    return names, np.asarray(width), np.asarray(height), np.asarray(terminal, dtype=bool)


def _read_placement(path: Path, names: list[str], width: np.ndarray, height: np.ndarray):
    index = {name: i for i, name in enumerate(names)}
    lower_left = np.full((len(names), 2), np.nan, dtype=np.float64)
    orientation = [""] * len(names)
    explicitly_fixed = np.zeros(len(names), dtype=bool)
    for line in _content_lines(path):
        fields = line.split()
        if len(fields) < 5 or fields[3] != ":":
            raise ValueError(f"Malformed placement record in {path}: {line}")
        try:
            node = index[fields[0]]
        except KeyError as exc:
            raise ValueError(f"{path}: placement names unknown object {fields[0]!r}") from exc
        lower_left[node] = (float(fields[1]), float(fields[2]))
        orientation[node] = fields[4]
        explicitly_fixed[node] = any(token.startswith("/FIXED") for token in fields[5:])
    if np.isnan(lower_left).any():
        missing = [names[i] for i in np.flatnonzero(np.isnan(lower_left).any(axis=1))[:5]]
        raise ValueError(f"{path}: missing placements, e.g. {missing}")
    unsupported = sorted({item for item in orientation if item != "N"})
    if unsupported:
        raise ValueError(
            "Only orientation N is currently supported because pin offsets must be rotated "
            f"for other Bookshelf orientations; found {unsupported}."
        )
    centres = lower_left + np.column_stack((width / 2.0, height / 2.0))
    return centres, tuple(orientation), explicitly_fixed


def _read_weights(path: Path) -> dict[str, float]:
    weights: dict[str, float] = {}
    for line in _content_lines(path):
        fields = line.split()
        if len(fields) != 2:
            raise ValueError(f"Malformed net weight record in {path}: {line}")
        weights[fields[0]] = float(fields[1])
    return weights


def _read_nets(path: Path, node_index: dict[str, int], weights: dict[str, float]):
    expected_nets = expected_pins = None
    net_names: list[str] = []
    net_start = [0]
    pin_node: list[int] = []
    pin_offset: list[tuple[float, float]] = []
    lines = iter(_content_lines(path))
    for line in lines:
        match = _COUNT_RE.match(line)
        if match:
            if match.group(1) == "NumNets":
                expected_nets = int(match.group(2))
            elif match.group(1) == "NumPins":
                expected_pins = int(match.group(2))
            continue
        fields = line.split()
        if len(fields) < 4 or fields[0] != "NetDegree" or fields[1] != ":":
            raise ValueError(f"Malformed net header in {path}: {line}")
        degree, net_name = int(fields[2]), fields[3]
        if net_name in weights and weights[net_name] < 0:
            raise ValueError(f"{path}: negative weight for net {net_name}")
        net_names.append(net_name)
        for _ in range(degree):
            try:
                pin_line = next(lines)
            except StopIteration as exc:
                raise ValueError(f"{path}: truncated net {net_name}") from exc
            pin_fields = pin_line.split()
            if len(pin_fields) < 5 or pin_fields[2] != ":":
                raise ValueError(f"Malformed pin record in {path}: {pin_line}")
            try:
                pin_node.append(node_index[pin_fields[0]])
            except KeyError as exc:
                raise ValueError(f"{path}: pin references unknown node {pin_fields[0]!r}") from exc
            pin_offset.append((float(pin_fields[3]), float(pin_fields[4])))
        net_start.append(len(pin_node))
    if expected_nets is not None and expected_nets != len(net_names):
        raise ValueError(f"{path}: NumNets={expected_nets}, read {len(net_names)}")
    if expected_pins is not None and expected_pins != len(pin_node):
        raise ValueError(f"{path}: NumPins={expected_pins}, read {len(pin_node)}")
    net_weight = np.asarray([weights.get(name, 1.0) for name in net_names], dtype=np.float64)
    return (
        tuple(net_names),
        np.asarray(net_start, dtype=np.int64),
        np.asarray(pin_node, dtype=np.int64),
        np.asarray(pin_offset, dtype=np.float64),
        net_weight,
    )


def _read_rows(path: Path) -> tuple[SiteRow, ...]:
    rows: list[SiteRow] = []
    current: dict[str, float] | None = None
    for line in _content_lines(path):
        if line.startswith("NumRows"):
            continue
        if line.startswith("CoreRow"):
            if current is not None:
                raise ValueError(f"{path}: nested CoreRow")
            current = {}
        elif line == "End":
            if current is None:
                raise ValueError(f"{path}: End without CoreRow")
            required = {"Coordinate", "Height", "Sitewidth", "SubrowOrigin", "NumSites"}
            if not required <= set(current):
                raise ValueError(f"{path}: incomplete CoreRow {current}")
            rows.append(
                SiteRow(
                    x=current["SubrowOrigin"], y=current["Coordinate"],
                    width=current["NumSites"] * current["Sitewidth"],
                    height=current["Height"], site_width=current["Sitewidth"],
                )
            )
            current = None
        elif current is not None and ":" in line:
            key, value = (part.strip() for part in line.split(":", 1))
            if key == "SubrowOrigin":
                fields = value.replace(":", " ").split()
                current["SubrowOrigin"] = float(fields[0])
                current["NumSites"] = float(fields[-1])
            elif key in {"Coordinate", "Height", "Sitewidth"}:
                current[key] = float(value)
    if current is not None:
        raise ValueError(f"{path}: unclosed CoreRow")
    if not rows:
        raise ValueError(f"{path}: no legal rows")
    return tuple(rows)


def _build_node_nets(node_count: int, pin_node: np.ndarray, net_start: np.ndarray):
    incident: list[list[int]] = [[] for _ in range(node_count)]
    for net in range(len(net_start) - 1):
        for node in np.unique(pin_node[net_start[net]:net_start[net + 1]]):
            incident[int(node)].append(net)
    offsets = np.zeros(node_count + 1, dtype=np.int64)
    offsets[1:] = np.cumsum([len(items) for items in incident], dtype=np.int64)
    nets = np.empty(int(offsets[-1]), dtype=np.int64)
    for node, items in enumerate(incident):
        nets[offsets[node]:offsets[node + 1]] = items
    return offsets, nets


def load_bookshelf(aux_path: str | Path) -> PlacementDB:
    """Load one Bookshelf RowBasedPlacement case into a validated ``PlacementDB``."""
    aux = Path(aux_path).resolve()
    files = _aux_files(aux)
    names, width, height, terminal = _read_nodes(files[".nodes"])
    centres, orientation, explicitly_fixed = _read_placement(files[".pl"], names, width, height)
    weights = _read_weights(files[".wts"])
    index = {name: i for i, name in enumerate(names)}
    net_names, net_start, pin_node, pin_offset, net_weight = _read_nets(files[".nets"], index, weights)
    rows = _read_rows(files[".scl"])
    node_net_start, node_nets = _build_node_nets(len(names), pin_node, net_start)
    return PlacementDB(
        source_aux=str(aux), node_names=tuple(names), net_names=net_names,
        width=width, height=height, fixed=terminal | explicitly_fixed,
        initial_centres=centres, orientation=orientation, pin_node=pin_node,
        pin_offset=pin_offset, net_start=net_start, net_weight=net_weight,
        node_net_start=node_net_start, node_nets=node_nets, rows=rows,
    )
