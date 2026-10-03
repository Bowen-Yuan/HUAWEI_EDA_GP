"""One-level hypergraph coarsening and safe coordinate expansion."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from gpplacer.model.types import PlacementDB
from gpplacer.multilevel.initialization import greedy_heavy_edge_groups


@dataclass(frozen=True, slots=True)
class Coarsening:
    """Mapping between one fine placement database and its coarse database."""

    fine: PlacementDB
    coarse: PlacementDB
    parent_of_fine: np.ndarray


def _node_nets(node_count: int, pin_node: np.ndarray, net_start: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    incident = [[] for _ in range(node_count)]
    for net in range(len(net_start) - 1):
        for node in np.unique(pin_node[net_start[net]:net_start[net + 1]]):
            incident[int(node)].append(net)
    offsets = np.zeros(node_count + 1, dtype=np.int64)
    offsets[1:] = np.cumsum([len(item) for item in incident], dtype=np.int64)
    values = np.empty(int(offsets[-1]), dtype=np.int64)
    for node, items in enumerate(incident):
        values[offsets[node]:offsets[node + 1]] = items
    return offsets, values


def coarsen_once(db: PlacementDB, degree_limit: int = 64) -> Coarsening:
    """Collapse greedily matched movable pairs and reconstruct a sparse hypergraph.

    Coarse dimensions preserve group area while keeping the common row height.
    Nets that entirely collapse inside one group are removed because they cannot
    contribute to coarse HPWL.
    """
    raw_group = greedy_heavy_edge_groups(db, degree_limit)
    _, parent = np.unique(raw_group, return_inverse=True)
    parent = parent.astype(np.int64)
    coarse_count = int(parent.max()) + 1
    # Group ids are arbitrary after greedy matching.  Sorting once makes child
    # access O(number of nodes), whereas repeated ``flatnonzero`` would turn a
    # million-node hierarchy build into an accidental quadratic operation.
    order = np.argsort(parent, kind="stable")
    group_sizes = np.bincount(parent, minlength=coarse_count)
    group_start = np.zeros(coarse_count + 1, dtype=np.int64)
    group_start[1:] = np.cumsum(group_sizes, dtype=np.int64)
    children = [order[group_start[group]:group_start[group + 1]] for group in range(coarse_count)]
    width = np.empty(coarse_count, dtype=np.float64)
    height = np.empty(coarse_count, dtype=np.float64)
    fixed = np.empty(coarse_count, dtype=bool)
    centres = np.empty((coarse_count, 2), dtype=np.float64)
    for coarse, members in enumerate(children):
        area = float(np.sum(db.width[members] * db.height[members]))
        group_height = max(db.row_height, float(np.max(db.height[members])))
        width[coarse] = area / group_height
        height[coarse] = group_height
        fixed[coarse] = bool(np.any(db.fixed[members]))
        centres[coarse] = np.mean(db.initial_centres[members], axis=0)
    net_names: list[str] = []
    net_start = [0]
    pin_node: list[int] = []
    pin_offset: list[tuple[float, float]] = []
    weights: list[float] = []
    for net, name in enumerate(db.net_names):
        members = parent[db.pin_node[db.net_start[net]:db.net_start[net + 1]]]
        unique_members = np.unique(members)
        if len(unique_members) < 2:
            continue
        net_names.append(name)
        weights.append(float(db.net_weight[net]))
        for member in unique_members:
            pin_node.append(int(member))
            pin_offset.append((0.0, 0.0))
        net_start.append(len(pin_node))
    pin_node_array = np.asarray(pin_node, dtype=np.int64)
    net_start_array = np.asarray(net_start, dtype=np.int64)
    node_net_start, node_nets = _node_nets(coarse_count, pin_node_array, net_start_array)
    coarse = PlacementDB(
        source_aux=f"{db.source_aux}#coarse", node_names=tuple(f"g{idx}" for idx in range(coarse_count)),
        net_names=tuple(net_names), width=width, height=height, fixed=fixed,
        initial_centres=centres, orientation=("N",) * coarse_count,
        pin_node=pin_node_array, pin_offset=np.asarray(pin_offset, dtype=np.float64),
        net_start=net_start_array, net_weight=np.asarray(weights, dtype=np.float64),
        node_net_start=node_net_start, node_nets=node_nets, rows=db.rows,
    )
    return Coarsening(fine=db, coarse=coarse, parent_of_fine=parent)


def uncoarsen_centres(coarsening: Coarsening, coarse_centres: np.ndarray, seed: int) -> np.ndarray:
    """Expand parent coordinates with bounded deterministic child perturbations."""
    fine = coarsening.fine
    centres = fine.initial_centres.copy()
    rng = np.random.default_rng(seed)
    parent = coarsening.parent_of_fine
    order = np.argsort(parent, kind="stable")
    group_sizes = np.bincount(parent, minlength=coarsening.coarse.node_count)
    group_start = np.zeros(coarsening.coarse.node_count + 1, dtype=np.int64)
    group_start[1:] = np.cumsum(group_sizes, dtype=np.int64)
    for parent in range(coarsening.coarse.node_count):
        members = order[group_start[parent]:group_start[parent + 1]]
        movable = members[~fine.fixed[members]]
        if len(movable) == 0:
            continue
        jitter = rng.uniform(-0.25 * fine.row_height, 0.25 * fine.row_height, size=(len(movable), 2))
        centres[movable] = coarse_centres[parent] + jitter
    centres[fine.fixed] = fine.initial_centres[fine.fixed]
    return centres
