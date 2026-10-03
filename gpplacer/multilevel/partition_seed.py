"""Rollback-safe recursive capacity-constrained hypergraph partition seed.

This module is deliberately separate from the existing initialisers.  It uses
heavy-edge clusters, recursively balanced regions, and inexpensive graph-cut
refinement to create a legal-capacity-oriented starting point without reading
movable coordinates from a Bookshelf ``.pl`` file.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from gpplacer.model.types import PlacementDB
from gpplacer.multilevel.initialization import greedy_heavy_edge_groups
from gpplacer.oracle.grid import DensityGrid


@dataclass(slots=True)
class _ClusterData:
    members: list[np.ndarray]
    area: np.ndarray
    target: np.ndarray
    adjacency: list[list[tuple[int, float]]]


def _clusters(db: PlacementDB, targets: np.ndarray, degree_limit: int) -> _ClusterData:
    """Build matched clusters and a sparse star expansion of the net hypergraph."""
    raw = greedy_heavy_edge_groups(db, degree_limit)
    _, group_of_node = np.unique(raw, return_inverse=True)
    group_count = int(group_of_node.max()) + 1
    order = np.argsort(group_of_node, kind="stable")
    sizes = np.bincount(group_of_node, minlength=group_count)
    starts = np.zeros(group_count + 1, dtype=np.int64)
    starts[1:] = np.cumsum(sizes)
    members: list[np.ndarray] = []
    area = np.zeros(group_count, dtype=np.float64)
    target = np.zeros((group_count, 2), dtype=np.float64)
    cell_area = db.width * db.height
    for group in range(group_count):
        group_nodes = order[starts[group]:starts[group + 1]]
        movable = group_nodes[db.movable[group_nodes]]
        members.append(movable)
        if len(movable):
            area[group] = float(cell_area[movable].sum())
            target[group] = targets[movable].mean(axis=0)
    adjacency: list[list[tuple[int, float]]] = [[] for _ in range(group_count)]
    for net in range(db.net_count):
        pins = db.pin_node[db.net_start[net]:db.net_start[net + 1]]
        groups = np.unique(group_of_node[pins[db.movable[pins]]])
        if len(groups) < 2 or len(groups) > degree_limit:
            continue
        root = int(groups[0])
        edge_weight = float(db.net_weight[net]) / (len(groups) - 1)
        for other in groups[1:]:
            other = int(other)
            adjacency[root].append((other, edge_weight))
            adjacency[other].append((root, edge_weight))
    return _ClusterData(members, area, target, adjacency)


def _region_capacity(grid: DensityGrid, bounds: tuple[float, float, float, float]) -> float:
    """Return target capacity of fine bins whose centres lie in a region."""
    left, right, bottom, top = bounds
    x = 0.5 * (grid.x_edges[:-1] + grid.x_edges[1:])
    y = 0.5 * (grid.y_edges[:-1] + grid.y_edges[1:])
    mask_x = (x >= left) & (x <= right)
    mask_y = (y >= bottom) & (y <= top)
    return float(grid.rho_target * grid.available_area[np.ix_(mask_y, mask_x)].sum())


def _split_groups(
    groups: np.ndarray, data: _ClusterData, axis: int, split: float,
    left_capacity: float, right_capacity: float, passes: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Capacity-balanced seed followed by bounded graph-cut improvement."""
    projection = data.target[groups, axis]
    order = groups[np.lexsort((groups, projection))]
    side = np.zeros(len(order), dtype=np.int8)
    index = {int(group): position for position, group in enumerate(order)}
    left_area = 0.0
    right_area = 0.0
    for position, group in enumerate(order):
        group_area = data.area[group]
        prefer_left = data.target[group, axis] <= split
        if prefer_left and left_area + group_area <= left_capacity:
            side[position] = 0; left_area += group_area
        elif right_area + group_area <= right_capacity:
            side[position] = 1; right_area += group_area
        elif left_area <= right_area:
            side[position] = 0; left_area += group_area
        else:
            side[position] = 1; right_area += group_area
    for _ in range(passes):
        changed = False
        for position, group in enumerate(order):
            current = int(side[position]); other = 1 - current
            group_area = data.area[group]
            if other == 0 and left_area + group_area > left_capacity:
                continue
            if other == 1 and right_area + group_area > right_capacity:
                continue
            same_weight = 0.0; other_weight = 0.0
            for neighbour, weight in data.adjacency[group]:
                neighbour_position = index.get(neighbour)
                if neighbour_position is None:
                    continue
                if side[neighbour_position] == current:
                    same_weight += weight
                else:
                    other_weight += weight
            # Moving is beneficial only when it reduces the local cut.
            if other_weight <= same_weight:
                continue
            if current == 0:
                left_area -= group_area; right_area += group_area
            else:
                right_area -= group_area; left_area += group_area
            side[position] = other
            changed = True
        if not changed:
            break
    return order[side == 0], order[side == 1]


def recursive_partition_seed(
    db: PlacementDB, grid: DensityGrid, targets: np.ndarray, *, leaf_bins: int,
    refine_passes: int, degree_limit: int,
) -> np.ndarray:
    """Place connectivity clusters through recursive balanced region partitioning."""
    if leaf_bins < 2 or leaf_bins & (leaf_bins - 1):
        raise ValueError("partition_leaf_bins must be a power of two >= 2")
    data = _clusters(db, targets, degree_limit)
    active = np.asarray([group for group, members in enumerate(data.members) if len(members)], dtype=np.int64)
    left, bottom, right, top = db.core_bounds
    leaves: list[tuple[np.ndarray, tuple[float, float, float, float]]] = []

    def partition(groups: np.ndarray, bounds: tuple[float, float, float, float], depth: int) -> None:
        region_left, region_right, region_bottom, region_top = bounds
        if depth == int(np.log2(leaf_bins)) or len(groups) <= 1:
            leaves.append((groups, bounds)); return
        split_x = region_right - region_left >= region_top - region_bottom
        if split_x:
            middle = 0.5 * (region_left + region_right)
            left_bounds = (region_left, middle, region_bottom, region_top)
            right_bounds = (middle, region_right, region_bottom, region_top)
        else:
            middle = 0.5 * (region_bottom + region_top)
            left_bounds = (region_left, region_right, region_bottom, middle)
            right_bounds = (region_left, region_right, middle, region_top)
        first, second = _split_groups(
            groups, data, 0 if split_x else 1, middle,
            _region_capacity(grid, left_bounds), _region_capacity(grid, right_bounds), refine_passes,
        )
        partition(first, left_bounds, depth + 1)
        partition(second, right_bounds, depth + 1)

    partition(active, (left, right, bottom, top), 0)
    result = targets.copy()
    fine_x = 0.5 * (grid.x_edges[:-1] + grid.x_edges[1:])
    fine_y = 0.5 * (grid.y_edges[:-1] + grid.y_edges[1:])
    xx, yy = np.meshgrid(fine_x, fine_y)
    remaining = grid.rho_target * grid.available_area
    for groups, (region_left, region_right, region_bottom, region_top) in leaves:
        local = np.flatnonzero(
            ((xx >= region_left) & (xx <= region_right) & (yy >= region_bottom) & (yy <= region_top)).ravel()
        )
        for group in groups:
            group_area = data.area[group]
            slack = remaining.ravel()[local]
            feasible = slack >= group_area
            if np.any(feasible):
                centres = np.column_stack((xx.ravel()[local], yy.ravel()[local]))
                distance = np.sum((centres - data.target[group]) ** 2, axis=1)
                distance[~feasible] = np.inf
                flat = int(local[np.argmin(distance)])
            else:
                flat = int(local[np.argmax(slack)])
            remaining.ravel()[flat] -= group_area
            iy, ix = divmod(flat, grid.nx)
            members = data.members[group]
            side = int(np.ceil(np.sqrt(len(members))))
            slots = np.arange(len(members))
            result[members, 0] = xx[iy, ix] + 0.55 * grid.bin_width * (((slots % side) + 0.5) / side - 0.5)
            result[members, 1] = yy[iy, ix] + 0.55 * grid.bin_height * (((slots // side) + 0.5) / side - 0.5)
    return result
