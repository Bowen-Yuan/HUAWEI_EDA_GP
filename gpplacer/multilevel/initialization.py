"""Deterministic multilevel initialisation driven by net connectivity.

The hierarchy is used only to remove gross initial disorder.  Fine placement is
always decided by the exact oracle in later stages.
"""

from __future__ import annotations

import numpy as np

from gpplacer.model.types import PlacementDB
from gpplacer.oracle.grid import DensityGrid


def capacity_aware_seed(db: PlacementDB, grid: DensityGrid, seed: int) -> np.ndarray:
    """Place movable centres throughout usable capacity with deterministic jitter."""
    centres = db.initial_centres.copy()
    rng = np.random.default_rng(seed)
    candidates = np.flatnonzero(grid.available_area.ravel() > 0.0)
    probabilities = grid.available_area.ravel()[candidates]
    probabilities = probabilities / probabilities.sum()
    bins = rng.choice(candidates, size=int(db.movable.sum()), p=probabilities)
    movable_nodes = np.flatnonzero(db.movable)
    iy, ix = np.divmod(bins, grid.nx)
    x = grid.x_edges[ix] + rng.random(len(ix)) * (grid.x_edges[ix + 1] - grid.x_edges[ix])
    y = grid.y_edges[iy] + rng.random(len(iy)) * (grid.y_edges[iy + 1] - grid.y_edges[iy])
    centres[movable_nodes, 0] = x
    centres[movable_nodes, 1] = y
    return centres


def connectivity_aware_seed(db: PlacementDB, rounds: int) -> np.ndarray:
    """Create a placement seed by propagating physical terminal anchors on nets.

    Movable cells never inherit their coordinate from the input ``.pl``.  They
    start at the core centre, while fixed terminals keep their physical
    locations.  Each round replaces a movable cell by a damped average of the
    centroids of its incident nets.  Consequently, cells sharing nets converge
    toward the same connectivity-defined cluster, and clusters connected to
    fixed terminals acquire an appropriate absolute location.

    Components without a fixed-terminal path remain near the core centre.  A
    later capacity-spreading stage is responsible for assigning such clusters
    distinct legal regions; this function deliberately does not introduce a
    random coordinate to break that symmetry.
    """
    if rounds < 0:
        raise ValueError("initialisation rounds must be non-negative")
    left, bottom, right, top = db.core_bounds
    centres = np.empty_like(db.initial_centres)
    centres[:, 0] = (left + right) / 2.0
    centres[:, 1] = (bottom + top) / 2.0
    centres[db.fixed] = db.initial_centres[db.fixed]
    net_degree = np.diff(db.net_start).astype(np.float64)
    pin_net = np.repeat(np.arange(db.net_count, dtype=np.int64), net_degree.astype(np.int64))
    node_degree = np.bincount(db.pin_node, minlength=db.node_count).astype(np.float64)
    movable = db.movable
    for _ in range(rounds):
        pin_positions = centres[db.pin_node] + db.pin_offset
        net_sum = np.add.reduceat(pin_positions, db.net_start[:-1], axis=0)
        net_centres = net_sum / net_degree[:, None]
        node_sum = np.zeros_like(centres)
        np.add.at(node_sum, db.pin_node, net_centres[pin_net])
        targets = node_sum / np.maximum(node_degree[:, None], 1.0)
        # Damping prevents a high-degree hyperedge from erasing all local
        # cluster structure in a single propagation round.
        centres[movable] = 0.35 * centres[movable] + 0.65 * targets[movable]
        centres[movable, 0] = np.clip(centres[movable, 0], left, right)
        centres[movable, 1] = np.clip(centres[movable, 1], bottom, top)
    return centres


def cluster_capacity_seed(
    db: PlacementDB, grid: DensityGrid, targets: np.ndarray, macro_bins: int,
) -> np.ndarray:
    """Assign heavy-edge clusters to capacity-limited macro and density bins.

    ``targets`` comes from connectivity propagation, not movable-cell `.pl`
    coordinates. Each matched cluster is kept intact during assignment and its
    members receive deterministic local slots in the selected density bin.
    """
    if macro_bins < 1:
        raise ValueError("cluster_macro_bins must be positive")
    left, bottom, right, top = db.core_bounds
    fine_capacity = grid.rho_target * grid.available_area
    fine_x = 0.5 * (grid.x_edges[:-1] + grid.x_edges[1:])
    fine_y = 0.5 * (grid.y_edges[:-1] + grid.y_edges[1:])
    xx, yy = np.meshgrid(fine_x, fine_y)
    mx = np.clip(((xx - left) / (right - left) * macro_bins).astype(np.int64), 0, macro_bins - 1)
    my = np.clip(((yy - bottom) / (top - bottom) * macro_bins).astype(np.int64), 0, macro_bins - 1)
    macro_of_fine = my * macro_bins + mx
    macro_count = macro_bins * macro_bins
    macro_remaining = np.bincount(
        macro_of_fine.ravel(), weights=fine_capacity.ravel(), minlength=macro_count,
    )
    macro_centres = np.empty((macro_count, 2), dtype=np.float64)
    for macro in range(macro_count):
        iy, ix = divmod(macro, macro_bins)
        macro_centres[macro] = (
            left + (ix + 0.5) * (right - left) / macro_bins,
            bottom + (iy + 0.5) * (top - bottom) / macro_bins,
        )

    raw = greedy_heavy_edge_groups(db)
    _, group_of_node = np.unique(raw, return_inverse=True)
    order = np.argsort(group_of_node, kind="stable")
    sizes = np.bincount(group_of_node)
    starts = np.zeros(len(sizes) + 1, dtype=np.int64)
    starts[1:] = np.cumsum(sizes)
    area = db.width * db.height
    groups: list[tuple[int, np.ndarray, float, np.ndarray]] = []
    for group in range(len(sizes)):
        members = order[starts[group]:starts[group + 1]]
        movable = members[db.movable[members]]
        if len(movable):
            groups.append((group, movable, float(area[movable].sum()), targets[movable].mean(axis=0)))
    groups.sort(key=lambda item: (-item[2], item[0]))

    assigned: list[list[tuple[np.ndarray, float, np.ndarray]]] = [[] for _ in range(macro_count)]
    for _, members, cluster_area, target in groups:
        feasible = macro_remaining >= cluster_area
        if np.any(feasible):
            distance = np.sum((macro_centres - target) ** 2, axis=1)
            distance[~feasible] = np.inf
            macro = int(np.argmin(distance))
        else:
            macro = int(np.argmax(macro_remaining))
        macro_remaining[macro] -= cluster_area
        assigned[macro].append((members, cluster_area, target))

    result = targets.copy()
    fine_remaining = fine_capacity.copy()
    for macro, clusters in enumerate(assigned):
        local_bins = np.flatnonzero(macro_of_fine.ravel() == macro)
        for members, cluster_area, target in clusters:
            remaining = fine_remaining.ravel()[local_bins]
            feasible = remaining >= cluster_area
            if np.any(feasible):
                local_centres = np.column_stack((xx.ravel()[local_bins], yy.ravel()[local_bins]))
                distance = np.sum((local_centres - target) ** 2, axis=1)
                distance[~feasible] = np.inf
                flat_bin = int(local_bins[np.argmin(distance)])
            else:
                flat_bin = int(local_bins[np.argmax(remaining)])
            fine_remaining.ravel()[flat_bin] -= cluster_area
            iy, ix = divmod(flat_bin, grid.nx)
            side = int(np.ceil(np.sqrt(len(members))))
            slots = np.arange(len(members))
            result[members, 0] = xx[iy, ix] + 0.55 * grid.bin_width * (((slots % side) + 0.5) / side - 0.5)
            result[members, 1] = yy[iy, ix] + 0.55 * grid.bin_height * (((slots // side) + 0.5) / side - 0.5)
    return result


def greedy_heavy_edge_groups(db: PlacementDB, degree_limit: int = 64) -> np.ndarray:
    """Return a compact group id per node using deterministic net-local matching.

    This retains the report's heavy-edge coarsening intent without materialising
    a dense hypergraph matrix.  Fixed nodes never join groups.
    """
    group = np.arange(db.node_count, dtype=np.int64)
    claimed = db.fixed.copy()
    next_group = 0
    for net in range(db.net_count):
        start, stop = int(db.net_start[net]), int(db.net_start[net + 1])
        if stop - start < 2 or stop - start > degree_limit:
            continue
        candidates = db.pin_node[start:stop]
        movable = candidates[~db.fixed[candidates]]
        if len(movable) < 2:
            continue
        first = int(movable[0])
        if claimed[first]:
            continue
        partner = -1
        for node in movable[1:]:
            if not claimed[node]:
                partner = int(node)
                break
        if partner >= 0:
            group[first] = next_group
            group[partner] = next_group
            claimed[first] = claimed[partner] = True
            next_group += 1
    for node in np.flatnonzero(~claimed):
        group[node] = next_group
        next_group += 1
    # Fixed-node group ids are unique and appear after movable groups.
    for node in np.flatnonzero(db.fixed):
        group[node] = next_group
        next_group += 1
    return group
