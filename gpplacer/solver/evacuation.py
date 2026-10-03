"""Deterministic capacity evacuation for bins trapped in a density flat region."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from gpplacer.model.types import PlacementDB
from gpplacer.oracle.grid import DensityGrid
from gpplacer.solver.steps import project_centres


@dataclass(frozen=True, slots=True)
class EvacuationStats:
    """Summary of the discrete feasibility correction before Explore."""

    moves: int
    remaining_centre_overflow: float


def _bin_indices(centres: np.ndarray, grid: DensityGrid) -> tuple[np.ndarray, np.ndarray]:
    """Return clipped x/y bin indices for every centre."""
    ix = np.searchsorted(grid.x_edges, centres[:, 0], side="right") - 1
    iy = np.searchsorted(grid.y_edges, centres[:, 1], side="right") - 1
    return np.clip(ix, 0, grid.nx - 1), np.clip(iy, 0, grid.ny - 1)


def _node_net_targets(db: PlacementDB, centres: np.ndarray) -> np.ndarray:
    """Return the mean incident-net pin centroid for each movable node."""
    pin_positions = centres[db.pin_node] + db.pin_offset
    net_sum = np.add.reduceat(pin_positions, db.net_start[:-1], axis=0)
    net_degree = np.diff(db.net_start).astype(np.float64)
    net_centres = net_sum / net_degree[:, None]
    pin_net = np.repeat(np.arange(db.net_count, dtype=np.int64), net_degree.astype(np.int64))
    node_sum = np.zeros_like(centres)
    np.add.at(node_sum, db.pin_node, net_centres[pin_net])
    node_degree = np.bincount(db.pin_node, minlength=db.node_count).astype(np.float64)
    return node_sum / np.maximum(node_degree[:, None], 1.0)


def _nearest_slack_bin(
    ix: int, iy: int, target: np.ndarray, area: float, remaining: np.ndarray,
    grid: DensityGrid, max_radius: int,
) -> tuple[int, int] | None:
    """Find a nearby bin with sufficient centre-assignment capacity.

    The search expands in square rings.  At the first ring with feasible bins,
    a net-target distance breaks ties, preserving connectivity where possible.
    """
    for radius in range(1, max_radius + 1):
        candidates: list[tuple[int, int]] = []
        x0, x1 = max(0, ix - radius), min(grid.nx - 1, ix + radius)
        y0, y1 = max(0, iy - radius), min(grid.ny - 1, iy + radius)
        for candidate_x in range(x0, x1 + 1):
            for candidate_y in (y0, y1):
                if remaining[candidate_y, candidate_x] >= area:
                    candidates.append((candidate_x, candidate_y))
        for candidate_y in range(y0 + 1, y1):
            for candidate_x in (x0, x1):
                if remaining[candidate_y, candidate_x] >= area:
                    candidates.append((candidate_x, candidate_y))
        if candidates:
            return min(
                candidates,
                key=lambda item: (
                    (0.5 * (grid.x_edges[item[0]] + grid.x_edges[item[0] + 1]) - target[0]) ** 2
                    + (0.5 * (grid.y_edges[item[1]] + grid.y_edges[item[1] + 1]) - target[1]) ** 2
                ),
            )
    return None


def evacuate_overfull_bins(
    db: PlacementDB, grid: DensityGrid, centres: np.ndarray, *, rounds: int,
    max_moves_per_round: int, max_radius: int,
) -> tuple[np.ndarray, EvacuationStats]:
    """Move low-connectivity cells out of centre-overfull bins before Explore.

    This is a discrete feasibility correction, not a smoothed density force. It
    uses a centre-bin load approximation to choose moves, then the exact Oracle
    remains the sole authority for objective evaluation and acceptance.
    """
    if rounds < 0 or max_moves_per_round < 0 or max_radius < 1:
        raise ValueError("invalid bin-evacuation configuration")
    result = centres.copy()
    movable = db.movable
    area = db.width * db.height
    capacity = grid.rho_target * grid.available_area
    node_degree = np.diff(db.node_net_start)
    total_moves = 0
    remaining = capacity.copy()
    for _ in range(rounds):
        ix, iy = _bin_indices(result, grid)
        flat_bins = iy * grid.nx + ix
        load = np.bincount(flat_bins[movable], weights=area[movable],
                           minlength=grid.nx * grid.ny).reshape(grid.ny, grid.nx)
        remaining = capacity - load
        overflow = np.maximum(-remaining, 0.0)
        source_bins = np.flatnonzero(overflow.ravel() > 0.0)
        if len(source_bins) == 0:
            break
        targets = _node_net_targets(db, result)
        moves_this_round = 0
        for source in source_bins[np.argsort(overflow.ravel()[source_bins])[::-1]]:
            if moves_this_round >= max_moves_per_round:
                break
            source_y, source_x = divmod(int(source), grid.nx)
            members = np.flatnonzero(movable & (flat_bins == source))
            # Moving a low-degree cell first normally harms fewer nets.
            members = members[np.argsort(node_degree[members], kind="stable")]
            for node in members:
                if remaining[source_y, source_x] >= 0.0:
                    break
                if moves_this_round >= max_moves_per_round:
                    break
                destination = _nearest_slack_bin(
                    source_x, source_y, targets[node], float(area[node]), remaining,
                    grid, max_radius,
                )
                if destination is None:
                    continue
                destination_x, destination_y = destination
                result[node, 0] = 0.5 * (grid.x_edges[destination_x] + grid.x_edges[destination_x + 1])
                result[node, 1] = 0.5 * (grid.y_edges[destination_y] + grid.y_edges[destination_y + 1])
                remaining[source_y, source_x] += area[node]
                remaining[destination_y, destination_x] -= area[node]
                moves_this_round += 1
        total_moves += moves_this_round
        if moves_this_round == 0:
            break
    result = project_centres(db, result)
    remaining_overflow = float(np.maximum(-remaining, 0.0).sum())
    return result, EvacuationStats(total_moves, remaining_overflow)
