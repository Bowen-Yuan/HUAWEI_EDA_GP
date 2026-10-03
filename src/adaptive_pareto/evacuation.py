"""Local axis-aligned strict-overflow evacuation candidates."""

from __future__ import annotations

import numpy as np

from gpplacer.model.types import PlacementDB
from gpplacer.oracle.evaluation import DualOracleResult
from gpplacer.oracle.grid import DensityGrid
from gpplacer.solver.steps import project_centres


def _centre_bins(centres: np.ndarray, grid: DensityGrid) -> tuple[np.ndarray, np.ndarray]:
    ix = np.searchsorted(grid.x_edges, centres[:, 0], side="right") - 1
    iy = np.searchsorted(grid.y_edges, centres[:, 1], side="right") - 1
    return np.clip(ix, 0, grid.nx - 1), np.clip(iy, 0, grid.ny - 1)


def _net_target(db: PlacementDB, centres: np.ndarray, node: int) -> np.ndarray:
    nets = db.node_nets[db.node_net_start[node]:db.node_net_start[node + 1]]
    if len(nets) == 0:
        return centres[node]
    targets = []
    for net in nets:
        pins = db.pin_node[db.net_start[net]:db.net_start[net + 1]]
        if len(pins):
            targets.append(centres[pins].mean(axis=0))
    return np.mean(np.asarray(targets), axis=0) if targets else centres[node]


def local_evacuation_candidate(
    db: PlacementDB, grid: DensityGrid, centres: np.ndarray, result: DualOracleResult,
    *, max_radius_bins: int, max_moves: int,
) -> tuple[np.ndarray | None, list[dict[str, float | int]]]:
    """Move cells along their nearest cardinal slack path, never to blocked bins.

    This is intentionally a proposal generator.  The caller must evaluate the
    resulting placement with the exact strict oracle and may split/reject it.
    """
    capacity = grid.rho_target * grid.available_area
    slack = capacity - result.strict.occupancy
    active = result.strict.active_bins
    sources = np.flatnonzero(active.ravel())
    if not len(sources):
        return None, []
    ix, iy = _centre_bins(centres, grid)
    flat = iy * grid.nx + ix
    area = db.width * db.height
    remaining = slack.copy()
    candidate = centres.copy()
    moves: list[dict[str, float | int]] = []
    # Prioritise the largest strict excess regions.
    excess = np.maximum(result.strict.occupancy - capacity, 0.0)
    source_order = sources[np.argsort(excess.ravel()[sources])[::-1]]
    for source in source_order:
        if len(moves) >= max_moves:
            break
        sy, sx = divmod(int(source), grid.nx)
        members = np.flatnonzero(db.movable & (flat == source))
        if not len(members):
            continue
        degree = np.diff(db.node_net_start)[members]
        members = members[np.argsort(degree * np.maximum(area[members], 1.0), kind="stable")]
        for node in members:
            if len(moves) >= max_moves:
                break
            targets: list[tuple[float, int, int, int]] = []
            # Search only cardinal rays: this is the short-axis escape policy.
            for dx, dy in ((-1, 0), (1, 0), (0, -1), (0, 1)):
                for radius in range(1, max_radius_bins + 1):
                    tx, ty = sx + dx * radius, sy + dy * radius
                    if tx < 0 or tx >= grid.nx or ty < 0 or ty >= grid.ny:
                        break
                    if grid.available_area[ty, tx] <= 0.0:
                        continue
                    if active[ty, tx] or remaining[ty, tx] < area[node]:
                        continue
                    centre = np.array([
                        .5 * (grid.x_edges[tx] + grid.x_edges[tx + 1]),
                        .5 * (grid.y_edges[ty] + grid.y_edges[ty + 1]),
                    ])
                    target = _net_target(db, centres, int(node))
                    hpwl_estimate = float(np.abs(centre - target).sum() - np.abs(centres[node] - target).sum())
                    risk = area[node] / max(remaining[ty, tx], 1.0)
                    score = float(radius) + 2.0 * max(hpwl_estimate, 0.0) + 4.0 * risk
                    targets.append((score, tx, ty, radius))
                    break
            if not targets:
                continue
            _, tx, ty, radius = min(targets, key=lambda item: item[0])
            old = candidate[node].copy()
            candidate[node] = (
                .5 * (grid.x_edges[tx] + grid.x_edges[tx + 1]),
                .5 * (grid.y_edges[ty] + grid.y_edges[ty + 1]),
            )
            remaining[ty, tx] -= area[node]
            moves.append({
                "node": int(node), "source_bin": int(source), "target_bin": int(ty * grid.nx + tx),
                "radius_bins": int(radius), "dx": float(candidate[node, 0] - old[0]),
                "dy": float(candidate[node, 1] - old[1]),
            })
    return (project_centres(db, candidate) if moves else None), moves
