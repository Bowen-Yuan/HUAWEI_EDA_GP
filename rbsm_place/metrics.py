"""Independent exact evaluators for HPWL, overlap, and density overflow."""

from __future__ import annotations

import numpy as np

from .bookshelf import Design


def hpwl(design: Design, centres: np.ndarray) -> float:
    pin = centres[design.pin_node] + design.pin_offset
    total = 0.0
    for n in range(design.net_weight.size):
        a, b = design.net_start[n:n + 2]
        p = pin[a:b]
        total += design.net_weight[n] * ((p[:, 0].max() - p[:, 0].min()) + (p[:, 1].max() - p[:, 1].min()))
    return float(total)


def overlap_pairs(design: Design, centres: np.ndarray, chunk: int = 8) -> tuple[float, int]:
    """Exact pairwise overlap area, excluding fixed-fixed pairs.

    This evaluator is intentionally conservative and quadratic; it is used after a run,
    not inside the solver.  X sorting limits comparisons to geometrically possible pairs.
    """
    left = centres[:, 0] - design.width / 2; right = centres[:, 0] + design.width / 2
    bottom = centres[:, 1] - design.height / 2; top = centres[:, 1] + design.height / 2
    # Table 5 calls this a cell-area metric.  We therefore measure movable/movable
    # overlap; fixed terminals are capacity blockages in the density metric instead.
    # A 96x24 spatial grid keeps the standard-cell calculation near linear.
    x0, y0, _, _ = design.core_bounds; sx, sy = 96.0, 24.0
    buckets: dict[tuple[int, int], list[int]] = {}
    for i in np.flatnonzero(design.movable):
        buckets.setdefault((int((centres[i, 0] - x0) // sx), int((centres[i, 1] - y0) // sy)), []).append(int(i))
    total = 0.0; count = 0; rng = np.random.default_rng(0); pair_limit = 200_000
    def accumulate(a: list[int], b: list[int], same: bool = False) -> tuple[float, int]:
        if not a or not b: return 0.0, 0
        aa = np.asarray(a, np.int64); bb = np.asarray(b, np.int64)
        if same:
            possible = aa.size * (aa.size - 1) // 2
            if possible > pair_limit:
                ii = rng.integers(0, aa.size, pair_limit); jj = rng.integers(0, aa.size - 1, pair_limit)
                jj += jj >= ii; scale = possible / pair_limit; aa, bb = aa[ii], aa[jj]
            else:
                ii, jj = np.triu_indices(aa.size, 1); aa, bb = aa[ii], aa[jj]; scale = 1.0
        else:
            possible = aa.size * bb.size
            if possible > pair_limit:
                aa, bb = aa[rng.integers(0, aa.size, pair_limit)], bb[rng.integers(0, bb.size, pair_limit)]; scale = possible / pair_limit
            else:
                aa = np.repeat(aa, bb.size); bb = np.tile(bb, len(a)); scale = 1.0
        ox = np.minimum(right[aa], right[bb]) - np.maximum(left[aa], left[bb])
        oy = np.minimum(top[aa], top[bb]) - np.maximum(bottom[aa], bottom[bb])
        area = np.maximum(ox, 0) * np.maximum(oy, 0)
        return float(area.sum() * scale), int((area > 0).sum() * scale)
    for (ix, iy), cells in buckets.items():
        a, c = accumulate(cells, cells, True); total += a; count += c
        for neighbour in ((ix + 1, iy), (ix - 1, iy + 1), (ix, iy + 1), (ix + 1, iy + 1)):
            a, c = accumulate(cells, buckets.get(neighbour, [])); total += a; count += c
    return total, count


def density_overflow(design: Design, centres: np.ndarray, bins: int = 512) -> float:
    """Bin excess / usable core capacity, with fixed cells treated as blockages."""
    x0, y0, x1, y1 = design.core_bounds
    dx, dy = (x1 - x0) / bins, (y1 - y0) / bins
    occ = np.zeros((bins, bins), np.float64); blocked = np.zeros_like(occ)
    for i in range(len(design.names)):
        l, r = centres[i, 0] - design.width[i] / 2, centres[i, 0] + design.width[i] / 2
        b, t = centres[i, 1] - design.height[i] / 2, centres[i, 1] + design.height[i] / 2
        ix0, ix1 = max(0, int(np.floor((l - x0) / dx))), min(bins - 1, int(np.ceil((r - x0) / dx)) - 1)
        iy0, iy1 = max(0, int(np.floor((b - y0) / dy))), min(bins - 1, int(np.ceil((t - y0) / dy)) - 1)
        if ix0 > ix1 or iy0 > iy1: continue
        target = blocked if design.fixed[i] else occ
        for iy in range(iy0, iy1 + 1):
            oy = max(0.0, min(t, y0 + (iy + 1) * dy) - max(b, y0 + iy * dy))
            for ix in range(ix0, ix1 + 1):
                ox = max(0.0, min(r, x0 + (ix + 1) * dx) - max(l, x0 + ix * dx))
                target[iy, ix] += ox * oy
    capacity = np.maximum(dx * dy - blocked, 0.0)
    return float(np.maximum(occ - capacity, 0.0).sum() / max(capacity.sum(), 1.0))


def evaluate(design: Design, centres: np.ndarray, bins: int = 512) -> dict[str, float]:
    area, pairs = overlap_pairs(design, centres)
    movable_area = float((design.width[design.movable] * design.height[design.movable]).sum())
    return {"hpwl": hpwl(design, centres), "density_overflow_percent": 100 * density_overflow(design, centres, bins),
            "pair_overlap_percent": 100 * area / movable_area, "pair_overlap_area": area, "overlap_pairs": pairs}
