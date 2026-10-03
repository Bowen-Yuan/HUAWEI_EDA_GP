"""PyTorch implementation of the paper's random-batch splitting method (RBSM).

The published pseudocode has conflicting f1/f2 names and an ascent sign.  This module
uses equations (8), (14), and (18): HPWL descent followed by local penalty descent.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import csv
import json
import math
import time
from typing import Union

import numpy as np
import torch

from .bookshelf import Design, read_placement, write_placement
from .metrics import evaluate


@dataclass
class SolverConfig:
    mode: str = "resolved"
    epochs: int = 85
    inner_steps: int = 25
    lr0: float = 0.1
    gamma0: float = 1000.0
    mean_field_alpha: float = 5.0
    batch_fraction: float = 0.20
    temperature: float = 50000.0
    initialization: str = "bookshelf"
    warm_start_path: str = ""
    pair_refresh_steps: int = 5
    max_pairs: int = 1_500_000
    max_displacement: float = 12.0
    noise_scale: float = 0.2
    noise_power: float = 3.0
    importance_weighted: bool = True
    device: str = "auto"
    density_bins: int = 512
    seed: int = 2026


def config_from_json(path: Union[str, Path]) -> SolverConfig:
    values = json.loads(Path(path).read_text(encoding="utf-8"))
    return SolverConfig(**values)


def _net_ids(design: Design) -> np.ndarray:
    return np.repeat(np.arange(design.net_weight.size, dtype=np.int64), np.diff(design.net_start))


def _approx_degree(design: Design) -> np.ndarray:
    """Projected-net degree approximation, stable enough for sampling calibration."""
    net_id = _net_ids(design)
    node_incidence = np.bincount(design.pin_node, minlength=len(design.names))
    return np.bincount(net_id, weights=np.maximum(node_incidence[design.pin_node] - 1, 0),
                       minlength=design.net_weight.size).astype(np.float64)


def _connectivity_seed(design: Design, seed: int) -> np.ndarray:
    """A deterministic warm start used only by the resolved mode/tuning experiments."""
    rng = np.random.default_rng(seed)
    x0, y0, x1, y1 = design.core_bounds
    p = np.empty_like(design.centres)
    p[:, 0], p[:, 1] = (x0 + x1) / 2, (y0 + y1) / 2
    movable = design.movable
    p[design.fixed] = design.centres[design.fixed]
    net_id = _net_ids(design); pin_node = design.pin_node
    counts = np.bincount(net_id, minlength=design.net_weight.size).astype(np.float64)
    for _ in range(5):
        pin = p[pin_node] + design.pin_offset
        sums = np.zeros((design.net_weight.size, 2), np.float64)
        np.add.at(sums, net_id, pin); net_centre = sums / counts[:, None]
        node_sum = np.zeros_like(p); node_count = np.zeros(len(design.names), np.float64)
        np.add.at(node_sum, pin_node, net_centre[net_id] - design.pin_offset)
        np.add.at(node_count, pin_node, 1)
        target = node_sum / np.maximum(node_count[:, None], 1.0)
        p[movable] = 0.35 * p[movable] + 0.65 * target[movable]
        p[movable, 0] = np.clip(p[movable, 0], x0 + design.width[movable] / 2, x1 - design.width[movable] / 2)
        p[movable, 1] = np.clip(p[movable, 1], y0 + design.height[movable] / 2, y1 - design.height[movable] / 2)
    # Capacity spreading preserves each propagated target's coarse region while
    # preventing unanchored components from collapsing at the core centre.
    macro = 32; bw, bh = (x1 - x0) / macro, (y1 - y0) / macro
    movable_ids = np.flatnonzero(movable); area = design.width * design.height
    capacity = np.full(macro * macro, 0.8 * bw * bh, np.float64)
    target_bin_x = np.clip(((p[movable_ids, 0] - x0) / bw).astype(np.int64), 0, macro - 1)
    target_bin_y = np.clip(((p[movable_ids, 1] - y0) / bh).astype(np.int64), 0, macro - 1)
    target_bins = target_bin_y * macro + target_bin_x
    order = movable_ids[np.argsort(area[movable_ids])[::-1]]
    target_of = np.empty(len(design.names), np.int64); target_of[movable_ids] = target_bins
    for node in order:
        wanted = int(target_of[node])
        if capacity[wanted] < area[node]:
            wx, wy = wanted % macro, wanted // macro
            candidates = np.flatnonzero(capacity >= area[node])
            if candidates.size:
                cx, cy = candidates % macro, candidates // macro
                wanted = int(candidates[np.argmin((cx - wx) ** 2 + (cy - wy) ** 2)])
            else:
                wanted = int(np.argmax(capacity))
        capacity[wanted] -= area[node]
        ix, iy = wanted % macro, wanted // macro
        p[node, 0] = x0 + (ix + rng.random()) * bw
        p[node, 1] = y0 + (iy + rng.random()) * bh
    p[movable, 0] = np.clip(p[movable, 0], x0 + design.width[movable] / 2, x1 - design.width[movable] / 2)
    p[movable, 1] = np.clip(p[movable, 1], y0 + design.height[movable] / 2, y1 - design.height[movable] / 2)
    return p


def initial_positions(design: Design, config: SolverConfig) -> np.ndarray:
    if config.initialization == "warm_start":
        if not config.warm_start_path: raise ValueError("warm_start_path is required for warm_start")
        return read_placement(design, config.warm_start_path)
    if config.initialization == "bookshelf": return design.centres.copy()
    if config.initialization == "connectivity": return _connectivity_seed(design, config.seed)
    if config.initialization != "random": raise ValueError(f"Unknown initialization {config.initialization}")
    rng = np.random.default_rng(config.seed); x0, y0, x1, y1 = design.core_bounds
    p = design.centres.copy(); movable = design.movable
    p[movable, 0] = rng.uniform(x0 + design.width[movable] / 2, x1 - design.width[movable] / 2)
    p[movable, 1] = rng.uniform(y0 + design.height[movable] / 2, y1 - design.height[movable] / 2)
    return p


def candidate_pairs(design: Design, centres: np.ndarray, max_pairs: int, rng_seed: int = 0) -> np.ndarray:
    """Find exact-overlap candidates with a grid; fixed macros span every covered bin."""
    movable = design.movable
    sx = max(96.0, float(design.width[movable].max()) * 1.25)
    sy = max(24.0, float(design.height[movable].max()) * 2.0)
    x0, y0, x1, y1 = design.core_bounds
    nx = max(1, int(math.ceil((x1 - x0) / sx)))
    ny = max(1, int(math.ceil((y1 - y0) / sy)))
    buckets: dict[tuple[int, int], list[int]] = {}
    for i in np.flatnonzero(movable):
        key = (int(math.floor((centres[i, 0] - x0) / sx)), int(math.floor((centres[i, 1] - y0) / sy)))
        buckets.setdefault(key, []).append(int(i))
    # Fixed macros are put in all cells they cover; fixed-fixed candidates are ignored.
    for i in np.flatnonzero(design.fixed):
        l, r = centres[i, 0] - design.width[i] / 2, centres[i, 0] + design.width[i] / 2
        b, t = centres[i, 1] - design.height[i] / 2, centres[i, 1] + design.height[i] / 2
        ix0, ix1 = max(0, int(math.floor((l - x0) / sx))), min(nx - 1, int(math.floor((r - x0) / sx)))
        iy0, iy1 = max(0, int(math.floor((b - y0) / sy))), min(ny - 1, int(math.floor((t - y0) / sy)))
        for iy in range(iy0, iy1 + 1):
            for ix in range(ix0, ix1 + 1):
                buckets.setdefault((ix, iy), []).append(int(i))
    pairs: list[tuple[int, int]] = []; rng = np.random.default_rng(rng_seed)
    seen: set[tuple[int, int]] = set()
    for key, items in buckets.items():
        remaining = max_pairs - len(pairs)
        if len(items) * (len(items) - 1) // 2 > remaining:
            # A collapsed seed can put every standard cell in one bin.  Taking
            # the lexicographic prefix would update only the first cell IDs;
            # uniformly sampled local pairs preserve stochastic coverage.
            item = np.asarray(items, np.int64); added = 0
            while added < remaining:
                a = int(item[rng.integers(item.size)]); b = int(item[rng.integers(item.size)])
                if a == b or (design.fixed[a] and design.fixed[b]): continue
                pairs.append((a, b) if a < b else (b, a)); added += 1
            return np.asarray(pairs, np.int64)
        for a in range(len(items)):
            i = items[a]
            for j in items[a + 1:]:
                if design.fixed[i] and design.fixed[j]: continue
                pair = (i, j) if i < j else (j, i)
                if pair not in seen:
                    seen.add(pair); pairs.append(pair)
                    if len(pairs) >= max_pairs: return np.asarray(pairs, np.int64)
    return np.asarray(pairs, np.int64).reshape((-1, 2))


class RBSMSolver:
    def __init__(self, design: Design, config: SolverConfig):
        self.design, self.config = design, config
        if config.device == "auto": config.device = "cuda" if torch.cuda.is_available() else "cpu"
        self.device = torch.device(config.device)
        if self.device.type == "cuda" and not torch.cuda.is_available(): raise RuntimeError("CUDA requested but unavailable")
        torch.manual_seed(config.seed); np.random.seed(config.seed)
        self.dtype = torch.float32
        self.width = torch.as_tensor(design.width, device=self.device, dtype=self.dtype)
        self.height = torch.as_tensor(design.height, device=self.device, dtype=self.dtype)
        self.fixed = torch.as_tensor(design.fixed, device=self.device)
        self.movable = ~self.fixed
        self.pin_node = torch.as_tensor(design.pin_node, device=self.device)
        self.pin_offset = torch.as_tensor(design.pin_offset, device=self.device, dtype=self.dtype)
        self.net_start = torch.as_tensor(design.net_start, device=self.device)
        self.net_weight = torch.as_tensor(design.net_weight, device=self.device, dtype=self.dtype)
        self.degree = _approx_degree(design)
        logits = torch.as_tensor(self.degree / max(config.temperature, 1.0), device=self.device, dtype=self.dtype)
        self.prob = torch.softmax(logits - logits.max(), 0)
        self.positions = torch.as_tensor(initial_positions(design, config), device=self.device, dtype=self.dtype)
        self.fixed_positions = self.positions[self.fixed].clone()
        self.pair_gamma: dict[int, float] = {}
        self.global_step = 0
        self.pair_refresh = 0

    def _batch_hpwl(self, pos: torch.Tensor) -> torch.Tensor:
        m = self.net_weight.numel(); b = max(1, min(m, int(math.ceil(m * self.config.batch_fraction))))
        nets = torch.multinomial(self.prob, b, replacement=True)
        starts, ends = self.net_start[nets], self.net_start[nets + 1]
        lengths = ends - starts; total = int(lengths.sum().item())
        offsets = torch.cumsum(lengths, 0) - lengths
        local = torch.repeat_interleave(torch.arange(b, device=self.device), lengths)
        idx = torch.repeat_interleave(starts - offsets, lengths) + torch.arange(total, device=self.device)
        pins = pos[self.pin_node[idx]] + self.pin_offset[idx]
        inf = torch.full((b,), float("inf"), device=self.device, dtype=self.dtype)
        ninf = torch.full((b,), float("-inf"), device=self.device, dtype=self.dtype)
        mnx = inf.scatter_reduce(0, local, pins[:, 0], reduce="amin", include_self=True)
        mxx = ninf.scatter_reduce(0, local, pins[:, 0], reduce="amax", include_self=True)
        mny = inf.scatter_reduce(0, local, pins[:, 1], reduce="amin", include_self=True)
        mxy = ninf.scatter_reduce(0, local, pins[:, 1], reduce="amax", include_self=True)
        each = self.net_weight[nets] * (mxx - mnx + mxy - mny)
        return (each / self.prob[nets]).mean() if self.config.importance_weighted else each.mean() * m

    def _full_hpwl_grad(self) -> torch.Tensor:
        pos = self.positions.detach().clone().requires_grad_(True)
        pin = pos[self.pin_node] + self.pin_offset
        net_id = torch.repeat_interleave(torch.arange(self.net_weight.numel(), device=self.device), self.net_start[1:] - self.net_start[:-1])
        m = self.net_weight.numel(); inf = torch.full((m,), float("inf"), device=self.device, dtype=self.dtype); ninf = -inf
        max_x = ninf.scatter_reduce(0, net_id, pin[:, 0], reduce="amax", include_self=True)
        min_x = inf.scatter_reduce(0, net_id, pin[:, 0], reduce="amin", include_self=True)
        max_y = ninf.scatter_reduce(0, net_id, pin[:, 1], reduce="amax", include_self=True)
        min_y = inf.scatter_reduce(0, net_id, pin[:, 1], reduce="amin", include_self=True)
        loss = (self.net_weight * ((max_x - min_x) + (max_y - min_y))).sum()
        loss.backward(); return pos.grad.detach()

    def _step(self, loss_fn, lr: float) -> float:
        pos = self.positions.detach().clone().requires_grad_(True)
        loss = loss_fn(pos)
        if not torch.isfinite(loss): raise FloatingPointError("Non-finite RBSM loss")
        loss.backward(); grad = pos.grad
        self.global_step += 1
        eps = self.config.noise_scale / (self.global_step ** self.config.noise_power)
        noise = torch.randn_like(grad)
        direction = grad + eps * torch.linalg.vector_norm(grad) * noise
        with torch.no_grad():
            delta = lr * direction
            # The paper's lr=0.1 and gamma0=1000 can move a standard cell by
            # hundreds of sites in one update.  The resolved mode caps the physical
            # move; paper_equations keeps the raw step for an explicit ablation.
            if self.config.mode == "resolved" and self.config.max_displacement > 0:
                norm = torch.linalg.vector_norm(delta, dim=1, keepdim=True).clamp_min(1.0e-12)
                delta = delta * torch.clamp(self.config.max_displacement / norm, max=1.0)
            next_pos = pos - delta
            next_pos[self.fixed] = self.fixed_positions
        self.positions = next_pos.detach()
        return float(loss.detach().cpu())

    def _pair_weights(self, pairs: torch.Tensor, hpwl_grad: torch.Tensor) -> torch.Tensor:
        """Equation (16) with a sparse, nondecreasing active-pair cache."""
        if pairs.numel() == 0: return torch.empty(0, device=self.device, dtype=self.dtype)
        i, j = pairs[:, 0], pairs[:, 1]
        sensitivity = torch.maximum(hpwl_grad[i].abs().amax(1), hpwl_grad[j].abs().amax(1))
        rx, ry = (self.width[i] + self.width[j]) / 2, (self.height[i] + self.height[j]) / 2
        # The active direction is not known until the loss is evaluated.  Using the
        # smaller slope is conservative and avoids a zero divisor at a hat kink.
        slope = torch.minimum(1 / rx, 1 / ry)
        current = torch.maximum(torch.full_like(sensitivity, self.config.gamma0), torch.ceil(sensitivity / slope))
        keys = (pairs[:, 0].long() * len(self.design.names) + pairs[:, 1].long()).detach().cpu().numpy()
        values = current.detach().cpu().numpy()
        merged: list[float] = []
        for k, value in zip(keys, values):
            old = self.pair_gamma.get(int(k), self.config.gamma0)
            value = max(old, float(value)); self.pair_gamma[int(k)] = value
            merged.append(value)
        # The cache deliberately retains only pairs that have actually become local.
        merged = np.asarray(merged, dtype=np.float32)
        return torch.as_tensor(merged, device=self.device, dtype=self.dtype)

    def _penalty_loss(self, pairs: torch.Tensor, weights: torch.Tensor):
        gamma0 = self.config.gamma0
        def loss(pos: torch.Tensor) -> torch.Tensor:
            boundary = torch.relu(self.width / 2 - pos[:, 0]) + torch.relu(pos[:, 0] - (self._x1 - self.width / 2)) + torch.relu(self.height / 2 - pos[:, 1]) + torch.relu(pos[:, 1] - (self._y1 - self.height / 2))
            boundary = (boundary[self.movable] * gamma0).sum()
            if pairs.numel() == 0: return boundary
            i, j = pairs[:, 0], pairs[:, 1]
            dx, dy = pos[i, 0] - pos[j, 0], pos[i, 1] - pos[j, 1]
            rx, ry = (self.width[i] + self.width[j]) / 2, (self.height[i] + self.height[j]) / 2
            hx, hy = torch.relu(1 - dx.abs() / rx), torch.relu(1 - dy.abs() / ry)
            hat = torch.minimum(hx, hy)
            return boundary + (weights * hat).sum()
        return loss

    def run(self, output_dir: Union[str, Path]) -> dict[str, float]:
        output = Path(output_dir); output.mkdir(parents=True, exist_ok=True)
        self._x0, self._y0, self._x1, self._y1 = self.design.core_bounds
        log_path = output / "epochs.csv"
        with log_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=["epoch", "lr", "hpwl_loss", "penalty_loss", "pairs", "elapsed_seconds"])
            writer.writeheader(); start = time.perf_counter()
            for epoch in range(self.config.epochs):
                lr = self.config.lr0 * (1 + math.cos(math.pi * epoch / max(self.config.epochs, 1))) / 2
                hpwl_grad = self._full_hpwl_grad()
                hpwl_loss = penalty_loss = 0.0; pairs = torch.empty((0, 2), device=self.device, dtype=torch.long); pair_weights = torch.empty(0, device=self.device)
                for step in range(self.config.inner_steps):
                    hpwl_loss += self._step(lambda p: self._batch_hpwl(p) + self.config.mean_field_alpha * ((p[self.movable] - p[self.movable].mean(0)) ** 2).sum(), lr)
                    if step % self.config.pair_refresh_steps == 0:
                        cpu_pos = self.positions.detach().cpu().numpy().astype(np.float64)
                        pairs = torch.as_tensor(candidate_pairs(self.design, cpu_pos, self.config.max_pairs, self.config.seed + self.pair_refresh), device=self.device, dtype=torch.long)
                        self.pair_refresh += 1
                        pair_weights = self._pair_weights(pairs, hpwl_grad)
                    penalty_loss += self._step(self._penalty_loss(pairs, pair_weights), lr)
                writer.writerow({"epoch": epoch + 1, "lr": lr, "hpwl_loss": hpwl_loss / self.config.inner_steps, "penalty_loss": penalty_loss / self.config.inner_steps, "pairs": int(pairs.shape[0]), "elapsed_seconds": time.perf_counter() - start}); handle.flush()
        centres = self.positions.detach().cpu().numpy().astype(np.float64)
        write_placement(self.design, centres, output / "solution.pl")
        metrics = evaluate(self.design, centres, self.config.density_bins)
        metrics.update({"epochs": self.config.epochs, "seed": self.config.seed, "device": str(self.device), "elapsed_seconds": time.perf_counter() - start})
        (output / "summary.json").write_text(json.dumps({"config": asdict(self.config), "metrics": metrics}, indent=2), encoding="utf-8")
        return metrics
