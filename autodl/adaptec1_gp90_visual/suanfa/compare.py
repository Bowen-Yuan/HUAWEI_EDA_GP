# -*- coding: utf-8 -*-
import os
import csv
import math
import random
import time
from pathlib import Path
import numpy as np
from typing import Optional
import torch.profiler as tprof
import torch
from torch import nn, optim
import time, functools, torch
from torch import amp
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
PROJECT_DIR = Path(__file__).resolve().parent.parent
ISPD2005_CASE_DIR = PROJECT_DIR / "ispd2005" / "adaptec1"
ENABLE_MEASURE = False  # 需要计时就改 True

class Measure:
    def __init__(self, name, cuda=False):
        self.name = name
        self.cuda = cuda
    def __enter__(self):
        if not ENABLE_MEASURE:
            # 返回 self，但啥也不做
            return self
        if self.cuda and torch.cuda.is_available():
            torch.cuda.synchronize()
            torch.cuda.reset_peak_memory_stats()
            self.mem0 = torch.cuda.memory_allocated()
        self.t0 = time.perf_counter()
        return self
    def __exit__(self, exc_type, exc, tb):
        if not ENABLE_MEASURE:
            return False
        t1 = time.perf_counter()
        if self.cuda and torch.cuda.is_available():
            torch.cuda.synchronize()
            mem1 = torch.cuda.memory_allocated()
            peak = torch.cuda.max_memory_allocated()
            print(f"[{self.name}] {1000*(t1-self.t0):.2f} ms | cuda_alloc Δ={ (mem1 - self.mem0)/1e6:.1f} MB | cuda_peak={peak/1e6:.1f} MB")
        else:
            print(f"[{self.name}] {1000*(t1-self.t0):.2f} ms")

def measure_fn(name=None, cuda=False):
    def deco(fn):
        tag = name or fn.__name__
        @functools.wraps(fn)
        def wrapper(*args, **kwargs):
            with Measure(tag, cuda=cuda):
                return fn(*args, **kwargs)
        return wrapper
    return deco
# =========================
# Data readers
# =========================
def read_scl_bounds(scl_file: str):
    """
    读取 ISPD2005 .scl 文件，返回布局区域宽度 W、高度 H 和边界 (x_min, x_max, y_min, y_max)
    """
    x_min = float("inf")
    x_max = float("-inf")
    y_min = float("inf")
    y_max = float("-inf")
    sitewidth = 1.0  # 默认值，会在文件里覆盖

    with open(scl_file, "r") as f:
        cur = {}
        for raw in f:
            line = raw.strip()
            if not line or line.startswith(("UCLA", "#", "NumRows")):
                continue

            if line.startswith("Sitewidth"):
                sitewidth = float(line.split(":")[-1])
            elif line.startswith("Coordinate"):
                cur["y"] = float(line.split(":")[-1])
            elif line.startswith("Height"):
                cur["h"] = float(line.split(":")[-1])
            elif line.startswith("SubrowOrigin"):
                parts = line.replace("\t", " ").split()
                # e.g. SubrowOrigin  :  459    NumSites  :  10692
                x0 = float(parts[2])
                nsites = int(parts[-1])
                cur["x0"] = x0
                cur["x1"] = x0 + nsites * sitewidth
            elif line.startswith("End"):
                if all(k in cur for k in ("y", "h", "x0", "x1")):
                    y0 = cur["y"]
                    y1 = cur["y"] + cur["h"]
                    x_min = min(x_min, cur["x0"])
                    x_max = max(x_max, cur["x1"])
                    y_min = min(y_min, y0)
                    y_max = max(y_max, y1)
                cur = {}

    W = x_max - x_min
    H = y_max - y_min
    return W, H, (x_min, x_max, y_min, y_max)




class Nodes:
    """Read ISPD2005 'nodes' (aka 'nodes/nodes.txt')"""
    def __init__(self, node_file: str):
        self.file_name = node_file

    def read_nodes(self):
        blocks = []       # [w, h] (terminal -> [0,0] for overlap-free)
        block_dict = {}   # name -> global index (movable first, then terminals)
        sb_count = 0      # movable blocks count
        p_count = 0       # terminals count
        area = 0.0

        with open(self.file_name, 'r') as f:
            for raw in f:
                line = raw.strip()
                if not line or line.startswith(('UCLA', '#', 'NumNodes', 'NumTerminals')):
                    continue
                parts = line.split()
                name = parts[0]
                # ISPD2005 nodes: movable lines are "name w h"
                # terminals lines are typically "name w h terminal"
                w = float(parts[1]); h = float(parts[2])
                is_terminal = (len(parts) >= 4 and parts[-1].lower().startswith('terminal'))
                if is_terminal:
                    block_dict[name] = sb_count + p_count
                    p_count += 1
                    blocks.append([0.0, 0.0])       # terminals have no size in our overlap model
                else:
                    block_dict[name] = sb_count
                    sb_count += 1
                    blocks.append([w, h])
                    area += w * h

        return np.asarray(blocks, dtype=np.float32), sb_count, block_dict, area


class Nets:
    """Read ISPD2005 'nets'."""
    def __init__(self, net_file: str, block_dict: dict):
        self.file_name = net_file
        self.block_dict = block_dict

    def read_nets(self):
        nets_list = []
        max_deg = 0
        num_nets = 0
        cur = None

        with open(self.file_name, 'r') as f:
            for raw in f:
                line = raw.strip()
                if not line or line.startswith('#') or line.startswith('UCLA'):
                    continue
                if line.startswith('NumNets'):
                    num_nets = int(line.split(':')[-1])
                    continue
                if line.startswith('NumPins'):
                    continue
                if line.startswith('NetDegree'):
                    if cur:
                        nets_list.append(cur)
                    # e.g. "NetDegree : 4   n0"
                    deg = int(line.split(':')[1].split()[0])
                    cur = {'degree': deg, 'cells': []}
                    max_deg = max(max_deg, deg)
                    continue

                # pin line: "name I|O : offsetX offsetY"
                name = line.split()[0]
                if cur is not None:
                    cur['cells'].append(name)

            if cur:
                nets_list.append(cur)

        nets = -np.ones((len(nets_list), max_deg), dtype=int)
        for i, net in enumerate(nets_list):
            for j, name in enumerate(net['cells']):
                if name in self.block_dict:
                    nets[i, j] = self.block_dict[name]

        return nets, num_nets, nets_list


class P_Position(nn.Module):
    """Read ISPD2005 'pl' (fixed terminals positions)."""
    def __init__(self, sb_count: int, block_dict: dict, plfile: str):
        super().__init__()
        num_terms = len(block_dict) - sb_count
        self.fixed_positions = torch.zeros((num_terms, 2), dtype=torch.float32, requires_grad=False)

        # terminal's index in all_positions is (sb_count + term_idx)
        # where term_idx = block_dict[name] - sb_count
        with open(plfile, 'r') as f:
            for raw in f:
                line = raw.strip()
                if not line or line.startswith(('UCLA', '#')):
                    continue
                parts = line.split()
                name = parts[0]
                if name not in block_dict:
                    continue
                gidx = block_dict[name]
                t_ofs = gidx - sb_count
                if 0 <= t_ofs < num_terms:
                    x = float(parts[1]); y = float(parts[2])
                    self.fixed_positions[t_ofs] = torch.tensor([x, y], dtype=torch.float32)


# =========================
# Placement variables
# =========================
class Position(nn.Module):
    def __init__(self, sb_count: int, W: float, H: float):
        super().__init__()
        self.positions = nn.Parameter(torch.empty((sb_count, 2), dtype=torch.float32))
        with torch.no_grad():
            nn.init.uniform_(self.positions[:, 0], 0.15 * W, 0.85 * W)
            nn.init.uniform_(self.positions[:, 1], 0.15 * H, 0.85 * H)

    def forward(self):
        return self.positions


# =========================
# Losses (HPWL / Boundary / Overlap with spatial hashing)
# =========================
class loss_function:
    def __init__(self, gamma1, gamma2, W, H, bin_size: Optional[float] = None,
                 max_pairs: Optional[int] = None, use_hat: bool = True):
        self.gamma1 = gamma1
        self.gamma2 = gamma2
        self.W = W
        self.H = H
        self.bin_size = bin_size
        self.max_pairs = max_pairs
        self.use_hat = use_hat

    # ---- hat kernel ----
    @staticmethod
    def hat_function(dx, dy, r, t, eps=1e-12):
        abs_x = dx.abs()
        abs_y = dy.abs()
        r = torch.clamp(r, min=eps)
        t = torch.clamp(t, min=eps)
        outside = (abs_x > r) | (abs_y > t)
        x_dom = (abs_y / t) <= (abs_x / r)
        zero = dx.new_zeros(())
        one = dx.new_ones(())
        val_x = one - abs_x / r
        val_y = one - abs_y / t
        inside_val = torch.where(x_dom, val_x, val_y)
        return torch.where(outside, zero, inside_val)

    # ---- HPWL ----
    def calculate_distance_loss(self, nets, batch_idx, positions):
        device = positions.device

        # 1) 从 nets 上做索引，再把结果搬到 positions.device
        if isinstance(nets, torch.Tensor):
            # nets 可能在 CPU；用 nets.device 做索引
            if batch_idx.device != nets.device:
                batch_idx_local = batch_idx.to(nets.device)
            else:
                batch_idx_local = batch_idx
            batch_cells = nets.index_select(0, batch_idx_local).to(device=device)
        else:
            # 万一 nets 还是 numpy，就把 index 放到 CPU 再转 tensor
            batch_cells = torch.as_tensor(
                nets[batch_idx.cpu().numpy()], dtype=torch.long, device=device
            )

        valid_mask = (batch_cells != -1)
        safe_idx = torch.clamp(batch_cells, 0, positions.shape[0]-1)
        all_pos = positions[safe_idx]  
        inf  = torch.tensor(float('inf'),  device=device, dtype=positions.dtype)
        ninf = torch.tensor(float('-inf'), device=device, dtype=positions.dtype)
        min_pos = torch.where(valid_mask.unsqueeze(-1), all_pos, inf).amin(dim=1)
        max_pos = torch.where(valid_mask.unsqueeze(-1), all_pos, ninf).amax(dim=1)
        return (max_pos - min_pos).sum(dim=1).sum()




    # ---- Boundary ----
    def calculate_boundary_loss1(self, batch_idx, positions, blocks):
        device = positions.device
        p = positions[batch_idx]
        b = torch.as_tensor(blocks[batch_idx], dtype=positions.dtype, device=device)
        avg = positions.mean(0)
        r = torch.linalg.norm(avg - p, dim=1)
        w2, h2 = b[:,0]*0.5, b[:,1]*0.5
        left   = torch.relu(-(p[:,0]-w2))
        right  = torch.relu( p[:,0]-(self.W - w2))
        bottom = torch.relu(-(p[:,1]-h2))
        top    = torch.relu( p[:,1]-(self.H - h2))
        base = left + right + bottom + top
        if torch.is_tensor(self.gamma1):
            g1 = self.gamma1.to(device)[batch_idx]
        else:
            g1 = torch.as_tensor(self.gamma1, dtype=base.dtype, device=device)
        return (g1 * base + p.new_tensor(5.0) * r).sum()

    # ---- Spatial hashing helpers ----
    @staticmethod
    def _build_bins(pos, sizes, bin_size):
        K = pos.shape[0]
        hx = sizes[:,0]*0.5
        hy = sizes[:,1]*0.5
        bx0 = torch.floor((pos[:,0]-hx)/bin_size).to(torch.int32)
        bx1 = torch.floor((pos[:,0]+hx)/bin_size).to(torch.int32)
        by0 = torch.floor((pos[:,1]-hy)/bin_size).to(torch.int32)
        by1 = torch.floor((pos[:,1]+hy)/bin_size).to(torch.int32)
        bins = {}
        for i in range(K):
            for bx in range(int(bx0[i]), int(bx1[i])+1):
                for by in range(int(by0[i]), int(by1[i])+1):
                    bins.setdefault((bx,by), []).append(i)
        aabb_bins = torch.stack([bx0,bx1,by0,by1], dim=1)  # [K,4]
        return bins, aabb_bins

    def _gen_pairs(self, pos, sizes, scope, batch_idx_local=None, halo_bins=1, bin_size=None):
        if bin_size is None:
            med_w = torch.median(sizes[:,0])
            med_h = torch.median(sizes[:,1])
            bin_size = 1.5 * torch.max(med_w, med_h).clamp(min=1e-6)
        else:
            bin_size = pos.new_tensor(bin_size).clamp(min=1e-6)

        bins, aabb_bins = self._build_bins(pos, sizes, bin_size)
        pair_set = set()

        if scope in ('intra', 'global'):
            for _, idxs in bins.items():
                if len(idxs) <= 1:
                    continue
                idxs.sort()
                for a in range(len(idxs)-1):
                    ia = idxs[a]
                    for b in range(a+1, len(idxs)):
                        ib = idxs[b]
                        pair_set.add((ia, ib))

        elif scope == 'batch_halo':
            assert batch_idx_local is not None, "batch_halo requires batch_idx_local"
            anchor_mask = torch.zeros(pos.shape[0], dtype=torch.bool, device='cpu')
            anchor_mask[batch_idx_local.cpu()] = True
            aabb = aabb_bins.cpu()
            for ia in batch_idx_local.cpu().tolist():
                bx0, bx1, by0, by1 = aabb[ia].tolist()
                for bx in range(bx0 - halo_bins, bx1 + halo_bins + 1):
                    for by in range(by0 - halo_bins, by1 + halo_bins + 1):
                        idxs = bins.get((bx,by))
                        if not idxs:
                            continue
                        for ib in idxs:
                            if ib == ia:
                                continue
                            i, j = (ia, ib) if ia < ib else (ib, ia)
                            if anchor_mask[i] or anchor_mask[j]:
                                pair_set.add((i, j))
        else:
            raise ValueError(f"unknown scope: {scope}")

        if not pair_set:
            empty = torch.empty((0,), dtype=torch.long, device=pos.device)
            return empty, empty

        pairs = torch.tensor(list(pair_set), dtype=torch.long, device=pos.device)
        if (self.max_pairs is not None) and (pairs.shape[0] > self.max_pairs):
            perm = torch.randperm(pairs.shape[0], device=pos.device)
            pairs = pairs[perm[:self.max_pairs]]
        return pairs[:,0], pairs[:,1]

    # ---- Overlap (scope-based, no per-pair weights) ----
    def calculate_overlap_loss(self, batch_idx, positions, blocks, scope='intra', halo_bins=1):
        device = positions.device
        dtype  = positions.dtype
        blocks_t = torch.as_tensor(blocks, dtype=dtype, device=device)

        if scope == 'intra':
            pos = positions[batch_idx]
            siz = blocks_t[batch_idx]
            if pos.shape[0] <= 1:
                return pos.new_zeros(())
            i_idx, j_idx = self._gen_pairs(pos, siz, scope='intra', bin_size=self.bin_size)
            if i_idx.numel() == 0:
                return pos.new_zeros(())
            pi, pj = pos[i_idx], pos[j_idx]
            si, sj = siz[i_idx], siz[j_idx]

        elif scope == 'batch_halo':
            pos_all = positions
            siz_all = blocks_t
            if positions.shape[0] <= 1 or len(batch_idx) == 0:
                return positions.new_zeros(())
            i_idx, j_idx = self._gen_pairs(pos_all, siz_all, scope='batch_halo',
                                           batch_idx_local=torch.as_tensor(batch_idx, dtype=torch.long, device=device),
                                           halo_bins=halo_bins, bin_size=self.bin_size)
            if i_idx.numel() == 0:
                return positions.new_zeros(())
            pi, pj = pos_all[i_idx], pos_all[j_idx]
            si, sj = siz_all[i_idx], siz_all[j_idx]

        elif scope == 'global':
            pos_all = positions
            siz_all = blocks_t
            if positions.shape[0] <= 1:
                return positions.new_zeros(())
            i_idx, j_idx = self._gen_pairs(pos_all, siz_all, scope='global', bin_size=self.bin_size)
            if i_idx.numel() == 0:
                return positions.new_zeros(())
            pi, pj = pos_all[i_idx], pos_all[j_idx]
            si, sj = siz_all[i_idx], siz_all[j_idx]
        else:
            raise ValueError(f"unknown scope: {scope}")

        if self.use_hat:
            dx = pi[:,0] - pj[:,0]
            dy = pi[:,1] - pj[:,1]
            r  = 0.5 * (si[:,0] + sj[:,0])
            t  = 0.5 * (si[:,1] + sj[:,1])
            penalties = self.hat_function(dx, dy, r, t)
        else:
            hx_sum = 0.5 * (si[:,0] + sj[:,0])
            hy_sum = 0.5 * (si[:,1] + sj[:,1])
            ovx = torch.clamp(hx_sum - (pi[:,0]-pj[:,0]).abs(), min=0)
            ovy = torch.clamp(hy_sum - (pi[:,1]-pj[:,1]).abs(), min=0)
            penalties = ovx * ovy

        # gamma2: scalar or per-node -> average; otherwise use scalar fallback
        if torch.is_tensor(self.gamma2):
            g2 = self.gamma2.to(device)
            if g2.numel() == 1:
                weight = g2
            else:
                # for simplicity
                weight = penalties.new_tensor(1.0)
        else:
            weight = penalties.new_tensor(float(self.gamma2))

        return (weight * penalties).sum()

    # ---- Overlap (pair list + optional weights) ----
    def calculate_overlap_loss_pairs(self, positions, blocks, i_idx, j_idx,
                                     pair_weights: Optional[torch.Tensor] = None,
                                     use_hat: Optional[bool] = None):
        if i_idx.numel() == 0:
            return positions.new_zeros(())
        device = positions.device
        dtype  = positions.dtype
        use_hat = self.use_hat if use_hat is None else use_hat

        siz = torch.as_tensor(blocks, dtype=dtype, device=device)
        pi, pj = positions[i_idx], positions[j_idx]
        si, sj = siz[i_idx], siz[j_idx]

        if use_hat:
            dx = pi[:,0]-pj[:,0]; dy = pi[:,1]-pj[:,1]
            r  = 0.5*(si[:,0]+sj[:,0]); t = 0.5*(si[:,1]+sj[:,1])
            penalties = self.hat_function(dx, dy, r, t)
        else:
            hx = 0.5*(si[:,0]+sj[:,0]); hy = 0.5*(si[:,1]+sj[:,1])
            ovx = (hx - (pi[:,0]-pj[:,0]).abs()).clamp_min(0)
            ovy = (hy - (pi[:,1]-pj[:,1]).abs()).clamp_min(0)
            penalties = ovx * ovy

        if pair_weights is None:
            w = penalties.new_tensor(1.0)
        else:
            w = pair_weights.to(device=device, dtype=dtype)
        return (penalties * w).sum()


# =========================
# Spatial hash for gamma (module-level helper)
# =========================
def _spatial_hash_pairs(
    pos: torch.Tensor,
    sizes: torch.Tensor,
    bin_size: Optional[float] = None,
    halo_bins: int = 0,
    BIN_CAP: int = 128,              # 每个 bin 内最多保留的 cell 数
    PAIRS_PER_BIN_CAP: int = 2000,   # 每个 bin 最多抽样的 pair 数
    MAX_PAIRS_GLOBAL: int = 300_000  # 全局最多 pair 数
):
    """
    受控版 pair 生成：每个 bin 限流 + 抽样；避免 set 的 O(n^2) 内存爆炸
    返回: (i_idx, j_idx) 都是 1D LongTensor（位于 pos.device）
    """
    device = pos.device
    N = pos.shape[0]
    if N <= 1:
        e = torch.empty(0, dtype=torch.long, device=device); return e, e

    if bin_size is None:
        med = torch.max(torch.median(sizes[:,0]), torch.median(sizes[:,1])).clamp_min(1e-6)
        bin_size = float(med) * 1.5
    else:
        bin_size = float(bin_size)

    # 先在 CPU 上做散列/采样，最后一次性搬回 device，避免挤占 GPU 显存
    hx = (sizes[:,0]*0.5).cpu().numpy()
    hy = (sizes[:,1]*0.5).cpu().numpy()
    px = pos[:,0].detach().cpu().numpy()
    py = pos[:,1].detach().cpu().numpy()

    bx0 = np.floor((px - hx) / bin_size).astype(np.int32)
    bx1 = np.floor((px + hx) / bin_size).astype(np.int32)
    by0 = np.floor((py - hy) / bin_size).astype(np.int32)
    by1 = np.floor((py + hy) / bin_size).astype(np.int32)

    bins = {}
    for i in range(N):
        for bx in range(int(bx0[i]) - halo_bins, int(bx1[i]) + 1 + halo_bins):
            for by in range(int(by0[i]) - halo_bins, int(by1[i]) + 1 + halo_bins):
                bins.setdefault((bx, by), []).append(i)

    import random
    pairs_cpu = []   # 收集在 CPU
    total_cap = MAX_PAIRS_GLOBAL

    for _, idxs in bins.items():
        if len(idxs) <= 1:
            continue
        # —— 每 bin 限流（随机取 BIN_CAP 个）
        if len(idxs) > BIN_CAP:
            idxs = random.sample(idxs, BIN_CAP)
        m = len(idxs)
        max_pairs_bin = min(PAIRS_PER_BIN_CAP, m*(m-1)//2)
        if max_pairs_bin <= 0:
            continue

        # 用 numpy/torch 向量化在 CPU 上一次性采样若干对，然后在 bin 内去重
        # 过采样一倍提高去重后命中率
        k = min(2*max_pairs_bin, m*(m-1)//2)
        a = np.random.randint(0, m, size=(k,), dtype=np.int64)
        b = np.random.randint(0, m, size=(k,), dtype=np.int64)
        mask = (a != b)
        a, b = a[mask], b[mask]
        ia = np.array(idxs, dtype=np.int64)[a]
        ib = np.array(idxs, dtype=np.int64)[b]
        # 规范化 (i<j)
        swap = ia > ib
        ia_sw = ia.copy(); ib_sw = ib.copy()
        ia_sw[swap] = ib[swap]
        ib_sw[swap] = ia[swap]
        ia, ib = ia_sw, ib_sw
        if ia.size == 0:
            continue
        # bin 内去重
        key = ia.astype(np.int64) * int(N) + ib.astype(np.int64)
        _, uniq_idx = np.unique(key, return_index=True)
        ia = ia[uniq_idx][:max_pairs_bin]
        ib = ib[uniq_idx][:max_pairs_bin]

        # 追加到全局池；到上限就停
        for x, y in zip(ia.tolist(), ib.tolist()):
            pairs_cpu.append((x, y))
            if len(pairs_cpu) >= total_cap:
                break
        if len(pairs_cpu) >= total_cap:
            break

    if not pairs_cpu:
        e = torch.empty(0, dtype=torch.long, device=device); return e, e

    pairs = torch.tensor(pairs_cpu, dtype=torch.long, device=device)
    i_idx, j_idx = pairs[:,0], pairs[:,1]

    # 可选：全局再去一次重（这时数量已经受限，代价很小）
   # 可选：全局再去一次重（兼容 torch 1.x，无 return_index）
    key = i_idx * int(N) + j_idx
    sorted_key, order = torch.sort(key)
    keep = torch.ones_like(sorted_key, dtype=torch.bool)
    keep[1:] = sorted_key[1:] != sorted_key[:-1]
    keep_idx = order[keep]
    i_idx = i_idx[keep_idx]
    j_idx = j_idx[keep_idx]
    # 调试输出（可留可去）
    # print(f"[pairs] generated {i_idx.numel()} pairs "
    #     f"(cap={MAX_PAIRS_GLOBAL}, bin_cap={BIN_CAP}, per_bin={PAIRS_PER_BIN_CAP})")
    return i_idx, j_idx



# =========================
# Trainer
# =========================
class TrainModel:
    def __init__(self, model: Position, net_file: str, node_file: str, pl_file: str,
                 W: float, H: float, max_epoch: int, batch_size: int, max_step: int):
        self.model = model
        self.W, self.H = W, H
        self.max_epoch = max_epoch
        self.batch_size = batch_size
        self.max_setp = max_step
        self.loss_history = []
        # Nodes / Nets
        self.blocks, self.sb_count, self.block_dict, self.area = Nodes(node_file).read_nodes()
        self.block1 = self.blocks[:self.sb_count]   # movable sizes only
        self.nets, self.num_nets, self.nets_list = Nets(net_file, self.block_dict).read_nets()

        # Fixed terminals
        p_pos = P_Position(self.sb_count, self.block_dict, pl_file)
        self.fixed_positions = p_pos.fixed_positions  # [P,2]

        # Common indices
        self.full_idx_t = torch.arange(self.sb_count, dtype=torch.long)
        self.net_idx_t  = torch.arange(self.num_nets, dtype=torch.long)

        # Grid size (for density metric / batching if needed)
        med_w = np.median(self.block1[:,0]) if self.block1.size else 1.0
        med_h = np.median(self.block1[:,1]) if self.block1.size else 1.0
        self.grid_size = float(max(med_w, med_h)) * 2.0
        # __init__ 里加
        self.nets_t = torch.from_numpy(self.nets).long()  # 稠密矩阵版
        self.net_idx_t = torch.arange(self.num_nets, dtype=torch.long)  # 已有，同步一下
        # ---- HPWL 分段归约的预处理（一次性做）----
        nets_np = self.nets  # numpy
        num_nets, max_deg = nets_np.shape
        pin2net, pin2node = [], []
        for n in range(num_nets):
            row = nets_np[n]
            for j in range(max_deg):
                u = row[j]
                if u >= 0:
                    pin2net.append(n)
                    pin2node.append(u)
        # 先放 CPU，首次用到时再搬设备
        self.pin2net_cpu  = torch.tensor(pin2net,  dtype=torch.long)
        self.pin2node_cpu = torch.tensor(pin2node, dtype=torch.long)
        self.num_pins = len(pin2net)
        # 批内 net 映射缓冲（复用，避免反复分配）
        self.net_id_buffer = torch.full((self.num_nets,), -1, dtype=torch.long)

    def _ensure_pinmap_on_device(self, device):
        # 首次使用时把映射搬到同设备
        if not hasattr(self, "pin2net") or self.pin2net.device != device:
            self.pin2net  = self.pin2net_cpu.to(device, non_blocking=True)
            self.pin2node = self.pin2node_cpu.to(device, non_blocking=True)
            self.net_id_buffer = self.net_id_buffer.to(device, non_blocking=True)

    def hpwl_segmented(self, positions: torch.Tensor, fixed_pos: torch.Tensor, batch_nets: torch.Tensor):
        """分段归约版 HPWL（不构造 [B,max_deg,2] 大张量）"""
        device = positions.device
        dtype  = positions.dtype
        self._ensure_pinmap_on_device(device)

        Nmov = positions.shape[0]
        Nfix = fixed_pos.shape[0]

        # 节点坐标缓冲：movable 参与反传，fixed 切断梯度
        node_x = torch.empty(Nmov + Nfix, dtype=dtype, device=device)
        node_y = torch.empty_like(node_x)
        node_x[:Nmov] = positions[:, 0]
        node_y[:Nmov] = positions[:, 1]
        with torch.no_grad():
            node_x[Nmov:] = fixed_pos[:, 0]
            node_y[Nmov:] = fixed_pos[:, 1]

        # 把 batch 的 net id 压到 [0..B-1]，用于 scatter_reduce
        bid = self.net_id_buffer
        bid.fill_(-1)
        bid[batch_nets] = torch.arange(batch_nets.numel(), device=device, dtype=torch.long)

        pin_bid = bid[self.pin2net]            # [P] -> [-1..B-1]
        mask = pin_bid >= 0
        if not torch.any(mask):
            return positions.new_zeros(())

        pin_bid = pin_bid[mask]                # [P_b]
        pin_nid = self.pin2node[mask]          # [P_b]
        px = node_x[pin_nid]
        py = node_y[pin_nid]

        B = batch_nets.numel()
        INF  = torch.tensor(float('inf'),  device=device, dtype=dtype)
        NINF = torch.tensor(float('-inf'), device=device, dtype=dtype)
        xmin = torch.full((B,), INF,  dtype=dtype, device=device)
        xmax = torch.full((B,), NINF, dtype=dtype, device=device)
        ymin = xmin.clone()
        ymax = xmax.clone()

        # 需要 PyTorch>=2.0 的 scatter_reduce_
        xmin.scatter_reduce_(0, pin_bid, px, reduce='amin')
        xmax.scatter_reduce_(0, pin_bid, px, reduce='amax')
        ymin.scatter_reduce_(0, pin_bid, py, reduce='amin')
        ymax.scatter_reduce_(0, pin_bid, py, reduce='amax')

        return (xmax - xmin + ymax - ymin).sum()


    # γ (sparse) —— dis_grad / jacobian → max → ceil
    def calculate_gamma_sparse(self, positions: torch.Tensor):
        device = positions.device
        dtype  = positions.dtype
        N = self.sb_count

        fixed_pos = self.fixed_positions.to(device=device, dtype=dtype)
        sizes = torch.as_tensor(self.block1, dtype=dtype, device=device)
        w2, h2 = sizes[:,0]*0.5, sizes[:,1]*0.5

        # === 1) 采样部分 nets 估计 HPWL 梯度 ===
        pos_var = positions.clone().detach().requires_grad_(True)
        all_pos = torch.cat([pos_var, fixed_pos], dim=0)
        LossTmp = loss_function(1.0, 1.0, self.W, self.H)

        ALL = self.net_idx_t.to(all_pos.device)
        SAMPLE = min(1000, ALL.numel())
        bi = ALL[torch.randperm(ALL.numel(), device=ALL.device)[:SAMPLE]]

        distance_loss = LossTmp.calculate_distance_loss(self.nets_t, bi, all_pos)

        # 先反传，再打印（避免任何奇怪的交互）
        dis_grad = torch.autograd.grad(distance_loss, pos_var, retain_graph=False, create_graph=False)[0]
        if dis_grad.is_cuda: torch.cuda.synchronize()
        print(f"[gamma] sample_nets={SAMPLE}, loss={float(distance_loss):.3e}, ||dis_grad||={dis_grad.norm().item():.3e}")

        # 从这里起全部用“断图”的 pos_ng
        pos_ng = pos_var.detach()

        # === 2) gamma1（边界） ===
        x, y = pos_ng[:,0], pos_ng[:,1]
        left   = (-(x - w2) > 0); right = (x - (self.W - w2) > 0)
        bottom = (-(y - h2) > 0); top   = (y - (self.H - h2) > 0)
        g1x = torch.zeros(N, dtype=dtype, device=device)
        g1y = torch.zeros(N, dtype=dtype, device=device)
        g1x[left]  = dis_grad[left,0].abs()
        g1x[right] = torch.maximum(g1x[right], dis_grad[right,0].abs())
        g1y[bottom]= dis_grad[bottom,1].abs()
        g1y[top]   = torch.maximum(g1y[top],   dis_grad[top,1].abs())
        gamma1_new = torch.ceil(torch.maximum(g1x, g1y)).detach()

        # === 3) gamma2（只对可能重叠对儿） ===
        coarse_bin = float(torch.max(torch.median(sizes[:,0]), torch.median(sizes[:,1])))*4.0
        i_idx, j_idx = _spatial_hash_pairs(pos_ng, sizes, bin_size=coarse_bin, halo_bins=0)

        if i_idx.numel() == 0:
            eL = torch.empty(0, dtype=torch.long, device=device)
            eF = torch.empty(0, dtype=dtype,      device=device)
            return gamma1_new, (eL, eL, eF)

        pi, pj = pos_ng[i_idx], pos_ng[j_idx]
        si, sj = sizes[i_idx], sizes[j_idx]
        dx = pi[:,0]-pj[:,0]; dy = pi[:,1]-pj[:,1]
        r  = 0.5*(si[:,0]+sj[:,0]); t = 0.5*(si[:,1]+sj[:,1])
        absx, absy = dx.abs(), dy.abs()
        inside = (absx <= r) & (absy <= t)
        if not inside.any():
            eL = torch.empty(0, dtype=torch.long, device=device)
            eF = torch.empty(0, dtype=dtype,      device=device)
            return gamma1_new, (eL, eL, eF)

        m = inside
        xdom = (absy[m]/t[m].clamp_min(1e-12)) <= (absx[m]/r[m].clamp_min(1e-12))
        sgnx = dx[m].sign(); sgny = dy[m].sign()
        rr = r[m].clamp_min(1e-12); tt = t[m].clamp_min(1e-12)

        dpi_x = torch.zeros_like(rr); dpi_y = torch.zeros_like(tt)
        dpj_x = torch.zeros_like(rr); dpj_y = torch.zeros_like(tt)
        if xdom.any():
            dpi_x[xdom] = - sgnx[xdom] / rr[xdom]
            dpj_x[xdom] = + sgnx[xdom] / rr[xdom]
        ydm = ~xdom
        if ydm.any():
            dpi_y[ydm] = - sgny[ydm] / tt[ydm]
            dpj_y[ydm] = + sgny[ydm] / tt[ydm]

        dpi = torch.sqrt(dpi_x**2 + dpi_y**2).clamp_min(1e-12)
        dpj = torch.sqrt(dpj_x**2 + dpj_y**2).clamp_min(1e-12)
        gi = torch.linalg.norm(dis_grad[i_idx[m]], dim=1)
        gj = torch.linalg.norm(dis_grad[j_idx[m]], dim=1)
        g2 = torch.ceil(torch.maximum(gi/dpi, gj/dpj)).detach()

        return gamma1_new, (i_idx[m].detach(), j_idx[m].detach(), g2)



    def update_gradient(self, positions: torch.Tensor, loss: torch.Tensor, optimizer,
                        epoch: int, noise_scale_base: float = 0.0,
                        clip_norm: Optional[float] = None):
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        if clip_norm is not None:
            torch.nn.utils.clip_grad_norm_([positions], max_norm=clip_norm)
        with torch.no_grad():
            eps = noise_scale_base / ((epoch + 1) ** 2)
            g = positions.grad
            g_norm = g.norm()
            if g_norm > 0:
                noise = torch.randn_like(g) * (eps * g_norm)
                g.add_(noise)
        optimizer.step()

    def approx_overlap_metric(self, positions: torch.Tensor, blocks_t: torch.Tensor, bins: int = 64):
        device = positions.device; dtype = positions.dtype
        bw = self.W / bins; bh = self.H / bins
        gx = (positions[:,0] / bw).floor().clamp(0, bins-1).long()
        gy = (positions[:,1] / bh).floor().clamp(0, bins-1).long()
        key = gy * bins + gx
        mass = (blocks_t[:,0] * blocks_t[:,1]) / (bw * bh + 1e-12)
        density = torch.zeros(bins*bins, dtype=dtype, device=device)
        density.index_add_(0, key, mass)
        overflow = (density - 1.0).clamp_min(0)
        return overflow.sum() * (bw * bh)

    def finalize_training(self, start_time, all_positions, positions, name, gamma1, gamma2_triplet):
        end_time = time.time()
        total_training_time = end_time - start_time
        print(f"[{name}] time: {total_training_time:.2f}s")
        if self.loss_history:
            print(f"[{name}] Final Loss = {self.loss_history[-1]:.6f}")
        # save results
        np.save(f'{name}_placement.npy', positions.detach().cpu().numpy())
        print(f"[{name}] saved positions to {name}_placement.npy")


    def train_RBSM(self, optimizer_distance, optimizer_penalty):
        Loss_eval_hat  = loss_function(gamma1=1.0, gamma2=1.0, W=self.W, H=self.H, use_hat=True)
        Loss_eval_area = loss_function(gamma1=1.0, gamma2=1.0, W=self.W, H=self.H, use_hat=False)

        device = next(self.model.parameters()).device
        # 统一迁移
        if isinstance(self.nets_t, torch.Tensor):
            self.nets_t = self.nets_t.to(device, non_blocking=True)
        self.fixed_positions = self.fixed_positions.to(device)
        start_time = time.time()                      
        #self.model().to(device)
        positions = self.model()                      
        # positions = self.model().to(device)       
        scaler = amp.GradScaler("cuda", enabled=torch.cuda.is_available())

        output_dir = PROJECT_DIR / "output"
        output_dir.mkdir(parents=True, exist_ok=True)
        metrics_path = output_dir / "convergence_metrics.csv"
        metric_fields = [
            "k", "step", "global_step", "hpwl_full", "overflow_area",
            "overflow_region_pct", "overflow_cell_pct", "distance_lr",
            "penalty_lr", "elapsed_sec",
        ]
        with metrics_path.open("w", newline="", encoding="utf-8") as stream:
            csv.DictWriter(stream, fieldnames=metric_fields).writeheader()

        def append_metrics(k_value, step_value, global_step, hpwl_value, overflow_area):
            row = {
                "k": k_value,
                "step": step_value,
                "global_step": global_step,
                "hpwl_full": hpwl_value,
                "overflow_area": overflow_area,
                "overflow_region_pct": 100.0 * overflow_area / max(self.W * self.H, 1e-12),
                "overflow_cell_pct": 100.0 * overflow_area / max(self.area, 1e-12),
                "distance_lr": optimizer_distance.param_groups[0]["lr"],
                "penalty_lr": optimizer_penalty.param_groups[0]["lr"],
                "elapsed_sec": time.time() - start_time,
            }
            with metrics_path.open("a", newline="", encoding="utf-8") as stream:
                csv.DictWriter(stream, fieldnames=metric_fields).writerow(row)

        blocks_t   = torch.as_tensor(self.block1, dtype=positions.dtype, device=device)
        fixed_pos  = self.fixed_positions.to(device=device, dtype=positions.dtype)
        nets_t     = torch.as_tensor(self.nets, dtype=torch.long, device=device)
        net_idx_full = torch.arange(self.num_nets, device=device)   # 如未用可删
        full_idx   = torch.arange(self.sb_count, device=device)
        # ==== 全局评估器（不带 γ2 权重） + 评估缓冲 ====
        Loss_eval_hat  = loss_function(gamma1=1.0, gamma2=1.0, W=self.W, H=self.H, use_hat=True)   # HPWL/hat 重叠（未加权）
        Loss_eval_area = loss_function(gamma1=1.0, gamma2=1.0, W=self.W, H=self.H, use_hat=False)  # 几何面积（未加权）

        # 复用 all_positions 的缓冲，避免每次 torch.cat 产生大临时
        all_positions_buf = torch.empty(
            (self.sb_count + self.fixed_positions.shape[0], 2), 
            dtype=positions.dtype, device=device
        )

        # HPWL 分块大小（显存 8GB 建议 4096~8192；不够就再小）
        EVAL_HPWL_CHUNK = 2048

        @torch.no_grad()
        def eval_hpwl_full_chunked():
            """全网 HPWL（精确），GPU 分块累加，防 OOM"""
            # 填充缓冲
            all_positions_buf[:self.sb_count].copy_(positions)
            all_positions_buf[self.sb_count:].copy_(fixed_pos)
            total = 0.0
            N = self.num_nets
            for s in range(0, N, EVAL_HPWL_CHUNK):
                bi = self.net_idx_t[s:s+EVAL_HPWL_CHUNK].to(device)
                val = Loss_eval_hat.calculate_distance_loss(self.nets_t, bi, all_positions_buf)
                total += float(val)
            return total

        @torch.no_grad()
        def eval_overlap_global_unweighted():
            """全局未加权重叠：用空间哈希抽样对，显存友好"""
            i_idx, j_idx = _spatial_hash_pairs(
                positions.detach(),
                torch.as_tensor(self.block1, dtype=positions.dtype, device=device),
                bin_size=None, halo_bins=0,
                BIN_CAP=128, PAIRS_PER_BIN_CAP=2000, MAX_PAIRS_GLOBAL=300_000
            )
            if i_idx.numel() == 0:
                return 0.0, 0.0
            ovl_hat  = Loss_eval_hat.calculate_overlap_loss_pairs(
                positions, blocks_t, i_idx, j_idx, pair_weights=None, use_hat=True)
            ovl_area = Loss_eval_area.calculate_overlap_loss_pairs(
                positions, blocks_t, i_idx, j_idx, pair_weights=None, use_hat=False)
            return float(ovl_hat), float(ovl_area)

        # Full HPWL is expensive. Use 5 or 25 for a faster, coarser curve.
        EVAL_EVERY_STEP = max(1, int(os.environ.get("RBSM_EVAL_EVERY_STEP", "1")))
        ENABLE_EARLY_STOP = os.environ.get("RBSM_ENABLE_EARLY_STOP", "0") == "1"
        RUN_OUTER_ROUNDS = max(1, int(os.environ.get("RBSM_OUTER_ROUNDS", "90")))

        snapshot_dir = output_dir / "layout_snapshots"
        snapshot_dir.mkdir(parents=True, exist_ok=True)
        snapshot_bins = 64
        snapshot_every_outer = max(1, int(os.environ.get("RBSM_SNAPSHOT_EVERY_OUTER", "1")))
        snapshot_max_points = max(1, int(os.environ.get("RBSM_SNAPSHOT_MAX_POINTS", "100000")))
        snapshot_count = min(self.sb_count, snapshot_max_points)
        snapshot_ids_np = np.sort(
            np.random.default_rng(0).choice(self.sb_count, snapshot_count, replace=False)
        )
        snapshot_ids = torch.as_tensor(snapshot_ids_np, dtype=torch.long, device=device)
        np.savez(
            snapshot_dir / "metadata.npz",
            sample_indices=snapshot_ids_np,
            sample_sizes=self.block1[snapshot_ids_np],
            fixed_positions=fixed_pos.detach().cpu().numpy(),
            width=np.asarray(self.W),
            height=np.asarray(self.H),
            bins=np.asarray(snapshot_bins),
        )

        @torch.no_grad()
        def save_layout_snapshot(outer_round, global_step, hpwl_value, overflow_area):
            bw = self.W / snapshot_bins
            bh = self.H / snapshot_bins
            gx = (positions[:, 0] / bw).floor().clamp(0, snapshot_bins - 1).long()
            gy = (positions[:, 1] / bh).floor().clamp(0, snapshot_bins - 1).long()
            density = torch.zeros(snapshot_bins * snapshot_bins, dtype=positions.dtype, device=device)
            density.index_add_(0, gy * snapshot_bins + gx, blocks_t[:, 0] * blocks_t[:, 1] / (bw * bh))

            half_w = blocks_t[:, 0] * 0.5
            half_h = blocks_t[:, 1] * 0.5
            outside = (
                (positions[:, 0] - half_w < 0)
                | (positions[:, 0] + half_w > self.W)
                | (positions[:, 1] - half_h < 0)
                | (positions[:, 1] + half_h > self.H)
            )
            np.savez(
                snapshot_dir / f"iteration_{outer_round:03d}.npz",
                positions=positions.index_select(0, snapshot_ids).detach().cpu().numpy(),
                density=density.reshape(snapshot_bins, snapshot_bins).detach().cpu().numpy(),
                outer_round=np.asarray(outer_round),
                global_step=np.asarray(global_step),
                hpwl_full=np.asarray(hpwl_value),
                overflow_region_pct=np.asarray(
                    100.0 * overflow_area / max(self.W * self.H, 1e-12)
                ),
                outside_pct=np.asarray(100.0 * float(outside.float().mean().cpu())),
            )


        # gamma
        gamma1 = torch.full((self.sb_count,), 0.0, dtype=positions.dtype, device=device)
        g2_i = torch.empty(0, dtype=torch.long,  device=device)
        g2_j = torch.empty(0, dtype=torch.long,  device=device)
        g2_v = torch.empty(0, dtype=positions.dtype, device=device)
        max_pairs_keep = int(2e6)
        if g2_v.numel() > max_pairs_keep:
            topk = torch.topk(g2_v, k=max_pairs_keep, sorted=False)
            g2_i, g2_j, g2_v = g2_i[topk.indices], g2_j[topk.indices], topk.values

        # LR schedulers
        switch_step = int(self.max_setp * 0.9)
        def lr_lambda_distance(k):
            base = 0.5*(1 + math.cos(math.pi*k/self.max_setp))
            return base if k < switch_step else base / math.sqrt(k+1)
        def lr_lambda_penalty(k):
            base = (1 + math.cos(math.pi*k/self.max_setp))
            return base if k < switch_step else base / math.sqrt(k+1)
        scheduler_distance = torch.optim.lr_scheduler.LambdaLR(optimizer_distance, lr_lambda=lr_lambda_distance)
        scheduler_penalty  = torch.optim.lr_scheduler.LambdaLR(optimizer_penalty,  lr_lambda=lr_lambda_penalty)

        # 均匀抽样网
        probs = torch.ones(self.num_nets, dtype=torch.float32, device=device) / max(1, self.num_nets)

        Loss = loss_function(gamma1=1.0, gamma2=1.0, W=self.W, H=self.H,
                            bin_size=None, max_pairs=int(2e6), use_hat=True)

        # monitors
        last_lengths = []
        max_convergences = 3
        convergence_threshold = 1e-5
        self.overlap_threshold = self.area * 0.01
        print(f"overlap_threshold={self.overlap_threshold:.3f}")

        k = 1
        stop_training = False
        batch_size = int(math.ceil(self.num_nets / 25))
        print(f"[train] positions={positions.device}, nets_t={'tensor:'+str(nets_t.device) if isinstance(nets_t, torch.Tensor) else 'numpy'}")

        logdir = "./prof_log"  # 日志目录
        with tprof.profile(
            activities=[tprof.ProfilerActivity.CPU, tprof.ProfilerActivity.CUDA],
            schedule=tprof.schedule(wait=0, warmup=0, active=25, repeat=1),  # 先等1步、热身1步、记录3步，重复2轮
            on_trace_ready=tprof.tensorboard_trace_handler(logdir),         # 写入 TensorBoard
            record_shapes=True,                # 记录张量shape
            profile_memory=True,               # 记录内存分配
            with_stack=False,                  # 需要堆栈可改 True（日志更大）
        ) as prof:
            with tprof.record_function("Full_value"):
                                hpwl_full = eval_hpwl_full_chunked()                     # 全网 HPWL（精确、不采样）
                                ovl_hat_global, ovl_area_global = eval_overlap_global_unweighted()  # 未加权重叠（hat & 几何）
            overflow_area = float(self.approx_overlap_metric(positions, blocks_t).detach().cpu())
            overflow_region_pct = 100.0 * overflow_area / max(self.W * self.H, 1e-12)
            append_metrics(0, -1, 0, hpwl_full, overflow_area)
            save_layout_snapshot(0, 0, hpwl_full, overflow_area)

            print(
                f"HPWL_full={hpwl_full:.3e}  "
                f"OVLhat_global={ovl_hat_global:.3e}  "
                f"OVLarea_global={ovl_area_global:.3e}  "
                f"Overflow={overflow_region_pct:.3f}%"
            )
            while True:
                # --- gamma update (sparse) ---
                with tprof.record_function("Gamma_update"):
                    g1_new, (gi_new, gj_new, gv_new) = self.calculate_gamma_sparse(positions)
                    #gamma1 = torch.maximum(gamma1, 5.0 * g1_new)
                    gamma1 = torch.minimum(gamma1, 0.0 * g1_new)

                if gv_new.numel() > 0:
                    if g2_v.numel() == 0:
                        g2_i, g2_j, g2_v = gi_new, gj_new, 100.0 * gv_new
                    else:
                        key_old = g2_i * self.sb_count + g2_j
                        key_new = gi_new * self.sb_count + gj_new
                        key = torch.cat([key_old, key_new], 0)
                        vi  = torch.cat([g2_v, 100.0 * gv_new], 0)
                        order = torch.argsort(key)
                        key = key[order]; vi = vi[order]
                        keep = torch.ones_like(key, dtype=torch.bool)
                        keep[:-1] = key[1:] != key[:-1]
                        key = key[keep]; vi = vi[keep]
                        g2_i = key // self.sb_count
                        g2_j = key %  self.sb_count
                        g2_v = vi
                

                # --- inner steps ---
                for step in range(25):
                    # distance phase：固定小批
                    N = self.num_nets
                    BATCH = min(batch_size, N, 10000)
                    batch_idx = torch.randperm(N, device=device)[:BATCH]

                    accum = 4                                   # 可按显存/速度调
                    optimizer_distance.zero_grad(set_to_none=True)
                    batch_idx_chunks = batch_idx.chunk(accum)   # 把这次的 nets 均分成 accum 份
                    #with tprof.record_function("HPWL_loss"):
                    with Measure("HPWL_forward", cuda=True):
                        # all_positions = torch.cat([positions, fixed_pos], dim=0)
                        # distance_loss = Loss.calculate_distance_loss(nets_t, batch_idx, all_positions)
                        # distance_loss = self.hpwl_segmented(positions, fixed_pos, batch_idx)
                        # ----- HPWL phase: AMP + 梯度累积 ----
                        distance_loss = positions.new_zeros(())
                        for t in range(accum):
                            bi_t = batch_idx_chunks[t]
                            with amp.autocast("cuda", dtype=torch.float16):   # 或 torch.bfloat16
                                loss_t = self.hpwl_segmented(positions, fixed_pos, bi_t) / accum   
                            distance_loss = distance_loss + loss_t

                    with Measure("HPWL_backward+step", cuda=True):
                        #self.update_gradient(positions, distance_loss, optimizer_distance, k)
                        scaler.scale(distance_loss).backward()
                        scaler.step(optimizer_distance)
                        scaler.update()
                        if torch.cuda.is_available(): torch.cuda.synchronize()

                    # penalty phase (boundary + overlap)
                    #with tprof.record_function("Penalty"):
                    with Measure("Boundary_forward", cuda=True):
                        Loss.gamma1 = gamma1
                        optimizer_penalty.zero_grad(set_to_none=True)
                        boundary_loss = Loss.calculate_boundary_loss1(full_idx, positions, blocks_t)
                        if torch.cuda.is_available(): torch.cuda.synchronize()
                    with Measure("Overlap_forward", cuda=True):
                        if g2_v.numel() > 0:
                            overlap_loss = Loss.calculate_overlap_loss_pairs(
                                positions, blocks_t, i_idx=g2_i, j_idx=g2_j, pair_weights=g2_v, use_hat=True
                            )
                            K = min(100000, g2_v.numel())
                            val, idx = torch.topk(g2_v, k=K, sorted=False)
                            g2_i, g2_j, g2_v = g2_i[idx], g2_j[idx], val
                        else:
                            overlap_loss = Loss.calculate_overlap_loss(
                                full_idx, positions, blocks_t, scope='batch_halo', halo_bins=1
                            )
                        if torch.cuda.is_available(): torch.cuda.synchronize()
                    with Measure("Penalty_backward+step", cuda=True):
                        penalty_loss = boundary_loss + overlap_loss
                        self.update_gradient(positions, penalty_loss, optimizer_penalty, k)
                        if torch.cuda.is_available(): torch.cuda.synchronize()
                    
                    with Measure("Value", cuda=True):
                        # === 日志（每 5 步）===
                        if step % 5 == 0:
                            hpwl_v    = float(distance_loss.detach().cpu())
                            bound_v   = float(boundary_loss.detach().cpu())
                            overlap_v = float(overlap_loss.detach().cpu())
                            pairs_v   = int(g2_v.numel())
                            msg = (f"[k={k:02d} step={step:02d}] "
                                f"HPWL={hpwl_v:.3e}  BND={bound_v:.3e}  OVL={overlap_v:.3e}  pairs={pairs_v}")
                            if torch.cuda.is_available():
                                alloc = torch.cuda.memory_allocated() / 1e9
                                rsrv  = torch.cuda.memory_reserved() / 1e9
                                msg += f"  [cuda mem alloc={alloc:.2f}G reserved={rsrv:.2f}G]"
                            print(msg)
                            with torch.no_grad():
                                total = distance_loss + penalty_loss
                                self.loss_history.append(float(total.detach().cpu()))
                        # === 早停 ===
                        with torch.no_grad():
                            last_lengths.append(float(distance_loss.detach().cpu()))
                            if len(last_lengths) > max_convergences:
                                last_lengths = last_lengths[-max_convergences:]
                            approx_ov_val = float(self.approx_overlap_metric(positions, blocks_t).detach().cpu())

                            conv = False
                            if len(last_lengths) >= max_convergences:
                                ll = np.array(last_lengths, dtype=np.float64)
                                diffs = np.abs(np.diff(ll))
                                ref   = ll[:-1] + 1e-12
                                conv  = np.all(diffs < (convergence_threshold * ref))

                            if ENABLE_EARLY_STOP and conv and approx_ov_val < self.overlap_threshold:
                                print("Early stop: small HPWL change & acceptable overlap.")
                                stop_training = True
                    # === 全局评估：每步都打印（HPWL 全网+不带权重重叠） ===
                    with Measure("Global value", cuda=True):
                        # ==== 每25步全局评估（未加权） ====
                        if (step % EVAL_EVERY_STEP) == 0:
                            with tprof.record_function("Full_value"):
                                hpwl_full = eval_hpwl_full_chunked()                     # 全网 HPWL（精确、不采样）
                                if step == 0:
                                    ovl_hat_global, ovl_area_global = eval_overlap_global_unweighted()

                            overflow_area = approx_ov_val
                            overflow_region_pct = 100.0 * overflow_area / max(self.W * self.H, 1e-12)
                            global_step = (k - 1) * 25 + step + 1
                            append_metrics(k, step, global_step, hpwl_full, overflow_area)

                            msg = (
                                f"[k={k:02d} step={step:02d}] "
                                f"HPWL_full={hpwl_full:.3e}  "
                                f"Overflow={overflow_region_pct:.3f}%"
                            )
                            if step == 0:
                                msg += (
                                    f"  OVLhat_global={ovl_hat_global:.3e}"
                                    f"  OVLarea_global={ovl_area_global:.3e}"
                                )
                            print(msg)

                    prof.step()
                    if stop_training:
                        break

                snapshot_global_step = (k - 1) * 25 + step + 1
                if (step % EVAL_EVERY_STEP) == 0:
                    snapshot_hpwl = hpwl_full
                    snapshot_overflow = approx_ov_val
                else:
                    snapshot_hpwl = eval_hpwl_full_chunked()
                    snapshot_overflow = float(
                        self.approx_overlap_metric(positions, blocks_t).detach().cpu()
                    )
                if (k % snapshot_every_outer) == 0 or stop_training or k == RUN_OUTER_ROUNDS:
                    save_layout_snapshot(
                        k, snapshot_global_step, snapshot_hpwl, snapshot_overflow
                    )
                scheduler_distance.step()
                scheduler_penalty.step()
                k += 1

                # ====  最终全局评估（未加权） ====
                if stop_training or k > RUN_OUTER_ROUNDS:
                    final_hpwl = eval_hpwl_full_chunked()
                    _, final_ovl_area = eval_overlap_global_unweighted()
                    final_overflow_area = float(self.approx_overlap_metric(positions, blocks_t).detach().cpu())
                    final_overflow_pct = 100.0 * final_overflow_area / max(self.W * self.H, 1e-12)
                    print(
                        f"[final] HPWL_full={final_hpwl:.3e}  "
                        f"OVLarea_global={final_ovl_area:.3e}  "
                        f"Overflow={final_overflow_pct:.3f}%  "
                        f"metrics={metrics_path}"
                    )

                    all_positions = torch.cat([positions, fixed_pos], dim=0)
                    self.finalize_training(start_time, all_positions, positions, 'RBSM_sparse',
                                        gamma1, (g2_i, g2_j, g2_v))
                    break

# =========================
# Pipeline / main
# =========================
class TrainingPipeline:
    def __init__(self, node_file, net_file, pl_file, max_epoch, batch_size, max_step):
        # read once to get sb_count
        blocks, sb_count, _, _ = Nodes(node_file).read_nodes()
        self.sb_count = sb_count
        self.node_file = node_file
        self.net_file  = net_file
        self.pl_file   = pl_file
        self.max_epoch = max_epoch
        self.batch_size= batch_size
        self.max_step  = max_step
        scl_file = ISPD2005_CASE_DIR / "adaptec1.scl"
        W, H, bounds = read_scl_bounds(scl_file)
        print(f"W={W}, H={H}")
        self.model = Position(sb_count, W, H)
        self.trainer = TrainModel(self.model, net_file, node_file, pl_file, W, H, max_epoch, batch_size, max_step)

    def train_rbsm(self, lr_distance: float, lr_penalty: float):
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
        self.model.to(device)  # 先搬到目标设备

        optimizer_distance = optim.SGD([self.model.positions], lr=lr_distance)
        optimizer_penalty  = optim.SGD([self.model.positions], lr=lr_penalty)

        self.trainer.train_RBSM(optimizer_distance, optimizer_penalty)

def main():
    seed = int(os.environ.get("RBSM_SEED", "2026"))
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    print(f"RBSM_SEED={seed}")

    # ========== 修改为 ISPD2005 文件路径 ==========
    node_file = ISPD2005_CASE_DIR / "adaptec1.nodes"
    net_file  = ISPD2005_CASE_DIR / "adaptec1.nets"
    pl_file   = ISPD2005_CASE_DIR / "adaptec1.pl"
    # ===============================================

    max_epoch = 40
    batch_size = 20
    # Keep the original 100-round learning-rate horizon. RBSM_OUTER_ROUNDS
    # independently truncates this diagnostic run after round 90.
    max_step = 100

    pipeline = TrainingPipeline(node_file, net_file, pl_file,  max_epoch, batch_size, max_step)
    pipeline.train_rbsm(lr_distance=5, lr_penalty=5)
    #device = 'cuda' if torch.cuda.is_available() else 'cpu'
    #pipeline.model.to(device)


if __name__ == "__main__":
    main()
