# HUAWEI EDA GP：非光滑协议下的修复、重启动与 DREAMPlace 思想迁移实施方案 V8

> 仓库：`Bowen-Yuan/HUAWEI_EDA_GP`  
> 分支：`nonsmooth-gp-v1`  
> 方案审计基线：`6fc2068c76744d346f04102767b75f8f4366a774`  
> 日期：2026-09-13  
> 用途：交给本地模型执行代码修复、契约测试、数值实验、结果分析和 Git 收尾。  
>
> **本文件是执行任务书，不是完成证明。执行前必须重新读取 `AGENTS.md` 和 `you should know about/` 全部 Markdown，重新确认当前 branch / HEAD / dirty worktree。若 HEAD 与本文件不同，以当前源码、CMake、测试和 fresh exact audit 为最终权威。**

---

# 0. 任务总目标

本轮一次完成三个任务：

1. **修复当前项目已经审计出的 bug、数值风险和不合理实现**，但暂时不处理 Linux/macOS SHA-256 跨平台问题。
2. **设计并执行“重启动 restart”数值实验**，检验保留布局、重置优化器/控制器状态后继续优化是否能带来增益；包括 H221-style 等预算因果对照，以及 V7 全链连续重启的实用实验。
3. **从 DREAMPlace 迁移两个与本项目最兼容、且不破坏非光滑协议的思想，分别做成新 module 并做数值实验：**
   - `spectral_bb_gp`：Barzilai–Borwein / spectral step-size，自适应估计步长；
   - `overflow_epsilon_continuation_gp`：把 DREAMPlace 的 overflow-driven smoothing continuation 思想映射为 **epsilon-active 集合半径 continuation**，仅改变非光滑方向 oracle，不改变 exact objective。

本轮禁止把 DREAMPlace 的 weighted-average / log-sum-exp wirelength、electrostatic potential、DCT/Poisson 等平滑 surrogate 重新引入正式目标。

---

# 1. 不可破坏的协议

## 1.1 Exact HPWL 是唯一权威线长

继续使用：

\[
W(x)=\sum_{e\in E}w_e\left[
\max_{p\in e}X_p-\min_{p\in e}X_p+
\max_{p\in e}Y_p-\min_{p\in e}Y_p
\right].
\]

其中：

\[
X_p=x_{v(p)}+\Delta x_p,\qquad
Y_p=y_{v(p)}+\Delta y_p.
\]

要求：

- pin offset 不得丢失；
- net weight 不得丢失；
- epsilon-active 只允许改变方向选择，不得改变报告的 exact HPWL；
- 禁止 WA/LSE 替代最终 HPWL。

## 1.2 Exact rectangle/bin overlap 是唯一权威 density

继续使用：

\[
\rho_b(x)=\frac{\operatorname{occ}_b(x)}{A_b},
\]

\[
D(x)=\sum_b\frac12 A_b[\rho_b(x)-\rho_t]_+^2,
\]

\[
O(x)=
\frac{\sum_b A_b[\rho_b(x)-\rho_t]_+}
{A_{\mathrm{movable}}}.
\]

保持：

- fixed macro 消耗容量且不移动；
- `terminal_NI` 不消耗物理容量；
- canonical 正式 evaluator 为 `512 × 512 / target_density = 1.0`；
- 对外 overflow 使用百分比；
- overflow < 7% 不等于 legalized。

## 1.3 Direction 与 metric 必须分离

允许：

- exact subgradient；
- epsilon-active direction；
- active-face ensemble；
- finite-radius exact secant；
- BB/spectral step-size；
- trust radius；
- coarse capacity candidate；
- exact backtracking；
- exact-audited relocation。

禁止：

- 将平滑 surrogate 的数值作为 final metric；
- 为了性能跳过 stage boundary fresh audit；
- 将历史 DCT/Poisson 重新包装成正式 non-smooth module。

## 1.4 Pure-descent module 的特殊规则

现有五个 non-local pure-descent module 若继续保持 pure-descent 定义，则 exact audit 只能用于：

- telemetry；
- lambda feedback；
- 研究记录。

不得偷偷增加：

- candidate accept/reject；
- rollback；
- best-feasible restore；
- hard overflow gate。

本轮新建的两个 DREAMPlace-inspired module 第一版也采用 **pure-descent default**，避免同时引入新方向、步长和 acceptance 三个变量。若后续需要 exact acceptance，单独做下一轮消融。

---

# 2. 执行前固定检查

本地模型开始修改前执行：

```powershell
Set-Location 'D:\codex_project\HUAWEI_EDA\non-smooth'

git branch --show-current
git rev-parse HEAD
git status --short
git log -5 --oneline

Get-Content .\AGENTS.md
# 然后完整阅读 you should know about 下全部 Markdown
```

必须确认：

```text
branch = nonsmooth-gp-v1
```

但如果 HEAD 已不是 `6fc2068c...`，不要强行回退；先审计新 HEAD 是否已完成本方案中的部分任务。

禁止：

```text
git reset --hard
git checkout -- <user-file>
```

不得覆盖用户已有未提交改动。

---

# 3. 本轮建议提交结构

不要把所有改动堆成一个巨大提交。推荐四个逻辑提交，每个提交同时更新对应 handoff 文档和 `CHANGELOG.md`。

## Commit A — correctness / audit repair

目标：

- 修 finite-radius；
- deterministic canonical density audit；
- nonlocal HPWL context；
- objective/local-query accounting；
- CSV/global iteration；
- canonical `audit` 对齐；
- weighted preconditioner；
- high-degree direction；
- tests。

建议 commit message：

```text
fix: harden exact nonsmooth audit and local oracle contracts
```

## Commit B — restart experiment

只增加必要配置、pipeline、结果摘要，不再改数学 kernel，除非 Commit A 暴露 bug。

建议：

```text
exp: evaluate optimizer restart on h221 and v7 chains
```

## Commit C — DREAMPlace spectral BB module

```text
feat: add nonsmooth spectral BB placement stage
```

## Commit D — DREAMPlace overflow-epsilon continuation module

```text
feat: add overflow-driven epsilon continuation stage
```

如果两个新 module 共用一个很小的 shared utility，可以在 Commit C 先加入通用 primitive，但不能建立第二套 optimizer hierarchy。

---

# 4. Task 1：修复当前项目所有已识别问题

下面按优先级执行。跨平台 SHA-256 **明确不在本轮范围**。

---

# 5. Fix A1：finite-radius oracle 必须只返回真正下降的 secant

当前 `finite_radius_oracle_gp` 会在 `+/-` 两个 probe 都使 density energy 变差时，仍选择“较不差”的那个 probe。

这会产生错误梯度。

## 5.1 正确数学定义

对方向：

\[
u\in\{+e_x,-e_x,+e_y,-e_y\}
\]

定义：

\[
\Delta D_{i,u}(\delta)=D(x+\delta u)-D(x).
\]

只允许：

\[
\Delta D_{i,u}<-\tau_D
\]

进入候选集合。

若某一 axis 没有任何下降方向，则：

\[
g^{oracle}_{i,axis}=0.
\]

若选中某个下降方向 \(u^\star\)，定义：

\[
s_{i,u^\star}
=
\frac{\Delta D_{i,u^\star}}{\delta_{\mathrm{actual}}},
\]

\[
g_i^{oracle}=s_{i,u^\star}u^\star.
\]

因为 \(s<0\)，统一更新：

\[
x\leftarrow x-\eta g
\]

会朝 probe 的下降方向移动。

## 5.2 必须使用实际位移长度

probe 经过 die clamp 后，实际位移可能小于名义 `delta`。

因此：

```cpp
actual_delta = abs(trial_coordinate - original_coordinate);
```

若：

```text
actual_delta <= coordinate_tolerance
```

则该 probe 无效，禁止除以名义 `delta`。

## 5.3 参数校验

增加：

```text
probe_count >= 0
radius_bin_scale > 0
descent_tolerance >= 0
```

建议默认：

```json
{
  "probe_count": 32,
  "radius_bin_scale": 1.0,
  "descent_tolerance": 1e-12
}
```

`descent_tolerance` 是 density-energy 数值容差，不是 feasibility tolerance。

---

# 6. Fix A2：finite-radius 不得继续复制整个 Database 做 full density probe

当前最浪费的路径是：

```text
Database trial = db0
→ move one node
→ construct ExactOverlapDensity
→ full evaluate
```

项目已有：

```cpp
ExactOverlapDensity::evaluate_move(...)
```

它已经能返回 touched-bin 的 exact：

- area changes；
- overflow-area delta；
- energy delta。

因此 finite-radius 必须改为：

```cpp
const DensityMove move =
    density.evaluate_move(node_id, new_x, new_y);

const Real delta_energy = move.energy_delta;
```

**probe 阶段不得 commit**。

只在统一 optimizer 真正更新 live layout 后，下一轮重新构造/刷新 occupancy。

目标复杂度从：

\[
O(\text{probe count}\times N_{\mathrm{movable}})
\]

降为接近：

\[
O(\text{probe count}\times \text{touched bins}).
\]

## 6.1 重要状态条件

`evaluate_move()` 依赖 `ExactOverlapDensity::occupancy_` 对应当前 live layout。

因此在 callback 进入前，common engine 必须已经对当前布局调用 density evaluation / rebuild。

新增契约测试防止未来调用顺序破坏这一点。

---

# 7. Fix A3：补全 local oracle 成本统计

当前 `objective_evaluations` 只统计 full exact HPWL/density evaluation，finite-radius 内部 probe 没有被记录。

不能继续把 128 次 local probe 和 0 次 local probe 写成同一个成本。

## 7.1 扩展统计结构

建议新增：

```cpp
struct AuxiliaryQueryStats {
    int local_density_queries = 0;
    int local_hpwl_queries = 0;
};
```

将：

```cpp
using AuxiliaryGradient = std::function<void(...)>;
```

调整为：

```cpp
using AuxiliaryGradient =
    std::function<AuxiliaryQueryStats(...)>;
```

五个现有 nonlocal module：

- 无 local exact query 的返回 `{0, 0}`；
- finite-radius 每次 `evaluate_move()` 增加 `local_density_queries`；
- 未来 wire-aware finite oracle 可增加 `local_hpwl_queries`。

## 7.2 StageStats

在保持旧 aggregate 初始化可编译的前提下，追加：

```cpp
int local_density_queries = 0;
int local_hpwl_queries = 0;
```

`objective_evaluations` 继续表示 full-layout exact evaluations，不要把 local query 混进去。

## 7.3 trajectory / experiment.md

在 stage row 追加：

```text
local_density_queries
local_hpwl_queries
```

`experiment.md` 总结也累计输出。

这样以后可以比较：

```text
quality / wall_seconds / full exact evals / local delta queries
```

而不是只有一个模糊的 `objective_evaluations`。

---

# 8. Fix A4：nonlocal engine 计算了 exact HPWL，却把值丢掉

当前逻辑中 exact HPWL evaluation 返回值没有保存，导致：

```text
LambdaController.initialize(... hpwl = 0)
LambdaController.update(... hpwl = 0)
AuxiliaryContext.hpwl = 0
```

修复为：

```cpp
const Real exact_hpwl =
    hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
```

之后传入：

```cpp
controller.initialize(..., exact_hpwl, metrics.overflow);
controller.update(iteration, exact_hpwl, metrics.overflow);

AuxiliaryContext ctx{
    iteration,
    lambda,
    wx,
    wy,
    dx,
    dy,
    exact_hpwl,
    metrics
};
```

增加 synthetic test：

```text
callback 读取 ctx.hpwl
assert ctx.hpwl == independently evaluated exact HPWL
```

---

# 9. Fix A5：canonical density audit 必须 deterministic

当前 movable occupancy 使用 OpenMP atomic：

```cpp
#pragma omp atomic update
occupancy_[index] += area;
```

数学定义没错，但浮点加法顺序不确定。

当 V7 已经到：

```text
overflow ≈ 6.9999995%
```

时，canonical feasibility 不应依赖 thread interleaving。

## 9.1 不要求把所有搜索计算都改成串行

保留两种路径：

```cpp
enum class DensityAccumulationMode {
    FastParallel,
    DeterministicCanonical
};
```

或者等价 API。

### FastParallel

用于：

- 搜索方向；
- 非关键 telemetry；
- 允许微小浮点差异的内部计算。

### DeterministicCanonical

用于：

- `exact_audit()`；
- final-selected；
- candidate acceptance；
- best-feasible selection；
- overflow-cap recovery 的精确决策；
- 任何会决定“是否满足 7%”的地方。

## 9.2 第一版 deterministic 实现优先正确性

最简单可靠的实现：

```text
fixed occupancy 按固定 node 顺序
movable occupancy 按 db.movable_ids 固定顺序
每个 rectangle 内 bin loop 固定 y/x 顺序
最终 bin reduction 固定 index 顺序
```

即 canonical metrics-only rebuild 可以先使用单线程固定顺序。

不要为了优化过早引入复杂并行 deterministic reduction。

先测 runtime。

如果 canonical audit runtime 过高，再考虑：

- deterministic contribution buffer；
- stable sort by `(bin_id, node_id)`；
- fixed-order segmented reduction。

但第一版不要求。

## 9.3 exact local recovery 的 occupancy

凡是后续通过 `evaluate_move()` / `commit_move()` 做 exact accept/reject 的模块，在开始该模块前必须从 deterministic occupancy rebuild 起步。

推荐提供：

```cpp
density.rebuild_occupancy(DensityAccumulationMode::DeterministicCanonical);
```

避免 recovery 的 local delta 基于 nondeterministic baseline occupancy。

---

# 10. Fix A6：trajectory.csv schema 错位和 global_iteration

## 10.1 修复列数

当前 iteration row 在 `hpwl` 后插入过多空列。

必须建立单一 schema，并用单元测试验证：

```text
每一行 split(',').size() == header_column_count
```

建议 header：

```text
record_type
stage_index
module
global_iteration
stage_iteration
hpwl
hpwl_before
hpwl_after
overflow_percent_before
overflow_percent_after
density_energy
max_density
iterations
accepted
rejected
objective_evaluations
local_density_queries
local_hpwl_queries
wall_seconds
```

iteration row：

- `hpwl` 有值；
- `hpwl_before/hpwl_after` 留空；
- `overflow_percent_before/after` 可统一写 post-update 值，或定义一个 `overflow_percent` 新字段；本轮为了兼容，保留当前两个字段但确保列位置正确；
- query/evaluation 使用 **累计值**，不是硬编码每轮 `"4"`。

## 10.2 global_iteration

现在每个 stage 的 `global_iteration` 都等于 `stage_iteration`，会重置。

runner 维护：

```cpp
int global_iteration_counter = 0;
```

每次收到 iteration telemetry：

```cpp
tagged["stage_iteration"] = stage_iteration;
tagged["global_iteration"] = global_iteration_counter++;
```

如果某 module 一轮发多条 telemetry，则应由 module 明确传 `iteration_event_kind`；当前模块是一轮一条，可以直接累计。

## 10.3 nonlocal module 名称

当前 common wrapper 把 telemetry module 写成：

```text
nonlocal
```

应改为具体 registry module：

```text
finite_radius_oracle_gp
density_charge_gp
...
```

建议：

```cpp
nonlocal::run(context, config, "finite_radius_oracle_gp", aux);
```

---

# 11. Fix A7：`nsgp audit` 与 pipeline canonical loading 对齐

当前 pipeline：

```text
read Bookshelf
→ load placement
→ clamp movable
→ exact audit
```

而 `audit` 命令可能直接 audit 未 clamp 的 placement。

H221 checkpoint 已经证明这会产生：

```text
pipeline canonical overflow ≈ 7.8193%
raw audit artifact ≈ 11.52%
```

因此默认 `audit` 必须：

```cpp
load_bookshelf_placement(...)
clamp_movable(db)
exact_audit(...)
```

并输出：

```text
audit_mode=canonical_clamped
```

可选地增加诊断：

```text
--raw-unclamped-audit
```

但若加入，必须明确输出：

```text
NONCANONICAL diagnostic
```

如果不需要，不要增加这个 flag。

## 11.1 audit density config

`audit` 默认使用 canonical：

```text
512 × 512
target_density = 1.0
```

同时允许：

```text
--pipeline <pipeline.json>
```

读取其中 `density_grid`，用于审计非 canonical smoke 结果。

不传 `--pipeline` 时仍使用 canonical，不改变常规命令。

---

# 12. Fix A8：preconditioner 使用 incident pin/net weight，而不是纯 pin count

当前：

\[
P_i =
\max\{1,d_i+\lambda A_i\}.
\]

改为：

\[
P_i =
\max\left\{
1,
\sum_{p\in i} w_{e(p)}
+
\lambda A_i
\right\}.
\]

这与 DREAMPlace 的 `sum_pin_weights_in_nodes + density_weight * node_area` 思想一致，同时仍然是 exact non-smooth direction preconditioning。

## 12.1 Database 新字段

建议：

```cpp
std::vector<Real> node_pin_weight_sum;
```

在所有 net weight 读取完成后：

```cpp
node_pin_weight_sum.assign(nodes.size(), 0.0);

for (const Net& net : nets) {
    for (pin in net) {
        node_pin_weight_sum[pin.node] += net.weight;
    }
}
```

对于全部 `net.weight=1` 的 case，数值应尽量接近旧 `node_pin_count`。

## 12.2 backward-compatible ablation

为了判断这一变化是否影响 V7，`PlaceConfig` 可提供：

```text
preconditioner = "weighted_pin_sum"   # 新默认
preconditioner = "pin_count"          # 历史回归
```

但不要长期保留复杂 policy tree；若验证 weighted 明显更合理且不破坏历史 replay，可在后续清理旧 mode。

---

# 13. Fix A9：high-degree net 不应在 gradient 中完全消失

当前超过 `degree_limit` 的 net：

- exact HPWL 继续计入；
- gradient 直接置零。

这会导致 optimizer 对大网完全失明。

## 13.1 新策略

保留历史模式：

```text
high_degree_mode = "ignore"
```

新增推荐模式：

```text
high_degree_mode = "exact_extrema"
```

对于：

```text
degree > degree_limit
```

不再使用 epsilon 邻域 active set，而只使用 **epsilon=0 exact extremal-face subgradient**：

- min-x ties 平分 `-w_e`；
- max-x ties 平分 `+w_e`；
- y 同理。

这仍然是 HPWL 的合法非光滑次梯度，不引入 smooth surrogate。

## 13.2 为什么不用 taper smooth weighting

本轮不要用：

```text
(limit / degree)^p
```

之类人为连续权重作为第一实现。

`exact_extrema` 更符合当前协议，且数学解释更清楚。

## 13.3 兼容性

历史 V7 control 第一轮仍跑：

```text
high_degree_mode = ignore
```

复现修复后 baseline。

之后单独做：

```text
ignore vs exact_extrema
```

小消融。

不要把 high-degree 新策略与 restart / BB / epsilon continuation 同时改变。

---

# 14. Fix A10：扩大 unit tests

现有 `test_numeric_contracts.cpp` 覆盖不足。

新增/扩展以下测试。

## 14.1 Density delta contract

随机或手工小布局：

\[
D(x')-D(x)
=
\texttt{evaluate_move.energy_delta}
\]

\[
O_{\mathrm{area}}(x')-O_{\mathrm{area}}(x)
=
\texttt{evaluate_move.overflow_area_delta}
\]

要求 full deterministic reevaluation 与 incremental delta 在严格容差内一致。

## 14.2 Group move contract

同理验证：

```cpp
evaluate_group_move(...)
```

## 14.3 HPWL local affected-net delta

如果将 `affected_hpwl_delta` 抽成 shared exact local primitive，则验证：

```text
local delta == full exact HPWL(after)-full exact HPWL(before)
```

包括：

- pin offsets；
- fixed node；
- weighted net；
- tie extrema；
- high-degree exact-extrema。

## 14.4 deterministic audit

同一布局：

```text
threads = 1
threads = 2
threads = 8
```

canonical audit 要求：

```text
HPWL identical
density_energy bitwise identical or explicitly exact-equal
overflow identical
max_density identical
```

如果 canonical path 是串行，则该测试应天然通过。

## 14.5 finite-radius contracts

必须覆盖：

### FR1 — both probes worse

```text
+x deltaD > 0
-x deltaD > 0
→ auxiliary x = 0
```

### FR2 — one direction improves

```text
+x deltaD < 0
→ gradient sign must make x -= delta move right
```

### FR3 — negative direction improves

验证符号相反。

### FR4 — boundary clamp

actual displacement < nominal probe radius 时，secant denominator 使用 actual displacement。

### FR5 — evaluate_move vs old full probe

在 synthetic case 上两种计算得到相同 `delta_energy`。

## 14.6 CSV schema test

构造小型 `ExperimentLog`，写 stage + iteration row，验证：

- 所有 row 列数相同；
- overflow 位于正确列；
- local query 字段正确；
- `global_iteration` 单调。

## 14.7 audit clamp contract

创建超 die 的 movable placement：

```text
pipeline loading audit == nsgp audit canonical result
```

---

# 15. Task 1 完成判据

Commit A 只有满足全部条件才算完成：

```text
[ ] cmake configure/build 成功
[ ] ctest 全通过
[ ] repaired V7 control 可运行
[ ] repaired H221 control 可运行
[ ] deterministic canonical audit thread-independent
[ ] finite-radius 不再 full-copy Database
[ ] no-descent finite probe 输出 0
[ ] ctx.hpwl 非 0 且等于 exact HPWL
[ ] trajectory 每行列数一致
[ ] global_iteration 跨 stage 单调
[ ] audit 默认 clamp
[ ] objective/local query accounting 可区分
[ ] weighted preconditioner contract 测试通过
[ ] high-degree exact-extrema contract 测试通过
[ ] retention contract 通过
```

---

# 16. 修复后必须建立新的 control baseline

由于 deterministic audit 和若干方向工程问题修复后，历史数值只能作为参考，不能作为当前严格对照。

## 16.1 H221 canonical input

使用历史：

```text
h221_dct_poisson_512/best.pl
```

期望 SHA-256：

```text
9d0bd99966c3d679780ea73f0c28d9b3eb4163661fe679dfd85d7d32277e79b2
```

加载 + clamp 后历史 exact 参考：

```text
HPWL ≈ 109,837,622.009
overflow ≈ 7.8193%
```

执行前必须重新计算 hash，不要只相信本文。

## 16.2 repaired V7 control

创建明确、冻结的参数文件，不再依赖临时实验 JSON：

```text
modules/exact_recovery/params/v8_v7_cap15.json
modules/exact_joint_gp/params/v8_v7_bridge_control.json
modules/exact_joint_gp/params/v8_v7_retighten_control.json
modules/exact_recovery/params/v8_v7_polish4.json
modules/exact_recovery/params/v8_v7_polish1.json

framework/params/pipelines/v8_v7_h221_control.json
```

历史 V7 结构：

```text
H221 best.pl
→ exact_recovery cap15
→ exact_joint_gp bridge
→ exact_joint_gp retighten
→ exact_recovery polish
→ exact_recovery polish
```

历史关键设置：

```text
bridge:
  iterations = 120
  optimizer = Adam
  lr = 0.002
  maximum_delta = 4.0
  density_weight_scale = 0.5
  net_batch.weight = 1
  net_batch.degree_limit = 32
  hpwl epsilon = 125
  active_power = 4
  selector cap = 0.15
  lambda policy = dreamplace

retighten:
  iterations = 180
  optimizer = Adam
  lr = 0.001
  maximum_delta = 4.0
  density_weight_scale = 4.75
  net_batch.weight = 1
  net_batch.degree_limit = 32
  hpwl epsilon = 125
  active_power = 4
  selector cap = 0.07
  lambda policy = dreamplace
```

H252 / polish recovery 参数以历史脚本和当前 V7 证据为准。

历史旧结果：

```text
87,599,730.55 @ 6.9999995%
```

只作为 regression reference。

**Commit A 后第一次重新执行得到的结果，定义为 `V8_REPAIRED_V7_CONTROL`。**

---

# 17. Task 2：Restart 数值实验

## 17.1 Restart 定义

本轮 restart 不改变布局：

\[
x^{restart}_0 = x^{previous}_{last}.
\]

只重置 solver state：

```text
optimizer moments
optimizer age
lambda controller state（full-stage restart 时）
direction memory / bundle（若有）
step controller streaks
```

禁止：

```text
rollback 到旧 best
重新随机初始化位置
重新读取 raw Bookshelf
使用 hidden checkpoint
```

模块边界本来就要求隐含 solver state 重置，所以最简单的 full restart 是在 pipeline 中重复同一 module。

---

# 18. Restart Experiment R1：H221-style 等预算因果实验

这是判断“restart 本身是否有效”的**主因果实验**。

因为直接跑 V7 ×3 会同时增加预算和重复 pressure schedule，不能单独归因于 restart。

## 18.1 输入

同一 H221 checkpoint：

```text
SHA = 9d0bd999...
canonical 512×512/1.0
```

threads 固定：

```text
threads = 1
```

用于最大限度对齐历史 V5 H221 screen。

## 18.2 配置

以当前 `h221_trajectory_like` 思想为基础，冻结一份 repaired config：

```text
optimizer = AMSGrad 或历史最优配置
learning_rate = 0.002
beta1 = 0.90
beta2 = 0.99
trajectory lambda
batch acceptance = off
exact HPWL/density
```

先做两个主要 variant。

### R1-CONTINUOUS

```text
exact_joint_gp(iterations=150)
```

一个 stage，optimizer/lambda state 连续 150 步。

### R1-RESTART-3x50

```text
exact_joint_gp(iterations=50)
→ exact_joint_gp(iterations=50)
→ exact_joint_gp(iterations=50)
```

每 50 步保留布局但自然重置：

```text
optimizer + lambda controller
```

总 GP budget 完全相同：

```text
150 iterations
```

这是真正公平的 restart test。

## 18.3 可选二阶段诊断：optimizer-only restart

只有 R1-RESTART-3x50 与 continuous 出现 material 差异后才实现。

新增最小参数：

```json
"optimizer": {
  "restart_iterations": [50, 100]
}
```

要求：

- 只 reset optimizer；
- lambda controller 连续；
- 默认空数组，不改变旧行为；
- restart event 写 trajectory。

得到：

### R1-OPT-ONLY

```text
single stage 150
optimizer restart at 50,100
lambda continuous
```

这样可以分辨：

```text
收益来自 optimizer moment reset
还是来自 lambda trajectory 重新开始
```

如果第一阶段没有信号，不要增加这个功能。

---

# 19. Restart Experiment R2：V7 full-chain 连续重启

这是实用实验，不是严格因果实验。

## 19.1 Variant

### R2-V7x1

```text
[V7 five-stage chain] ×1
```

即 repaired control。

### R2-V7x2

```text
[V7 five-stage chain]
→ [V7 five-stage chain]
```

第二轮直接接第一轮最后布局，所有 module state 自然重置。

### R2-V7x3

```text
[V7 five-stage chain] ×3
```

总 GP budget：

```text
x1 = 300
x2 = 600
x3 = 900
```

重点看：

```text
每一轮 V7 restart 的 marginal feasible HPWL gain
```

而不是简单说“900 步比 300 步好”。

## 19.2 每 cycle 记录

```text
cycle index
cycle input HPWL
cycle input overflow
after cap15
after bridge
after retighten
after polish1
after polish2
cycle final HPWL
cycle final overflow
cycle runtime
```

计算：

\[
\Delta W_c=W_{c,\mathrm{final}}-W_{c,\mathrm{input}}.
\]

若：

```text
cycle2 / cycle3 仍有明显负 ΔW 且 final overflow <7%
```

说明 V7 作为 restart macro-cycle 仍有继续价值。

若：

```text
cycle2 很小
cycle3 ≈ 0 或恶化
```

说明第一轮已基本耗尽该 schedule。

---

# 20. Restart Experiment R3：V7 前半段 / H221 前半段的短 restart screen

为了对应“连续运行三个 V7 / H221 前半部分”的直觉，增加一个低成本 screen。

## 20.1 V7 front-half

定义：

```text
cap15 recovery
→ bridge 120
```

比较：

### R3-A

```text
cap15
→ bridge 360 continuous
```

### R3-B

```text
cap15
→ bridge120
→ bridge120
→ bridge120
```

总 bridge budget 360。

注意：R3-B 每 120 步重置 optimizer/lambda，但不重复 cap15。

这是比 “完整 V7 ×3” 更干净的 low-pressure restart 测试。

## 20.2 H221 trajectory front-half

若 R1 使用 150 步，则额外看前 75 步：

```text
75 continuous
vs
3 × 25 restart
```

用于判断 restart 收益是否只出现在后期。

---

# 21. Restart 结论判据

以最终 canonical feasible 结果为主。

定义：

### Positive

相同 iteration budget 下：

\[
W_{restart}<W_{continuous}
\]

且：

```text
final overflow < 7%
```

并达到至少：

```text
0.10% HPWL improvement
```

### Strong positive

```text
>= 0.30% HPWL improvement
```

### Neutral

```text
|ΔHPWL| < 0.10%
```

### Negative

restart 更差，或导致 final infeasible。

同时报告：

```text
best-feasible
last
stage final
```

不要只挑中间最好点。

---

# 22. Task 3A：新增 DREAMPlace-inspired `spectral_bb_gp`

## 22.1 迁移的不是平滑 objective，而是步长估计思想

DREAMPlace Nesterov optimizer 中值得迁移的是：

\[
s_k=x_k-x_{k-1},
\]

\[
y_k=g_k-g_{k-1},
\]

以及 BB / local Lipschitz 型步长：

\[
\alpha_{BB1}
=
\frac{s_k^Ts_k}{s_k^Ty_k},
\]

\[
\alpha_{BB2}
=
\frac{s_k^Ty_k}{y_k^Ty_k},
\]

\[
\alpha_L
=
\frac{\|s_k\|_2}{\|y_k\|_2}.
\]

本项目只迁移这一 **step-size adaptation**。

不迁移：

```text
WA wirelength
LSE wirelength
electric potential objective
smooth Nesterov objective model
```

---

# 23. `spectral_bb_gp` 模块定位

新增：

```text
modules/spectral_bb_gp/
├─ code/
│  ├─ module.cpp
│  ├─ spectral_bb_gp.cpp
│  └─ spectral_bb_gp.hpp
└─ params/
   ├─ smoke.json
   ├─ h221_primary.json
   └─ v7_bridge_primary.json
```

CMake：

```text
target_sources(nsgp ...)
```

如 core 搜索实现独立可测试，可加入 `nsgp_numeric_kernel`。

registry：

```cpp
modules::register_spectral_bb_gp(registry);
```

---

# 24. `spectral_bb_gp` 的方向必须保持 non-smooth exact

每轮：

1. exact HPWL epsilon-active direction；
2. exact overlap density direction；
3. shared LambdaController；
4. weighted pin/area preconditioner；
5. spectral step-size；
6. unconditional update；
7. die clamp；
8. post-update exact telemetry。

方向：

\[
g_k
=
P^{-1}
\left(
g_W^{(\varepsilon_W)}
+
\lambda_k g_D^{(\varepsilon_D)}
\right).
\]

其中 exact HPWL/density metric 仍用 epsilon=0 canonical evaluator 报告。

---

# 25. BB 在非光滑问题中的保护

不能直接照搬 smooth BB。

## 25.1 有效性条件

计算：

\[
s^Ty.
\]

只有：

```text
sTy > curvature_epsilon
yTy > gradient_delta_epsilon
all finite
```

才使用 BB。

否则：

```text
fallback to base step
```

## 25.2 推荐 step

第一版推荐 BB2：

\[
\alpha_{BB2}
=
\frac{s^Ty}{y^Ty}.
\]

因为通常更保守。

再计算：

\[
\alpha_L=\frac{\|s\|}{\|y\|}.
\]

使用：

\[
\alpha_k
=
\operatorname{clip}
\left(
\min(\alpha_{BB2},c_L\alpha_L),
\alpha_{\min},
\alpha_{\max}
\right).
\]

建议：

```text
c_L = 1.0
```

同时增加 growth limiter：

\[
\alpha_k\le g_{\max}\alpha_{k-1}
\]

例如：

```text
max_step_growth = 2.0
```

避免 active-set 切换后突然放大。

## 25.3 trust displacement

最终仍限制：

```text
maximum_delta_bins
```

BB 只能调整 scalar learning rate，不能绕过 per-coordinate displacement cap。

---

# 26. BB restart / fallback 规则

出现任一情况：

```text
sTy <= eps
alpha non-finite
alpha < alpha_min
alpha > raw_sanity_max
||y|| extremely small
```

则：

```text
spectral state reset
alpha = base_learning_rate
```

注意：

**这不是 candidate rollback。**

只是 step estimator restart，live layout 不回滚，符合 pure-descent。

记录：

```text
bb_valid
bb_restart
alpha_bb
alpha_lip
alpha_used
sTy
yTy
```

---

# 27. `spectral_bb_gp` 第一轮 JSON

示意：

```json
{
  "iterations": 100,
  "hpwl_direction": {
    "epsilon": 125.0,
    "active_power": 4.0,
    "degree_limit": 100,
    "high_degree_mode": "exact_extrema"
  },
  "density_direction": {
    "epsilon": 0.0,
    "active_power": 1.0
  },
  "lambda": {
    "policy": "trajectory",
    "density_weight_scale": 1.0,
    "update_interval": 2,
    "stop_overflow": 0.07
  },
  "spectral_step": {
    "base_learning_rate": 0.002,
    "method": "bb2",
    "alpha_min_ratio": 0.125,
    "alpha_max_ratio": 8.0,
    "max_step_growth": 2.0,
    "curvature_epsilon": 1e-14
  },
  "maximum_delta_bins": 0.25,
  "preconditioner": "weighted_pin_sum"
}
```

具体单位必须在 module README / handoff 文档说明。

---

# 28. `spectral_bb_gp` 契约测试

新增：

```text
tests/test_dreamplace_inspired_contracts.cpp
```

至少：

### BB1 positive curvature

给定手工 `s,y`，验证公式。

### BB2 negative curvature

```text
sTy <= 0
→ fallback
```

### BB3 zero y

不得 NaN/Inf。

### BB4 clipping

`alpha_used` 始终在 bounds。

### BB5 growth limit

相邻 step 不超过配置倍数。

### BB6 no rollback

即使 post-update exact metric 变差，pure-descent 第一版仍保留更新。

### BB7 exact metric invariance

改变 BB 参数只改变位置轨迹，不改变 evaluator 定义。

---

# 29. `spectral_bb_gp` 数值实验

## E-BB-1：H221 screen

输入同 H221。

对照总 budget 100：

```text
BB0: exact_joint_gp + SGD / normalized-SGD fixed step
BB1: spectral_bb_gp BB2
BB2: spectral_bb_gp BB1
```

为了隔离 step-size，方向、lambda、preconditioner 必须一致。

若 fixed-step SGD 明显太弱，再增加：

```text
historical Adam repaired control
```

作为 performance reference，但不能用 Adam vs BB 得出纯 step-size 因果结论。

## E-BB-2：300-step H221

只有 100-step BB2 出现正信号再运行：

```text
300 continuous
```

和 repaired V6 / V7 参考比较。

## E-BB-3：V7 bridge replacement

只替换 V7 的 bridge：

```text
cap15
→ spectral_bb_gp bridge120
→ original retighten180
→ polish×2
```

其它全部不变。

若正向，再测试：

```text
bridge + retighten 都换 spectral_bb_gp
```

不要一开始两个阶段都替换，否则无法定位收益来源。

---

# 30. Task 3B：新增 DREAMPlace-inspired `overflow_epsilon_continuation_gp`

## 30.1 迁移思想

DREAMPlace 会根据 overflow 动态改变 wirelength smoothing `gamma`。

本项目不能使用 smooth WA/LSE。

因此迁移为：

> overflow 高时扩大 epsilon-active 邻域，获得更宽的非光滑 active-face 信息；接近目标 overflow 时逐步缩小 epsilon，回到精确 active face。

改变的是：

```text
direction oracle radius
```

不改变：

```text
exact HPWL
exact density
final overflow
```

---

# 31. 模块结构

新增：

```text
modules/overflow_epsilon_continuation_gp/
├─ code/
│  ├─ module.cpp
│  ├─ overflow_epsilon_continuation_gp.cpp
│  └─ overflow_epsilon_continuation_gp.hpp
└─ params/
   ├─ smoke.json
   ├─ h221_primary.json
   └─ v7_primary.json
```

registry：

```cpp
register_overflow_epsilon_continuation_gp(...)
```

---

# 32. epsilon schedule

定义：

```text
overflow_low
overflow_high
```

例如：

```text
overflow_low  = 0.07
overflow_high = 0.15
```

标准化：

\[
q_k
=
\operatorname{clip}
\left(
\frac{O_k-O_{low}}
{O_{high}-O_{low}},
0,1
\right).
\]

使用 smoothstep 仅作为**参数调度函数**：

\[
h(q)=q^2(3-2q).
\]

这不是 objective surrogate，因此允许。

定义：

\[
\varepsilon_W(k)
=
\varepsilon_{W,\min}
+
h(q_k)
(\varepsilon_{W,\max}-\varepsilon_{W,\min}),
\]

\[
\varepsilon_D(k)
=
\varepsilon_{D,\min}
+
h(q_k)
(\varepsilon_{D,\max}-\varepsilon_{D,\min}).
\]

注意：

- epsilon 只进入 active-set direction；
- exact HPWL / overflow 每轮仍 fresh exact 计算；
- final audit 完全不使用 epsilon。

---

# 33. epsilon 单位

为避免当前代码中 absolute coordinate / bin-scale 混淆，新 module 参数显式使用：

```text
*_epsilon_bins
```

实际：

\[
\varepsilon
=
\text{epsilon_bins}
\times
\min(\text{bin_width},\text{bin_height}).
\]

如果要复现历史 `hpwl_epsilon=125`，提供显式：

```text
epsilon_unit = "absolute"
```

但第一版最好统一 bins，并在实验中重新校准。

---

# 34. 推荐第一组 epsilon

不要大矩阵搜索。

H221 初筛只做三档：

### EC0 fixed exact-like control

```text
wire eps bins = fixed baseline
density eps bins = 0
```

### EC1 conservative continuation

```text
wire: 0.25 → 1.0 bins
density: 0 → 0.5 bins
```

### EC2 wider continuation

```text
wire: 0.25 → 2.0 bins
density: 0 → 1.0 bins
```

方向上：

```text
overflow high → epsilon large
overflow near 7% → epsilon small
```

不要引入 10+ 个网格参数。

---

# 35. continuation 与 optimizer moment

epsilon 大幅改变时，Adam/AMSGrad 的旧 moment 可能不再代表当前 active face。

因此增加显式、可记录的 restart 规则：

当：

\[
\frac{\varepsilon_{old}+\epsilon_0}
{\varepsilon_{new}+\epsilon_0}
\]

跨越配置阈值，例如 2 倍，或 schedule 从一个离散 band 进入下一 band 时：

```text
optimizer reset
```

第一版推荐更简单的 band：

```text
q > 0.66      → wide
0.33<q<=0.66 → medium
q<=0.33      → narrow
```

band 改变：

```text
reset optimizer
```

lambda controller 不重置。

记录：

```text
epsilon_wire
epsilon_density
epsilon_band
epsilon_restart_count
```

---

# 36. `overflow_epsilon_continuation_gp` 第一轮 JSON

示意：

```json
{
  "iterations": 100,
  "optimizer": {
    "name": "amsgrad",
    "learning_rate": 0.002,
    "maximum_delta": 0.25,
    "beta1": 0.90,
    "beta2": 0.99
  },
  "epsilon_continuation": {
    "overflow_low": 0.07,
    "overflow_high": 0.15,
    "wire_min_bins": 0.25,
    "wire_max_bins": 1.0,
    "density_min_bins": 0.0,
    "density_max_bins": 0.5,
    "schedule": "smoothstep",
    "restart_on_band_change": true
  },
  "lambda": {
    "policy": "trajectory",
    "update_interval": 2,
    "stop_overflow": 0.07
  },
  "preconditioner": "weighted_pin_sum"
}
```

---

# 37. epsilon continuation 契约测试

### EC-C1 monotonicity

overflow 增大：

```text
epsilon 不得减小
```

### EC-C2 endpoint

```text
O <= low → epsilon = min
O >= high → epsilon = max
```

### EC-C3 metric invariance

同一 placement：

```text
HPWL(eps=0) == HPWL(eps>0) returned exact value
```

仅 gradient 不同。

density metric 同理。

### EC-C4 restart event

跨 band 只 reset optimizer，不移动/回滚布局。

### EC-C5 final exact

stage final metrics 必须来自 canonical epsilon=0 deterministic audit。

---

# 38. epsilon continuation 数值实验

## E-EC-1：H221 100-step screen

固定：

```text
optimizer = AMSGrad
lr = 0.002
max delta = repaired control
lambda identical
preconditioner identical
```

只改变：

```text
epsilon schedule
```

比较：

```text
EC0 fixed
EC1 conservative
EC2 wider
```

晋级规则：

```text
best feasible HPWL 改善 >= 0.10%
或相同 HPWL 下 overflow recovery 明显更快
```

## E-EC-2：H221 300-step

只跑 E-EC-1 最佳者。

## E-EC-3：V7 bridge replacement

```text
cap15
→ overflow_epsilon_continuation_gp(bridge120)
→ original retighten180
→ polish×2
```

若正向，再替换 retighten。

---

# 39. 两个 DREAMPlace module 的 fusion 实验

只有 BB 和 epsilon continuation **各自独立实验都为正**，才创建第三阶段组合。

不要直接新建第三 module。

可以先在：

```text
overflow_epsilon_continuation_gp
```

中选择：

```text
step_policy = spectral_bb
```

前提是 spectral step primitive 已抽成小型 shared utility。

或者反过来让 `spectral_bb_gp` 支持 dynamic epsilon schedule。

但第一轮必须先独立，避免无法归因。

Fusion 只跑：

```text
H221 100
H221 300
V7 bridge replacement
```

不做大网格。

---

# 40. 实验矩阵总表

## Phase A — 修复回归

| ID | 输入 | 变化 | Budget | 目的 |
|---|---|---|---:|---|
| A0 | H221 | repaired H221 control | 100 | 建新 baseline |
| A1 | H221 | repaired V7 control | 300 GP + recovery | 重建 V7 baseline |
| A2 | H375 | finite-radius smoke | 1/20 | 验证修复健康度 |

## Phase R — Restart

| ID | 输入 | Pipeline | Budget |
|---|---|---|---:|
| R1C | H221 | exact_joint 150 continuous | 150 |
| R1R | H221 | exact_joint 50×3 | 150 |
| R3C | H221/H252 | bridge 360 continuous | 360 |
| R3R | H221/H252 | bridge120×3 | 360 |
| R2-1 | H221 | V7×1 | 300 |
| R2-2 | H221 | V7×2 | 600 |
| R2-3 | H221 | V7×3 | 900 |

## Phase BB

| ID | 输入 | Variant | Budget |
|---|---|---|---:|
| BB0 | H221 | fixed-step matched control | 100 |
| BB1 | H221 | BB2 | 100 |
| BB2 | H221 | BB1 | 100 |
| BB3 | H221 | best BB | 300 |
| BB4 | H221 | V7 bridge→BB | V7 budget |

## Phase EC

| ID | 输入 | Variant | Budget |
|---|---|---|---:|
| EC0 | H221 | fixed epsilon | 100 |
| EC1 | H221 | conservative | 100 |
| EC2 | H221 | wide | 100 |
| EC3 | H221 | best EC | 300 |
| EC4 | H221 | V7 bridge→EC | V7 budget |

## Phase FUSION

只有 BB/EC 各自成功才运行。

---

# 41. 所有实验固定变量

同一实验组必须固定：

```text
case
input checkpoint SHA-256
canonical grid 512×512
target density 1.0
seed
threads
iteration budget
retention
exact evaluator implementation
Git commit
high-degree mode
preconditioner mode
lambda policy（除非它就是实验变量）
```

特别注意：

Commit A 修复后，不要把新实验直接与旧 commit 数值做强因果比较。

必须在同一 commit 下重新运行 control。

---

# 42. 结果记录

每个 run 仍默认只保留：

```text
params.json
experiment.md
trajectory.csv
```

不要默认保存 `.pl`。

对于 pipeline 中重复 V7 stage：

不需要中间 placement 文件，直接使用内存 `Database` 交接。

只有必须人工复核某个关键阶段时，才使用：

```text
--save-stage <module>
```

但重复同名 module 存在命名冲突风险，因此本轮最好不要依赖 `--save-stage` 做 restart 实验。

若后续需要保存，先修 `--save-stage` 支持 `stage_index:name`，不要覆盖同名 stage；此项不是本轮必须。

---

# 43. 需要新增的 trajectory telemetry

## 通用

```text
global_iteration
stage_iteration
exact_hpwl
overflow_percent
density_energy
max_density
lambda
full_objective_evaluations
local_density_queries
local_hpwl_queries
step_rms
step_max
fraction_coordinates_clipped
```

## restart

在 stage boundary `experiment.md` 记录：

```text
solver_state_reset = optimizer,lambda,...
restart_cycle
```

## BB

```text
bb_valid
bb_restart
alpha_bb1
alpha_bb2
alpha_lip
alpha_used
sTy
yTy
```

## epsilon continuation

```text
epsilon_wire
epsilon_density
epsilon_band
epsilon_restart
```

不要为了 telemetry 新建独立大文件。

---

# 44. 分析脚本

允许新增轻量脚本：

```text
scripts/analyze_v8_experiments.py
```

只读取三文件结果，不修改 experiment。

输出到：

```text
framework/results/analysis/<date>_v8_.../
```

该目录默认不入 Git。

建议生成：

```text
runs.csv
stage_boundaries.csv
per_iteration_metrics.csv
```

以及图：

```text
HPWL vs iteration
overflow vs iteration
HPWL vs overflow
step size vs iteration
lambda vs iteration
epsilon vs iteration
```

图不是完成任务的必要证据；CSV + exact summary 更重要。

---

# 45. 性能测量

Fix A2 会显著改变 finite-radius 成本。

必须记录修复前后的：

```text
wall_seconds / iteration
full exact evaluations
local density queries
```

至少做：

```text
probe_count = 4
probe_count = 32
```

若 `evaluate_move()` 版 32 probe 仍比旧 4 probe 更快或同量级，说明重构成功。

不允许用“看起来快了”代替测量。

---

# 46. 代码所有权和文件清单

预计修改/新增：

## framework

```text
framework/code/main.cpp
framework/code/microkernel.hpp
framework/code/microkernel.cpp

framework/kernel/include/epsilon_active/types.hpp
framework/kernel/include/epsilon_active/density.hpp
framework/kernel/include/epsilon_active/hpwl.hpp
framework/kernel/include/epsilon_active/nonlocal_descent.hpp

framework/kernel/src/bookshelf.cpp
framework/kernel/src/density.cpp
framework/kernel/src/hpwl.cpp
framework/kernel/src/nonlocal_descent.cpp
```

只改实际需要文件。

## existing modules

```text
modules/finite_radius_oracle_gp/code/module.cpp
modules/nonlocal_common.hpp
modules/exact_joint_gp/code/placer.cpp
modules/exact_joint_gp/code/module.cpp
```

若 high-degree/preconditioner 已完全在 kernel 层解决，则 module 只解析参数。

## new modules

```text
modules/spectral_bb_gp/...
modules/overflow_epsilon_continuation_gp/...
```

## tests

```text
tests/test_numeric_contracts.cpp
tests/test_nonlocal_oracle_contracts.cpp
tests/test_experiment_log_contracts.cpp       # 可选独立
tests/test_dreamplace_inspired_contracts.cpp
```

## pipelines / params

新增 V8 control、restart、BB、EC 配置。

---

# 47. CMake 要求

新增 module 必须：

```text
只增加必要 source
不创建第二套 library hierarchy
不复制 hpwl/density/optimizer 源码
```

CTest 新增：

```text
nsgp_experiment_log_contracts       # 若独立
nsgp_dreamplace_inspired_contracts
```

运行：

```powershell
cmake -S . -B build
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

---

# 48. 文档同步要求

Commit A 若改变：

- evaluator deterministic mode；
- StageStats；
- trajectory schema；
- audit CLI；
- preconditioner；
- high-degree policy；

至少同步：

```text
you should know about/00_START_HERE.md          # audit CLI 若变化
you should know about/01_PROJECT_REQUIREMENTS.md # 仅需补 deterministic canonical contract 时
you should know about/02_ARCHITECTURE_AND_CODE_MAP.md
you should know about/03_MODULE_CATALOG.md
you should know about/04_EXPERIMENT_RULES.md
you should know about/05_CURRENT_STATE.md
you should know about/07_FILE_INDEX.md
you should know about/CHANGELOG.md
```

Commit C/D 新 module 必须更新：

```text
02
03
05
07
CHANGELOG
```

本方案本身建议放入：

```text
plan/NONSMOOTH_GP_REPAIR_RESTART_DREAMPLACE_INTEGRATION_PLAN_V8.md
```

---

# 49. 本轮明确不做

以下不要顺手扩大范围：

```text
Linux/macOS SHA-256
完整 GPU/CUDA 化
PyTorch 重写
legalization
detailed placement
routing congestion
timing-driven net weighting
WA/LSE objective
DCT/Poisson 正式化
大型 min-cost-flow
八 case 全量 benchmark
大量 seed sweep
大规模超参数网格
```

只有 adaptec1 单 case 获得稳定正信号后，才扩展。

---

# 50. 最终执行顺序

本地模型严格按下列顺序做。

## Step 1 — 只读审计

```text
read protocols
git status
read current source
确认哪些旧 plan 已经实现
```

## Step 2 — Commit A correctness repair

按 A1→A10 修。

## Step 3 — Build + CTest

任何测试失败先修，不跑大实验。

## Step 4 — smoke

```text
finite radius synthetic
finite radius 1-round
H221 1-round
V7 pipeline structure smoke
```

## Step 5 — repaired controls

建立：

```text
V8_REPAIRED_H221_CONTROL
V8_REPAIRED_V7_CONTROL
```

## Step 6 — restart R1

先做等预算：

```text
150 continuous
vs
50×3
```

这是 restart 是否值得继续的首要判断。

## Step 7 — restart R2/R3

若 R1 正向，跑：

```text
bridge120×3 vs bridge360
V7×1/×2/×3
```

若 R1 明显负向，V7×2 仍可做一次实用验证，但不要直接跑 ×3 浪费预算。

## Step 8 — spectral BB module

先 unit tests，再 100-step H221 screen。

## Step 9 — epsilon continuation module

先 unit tests，再 100-step H221 screen。

## Step 10 — 300-step / V7 replacement

只有 100-step 晋级者运行。

## Step 11 — fusion

只有两个 module 各自独立为正才运行。

## Step 12 — docs + Git

每个逻辑提交同步 handoff 文档。

---

# 51. 最终报告模板

任务执行完成后，本地模型最终回答必须至少给出：

```text
Branch:
HEAD commits:

Protocol:
- exact HPWL unchanged
- exact density unchanged
- no smooth surrogate introduced

Fixes:
- finite-radius descent gate
- incremental probe
- deterministic audit
- HPWL context
- query accounting
- trajectory schema/global iteration
- audit clamp
- weighted preconditioner
- high-degree exact extrema
- tests

Tests:
- cmake build
- ctest ...
- retention ...

Repaired controls:
H221: ...
V7: ...

Restart:
150 continuous: ...
3×50 restart: ...
bridge continuous/restart: ...
V7×1: ...
V7×2: ...
V7×3: ...

spectral_bb_gp:
best result: ...
vs matched control: ...

overflow_epsilon_continuation_gp:
best result: ...
vs matched control: ...

Fusion:
run / not run and why

Performance:
finite-radius old/new runtime and query counts

Remaining limitations:
- adaptec1 only
- one checkpoint
- no legalization
- cross-platform SHA still pending
```

失败结果同样必须报告，禁止只写成功项。

---

# 52. 研究判断原则

这轮最重要的不是“必须打到某个数字”，而是把三个问题做成可证伪实验：

## Q1 — 修复后，当前负结果有多少是实现问题？

finite-radius / deterministic audit / accounting 修复后重新测。

## Q2 — restart 是否有独立价值？

以：

```text
150 continuous vs 3×50
```

这个等预算对照回答，不用 V7×3 的更大预算替代因果证据。

## Q3 — DREAMPlace 的哪些机制能在 exact non-smooth 框架中保留价值？

分别检验：

```text
spectral step-size adaptation
overflow-driven active-set continuation
```

而不是把 DREAMPlace 的平滑目标搬进来。

如果最终结论是：

```text
restart neutral
BB negative
epsilon continuation positive
```

也是完整且有价值的研究结论。

---

# 53. 最终红线

本地模型在整个执行过程中必须一直满足：

\[
\boxed{
\text{exact objective}
+
\text{non-smooth directions}
+
\text{reproducible canonical audit}
}
\]

任何时候如果一个新实现需要：

```text
WA/LSE
electric potential objective
DCT/Poisson objective
smooth envelope replacing exact metrics
```

才能工作，就停止该实现，不得为了结果越过协议。

DREAMPlace 在本轮只是提供：

\[
\boxed{
\text{step-size adaptation}
+
\text{continuation control idea}
}
\]

不是提供新的最终目标函数。

---

# 54. 验收清单

## Task 1

- [ ] finite-radius 无下降方向返回 0
- [ ] finite-radius 使用 `evaluate_move`
- [ ] actual displacement denominator
- [ ] query accounting
- [ ] nonlocal exact HPWL context
- [ ] deterministic canonical density
- [ ] deterministic feasibility decisions
- [ ] CSV schema 修复
- [ ] global iteration 修复
- [ ] concrete nonlocal module telemetry name
- [ ] canonical audit clamp
- [ ] optional pipeline density config for audit
- [ ] weighted pin sum preconditioner
- [ ] high-degree exact-extrema mode
- [ ] expanded tests
- [ ] cross-platform SHA 未动

## Task 2

- [ ] repaired V7 baseline
- [ ] 150 continuous
- [ ] 3×50 restart
- [ ] bridge360 continuous
- [ ] bridge120×3
- [ ] V7×1
- [ ] V7×2
- [ ] V7×3 only if justified
- [ ] restart marginal gain analysis

## Task 3

- [ ] `spectral_bb_gp` module
- [ ] BB contract tests
- [ ] BB 100-step screen
- [ ] BB promoted experiment if positive
- [ ] `overflow_epsilon_continuation_gp` module
- [ ] epsilon schedule contract tests
- [ ] EC 100-step screen
- [ ] EC promoted experiment if positive
- [ ] fusion only after independent positive evidence

## Project handoff

- [ ] relevant `you should know about/*.md` updated
- [ ] `CHANGELOG.md` updated in same commits
- [ ] `git diff --check`
- [ ] no dataset/binary/ordinary experiment directories committed
- [ ] push `origin/nonsmooth-gp-v1`
