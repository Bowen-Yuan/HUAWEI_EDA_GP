# nonsmooth-gp-v1：低产物实验运行 + 全局视野策略 + 优化器数值实验编码方案

> 目标：直接交给本地 AI，在 `Bowen-Yuan/HUAWEI_EDA_GP` 的 `nonsmooth-gp-v1` 分支上继续开发。
>
> 本轮有三个目标：
>
> 1. **先以当前分支真实代码为准做只读扫描**，不要凭旧方案重建项目。
> 2. 将实验产物策略改为 **metrics-only**：实验完成后默认不保留任何 `.pl`、snapshot、候选 dump、缓存或其他中间产物，只保留简易参数、实验说明和收敛指标。
> 3. 在该基础上实现“全局视野”非光滑策略与 optimizer / step-size 数值实验，默认从一个已存在的 `adaptec1`、exact overflow≈7% 的外部 checkpoint 开始。

---

# 0. 最高优先级：本地 AI 先扫描当前 `nonsmooth-gp-v1`

远程规划文件中的文件名只能当参考。真正修改前，本地 AI 必须以当前工作树为准。

执行：

```bash
git checkout nonsmooth-gp-v1
git branch --show-current
git rev-parse HEAD
git log -n 15 --date=iso --oneline
git status --porcelain
```

然后定位当前代码：

```bash
rg -n "ModuleResult|RunContext|Pipeline|pipeline|module_registry|selected\.pl|snapshot|write_bookshelf_placement|write.*placement" .
rg -n "ExactHpwl|ExactOverlap|overflow|target_density|density_energy|evaluate_move|commit_move" .
rg -n "Optimizer|Adam|AMSGrad|AdaGrad|HeavyBall|SGD|learning_rate|maximum_delta|backtrack|trust" .
rg -n "results/|result_dir|output_dir|checkpoint|manifest|trajectory" .
```

必须先产出一个**不进入实验结果目录**的临时控制台摘要，确认：

- 当前 pipeline runner 的真实文件；
- module API 的真实文件；
- placement 写盘入口；
- snapshot 写盘入口；
- experiment/run directory 创建逻辑；
- 当前 optimizer 实现；
- 当前 `exact_joint_gp` 或等价连续 GP 模块；
- 当前 exact auditor/evaluator；
- 当前 `adaptec1` canonical evaluator grid、`target_density`；
- 当前可用的 overflow≈7% 外部 checkpoint。

如果当前分支已有同等职责的类/函数，**直接修改现有实现，不新造平行框架**。

---

# 1. 本轮不改变的数值契约

所有实验继续使用统一权威 evaluator：

```text
HPWL = exact weighted pin-offset max-min HPWL

rho_b =
    exact physical rectangle/bin overlap occupancy
    / bin_area

Overflow =
    sum_b bin_area * max(rho_b - rho_t, 0)
    / total_movable_area
```

正式实验的优化方向可以变，但权威评价不能变。

禁止：

- WA / LSE 替代最终 HPWL；
- smooth density surrogate 替代最终 exact rectangle/bin overlap；
- Moreau envelope 替代最终目标；
- DCT/Poisson 进入正式 nonsmooth pipeline。

允许：

- epsilon-active；
- 多 active-face 聚合；
- temporal bundle / direction memory；
- exact-audited relocation；
- exact backtracking / trust radius；
- coarse capacity map 作为**候选生成**；
- 临时 overflow corridor，但最终结果仍由 exact final cap 选择。

---

# 2. 第一项代码修改：实验结果改成 metrics-only

## 2.1 默认永久保留内容

每次 experiment 完成后，默认只允许存在：

```text
<experiment_dir>/
├─ params.json
├─ experiment.md
└─ trajectory.csv
```

其中：

### `params.json`

保存完整但简洁的实验配置：

- case；
- external input checkpoint path；
- input SHA256；
- module chain；
- 每个 module 参数；
- optimizer；
- step policy；
- acceptance policy；
- threads；
- seed；
- exact evaluator 参数；
- retention policy；
- git branch / commit。

### `experiment.md`

保存面向人的简易说明：

- 实验目的；
- 输入 checkpoint 来源；
- 模块调用顺序；
- 起始 exact HPWL / overflow；
- best-feasible exact HPWL / overflow；
- last exact HPWL / overflow；
- 总 runtime；
- 每模块 runtime；
- accepted/rejected step 数；
- 是否发生 NaN/异常；
- 是否显式要求保留某个 stage placement；
- 简短观察。

### `trajectory.csv`

用于画收敛曲线，不保存布局。

推荐固定字段：

```text
global_iter,
stage_index,
module,
module_iter,
hpwl,
overflow,
density_energy,
max_density,
best_feasible_hpwl,
best_feasible_overflow,
lambda,
step_size,
trust_radius,
accepted,
wall_seconds,
direction_cos_prev,
bundle_norm_ratio
```

某个 module 不使用的字段写空值。

**不要默认保存单独 summary.json / manifest.json / final.pl。**
最终标量总结直接写进 `experiment.md`，机器读取时也可由 `trajectory.csv` 最后一行 / best-feasible 列恢复。

---

# 3. 默认禁止持久化的内容

以下默认一律不能出现在 experiment 完成后的目录：

```text
*.pl
snapshots/
checkpoint/
candidate_dump/
gradient_dump/
density_map_dump/
price_field_dump/
bundle_history_dump/
optimizer_state/
occupancy_cache/
debug/
*.bin
*.npy
*.raw
```

也不要复制：

- 原始 benchmark；
- external input checkpoint；
- 每个 stage 的 input/output placement；
- trajectory 对应的每轮坐标。

---

# 4. `RetentionPolicy`：把“不保存”做成框架默认行为

在当前 `RunContext` 或等价对象中加入一个很小的配置结构，不要建立大型 artifact manager。

接口草案：

```cpp
struct RetentionPolicy {
    bool keep_final_placement = false;
    bool keep_snapshots = false;
    bool keep_debug_artifacts = false;

    // 默认空。只有用户明确指定 stage id 时才保留该 stage 的 audited .pl
    std::vector<std::string> keep_stage_placements;

    bool cleanup_on_success = true;
    bool cleanup_on_failure = true;
};
```

默认配置：

```json
{
  "retention": {
    "mode": "metrics_only",
    "keep_final_placement": false,
    "keep_stage_placements": [],
    "keep_snapshots": false,
    "keep_debug_artifacts": false,
    "cleanup_on_success": true,
    "cleanup_on_failure": true
  }
}
```

如果用户明确要求保存某个实验中某个模块：

```json
{
  "retention": {
    "mode": "metrics_only",
    "keep_final_placement": false,
    "keep_stage_placements": ["capacity_slack"],
    "keep_snapshots": false
  }
}
```

此时只允许增加：

```text
<experiment_dir>/saved/capacity_slack.pl
```

不要顺便保留前后其他 stage。

建议同时支持 CLI：

```text
--save-stage <stage_id>
```

该选项可以重复出现。

---

# 5. 临时 PL 的生命周期

## 5.1 第一选择：同一进程直接传 `Layout`

若当前 pipeline 已经是：

```cpp
ModuleResult run(
    const Problem&,
    const Layout&,
    ...);
```

则 stage 间直接：

```cpp
Layout current = input;

for (stage : pipeline) {
    ModuleResult result = run_module(problem, current, ...);
    current = std::move(result.layout);
}
```

这种情况下 **完全不需要 `.pl` 作为 stage handoff**。

`.pl` 只保留为：

- 外部输入格式；
- 用户明确要求保存时的输出格式；
- 某个旧 module 仍必须通过文件启动时的临时兼容格式。

## 5.2 若某个现有模块仍依赖 `.pl` 路径

建立一个**临时 workspace**：

```text
<system-temp>/nonsmooth-gp/<run_id>/
```

运行中可出现：

```text
stage_000_input.pl
stage_000_output.pl
stage_001_output.pl
```

但是：

1. 这些文件绝不放进正式 experiment_dir；
2. 下一个 stage 读取成功后，立即删除已无用的上一个 `.pl`；
3. experiment 完成后删除整个 temp workspace；
4. experiment 失败也默认删除；
5. 外部用户提供的输入 `.pl` 绝不能删除。

必须实现路径安全检查：

```cpp
bool is_under_temp_workspace(path);
```

只有 `true` 才允许自动 `remove()`。

---

# 6. `ExperimentWorkspace`：最小实现

建议在现有 runner 附近增加：

```cpp
class ExperimentWorkspace {
public:
    ExperimentWorkspace(std::filesystem::path temp_root,
                        std::string run_id);
    ~ExperimentWorkspace();

    std::filesystem::path temp_path(std::string_view name) const;
    void cleanup_noexcept();

private:
    std::filesystem::path root_;
    bool cleanup_enabled_ = true;
};
```

规则：

- constructor 创建唯一 temp dir；
- destructor 尽力递归清理；
- success/failure 都清理；
- 所有 module 临时 artifact 必须在该 root 下；
- 正式结果目录只允许写三个 metrics-only 文件；
- 若用户显式 `--save-stage`，runner 在 fresh audit 后把指定 stage 写入 `saved/`，它不是 temp 文件。

---

# 7. 必须清查所有直接写盘点

本地 AI 修改时必须执行：

```bash
rg -n "write_bookshelf_placement|ofstream|snapshot|selected\.pl|global\.pl|last\.pl|best\.pl" framework modules
```

处理原则：

### A. 权威 framework writer

保留 `write_bookshelf_placement()` 本身，它是 I/O primitive。

### B. module 内部无条件写 snapshot

删除或改为：

```cpp
if (context.retention.keep_snapshots) {
    ...
}
```

但默认 false。

### C. module 内部无条件写 `selected.pl`

删除。
module 应返回 `Layout`，不应该决定永久文件。

### D. runner 无条件写 final `.pl`

删除。
只有：

```cpp
retention.keep_final_placement == true
```

时才写。

### E. 中间 module 需要路径 handoff

写到 temp workspace，不写正式结果目录。

---

# 8. 内存中的中间大对象也要及时释放

每个 module 结束后，不允许 `ModuleResult` 带回：

- full gradient history；
- occupancy scratch；
- candidate list；
- density map；
- bundle work matrix；
- optimizer scratch arrays。

`ModuleResult` 只保留：

```cpp
struct ModuleResult {
    Layout layout;
    ModuleStats stats;   // 标量/小数组
};
```

模块内部大对象必须限制在 `run()` 局部 scope。

对明显大的 `std::vector`，stage 边界若确实需要立即归还内存，可使用：

```cpp
std::vector<Real>().swap(buffer);
```

但不要在每个 iteration 做 `shrink_to_fit()`。

Temporal bundle 属于 module 内状态，module 结束后必须析构。

---

# 9. 结果生命周期测试

增加小型测试，不需要复杂测试框架。

## T1：metrics-only success

运行一个 1～2 iter smoke pipeline。

结束后 assert：

```text
params.json exists
experiment.md exists
trajectory.csv exists
```

并 assert：

```text
find experiment_dir -name "*.pl" -> empty
find experiment_dir -name "snapshots" -> empty
```

## T2：failure cleanup

人为让第二 stage 返回 failure。

结束后：

- experiment.md 记录 failure；
- trajectory.csv 保留已有标量；
- temp workspace 已删除；
- 没有 `.pl` 泄漏。

## T3：explicit stage save

设置：

```text
--save-stage stage_b
```

结束后只允许：

```text
saved/stage_b.pl
```

## T4：external input safety

external input checkpoint 在 run 前后 SHA256 相同且文件仍存在。

## T5：handoff correctness

同一个两 stage pipeline：

- 内存 `Layout` 直接交接；
- `.pl` temp round-trip 交接；

fresh audit 结果应在数值容差内一致。

---

# 10. 实验入口：固定一个已有 7% checkpoint

本轮不要重复前面的 spreading 链。

从一个**已经存在**的 `adaptec1` checkpoint 开始：

```text
exact overflow ≈ 0.07
```

该文件是 external input，runner：

1. 只读它；
2. 计算 SHA256；
3. fresh audit；
4. 在 `experiment.md` 写清来源；
5. 永远不复制到 experiment_dir；
6. 永远不自动删除它。

配置示例：

```json
{
  "case": "adaptec1",
  "input_placement": "D:/.../known_ov07_checkpoint.pl",
  "input_policy": "external_readonly",
  "optimizer_state": "reset"
}
```

公平比较时所有实验必须使用同一 checkpoint hash。

---

# 11. 全局策略 A：`global_view_gp`

## 11.1 目标

解决：

> 当前 subgradient / epsilon-active 仍然主要描述当前 active face；希望同时看到多个 plausible active sets，并利用过去若干 accepted iterations 的方向关系。

推荐不是复制一份新的完整 GP，而是尽可能复用当前连续 GP 主循环：

```text
exact density direction
+ exact HPWL direction
+ lambda controller
+ optimizer
+ exact acceptance
```

只替换/扩展 **HPWL direction provider** 和 **temporal direction composition**。

若当前代码已有等价 strategy enum，直接扩展 enum。

---

# 12. Active Ensemble：多 epsilon active-face 视野

同一个 Layout 上计算：

```text
epsilon scales = [0.0, 0.25, 0.5, 1.0] × bin_size
```

每个 epsilon：

- exact HPWL 数值不变；
- 只改变搜索方向中包含的 active/near-active pins。

得到：

```text
g_w[0], g_w[1], ..., g_w[K-1]
```

每个先做一致尺度归一：

```text
g_k <- g_k / max(rms(g_k), eps)
```

支持两种 aggregation。

## Control：weighted mean

```text
g = sum_k w_k g_k
```

## 推荐：minimum-norm convex combination

求：

```text
min_alpha || sum_k alpha_k g_k ||^2
s.t. alpha_k >= 0
     sum alpha_k = 1
```

K 只取 3～4。

Gram matrix：

```cpp
Q[i][j] = dot(g[i], g[j]);
```

解：

```text
min alpha^T Q alpha
on probability simplex
```

不用外部 QP 库。

实现一个小型 projected-gradient + simplex projection：

```cpp
std::vector<Real> solve_min_norm_simplex(
    const SmallMatrix& gram,
    int iterations);
```

默认 20 次小迭代即可。

配置：

```json
{
  "active_ensemble": {
    "enabled": true,
    "epsilon_unit": "bin",
    "epsilon_scales": [0.0, 0.25, 0.5, 1.0],
    "aggregation": "min_norm_convex",
    "normalize_each": "rms",
    "qp_iterations": 20
  }
}
```

第一轮如果 runtime 过高，先用：

```text
[0.0, 0.25, 0.5]
```

---

# 13. Temporal Bundle：跨迭代“时间视野”

## 13.1 保存什么

不要保存每个 epsilon 的历史。

每个 accepted iteration 只保存最终组合后的 joint search gradient：

```text
g_t = normalized(
    g_hpwl_ensemble
    + lambda * g_density
)
```

保存最近：

```text
K = 4
```

个 accepted gradients。

为了省内存，history 用 `float`：

```cpp
std::deque<std::vector<float>> history;
```

当前计算仍用项目的 `Real`。

---

# 14. Temporal bundle 聚合

从：

```text
g_t, g_{t-1}, ..., g_{t-K+1}
```

构造 K×K Gram matrix，解同样的 simplex min-norm：

```text
g_bundle = sum alpha_i g_i
```

最终：

```text
g_search =
    (1 - bundle_mix) * g_current
    + bundle_mix * g_bundle
```

首轮：

```text
bundle_mix = 0.5 或 0.75
K = 4
```

不要一开始 `mix=1`。

配置：

```json
{
  "temporal_bundle": {
    "enabled": true,
    "history": 4,
    "storage": "float32",
    "mix": 0.75,
    "qp_iterations": 20,
    "max_history_mb": 256,
    "reset_lambda_relative_change": 0.20,
    "reset_after_rejections": 4
  }
}
```

---

# 15. Bundle reset

以下任一发生时清空 history：

- lambda 相对变化 > 20%；
- density grid 改变；
- target density 改变；
- 连续 4 次 exact reject；
- external checkpoint 启动；
- NaN / Inf；
- optimizer 被 reset；
- module restart。

trajectory 额外记录：

```text
direction_cos_prev
bundle_norm_ratio =
    ||g_bundle|| / max(||g_current||, eps)
```

这样以后可以判断“是否真的减少 face oscillation”。

不要保存 bundle gradient 本身。

---

# 16. `global_view_gp` 主循环草案

```cpp
Layout x = input;
Layout best_feasible = input;
Metrics best_metrics = fresh_audit(input);

for (int iter = 0; iter < cfg.iterations; ++iter) {

    // 1. exact density metric + nonsmooth direction
    DensityEval d = density.evaluate_with_direction(x, ...);

    // 2. exact HPWL + multi-epsilon direction
    ActiveEnsembleResult w = active_ensemble.compute(x, ...);

    // 3. lambda
    Real lambda = lambda_controller.update(...);

    // 4. current joint direction
    Vector g_current = combine_and_normalize(
        w.direction, d.direction, lambda);

    // 5. temporal view
    Vector g_search = temporal_bundle.compose(
        g_current, lambda);

    // 6. optimizer proposal
    Delta proposal = optimizer.compute_delta(
        g_search, step_policy.current_scale(), ...);

    // 7. exact audited acceptance
    StepResult step = step_policy.try_step(
        x, proposal, exact_auditor, acceptance_cfg);

    if (step.accepted) {
        x = std::move(step.layout);
        temporal_bundle.on_accept(g_current, lambda);
    } else {
        temporal_bundle.on_reject();
    }

    // 8. maintain best final-feasible in MEMORY ONLY
    Metrics m = step.metrics_or_fresh_audit();

    if (m.overflow <= cfg.final_overflow_cap &&
        m.hpwl < best_metrics.hpwl) {
        best_feasible = x;
        best_metrics = m;
    }

    trajectory_writer.write_scalar_row(...);
}

return ModuleResult{
    .layout = std::move(best_feasible),
    .stats = ...
};
```

注意：

```text
best_feasible
```

只存在内存，module 结束后交给下一个 stage。

默认绝不写 `selected.pl`。

---

# 17. 全局策略 B：`global_capacity_transport`

## 17.1 假设

从正好约 7% overflow 的布局开始时：

```text
HPWL recovery move
```

很容易因为稍微增加局部密度而被 7% cap 卡住。

因此先试：

```text
7%
→ exact capacity transport 降到 6.5% / 6.0%
→ 再用连续方法在 7% cap 下恢复 HPWL
```

这可以区分：

- 是局部方向不好；
- 还是 7% 边界根本没有可移动 slack。

---

# 18. capacity map

使用较粗 grid：

```text
64×64 或 128×128
```

但 occupancy 仍由 exact rectangle/bin overlap 构造。

对 coarse bin：

```text
excess_b =
    max(occupancy_b - rho_t * area_b, 0)

capacity_b =
    max(rho_t * area_b - occupancy_b, 0)
```

coarse grid 只用于候选生成。

正式评价仍用 canonical exact evaluator。

---

# 19. source → destination 搜索

第一版不要实现复杂 optimal transport。

对 source bin 按 excess 降序：

```text
radius = 1, 2, 4, 8, 16...
```

寻找有 capacity 的 destination。

destination score：

```text
cost =
    distance_weight * normalized_distance
    + hpwl_weight * positive_hpwl_delta_estimate
```

每个 source 只尝试有限数量 cell。

cell score：

```text
score =
    relieved_excess_area
    / (1 + positive_hpwl_delta / hpwl_scale)
```

优先：

- 释放 excess 多；
- HPWL 损失低；
- destination 容量充分；
- 位移相对短。

第一版：

```text
move_macros = false
```

---

# 20. exact-audited move

candidate anchors：

```text
nearest feasible point in destination bin
bin center
HPWL-guided point clamped to destination
```

每个 candidate 必须通过已有 incremental/fresh exact oracle：

```text
exact overflow delta
exact HPWL delta
```

接受：

```text
new_overflow < current_overflow
AND hpwl_damage <= configured budget
```

达到目标 overflow 即停止。

默认参数起点：

```json
{
  "coarse_bins": [64, 64],
  "target_overflow": 0.060,
  "max_sweeps": 4,
  "move_macros": false,
  "search_radii": [1, 2, 4, 8, 16],
  "candidate_destinations": 8,
  "candidate_cells_per_source": 64,
  "exact_acceptance": true
}
```

不要保存：

- source map；
- destination map；
- candidate list；
- 每轮 `.pl`。

只在 trajectory 写：

```text
sweep
hpwl
overflow
accepted_moves
wall_seconds
```

---

# 21. optimizer 数值实验：原则

必须把：

```text
direction
optimizer
step policy
acceptance
```

拆开。

第一轮 optimizer 对比时：

```text
direction 完全固定
step policy 完全固定
acceptance 完全固定
input checkpoint 完全固定
```

只换 optimizer。

所有从 7% checkpoint 开始的实验：

```text
optimizer_state = reset
```

不要继承原始 checkpoint 之前的 Adam moments。

---

# 22. 第一轮 optimizer

优先复用当前分支已经实现的：

```text
Adam
AMSGrad
AdaGrad
HeavyBall
SGD
```

若缺少再补。

新增两个小型 experiment optimizer：

## Normalized SGD

```text
d = -g / max(rms(g), eps)
delta = step_size * d
```

目的：

> 判断 Adam 的优势是 adaptive-moment 本身，还是主要来自稳定 gradient scale。

## Dual Averaging

```text
s_t = s_{t-1} + g_t
alpha_t = alpha0 / sqrt(t + 1)
x_t = Project(
    x0 - alpha_t * s_t / scale_t
)
```

必须有：

- region clamp；
- max displacement；
- exact acceptance。

这里只作为非凸数值实验，不声明理论保证。

---

# 23. optimizer API

如果当前已有：

```cpp
class Optimizer {
public:
    virtual void reset(size_t dims) = 0;
    virtual void compute_delta(...) = 0;
};
```

继续沿用。

不要新建第二套 optimizer hierarchy。

必要时补：

```cpp
enum class OptimizerKind {
    Adam,
    AMSGrad,
    AdaGrad,
    HeavyBall,
    SGD,
    NormalizedSGD,
    DualAveraging
};
```

如果当前项目已经用 factory/enum，则扩展现有 enum。

---

# 24. step policy 实验

很多非光滑实验差异可能主要来自“走多远”，所以 optimizer 和 step 必须单独研究。

第一轮保留：

```text
current_control
1/sqrt(t)
exact_backtracking
exact_trust_radius
```

第二轮再试：

```text
BB spectral step
```

---

# 25. exact backtracking

optimizer 产生原 proposal：

```text
scale = 1
```

依次：

```text
1
1/2
1/4
1/8
...
```

每个 candidate fresh/incremental exact audit。

配置：

```json
{
  "type": "exact_backtracking",
  "max_backtracks": 8,
  "shrink": 0.5
}
```

---

# 26. exact trust radius

维护最大 cell displacement：

```text
R_t
```

proposal 先 clip：

```text
||delta_i|| <= R_t
```

规则：

```text
reject        -> R *= 0.5
连续3次有效accept -> R *= 1.25
```

限制：

```text
R_min <= R <= R_max
```

建议用 bin 为单位：

```json
{
  "type": "exact_trust_radius",
  "initial_radius_bins": 0.25,
  "min_radius_bins": 0.01,
  "max_radius_bins": 2.0,
  "grow": 1.25,
  "shrink": 0.5,
  "grow_after_successes": 3
}
```

这是本轮非常值得优先试的机制。

---

# 27. BB spectral step：第二轮

只给 SGD / normalized SGD。

```text
s = x_t - x_{t-1}
y = g_t - g_{t-1}

BB1 = (s^T s) / (s^T y)
BB2 = (s^T y) / (y^T y)
```

必须 guard：

```text
if sTy <= curvature_eps:
    fallback

alpha = clamp(alpha, alpha_min, alpha_max)
```

nonsmooth 下要配 exact/nonmonotone guard，不能裸跑。

---

# 28. acceptance 实验

## Control：strict 7% cap

```text
candidate overflow <= 0.07 + tol
AND candidate HPWL < current HPWL
```

用于最干净比较。

## 第二阶段：shrinking corridor

只有 strict-cap 实验出现“大量 reject / 基本无法走动”时才试：

```text
8% -> 7%
```

例如：

```text
cap(t) =
    0.07
    + 0.01 * (1 - t/T)^1.5
```

但是 best-feasible 永远只接受：

```text
overflow <= 0.07
```

如果最后没有找到新的 7% feasible point：

```text
return original input layout
```

而不是返回 7.5% last iterate。

---

# 29. 实验顺序

不要全排列。

## Phase 0：锁定输入

只做一次：

```text
same adaptec1 external checkpoint
same SHA256
fresh audit
same threads
same seed
same build
```

在每个 experiment.md 引用该 input hash。

## Phase 1：optimizer-only

固定当前 control direction + current control step + strict cap。

```text
O0 Adam
O1 AMSGrad
O2 AdaGrad
O3 HeavyBall
O4 NormalizedSGD
O5 DualAveraging
```

先：

```text
50 iterations
```

若成本可接受，再 top 3：

```text
100 / 200 iterations
```

比较：

```text
best-feasible HPWL
best-feasible overflow
last HPWL
last overflow
runtime
accepted/rejected
```

## Phase 2：step policy

取 Phase 1 最好的 1～2 个 optimizer。

例如：

```text
Adam + control
Adam + exact backtracking
Adam + exact trust radius

best_non_adam + control
best_non_adam + exact backtracking
best_non_adam + exact trust radius
```

## Phase 3：global view

固定 Phase 2 最优 optimizer+step。

```text
G0 current epsilon-active
G1 ActiveEnsemble [0, .25, .5]
G2 ActiveEnsemble [0, .25, .5, 1.0]
G3 TemporalBundle K=4
G4 TemporalBundle K=8
G5 ActiveEnsemble + TemporalBundle K=4
```

先 strict cap。

如果大量 reject，再对 top 2 试 corridor。

## Phase 4：capacity slack

```text
C0 7% -> best continuous method

C1 7%
   -> capacity target 6.5%
   -> best continuous method, final cap 7%

C2 7%
   -> capacity target 6.0%
   -> best continuous method, final cap 7%

C3 7%
   -> capacity target 5.5%
   -> best continuous method, final cap 7%
```

---

# 30. 每个实验的持久化结果示例

```text
framework/results/experiments/
└─ 20260907_adaptec1_G3_bundle_k4/
   ├─ params.json
   ├─ experiment.md
   └─ trajectory.csv
```

不要有：

```text
selected.pl
last.pl
best.pl
snapshots/
stage_*.pl
```

运行中如果不得不用 temp PL：

```text
%TEMP%/nonsmooth-gp/20260907_adaptec1_G3_bundle_k4/
```

实验结束后该目录必须不存在。

---

# 31. `experiment.md` 模板

```markdown
# Experiment: G3 TemporalBundle K=4

## Purpose
Compare temporal bundle against current epsilon-active control from the same 7% checkpoint.

## Provenance
- branch: nonsmooth-gp-v1
- commit: <hash>
- case: adaptec1
- input checkpoint: <external path>
- input SHA256: <hash>
- threads: <n>
- seed: <seed>

## Module chain
1. global_view_gp

## Key parameters
- optimizer: Adam
- step policy: exact_trust_radius
- active ensemble: false
- temporal bundle: K=4, mix=0.75
- final overflow cap: 0.07

## Exact metrics
- start HPWL:
- start overflow:
- best-feasible HPWL:
- best-feasible overflow:
- last HPWL:
- last overflow:

## Runtime
- total:
- module global_view_gp:

## Search statistics
- accepted:
- rejected:
- bundle resets:

## Retention
- mode: metrics_only
- saved stage placement: none

## Observation
<2~5 lines>
```

---

# 32. `trajectory.csv` 记录频率

不要为了省文件而只记录最后一点。

CSV 很小，建议：

```text
每个 accepted/rejected optimization iteration 1 row
```

如果某 module 一次 sweep 内 candidate 很多，不要每 candidate 写一行。

例如 capacity transport：

```text
每 sweep 或每固定 N 个 committed moves 写一行
```

保证：

- 能画 HPWL vs iteration；
- 能画 overflow vs iteration；
- 能画 HPWL vs wall time；
- 能比较 trust radius；
- 能观察 face oscillation。

---

# 33. 代码层面的最小修改顺序

本地 AI 按以下顺序提交/验证。

## Commit A：只改 retention，不改算法

1. 扫描所有 placement/snapshot 写盘点。
2. 加 `RetentionPolicy`。
3. 加 temp `ExperimentWorkspace`。
4. runner 默认 metrics-only。
5. stage 之间优先传内存 Layout。
6. temp `.pl` 自动清理。
7. experiment_dir 只剩 3 类文件。
8. 加 artifact lifecycle tests。
9. 跑 adaptec1 极短 smoke，确认 exact metrics 与改动前相同。

**这一 commit 禁止改 optimizer / direction。**

## Commit B：optimizer lab

1. 抽出/确认 optimizer kind 配置。
2. 补 NormalizedSGD。
3. 补 DualAveraging。
4. step policy 独立配置。
5. exact backtracking。
6. exact trust radius。
7. trajectory 增加 step/trust 字段。
8. 跑 Phase 1/2 quick experiments。

## Commit C：Active Ensemble

1. 复用 exact HPWL oracle。
2. 加 multi-epsilon direction。
3. 加 small simplex QP。
4. 不改变 exact HPWL metric。
5. 跑 G0/G1/G2。

## Commit D：Temporal Bundle

1. history K=4 float32。
2. min-norm convex aggregation。
3. reset rules。
4. trajectory diagnostics。
5. 跑 G3/G4/G5。

## Commit E：Global Capacity Transport

1. coarse exact occupancy。
2. source/destination list。
3. candidate anchor。
4. exact move audit。
5. 不保存 map/candidate dump。
6. 跑 C0～C3。

---

# 34. 每次 commit 的 Definition of Done

## Retention DoD

- 默认实验完成后没有 `.pl`；
- 没有 snapshots；
- 没有 candidate/gradient/cache dump；
- temp workspace 被删除；
- external input 不被修改；
- 只有 `params.json / experiment.md / trajectory.csv`；
- 显式 `--save-stage X` 时只保存 X。

## Global-view DoD

- exact HPWL 返回值与 control 一致；
- exact overflow evaluator 不变；
- multi-epsilon 只改变 direction；
- bundle memory 有上限；
- module 结束后 bundle memory 释放；
- no PL/snapshot retained。

## Optimizer DoD

- 同一个 input SHA256；
- optimizer state reset；
- 被研究项之外参数固定；
- trajectory 能比较 convergence；
- runtime 记录；
- no PL retained。

## Capacity DoD

- coarse map 不作为最终 evaluator；
- 每个 accepted move exact-audited；
- 最终用 canonical exact evaluator；
- no map dump；
- no candidate dump；
- no PL retained。

---

# 35. 首轮推荐的最小实验集合

如果只想尽快判断方向，不要跑几十组。

先跑：

```text
E0 Adam control
E1 AMSGrad control
E2 NormalizedSGD control
E3 Adam + exact trust radius
E4 current direction + TemporalBundle K=4
E5 ActiveEnsemble [0,.25,.5] + TemporalBundle K=4
E6 capacity 6.0% -> E5-like recovery to 7%
```

全部：

```text
case = adaptec1
same input checkpoint
same hash
same threads
same seed
quick iterations = 50
metrics_only
```

看：

```text
best-feasible HPWL
overflow
runtime
accepted ratio
HPWL-vs-time
overflow-vs-time
direction_cos_prev
```

只有明显有潜力的方案再跑 100/200 iterations。

---

# 36. 研究解释优先级

完成首轮后按下面的问题解释结果，而不是只挑最低 HPWL：

1. **optimizer sensitivity 大吗？**
   - 如果 Adam / AMSGrad / normalized SGD 差很大，说明 scale/momentum 本身重要。

2. **trust radius 是否显著减少 reject？**
   - 若是，问题可能主要在 nonsmooth step control。

3. **TemporalBundle 是否提高 direction cosine 稳定性？**
   - 若是但 HPWL 没改善，说明稳定≠有用，需要调整 mix/acceptance。

4. **ActiveEnsemble 是否值得额外 HPWL direction 计算成本？**
   - 要用 HPWL-vs-wall-time 比，而不只看 iteration。

5. **先造 6% slack 后是否能换回更好 HPWL？**
   - 如果是，说明 7% boundary geometry 是当前局部 recovery 的主要限制之一。

---

# 37. 明确不做的事情

本轮不要：

- 跑完整八 case；
- baseline 比较；
- 大规模参数网格；
- 同时常驻多个 branch layout；
- 默认保存 best/final PL；
- 保存 optimizer state；
- 保存 bundle history；
- 保存 density map；
- 保存 snapshot；
- 引入通用 artifact database；
- 引入复杂 plugin architecture；
- 把 DCT/Poisson 放回正式路线；
- 为了“全局视野”改变 exact objective。

---

# 38. 最终执行原则

本轮代码应让日常研究变成：

```text
一个 external 7% checkpoint
        ↓
选择实验 module / optimizer / step
        ↓
顺序在内存中传 Layout
        ↓
必要时 temp PL 兼容旧 module
        ↓
每轮只记 scalar trajectory
        ↓
fresh exact audit
        ↓
写 params.json + experiment.md + trajectory.csv
        ↓
删除所有 temp PL / snapshot / cache
```

只有用户明确说：

> “保留这个实验的 stage X 布局”

才产生：

```text
saved/stage_X.pl
```

这应成为框架默认且难以误用的行为，而不是靠研究者每次手动删除文件。
