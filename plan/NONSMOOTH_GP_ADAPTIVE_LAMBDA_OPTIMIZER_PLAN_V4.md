# nonsmooth-gp-v1：Adaptive Lambda GP + 可组装 Optimizer / Step Policy 修改方案 V4

> 面向本地 AI 的执行方案。  
> 仓库：`Bowen-Yuan/HUAWEI_EDA_GP`  
> 目标分支：`nonsmooth-gp-v1`  
> 本方案编写时远端 HEAD：`cca2e8b5d11ffe7ae5b7aed07a7b48984bfccd3d`  
> 日期：2026-09-07  
>
> **最高优先级：本方案是执行蓝图，不替代当前工作树。真正修改前必须重新读取 `AGENTS.md` 与 `you should know about/` 全部 Markdown，并检查 branch / HEAD / dirty 状态。若 HEAD 已变化，以当前源码和 fresh exact audit 为最终事实。**

---

## 0. 本轮目标

本轮只解决两个研究问题，不修改现有优化模块的算法行为：

1. **新增一个独立的 `adaptive_lambda_gp` stage**：
   - 从现有约 7% overflow 的 `adaptec1` external checkpoint 出发；
   - 不再把 7% 作为每一步的硬接受上限；
   - 允许短期甚至中等时间进入更高 overflow 区域，以换取 HPWL 探索；
   - 使用动态 `lambda` 和逐渐收紧的 overflow funnel/corridor，把搜索重新拉回 7%；
   - 模块内部始终保存一个**仅内存存在**的 `best_feasible` 布局；
   - stage 结束时输出的布局必须恢复为 canonical exact `overflow <= 7%` 的最佳 feasible 布局。

2. **建立轻量、可组装的 optimizer + step policy 机制**：
   - 复用当前已有 `ea::Optimizer` 工厂，不创建第二套 optimizer hierarchy；
   - 新增一个很小的、独立于 optimizer 的 step controller；
   - 新模块可以用 JSON 独立选择 optimizer 与 step policy；
   - 第一轮直接公平比较现有 7 种 optimizer，不先增加新的 optimizer 算法。

本轮不是重构 `global_view_gp`，也不是重写 `exact_joint_gp`。现有模块必须保持行为不变，以便作为 control/baseline。

---

# 1. 当前事实与问题定义

## 1.1 当前权威输入

正式数值实验从当前 handoff 文档记录的 `adaptec1` checkpoint 出发：

```text
D:\codex_project\HUAWEI_EDA\epsilon-active\experiments\
h375_a1_surplus_recovery10_exchange\run\global.pl
```

必须在实验前重新计算并核对：

```text
SHA-256:
ba0fd17e9111b5bcbd13ed746ab9857106d962652552d26076c6e82007551e6b
```

canonical evaluator：

```text
bins_x = 512
bins_y = 512
target_density = 1.0
```

当前记录的 exact 初始指标：

```text
HPWL     = 85,999,318.311947
Overflow = 6.99999030053%
```

输入 checkpoint 只读使用；实验前后 SHA-256 必须完全一致。

---

## 1.2 当前 `global_view_gp` 的核心限制

当前搜索大致采用：

```cpp
cap = initial_overflow + tiny_tolerance;

accept =
    finite(candidate)
    && candidate.overflow <= cap
    && candidate.hpwl < current.hpwl;
```

这意味着从约 7% checkpoint 启动后，搜索几乎被限制在一条非常薄的 feasible boundary 上。

已有 50-iteration V3 control / active-bundle / backtracking 证据中存在“零 accepted step”。因此本轮的研究问题不是继续微调固定 `lambda=0.2`，而是：

> **如果允许搜索暂时离开 7% feasible boundary，并通过动态 density penalty 在后半程主动回收 feasibility，是否能找到更低的最终 feasible HPWL？**

---

# 2. 不允许改变的数学与工程契约

以下契约与当前项目完全一致，本轮不得修改。

## 2.1 HPWL

```text
pin_x = node_center_x + pin_offset_x
pin_y = node_center_y + pin_offset_y

HPWL_e =
    max(pin_x) - min(pin_x)
  + max(pin_y) - min(pin_y)

HPWL = sum_e net_weight_e * HPWL_e
```

方向可以使用 epsilon-active / active-face 等非光滑方向，但报告与接受使用的 HPWL 必须是 exact max-minus-min。

## 2.2 Density / overflow

```text
occupancy_b =
    movable rectangle 与 bin 的精确相交面积
  + physical fixed macro 与 bin 的精确相交面积

rho_b = occupancy_b / bin_area

density_energy =
    sum_b 0.5 * bin_area
    * max(rho_b - target_density, 0)^2

overflow_ratio =
    sum_b bin_area
    * max(rho_b - target_density, 0)
    / total_movable_area
```

仍然满足：

- fixed macro 不移动并消耗容量；
- `terminal_NI` 不消耗物理容量；
- 内存 `Node.x/y` 是中心坐标；
- Bookshelf `.pl` 是左下角；
- 所有候选接受都使用 canonical fresh exact audit；
- overflow 低于 7% 不代表 legalized；
- 本项目仍不做 legalization / detailed placement。

## 2.3 Surrogate 使用边界

允许：

- epsilon-active direction；
- active-face ensemble；
- gradient normalization；
- dynamic lambda；
- trust radius；
- exact backtracking；
- temporary overflow corridor。

禁止：

- WA / LSE 替代最终 HPWL；
- DCT / Poisson / electrostatic density 替代最终 density；
- smooth merit 冒充最终 exact evaluator；
- 因性能跳过 candidate / stage boundary fresh exact audit。

---

# 3. 新模块：`adaptive_lambda_gp`

建议新增：

```text
modules/adaptive_lambda_gp/
├─ code/
│  ├─ adaptive_lambda_gp.hpp
│  ├─ adaptive_lambda_gp.cpp
│  └─ module.cpp
├─ params/
│  ├─ smoke.json
│  ├─ adam_trust.json
│  ├─ optimizer_screen.json        # 仅作为模板，不在一次 run 内循环 optimizer
│  └─ README.md                    # 可选，说明参数单位
└─ results/.gitkeep
```

不要修改：

```text
modules/global_view_gp/code/*
modules/exact_joint_gp/code/*
modules/exact_recovery/code/*
```

新模块只复用：

```text
ea::Database
ea::ExactHpwl
ea::ExactOverlapDensity
ea::Optimizer / make_optimizer()
microkernel StageContext / StageStats
```

---

# 4. 新模块的搜索状态

建议新模块维护以下轻量状态：

```cpp
struct AdaptiveSearchState {
    ExactMetrics initial;
    ExactMetrics current;
    ExactMetrics best_any;
    ExactMetrics best_feasible;

    std::vector<Real> best_feasible_positions; // movable only

    Real lambda;
    Real lambda_integral;
    Real previous_overflow;

    int accepted;
    int rejected;
    int objective_evaluations;
    int reject_streak;
    int optimizer_resets;

    Real maximum_observed_overflow;
    int first_return_to_feasible_iteration;
};
```

禁止保存：

- occupancy cache 跨 stage；
- optimizer state 到磁盘；
- bundle history 到磁盘；
- intermediate `.pl`；
- snapshot。

`best_feasible_positions` 只存在于当前进程内。

因为输入本身已经约 7% feasible，所以从 iteration 0 起一定存在 fallback feasible placement。

---

# 5. Direction 设计：先保持简单，避免一次引入过多变量

第一版 `adaptive_lambda_gp` 不建议同时引入 temporal bundle。

默认只使用两类 direction：

```text
wire_direction
density_direction
```

## 5.1 Wire direction

建议保留可选 active ensemble：

```text
epsilon/bin_scale = {0, 0.25, 0.5, 1.0}
```

每个 epsilon 得到 exact HPWL active direction，然后：

1. 对每个 direction 做 RMS normalization；
2. 在 simplex 上求 minimum-norm convex combination；
3. 得到 `wire_direction`。

smoke 可以先只用：

```text
epsilon = 0
```

正式默认建议：

```text
active_ensemble = true
```

## 5.2 Density direction

调用 canonical：

```cpp
ExactOverlapDensity::evaluate(
    epsilon = 0.0,
    active_power = 1.0,
    &gx,
    &gy
);
```

得到 exact-overlap density direction。

**重要变化：新模块中 wire 与 density direction 分别独立做 RMS normalization。**

因此组合：

```text
g = normalized_wire
  + lambda * normalized_density
```

最后再对 `g` 做一次 RMS normalization。

这样 `lambda` 成为近似无量纲的“方向权重”，不再高度依赖 wire / density 原始梯度量纲，便于跨 optimizer 做公平实验。

这只改变搜索方向，不改变 exact metric。

---

# 6. Dynamic Lambda：Funnel + Log-domain PI(D) Controller

本轮不直接复用或修改现有 `LambdaController`。

原因：

- 现有 controller 已服务于 `exact_joint_gp`；
- 它的现有语义不应因为本轮实验被改变；
- 新模块需要明确支持“先主动放宽 feasibility，再收紧”的 funnel 语义。

因此在 `adaptive_lambda_gp.cpp` 内实现一个小型私有：

```cpp
class FunnelLambdaController;
```

不建立通用 controller hierarchy。

---

## 6.1 Overflow funnel

设最终目标：

```text
O_final = 7.0%
```

默认建议：

```text
O_explore = 15.0%
O_hard    = 16.0%
```

15% 不是新的可行性标准，只是早期搜索 corridor。选择 16% hard ceiling 是为了给 15% corridor 留少量 exact/backtracking 数值裕量，同时仍阻止无界 density 恶化。

其中：

- `O_explore`：允许 HPWL 探索的动态 corridor 初始高度；
- `O_hard`：任何时期都不能超过的 safety ceiling；
- `O_final`：最终 stage 输出必须满足的 canonical exact overflow。

设进度：

```text
p = iteration / total_iterations
```

默认 schedule：

### Phase A — exploration hold

```text
0 <= p <= 0.20
corridor = 15.0%
```

允许搜索在 7–15% 区间中寻找更低 HPWL。

### Phase B — smooth contraction

```text
0.20 < p < 0.75
```

定义：

```text
u = (p - 0.20) / (0.75 - 0.20)
s = u*u*(3 - 2*u)       # smoothstep

corridor =
    O_final
  + (O_explore - O_final) * (1 - s)
```

因此 corridor 从 15% 平滑下降至 7%。

### Phase C — final feasible lock

```text
p >= 0.75
corridor = 7.0%
```

最后 25% iteration 不再放宽目标。

注意：

`corridor` 是搜索接受漏斗，不是 evaluator target density，也不能修改 canonical density grid。

---

## 6.2 Lambda 更新

建议 lambda 初值：

```text
lambda_initial = 0.25
lambda_min     = 0.03
lambda_max     = 16.0
update_interval = 5
```

由于 wire/density direction 已分别 normalized，lambda 可直接作为相对权重。

每 `update_interval` 次迭代更新一次。

定义百分比尺度：

```text
overflow_scale = 1.0% = 0.01 ratio
```

误差：

```text
e =
    (current_overflow - corridor)
    / 0.01
```

趋势：

```text
trend =
    (current_overflow - previous_update_overflow)
    / 0.01
```

积分项：

```text
I = clamp(0.8 * I + e, -4, 4)
```

建议初始控制参数：

```text
kp = 0.20
ki = 0.02
kd = 0.10
```

log-domain update：

```text
log_step =
    clamp(
        kp * e
      + ki * I
      + kd * trend,
      -0.35,
       0.55
    )

lambda =
    clamp(
        lambda * exp(log_step),
        lambda_min,
        lambda_max
    )
```

解释：

- overflow 低于 corridor 时，lambda 可以缓慢下降；
- overflow 高于 corridor 时，lambda 指数式增大；
- overflow 仍在继续恶化时，trend 项进一步提高 lambda；
- 所有变化有 log-step 上限，不允许 controller 一次爆炸数个数量级。

---

## 6.3 Final lock boost

进入：

```text
p >= 0.75
```

且：

```text
overflow > 7%
```

时增加恢复压力：

```text
extra =
    clamp(
        0.10 * (overflow - 0.07) / 0.01,
        0.0,
        0.25
    )

lambda *= exp(extra)
```

仍然 clamp 到 `lambda_max`。

不要直接把 lambda 设置成一个巨大常数。

---

# 7. Exact Funnel Acceptance

这是本轮与当前 strict-cap 模块最关键的行为差异。

所有 candidate 先进行 canonical fresh exact audit。

---

## 7.1 全局拒绝条件

以下任一成立直接 reject：

```text
NaN / Inf
candidate overflow > O_hard
position 超出 placement boundary（clamp 后仍异常）
```

默认：

```text
O_hard = 16%
```

这是安全红线，不是最终目标。

由于 15% 相比初始约 7% 是明显更宽的探索区，第一版必须额外记录 `maximum_observed_overflow` 和 `iterations_above_10_percent`。若搜索长期停留在 12–15% 而不能随 funnel 收紧回落，应视为 controller/step policy 失败，而不是继续放宽 hard ceiling。

---

## 7.2 当前状态在 corridor 内：探索模式

如果：

```text
current_overflow <= corridor
```

则允许：

```text
candidate_overflow <= corridor
AND
candidate_hpwl < current_hpwl - hpwl_tolerance
```

也就是说：

> 只要 candidate 仍在当前动态 funnel 中，就允许为了 HPWL 改善把 overflow 从 7% 提高到 9%、12%、甚至接近 15%。

这正是当前 strict 7% acceptance 无法完成的探索。

---

## 7.3 当前状态在 corridor 外：恢复模式

当 corridor 随时间下降后，当前布局可能暂时位于 corridor 之外。

这时不要求一步回到 7%，而允许连续恢复：

```text
candidate_overflow
    < current_overflow - overflow_improvement_epsilon
```

同时限制 HPWL 代价：

```text
candidate_hpwl
    <= current_hpwl
       * (1 + recovery_step_hpwl_budget)
```

并增加一个累计 guard：

```text
candidate_hpwl
    <= initial_hpwl
       * (1 + recovery_total_hpwl_budget)
```

建议初值：

```text
overflow_improvement_epsilon_percent = 0.002%
recovery_step_hpwl_budget_percent    = 0.05%
recovery_total_hpwl_budget_percent   = 0.30%
```

这些是第一轮实验起点，不是固定真值。

---

## 7.4 Backtracking

所有 optimizer / step policy 共用同一个 exact backtracking acceptance loop：

```text
scale = 1, 1/2, 1/4, ...
max_backtracks = 9
```

每个 trial 都 fresh exact audit。

这样 optimizer 对比时 acceptance contract 完全一致。

---

# 8. Best-feasible 保底与 stage 输出语义

这是新模块必须实现的安全机制。

每次得到 canonical exact：

```text
overflow <= 7.0% + numerical_tolerance
```

且：

```text
hpwl < best_feasible_hpwl
```

时：

```text
best_feasible_metrics = candidate
best_feasible_positions = current movable positions
```

stage 搜索结束后：

```cpp
restore(best_feasible_positions);
clamp_movable(db);
fresh_exact_audit();
```

无论 `last` 是否 feasible，都输出 `best_feasible`。

因此：

```text
中间状态：允许 >7%
stage 边界：必须 <=7%
```

如果最后 25% 恢复失败，最坏情况回到 iteration 0 的输入 feasible checkpoint，而不是把 >7% 状态交给后续 module。

不允许为了 fallback 写磁盘 checkpoint。

---

# 9. Optimizer 架构：复用现有 `ea::Optimizer`

当前 framework 已经支持：

```text
adam
amsgrad
adagrad
heavy-ball
sgd
normalized-sgd
dual-averaging
```

本轮不要重新定义：

```cpp
class NewOptimizerBase ...
```

不要在 `adaptive_lambda_gp` 中复制 Adam 实现。

直接继续使用：

```cpp
ea::parse_optimizer(...)
ea::make_optimizer(...)
optimizer->reset(...)
optimizer->compute_delta(...)
```

---

# 10. 新增轻量 Step Policy

当前 optimizer API 已经把：

```text
gradient
learning_rate
maximum_delta
```

分开，因此只需要增加一个很小的 step controller。

建议新增：

```text
framework/kernel/include/epsilon_active/step_policy.hpp
framework/kernel/src/step_policy.cpp
```

接口保持很小：

```cpp
enum class StepPolicyKind {
    Constant,
    CosineDecay,
    TrustRadius
};

struct StepDecision {
    Real learning_rate;
    Real maximum_delta;
};

struct StepObservation {
    bool accepted;
    int backtracks;
};

class StepController {
public:
    virtual ~StepController() = default;

    virtual void reset(
        int total_iterations,
        Real bin_size
    ) = 0;

    virtual StepDecision propose(
        int iteration
    ) = 0;

    virtual void observe(
        const StepObservation&
    ) = 0;
};
```

提供：

```cpp
parse_step_policy(...)
make_step_controller(...)
```

不要引入 DI / plugin / registry framework。

---

## 10.1 Constant

```text
learning_rate = lr0
maximum_delta = maximum_delta_bins * bin_size
```

它是最简单 control。

---

## 10.2 CosineDecay

```text
q = iteration / total_iterations

lr =
    lr_min
  + 0.5 * (lr0 - lr_min)
  * (1 + cos(pi*q))
```

默认：

```text
lr_min_ratio = 0.10
```

`maximum_delta` 保持固定。

用途：测试“后期自然缩步”是否比 trust reset 更稳定。

---

## 10.3 TrustRadius

初始：

```text
radius = trust_radius_bins * bin_size
```

建议：

```text
3 consecutive accepted:
    radius *= 1.25

rejected:
    radius *= 0.5
```

clamp：

```text
radius_min_bins = 0.02
radius_max_bins = 2.0
```

`learning_rate` 仍由 optimizer 使用，`maximum_delta = radius`。

backtracking 仍由统一 acceptance loop 实现，不放进 optimizer。

---

# 11. Optimizer state reset 规则

dynamic lambda 会改变 direction composition，动量类 optimizer 可能保留过时方向。

第一版使用统一、轻量 reset 规则：

```text
reject_streak >= 4
    -> optimizer.reset()

lambda_new / lambda_old >= 4
    -> optimizer.reset()

进入 final feasible lock
    -> reset once
```

同时 reset：

```text
step controller 的 accept/reject streak
```

但不重置 lambda controller。

记录：

```text
optimizer_resets
```

不要每次 lambda update 都 reset，否则 Adam / AdaGrad 无法积累状态。

---

# 12. 建议 JSON schema

示例：

```json
{
  "iterations": 200,

  "direction": {
    "active_ensemble": true,
    "epsilon_bin_scales": [0.0, 0.25, 0.5, 1.0],
    "active_power": 1.0
  },

  "optimizer": {
    "name": "adam",
    "learning_rate": 0.02,
    "beta1": 0.9,
    "beta2": 0.999,
    "momentum": 0.9,
    "numerical_epsilon": 1e-8
  },

  "step_policy": {
    "name": "trust",
    "maximum_delta_bins": 0.5,
    "trust_radius_bins": 1.0,
    "radius_min_bins": 0.02,
    "radius_max_bins": 2.0,
    "grow_factor": 1.25,
    "shrink_factor": 0.5,
    "grow_after_accepts": 3,
    "max_backtracks": 9
  },

  "lambda_policy": {
    "name": "funnel_pid",
    "initial": 0.25,
    "minimum": 0.03,
    "maximum": 16.0,
    "update_interval": 5,

    "kp": 0.20,
    "ki": 0.02,
    "kd": 0.10,

    "explore_overflow_percent": 15.0,
    "final_overflow_percent": 7.0,
    "hard_overflow_percent": 16.0,

    "hold_fraction": 0.20,
    "contract_end_fraction": 0.75
  },

  "acceptance": {
    "policy": "exact_funnel",
    "hpwl_tolerance_relative": 1e-12,
    "overflow_improvement_epsilon_percent": 0.002,
    "recovery_step_hpwl_budget_percent": 0.05,
    "recovery_total_hpwl_budget_percent": 0.30
  },

  "reset": {
    "reject_streak": 4,
    "lambda_ratio": 4.0,
    "reset_on_final_lock": true
  }
}
```

所有字段必须：

- 有默认值；
- 有边界校验；
- 自动进入 pipeline `module_chain.parameters` provenance；
- overflow 对外字段明确使用 `percent`。

---

# 13. CMake / registry 修改

允许修改 framework 注册与构建，但不修改现有算法模块。

## 13.1 CMake

在：

```text
NUMERIC_KERNEL_SOURCES
```

增加：

```text
framework/kernel/src/step_policy.cpp
```

在 `nsgp` sources 增加：

```text
modules/adaptive_lambda_gp/code/adaptive_lambda_gp.cpp
modules/adaptive_lambda_gp/code/module.cpp
```

新模块私有算法不放进另一个已有 module。

---

## 13.2 microkernel

在：

```text
framework/code/microkernel.hpp
```

只新增声明：

```cpp
void register_adaptive_lambda_gp(ModuleRegistry&);
```

在：

```text
framework/code/microkernel.cpp
```

只新增：

```cpp
modules::register_adaptive_lambda_gp(registry);
```

不要改变已有 module 注册顺序或实现。

---

# 14. Pipeline 配置

新增：

```text
framework/params/pipelines/adaptive_lambda_smoke.json
framework/params/pipelines/adaptive_lambda_gp.json
```

正式单 stage baseline：

```text
external H375 checkpoint
    ↓
adaptive_lambda_gp
    ↓
best feasible in-memory output
```

后续验证组合能力时再测试：

```text
external H375
    ↓
adaptive_lambda_gp
    ↓
exact_recovery
    ↓
equal_shape_swap
```

第一轮研究不要一开始就把 recovery/swap 加进主对照，否则无法判断 improvement 来自哪个 stage。

---

# 15. 测试要求

## 15.1 Step policy unit contracts

至少验证：

### Constant

```text
iteration 不改变 lr
maximum_delta 正确转换 bins -> placement units
```

### Cosine

```text
lr 单调不增
lr(0) = lr0
lr(last) 接近 lr_min
```

### Trust

```text
连续 accepted -> radius 增大
reject -> radius 缩小
radius 始终处于 min/max 内
```

---

## 15.2 Lambda controller synthetic contracts

不依赖 benchmark。

构造 synthetic overflow trajectory：

### Case A — below corridor

```text
overflow = 7.2%
corridor = 15%
```

要求：

```text
lambda 不应快速增大
```

### Case B — above corridor

```text
overflow = 8.5%
corridor = 7.5%
```

要求：

```text
lambda 增大
```

### Case C — worsening trend

```text
8.0% -> 8.5%
```

要求：

```text
lambda increase > 同误差但下降趋势的 increase
```

### Case D — final lock

```text
p > 0.75
overflow > 7%
```

要求：

```text
lambda 有额外 boost
```

### Case E — bounds

长时间极端误差后：

```text
lambda_min <= lambda <= lambda_max
```

---

## 15.3 Acceptance contracts

至少构造以下 exact-metric-only 小例：

```text
current:   HPWL 100, overflow 7.0
corridor: 9.0
candidate: HPWL  99, overflow 8.5
=> accept
```

```text
current:   HPWL 100, overflow 8.5
corridor: 8.0
candidate: HPWL 100.03, overflow 8.2
=> recovery accept（在预算内）
```

```text
candidate overflow > hard ceiling
=> reject
```

```text
candidate overflow worsens while current already outside corridor
=> reject
```

---

## 15.4 Best-feasible restore contract

synthetic 搜索结束在：

```text
last overflow = 8%
```

但历史 best feasible：

```text
overflow = 6.99%
```

要求：

```text
stage return 前 restore best feasible
fresh audit 后 <= 7%
```

---

## 15.5 Optimizer × Step smoke

对全部：

```text
7 optimizers
×
3 step policies
```

使用 synthetic normalized gradient 验证：

```text
finite
bounded delta
reset 后状态可重复
```

不要在 unit test 中判断哪个 optimizer “更优”。

---

# 16. 数值实验：从固定 SHA-256 checkpoint 开始

## 16.1 每次实验的固定条件

第一轮算法比较全部固定：

```text
case          = adaptec1
input SHA     = ba0fd17e...
grid          = 512 x 512
target        = 1.0
seed          = 219
threads       = 1
retention     = metrics-only
max backtrack = 9
input         = identical external checkpoint
```

先用 `threads=1` 做算法比较，减少 OpenMP reduction 顺序导致的微小数值差异。

只有选出最终候选后，再用相同参数做 `threads=8` runtime smoke。

---

# 17. Experiment A：验证“动态 lambda + overflow exploration”本身

## A0 — 现有 control

使用现有、未修改的：

```text
global_view_gp
optimizer = adam
step = trust
strict 7% cap
iterations = 50
```

已有历史证据是 50 轮存在零 accepted step，但本地 AI仍应在当前 HEAD / 当前输入 SHA 上重新跑一次 control，作为这次实验的同环境证据。

记录：

```text
initial HPWL / overflow
final HPWL / overflow
accepted / rejected
objective evaluations
runtime
```

---

## A1 — 新模块 control

```text
adaptive_lambda_gp
optimizer = adam
step = trust
iterations = 50
explore corridor = 15%
final target = 7%
```

目的不是直接证明最终最优，而是回答：

1. 是否出现当前 strict-cap 下不存在的 accepted move？
2. 是否实际进入过 >7% 区域，以及是否利用了 10–15% 的更宽探索区？
3. lambda 是否在 funnel 收紧后明显增大？
4. 是否重新回到 <=7%？
5. best-feasible HPWL 是否优于输入？

---

## A2 — 禁止探索 ablation

仍然使用新模块，但设：

```text
explore_overflow_percent = 7.0
hard_overflow_percent = 7.0 + tiny tolerance
```

其他完全相同。

比较：

```text
A1 vs A2
```

这是最重要的消融之一。

如果 A1 优于 A2，才能说明收益确实与“允许中期 infeasible exploration”有关，而不是单纯新代码路径造成。

---

## A3 — 固定 lambda ablation

新模块保持 funnel acceptance，但：

```text
lambda_policy = fixed
lambda = 0.25
```

比较：

```text
A1 dynamic lambda
vs
A3 fixed lambda
```

用于区分：

```text
收益来自 overflow corridor
还是
收益来自真正动态 lambda
```

---

# 18. Experiment B：Optimizer screening

只有 A1 能可靠返回 feasible 后再执行。

第一轮固定：

```text
dynamic lambda = identical
funnel          = identical
step policy     = trust
iterations      = 20
learning rate   = 0.02
threads         = 1
```

依次：

```text
B1 adam
B2 amsgrad
B3 adagrad
B4 heavy-ball
B5 sgd
B6 normalized-sgd
B7 dual-averaging
```

这只是 smoke screen。

淘汰条件：

```text
NaN / Inf
完全无法产生有效 proposal
反复撞 hard overflow ceiling
20 轮内全部 reject 且无恢复迹象
```

不要仅凭 20 轮 HPWL 排名直接下结论。

---

# 19. Experiment C：公平 learning-rate screening

从 B 中保留最多 4 个 optimizer。

每个 optimizer 使用相同三点预算：

```text
learning_rate ∈ {0.005, 0.02, 0.08}
iterations = 50
step = trust
```

每个 optimizer 都获得完全相同的 3 × 50 evaluation budget。

选择该 optimizer 自己的最佳 feasible setting。

主要排序：

```text
1. final selected overflow <= 7%
2. best feasible exact HPWL
3. accepted / rejected
4. objective evaluations
5. runtime
```

辅助观察：

```text
maximum overflow excursion
进入 >7% 的 iteration 数
首次重新回到 <=7% 的 iteration
lambda_initial / lambda_max / lambda_final
optimizer resets
```

---

# 20. Experiment D：Step policy screening

取 Experiment C 的前 2 个 optimizer。

固定各自最佳 learning rate，比较：

```text
constant
cosine
trust
```

每个：

```text
iterations = 100
```

同一 optimizer 内只改变 step policy。

不要同时改变 lambda 参数。

---

# 21. Experiment E：正式长跑

最终选：

```text
top 2 optimizer/step combinations
+
Adam/trust dynamic-lambda reference
```

执行：

```text
iterations = 200
threads = 1
```

如果 200 iteration 仍持续改善，可追加：

```text
iterations = 500
```

但 500 不作为第一轮默认预算。

---

# 22. 结果评价准则

## 22.1 必须满足

正式结果必须：

```text
final selected exact overflow <= 7%
input SHA-256 前后不变
无 NaN / Inf
stage boundary fresh exact audit 通过
metrics-only retention 通过
```

## 22.2 Primary metric

```text
best feasible exact HPWL
```

不能用：

```text
best-any HPWL
```

冒充结果，因为 best-any 可能发生在 >7% overflow。

## 22.3 必须同时报告

```text
initial HPWL
initial overflow %

best-any HPWL
best-any overflow %

best-feasible HPWL
best-feasible overflow %

last-before-restore HPWL
last-before-restore overflow %

final-selected HPWL
final-selected overflow %

maximum overflow excursion %
accepted
rejected
objective evaluations
optimizer resets
runtime
```

## 22.4 “成功”的表述

如果：

```text
best-feasible HPWL < control
```

可以说：

```text
在当前固定 checkpoint / evaluator / budget 下观察到 exact feasible HPWL 改善。
```

如果相对改善：

```text
>= 0.01%
```

可额外标记为“material improvement”。

如果只有几十或几百 HPWL 的微小差异，必须保留精确数值，不夸大。

---

# 23. 关于不同 optimizer 的研究假设

这些只是需要证伪的假设，不能在实验前写成结论。

## Adam

优点：

```text
已有项目经验最多
适合作 control
```

风险：

```text
dynamic lambda 变化后 moment 可能滞后
```

## AMSGrad

假设：

```text
对 penalty 权重增大阶段可能比 Adam 更稳定，
因为 second moment 不会回落。
```

## AdaGrad

假设：

```text
与非光滑 subgradient 相容性较好；
后期自然缩步可能有利于 feasibility recovery。
```

风险：

```text
早期探索后累计梯度可能导致后期过早冻结。
```

## NormalizedSGD

假设：

```text
因为新模块已经对 wire/density/combined direction 做 normalization，
它可能成为非常强的轻量 control。
```

## DualAveraging

假设：

```text
对长期一致的 subgradient 方向可能有效，
值得作为非光滑方法重点观察。
```

风险：

```text
当前问题非凸且 direction regime 随 lambda 改变，
历史累计可能成为负担。
```

## HeavyBall

风险优先：

```text
momentum 可能在 overflow funnel 收紧时产生明显 overshoot。
```

必须配 exact backtracking / trust。

## SGD

作为最简单基线保留。

---

# 24. 第一轮不要新增 RMSProp / Lion 等 optimizer

原因：

当前项目已经有 7 个 optimizer，但尚没有在**相同 dynamic-lambda + 相同 exact acceptance**条件下完成公平比较。

因此第一轮先回答：

> 现有 optimizer 是否已经足够？

只有在 Experiment C/D 证明现有方法存在明确缺口后，再单独提出第二次改动增加 RMSProp/Lion 等。

不要把“新增 optimizer 算法”和“动态 lambda 模块”混在同一个研究变量中。

---

# 25. 实验产物

遵守现有 metrics-only contract。

正常 experiment directory 仍然只允许：

```text
params.json
experiment.md
trajectory.csv
```

默认不保存：

```text
.pl
snapshot
checkpoint
gradient
occupancy
optimizer state
bundle history
```

用户明确要求时才：

```text
--save-stage adaptive_lambda_gp
```

并且只保存：

```text
saved/adaptive_lambda_gp.pl
```

---

# 26. 关于 iteration-level lambda 轨迹

第一版**不要为了记录 lambda 轨迹而重构整个 microkernel experiment schema**。

正式 pipeline 仍以当前三文件 stage-level记录为主。

开发 smoke 时允许：

```text
-- verbose console / stdout
```

打印轻量 iteration telemetry：

```text
iter
hpwl
overflow_percent
corridor_percent
lambda
step
accepted
backtracks
```

但不要生成第四个长期文件。

如果新算法取得正结果，再单独设计一个“generic stage telemetry sink”任务，避免本轮扩大架构修改范围。

---

# 27. 推荐实现顺序

## Step 1 — 只读确认

执行：

```powershell
git branch --show-current
git rev-parse HEAD
git status --short
git log -5 --oneline
```

然后重新读取：

```text
AGENTS.md
you should know about/*.md
```

核对 external checkpoint SHA。

---

## Step 2 — Step policy

先加入：

```text
step_policy.hpp
step_policy.cpp
```

写 unit contract。

不碰当前 modules。

---

## Step 3 — adaptive lambda controller

只写纯逻辑和 synthetic tests。

先不接真实 placement。

---

## Step 4 — 新 stage

接：

```text
ExactHpwl
ExactOverlapDensity
Optimizer
StepController
FunnelLambdaController
exact funnel acceptance
best-feasible restore
```

---

## Step 5 — registry / CMake

只增加新 module 注册，不修改已有 stage 实现。

---

## Step 6 — smoke

```text
1 iteration
5 iterations
20 iterations
```

检查：

```text
fresh audit
no NaN
best-feasible restore
retention
external hash unchanged
```

---

## Step 7 — A0/A1/A2/A3

先验证 dynamic-lambda 研究问题。

如果 A1 无法安全返回 7%，不要继续 optimizer sweep。

---

## Step 8 — B/C/D/E

逐级筛 optimizer / learning rate / step policy。

不要一次跑全部大矩阵。

---

# 28. 文档同步要求

因为会新增 module 和 shared step-policy primitive，同一个最终 commit 至少更新：

```text
you should know about/02_ARCHITECTURE_AND_CODE_MAP.md
you should know about/03_MODULE_CATALOG.md
you should know about/05_CURRENT_STATE.md
you should know about/07_FILE_INDEX.md
you should know about/CHANGELOG.md
```

如果本轮没有改变：

```text
exact metric
overflow 定义
retention schema
```

则不需要为了形式修改 `01_PROJECT_REQUIREMENTS.md` 或 `04_EXPERIMENT_RULES.md`。

若最终实际实现改变 trajectory schema，则必须同步更新 `04_EXPERIMENT_RULES.md`。

---

# 29. Git 要求

完成后：

```powershell
git diff --check
git status --short
ctest --test-dir build -C Release --output-on-failure
```

运行所需 smoke / experiment / retention tests。

确认没有进入 Git：

```text
dataset
external checkpoint
build/
*.exe
normal experiment directories
*.pl
snapshot
cache
```

建议一个提交完成本轮：

```text
feat: add exploratory adaptive-lambda GP stage
```

commit 内必须同时包含：

```text
source
params
tests
handoff docs
CHANGELOG
```

随后推送：

```text
origin/nonsmooth-gp-v1
```

---

# 30. 本轮明确不做

禁止顺手进行：

- 修改 `global_view_gp` acceptance；
- 修改 `exact_joint_gp` 原有 lambda controller；
- 删除历史实验；
- 重写 microkernel；
- 新增 plugin framework；
- 新增 artifact manager；
- 把 optimizer 重新复制进 module；
- 添加 smooth HPWL / smooth density；
- 把 temporary >7% 状态当成 final result；
- 为了得到好结果取消 final exact <=7% 约束。

---

# 31. 预期研究结论形式

本轮结束后应能回答三个独立问题：

## Q1

```text
允许 7% 以上的中期探索，
是否比 strict 7% cap 找到更低的最终 feasible HPWL？
```

由：

```text
A0 / A1 / A2
```

回答。

## Q2

```text
收益主要来自 overflow funnel，
还是 dynamic lambda 本身？
```

由：

```text
A1 / A3
```

回答。

## Q3

```text
在完全相同的 dynamic-lambda / exact acceptance 下，
哪一种 optimizer + step policy 最适合当前 non-smooth GP？
```

由：

```text
B / C / D / E
```

回答。

不要把三个问题混在一次实验中。

---

# 32. 当前可用的数值证据与尚未执行的部分

截至本方案编写时，仓库已有的真实数值证据包括：

```text
external H375 checkpoint:
HPWL     85,999,318.311947
overflow 6.99999030053%
```

以及当前 strict 7% global-view 路线已有 50-iteration 零 accepted-step 的历史记录。

**本方案中的 `adaptive_lambda_gp` 尚未实现，因此不能在不写代码的前提下伪造 A1/B/C/D/E benchmark 结果。**

本地 AI 完成实现后，必须从上述固定 SHA checkpoint 真实执行 Experiment A→E，并把 exact 数值写回：

```text
you should know about/05_CURRENT_STATE.md
you should know about/CHANGELOG.md
```

失败、零 accepted、无法回到 7%、optimizer 发散等负结果同样必须记录。

---

# 33. 最小成功标准

如果第一轮只做到以下几点，也算实现正确：

1. 新 module 完全独立，旧 module 行为不变；
2. starting SHA checkpoint 不被修改；
3. 中间可以 exact overflow >7%；
4. hard safety ceiling 生效；
5. lambda 会随 funnel violation 动态增减；
6. 末态无条件 restore best feasible；
7. stage boundary exact overflow <=7%；
8. 7 个 optimizer 与 3 个 step policy 都可通过 JSON 选择；
9. current numeric / retention contracts 全部通过；
10. A0/A1/A2/A3 得到可复现 exact 数值。

只有满足以上之后，才进入大规模 optimizer 排名。
