# h221 → 86M：`exact_joint_gp` 连续次梯度调参方案 V6

> 仓库：`D:\codex_project\HUAWEI_EDA\non-smooth`  
> 分支：`nonsmooth-gp-v1`  
> 方案基线提交：`0fbf1b9863420601e1442814e80eac242fc27b41`  
> 日期：2026-09-07  
> 性质：供后续低成本模型执行和持续填写的实验蓝图；本文不是完成证明。

## 0. 一句话目标

从历史中间结果 `h221_dct_poisson_512/best.pl` 出发，只使用 `exact_joint_gp` 的连续、无条件次梯度更新，在 300 轮预算内寻找：

```text
exact HPWL <= 86,000,000
canonical exact overflow < 7.0000%
```

整个搜索过程不使用 overflow 硬可行域、不使用 candidate accept/reject、不使用 exact backtracking；允许轨迹暂时超过 7%，但最终用于排名的状态必须满足 canonical exact overflow `< 7%`。

---

## 1. 不可变实验契约

### 1.1 输入与权威审计

```text
case: adaptec1
checkpoint:
D:\codex_project\HUAWEI_EDA\non-smooth\framework\experiment_logs\
20260906_233448_h375_raw_replay_local\h221_dct_poisson_512\best.pl

SHA-256:
9d0bd99966c3d679780ea73f0c28d9b3eb4163661fe679dfd85d7d32277e79b2

canonical evaluator:
bins_x = 512
bins_y = 512
target_density = 1.0
threads = 1
seed = 219（记录用；当前 stage 为确定性）
```

该 `.pl` 是历史 DCT/Poisson 平滑阶段的输出，只能作为只读起点；后续评价全部使用当前项目的 exact HPWL 与 exact rectangle/bin overlap。

加载 checkpoint 后必须执行 `clamp_movable`。当前权威起始值为：

| 指标 | 数值 |
|---|---:|
| exact HPWL | 109,837,622.009 |
| exact overflow | 7.8193002333% |

原始 `.pl` 中有 480 个可动节点位于 die 外。直接运行不带 clamp 语义的 `audit` 会得到约 11.52% 的边界堆积伪影，不能作为本方案起点。

### 1.2 “丝滑下降”的具体含义

所有正式 run 必须满足：

- `module = exact_joint_gp`；
- `exact_batch_acceptance.enabled = false`；
- 每一步只要数值有限就无条件应用，不因 HPWL 或 overflow 变差而回滚；
- 不设置 overflow hard ceiling，不做 feasibility gate，不做 line search/backtracking；
- optimizer state 在一个 300 轮 run 内连续保留；
- 只允许 die-boundary clamp 和 optimizer 的有限步长裁剪；它们是数值定义，不是候选接受机制；
- `selector` 只在 run 结束时从内存轨迹中选取最低 HPWL 的 `<7%` 状态，不参与逐步更新。

若出现 NaN/Inf，应将 run 记为失败并停止该 run；不得悄悄恢复旧布局后宣称成功。

### 1.3 公平性

除正在研究的变量外，以下项固定：输入 hash、case、512×512/1.0 evaluator、threads=1、300 iterations、HPWL direction、无接受机制、retention 和日志频率。每个正式 run 均从原始 h221 checkpoint 重新加载，禁止从上一个调参 run 的末态续跑。

---

## 2. 已有 100 轮基线及其含义

共同参数：`trajectory` λ、`density_weight_scale=0.5`、`update_interval=5`、`trajectory_horizon=100`、`lr=0.002`、`maximum_delta=4`。

| optimizer | best feasible HPWL | overflow | 首次 `<7%` | λ 峰值 |
|---|---:|---:|---:|---:|
| AMSGrad | **99,136,915.82** | 5.4217% | 54 | 51 |
| Adam | 99,960,810.50 | 5.5886% | 56 | 51 |
| AdaGrad | 101,850,109.32 | 6.9879% | 75 | 130 |
| Dual Averaging | 105,949,370.15 | 6.9907% | 96 | 131 |
| SGD | 108,191,594.03 | 5.2690% | 10 | 3.2 |
| Heavy Ball | 108,228,173.46 | 5.2480% | 11 | 3.3 |
| Normalized SGD | 108,985,009.69 | 4.9752% | 21 | 3.4 |

目标要求从起点下降 23,837,622.009（21.7026%）。当前 AMSGrad 已下降 10,700,706.192（9.7423%），但距离 86M 仍有 13,136,915.817（相对当前 13.2513%）。因此直接把原配置从 100 轮延长到 300 轮只能作为 control，不能视为有根据的主方案。

现有曲线还给出两个调参信号：

1. AMSGrad/Adam 的第一次步进把 overflow 推到约 19%，峰值达到 22%–23%，说明 `lr=0.002` 与当前裁剪允许很强的初始探索。
2. AMSGrad 在有效 λ 约 12 时首次回到 7% 以下，随后 λ 继续升到约 51，overflow 被压到 5.42%。这部分容量裕量没有直接服务于 86M 目标，后续应让 λ 围绕目标边界调节，而不是单向持续增大。

---

## 3. 目标驱动的 λ 设计

### 3.1 首选：用 `ratio` 策略把 86M/7% 编成平衡点

`ratio` 控制器比较：

```text
hpwl_ratio     = current_hpwl / hpwl_baseline
overflow_ratio = current_overflow / overflow_baseline
```

取：

```text
hpwl_baseline = 109,837,622.009
target_hpwl_ratio = 86,000,000 / 109,837,622.009
                  = 0.7829739795
```

为了让 `HPWL=86M, overflow=7%` 附近对应两个 ratio 相等：

```text
overflow_baseline* = 0.07 / 0.7829739795
                   = 0.0894027156
```

因此 λ 主扫描中心不是随意的 `0.07`，而是：

```json
"lambda": {
  "policy": "ratio",
  "density_weight_scale": 0.5,
  "update_interval": 5,
  "hpwl_baseline": 109837622.009,
  "overflow_baseline": 0.0894027156,
  "stop_overflow": 0.0695,
  "control_min": 1e-6,
  "control_max": 1e12
}
```

`stop_overflow=0.0695` 仅用于 λ 反馈留出数值余量；stage 结束时的 selector 使用 `overflow_cap=0.0699`。二者都不构成逐步接受条件。

中心值左右只做小范围扫描：

```text
overflow_baseline ∈ {0.0850, 0.0894027156, 0.0940}
```

解释：较小值更重视 density，较大值更愿意贴近 7% 换取 HPWL。`0.0894027156` 是与 86M 目标直接对应的 control。

### 3.2 对照：保留 `trajectory` 策略但拉长回收时间

300 轮下使用：

```text
trajectory_horizon ∈ {180, 240, 270}
density_weight_scale = 0.5
update_interval = 5
stop_overflow = 0.0695
```

`horizon=240` 为中心：前段允许线长探索，后 60 轮用于贴回目标带。`horizon=300` 不作为首选，因为没有稳定尾段；`horizon=100` 仅作为原机制延长 control。

可选诊断：在 `density_weight_scale=0.5` 时，已观测首轮 `lambda_base≈0.622172`。若 trajectory 继续把 overflow 压到明显低于 6%，可测试有效 λ 上限 `{12, 20, 40}`，对应近似：

```text
control_max = lambda_effective_max / 0.622172
            ≈ {19.3, 32.2, 64.3}
```

这只是 λ 数值范围消融，不是 placement 可行域，也不是 candidate gate。每次改变 `density_weight_scale` 后必须重新读取首轮 `lambda_base`，不能照抄上述 control 值。

---

## 4. 固定默认配置

除每阶段明确改变的字段外，默认配置如下：

```json
{
  "iterations": 300,
  "hpwl_direction": {
    "epsilon": 125.0,
    "active_power": 4.0,
    "degree_limit": 100
  },
  "optimizer": {
    "name": "amsgrad",
    "learning_rate": 0.002,
    "maximum_delta": 4.0,
    "beta1": 0.90,
    "beta2": 0.99
  },
  "lambda": {
    "policy": "ratio",
    "density_weight_scale": 0.5,
    "update_interval": 5,
    "hpwl_baseline": 109837622.009,
    "overflow_baseline": 0.0894027156,
    "stop_overflow": 0.0695,
    "control_min": 1e-6,
    "control_max": 1e12
  },
  "exact_batch_acceptance": {
    "enabled": false,
    "max_backtracks": 0
  },
  "selector": {
    "mode": "lowest_hpwl_under_overflow_cap",
    "overflow_cap": 0.0699
  },
  "log_every": 1,
  "adaptive_epsilon": false
}
```

注意：当前 `exact_joint_gp` 适配器支持 optimizer 的 `name/learning_rate/maximum_delta/beta1/beta2/momentum`，以及 λ 的 `policy/density_weight_scale/update_interval/trajectory_horizon/hpwl_baseline/overflow_baseline/control_min/control_max/stop_overflow`。不要把未接线字段写入 JSON 后假定其生效。

---

## 5. 分阶段实验矩阵

所有正式对照均运行 300 轮。每一阶段只晋级满足数值有限且曾达到 `<7%` 的候选；排名先看 best feasible HPWL，再看目标是否命中。

### A. 复现与 300 轮 control（2 runs）

| ID | optimizer | λ policy | horizon | lr | 目的 |
|---|---|---|---:|---:|---|
| A0 | AMSGrad | trajectory | 100 | 0.002 | 原 100 轮配置的直接延长 |
| A1 | AMSGrad | trajectory | 240 | 0.002 | 与 300 轮预算匹配的轨迹 |

若 A0 前 100 轮不能复现既有曲线到打印精度，停止调参并检查输入 hash、clamp、代码 commit 和参数解析。

### B. λ 主筛选（9 runs）

先只改变一个 λ 变量：

| ID | policy | overflow baseline / horizon | density scale | interval |
|---|---|---:|---:|---:|
| B0 | ratio | 0.0894027156 | 0.5 | 5 |
| B1 | ratio | 0.0850 | 0.5 | 5 |
| B2 | ratio | 0.0940 | 0.5 | 5 |
| B3 | ratio | B0 最优值 | 0.25 | 5 |
| B4 | ratio | B0 最优值 | 1.0 | 5 |
| B5 | ratio | B0 最优值 | B0/B3/B4 最优值 | 2 |
| B6 | ratio | B0 最优值 | B0/B3/B4 最优值 | 10 |
| B7 | trajectory | 180 | 0.5 | 5 |
| B8 | trajectory | 270 | 0.5 | 5 |

晋级 B 阶段 best feasible HPWL 最低的 2 个配置。临时 peak overflow 只作为诊断，不作为拒绝条件；若出现 NaN/Inf，则该 run 明确失败。

### C. 自适应 optimizer 与步幅筛选（最多 14 runs）

AMSGrad/Adam 本身已经提供逐坐标自适应步长，先利用现有实现，不急于新增控制器。

在 B 阶段前 2 个 λ 配置上扫描：

```text
optimizer ∈ {amsgrad, adam}
learning_rate ∈ {0.0005, 0.001, 0.002, 0.004}
maximum_delta = 4.0
```

这是 2 个 λ 配置分别配其原晋级 optimizer；若两个晋级项相同 optimizer，则仍各做 4 个 lr，共 8 runs。随后在最佳 `λ × optimizer × lr` 上扫描：

```text
maximum_delta ∈ {0.5, 1.0, 2.0, 4.0}
```

`maximum_delta<1` 会实质裁剪 Adam/AMSGrad 的初始坐标步，可检验“减少最初 19%–23% overflow 冲击是否给后续 HPWL 留出更多有效轮次”。最后仅对当前最佳配置补两个 optimizer control：

```text
optimizer ∈ {adagrad, dual-averaging}
```

不要重新全扫 7 个 optimizer；100 轮证据已经排除了大部分低收益组合。

### D. 动量参数精调（5 runs）

只在 C 阶段最佳配置上做：

| ID | beta1 | beta2 |
|---|---:|---:|
| D0 | 0.90 | 0.99 |
| D1 | 0.80 | 0.99 |
| D2 | 0.95 | 0.99 |
| D3 | 0.90 | 0.95 |
| D4 | 0.90 | 0.999 |

若最佳 optimizer 不是 Adam/AMSGrad，则跳过本阶段，改为只调该 optimizer 已接线的参数。

### E. 可选的平滑全局步长机制（3 runs，只有前四阶段停滞时才实现）

触发条件：A–D 均完成，最佳 `<7%` HPWL 仍高于 90M，且最后 50 轮 HPWL 下降不足 0.2%。

在 `exact_joint_gp` 内增加一个很小的 learning-rate schedule 参数面，不引入接受机制：

```text
constant
cosine_floor
overflow_ema
```

建议公式：

```text
cosine_floor:
lr(t) = lr_min + 0.5 * (lr0 - lr_min) * (1 + cos(pi*t/(T-1)))
lr_min = 0.25 * lr0

overflow_ema:
o_bar(t) = 0.9 * o_bar(t-1) + 0.1 * overflow(t)
lr(t) = lr0 / (1 + k * softplus((o_bar(t)-0.07)/0.01))
k ∈ {0.5, 1.0}
```

该机制只连续改变步长，不能拒绝、回滚或重试任何一步。正式比较保持同一 λ、optimizer、300 轮预算，并记录每轮实际 `learning_rate`。修改代码后必须补配置解析和 schedule 单元测试，并按交接规则更新架构/模块/状态文档与 `CHANGELOG.md`。

### F. 命中后的确认

一旦出现 `HPWL<=86M && overflow<7%`：

1. 使用完全相同参数、同一 commit、同一 checkpoint 再运行一次，验证确定性；
2. fresh exact audit 最终内存选择；
3. 再跑一个 `selector overflow_cap=0.0695` 的安全裕量版本；
4. 默认仍不保存 `.pl`。如果之后确需 checkpoint，等待用户明确要求，再用相同参数单独复跑并显式 `--save-stage exact_joint_gp`。

命中目标后停止低优先级宽扫，把算力用于确认和结果整理。

---

## 6. 排名、停止与失败判据

### 6.1 主判据

```text
success = best_feasible_exact_hpwl <= 86,000,000
          && best_feasible_exact_overflow_percent < 7.0000
```

严格使用 `<7%`，不是 `<=7%` 的打印近似。推荐 selector cap 为 6.99%，λ target 为 6.95%。

### 6.2 词典序排名

1. 是否命中 86M/<7%；
2. 是否存在 `<7%` 状态；
3. best feasible exact HPWL；
4. best feasible overflow；
5. 最后 50 轮 feasible HPWL 斜率；
6. runtime。

best-any HPWL 若 overflow≥7%，只能用于说明探索深度，不能作为最终结果。

### 6.3 诊断指标

每个 run 还应计算：

- first iteration `<7%`；
- best-any HPWL 及其 overflow；
- best-feasible iteration；
- peak overflow；
- `sum_t max(overflow_percent_t-7,0)`（不可行区面积）；
- λ 初值、峰值、末值；
- 前 20、100、200、300 轮 HPWL；
- 最后 50 轮 HPWL 线性斜率；
- 最大实际坐标步和裁剪比例（若代码已有/新增遥测）；
- NaN/Inf、异常退出和输入 hash 是否变化。

### 6.4 决策门槛

- 若 B 阶段所有 ratio run 都不能回到 `<7%`，先降低 `overflow_baseline` 或提高 `density_weight_scale`，不要进入 lr 扩张。
- 若 run 很早降到 `<6%` 且 HPWL 平台，优先减小 density scale、增大 overflow baseline 或降低 λ 上限。
- 若 peak overflow 很高且大部分预算都在回收，优先减小 lr 或 `maximum_delta`，而不是继续增大 λ。
- 若 300 轮末尾 HPWL 仍明显下降，下一轮优先增加预算到 500；不要在本轮偷偷改变 300 轮公平预算。
- 若 A–D 最佳仍高于 90M 且已平台，进入 E；若仍不改善，应判断当前方向/预条件器是否构成结构性瓶颈，而不是无限细扫参数。

---

## 7. 记录与产物保留

### 7.1 禁止保存中间布局

正式命令不得带 `--save-stage`。保持：

```text
snapshot_every = 0
output_dir = empty
```

禁止保存：中间/最终 `.pl`、snapshot、optimizer state、gradient dump、density map、occupancy cache、bundle history、candidate dump。实验目录仍只能有：

```text
params.json
experiment.md
trajectory.csv
```

这里的“不要保存中间结果”指不保存布局和大对象；逐轮标量是收敛证据，必须记录。

### 7.2 逐轮标量记录

当前 stage-boundary `trajectory.csv` 不含内部 300 轮轨迹。执行器应采用以下轻量方式：

1. `log_every=1`；
2. stdout 只写到 `%TEMP%\nonsmooth-gp\<run_id>\stdout.log`；
3. 每个 run 完成后立即解析 `[GP]` 行；
4. 把标量追加到一个合并 CSV；
5. 删除临时 stdout 和临时 JSON 工作区；
6. 输入 checkpoint 前后重新计算 SHA-256。

若为 E 阶段改代码，建议 stdout 每轮至少包含：

```text
iteration, exact_hpwl, overflow_percent, density_energy, max_density,
lambda_base, lambda_control, lambda_effective, learning_rate,
optimizer, maximum_delta
```

汇总目录建议：

```text
framework/results/analysis/YYYYMMDD_h221_86m_tuning/
├─ README.md                       # 假设、结论、失败说明
├─ runs.csv                        # 每个 run 一行
├─ per_iteration_metrics.csv       # 仅标量轨迹
├─ convergence_hpwl_overflow.png
├─ lambda_step_diagnostics.png
└─ hpwl_overflow_phase_portrait.png
```

这些是本地分析产物，默认不提交批量轨迹或图片。长期文档只提交本方案、关键结果表和必要的交接摘要。

### 7.3 `runs.csv` 字段

```text
run_id,commit,input_sha256,optimizer,lr,maximum_delta,beta1,beta2,
lambda_policy,density_weight_scale,update_interval,trajectory_horizon,
hpwl_baseline,overflow_baseline,control_min,control_max,
best_any_hpwl,best_any_overflow_percent,
best_feasible_hpwl,best_feasible_overflow_percent,best_feasible_iteration,
first_under_7_iteration,peak_overflow_percent,infeasible_area,
lambda_initial,lambda_peak,lambda_final,last50_hpwl_slope,
runtime_seconds,status,target_hit,notes
```

不得把缺失值填成 0；使用空值或 `NA`。失败 run 也保留一行并记录原因。

---

## 8. 收敛图布置

### 图 1：主收敛图（2×1，共享横轴）

- 上图：`exact HPWL / 1e6` 对 iteration；加 86M 水平虚线。
- 下图：`overflow_percent` 对 iteration；加 7% 虚线，并浅色标出 `overflow>7%` 区域。
- 主图只画 A0、目标校准 B0 和当前各阶段 top-1，避免 30 条线挤在一起。
- 同一 run 在两张子图使用完全相同颜色；实线表示最终可行，虚线表示全程未找到可行点。

### 图 2：λ 与步长诊断（2×1）

- 上图：`lambda_effective`，y 轴用 log scale。
- 下图：实际 `learning_rate`；常数步长也必须画出，便于与 E 阶段公平比较。
- 在 first-under-7 和 best-feasible iteration 处画竖线。

### 图 3：HPWL–overflow 相图

- x 轴 `overflow_percent`，y 轴 `HPWL/1e6`；
- 用颜色深浅表示 iteration 0→299；
- 标出起点、best-any、best-feasible 和末态；
- 阴影标出 `overflow<7% && HPWL<=86M` 目标区域。

相图比单独两条时间曲线更容易判断：算法是在“低 HPWL/高 overflow”与“高 HPWL/低 overflow”之间往返，还是整体向目标左下角推进。

### 附录图

所有 run 可用低透明度画在 appendix 图中；不要把完整 sweep 全塞进主图。PNG 建议 180–220 dpi，标题注明 checkpoint hash 前 8 位、commit、300 iterations、512×512/1.0。

---

## 9. 可执行命令模板

执行器可在被 Git 忽略的 `build\h221_86m_tuning\` 下生成每个 run 的 config 和 pipeline；完整解析参数会自动写入该 run 的 `params.json`。

```powershell
Set-Location 'D:\codex_project\HUAWEI_EDA\non-smooth'

$Nsgp = 'D:\codex_project\HUAWEI_EDA\non-smooth\build\Release\nsgp.exe'
$Checkpoint = 'D:\codex_project\HUAWEI_EDA\non-smooth\framework\experiment_logs\20260906_233448_h375_raw_replay_local\h221_dct_poisson_512\best.pl'
$Pipeline = 'D:\codex_project\HUAWEI_EDA\non-smooth\build\h221_86m_tuning\pipeline_B0.json'
$RunId = '20260908_adaptec1_h22186_B0_ratio_target'

& $Nsgp run `
  --case adaptec1 `
  --placement $Checkpoint `
  --pipeline $Pipeline `
  --threads 1 `
  --run-id $RunId
```

不得加入 `--save-stage`。若本机构建产物为 `build\nsgp.exe`，先核对其 commit 元数据再替换 `$Nsgp`，不要混用旧二进制。

每一批实验前：

```powershell
git branch --show-current
git rev-parse HEAD
git status --short
Get-FileHash -Algorithm SHA256 -LiteralPath $Checkpoint
ctest --test-dir build -C Release --output-on-failure
```

每个 run 后运行 retention contract，并核对实验目录没有 `.pl` 或 snapshot。

---

## 10. 活结果台账（后续执行器必须回填）

### 10.1 批次信息

| 字段 | 值 |
|---|---|
| 执行日期 | 2026-09-08 |
| Git commit / dirty | 基线 0fbf1b9（dirty：07/CHANGELOG 交接文档与 V6 方案文件） |
| 输入 SHA-256 | `9d0bd99966c3d679780ea73f0c28d9b3eb4163661fe679dfd85d7d32277e79b2`（每批前复核，全程未变） |
| binary 路径 | `build/nsgp.exe`（`scripts/build_manual_gcc.ps1` 手工构建，本机无 CMake；宏元数据 commit=0fbf1b9） |
| CTest | CTest 不可用（无 CMake）；等价物 `build/nsgp_tests.exe` 与 `build/nsgp_adaptive_tests.exe` 全部通过 |
| retention contract | 30 个实验目录（29 调参 + 1 smoke）全部通过 |

### 10.2 结果表

| ID | optimizer | lr | max Δ | λ policy/关键参数 | best-any HPWL @ ovf | best feasible HPWL @ ovf | first `<7%` | peak ovf | λ peak | last-50 slope | runtime | 状态/结论 |
|---|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---|
| A0 | AMSGrad | .002 | 4 | traj h=100 | 98,900,700 @ 4.785% | 98,900,721.01 @ 4.7852% | 54 | 22.42% | 412.2 | -326.0 | 53.7s | 复现通过：前 100 轮与 V5 基线一致到打印精度 |
| A1 | AMSGrad | .002 | 4 | traj h=240 | 98,387,100 @ 4.764% | 98,387,144.68 @ 4.7642% | 54 | 22.42% | 81.5 | -1619.7 | 52.7s | 优于 A0；4.8% 深可行域，非 86M 路径 |
| B0 | AMSGrad | .002 | 4 | ratio ob=.0894027156 | 94,878,500 @ 7.125% | 101,458,303.23 @ 6.9229% | 42 | 22.42% | 7.9 | -8732.6 | 50.9s | 可行但 HPWL 差 |
| B1 | AMSGrad | .002 | 4 | ratio ob=.0850 | 95,051,900 @ 7.019% | 101,036,898.48 @ 6.9883% | 41 | 22.42% | 8.9 | -6478.8 | 39.5s | B0-B2 中最优 ob |
| B2 | AMSGrad | .002 | 4 | ratio ob=.0940 | 94,878,500 @ 7.125% | 101,458,303.23 @ 6.9229% | 42 | 22.42% | 7.9 | -8732.6 | 17.4s | 与 B0 逐位一致：interval 5 时早期 1.25 增长项主导，ob 不分化 |
| B3 | AMSGrad | .002 | 4 | ratio ob=.0850 ds=.25 | 97,078,000 @ 9.033% | 97,315,026.07 @ 7.3081% | 从未 | 24.40% | 33230.9 | +69.2 | 51.8s | 无 <7% 状态，失败 |
| B4 | AMSGrad | .002 | 4 | ratio ob=.0850 ds=1.0 | 95,152,200 @ 7.272% | 101,716,342.08 @ 6.9831% | 49 | 21.02% | 4.7 | -9253.0 | 52.6s | ds=1.0 过强 |
| B5 | AMSGrad | .002 | 4 | ratio ob=.0850 itv=2 | 94,890,000 @ 7.693% | **95,827,835.14 @ 6.9621%** | 21 | 22.42% | 16.7 | -9502.1 | 16.2s | interval 2 激活 ratio 策略，晋级 |
| B6 | AMSGrad | .002 | 4 | ratio ob=.0850 itv=10 | 94,409,000 @ 7.350% | 96,123,329.43 @ 7.1310% | 从未 | 22.42% | 16.6 | -5672.3 | 16.0s | 无 <7% 状态 |
| B7 | AMSGrad | .002 | 4 | traj h=180 | 98,647,600 @ 4.739% | 98,647,648.41 @ 4.7390% | 54 | 22.42% | 192.5 | -572.9 | 16.5s | 与 traj 系一致 |
| B8 | AMSGrad | .002 | 4 | traj h=270 | 98,262,500 @ 4.811% | **98,262,496.79 @ 4.8107%** | 54 | 22.42% | 53.6 | -2668.2 | 16.4s | traj 系最优，晋级 |
| C_B5_lr0005 | AMSGrad | .0005 | 4 | ratio ob=.0850 itv=2 | 97,593,800 @ 8.335% | 97,837,619.00 @ 6.9776% | 14 | 13.71% | 18.8 | -9479.2 | 16.8s | lr 过小 |
| C_B5_lr001 | AMSGrad | .001 | 4 | ratio ob=.0850 itv=2 | 96,223,000 @ 8.037% | 96,491,160.77 @ 6.9431% | 16 | 18.87% | 14.7 | -10382.4 | 38.0s | 次于 lr .002 |
| C_B5_lr004 | AMSGrad | .004 | 4 | ratio ob=.0850 itv=2 | 97,892,900 @ 8.971% | 99,386,752.93 @ 6.9810% | 36 | 23.90% | 171.1 | -23946.6 | 18.2s | lr 过大 |
| C_B8_lr0005 | AMSGrad | .0005 | 4 | traj h=270 | 100,592,000 @ 4.821% | 100,592,000.66 @ 4.8207% | 41 | 13.71% | 13.5 | -3732.7 | 16.5s | lr 过小 |
| C_B8_lr001 | AMSGrad | .001 | 4 | traj h=270 | 99,285,400 @ 4.781% | 99,285,435.05 @ 4.7814% | 45 | 18.87% | 21.9 | -2978.9 | 16.4s | 次于 lr .002 |
| C_B8_lr004 | AMSGrad | .004 | 4 | traj h=270 | 98,023,500 @ 4.794% | 98,023,517.80 @ 4.7945% | 58 | 23.90% | 26.8 | -14847.5 | 17.7s | traj 系最佳 lr，整体仍逊 B5 |
| C_B5_md05 | AMSGrad | .002 | 0.5 | ratio ob=.0850 itv=2 | 94,526,000 @ 7.791% | 94,746,041.05 @ 6.9084% | 19 | 19.46% | 13.8 | -9966.5 | 49.8s | md 单调趋势成立 |
| C_B5_md1 | AMSGrad | .002 | 1.0 | ratio ob=.0850 itv=2 | 95,109,700 @ 7.388% | 95,226,587.83 @ 6.9014% | 21 | 22.42% | 16.6 | -7694.6 | 27.4s | - |
| C_B5_md2 | AMSGrad | .002 | 2.0 | ratio ob=.0850 itv=2 | 95,120,900 @ 7.374% | 95,237,671.77 @ 6.9001% | 21 | 22.42% | 16.6 | -7491.5 | 25.4s | - |
| C_B5_md025 | AMSGrad | .002 | 0.25 | ratio ob=.0850 itv=2 | 93,983,000 @ 7.704% | **94,061,868.14 @ 6.7474%** | 18 | 13.79% | 12.3 | -12820.1 | 48.9s | **全程最佳** |
| C_B5_adagrad | AdaGrad | .002 | 0.5 | ratio ob=.0850 itv=2 | 99,700,500 @ 7.286% | 103,492,172.23 @ 6.9371% | 20 | 19.80% | 11.1 | -7855.5 | 46.7s | 远逊 amsgrad |
| C_B5_dualavg | DualAvg | .002 | 0.5 | ratio ob=.0850 itv=2 | 102,705,000 @ 16.370% | 107,167,246.43 @ 6.9609% | 90 | 20.32% | 4.5e11 | +569347 | 50.4s | λ 失控发散 |
| D1 | AMSGrad | .002 | 0.25 | β1=.80 | 95,561,800 @ 8.387% | 96,008,149.87 @ 6.5556% | 18 | 13.76% | 14.8 | -3477.7 | 50.0s | 劣于基线动量 |
| D2 | AMSGrad | .002 | 0.25 | β1=.95 | 94,909,600 @ 9.673% | 97,065,074.84 @ 6.9293% | 20 | 13.82% | 82.5 | -15336.4 | 49.4s | 劣于基线动量 |
| D3 | AMSGrad | .002 | 0.25 | β2=.95 | 93,896,000 @ 8.015% | 94,213,110.62 @ 6.8036% | 17 | 13.79% | 12.0 | -11414.7 | 48.0s | 接近但劣于基线 |
| D4 | AMSGrad | .002 | 0.25 | β2=.999 | 93,929,800 @ 8.259% | 94,240,008.14 @ 6.8780% | 18 | 13.79% | 15.6 | -8630.0 | 48.1s | 接近但劣于基线 |
| D_probe_md0125 | AMSGrad | .002 | 0.125 | ratio ob=.0850 itv=2（诊断扩展） | 93,690,700 @ 8.205% | 94,188,027.39 @ 6.8695% | 17 | 9.78% | 15.4 | -12186.8 | 47.9s | md 趋势在 0.25 反转 |

后续每完成一个 run 就追加一行；不得只保留胜者。完成一个阶段后，在表下写三句话：观察、证据、下一阶段晋级理由。

### 10.3 最终结论模板

```text
最佳参数：exact_joint_gp + AMSGrad(lr=0.002, maximum_delta=0.25, beta 0.90/0.99)
  + ratio λ（hpwl_baseline=109837622.009, overflow_baseline=0.085,
  density_weight_scale=0.5, update_interval=2, stop=0.0695），selector cap
  0.0699，300 轮，threads=1。
最佳 exact HPWL：94,061,868.141472（run C_B5_md025）
对应 exact overflow：6.7474155990704%
是否命中 86M/<7%：否（差 8,061,868.14，需再降 8.57%）
相对 h221 起点改善：-15,775,753.87（-14.36%）
相对 100 轮 AMSGrad 基线改善：-5,075,047.68（-5.12%）
首次进入 <7% 的 iteration：18
best-feasible iteration：全程窗口最低 HPWL（约第 299 轮区域）
关键 λ/步长行为：λ 峰值 12.3、围绕目标带小幅调节；maximum_delta=0.25 把首轮
  overflow 冲击从 22.4% 压到 13.8%，是本轮最大单项收益（约 -108 万 HPWL）；
  md=0.125 时趋势反转。
失败或局限：未命中 86M；ratio 策略在 interval=5 时早期不分化（B0/B2 逐位一致）；
  trajectory 系把 overflow 压到 4.7-4.8% 但距 86M 更远；adagrad/dual-averaging
  明显更差（dual-averaging λ 发散到 4.5e11）；单 case、单 checkpoint、单确定性
  轨迹；E 阶段触发条件未满足（最后 50 轮下降 0.647% > 0.2%），未实现 schedule。
是否保存 placement：否（除非用户之后明确要求）
```

---

## 11. 执行顺序清单

- [x] 重读 `AGENTS.md` 和 `you should know about/` 全部 Markdown。
- [x] 核对 branch、HEAD、dirty 状态、binary commit 和 checkpoint SHA-256。
- [x] 编译并通过 CTest（本机无 CMake，等价运行两个测试可执行）。
- [x] 1 轮 smoke：确认 ratio 参数被解析、无 acceptance、三文件 retention。
- [x] 执行 A0/A1，验证已有前 100 轨迹（A0 iter54/99 与 V5 一致到打印精度）。
- [x] 执行 B0–B8，晋级 top-2 λ 配置（B5、B8）。
- [x] 执行 C，先 lr，再 maximum delta，再两个 optimizer control。
- [x] 仅在最佳配置上执行 D（含 md=0.125 诊断探针）。
- [x] E 触发条件未满足（最后 50 轮下降 0.647% > 0.2%），按方案不实现不执行。
- [x] 每个 run 都保留失败记录、合并逐轮标量、清理临时日志。
- [x] 绘制三类收敛图，主图只放 control 与各阶段优胜者。
- [x] 回填本文结果台账；关键新结果同步到 `05_CURRENT_STATE.md` 和 `CHANGELOG.md`。
- [x] 不保存任何中间或最终 placement（未命中目标，无需等待保存决定）。

## 12. 本方案不允许宣称的内容

- `<7%` overflow 不等于 legalized。
- h221 是历史 smooth-surrogate 链的中间结果，不是新的最终评价方法。
- 单 case、单确定性轨迹不能证明参数对八个 ISPD 2005 case 通用。
- best-any 的低 HPWL 不能替代 best-feasible。
- 未达到 86M 时不能用“仍在下降”冒充达标；应如实报告差距和斜率。
- 目标命中后仍需 fresh exact audit 和确定性复跑，不能只引用 stdout 的低精度打印值。
