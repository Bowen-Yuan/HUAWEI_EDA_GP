# 05 — 当前状态、数值结果与已知缺口

记录日期：2026-09-08（Asia/Shanghai）。架构：registry-driven microkernel。接手时以 `git rev-parse HEAD` 为准。

## 1. 当前架构状态

已完成：

- `framework/kernel` 只保存 Bookshelf I/O、统一运行状态、exact HPWL/density 和 optimizer。
- `framework/code/microkernel.*` 提供注册/调用、fresh audit、threads、SHA-256 和三文件实验记录。
- 普通策略的权威源码和 `module.cpp` 位于 `modules/<stage>/code`，CMake 直接编译这些文件。
- 原 `framework/legacy` 权威目录和逐字节重复的策略源码已经移除。
- 未接线的第二套 `Problem/Layout/ModuleResult` 骨架已移除。
- `nsgp run` 和 `nsgp lab` 默认都是 metrics-only；`--save-stage` 是唯一普通 placement 保留入口。
- `global_place` 收到空 output path 时不写 metrics、snapshot 或 `last.pl`。
- `lab --threads` 会实际设置 OpenMP；pipeline 在入口统一设置 threads。
- pipeline params 记录完整解析后 module chain；lab params/summary 记录输入 hash、optimizer、step、acceptance、best/last、runtime 和搜索统计。
- `global_view_gp` 已作为普通 registry stage 接入 JSON pipeline；`nsgp lab` 与该 stage 复用同一搜索实现。
- retention 测试覆盖正常三文件、显式 stage save、外部 hash 和 workspace cleanup；数值测试覆盖 placement round trip。
- pipeline 和 lab 在搜索异常时也保留三文件 metrics-only 失败摘要；无已完成阶段时允许只有 CSV 表头。
- CMake 构建会记录 Git remote/branch/commit 和配置时 tracked-worktree dirty 状态；pipeline 可透传顶层 `experiment` 研究元数据。
- 新增共享 step-policy 原语 `ea::StepController`（constant/cosine/trust）与探索型 stage
  `adaptive_lambda_gp`（funnel lambda + exact funnel acceptance + 内存内 best-feasible restore）；
  契约测试 `nsgp_adaptive_lambda_contracts` 覆盖 step policy、lambda controller 合成轨迹、
  acceptance 四象限、restore 与 7×3 optimizer×step smoke。
- 本机无 CMake 时可用 `scripts/build_manual_gcc.ps1` 手工构建（目标与 CMakeLists 一一对应）。

当前注册模块：

```text
layout_init
hpwl_adam
exact_joint_gp
exact_recovery
surplus_bisection
equal_shape_swap
 global_capacity_transport
 density_coordinate
 global_view_gp
adaptive_lambda_gp
```

历史 DCT/Poisson 和 historical exact replay 保持隔离。

## 2. 权威 `adaptec1` 输入 checkpoint

```text
D:\codex_project\HUAWEI_EDA\epsilon-active\experiments\h375_a1_surplus_recovery10_exchange\run\global.pl
SHA-256: ba0fd17e9111b5bcbd13ed746ab9857106d962652552d26076c6e82007551e6b
```

canonical evaluator：512×512，target density 1.0。

| HPWL | Overflow |
|---:|---:|
| 85,999,318.311947 | 6.99999030053% |

文件只读使用，不复制、不清理、不提交。

## 3. 历史 V3 实验

| Run | Iterations | Best feasible HPWL | Best overflow | 观察 |
|---|---:|---:|---:|---|
| `v3_e0_adam_control_20260907` | 50 | 85,999,318.311947 | 6.99999030053% | 零接受步 |
| `v3_e5_active_bundle_20260907` | 50 | 85,999,318.311947 | 6.99999030053% | 零接受步 |
| `v3_e5_active_bundle_backtrack_20260907` | 50 | 85,999,318.311947 | 6.99999030053% | exact backtracking 后仍零接受 |
| `v3_e6_capacity_active_bundle_20260907` | 50 | 86,390,926.507477 | 7.0000002998563% | transport 24 move，未达到 6% |

历史 raw-to-H375 本地复现：HPWL `88,664,604.96`，overflow `6.99999968364%`；它与外部 checkpoint 不是同一个结果。

## 3.1 V4 adaptive-lambda 实验（同一 checkpoint、512×512/1.0、threads=1、seed 219）

A 系列（50 轮；A0 为未改动的 `global_view_gp` strict-cap 对照）：

| Run | Accepted/Rejected | best-feasible HPWL | best-any（overflow） | 结论 |
|---|---|---:|---:|---|
| `a0_gv_control50` | 12/38 | 85,999,314.800417（7.0000%） | 同左 | strict cap 只能沿容差边界滑动 −3.51 |
| `a1_alambda_dynamic` | 30/20 | 85,999,318.311945（=输入） | 85,878,982.41（8.228%） | 探索 −0.14%，50 轮内未回收 |
| `a2_alambda_noexplore` | 0/50 | 85,999,318.311945（=输入） | 同左 | corridor 关闭后完全零接受 |
| `a3_alambda_fixed` | 31/19 | 85,999,318.311945（=输入） | 85,879,071.25（8.114%） | 与 A1 几乎相同 |

B 系列（20 轮 screen，7 optimizer 全部通过基本健康检查）：dual-averaging 探索最强
（best-any 85,835,773 @ 8.472%，且是唯一在运行中把 λ 推高到 0.296 并尝试回收的）；
normalized-sgd 几乎不动（7.001%）。无 NaN、无 hard-ceiling 卡死。

C 系列（50 轮，4 optimizer × lr {0.005, 0.02, 0.08}）：更高 lr 探索更深
（dual-averaging/lr0.08 best-any 85,084,036 @ 13.186%，λ 到 10.79），但 50 轮内全部未能回到
≤7%。唯一一次 recovery-mode 接受出现在 sgd/lr0.08（11.42%→11.08%）。

D 系列（100 轮，adam 与 dual-averaging × 3 step policy，lr=0.08）：

| Run | best-feasible HPWL | best-feasible overflow | best-any | first return iter |
|---|---:|---:|---:|---:|
| `d_dual-averaging_constant` | **85,445,335.249799** | **6.8253943126%** | 84,995,882.11（12.653%） | 95 |
| `d_dual-averaging_trust` | 85,449,211.827444 | 6.8499194432% | 85,007,987.24（12.661%） | 95 |
| `d_dual-averaging_cosine` | 85,999,318.311945（=输入） | — | 85,067,886.80（13.417%） | 未回收 |
| `d_adam_constant` / `d_adam_trust` | 85,999,318.311945（=输入） | — | 85,424,004.74（12.306%） | 未回收 |
| `d_adam_cosine` | 85,999,318.311945（=输入） | — | 85,487,327.83（12.571%） | 未回收 |

E 系列（200 轮）：

| Run | final selected HPWL | overflow | best-any（overflow） | first return iter |
|---|---:|---:|---:|---:|
| `e_dual-averaging_constant` | **85,698,886.84451** | 6.8921080706% | 84,991,472.64（12.036%） | 154 |
| `e_dual-averaging_trust` | 85,718,118.441397 | 6.8962279280% | 85,008,190.20（12.021%） | 155 |
| `e_adam_trust` | 85,999,318.311945（=输入） | 6.9999903005% | 85,627,988.26（10.892%） | 未回收 |

关键数值结论（当前 checkpoint/evaluator/budget 下）：

1. 探索型 funnel 相比 strict cap 有 material feasible 改善：100 轮 dual-averaging/constant 达到
   best-feasible 85,445,335.249799 @ 6.8254%，相对输入 −553,983.06（−0.644%）。
2. 100 轮的走廊调度优于 200 轮：更深的探索（贴 15% 上限、>10% 停留 109 轮）换来的回归点更差
   （−0.349% vs −0.644%）。
3. Recovery 成功只发生在 dual-averaging 上；其 dual average 在 λ 增大后自发产生 overflow 下降方向。
   Adam 系探索正常但 recovery 全部失败；cosine step 使 recovery 窗口内步长过小，未回收。
4. 结构性限制：density direction 最小化 energy（平方超额），不保证降低 overflow_ratio；恢复期提案
   常表现为 energy 下降而 overflow 不降，从而被 exact_funnel 拒绝。

## 3.2 V5 h221 续跑实验（trajectory-lambda 纯下降，无逐步接受机制）

起点改为历史链中间 checkpoint `h221_dct_poisson_512/best.pl`（DCT/Poisson 同伦 512 精化，历史平滑
代理；SHA-256 `9d0bd999...`）。canonical 加载 + clamp 后 exact audit 为 HPWL 109,837,622.009 /
overflow 7.8193%（与 legacy 链一致）。注意：该 .pl 有 480 个可动节点共 23.3% 面积在 die 外；
无 clamp 的 `nsgp audit` CLI 会得到 11.52% 的伪影值，不是 pipeline 语义。

协议：`exact_joint_gp`，batch acceptance 关闭（每步无条件移动），lambda 用 `trajectory` 策略
（overflow 相对 smoothstep 目标轨迹的 PI 反馈，interval 5，horizon 100，stop 0.07），lr 0.002，
100 轮，threads=1。7 optimizer 对照（run id 前缀 `20260908_0100_adaptec1_h221traj_`）：

| Optimizer | best-feasible HPWL | overflow | 首次 ≤7% | λ 峰值 |
|---|---:|---:|---:|---:|
| amsgrad | **99,136,915.82** | 5.4217% | 54 | 51 |
| adam | 99,960,810.50 | 5.5886% | 56 | 51 |
| adagrad | 101,850,109.32 | 6.9879% | 75 | 130 |
| dual-averaging | 105,949,370.15 | 6.9907% | 96 | 131 |
| sgd | 108,191,594.03 | 5.2690% | 10 | 3.2 |
| heavy-ball | 108,228,173.46 | 5.2480% | 11 | 3.3 |
| normalized-sgd | 108,985,009.69 | 4.9752% | 21 | 3.4 |

结论：overflow 反馈 lambda 的纯下降机制可行——无任何逐步接受/回溯，7 个 optimizer 全部回到 ≤7%，
HPWL 改善 −1.6% 到 −9.7%；逐坐标自适应步长（amsgrad/adam）在该机制下显著占优，且 amsgrad 在
第 99 轮仍在下降。原始子梯度步长几乎不动布局（λ 无需超过 3.4）。曲线与数据在
`framework/results/analysis/20260908_h221traj_optimizer_screen/`（本地，不入库）。

### 3.3 V6 h221→86M 调参（exact_joint_gp 连续次梯度，300 轮）

目标 86M/<7% **未命中**。最佳 run `C_B5_md025`：AMSGrad(lr .002, md .25, β .9/.99) +
ratio λ(ob=.085, itv=2, ds=.5) → best-feasible **94,061,868.14 @ 6.7474%**，相对起点
−14.36%，相对 V5 100 轮基线 −5.12%。关键发现：

1. ratio λ 在 interval=5 时早期被无条件 1.25 增长项主导（B0/B2 逐位一致）；
   interval=2 才真正激活（B5 95.83M，比 traj 系最好好 240 万）。
2. maximum_delta（步长裁剪 = md×lr）单调有效：4→0.25 把首轮 overflow 冲击从
   22.4% 压到 13.8%，HPWL 95.83M→94.06M；md=0.125 反转（94.19M）。
3. trajectory λ 把 overflow 压到 4.7-4.8%（λ 峰值 81-412），容量裕量没有转化为
   HPWL，距 86M 更远——与 V4 的结论一致：direction 质量而非 λ 强度是瓶颈。
4. 动量微调全部劣于 (0.90, 0.99)；adagrad/dual-averaging 远差（后者 λ 发散）。
5. A0 前 100 轮复现 V5 基线到打印精度（first<7% 同为 iter 54）。

证据与逐轮标量：`framework/results/analysis/20260908_h221_86m_tuning/`
（runs.csv、per_iteration_metrics.csv、三张收敛图，本地不入库）；台账回填在
`plan/NONSMOOTH_GP_H221_86M_TUNING_PLAN_V6.md` §10。30 个实验目录全部通过
retention 契约，checkpoint SHA-256 全程未变。E 阶段触发条件未满足（最后 50 轮
仍下降 0.647% > 0.2%），未实现 schedule；证据支持的下一步是最佳配置延长到 500 轮。

## 4. 本次微内核回归证据

### 4.1 Non-local oracle implementation health smoke（H375 canonical input）

2026-09-08 已实现并注册五个独立 pure-descent oracle，使用 checkpoint SHA
`ba0fd17e...07551e6b`、512×512/1.0、threads=1、AMSGrad/lr .002、一轮。每个 run 都是 metrics-only，
无接受/拒绝/回滚，输入 hash 未被写入。末态仅用于健康检查：multiscale HPWL 89,408,103.56 / overflow
22.7080%；cut 88,791,145.10 / 22.8037%；transport 88,777,884.40 / 23.5175%；charge 88,431,118.63 /
28.4872%；finite-radius 88,924,114.36 / 23.6278%。这些大幅 excursion 说明无条件更新路径确实生效，
但也说明默认量纲/步幅尚未校准；它们不是 20/100 轮比较结果。

The common 100-round AMSGrad primary completed on the same input. Final (HPWL / overflow) was multiscale
139,548,754.76 / 14.8103%; cut 156,286,419.76 / 69.0618%; transport 224,641,758.21 / 24.2072%; charge
87,801,002.56 / 87.1644%; finite-radius (4 probes) 195,518,322.69 / 13.8906%. This is a negative result:
no initial oracle returned to 7%, therefore none proceeds to optimizer/LR or fusion ranking without scale repair.

- 模块权威构建 smoke：`layout_init → hpwl_adam → exact_joint_gp` 可运行。
- 默认 `run` 结果目录只有 `params.json`、`experiment.md`、`trajectory.csv`，仓库根无 `global_metrics.csv`/`last.pl` 泄漏。
- 显式 `--save-stage layout_init` 只增加 `saved/layout_init.pl`。
- `lab` 1 轮、2 threads：HPWL `85,999,279.708727`，overflow `6.9999980678266%`，1 accepted；输入 SHA-256 前后不变。
- invalid optimizer 异常 smoke：exit 1，`%TEMP%/nonsmooth-gp/<run_id>` 不残留。
- C++ numeric/placement round-trip contracts 通过。
- `global_view_gp` registry/lab 同参数回归（1 轮、2 threads、active ensemble + bundle）末态均为 HPWL `85,999,318.311948`、overflow `6.99999030053%`，该步被 exact 接受规则拒绝；两份结果均通过三文件 retention 检查。
- 未知模块失败 smoke：exit 1，仍生成 `params.json`、带百分比表头的 `trajectory.csv` 和 `status: failed` 的 `experiment.md`；retention 检查通过。

这些是结构和回归 smoke，不替代新的 50 轮公平算法对照。

## 5. 剩余缺口

0. 新 non-local modules 尚缺 planned per-iteration trajectory telemetry、synthetic module contracts、
20-round screen、100-round control/primary comparison和后续 optimizer/LR/fusion study；默认一轮仅是
结构 health smoke，不能用于效果排名。

1. `nsgp_numeric_kernel` 为减少链接复杂度仍同时编译共享 kernel 和模块算法源；物理所有权已正确。不要仅为目录形式改成多层 object libraries，除非测得编译收益。
2. Git 元数据在 CMake configure 时记录，因此 build 后又修改源码时不会自动刷新；手工 g++ 构建显示 `unknown`/`not_checked`。
3. SHA-256 使用 Windows CryptoAPI；Linux 尚无实现。
4. 依赖外部 benchmark 的 retention lifecycle smoke 未注册进 CTest；当前由 PowerShell 测试显式运行，避免硬编码本机数据路径。
5. global-view 仍使用固定 lambda，bundle reset 和 capacity transport 搜索尚未覆盖 V3 计划全部变体。
6. `reference_nonsmooth_chain.json` 使用 64×64/0.9，`smoke.json` 使用 32×32/0.9；它们不是 canonical 512×512/1.0 正式实验。`global_view_smoke.json` 使用 canonical evaluator，但只有 1 轮，仍不是正式对照。
7. `adaptive_lambda_gp` 的 recovery 依赖 lambda 长时间增大（D 系列成功案例在 it95 才首次回到 ≤7%），
   50 轮预算内 recovery 从未成功；constant 与 trust step 在 D 系列中等价。

## 6. 推荐下一步

V4 已回答 Q1（探索 + funnel 收益为正且 material）、Q2（50 轮内收益主要来自 corridor；dynamic lambda
是 recovery 的必要压力源但非探索收益来源）、Q3（dual-averaging + constant/trust 最优）。后续自然延伸：

1. 把 recovery 失败的根因（energy direction ≠ overflow direction）变成可检验改动：为
   `adaptive_lambda_gp` 增加 overflow 导出方向或 bin-price 方向，先 smoke 再做 100 轮对照。
2. 组合链验证：global_capacity_transport → adaptive_lambda_gp → exact_recovery/equal_shape_swap。
3. 若继续调参：在 100 轮预算附近做 hold/contract fraction 的两点消融（如 0.30/0.70），不要一次扫大矩阵。

## 7. 不能据此宣称

- 不能宣称 active ensemble/bundle 已优于 baseline。
- 不能宣称 capacity transport 达到 6%。
- 不能把 7% overflow 称为 legalized。
- 不能把 smoke 中的一步改善当作稳定统计结论。
- 不能把深探索处的 best-any HPWL 当作结果；正式结论只用 best-feasible。
- 不能宣称 dual-averaging 对所有 case/预算通用；证据仅限 adaptec1、该 checkpoint 与 20–200 轮预算。
