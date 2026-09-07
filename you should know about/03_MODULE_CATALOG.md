# 03 — 模块目录

## 模块调用关系

参考纯非光滑组合链：

```text
layout_init
  → hpwl_adam
  → exact_joint_gp
  → surplus_bisection
  → exact_recovery
  → equal_shape_swap
```

这只是可组合示例。模块可以重复、删除和重排，但每个边界都要 fresh exact audit，并记录输入来源和状态重置策略。

历史 H375 复现链：

```text
raw Bookshelf .pl
  → H219 HPWL-only Adam seed
  → H219 128×128 DCT/Poisson homotopy       [historical smooth]
  → H221 512×512 refine                     [historical smooth]
  → H252–H257 exact recovery to 7% band
  → H372 surplus-only recursive bisection
  → H372 exact recovery ×10
  → H375 equal-shape net-aware swap ×5
```

该完整链由 `scripts/run_h375_replay.ps1` 描述；不能把其中 DCT/Poisson 阶段重新命名后放进正式 challenge pipeline。

## `layout_init`

- 职责：创建起始布局；保留 fixed nodes。
- 模式：`raw` 或 `center_gaussian`。
- 参数：`modules/layout_init/params/raw.json`、`center_gaussian_219.json`。
- 当前执行位置：runner 只负责 `read_bookshelf`；`modules/layout_init/code/module.cpp` 处理 raw 或 `initialize_center_gaussian` 策略。
- 输出：内存中的 `ea::Database` 布局；只有显式 `--save-stage layout_init` 才保留 placement。

## `hpwl_adam`

- 职责：使用 epsilon-active exact HPWL subgradient 产生 seed；density weight 为零。
- 阶段适配器：`modules/hpwl_adam/code/module.cpp`；共享 exact HPWL/optimizer 位于 `framework/kernel`。
- 参数：`seed_h219_like.json`、`smoke.json`。
- 重要边界：epsilon/active power 只影响方向，不改变 exact HPWL 报告值。

## `exact_joint_gp`

- 职责：组合 exact HPWL 方向和 exact overlap 方向，以 lambda controller 和 optimizer 推进全局布局。
- 权威代码：`modules/exact_joint_gp/code/placer.cpp`、`lambda_controller.cpp`、`batch_acceptance.cpp`、`module.cpp`。
- 参数：`pure_coarse.json`、`bridge_h253_like.json`、`retighten_h254_like.json`、`smoke.json`。
- 可选能力：net-batch direction、density direction sharing、regional price、late-stage switch、bisection/coarse flow/transport/recovery/swap 后处理。
- `global_place` 只有显式非空 output directory 时才写内部指标/snapshot；微内核模块传空路径，正常 pipeline 不产生这些副作用。

## `global_view_gp`

- 位置：`modules/global_view_gp/code/global_view_lab.cpp`；已注册为普通 `StageFunction`，可在 JSON pipeline 中重复、删除、重排。
- 入口：组合实验使用 `nsgp run`；`nsgp lab` 只是复用同一搜索函数的单阶段便捷入口，不是第二套算法。
- 输入：外部 checkpoint、512×512/target 1.0 默认 evaluator、optimizer/step 参数。
- active ensemble：epsilon scale `{0, 0.25, 0.5, 1.0} × bin_size`。
- 聚合：各方向 RMS normalization 后做 simplex 上的 minimum-norm convex combination。
- temporal bundle：最多四个方向，bundle mix 默认 0.5；连续 reject 会清空历史并重置 optimizer。
- 接受：canonical exact audit、严格 overflow cap、HPWL 必须下降；最多九次二分回溯。
- 输出：best feasible 指标；布局默认不持久化。
- 1 轮同参数回归中，registry stage 与 lab 的末态 exact 指标一致到打印精度；历史 50 轮实验仍显示严格 7% 边界下 E0/E5 零接受步，详见当前状态文档。

## `adaptive_lambda_gp`

- 职责：从 canonical feasible checkpoint 出发，允许搜索临时进入 >7% overflow 区域换取 HPWL 探索，
  再用动态 lambda 与逐渐收紧的 overflow funnel 把搜索拉回 7%。
- 权威代码：`modules/adaptive_lambda_gp/code/adaptive_lambda_gp.hpp`（纯策略逻辑，header-only，可被
  synthetic 测试直接包含）、`adaptive_lambda_gp.cpp`（布局搜索）、`module.cpp`（registry 适配器与 JSON 解析）。
- 方向：wire（epsilon-active ensemble，各 RMS normalize 后 simplex 最小范数组合）与 exact overlap
  density direction 分别 RMS normalize，组合为 normalize(wire + lambda * density)，因此 lambda 是
  无量纲方向权重。
- Lambda：模块私有 `FunnelLambdaController`（不改动 `exact_joint_gp` 的 controller）。corridor 由
  explore 高度 smoothstep 收缩到 final（默认 hold 0.20、lock 0.75）；每 `update_interval` 轮做
  log-domain PI(D) 更新（百分比尺度误差、积分 clamp ±4、log-step clamp [-0.35, 0.55]）；final lock
  且 overflow > final 时额外乘 exp(clamp(0.10·误差, 0, 0.25))。`name: "fixed"` 表示 fixed-lambda 消融。
- Acceptance：`exact_funnel`——candidate 先做 canonical fresh exact audit；overflow > hard ceiling
  （默认 16%）必拒；current 在 corridor 内为探索模式（candidate ≤ corridor 且 HPWL 严格下降），否则为
  恢复模式（overflow 至少改善 0.002%、单步 HPWL 代价 ≤ 0.05%、累计 ≤ 初始 HPWL 的 0.30%）。所有
  optimizer/step policy 共用同一 9 次二分 exact backtracking。
- Best-feasible 保底：iteration 0 的 checkpoint 即为内存内 fallback；任何 canonical feasible
  （≤ final + 1e-9）且 HPWL 更好的状态都保存 movable 坐标；stage 结束无条件 restore + clamp + fresh
  audit。若输入本身在 final 目标下不可行，stage 显式抛错而不是输出 >7% 布局。
- Optimizer reset：reject_streak ≥ 4、lambda 比值 ≥ 4、进入 final lock 各触发一次；不重置 lambda
  controller。JSON 可选 7 种 optimizer 与 3 种 step policy（constant/cosine/trust）。
- 产物：metrics-only 三文件；stage 结束向 stdout 打印 initial/best-any/best-feasible/
  last-before-restore/final-selected、maximum excursion、iterations above 10%、first return iteration、
  optimizer resets 与 lambda 初值/峰值/末值；`verbose: true` 时附加逐轮 key=value 遥测（仅 stdout，
  不落盘）。输入 checkpoint 必须在 final target 下 feasible。

## `global_capacity_transport`（V3 lab 内部能力）

- 轻量 V3 版本位于 `modules/global_view_gp/code/global_view_lab.cpp` 的 `capacity_transport`。
- 候选：64×64 coarse occupancy 找最大 surplus source 和最低 occupancy destination。
- 限制：跳过宽/高超过 coarse bin 的宏；最多 24 pass。
- 接受：每个候选都在 canonical 512×512 exact evaluator 上重新审计；overflow 必须下降，HPWL 增幅不得超过 0.2%。
- 参数：`--capacity-transport --transport-target-percent <percent>`。
- 完整独立 pipeline stage 为 `global_capacity_transport`，代码在同名模块，支持 nearest/auction/Hilbert 和更复杂原子组策略。

## `exact_recovery`

- 职责：在 overflow cap 内通过 node move、net block、breakpoint、compact direction 等搜索 exact HPWL 改善。
- 权威代码：`modules/exact_recovery/code/recovery.cpp`、`compact_recovery.cpp`、`module.cpp`。
- 参数：`cap15_netblock.json`、`cap07_node.json`、`cap07_full.json`。
- 接口：`recover_hpwl_under_overflow(Database&, ExactOverlapDensity&, RecoveryConfig)`。
- 接受与统计必须使用 exact affected HPWL 和 exact density move/group move audit。

## `surplus_bisection`

- 职责：根据 exact capacity 将超额区域递归二分，并用位置、net 和容量启发式分配节点。
- 权威代码：`modules/surplus_bisection/code/bisection.cpp`、`module.cpp`。
- 参数：`h372_like.json`。
- 常用选项：`surplus_only`、position seeded、nearest-capacity leaf、HPWL-guided leaf。
- 它是候选布局模块，边界后仍须 canonical fresh audit。

## `equal_shape_swap`

- 职责：在同形单元之间做 occupancy-invariant 或 exact-audited 交换，使用空间与 net-aware 候选降低 HPWL。
- 权威代码：`modules/equal_shape_swap/code/swap_recovery.cpp`、`module.cpp`。
- 参数：`h375_like.json`。
- 扩展能力：permutation、assignment、density-guided swap 和有限 HPWL budget。
- equal-shape occupancy invariance 是重要回归测试点。

## 其他显式模块

- `global_capacity_transport`：`coarse_flow.cpp`、`transport.cpp` 与 `module.cpp`，提供 excess-to-capacity relocation、原子组、auction/Hilbert、identity exchange。
- `density_coordinate`：`density_coordinate.cpp` 与 `module.cpp`，提供 exact overlap coordinate descent。
- `compact_recovery.cpp`：support contraction。
- `batch_acceptance.cpp`：对 proposal batch 做 exact backtracking/接受。

这些能力已经进入 `ModuleRegistry`，可直接写入 pipeline。新模块同样必须新增清晰参数文件和调用链记录，不要只增加难以发现的 CLI flag。

## `historical_dct_poisson`

- 代码：`modules/historical_dct_poisson/code`。
- 独立 target：`nsgp_historical_homotopy`。
- 内容：独立 Bookshelf 数据结构/读取、FFT/DCT、electric density、homotopy main。
- 参数标记：`disabled.json` 明确 `uses_smooth_surrogate=true`。
- 用途：只复现 H219/H221 历史链。
- 禁止：进入正式 `challenge_nonsmooth` pipeline 或作为最终 evaluator。
