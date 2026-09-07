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
