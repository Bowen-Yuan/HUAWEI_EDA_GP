# 05 — 当前状态、数值结果与已知缺口

记录日期：2026-09-07（Asia/Shanghai）。架构：registry-driven microkernel。接手时以 `git rev-parse HEAD` 为准。

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
- retention 测试覆盖正常三文件、显式 stage save、外部 hash 和 workspace cleanup；数值测试覆盖 placement round trip。

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
global_view_gp       # 通过 lab 单阶段入口
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

## 4. 本次微内核回归证据

- 模块权威构建 smoke：`layout_init → hpwl_adam → exact_joint_gp` 可运行。
- 默认 `run` 结果目录只有 `params.json`、`experiment.md`、`trajectory.csv`，仓库根无 `global_metrics.csv`/`last.pl` 泄漏。
- 显式 `--save-stage layout_init` 只增加 `saved/layout_init.pl`。
- `lab` 1 轮、2 threads：HPWL `85,999,279.708727`，overflow `6.9999980678266%`，1 accepted；输入 SHA-256 前后不变。
- invalid optimizer 异常 smoke：exit 1，`%TEMP%/nonsmooth-gp/<run_id>` 不残留。
- C++ numeric/placement round-trip contracts 通过。

这些是结构和回归 smoke，不替代新的 50 轮公平算法对照。

## 5. 剩余缺口

1. `global_view_gp` 已物理归入模块目录，但仍通过专用 `lab` CLI，而不是普通 JSON pipeline 的 `StageFunction`；若需要与其他阶段任意重排，应提取为同一 registry adapter。
2. `nsgp_numeric_kernel` 为减少链接复杂度仍同时编译共享 kernel 和模块算法源；物理所有权已正确，但可进一步使用小型 object libraries，前提是确有编译收益。
3. Git branch/commit 在 CMake configure 时记录；手工 g++ 构建显示 `unknown`。dirty state 仍为 `not_checked`。
4. SHA-256 使用 Windows CryptoAPI；Linux 尚无实现。
5. failure experiment 当前保证临时 workspace 清理，但失败后的持久结果摘要还不是统一的三文件失败报告。
6. retention lifecycle smoke 尚未全部注册进 CTest/CTest fixture。
7. global-view 仍使用固定 lambda，bundle reset 和 capacity transport 搜索尚未覆盖 V3 计划全部变体。
8. `reference_nonsmooth_chain.json` 使用 64×64/0.9，`smoke.json` 使用 32×32/0.9；它们不是 canonical 512×512/1.0 正式实验。

## 6. 推荐下一步

除非用户另有任务，优先把 `global_view_gp` 提取成普通 registry stage，然后用 JSON pipeline 表达：

```text
external H375 checkpoint
→ global_capacity_transport
→ global_view_gp
→ exact_recovery / equal_shape_swap（可选）
```

完成后进行相同 checkpoint/hash、threads、seed、预算的 50 轮对照。不要为了“架构完整”再增加服务层、插件框架或复杂 artifact manager。

## 7. 不能据此宣称

- 不能宣称 active ensemble/bundle 已优于 baseline。
- 不能宣称 capacity transport 达到 6%。
- 不能把 7% overflow 称为 legalized。
- 不能把 smoke 中的一步改善当作稳定统计结论。
