# 02 — 微内核架构与代码地图

## 1. 真实构建与调用图

```text
Bookshelf / external checkpoint
                │
                ▼
framework/kernel
I/O + ea::Database + exact HPWL/density + optimizer
                │
                ▼
framework/code/microkernel
registry + stage context + fresh audit + three-file experiment log
                │
       pipeline module sequence
                │
      ┌─────────┼────────────────────────────┐
      ▼         ▼                            ▼
modules/layout_init   modules/exact_joint_gp   modules/.../module.cpp
      │         │                            │
      └─────────┴──── in-memory Database ────┘
                │
                ▼
params.json + experiment.md + trajectory.csv
```

顶层没有算法 `if/else` 链。`main.cpp` 只解析命令和 pipeline，通过 `ModuleRegistry` 按名字调用阶段。新增普通优化阶段时，主要改动应局限于：

```text
modules/<new_stage>/code/module.cpp
modules/<new_stage>/params/*.json
CMakeLists.txt                 # 加入一个 source
microkernel.cpp               # 注册一行
```

如果新阶段需要新的共享 exact 原语，才修改 `framework/kernel`。

## 2. framework 的最小职责

### `framework/kernel`

共享数值内核，不包含固定优化流程：

- `include/epsilon_active/*.hpp`：`ea::Database`、exact evaluator 和算法接口。
- `src/bookshelf.cpp`：Bookshelf I/O 和坐标转换。
- `src/hpwl.cpp`：exact weighted pin-offset HPWL 与 active direction。
- `src/density.cpp`：exact rectangle/bin occupancy、overflow、energy 和 incremental audit。
- `src/optimizer.cpp`：统一 optimizer 实现。
- `src/step_policy.cpp`：轻量 step controller（constant / cosine / trust），只决定 learning rate 和
  per-step maximum delta；acceptance 和 backtracking 仍在调用方，保证所有 optimizer 对比共享同一个
  exact 接受契约。

### `framework/code/microkernel.hpp/.cpp`

微内核只提供：

- `DensityConfig`、`ExactMetrics`、`StageContext`、`StageStats`；
- `ModuleRegistry`；
- `exact_audit`、movable clamp、OpenMP thread 设置；
- checkpoint SHA-256；
- `ExperimentLog` 三文件记录器；
- 默认模块的最小注册表。

### `framework/code/main.cpp`

- 命令：`run`、`batch`、`audit`、`list-modules`、`lab`。
- 读取 JSON pipeline 和每阶段 JSON 参数。
- 构造完整 module chain provenance。
- 在每个阶段前后 fresh exact audit。
- 默认只产生 `params.json`、`experiment.md`、`trajectory.csv`。
- 只有 `--save-stage MODULE` 才写 `saved/MODULE.pl`。

`global_view_gp` 同时提供普通 registry stage 和薄的单阶段 `lab` 便捷入口；两者复用同一个搜索实现，代码都位于模块目录，不在 framework 放置算法。

## 3. 模块目录是策略权威实现

| 模块 | 权威实现/适配器 | 作用 |
|---|---|---|
| `layout_init` | `code/module.cpp` | raw 或中心高斯初始化 |
| `hpwl_adam` | `code/module.cpp` | 配置共享 exact HPWL/optimizer，执行 HPWL-only seed |
| `exact_joint_gp` | `placer.cpp`、`lambda_controller.cpp`、`batch_acceptance.cpp`、`module.cpp` | exact joint GP 阶段 |
| `exact_recovery` | `recovery.cpp`、`compact_recovery.cpp`、`module.cpp` | cap 内 recovery |
| `surplus_bisection` | `bisection.cpp`、`module.cpp` | surplus-only 容量二分 |
| `equal_shape_swap` | `swap_recovery.cpp`、`module.cpp` | equal-shape / net-aware swap |
| `global_capacity_transport` | `coarse_flow.cpp`、`transport.cpp`、`module.cpp` | coarse 候选与 exact-audited transport |
| `density_coordinate` | `density_coordinate.cpp`、`module.cpp` | exact density coordinate search |
| `global_view_gp` | `global_view_lab.cpp/.hpp` | 普通 registry stage；active ensemble、bundle、backtracking/trust；兼容薄 `lab` 入口 |
| `adaptive_lambda_gp` | `adaptive_lambda_gp.cpp/.hpp`、`module.cpp` | 探索型 dynamic-lambda GP；overflow funnel + log-domain PI(D) lambda + exact funnel acceptance + 内存内 best-feasible restore |
| `historical_exact_replay` | `legacy_stage_main.cpp` | 历史 exact CLI 兼容入口 |
| `historical_dct_poisson` | 独立 Bookshelf/electric/spectral/homotopy 代码 | 只做历史 smooth 复现 |

旧的 `framework/legacy` 权威目录已移除。原来与 module 逐字节相同的七份策略副本已删除；原 `hpwl_adam` 下未参与构建的 HPWL/optimizer 镜像也已删除。Git 历史仍可恢复，但不要重新引入源码镜像。

## 4. CMake target

- `nsgp_numeric_kernel`：四个 framework 数值原语 + 各模块算法源；名称表示共享链接单元，不表示固定算法流程。
- `nsgp`：微内核 CLI、模块适配器和复用同一 stage 实现的 global-view lab 入口。
- `nsgp_tests`：exact 数值和 placement round-trip 契约。
- `nsgp_adaptive_tests`：step policy、funnel lambda controller、funnel acceptance、best-feasible
  restore 和 optimizer × step policy 组合的 synthetic 契约测试（不需要 benchmark 数据）。
- `nsgp_legacy_stage`：历史 exact 兼容目标，源码在 module 目录。
- `nsgp_historical_homotopy`：历史 DCT/Poisson 目标。

CMake 配置时将 Git remote/branch/commit 和当时的 tracked-worktree dirty 状态编译进 `nsgp` 的实验元数据；手工编译无法提供的字段明确记录为 `unknown`/`not_checked`。

## 5. 配置与结果

- `framework/params/defaults.json`：默认数据目录、seed、thread、canonical evaluator。
- `framework/params/cases.json`：默认 case 和八 case 清单。
- `framework/params/pipelines/*.json`：模块顺序；不是硬编码算法。
- `modules/<stage>/params/*.json`：阶段参数和实验示例。
- `framework/results/experiments/<run_id>`：统一实验目录，默认被 Git 忽略。
- `framework/experiment_logs`：历史重产物复现日志，不作为新实验模板。

## 6. 状态与依赖方向

```text
framework/kernel  ← modules  ← microkernel runner
```

- 模块之间通过同一个内存 `ea::Database` 顺序交接。
- 每个阶段自行构造短生命周期 evaluator/optimizer；不跨阶段偷偷共享 cache/moment。
- framework 不 include 某个具体策略实现。
- module 可以调用 kernel API；普通 module 不应依赖另一个 module 的私有实现。
- `exact_joint_gp` 链接其他低层算法函数是历史组合能力；新 pipeline 应优先把它们作为独立 stage 显式调用。

## 7. 数据和坐标

- benchmark base：`<dataset>/<case>/<case>`。
- 内存 `Node.x/y` 为中心，`.pl` 为左下角。
- fixed node 不移动；physical fixed macro 占容量；`terminal_NI` 不占容量。
- 数据集和外部 checkpoint 不进入 Git。
