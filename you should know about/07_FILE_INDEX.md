# 07 — 文件级索引

## 根目录

- `.gitignore`：忽略 build、exe、常规 results/logs 和临时 `.pl`。
- `AGENTS.md`：强制 AI 阅读并同步交接知识库。
- `CMakeLists.txt`：权威构建与源码所有权清单。
- `README.md`：简明构建、CLI 和微内核说明。
- `third_party/json.hpp`：vendored nlohmann JSON；不要手改内部实现。

## `framework/code`

- `main.cpp`：极薄 CLI/pipeline orchestration；不实现具体优化算法。
- `microkernel.hpp`：阶段上下文、统计、registry、audit、experiment log API。
- `microkernel.cpp`：上述最小机制的实现与默认模块注册。

原来未接线的 `common/problem/layout/metrics/module_api` 骨架已删除，避免第二套运行状态。

## `framework/kernel`

- `include/epsilon_active/types.hpp`：唯一运行时 `ea::Database` 和基础指标。
- `include/epsilon_active/*.hpp`：exact evaluator、optimizer 与各模块算法接口。
- `src/bookshelf.cpp`：Bookshelf I/O。
- `src/hpwl.cpp`：exact HPWL 与 active direction。
- `src/density.cpp`：exact density/overflow、node/group incremental audit。
- `src/optimizer.cpp`：Adam、AMSGrad、AdaGrad、HeavyBall、SGD、NormalizedSGD、DualAveraging。

`framework/legacy` 已移除，不再是构建权威。

## `modules`

每个普通 stage 以 `code/module.cpp` 连接微内核，以 `params/*.json` 表达实验变量：

- `layout_init`：初始化适配器和 raw/center-Gaussian 参数。
- `hpwl_adam`：HPWL-only stage 适配器；HPWL/optimizer 直接复用 kernel，没有源码镜像。
- `exact_joint_gp`：权威 placer、lambda、batch acceptance、stage adapter 和参数。
- `exact_recovery`：权威 recovery、compact recovery、stage adapter 和 cap 参数。
- `surplus_bisection`：权威 bisection、stage adapter 和 H372-like 参数。
- `equal_shape_swap`：权威 swap recovery、stage adapter 和 H375-like 参数。
- `global_capacity_transport`：权威 coarse flow、transport、stage adapter 和基础参数。
- `density_coordinate`：权威 coordinate search、stage adapter 和基础参数。
- `global_view_gp`：`global_view_lab.cpp/.hpp`、control/active-bundle 参数；目前为专用 lab stage。
- `historical_exact_replay`：历史 exact CLI 源，保留 target 名 `nsgp_legacy_stage` 以兼容脚本。
- `historical_dct_poisson`：独立 Bookshelf/electric/spectral/homotopy 历史 smooth 模块。

各 `results/.gitkeep` 只保留空目录，不在模块目录堆积普通实验产物。

## `framework/params`

- `defaults.json`：数据目录、threads、seed、canonical 512×512/1.0。
- `cases.json`：默认与完整 case 清单。
- `pipelines/smoke.json`：32×32 三阶段结构 smoke。
- `pipelines/reference_nonsmooth_chain.json`：64×64 六模块组合示例。

V3 global-view 示例参数已经归入 `modules/global_view_gp/params`，不再在 framework 复制一份。

## `scripts`

- `run_h375_replay.ps1`：历史 H219→H375 重产物复现脚本；不作为新实验模板。

## `tests`

- `test_numeric_contracts.cpp`：exact HPWL/density 小例及 placement round trip。
- `test_v3_retention_contract.ps1`：三文件、百分比字段、显式 save、外部 hash、workspace cleanup 检查。

## `plan`

- `NEW_PROJECT_IMPLEMENTATION_PLAN_V2.md`：项目数学契约、模块化目标和迁移设计。
- `NONSMOOTH_GP_RETENTION_GLOBAL_OPTIMIZER_PLAN_V3.md`：metrics-only、global-view、optimizer/step 和 capacity 实验方案。

计划是审查依据，不是完成证明；以当前 CMake、源码、测试和 exact 实验为事实。
