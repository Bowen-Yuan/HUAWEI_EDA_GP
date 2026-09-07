# 07 — 文件级索引

本索引覆盖受管项目文件的职责，不列出 `bin/build` 二进制、被忽略的普通实验内容以及 vendored `json.hpp` 的内部实现细节。

## 根目录

- `.gitignore`：忽略 build、exe、framework/module results、experiment logs 和临时 `.pl`。
- `AGENTS.md`：让支持仓库级指令的 AI 自动获知交接文档和同步更新义务。
- `CMakeLists.txt`：唯一权威构建清单。
- `README.md`：面向使用者的简要说明。
- `third_party/json.hpp`：nlohmann JSON 单头依赖；除非升级依赖，不要手工修改。

## `framework/code`

- `main.cpp`：主 CLI、JSON pipeline runner、module 静态分派、旧式结果 writer、audit。
- `experiment_lab.cpp`：V3 metrics-only 实验、SHA-256、active/bundle、backtracking/trust、capacity transport。
- `experiment_lab.hpp`：`run_global_view_lab` 声明。
- `common.hpp`：未来 `nsgp` 数据层的基本类型。
- `problem.hpp`：未来只读问题结构。
- `layout.hpp`：未来独立布局和 revision。
- `metrics.hpp`：未来统一指标结构。
- `module_api.hpp`：未来模块上下文、状态产物和返回值协议。

## `framework/legacy/src`（当前权威 core）

- `bookshelf.cpp`：Bookshelf 文件读取、placement 读写、中心初始化。
- `hpwl.cpp`：exact weighted pin-offset HPWL 和 epsilon-active direction。
- `density.cpp`：exact rectangle/bin occupancy、energy/overflow、direction、node/group move delta。
- `optimizer.cpp`：统一 stateful optimizer 实现和字符串解析。
- `lambda_controller.cpp`：density multiplier 初始化与更新策略。
- `batch_acceptance.cpp`：批量候选的 exact backtracking/接受。
- `placer.cpp`：旧 `global_place` 主循环、direction 组合、可选后处理、snapshot 和结果写盘。
- `recovery.cpp`：overflow cap 下 node/net-block exact recovery。
- `compact_recovery.cpp`：compact support contraction。
- `density_coordinate.cpp`：基于 exact overlap 的坐标下降。
- `bisection.cpp`：容量约束递归二分和 leaf assignment。
- `coarse_flow.cpp`：coarse capacity flow 候选生成与接受。
- `transport.cpp`：容量搬运、atomic/group/auction/Hilbert/identity 策略。
- `swap_recovery.cpp`：equal-shape、net-aware、assignment/permutation swap。
- `legacy_stage_main.cpp`：历史参数面完整的独立 CLI；不进入 core library。

## `framework/legacy/include/epsilon_active`

每个 `.hpp` 与同名 `.cpp` 对应；`types.hpp` 提供当前权威 `ea::Database` 和指标结构，`placer.hpp` 汇总组合 GP 的所有子模块配置。修改配置字段时必须同步 parser、默认值、记录器和测试。

## `framework/params`

- `defaults.json`：数据目录、thread、seed、canonical 512×512/1.0 和旧 artifacts 默认。
- `cases.json`：default case 和八 case 清单。
- `pipelines/smoke.json`：32×32 三阶段 smoke。
- `pipelines/reference_nonsmooth_chain.json`：64×64 六模块参考组合。
- `experiments/v3_e0_adam_control.json`：V3 E0 命令说明。
- `experiments/v3_e5_active_bundle.json`：V3 E5 命令说明。

## `modules`

每个模块都有 `params/`；`results/.gitkeep` 只保留空目录。`code/` 当前是研究可读镜像，不是独立 CMake library：

- `layout_init`：只有 raw/center Gaussian 参数，无独立 code 文件。
- `hpwl_adam`：`hpwl.cpp`、`optimizer.cpp` 镜像；后者当前落后于 core。
- `exact_joint_gp`：`placer.cpp`、`lambda_controller.cpp`、`batch_acceptance.cpp` 镜像及四组参数。
- `exact_recovery`：`recovery.cpp`、`compact_recovery.cpp` 镜像及 7%/15% 参数。
- `surplus_bisection`：`bisection.cpp` 镜像和 H372-like 参数。
- `equal_shape_swap`：`swap_recovery.cpp` 镜像和 H375-like 参数。
- `historical_dct_poisson`：独立的 `types.h`、Bookshelf reader、FFT/DCT spectral 实现、electric density 和 homotopy CLI；只用于历史复现。

## `scripts`

- `run_h375_replay.ps1`：按 H219 → H221 → H252–H257 → H372 → H375 顺序调用两个历史 executable；每阶段写 invocation、stdout/stderr 和 placements。它是重产物历史复现脚本，不是新实验模板。

## `tests`

- `test_legacy_numeric_contracts.cpp`：2-node toy DB，验证 weighted pin-offset HPWL 数值以及 density move 前后有限性。覆盖面较小。
- `test_v3_retention_contract.ps1`：验证指定 experiment 根目录只有三个 metrics 文件并存在 overflow percentage 列。

## `plan`

- `NEW_PROJECT_IMPLEMENTATION_PLAN_V2.md`：完整 V2 目标设计，涵盖数学契约、模块 API、pipeline、结果协议、测试矩阵和迁移顺序。
- `NONSMOOTH_GP_RETENTION_GLOBAL_OPTIMIZER_PLAN_V3.md`：V3 metrics-only、global view、optimizer/step、capacity transport 和实验矩阵方案。

计划是目标与审查依据。判断“是否已完成”时，必须对照当前 CMake、源码、测试和实际输出，不能只看计划标题或以前模型的总结。
