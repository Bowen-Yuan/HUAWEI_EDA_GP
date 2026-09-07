# 02 — 顶层架构与代码地图

## 1. 当前真实构建图

```text
Bookshelf dataset / external .pl
              │
              ▼
framework/legacy/src/bookshelf.cpp
              │  ea::Database（当前运行时权威数据结构）
              ├──────── ExactHpwl
              ├──────── ExactOverlapDensity
              └──────── placement/search modules
                              │
        ┌─────────────────────┴──────────────────────┐
        ▼                                            ▼
framework/code/main.cpp                  framework/code/experiment_lab.cpp
run / batch / audit                      lab：V3 metrics-only 实验
        │                                            │
旧式 stage .pl + manifest                params.json + experiment.md + trajectory.csv
```

`CMakeLists.txt` 将 `framework/legacy/src` 中 14 个源文件编译为 `nsgp_exact_core`，再链接给各个可执行文件。因此，当前权威算法实现位于 `framework/legacy`；目录名虽叫 legacy，但它不是可忽略的废弃代码。

## 2. 顶层目录

| 目录/文件 | 内容与职责 | 是否进入正常构建 |
|---|---|---|
| `CMakeLists.txt` | C++17、OpenMP、四个 target、Windows Advapi32 SHA-256 依赖 | 是，权威构建描述 |
| `AGENTS.md` | 强制后续模型先读并同步维护本交接目录 | 否 |
| `README.md` | 简短使用说明 | 否 |
| `framework/code` | CLI、V3 lab、未来轻量模块 API 骨架 | `main.cpp`、`experiment_lab.cpp` 进入 `nsgp` |
| `framework/legacy/include/epsilon_active` | 当前 `ea` API、配置和结果结构 | 是 |
| `framework/legacy/src` | 当前 exact evaluator、I/O、optimizer 和算法实现 | 是，组成 `nsgp_exact_core` |
| `framework/params` | 默认 case、默认 evaluator、pipeline 与 V3 实验示例配置 | 运行时/说明 |
| `framework/results` | 普通运行和 V3 实验结果；默认被 Git 忽略 | 否 |
| `framework/experiment_logs` | 历史 H375 完整复现日志；默认被 Git 忽略 | 否 |
| `modules` | 面向研究者的模块目录、参数和部分源码镜像 | 参数被 pipeline 使用；大多数 code 镜像不由 CMake 单独编译 |
| `plan` | V2 项目设计和 V3 retention/global-view 方案 | 否；目标文档，不等于现状 |
| `scripts` | 历史 H375 调用链复现脚本 | 显式运行才使用 |
| `tests` | exact 数值契约测试和 retention PowerShell 测试 | 部分由 CTest 使用 |
| `third_party/json.hpp` | vendored nlohmann JSON 单头文件 | 是 |
| `bin`、`build` | 本机构建产物 | 不提交 |
| `you should know about` | 本交接知识库 | 必须随架构同步提交 |

## 3. `framework/code`

### `main.cpp`

当前紧凑 CLI 和旧 pipeline runner：

- 命令：`run`、`batch`、`audit`、`list-modules`；遇到 `lab` 时转发给 `run_global_view_lab`。
- JSON pipeline 的模块名静态分派到 `layout_init`、`hpwl_adam`、`exact_joint_gp`、`exact_recovery`、`surplus_bisection`、`equal_shape_swap`。
- `challenge_nonsmooth` 会拒绝 `historical_dct_poisson`。
- 每个模块前后 fresh audit。
- 旧 `run` 当前仍在 `framework/results` 写 stage `selected.pl`、`stage_manifest.json`、`case_summary.csv` 和最终 `.pl`；不要误认为全项目都已 metrics-only。
- `hpwl_adam`/`exact_joint_gp` 通过 `ea::global_place` 执行，并使用 `_scratch` 输出目录，仍有旧式写盘副作用。

### `experiment_lab.cpp/.hpp`

V3 研究入口，直接复用 `ea` exact core：

- 默认从外部 H375 7% `adaptec1` checkpoint 读取布局并用 Windows CryptoAPI 计算 SHA-256。
- 实现 active ensemble、minimum-norm simplex 聚合、4 步 temporal bundle、exact backtracking、trust radius 和轻量 capacity transport。
- 同一进程内传递布局，`ExperimentWorkspace` 位于系统临时目录并在析构时清理。
- 正常实验只输出三个 metrics 文件；仅显式 `--save-stage NAME` 时写 `saved/NAME.pl`。
- 当前限制见 `05_CURRENT_STATE.md`：元数据尚未覆盖 V3 方案所有字段，threads 未实际设置 OpenMP，SHA-256 路径是 Windows-only。

### 架构骨架头文件

- `common.hpp`：`Real`、`NodeId`、`Rect`。
- `problem.hpp`：不可变问题描述 `Problem`、`NodeInfo`。
- `layout.hpp`：坐标、orientation 和 revision。
- `metrics.hpp`：轻量 `Metrics`。
- `module_api.hpp`：`RunContext`、`StateArtifact`、`ModuleStats`、`ModuleResult`。

这些结构目前不是 `ea::Database` 主运行路径的权威实现。若继续完成 V2 解耦，应逐步接线并一次只保留一个权威 evaluator；不要简单复制 `ea` 状态形成长期双轨。

## 4. `framework/legacy/include/epsilon_active`

| 文件 | 核心接口 |
|---|---|
| `types.hpp` | `Node`、`Pin`、`Net`、`Row`、`Database`、`DensityMetrics`、`IterationMetrics` |
| `bookshelf.hpp` | 读取 benchmark、加载/写出 `.pl`、中心高斯初始化 |
| `hpwl.hpp` | `ExactHpwl::evaluate` |
| `density.hpp` | `ExactOverlapDensity`，full evaluate、price direction、node/group incremental move audit |
| `optimizer.hpp` | Adam、AMSGrad、AdaGrad、HeavyBall、SGD、NormalizedSGD、DualAveraging |
| `lambda_controller.hpp` | Dreamplace、Trajectory、Ratio 三种 lambda policy |
| `placer.hpp` | `PlaceConfig`、`PlaceResult`、`global_place` |
| `batch_acceptance.hpp` | exact batch proposal acceptance/backtracking |
| `recovery.hpp` | overflow cap 下的 exact HPWL recovery |
| `compact_recovery.hpp` | support contraction |
| `density_coordinate.hpp` | exact overlap coordinate descent |
| `bisection.hpp` | 递归超图/容量二分 |
| `coarse_flow.hpp` | 粗粒度容量流 |
| `transport.hpp` | excess-to-capacity transport、auction/Hilbert/identity exchange |
| `swap_recovery.hpp` | equal-shape、net-aware、assignment/permutation recovery |

对应实现全部在 `framework/legacy/src`。`legacy_stage_main.cpp` 不进入 core library，而是单独构建成历史复现 CLI。

## 5. 配置层

- `framework/params/defaults.json`：默认数据目录、seed 219、1 thread、canonical 512×512、target density 1.0。注意旧 pipeline 文件可自行覆盖为 32/64 和 0.9。
- `framework/params/cases.json`：默认 `adaptec1` 与八 case 清单。
- `framework/params/pipelines/smoke.json`：32×32、三阶段快速 smoke。
- `framework/params/pipelines/reference_nonsmooth_chain.json`：64×64 的模块组合示例，不是最终算法或 canonical 实验配置。
- `framework/params/experiments/*.json`：V3 E0/E5 命令说明，目前不是被 `lab` 自动读取的完整配置文件。
- `modules/*/params/*.json`：各模块参数；pipeline 的相对路径从 pipeline 文件父目录解析。

## 6. 源码镜像风险

`modules/*/code` 多数文件与 `framework/legacy/src` 是相同副本，但 CMake 当前编译后者。盘点时以下镜像相同：HPWL、placer、lambda controller、batch acceptance、recovery、compact recovery、bisection、swap recovery。

`modules/hpwl_adam/code/optimizer.cpp` 已经与权威 `framework/legacy/src/optimizer.cpp` 不同：后者包含 V3 新增的 NormalizedSGD 和 DualAveraging。修改算法时应先改权威实现；若模块镜像仍被保留，应同步或明确删除镜像，不能让两者无说明漂移。

## 7. 数据与坐标边界

- Bookshelf benchmark base 传入形式是 `<dataset>/<case>/<case>`，由 `.aux` 继续解析 `.nodes/.nets/.pl/.scl/.wts`。
- 数据不能提交进本仓库。
- 移动单元必须 clamp 到 region，固定单元不得被移动。
- `.pl` round trip 必须保持中心/左下角转换和 fixed/terminal 标志。
