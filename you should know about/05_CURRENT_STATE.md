# 05 — 当前状态、数值结果与已知缺口

记录日期：2026-09-07（Asia/Shanghai）。记录基线：分支 `nonsmooth-gp-v1`，提交 `4830135`。接手时必须重新运行 Git 只读检查，因为这里的 commit 会随开发变化。

## 1. 已实现并进入当前构建

- Bookshelf I/O 和 `ea::Database`。
- exact weighted pin-offset HPWL。
- exact rectangle/bin density、overflow 和 density energy。
- node/group incremental density move audit。
- Adam、AMSGrad、AdaGrad、HeavyBall、SGD、NormalizedSGD、DualAveraging。
- exact joint GP、lambda controller、batch acceptance。
- recovery、compact recovery、density coordinate、bisection、coarse flow、transport、equal-shape swap。
- `run/batch/audit/list-modules` CLI 和 JSON pipeline。
- 历史 H375 exact/DCT-Poisson 分离回放工具。
- V3 `lab`：metrics-only、checkpoint SHA-256、active ensemble、temporal bundle、exact backtracking、trust radius、轻量 capacity transport。
- 一个 C++ exact 数值测试和一个 PowerShell retention 基本测试。

## 2. 权威 `adaptec1` 输入 checkpoint

默认外部路径：

```text
D:\codex_project\HUAWEI_EDA\epsilon-active\experiments\h375_a1_surplus_recovery10_exchange\run\global.pl
```

SHA-256：

```text
ba0fd17e9111b5bcbd13ed746ab9857106d962652552d26076c6e82007551e6b
```

canonical evaluator：512×512 bins，target density 1.0。

输入 exact 指标：

| HPWL | Overflow |
|---:|---:|
| 85,999,318.311947 | 6.99999030053% |

该文件属于 `epsilon-active` 外部实验目录，只读使用，不复制、不清理、不提交。

## 3. 已完成的 V3 本地实验

这些目录位于被 Git 忽略的 `framework/results/experiments`，不保证在其他主机存在。

| Run | Iterations | Best feasible HPWL | Best overflow | 轨迹 accepted 总数 | 观察 |
|---|---:|---:|---:|---:|---|
| `v3_e0_adam_control_20260907` | 50 | 85,999,318.311947 | 6.99999030053% | 0 | 严格 7% cap 下无可接受改善 |
| `v3_e5_active_bundle_20260907` | 50 | 85,999,318.311947 | 6.99999030053% | 0 | active+bundle 未产生接受步 |
| `v3_e5_active_bundle_backtrack_20260907` | 50 | 85,999,318.311947 | 6.99999030053% | 0 | 加 exact backtracking 后仍无接受步 |
| `v3_e6_capacity_active_bundle_20260907` | 50 | 86,390,926.507477 | 7.0000002998563% | 15 | transport 先完成 24 个移动，但未达到 6% 目标；后续在输入 cap 容差内恢复 |

注意：零接受步只证明当前参数化没有找到可接受下降，不证明 checkpoint 是局部或全局最优。E6 的 accepted 列来自 global-view trajectory，不包含先前 24 个 transport move。

历史 raw-to-H375 本地复现曾得到：HPWL `88,664,604.96`，overflow `6.99999968364%`。它与上面的外部 H375 checkpoint 不是同一个最终数值，不能混作同一 baseline。

## 4. 当前最重要的技术缺口

按优先级排列：

1. **双路径 retention**：`nsgp lab` 是 metrics-only，但 `nsgp run`/`global_place` 仍无条件或默认写 stage/final placement、manifest 和 `_scratch` 内容。V3 retention 尚未统一到 runner/module API。
2. **V3 元数据不完整**：当前 `params.json` 未完整记录 module chain、seed、git branch/commit、dirty state、learning rate、lambda、acceptance policy等；`experiment.md` 未在运行结束后补 final/best/last、runtime 和接受统计。
3. **threads 未生效**：`lab --threads` 当前只解析和记录，没有调用 OpenMP thread 配置；公平实验前应修复并验证。
4. **实验 JSON 未接线**：`framework/params/experiments/*.json` 只是命令示例，`lab` 不会读取其配置。
5. **Retention 测试不足**：目前只验证正常目录三个文件，尚缺 failure cleanup、explicit save、external hash safety、handoff exact equivalence。
6. **临时工作区用途有限**：`ExperimentWorkspace` 会创建/清理目录，但当前 lab 没有真正需要临时 PL 的模块；失败路径和安全边界尚无自动测试。
7. **Windows-only SHA-256**：`experiment_lab.cpp` 非 Windows 分支直接报错；若迁移 Linux，需要小型跨平台 SHA-256 或受控依赖。
8. **源码镜像漂移**：`modules/hpwl_adam/code/optimizer.cpp` 尚未同步 NormalizedSGD/DualAveraging。长期应明确生成/同步机制或移除冗余镜像。
9. **模块 API 骨架未接线**：`nsgp::Problem/Layout/ModuleResult` 与 `ea::Database` 并存，但真正 runner 仍用后者。继续 V2 时要避免形成永久双状态源。
10. **V3 全局视野实现仍是研究原型**：lambda 固定；bundle reset 只覆盖 reject streak，未覆盖 V3 计划列出的所有事件；capacity transport 目的 bin 搜索很粗，且 24 pass 没达到 6%。
11. **轨迹 schema 不完全统一**：lab trajectory 缺 stage/module/global iteration 等统一字段，旧 runner 使用另一份 CSV schema。
12. **自动化构建覆盖有限**：当前 CTest 只注册 `test_legacy_numeric_contracts.cpp`，缺模块级、round-trip、fixed macro、terminal_NI、group delta 和 pipeline smoke 等完整回归。

## 5. 建议下一个模型的起点

若用户没有另行指定，最有价值的下一步是先完成 retention/provenance 一致化：

1. 在统一 `RunContext` 中落地小型 `RetentionPolicy`。
2. 让 `run` 与 `lab` 共享一个实验记录器，而不是再建第三套 runner。
3. 补全 params/experiment 结束摘要和 Git/输入 hash。
4. 让 `--threads` 真正设置 OpenMP，并记录实际线程数。
5. 补齐 T1–T5 测试。
6. 只在以上基础可靠后继续 optimizer/global-view 数值实验。

## 6. 不能被当前结果支持的结论

- 不能宣称 active ensemble 或 bundle 优于/等于 baseline；当前只观察到零接受。
- 不能宣称 capacity transport 已达到 6% overflow；它没有达到。
- 不能宣称全项目已完成 metrics-only；只有 `lab` 满足基本目录形式。
- 不能宣称 H375 raw replay 与外部 checkpoint 数值一致。
- 不能宣称 overflow 7% 表示 legalized。
