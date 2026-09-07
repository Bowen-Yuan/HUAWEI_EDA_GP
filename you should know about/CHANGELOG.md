# Handoff changelog

这里记录影响后续接手者判断的架构、算法、实验协议和关键数值结果。它不是 Git commit 日志的复制品。

## 2026-09-07 — 建立多模型交接知识库

- Git：文档创建时基线为 `4830135`，分支 `nonsmooth-gp-v1`；本条最终 commit 见 Git 历史。
- Changed：建立接手入口、项目契约、真实架构地图、模块目录、实验/retention 规则、当前状态和维护清单。
- Evidence：交叉检查 CMake、91 个受管源码/配置路径、当前 V3 本地 trajectory、外部 H375 checkpoint 及源码镜像 hash。
- Contract impact：无算法或 evaluator 行为修改。
- Key finding：CMake 权威 core 位于 `framework/legacy`；`modules/*/code` 多数是镜像，且 optimizer 镜像已经漂移。只有 `nsgp lab` 满足基础 metrics-only 目录形式，旧 `run` 仍写 placement。
- Remaining gaps：见 `05_CURRENT_STATE.md` 的技术缺口列表。

## 2026-09-07 — V3 retained global-view experiments

- Git：`4830135` on `nonsmooth-gp-v1`。
- Changed：加入 `nsgp lab`、基本 metrics-only retention、外部 checkpoint SHA-256、active ensemble、temporal bundle、exact backtracking、trust radius、轻量 capacity transport、NormalizedSGD 和 DualAveraging。
- Numerical evidence：E0/E5/E5-backtrack 在 50 轮严格 7% cap 下均为零接受；E6 capacity+active+bundle 最佳记录为 HPWL `86390926.507477`、overflow `7.0000002998563%`。
- Contract impact：正式评价仍使用 exact weighted pin-offset HPWL 和 exact rectangle/bin overlap；对外 overflow 使用百分比。
- Remaining gaps：V3 metadata、thread 应用、retention T2–T5、旧 runner 统一仍未完成。

## 2026-09-06 — H375 historical replay and provenance

- Git：`6d2d0c7`。
- Changed：加入历史 H219/H221 DCT/Poisson 独立 target、H252–H375 exact 模块调用脚本和独立实验日志区域。
- Numerical evidence：本地 raw replay 最终 HPWL `88664604.96`、overflow `6.99999968364%`。
- Contract impact：历史 smooth surrogate 与正式 non-smooth pipeline 分离。

## 2026-09-06 — Initial modular project

- Git：`a64dc60`。
- Changed：建立 C++17/OpenMP 项目、exact core、静态模块 runner、JSON 参数、测试和空结果目录。
- Contract impact：确立 exact HPWL/density 与模块化研究主干。
