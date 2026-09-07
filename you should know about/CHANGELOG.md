# Handoff changelog

这里记录影响后续接手者判断的架构、算法、实验协议和关键数值结果。它不是 Git commit 日志的复制品。

## 2026-09-07 — 收敛为轻量微内核与真实 stage modules

- Git：基于 `00c0780` 开发，最终提交见 Git 历史，分支 `nonsmooth-gp-v1`。
- Changed：新增最小 `microkernel`（registry、stage context、fresh audit、threads、SHA-256、三文件记录）；`main.cpp` 改为 registry 调度；`run` 与 `lab` 默认统一 metrics-only；普通策略源码归入对应 module 并由 CMake 直接编译；共享 exact/I/O/optimizer 归入 `framework/kernel`；移除 `framework/legacy`、重复源码镜像和未接线第二状态骨架。
- Modules：新增显式 `global_capacity_transport`、`density_coordinate` stage；`global_view_gp` 的 lab 源移动到自己的 module。
- Provenance：pipeline 记录解析后的 module chain、每阶段参数、evaluator 和输入 hash；lab 补全 optimizer/step/acceptance、best/last、runtime、接受/拒绝/reset 统计；threads 实际应用到 OpenMP。
- Verification：模块权威构建通过；pipeline 和 lab smoke 通过；默认三文件、显式 save、外部 hash、workspace cleanup、placement round trip 均通过。
- Numerical evidence：2-thread `adaptec1` lab 1 轮得到 HPWL `85999279.708727`、overflow `6.9999980678266%`，1 accepted；该结果仅为 smoke。
- Contract impact：exact HPWL/density 公式不变；overflow 对外仍为百分比。
- Remaining gaps：`global_view_gp` 尚需变为普通 registry stage；Git dirty state、跨平台 SHA-256 和失败摘要仍需完善。

## 2026-09-07 — 建立多模型交接知识库

- Git：文档创建时基线为 `4830135`，分支 `nonsmooth-gp-v1`；本条最终 commit 见 Git 历史。
- Changed：建立接手入口、项目契约、真实架构地图、模块目录、实验/retention 规则、当前状态和维护清单。
- Evidence：交叉检查 CMake、91 个受管源码/配置路径、当前 V3 本地 trajectory、外部 H375 checkpoint 及源码镜像 hash。
- Contract impact：无算法或 evaluator 行为修改。
- Historical finding：当时 CMake 权威 core 位于 `framework/legacy`、module 多为镜像且只有 `lab` metrics-only；这些问题已由上面的微内核重构解决。
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
