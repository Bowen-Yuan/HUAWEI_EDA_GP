# Handoff changelog

这里记录影响后续接手者判断的架构、算法、实验协议和关键数值结果。它不是 Git commit 日志的复制品。

## 2026-09-14 — V8 restart 实验与 DREAMPlace-inspired 模块（spectral_bb_gp、overflow_epsilon_continuation_gp）

- Git：基于 V8 Commit A（e4aecb8，nonsmooth-gp-v1）；最终提交见 Git 历史。
- Changed：新增两个 DREAMPlace-inspired 阶段：
  `spectral_bb_gp`（Barzilai-Borwein 步长估计，BB1/BB2 公式 + curvature 门控 +
  clip + growth limiter + restart fallback；方向 oracle 与 exact joint GP 相同，
  仅标量学习率自适应；要求 SGD + 无 batch acceptance）与
  `overflow_epsilon_continuation_gp`（overflow 驱动的 epsilon-active 半径 continuation，
  band 切换可选 optimizer reset；exact 指标不变）。
  kernel 新增纯逻辑原语 `ea::SpectralStep`（spectral_step.hpp）与
  `ea::EpsilonContinuation`（epsilon_continuation.hpp）；
  `PlaceConfig` 新增 spectral_step/spectral/epsilon_continuation/iteration_hook；
  placer GP 循环集成两者并输出每轮 lr。契约测试
  `nsgp_dreamplace_inspired_contracts`（BB 公式/fallback/clip/growth/EC 单调/端点/band）。
- New V8 experiments（h221 checkpoint SHA 9d0bd999…，canonical 512×512/1.0，threads=1，metrics-only）：
  V8_REPAIRED_V7_CONTROL（cap15→bridge120→retighten180→polish×2）
  = **90,428,943.156329 @ 6.9999992385929%**（新 baseline）；
  R1 等预算 150：continuous **98,139,619.875507 @ 5.0737%** vs 3×50 restart
  102,847,407.27813 @ 3.1587%（HPWL 口径 restart 负向 −4.8%，两者 feasible）；
  R2 V7×2 = 87,477,937.531673 @ **7.4436%（infeasible）**，第二轮宏循环未能收回 feasibility；
  BB screen 100：BB0 fixed-SGD 108,334,816.66594 @ 5.0971% vs
  BB2-method 109,357,983.39731 @ 5.4527% vs **BB1-method 106,609,642.22966 @ 5.3065%（晋级）**；
  BB1-method 300 = **105,122,976.93855 @ 5.1116%**；
  EC：fixed-epsilon EC0 99,124,454.706286 @ 5.2284% 正常收敛，但
  widening continuation（EC1 4→16 bins、EC2 4→32 bins）overflow 发散至 43–56%，
  stage 结束回滚输入，**continuation widening 明确负向**；
  band-change optimizer reset 在 overflow 压力震荡期反复清零 moments（EC0 上实证有害）。
- Contract impact：exact HPWL/density 公式不变；BB 只改变标量步长；
  EC 只改变方向 oracle 半径；两者都是 pure-descent 无 rollback。
- Documentation updated：03、05、07 与本文件。
- Remaining gaps：Fusion 未运行（EC 负向未满足独立正向前提）；
  BB1-method 与 fixed-SGD 的对照在 100 步口径成立，300 步口径缺同长 fixed-SGD 对照；
  EC 的 hysteresis band 设计未实现；单 case/单 checkpoint 局限不变。

## 2026-09-14 — V8 correctness/audit repair（finite-radius、deterministic audit、CSV schema、preconditioner、high-degree）

- Git：基于 6fc2068（nonsmooth-gp-v1）执行 V8 方案 Commit A；最终提交见 Git 历史。
- Changed：`ExactOverlapDensity` 新增 `DensityAccumulationMode` 与
  `rebuild_occupancy(mode)`；`exact_audit` 走 DeterministicCanonical（单线程固定顺序），
  搜索路径保留 FastParallel。`ExactHpwl::evaluate` 新增 `HighDegreeMode`
  （Ignore=历史；ExactExtrema=epsilon=0 exact extremal-face subgradient，ties 平分权重）。
  finite-radius probe 重写为共享 `finite_radius_probe`：只接受
  `energy_delta < -descent_tolerance` 的方向、secant 分母用 clamp 后实际位移、
  probe 用 `evaluate_move` 增量计算且不 commit，返回 local query 统计。
  nonlocal 引擎把 exact HPWL 传入 LambdaController 与 `AuxiliaryContext.hpwl`
  （修复恒 0）；`AuxiliaryGradient` 回调返回 `AuxiliaryQueryStats`。
  `StageStats` 增加 local_density_queries/local_hpwl_queries；
  trajectory.csv 统一 19 列 schema，修复 iteration 行错位与 stage 行缺列，
  global_iteration 跨 stage 单调，telemetry module 记录具体模块名。
  `nsgp audit` 默认 clamp+canonical（audit_mode=canonical_clamped），
  `--pipeline` 可读非 canonical grid，`--raw-unclamped-audit` 显式 NONCANONICAL 诊断。
  `PlaceConfig.preconditioner`（默认 weighted_pin_sum）与
  `PlaceConfig.high_degree_mode`；Database 新增 `node_pin_weight_sum`。
- Tests：numeric contracts 增加 density delta/group move/deterministic audit/
  weighted-net HPWL delta/exact-extrema/pin weight sum；nonlocal oracle contracts
  增加 FR1–FR5、ctx.hpwl 真值契约；新增 `nsgp_experiment_log_contracts`
  （19 列 schema、列位置、global_iteration 单调、local query 记账）。
  四个测试目标全部通过；finite-radius 1-round 与 GP 1-round smoke 通过。
- Numerical evidence：h221 best.pl canonical clamped audit =
  HPWL 109,837,622.00901 / overflow 7.8193002333393%（与历史参考一致，
  raw audit artifact 差异因此消除）。
- Contract impact：exact HPWL/overflow 公式不变；deterministic audit 只固定累加顺序；
  high-degree exact_extrema 是合法 non-smooth subgradient，不引入 smooth surrogate。
- Documentation updated：05 与本文件。
- Remaining gaps：weighted preconditioner 与 exact_extrema 的历史 replay 消融、
  restart/BB/EC 数值实验见 V8 方案后续提交。

## 2026-09-08 — V7 H252-H257 四模块联合链：从 h221 命中 88M/7%

- Git：基于 e230425（nonsmooth-gp-v1），最终提交见 Git 历史。
- Changed：`exact_recovery` module adapter 新增可选 JSON 键 degree_limit/axis_separated/
  breakpoint_oracle/net_block_density_direction/net_block_degree_limit/net_block_max_nodes/
  net_block_max_blocks；`exact_joint_gp` adapter 新增可选 net_batch.weight/degree_limit。
  默认值全部等于原 RecoveryConfig/NetBatchConfig 默认，既有 pipeline 行为不变。
- Numerical evidence：五阶段链（cap15 recovery → bridge dws .5/cap .15 → retighten dws 4.75/
  cap .07 → cap07 polish ×2）从 h221 best.pl 出发，11 个 300 步 run 筛 optimizer/步长。
  adam lr .002/.001 final 87,599,730.55 @ 6.9999995%（命中 88M/7%，−20.3% vs 起点）；
  amsgrad 89.83M；sgd/heavy-ball/nsgd retighten 失败（~10%）；步长 .002 为 adam 峰值。
  与 V6 单阶段（94.06M）相比，分段压力调度收益 +7.3%。
- Contract impact：exact HPWL/overflow 语义不变；新增键只暴露既有 config 字段；实验
  metrics-only，输入 checkpoint hash 前后不变。
- Documentation updated：05 §3.4 与本文件。
- Remaining gaps：单 case/单 checkpoint/threads=16；未做 8 case 外推与 seed 稳定性。

## 2026-09-08 — Five independent non-local non-smooth pure-descent oracle stages (initial implementation)

- Git: pending on `nonsmooth-gp-v1`.
- Changed: added shared direction-normalization/trajectory-lambda/unconditional-update mechanics and registered
  multiscale-active, cut-pressure, transport, charge, and finite-radius exact-secant stages.  Each strategy owns
  its direction construction and its own smoke JSON; no acceptance gate, backtracking, rollback, best-feasible
  restore, or cross-module state was added.
- Numerical evidence: five canonical H375 one-step metrics-only health runs completed. The 100-round AMSGrad
  primary ended at 14.81%/69.06%/24.21%/87.16%/13.89% overflow for multiscale/cut/transport/charge/finite(4 probes),
  respectively: an explicit negative result, not a ranking.
- Contract impact: exact evaluator and final metrics are unchanged; pure-descent stages now explicitly use audit only
  for measurement/lambda feedback and hand off their last applied layout.
- Documentation updated: 01, 02, 03, 05, 07 and this file.
- Remaining gaps: per-iteration exact telemetry is now in trajectory.csv; synthetic strategy contracts, exact-only
  pure-descent control, calibrated optimizer/LR screens, and pipeline fusion remain.

## 2026-09-07 — h221 → 86M 连续下降调参方案 V6

- Git：`pending` on `nonsmooth-gp-v1`，方案基线 `0fbf1b9`。
- Changed：新增 `plan/NONSMOOTH_GP_H221_86M_TUNING_PLAN_V6.md`，规定从历史 h221 checkpoint
  进行 300 轮、无 hard feasibility gate、无 accept/reject/backtracking 的 `exact_joint_gp` 调参；以
  `ratio` λ 策略和由 86M/7% 推导的 `overflow_baseline=0.0894027156` 为主线，分阶段筛选 λ、
  AMSGrad/Adam 步幅与动量，并定义逐轮标量、三类收敛图和活结果台账。
- Numerical evidence：本条仅设计实验；引用的 100 轮 7-optimizer 基线来自当前 V5 记录，没有新增 run。
- Contract impact：exact evaluator 与 metrics-only retention 不变；明确禁止保存中间/最终 placement，
  只保留逐轮标量与汇总。
- Documentation updated：07 与本文件。
- Remaining gaps：A–F 调参批次尚未执行；86M/<7% 尚未达到。

## 2026-09-08 — V6 h221→86M 调参执行（A–D 完成，未命中 86M，相对起点 −14.36%）

- Git：基线 0fbf1b9（dirty 交接文档 + V6 方案），最终提交见 Git 历史，分支 nonsmooth-gp-v1。
- Changed：无源码改动（纯参数调优轮次）；新增本地执行器 build/h221_86m_tuning/run_tuning.ps1
  （不入库），生成 config/pipeline、解析 [GP] 逐轮标量、汇总 runs.csv 并逐 run 跑 retention。
- Protocol：exact_joint_gp 从 h221 best.pl（SHA-256 9d0bd999…，canonical 109,837,622.009 /
  7.8193%）出发，无接受/回溯/硬上限，512×512/1.0，threads=1，300 轮，metrics-only。
- Numerical evidence：A（A0 复现 V5 前 100 轮到打印精度；A1 h=240 得 98,387,144.68 @ 4.7642%）；
  B（B5 ratio ob=.085 itv=2 → 95,827,835.14 @ 6.9621% 晋级；interval=5 时 ratio 早期不分化，
  B0/B2 逐位一致）；C（lr .002 最佳；md 缩小单调有效，md=0.25 → **94,061,868.14 @ 6.7474%**
  全程最佳，md=0.125 反转；adagrad/dual-averaging 远差）；D（动量微调全部劣于 0.90/0.99）。
  30 个实验目录全部通过 retention，checkpoint hash 全程未变，无任何 placement 保存。
- Contract impact：exact evaluator 与 metrics-only retention 不变；E 阶段触发条件未满足
  （最佳 run 最后 50 轮仍下降 0.647% > 0.2%），未实现 schedule 机制。
- Documentation updated：V6 方案 §10 台账/§11 清单、05 §3.3、本文件。
- Remaining gaps：未命中 86M/<7%（差 8.57%）；证据支持的下一步是最佳配置延长到 500 轮；
  单 case、单 checkpoint、单确定性轨迹的局限不变。

## 2026-09-08 — Adaptive lambda GP stage 与可组装 step policy（V4 方案执行）

- Git：基于 `cca2e8b` 开发，最终提交见 Git 历史，分支 `nonsmooth-gp-v1`。
- Changed：新增共享 `ea::StepController`（constant/cosine/trust，`framework/kernel/src/step_policy.*`，
  复用 `ea::Optimizer`，未建第二套 optimizer hierarchy）；新增 `adaptive_lambda_gp` stage
  （overflow funnel + log-domain PI(D) lambda + exact funnel acceptance + 内存内 best-feasible restore；
  wire 与 density direction 分别 RMS normalize）。未改动任何既有模块的算法行为。
- Tests：新增 `tests/test_adaptive_lambda_contracts.cpp`（CTest 目标 `nsgp_adaptive_lambda_contracts`）：
  step policy 单调/bounds、lambda controller 合成 Case A–E、funnel acceptance 四象限、合成 3-cell 布局上
  的 explore-accept + restore、infeasible-input 抛错、7 optimizer × 3 step policy finite/bounded/reset
  可复现 smoke；`nsgp_tests` 全部通过。
- Build：新增 `scripts/build_manual_gcc.ps1`（无 CMake 机器的手工 MinGW 构建，目标与 CMakeLists
  一一对应；Git 元数据经 forced-include 头注入）。
- Numerical evidence（canonical checkpoint SHA ba0fd17e…，512×512/1.0，threads=1，seed 219，
  metrics-only，21 个实验目录全部通过 retention 契约，输入 hash 前后不变）：
  A0 strict-cap 对照 50 轮只沿容差边界移动 −3.51 HPWL；
  A1/A2/A3（50 轮）证明收益来自 corridor 而非新代码路径，但 50 轮内 recovery 未成功；
  B/C optimizer 筛选中 dual-averaging 探索最强且唯一尝试回收；
  **D 系列 100 轮 dual-averaging/constant 得到 best-feasible HPWL 85,445,335.249799 @ 6.8253943126%，
  相对输入 −553,983.06（−0.644%，material improvement）**；200 轮 E 系列改善反而较小（−0.349%）。
  Adam 系全部能探索但 recovery 失败；cosine step 未回收。
- Contract impact：exact HPWL/density 公式不变；overflow 仍以百分比输出；metrics-only retention 不变；
  stage 边界输出无条件为 canonical feasible（≤ final + 1e-9），输入在 final 目标下不可行时显式抛错。
- Documentation updated：02、03、05、07 与本文件。
- Remaining gaps：recovery 失败的结构性根因（density direction 最小化 energy 而非 overflow）待解决；
  SHA-256 仍为 Windows-only；constant 与 trust 等价说明 trust 机制在该 stage 未发挥作用；
  证据仅限 adaptec1 单 case，未做八 case 与 threads=8 复核。

## 2026-09-07 — 收敛为轻量微内核与真实 stage modules

- Git：基于 `00c0780` 开发，最终提交见 Git 历史，分支 `nonsmooth-gp-v1`。
- Changed：新增最小 `microkernel`（registry、stage context、fresh audit、threads、SHA-256、三文件记录）；`main.cpp` 改为 registry 调度；`run` 与 `lab` 默认统一 metrics-only；普通策略源码归入对应 module 并由 CMake 直接编译；共享 exact/I/O/optimizer 归入 `framework/kernel`；移除 `framework/legacy`、重复源码镜像和未接线第二状态骨架。
- Modules：新增显式 `global_capacity_transport`、`density_coordinate` stage；`global_view_gp` 已成为普通 registry stage，且与薄 `lab` 入口共用模块内搜索实现。
- Provenance：pipeline 记录开始时间、解析后的 module chain、每阶段参数、evaluator、输入 hash 和可选研究 metadata；lab 补全 optimizer/step/acceptance、best/last、runtime、接受/拒绝/reset 统计；CMake 写入 Git remote/branch/commit/dirty，threads 实际应用到 OpenMP。
- Verification：模块权威构建通过；pipeline 和 lab smoke 通过；默认三文件、失败三文件摘要、显式 save、两种 input schema 的外部 hash、workspace cleanup、placement round trip 均通过。
- Numerical evidence：早期 2-thread `adaptec1` lab 1 轮得到 HPWL `85999279.708727`、overflow `6.9999980678266%`，1 accepted；接入 registry 后用 active ensemble + bundle 同参数核对，registry 与 lab 末态均为 HPWL `85999318.311948`、overflow `6.99999030053%`、0 accepted。两者均仅为 smoke。
- Contract impact：exact HPWL/density 公式不变；overflow 对外仍为百分比。
- Remaining gaps：跨平台 SHA-256、构建后 Git 状态刷新、依赖外部数据的 CTest 接线和完整 50 轮公平对照仍需完善。

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
- Changed：`exact_joint_gp` 模块 JSON 参数面扩展（optimizer.name/beta、lambda.policy/update_interval/
  trajectory_horizon/stop_overflow、log_every、adaptive_epsilon），支持 trajectory-lambda 纯下降实验；
  新增 V5 h221 续跑 7-optimizer 对照（见 05_CURRENT_STATE 3.2）。
