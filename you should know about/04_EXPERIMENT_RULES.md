# 04 — 实验、调用链和产物规则

## 1. 实验目录原则

新的正式研究实验优先使用 `nsgp lab`，默认目录：

```text
framework/results/experiments/<run_id>/
├─ params.json
├─ experiment.md
└─ trajectory.csv
```

默认禁止保留：

```text
*.pl
snapshots/
checkpoint/
candidate_dump/
gradient_dump/
density_map_dump/
price_field_dump/
bundle_history_dump/
optimizer_state/
occupancy_cache/
debug/
*.bin
*.npy
*.raw
```

只有用户明确要求某个阶段布局时，才允许 `--save-stage <stage>` 生成 `saved/<stage>.pl`。外部输入 checkpoint 只记录绝对路径和 SHA-256，不复制到实验目录。

`nsgp run` 与 `nsgp lab` 都默认遵守 metrics-only。只有 `modules/historical_exact_replay` 和 `historical_dct_poisson` 的显式历史复现仍会产生较多 `.pl` 和日志；新算法实验不要沿用其重产物方式。

## 2. `run_id` 和实验名字

推荐：`YYYYMMDD_HHMMSS_<case>_<question>_<variant>`，例如：

```text
20260907_143000_adaptec1_optimizer_adam_control
20260907_143500_adaptec1_globalview_active_bundle4
```

名字表达实验问题和变量，不使用无语义的 `test2_final_new`。相同 `run_id` 不得覆盖；程序当前会拒绝已存在目录。

## 3. 必须记录的 provenance

每次实验必须能恢复以下链条：

- 时间、run ID、case、目的和假设；
- Git remote、branch、commit，以及是否有 dirty working tree；
- 数据集 benchmark base；
- 初始化方式：raw，或外部/前序 experiment checkpoint；
- 输入 checkpoint 绝对路径和 SHA-256；
- 父实验 ID；若父实验也有父节点，记录完整父调用链；
- 按顺序列出的 module chain；
- 每个 module 的参数、额外 CLI override、seed、threads；
- evaluator 的 bins、target density 和 exact 语义版本；
- optimizer、step policy、acceptance policy 和 retention policy；
- 起始、best-feasible、last、final-selected 的 exact metrics；
- 每模块 runtime、总 runtime、accepted/rejected、objective evaluations；
- NaN/Inf、异常退出、早停、bundle/optimizer reset 等关键事件；
- 是否显式保留布局，以及原因。

微内核 pipeline 会记录解析后的 module chain、参数、输入 hash（若有）、exact evaluator、threads 和 retention；global-view lab 还记录 optimizer、step 与 acceptance。Git dirty 状态目前仍为 `not_checked`。

## 4. 调用链格式

推荐在人类说明中使用明确箭头并附输入类型：

```text
external:H375/global.pl
  [sha256=...]
  → global_capacity_transport(target=6%, coarse=64, audit=512)
  → global_view_gp(active={0,.25,.5,1}, bundle=4, optimizer=Adam)
  → best_feasible (in-memory; no .pl retained)
```

如果读取前序实验：

```text
experiment:E17/saved/capacity_slack.pl
  ← parent E17
  ← external H375 checkpoint
```

不能只写“使用上次结果”。

## 5. 指标格式

- 内部 evaluator 返回 `overflow_ratio`。
- CSV 和 Markdown 面向人的字段优先命名为 `overflow_percent` 并写 `ratio × 100`。
- 如果因兼容保留名为 `overflow` 的字段，文档必须声明其单位，禁止同一表混用比例与百分比。
- HPWL 保存足够精度，不擅自四舍五入原始证据。
- 至少记录：HPWL、overflow percent、density energy、max density、best feasible、accepted、runtime。
- 所有 stage 边界、候选接受和最终选择必须基于 fresh exact audit。

## 6. 公平对照规则

对照实验必须固定：

- 相同 input checkpoint SHA-256；
- 相同 case、evaluator grid、target density；
- 相同 seed、threads、迭代/评估预算；
- 相同 acceptance 和 retention policy；
- 一次只改变需要研究的因素。

推荐先跑 1 轮/小迭代 smoke，再跑 50 轮对照。若 control 和 variant 都零接受，结论是“在当前方向/步长/可行带下未发现改进”，不是“新方法相同”或“算法已最优”。

## 7. 临时工作区与清理

- 首选同进程传 `Database/Layout`，不产生中间 `.pl`。
- 必须用路径交接时，放进 `%TEMP%\nonsmooth-gp\<run_id>`，并验证删除目标确实位于该目录下。
- 成功和失败都清理临时目录；清理失败应记录，但不能删除外部 checkpoint、数据集或仓库根目录。
- 不对用户命名的目录、`%USERPROFILE%`、仓库根或未解析变量执行递归删除。

## 8. 实验生命周期测试

至少覆盖：

1. metrics-only success：目录只有三个规定文件。
2. failure cleanup：异常退出不残留临时 `.pl`。
3. explicit stage save：只有明确指定时出现 `saved/<stage>.pl`。
4. external input safety：实验前后 checkpoint SHA-256 相同。
5. handoff correctness：内存交接与临时 `.pl` round trip 的 exact metrics 一致。

当前覆盖：PowerShell retention 测试验证 metrics-only、显式 stage save、外部 checkpoint 当前 hash 和临时 workspace 清理；C++ 数值测试验证 placement round trip 的坐标与 exact HPWL。异常 smoke 还验证错误退出后 workspace 不残留。后续仍应把这些 smoke 全部注册进自动 CTest。

## 9. Git 中保存什么

应提交：算法源码、参数模板、测试、方案、交接文档和有长期价值的小型结果摘要。

默认不提交：二进制、完整 trajectory 批量集合、普通实验目录、原始/中间 `.pl`、数据集、stdout 大日志、snapshot 和缓存。若某个关键 checkpoint 确有长期复现价值，应先明确用户意图，再使用专门存储或 Git LFS，并在本目录记录其生命周期。
