# 06 — 后续模型维护清单

## 修改前

- [ ] 完整阅读 `you should know about`。
- [ ] 确认当前 branch、HEAD、remote 和 `git status --short`。
- [ ] 检查用户已有未提交修改并避让。
- [ ] 读取任务涉及的当前源码、测试、参数与最近计划；不只读旧方案。
- [ ] 确认要修改的是 CMake 实际编译的权威文件，还是仅用于说明的 module 镜像。
- [ ] 写出本次可证伪的数值问题、baseline、唯一变量和成功/失败判据。
- [ ] 锁定 case、输入 checkpoint SHA-256、grid、target density、seed、threads 和预算。

## 修改算法时

- [ ] 保持 exact HPWL 和 exact density 评价语义不变。
- [ ] 搜索方向与报告指标在类型/命名上明确分离。
- [ ] fixed node 不移动，`terminal_NI` 容量语义不改变。
- [ ] 模块边界重置隐含状态，或显式记录 state artifact。
- [ ] 不建立重复 evaluator、重复 optimizer hierarchy 或新的平行 runner。
- [ ] 新参数进入 JSON/CLI 校验和实验记录，不复制源文件制造变体。
- [ ] 粗网格或 surrogate 只产生候选，最终接受使用 canonical fresh exact audit。

## 修改实验框架时

- [ ] 正常 experiment 默认只保留 `params.json`、`experiment.md`、`trajectory.csv`。
- [ ] 仅显式 `--save-stage` 才保留 placement。
- [ ] 输入 checkpoint 保存 path+SHA-256，不复制且前后 hash 不变。
- [ ] 记录完整父实验和 module 调用链。
- [ ] overflow 对外显示为百分比，字段名/单位无歧义。
- [ ] 记录 best-feasible、last、selected，而不是只留一个含糊的 final。
- [ ] 成功/失败都清理受控临时目录；删除前验证目标位于 `%TEMP%\nonsmooth-gp`。

## 验证

- [ ] 编译 `nsgp`、历史 target（若受影响）和 `nsgp_tests`。
- [ ] 运行 CTest 数值契约测试。
- [ ] 先做 1 轮 smoke，再做正式实验。
- [ ] 运行 retention contract。
- [ ] 对 evaluator 改动增加精确小例：pin offset、fixed occupancy、terminal_NI、rectangle intersection、node/group delta、round trip。
- [ ] 对 algorithm 改动报告 exact before/after、接受统计和失败情况。
- [ ] 检查没有意外 `.pl`、snapshot、cache、binary 或数据集进入 Git。

## 必须同步更新本目录的触发条件

发生以下任一变化，必须在同一个 commit 更新相应交接文件：

- 顶层目录、CMake target、入口 CLI 或默认命令变化；
- 权威 evaluator、公式、坐标或 fixed/terminal 语义变化；
- module 新增、删除、重命名、职责或参数变化；
- runner、retention、trajectory schema 或调用链格式变化；
- 数据集位置、canonical grid/target、默认 checkpoint/hash 变化；
- 新 baseline、关键正结果、关键负结果或被推翻的结论；
- 已知缺口被解决或出现新的风险；
- Git branch/remote/发布策略变化。

更新映射：

| 变化 | 至少更新 |
|---|---|
| 新人入口/命令 | `00_START_HERE.md` |
| 原则或数学语义 | `01_PROJECT_REQUIREMENTS.md` |
| 目录/构建/依赖 | `02_ARCHITECTURE_AND_CODE_MAP.md` |
| 模块/参数/调用关系 | `03_MODULE_CATALOG.md` |
| 实验 schema/retention | `04_EXPERIMENT_RULES.md` |
| commit、结果、缺口 | `05_CURRENT_STATE.md` |
| 所有实质变化 | `CHANGELOG.md` |

## Git 收尾

- [ ] `git diff --check` 无空白错误。
- [ ] `git status --short` 只包含预期文件。
- [ ] 文档中的文件路径和命令确实存在。
- [ ] 提交信息说明研究/架构目的。
- [ ] 推送 `origin nonsmooth-gp-v1`。
- [ ] 最终回复给出 commit、分支、测试、实验结果和仍存在的限制。

## `CHANGELOG.md` 条目模板

```markdown
## YYYY-MM-DD — <short title>

- Git: `<commit or pending>` on `<branch>`
- Changed: <code/config/architecture changes>
- Numerical evidence: <test/run IDs and exact outcomes>
- Contract impact: <none, or exact description>
- Documentation updated: <files>
- Remaining gaps: <honest limitations>
```
