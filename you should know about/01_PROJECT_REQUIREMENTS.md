# 01 — 项目要求与不可变契约

## 1. 项目目标

项目用于研究非光滑全局布局算法和数值行为。目标是快速提出、组合、验证算法模块，并通过 exact evaluator 判断结果；不是追求企业软件式抽象、服务化、GUI、完整物理设计流程或大规模历史兼容层。

工作优先级：

1. 数值正确性和评价语义一致。
2. 实验可复现、调用链清晰。
3. 模块可替换、可重排、可独立测试。
4. 代码和产物轻量，参数面保持可理解。
5. 性能优化必须有测量依据，不能以牺牲 exact audit 为代价。

## 2. 非光滑算法要求

正式 non-smooth pipeline 的权威目标和可行性判断必须来自 exact evaluator。允许 epsilon-active、active-face ensemble、bundle/direction memory、subgradient、exact-audited relocation、exact backtracking、trust radius、容量粗网格候选生成；粗糙或平滑量只能帮助生成搜索方向，不能冒充最终指标或最终接受条件。

正式路径禁止：

- 用 weighted-average 或 log-sum-exp 线长作为最终 HPWL；
- 用 electrostatic/DCT/Poisson/smoothed density 作为最终 overflow；
- 用 Moreau envelope 等光滑目标替代挑战目标；
- 将 `historical_dct_poisson` 混入 `challenge_nonsmooth` pipeline；
- 因性能原因跳过模块边界或最终的 fresh exact audit。

Pure-descent research stages may use the exact audit only for telemetry and
dynamic-lambda feedback.  They must state this explicitly and must not turn an
overflow/HPWL measurement into a candidate acceptance, rollback, or final
best-feasible restore rule.

历史 H219/H221 DCT/Poisson 仅用于复现已有调用链，必须在名称、日志和报告中明确标为 historical smooth surrogate。

## 3. 权威评价定义

对 net `e`：

```text
pin_x = node_center_x + pin_offset_x
pin_y = node_center_y + pin_offset_y

HPWL_e = max(pin_x) - min(pin_x)
       + max(pin_y) - min(pin_y)

HPWL = sum_e net_weight_e * HPWL_e
```

对 density bin `b`：

```text
occupancy_b = movable rectangle 与 bin 的精确相交面积
            + physical fixed macro 与 bin 的精确相交面积

rho_b = occupancy_b / bin_area

density_energy = sum_b 0.5 * bin_area
                 * max(rho_b - target_density, 0)^2

overflow_ratio = sum_b bin_area
                 * max(rho_b - target_density, 0)
                 / total_movable_area

overflow_percent = 100 * overflow_ratio
```

额外语义：

- 内存中 `Node.x/y` 是中心坐标；Bookshelf `.pl` 是左下角坐标，读写时必须转换。
- pin offset 相对中心坐标。
- fixed macro 占容量且位置不可改变。
- `terminal_NI` 不占物理容量。
- epsilon-active 只改变方向选择，不能改变上面定义的 exact HPWL 数值。
- 最终结果必须重新构造或刷新 evaluator 后审计，防止缓存与布局状态不一致。

## 4. 模块化要求

- 保持 `framework + modules` 的两层微内核结构：framework 只提供 I/O、权威 evaluator、通用 optimizer、模块注册/调用和实验记录；module 实现一个优化阶段及其参数。
- pipeline 只是模块序列，不是写死的新算法。模块应能删除、重复、重排，并通过显式 layout/checkpoint 交接。
- 模块之间默认传布局状态，不共享危险的 occupancy cache、optimizer moment 或隐含全局变量。
- 跨模块 solver state 必须显式声明、版本化并记录 hash；否则在模块边界重置。
- 新模块优先使用统一的 `Database`、`ExactHpwl`、`ExactOverlapDensity` 和 optimizer API，不得私自重新定义另一套 HPWL/overflow。
- 不建立大型 plugin system、DI container、artifact manager 或企业级配置框架。注册表、JSON 参数和小型结构即可。
- 当前运行状态统一使用 `ea::Database`，模块协议统一使用 `framework/code/microkernel.hpp`；不要再引入第二套 `Problem/Layout` 状态模型。

## 5. 数值实验要求

- 每次算法改动都应由可证伪的数值问题驱动，并有 control/baseline。
- 默认从 `adaptec1` 开始；只有单 case 结果可信后才扩大到八 case。
- 比较实验必须使用相同数据、相同输入 checkpoint SHA-256、相同 evaluator grid/target density、相同 seed、threads 和迭代预算。
- 报告 exact HPWL、overflow 百分比、density energy、max density、runtime、accepted/rejected 数以及 best-feasible 与 last 的区别。
- 不得只报告“看起来更好”；失败、零接受步、未达目标和数值异常也必须记录。
- 调参应放在 JSON 或明确 CLI 参数中，不通过复制源文件制造实验分支。

## 6. 灵活与轻量化要求

- 优先最小实现和小型可组合函数。
- 默认不保存大对象；布局、方向、occupancy 和 optimizer state 应在无后续用途时释放。
- 不迁移或复制旧项目数百个实验目录。
- 不因抽象而重复构建昂贵 evaluator 或复制整份数据库。
- 新参数必须有明确实验问题、默认值、边界校验和记录位置。
- 一项实验结束后，应能仅凭代码提交、参数、输入 hash 和指标记录理解并复现。

## 7. Git 和 GitHub 要求

- 工作分支默认是 `nonsmooth-gp-v1`，远端为 `origin`。
- 开始前检查分支、HEAD 和 dirty 状态；用户已有修改不能覆盖或丢弃。
- 一个提交解决一个可说明的研究/架构问题，并带上测试或实验证据。
- 提交源码、配置、测试、方案和交接文档；不提交数据集、编译产物、常规 `.pl`、snapshot、缓存或批量实验目录。
- `.gitignore` 已忽略 `build/`、`*.exe`、常规 `framework/results/**`、`framework/experiment_logs/**` 和 module results。
- 用户要求保存进度或任务完成时，应提交并推送到 GitHub，然后报告分支、commit hash、测试和实验结果。
- 禁止使用会丢失用户工作树的 `git reset --hard`、`git checkout --` 等操作，除非用户明确授权具体目标。
