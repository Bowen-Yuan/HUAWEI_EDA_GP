# ISPD2005 八数据集 Bundle 方法评测

本轮使用完全相同的 3000 步 bundle 配置运行 adaptec1-4 与 bigblue1-4。
HPWL 始终是严格非光滑的 `max-min`，只对静电密度项进行平滑建模。

## 统一设置

- Bundle：2 cuts，interval 5，prox 2.0，current mix 0.50。
- Lambda：控制基线 6%，恢复目标 6.48%，fine cooldown floor 0.25。
- 可行性筛选：10%。该值来自 DREAMPlace 官方代码的默认
  `stop_overflow=0.1`；论文中的 20%/10% 分别是拥塞优化启动条件和单轮
  膨胀上限，不是最终 GP 门槛。
- 可视化：每 50 步一个快照，共 61 帧。

## 结果

| 数据集 | 初始化 | HPWL (M) | Overflow | 求解器时间 (s) | 首次 <=10% |
|---|---|---:|---:|---:|---:|
| adaptec1 | ePlace-IP | 75.433 | 6.477% | 160.12 | 1241 |
| adaptec2 | supplied .pl | 115.883 | 6.486% | 108.53 | 1372 |
| adaptec3 | supplied .pl | 326.346 | 6.480% | 181.48 | 1351 |
| adaptec4 | supplied .pl | 317.243 | 6.520% | 202.84 | 1300 |
| bigblue1 | supplied .pl | 103.900 | 6.486% | 96.41 | 1427 |
| bigblue2 | supplied .pl | 213.835 | 6.410% | 224.84 | 1312 |
| bigblue3 | supplied .pl | 707.435 | 7.162% | 481.70 | 2162 |
| bigblue4 | supplied .pl | 1497.324 | 6.511% | 1007.55 | 1485 |

八组全部满足 10% 参考门槛。总求解器时间 2463.47 秒，批处理墙钟时间约
2506 秒。结果是恢复后的第 2596 步最优可行状态；GIF 末帧是恢复前的第
3000 步轨迹，所以两者 HPWL 会略有差异。

## 汇总材料

- [综合 PDF 报告](ispd2005_bundle_evaluation.pdf)
- [原始结果 CSV](../research/bundle_hpwl/data/ispd2005_bundle_results.csv)
- [总收敛曲线 PNG](../research/bundle_hpwl/experiments/ispd2005_eight_benchmarks/figures/fig_ispd2005_aggregate_convergence.png)
- [结果汇总图 PNG](../research/bundle_hpwl/experiments/ispd2005_eight_benchmarks/figures/fig_ispd2005_result_summary.png)
- [实验协议](../research/bundle_hpwl/experiments/ispd2005_eight_benchmarks/protocol.md)
- [详细分析](../research/bundle_hpwl/experiments/ispd2005_eight_benchmarks/analysis.md)

## 动图

- [adaptec1](../research/bundle_hpwl/experiments/ispd2005_eight_benchmarks/animations/adaptec1.gif)
- [adaptec2](../research/bundle_hpwl/experiments/ispd2005_eight_benchmarks/animations/adaptec2.gif)
- [adaptec3](../research/bundle_hpwl/experiments/ispd2005_eight_benchmarks/animations/adaptec3.gif)
- [adaptec4](../research/bundle_hpwl/experiments/ispd2005_eight_benchmarks/animations/adaptec4.gif)
- [bigblue1](../research/bundle_hpwl/experiments/ispd2005_eight_benchmarks/animations/bigblue1.gif)
- [bigblue2](../research/bundle_hpwl/experiments/ispd2005_eight_benchmarks/animations/bigblue2.gif)
- [bigblue3](../research/bundle_hpwl/experiments/ispd2005_eight_benchmarks/animations/bigblue3.gif)
- [bigblue4](../research/bundle_hpwl/experiments/ispd2005_eight_benchmarks/animations/bigblue4.gif)

## 结论

静电密度可行性可以跨八个实例迁移，但 HPWL 质量和 lambda 尺度不能直接
迁移。bigblue3 的 10% 可行点明显更晚，说明固定 lambda 上限/增长速度缺乏
规模不变性。下一步应统一初始化，并按每个实例的 HPWL 与密度梯度 RMS
归一化 lambda，再进行 bundle 与 no-bundle 的配对实验。
