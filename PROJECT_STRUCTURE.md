# GPPlacer 项目结构

`gpplacer` 是面向 ISPD2005 Bookshelf 数据的原生非光滑 Global Placement 原型。它直接计算 HPWL 的 max/min 形式与 cell-bin 精确矩形交叠，不使用 WA、LSE 或 Moreau envelope 平滑近似。

## 目录职责

| 目录 | 职责 |
| --- | --- |
| `gpplacer/io` | 读取 `.aux/.nodes/.nets/.wts/.pl/.scl`，将 Bookshelf 左下角坐标转换为内部中心坐标，并写回完整 `.pl`。 |
| `gpplacer/model` | 类型安全的数据模型、CSR 索引与统一配置。 |
| `gpplacer/oracle` | 原始 HPWL、精确 cell-bin overlap、density grid 容量和次梯度。 |
| `gpplacer/multilevel` | 确定性的重边匹配与容量感知初始布局。 |
| `gpplacer/solver` | Explore、受限记忆近端 Bundle、局部 Refine、阶段切换与恢复。 |
| `gpplacer/metrics` | 每轮标量日志、快照、运行环境和资源统计。 |
| `gpplacer/visualization` | 从已保存日志生成收敛、活跃集、资源、布局和 density 图。 |
| `gpplacer/benchmark` | 八案例 baseline 比例验收，不把缺失 baseline 伪造成结论。 |

## 坐标与密度口径

Bookshelf `.pl` 的 `(x, y)` 是对象左下角；优化变量为对象中心。每个 pin 坐标为：

```text
pin = centre_of_instance + pin_offset_from_nets
```

`.scl` 的合法行段决定 bin 的可用面积；固定对象面积从该容量中扣除。默认 `rho_target=0.8`，优化器使用：

```text
D_linear = sum_b available_area_b * max(rho_b - 0.8, 0)
```

同时记录正式难题5公式使用的平方超限诊断值，供消融和合规报告使用。该诊断值不改变用户确认的线性搜索模型。

## 主流程

```text
Bookshelf DB
  -> capacity-aware coarse seed
  -> preconditioned subgradient Explore
  -> limited-memory proximal Bundle
  -> optional local Refine
  -> complete .pl + structured logs + figures
```

每个候选点都由原始目标验收。出现 NaN、越界、异常目标或子问题失败时，求解器恢复最近 serious center 并缩小步长或近端范围。

## 常用命令

在 `D:\codex_project\HUAWEI_EDA` 下运行：

```powershell
$env:PYTHONPATH = "$PWD\code"
python -m gpplacer.cli inspect --aux dataset\adaptec1\adaptec1.aux
python -m gpplacer.cli solve --aux dataset\adaptec1\adaptec1.aux --max-iterations 20 --threads 1 --plot
python -m gpplacer.cli evaluate --aux dataset\adaptec1\adaptec1.aux --placement runs\adaptec1\<run_id>\solution.pl
python -m gpplacer.cli visualize-run --run-dir runs\adaptec1\<run_id>
python -m gpplacer.cli benchmark --manifest baseline_manifest.json
```

## 运行产物与可视化

每次运行写入 `runs/<case>/<UTC timestamp>/`：

- `metadata.json`：数据规模、配置、线程数、目标摘要和运行状态；
- `iterations.csv`：HPWL、overflow、阶段、serious/null 信息、活跃集变化与 RSS；
- `snapshots/*.npz`：初始/周期性/最佳位置和 density 状态；
- `figures/*.pdf`、`figures/*.png`：收敛、活跃集、步长、阶段、内存、Oracle 耗时、布局采样，以及初始/最佳 density 对比图；
- `solution.pl`：完整 Global Placement 输出。

图表从日志读取而不重复实现求解逻辑。数值图采用色盲友好配色，同时导出 PDF 与 300 DPI PNG；百万对象布局仅绘制确定性样本并使用栅格化。

## 可读性约定

- 公共函数和类都带类型标注、docstring 和明确的坐标单位。
- 数学内核与阶段控制、输入输出、可视化严格分离。
- 参数只出现在 `SolverConfig`；不在算法逻辑中引入案例特例或无说明常数。
- 每个性能内核应与小规模参考行为通过测试一致。

## 验收边界

当前数据只有 adaptec1 和 adaptec3，可用于开发、回归和冒烟。正式难题5要求八个 ISPD2005 case、同一 baseline、相同 CPU 并行数下比较 HPWL、overflow、runtime 和迭代数；`benchmark` 模块为此保留严格八案例 manifest 接口，缺少其余数据或 baseline 时不会输出“通过正式验收”。
