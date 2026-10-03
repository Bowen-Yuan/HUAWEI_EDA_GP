# 25% Overflow 上限下的交替优化实验

本实验使用独立脚本 [alternating_overflow_hpwl.py](alternating_overflow_hpwl.py)，不修改 `hardplacer` 的既有运行结果。输入仍为 adaptec1 的低 HPWL、高 overflow 历史 seed：HPWL 为 `55,403,559`，线性 capacity overflow 为 `51.2802%`。

命令：

```powershell
$env:PYTHONPATH = '.\code'
python .\code-hard\alternating_overflow_hpwl.py `
  --aux .\dataset\adaptec1\adaptec1.aux `
  --seed-placement .\runs\adaptec1\20260713T081941Z\solution.pl `
  --overflow-cap-percent 25 --cycles 3 `
  --output-root .\code-hard\runs_alternating_25
```

运行目录：[20260720T144806Z](runs_alternating_25/adaptec1/20260720T144806Z)。

## 结果

| 指标 | 初始 seed | 最终 |
| --- | ---: | ---: |
| HPWL | 55,403,559 | 1,020,438,325 |
| 线性 capacity overflow | 51.2802% | 23.4976% |
| 阶段记录数 | - | 400 |
| 求解循环时间 | - | 113.15 s |

严格交替顺序和阶段终点如下：

| 阶段 | HPWL 终点 | overflow 终点 |
| --- | ---: | ---: |
| `density_init` | 1,119,761,399 | 20.4129% |
| `hpwl_1` | 1,066,236,944 | 25.0000% |
| `density_1` | 1,066,334,101 | 23.5710% |
| `hpwl_2` | 1,041,755,733 | 25.0000% |
| `density_2` | 1,041,880,898 | 23.5839% |
| `hpwl_3` | 1,020,265,416 | 25.0000% |
| `density_3` | 1,020,438,325 | 23.4976% |

每个 HPWL 阶段都使用 `overflow <= 25%` 硬接受条件；每个后续 density 阶段只接受线性 overflow 下降，但以该轮 HPWL 终点的 `0.25%` 作为保护边界，因此会出现很小的 HPWL 回升。

![完整交替收敛曲线](runs_alternating_25/adaptec1/20260720T144806Z/figures/convergence.png)

![seed 与最终布局/密度](runs_alternating_25/adaptec1/20260720T144806Z/figures/placement_density.png)

## 与 6% 实验的解释

较早的 6% 交替实验最终为 HPWL `1.494B`、overflow `5.8677%`。本次 25% 实验的 HPWL `1.020B` 低约 31.7%，说明放宽 density 约束确实允许保留更多连接结构；但两次运行的交替轮数不同（本次 3 轮、此前 2 轮），不应将该差值解释为严格的单变量对照。

两组结果都仍远高于初始 55.4M HPWL，原因是第一个 density-only 批量疏散阶段必须打散中心拥塞。要进一步改善，需要采用更好的 `4%--25%` 初始解构造器，而不是仅延长交替轮数。

> 口径说明：以上 overflow 是当前共享 `gpplacer` oracle 的线性 capacity overflow。该 oracle 对完全被 fixed macro 遮挡且 `available_area=0` 的 bin 存在已知漏罚问题，因此这些数值不能作为最终物理可行性证明；修复该 oracle 后应重新运行本实验。
