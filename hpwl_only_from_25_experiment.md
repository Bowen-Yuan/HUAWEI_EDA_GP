# 从 25% 解出发的 HPWL-only 长跑实验

本实验与交替求解器分离。输入是 25% overflow 交替实验的最终 placement；overflow 仍由相同的 exact oracle 记录，但不参与方向、步长、信赖域、回溯或接受条件。

运行入口：[hpwl_only_from_placement.py](hpwl_only_from_placement.py)。

```powershell
$env:PYTHONPATH = '.\code'
python .\code-hard\hpwl_only_from_placement.py `
  --aux .\dataset\adaptec1\adaptec1.aux `
  --placement .\code-hard\runs_alternating_25\adaptec1\20260720T144806Z\solution.pl `
  --iterations 1000 --max-seconds 600 `
  --output-root .\code-hard\runs_hpwl_only_from_25
```

运行目录：[20260720T150001Z](runs_hpwl_only_from_25/adaptec1/20260720T150001Z)。

| 指标 | 起点（25% 实验最终解） | HPWL-only 最终解 |
| --- | ---: | ---: |
| HPWL | 1,020,438,325 | 50,147,498 |
| 线性 capacity overflow | 23.4976% | 36.1372% |
| 迭代数 | - | 1,000 |
| 接受 HPWL 步数 | - | 756 |
| 求解循环时间 | - | 148.97 s |

HPWL 降低约 95.1%，而 overflow 上升 12.64 个百分点。曲线在约第 500 轮后接近平台，说明从当前起点出发，纯 HPWL 次梯度已经成功重新聚集 cell 并接近其无约束低线长区域；这不是 HPWL 难以优化，而是与 overflow 约束的明显冲突。

![HPWL-only 收敛与 overflow 诊断](runs_hpwl_only_from_25/adaptec1/20260720T150001Z/figures/hpwl_only_convergence.png)

该最终解不能作为可行 placement 使用：overflow 没有任何硬约束，且当前共享 oracle 对完全 blocked bin 有已知漏罚问题。它的用途仅是量化“去掉 density 约束后，HPWL 还可以下降多少”。
