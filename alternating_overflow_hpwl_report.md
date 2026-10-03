# 交替式 Overflow / HPWL 独立开发说明

代码入口：[alternating_overflow_hpwl.py](alternating_overflow_hpwl.py)。它不导入 `hardplacer`，不修改已有 `runs`，并将输出固定写入 `code-hard/runs_alternating/`。

## 阶段顺序

执行顺序固定为 `density_init → hpwl_1 → density_1 → hpwl_2 → density_2 ...`。每个阶段的接受准则严格分离：

1. `density_init`：只接受精确线性 overflow 下降的候选，直至达到 6%。次梯度停滞时，使用不看 net 的批量疏散跨越 bin 台阶，并由精确 oracle 验收。
2. `hpwl_k`：只接受 HPWL 严格下降且 overflow 不超过 6% 的候选。
3. `density_k`：再次只优化 overflow，不以 HPWL 作目标；仅设置该轮 HPWL 终点的 0.25% 守护线，防止无界回退。

如果第一阶段不能把 overflow 降至 6%，程序直接报错。

## 已验证运行

种子是低 HPWL、高 overflow 的历史输入文件，仅作为输入而不被改写：

```powershell
$env:PYTHONPATH = '.\code'
python .\code-hard\alternating_overflow_hpwl.py `
  --aux .\dataset\adaptec1\adaptec1.aux `
  --seed-placement .\runs\adaptec1\20260713T081941Z\solution.pl --cycles 2
```

运行目录：[20260717T073740Z](runs_alternating/adaptec1/20260717T073740Z)。

| 指标 | 种子 | 最终 |
| --- | ---: | ---: |
| HPWL | 55,403,559 | 1,494,433,533 |
| 线性 capacity overflow | 51.2802% | 5.8677% |
| 求解循环 | - | 65.49 s，192 条阶段记录 |

阶段终点：`density_init` 到 5.9648%，`hpwl_1` 回到 6.0000%，`density_1` 到 5.8776%，`hpwl_2` 回到 6.0000%，`density_2` 到 5.8677%。收敛与布局可视化：[收敛曲线](runs_alternating/adaptec1/20260717T073740Z/figures/convergence.png)、[布局与密度](runs_alternating/adaptec1/20260717T073740Z/figures/placement_density.png)。

密度修复会显著损失初始 seed 的低 HPWL；这是将 51.28% 降到 6% 以下的代价。之后的 HPWL 阶段在 6% 约束下继续压线长，后续 density 阶段再小幅回压 overflow。不要把该单案例结果与 `hardplacer` 的容量均衡初始化结果混合比较。
