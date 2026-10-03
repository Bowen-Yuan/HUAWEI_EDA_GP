# HardPlacer：难题五独立 density-first 原型报告

**状态：已完成独立单案例验证；未完成八案例正式验收。** 本报告只引用 `code-hard/runs_density_first/`，不引用 `runs/`、`runs_anchor/`、`runs_density_repair/` 或 `dai-code/` 的既有结果。

## 结论

本轮采用“先构造低 overflow 安全锚点，再在 `overflow <= 6%` 的硬约束下压 HPWL”。6% 是上限，不会因为锚点恰好为 0% 而被错误收紧为严格零约束。

adaptec1 的独立运行结果：HPWL 从 **1,965,941,020** 降至 **1,776,703,581**（下降 **9.63%**）；最终线性 capacity overflow 为 **2.7303%**，始终小于 6% 上限。求解循环耗时 11.57 s，共 160 次迭代。

## overflow 口径核对

题目、文献和现有代码都基于精确 cell-bin overlap，但汇总指标不完全一致，不能直接拿数值横比。

| 指标 | 定义 | 本实现用途 |
| --- | --- | --- |
| 线性 capacity overflow (%) | `100 * Σ max(occ_b - rho_t * available_b, 0) / Σ available_b` | 本次用户指定的 `<= 6%` 硬约束。`available_b` 扣除固定块与 blockage。 |
| 现有 `density_official` | `Σ available_b * max(occ_b/available_b - rho_t, 0)^2` | 保留的 capacity 加权平方诊断。 |
| 题目形式的平方密度 | `Σ max(rho_b-rho_t, 0)^2`，`rho_b=overlap_area/geometric_bin_area` | 记录为 `challenge_square_density`，包含固定 instance。 |

因此，题目形式的平方密度不是百分比，也不与线性 capacity overflow 同量纲。本次运行中它由 163.436 增至 629.658：这是允许在 6% 的可行余量内降低 HPWL 的结果，默认不将其当作第二个硬约束。若需要更保守的双指标模式，可在配置中设置 `enforce_challenge_square_budget: true`。

此前观察到的很大 overflow 还存在输入起点问题：Bookshelf 原始 `.pl` 中 adaptec1 的 movable cell 坐标整体在 core 外；直接投影到 core 边缘会造成非物理的拥塞。本实现使用独立的容量均衡 seed，不读取任何历史 `.pl` 结果。

## 算法

1. **安全锚点**：按每个 bin 的剩余目标容量分配 movable cell；若给定初始解超出上限，再以 overfull bin 到 slack bin 的搬移、面积相近交换修复。未达到目标上限时直接报错，不输出伪可行解。
2. **约束压缩**：使用精确 HPWL 次梯度及节点连接度预条件；每个分量由独立信赖域裁剪。只有 HPWL 严格下降且线性 overflow 不超过 6% 的候选会被接受。
3. **停滞处理**：保留固定数量的次梯度切平面做轻量 bundle 聚合；必要时短暂放宽探测，再通过密度修复回到 6% 以内才接受。

没有使用 WA、LSE、Moreau 或密度平滑近似。

## 可复现实验

```powershell
python -m pip install -e .\code -e .\code-hard
python -m hardplacer.cli solve --aux dataset\adaptec1\adaptec1.aux `
  --config code-hard\configs\adaptec1_smoke.json --max-iterations 160 --plot
```

运行目录：[20260717T072328Z](runs_density_first/adaptec1/20260717T072328Z)

| 指标 | 初始 / 锚点 | 最终 | 约束检查 |
| --- | ---: | ---: | --- |
| HPWL | 1,965,941,020 | 1,776,703,581 | 下降 189,237,439（9.63%） |
| 线性 overflow | 0.0000% | 2.7303% | 小于 6.0000% |
| 题目形式平方密度 | 163.436 | 629.658 | 仅诊断，默认不约束 |
| 迭代 / 求解循环 | - | 160 / 11.57 s | 单案例结果 |

![收敛曲线](runs_density_first/adaptec1/20260717T072328Z/figures/convergence.png)

![布局与密度分布](runs_density_first/adaptec1/20260717T072328Z/figures/placement_and_density.png)

## 验证边界

- 单元测试覆盖：6% 预算不被锚点错误收紧、可选平方诊断硬约束、固定 cell 不动、分量信赖域边界和数值有限。
- 当前工作区仅具备 adaptec1/adaptec3，且没有指定 baseline，不能声称已达到八案例平均 HPWL、runtime 或官方 overflow 的竞赛验收指标。

相关文件：[算法实现](hardplacer/solver.py)、[配置](configs/adaptec1_smoke.json)、[测试](tests/test_safe_anchor.py)。
