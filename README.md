# HardPlacer

`hardplacer` 是难题五的独立求解控制层：它复用工作区中 `gpplacer` 的 Bookshelf 解析和精确 Oracle，但不修改原项目。

## 方法

1. 先以独立的 density-first 容量搬移和面积相近交换建立低 overflow 安全锚点；此阶段只以精确 overflow 为目标，不引用 `runs/` 中的任何历史解。
2. 达到用户目标（默认示例为 `overflow <= 6%`）后，将该目标作为硬预算，再以预条件 HPWL 次梯度压缩线长；若锚点优于 6%，允许在不超过 6% 的空间内换取更低 HPWL。
3. 以每个 cell 独立的信赖域限制位移；拒绝的 cell 收缩，成功 cell 扩展。
4. 停滞时启用最多四个切平面的有限内存 Bundle；仍停滞时短暂加压，再通过搬移/交换修回复原锚点预算。

常规阶段以显式 `target_overflow_percent` 为硬上限；未配置该值时才锁定实际锚点 overflow。PDF 形式的平方密度默认仅记录诊断，可通过 `enforce_challenge_square_budget: true` 额外收紧。

## 安装与运行

在仓库根目录安装两个本地包：

```powershell
python -m pip install -e .\code -e .\code-hard
python -m hardplacer.cli solve --aux dataset\adaptec1\adaptec1.aux --config code-hard\configs\adaptec1_smoke.json --plot
```

若确实需要从外部候选解开始（该路径与默认独立开发运行分开记录）：

```powershell
python -m hardplacer.cli solve --aux dataset\adaptec1\adaptec1.aux --initial-placement runs\adaptec1\20260713T081941Z\solution.pl --target-overflow-percent 6 --plot
```

每个运行目录含 `solution.pl`、`iterations.csv`、`metadata.json`、`snapshots.npz` 和 `figures/`。图表为 PNG（300 DPI）及 PDF。

## 验收边界

本地仅有 adaptec1/adaptec3。正式八案例比较仍使用同线程数、同 baseline 的 manifest；没有其余六例和指定 baseline 时，任何本地结果都不代表正式验收通过。
