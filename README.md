# adaptec1 双向交叉求解实验

本目录把两套求解流程通过完整的 Bookshelf `.pl` 文件连接起来：

1. 用 `adaptive_pareto_soft` 的 `capacity_aware_seed` 生成初始解，再交给
   `combine/Alg/Alg/nsp_placer` 求解。
2. 把 `dataset/adaptec1/adaptec1.pl` 原始默认初始解交给
   `adaptive_pareto_soft` 求解。

两组输入、输出均由 `adaptive_pareto_soft` 的同一个精确评估器复算 HPWL、
legacy overflow 和 strict overflow，结果写入 `combine/results/<run>/summary.json`
和 `REPORT.md`。

注意：Adaptive 求解器在第 0 步会把 movable cell 投影到全局合法行边界。
因此第二组实验同时报告 `raw_input`（原始 `.pl`）与 `input`（投影后真正进入
迭代的解），求解增量以 `input` 为基线。

## 运行

先编译 `combine/Alg/Alg`，然后从仓库根目录执行：

```powershell
python combine/run_cross_experiments.py --solver combine/Alg/Alg/nsp_placer.exe
```

当前目录已生成并验证 Windows 二进制 `combine/Alg/Alg/nsp_placer.exe`。

Linux/macOS 下二进制通常为 `combine/Alg/Alg/nsp_placer`。脚本会自动寻找这
两个默认位置，也可以显式传入 `--solver`。快速单边验证可使用
`--skip-alg` 或 `--skip-adaptive`。
用 `--run-root <已有目录>` 可以补跑其中一边；被跳过但已经完成的结果会保留。

`combine/Alg` 新增参数：

- `--init-pl <file>`：加载完整的外部 Bookshelf 初始解；优先级高于 `--init`。
- `--prelegal-output <file>`：保存 legalizer 之前的全局放置结果，便于区分
  全局求解器和 legalizer 各自的影响。
- `--output <file>`：指定输出位置，避免在原始数据集旁写文件。

Adaptive 实验预算定义在 `combine/configs/a1_adaptive_from_default.json`；默认
最多 180 秒或 100 步，并关闭可视化，以便把时间集中在求解上。
