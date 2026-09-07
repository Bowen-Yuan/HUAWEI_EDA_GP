# 00 — 开始工作前必须知道

## 固定位置

- 工作目录：`D:\codex_project\HUAWEI_EDA\non-smooth`
- 数据集默认目录：`D:\codex_project\HUAWEI_EDA\alg-electronic\ispd2005`
- GitHub：`https://github.com/Bowen-Yuan/HUAWEI_EDA_GP.git`
- 当前开发分支：`nonsmooth-gp-v1`
- 架构基线：registry-driven microkernel；实际提交以 `git rev-parse HEAD` 为准。
- 默认开发 case：`adaptec1`
- 完整 case：`adaptec1..4`、`bigblue1..4`；除非明确要求，不要日常全量运行八个 case。

数据集和外部 checkpoint 不属于仓库，不能复制、修改、删除或提交。

## 每次接手的第一组只读命令

```powershell
Set-Location 'D:\codex_project\HUAWEI_EDA\non-smooth'
git branch --show-current
git rev-parse HEAD
git status --short
git log -5 --oneline
rg -n "selected\.pl|snapshot|write_bookshelf_placement|output_dir|trajectory" framework modules scripts
```

先确认真实工作树，保留用户已有的未提交修改。不要根据旧对话或规划文件重建一套平行框架。

## 当前最短开发路径

1. 先读本目录全部 Markdown。
2. 根据任务定位权威实现：共享数值内核在 `framework/kernel`；优化阶段的权威代码在 `modules/<stage>/code`。
3. 只改完成任务所需的最小范围。
4. 编译并先跑数值契约测试或 1 轮 smoke。
5. 正式实验默认用 `adaptec1`，从相同输入 checkpoint 和相同 SHA-256 开始。
6. 检查实验目录的 retention 合约和调用链记录。
7. 同步更新本目录及 `CHANGELOG.md`。
8. 检查 `git diff`、提交，并推送到 `origin/nonsmooth-gp-v1`；不得提交数据集、二进制或常规实验产物。

## 当前构建入口

推荐 CMake：

```powershell
cmake -S . -B build
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

本机历史上也使用过 `D:\MingGW\ucrt64\bin\g++.exe` 手工编译；这种方式只适合诊断，不应取代 `CMakeLists.txt` 作为权威构建描述。

主要目标：

- `nsgp`：微内核 pipeline、audit、batch 和 `global_view_gp` 的 `lab` 入口。
- `nsgp_tests`：exact HPWL/density 小型数值契约测试。
- `nsgp_legacy_stage`：位于 `modules/historical_exact_replay` 的兼容目标，会产生较多中间文件。
- `nsgp_historical_homotopy`：历史 DCT/Poisson 复现工具，不能进入正式 non-smooth pipeline。

## 最常用命令

```powershell
# 查看模块
build\Release\nsgp.exe list-modules

# 权威指标审计；输出 overflow_percent
build\Release\nsgp.exe audit --case adaptec1 --placement <placement.pl>

# 模块组合 pipeline；默认同样是 metrics-only
build\Release\nsgp.exe run --case adaptec1 --pipeline framework\params\pipelines\smoke.json --threads 1

# V3 metrics-only 实验入口
build\Release\nsgp.exe lab --case adaptec1 --iterations 50 --optimizer adam --active-ensemble --bundle

# retention 结果结构检查
powershell -ExecutionPolicy Bypass -File tests\test_v3_retention_contract.ps1 -ExperimentDirectory <experiment_dir>
```

## 红线摘要

- 不得用 WA、LSE、DCT/Poisson 或其他平滑 surrogate 替代最终 exact HPWL / exact density 评价。
- 内部 overflow 是比例，面向人的输出必须明确显示为百分比；例如 `0.07` 显示为 `7%`。
- overflow 达标不等于 legalized。本项目当前不做 row/site snapping、overlap removal、legalization 或 detailed placement。
- 固定宏消耗容量且不能移动；`terminal_NI` 不消耗物理容量。
- 默认只保存数值实验所必需的指标和说明，不堆积 placement、snapshot、缓存和调试 dump。
- 每次实验必须能回答：从哪里开始、经过哪些模块、每步参数是什么、读取了哪个前序结果、最终 exact 指标是多少。
