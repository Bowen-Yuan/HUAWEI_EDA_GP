# epsilon-active 模块化重构规划输入

> 用途：将本文与后续必要的核心代码一并交给高强度推理模型，用于规划一个**新的、轻量、灵活、面向数值优化实验**的代码项目。本文描述当前项目的核心数学契约、主要接口、实际阶段调用链、已知耦合与重构要求；它不是详细重构设计，也不要求在原仓库内直接改造。
>
> 本文基于 2026-09-06 的本地工作树。该工作树包含未提交代码和大量未跟踪实验产物，因此应把本文理解为“当前可运行研究状态”的摘要，而不是某个 Git commit 的精确说明。

## 0. 给网页版 GPT 的最高优先级指令

请直接根据本文输出一份**可以落地编写代码的新项目方案**。方案必须首先服从以下四项要求；如果后文中的现状描述、历史做法或一般软件工程惯例与它们冲突，以本节为准。

### 0.1 必须研究和实现非光滑优化，不得把问题改写成平滑 GP

新项目的主目标是从原始 nonsmooth Global Placement 模型出发，直接处理：

```text
min W(x, y) + lambda * D(x, y)

W = exact weighted pin-offset HPWL（max-min）
D = exact rectangle/bin overlap 形成的分段非光滑 density penalty
```

硬约束如下：

- 不得使用 weighted-average wirelength、log-sum-exp、Moreau envelope 等光滑近似来替代原始 HPWL；
- 不得用平滑 density surrogate 替代最终的 exact rectangle/bin overlap 评价；
- 可以研究 subgradient、epsilon-active piece selection、bundle/proximal 思想、exact breakpoint/coordinate search、exact-audited discrete move、capacity assignment 等非光滑求解手段；
- 任何内部方向、辅助场或候选生成器都必须与权威 exact objective/metrics 区分，所有 checkpoint 必须重新进行 exact HPWL 与 exact overflow 审计；
- 历史 H219/H221 使用的 DCT/Poisson electrostatic homotopy 只能作为明确标注的历史/可选实验模块，**不能默认被视为满足“不使用光滑化近似”要求的主方案**。若规划者希望保留它，必须单独论证“只作为最终 `mu=0` 前的搜索辅助”是否符合题意；否则应给出纯非光滑替代路径。不要让该历史模块弱化新项目的非光滑定位。

上述要求来自 `D:\codex_project\HUAWEI_EDA\第148期 - EDA专题第四期.pdf` 第五题“面向大规模 Global Placement 的非光滑优化建模与求解”。该题要求兼顾结构保真性、数值稳定性和大规模可扩展性，并在 ISPD2005 的 8 个 case 上**同时**满足：

1. 不使用光滑化近似；
2. 无 case 崩溃；
3. 平均 HPWL 相对选定 baseline 劣化小于 1.0%；
4. overflow 达到与 baseline 同量级；
5. 相同 CPU 并行数下，runtime 不超过 baseline 的 2 倍；
6. 题面第（6）项原文为“优化的迭代次数 5倍”。原页没有完整写出比较关系，正式方案必须把它标为待澄清项，并给出一个明确、可配置、可审计的迭代预算解释，而不能悄悄猜测原意。

验证对象就是本地数据集中的 8 个目录：`adaptec1`–`adaptec4`、`bigblue1`–`bigblue4`。方案必须说明 baseline 的选择、相同并行度的计时方法，以及如何一次性汇总八个 case 的全部门槛；不能只针对 adaptec1 设计。

### 0.2 必须轻量、灵活，直接服务数值实验

这不是企业软件或通用优化平台。请设计一个研究者可以直接阅读、修改、编译、运行的简单数值实验项目：

- 一个薄的顶层架构负责数据、权威指标和顺序调度；
- 每种优化算法是一个统一输入/输出的小模块；
- 用简单配置文件排列组合模块；
- 优先少量清楚的数据结构和普通函数；若函数表/简单 factory 足够，就不要设计复杂插件系统；
- 不要引入微服务、RPC、数据库、依赖注入、大型 registry、工作流引擎或过度抽象的 DSL；
- 不要为了“软件工程完整性”制造许多空层、adapter、manager、service；
- 第一版只需直接实现能复现和继续研究核心数值链的方案。

### 0.3 数据必须直接从指定目录读取

默认且首要的数据根目录固定为：

```text
D:\codex_project\HUAWEI_EDA\alg-electronic\ispd2005
```

每个 case 从以下形式读取 Bookshelf 文件：

```text
<dataset_root>/<case>/<case>.aux
<dataset_root>/<case>/<case>.nodes
<dataset_root>/<case>/<case>.nets
<dataset_root>/<case>/<case>.pl
<dataset_root>/<case>/<case>.scl
<dataset_root>/<case>/<case>.wts    # 若存在
```

该目录已包含 `adaptec1`–`adaptec4` 和 `bigblue1`–`bigblue4`。新项目不要复制这批数据，不要把数据塞进代码仓库，也不要默认从旧 `epsilon-active/experiments` 目录寻找 benchmark。可以允许命令行覆盖 `dataset_root`，但默认配置和示例必须直接使用上述绝对路径。

### 0.4 目录只需要“顶层架构 + 模块”，每个模块明确三分

首选的最小物理结构是：

```text
new-project/
├─ framework/                     # 顶层架构：Problem/Layout、Bookshelf I/O、exact metrics、pipeline runner
│  ├─ code/
│  ├─ params/                     # 全局默认值、八 case 清单、pipeline 配置
│  └─ results/                    # 顶层汇总；默认不提交大文件
└─ modules/
   ├─ hpwl_adam/
   │  ├─ code/                    # 该模块的实现
   │  ├─ params/                  # 该模块可提交 Git 的参数 preset
   │  └─ results/                 # 该模块的中间布局、CSV、日志、快照
   ├─ exact_recovery/
   │  ├─ code/
   │  ├─ params/
   │  └─ results/
   └─ ...
```

目录名允许做很小的调整，但两层逻辑不能变：一个简单顶层架构文件夹；一个专门的模块文件夹；并且**每个模块都清楚分成 `code/`、`params/`、`results/` 三部分**。`code/` 与 `params/` 应可直接上传 GitHub；`results/` 目录结构可以保留，但大体积内容默认 `.gitignore`。

## 1. 结论：这个方向可行

当前项目已经存在若干可复用的算法边界：Bookshelf 读写、精确 HPWL、精确 rectangle/bin overlap、优化器、lambda 控制、递归二分、exact recovery、equal-shape swap 等都有相对独立的 C++ 接口。真正造成庞大和繁杂的主要原因是：

1. 一个 `PlaceConfig` 聚合了几乎所有历史实验开关；`src/main.cpp` 约有 210 个 CLI 解析分支。
2. `global_place()` 内部写死了许多算子的顺序，而主要成果链又通过 PowerShell 启动多个进程、传递不同语义的 `.pl` checkpoint 来绕开这个固定顺序。
3. 算法代码、参数组合、实验 checkpoint、CSV、图片、PDF 和研究记录长期共处一个仓库；当前 `experiments/` 下约有 477 个一级实验目录。
4. 主链横跨当前 exact core 与一个归档的 DCT/Poisson homotopy 可执行文件，源码和依赖边界并不统一。

因此，重构成“共享问题/指标核心 + 统一模块协议 + 可声明组合的 pipeline + 仓库外运行产物”是自然且低风险的方向。关键不是建立大型软件框架，而是保住数值契约、实验可复现性和快速替换阶段的能力。

## 2. 项目性质与不可误解的边界

### 2.1 项目是什么

`epsilon-active` 是 Huawei EDA Challenge 5 / ISPD 2005 Bookshelf benchmark 上的全局布局数值优化研究代码。主要研究问题是：如何在精确非光滑 HPWL 与精确网格重叠密度的约束下得到低 HPWL、低 overflow 的连续坐标布局。

它输出 global placement，不是完整 P&R 软件。核心实现使用 C++17 和 OpenMP，线程数限制为 1–40。

### 2.2 当前严格评价模型

对 net `e`，pin 坐标包含 Bookshelf pin offset：

```text
HPWL_e = max_i(x_i + offset_x_i) - min_i(x_i + offset_x_i)
       + max_i(y_i + offset_y_i) - min_i(y_i + offset_y_i)

HPWL = sum_e weight_e * HPWL_e
```

对 bin `b`：

```text
occupancy_b = movable rectangles 与 bin 的精确相交面积
            + physical fixed macros 与 bin 的精确相交面积

rho_b = occupancy_b / bin_area

density_energy = sum_b 0.5 * bin_area * max(rho_b - target_density, 0)^2

overflow = sum_b bin_area * max(rho_b - target_density, 0)
           / total_movable_area
```

必须保留的语义：

- fixed macro 消耗 bin 容量且不能移动；
- `terminal_NI` 是非物理端口，不消耗 bin 容量；
- epsilon-active 只选择非光滑函数的搜索方向，不改变精确 HPWL、精确 occupancy 或报告指标；
- overflow 低于阈值不等于几何合法；
- 当前项目不做 row/site snapping、overlap removal、legalization、detailed placement 或 legality gate。

### 2.3 DCT/Poisson 辅助项的边界

当前仓库根目录下的 `epsilon_active_core` 声称并确实以 exact HPWL + exact overlap 为核心，不包含 Poisson/DCT 路径。但成果主链中的 H219/H221 会调用：

```text
experiments/h72_homotopy_electrostatic/homotopy_main.cpp
experiments/h72_homotopy_electrostatic/homotopy.exe
```

该归档程序把 exact HPWL、exact overlap 与一个阶段性 DCT/Poisson electrostatic 辅助项结合，并强制最后一个阶段 `mu=0`。它通过 `electric.h`、`bookshelf.h` 使用相邻 `dreamplace-cpp` 的实现；这些依赖不在当前 `epsilon-active` 核心库内。

检查时 `homotopy_main.cpp` 已被 Git 跟踪，但 `homotopy.exe` 与当前复现脚本 `scripts/run_h375_transfer.ps1` 尚未被跟踪；后者是当前主链最直接的事实来源，却还不是稳定的 commit provenance。

因此，新的架构若迁移“DCT/Poisson homotopy”，必须把它视为一个明确的**历史复现/可选对照模块**，并明确其依赖与最终 exact 审计，不能继续让它以归档二进制的形式隐含在实验目录里。面向第五题的默认主 pipeline 必须是非光滑路线；是否允许在主路线中使用最终会归零的光滑辅助项，需要规划者明确判定，不能通过改名规避“不使用光滑化近似”的限制。H252 及之后的历史主链回到当前 exact core。

## 3. 当前核心数据与 I/O 契约

核心类型位于 `include/epsilon_active/types.hpp`。

### 3.1 `Database`

`Database` 同时保存静态问题和可变布局状态：

- 静态拓扑：`nodes`、`pins`、`nets`、`rows`、node-to-pin 索引、net weight；
- 分类索引：`movable_ids`、`fixed_ids`；
- 版图边界：`xl, yl, xh, yh`；
- 面积与来源：`movable_area`、`benchmark_base`、`raw_pl_path`；
- 可变状态：每个 `Node` 中的 `x, y, orientation`。

内存中的 `Node.x/y` 是单元中心坐标。Bookshelf `.pl` 文件中的坐标是左下角；读取时加半宽/半高，写出时再减半宽/半高。

新项目的默认 `benchmark_base` 必须由
`D:\codex_project\HUAWEI_EDA\alg-electronic\ispd2005\<case>\<case>`
构造，不再沿用历史实验中出现过的其他 dataset root。

### 3.2 Bookshelf 接口

`include/epsilon_active/bookshelf.hpp`：

```cpp
Database read_bookshelf(const std::filesystem::path& benchmark);
void load_bookshelf_placement(Database& db, const std::filesystem::path& path);
void initialize_center_gaussian(Database& db, std::uint64_t seed, Real sigma_ratio);
void write_bookshelf_placement(const Database& db, const std::filesystem::path& path);
```

`read_bookshelf()` 读取 `.nodes/.nets/.pl/.scl` 和可选 `.wts`，建立静态问题与原始坐标。`load_bookshelf_placement()` 只覆盖节点位置/方向，因而当前阶段间的实际传递格式是完整 Bookshelf `.pl`。

### 3.3 当前 checkpoint 名称并非同义词

- `global.pl`：当前可执行程序最终选择的输出；若存在 overflow 可行点，选最低 HPWL 可行点，否则选最低 overflow 点。
- `last.pl`：主 GP 轨迹的最后状态，专门用于阶段 continuation；它与最终选择器可能不同。
- `best.pl`：归档 homotopy 程序在最终 `mu=0` 阶段选择出的 checkpoint。
- `snapshots/*.pl`：迭代快照。

历史 H222 曾发现 H221 同名 metrics row 与 snapshot 坐标处在 optimizer update 的两侧，导致重放审计不一致；后来使用 `best.pl` 修正。新项目应让布局状态、迭代编号、参数摘要和由该布局重新计算的指标形成一个不可混淆的模块输出记录，而不是依赖文件名约定。

## 4. 当前核心代码接口与实现职责

以下是规划新项目时应优先阅读和迁移的核心，而不是所有历史实验代码。

### 4.1 权威指标/oracle

#### `ExactHpwl` — `include/epsilon_active/hpwl.hpp`, `src/hpwl.cpp`

```cpp
Real ExactHpwl::evaluate(
    Real epsilon,
    Real active_power,
    int degree_limit,
    std::vector<Real>* grad_x,
    std::vector<Real>* grad_y);
```

- 始终返回精确 weighted pin-offset HPWL。
- 若请求 gradient，则对每条 net 的 min/max 邻域生成归一化 epsilon-active 方向，再按 node 聚合。
- `degree_limit` 只限制哪些 net 贡献搜索方向，不影响返回的精确 HPWL。

#### `ExactOverlapDensity` — `include/epsilon_active/density.hpp`, `src/density.cpp`

```cpp
DensityMetrics evaluate(...);
DensityMetrics evaluate_with_prices(...);
DensityMove evaluate_move(int node_id, Real new_x, Real new_y) const;
DensityMove evaluate_group_move(const std::vector<DensityNodeMove>& moves) const;
void commit_move(const DensityMove& move);
```

- 构造时建立 bin 几何和 fixed occupancy。
- `evaluate()` 从 fixed occupancy 开始重建 movable occupancy，计算 exact energy/overflow/max density；可同时生成 overlap active-piece 方向。
- `evaluate_with_prices()` 只用 price 调制搜索方向，报告指标仍来自无权重 exact oracle。
- `evaluate_move/evaluate_group_move/commit_move` 为 recovery、swap、transport 等离散候选提供增量式精确容量审计。

这两个 oracle 应成为新顶层核心的权威评价入口，优化模块不能各自定义不同版本的 HPWL/overflow。

### 4.2 通用一阶优化部件

#### `Optimizer` — `include/epsilon_active/optimizer.hpp`

```cpp
class Optimizer {
public:
    virtual void reset(std::size_t dimensions) = 0;
    virtual void compute_delta(
        const std::vector<Real>& gradient,
        Real learning_rate,
        Real maximum_delta,
        std::vector<Real>& delta) = 0;
};
```

实现包含 Adam、AMSGrad、AdaGrad、HeavyBall、SGD。这里已经是较好的统一接口。

#### `LambdaController` — `include/epsilon_active/lambda_controller.hpp`

根据初始 HPWL/density gradient scale 初始化密度 multiplier，并支持 DreamPlace、trajectory、ratio 三类更新策略。实验已经表明，初始化或 grid resolution 变化后直接继承无量纲控制值会失真，因此不同模块/分辨率交接时必须明确是重算、映射还是继承控制状态。

#### `accept_exact_batch()` — `include/epsilon_active/batch_acceptance.hpp`

对候选位移做 exact HPWL/overflow 回溯审计。可允许可行前的 overflow 改善，也可施加 overflow cap 与相对 HPWL 损伤预算。

### 4.3 当前算法阶段接口

这些函数几乎都采用同一种隐式模式：传入并原地修改 `Database`，共享一个与该 `Database` 绑定的 `ExactOverlapDensity`，返回各自的 `Stats`。

| 阶段 | 当前公开函数 | 核心实现含义 |
|---|---|---|
| recursive bisection | `recursive_hypergraph_bisection(Database&, ExactOverlapDensity&, BisectionConfig)` | 按 region 容量递归切分 hypergraph；支持 position seed、FM、leaf capacity assignment、surplus-only、HPWL-guided leaf 等 |
| coarse flow | `coarse_capacity_flow(...)` | 在粗网格上从过载源向容量目的地规划 flow，可尝试 exact anchor、connected commodity 和 joint assignment |
| transport | `transport_excess_to_capacity(...)` | 多种单元/分组容量转移、auction/Hilbert/nearest destination、atomic group、identity exchange |
| density coordinate | `coordinate_descent_overlap(...)` | exact overlap 坐标/断点下降，可含 node、net block、容量 assignment 和层次 cluster 候选 |
| exact recovery | `recover_hpwl_under_overflow(...)` | 在 exact overflow cap 下进行 HPWL recovery；支持 axis-separated breakpoint、compact direction、node move、net-rigid block、density direction 与 contraction |
| compact recovery | `compact_support_contraction(...)` | 在 overflow guard 下尝试缩小布局支撑区域 |
| swap recovery | `recover_hpwl_with_equal_shape_swaps(...)` | equal-shape 精确交换、全局同形排列、anchor assignment；可用 net-aware neighborhood，等形交换保持 occupancy 不变 |

主要实现文件分别为 `src/bisection.cpp`、`src/coarse_flow.cpp`、`src/transport.cpp`、`src/density_coordinate.cpp`、`src/recovery.cpp`、`src/compact_recovery.cpp`、`src/swap_recovery.cpp`。

### 4.4 当前总控 `global_place()`

接口位于 `include/epsilon_active/placer.hpp`：

```cpp
PlaceResult global_place(Database& db, const PlaceConfig& config);
```

`PlaceConfig` 直接嵌套所有阶段 config 与主 GP 参数；`PlaceResult` 同样嵌套所有阶段 stats。`src/placer.cpp` 内部顺序固定为：

```text
recursive bisection
→ coarse capacity flow
→ transport
→ density coordinate descent
→ exact recovery
→ compact recovery
→ equal-shape/swap recovery
→ post recovery
→ post swap recovery
→ anchor assignment
→ 主 GP loop
→ 写 last.pl
→ 选择最低 HPWL feasible checkpoint；若无 feasible，则选择最低 overflow checkpoint
```

主 GP loop 每轮：

```text
exact overlap metrics + density search direction
→ 可选 net-share / regional price
→ exact HPWL + epsilon-active wire direction
→ lambda 初始化或更新
→ 组合并预条件梯度
→ optimizer.compute_delta()
→ 可选 exact batch backtracking，否则直接更新并 clamp 到区域
→ 写 metrics/checkpoint
```

注意：当前很多“后处理”函数实际上排在主 GP loop 之前。成果主链通常令 `--iterations 1`，在单独进程中只打开一个算子，使该进程近似成为一个阶段模块。这正是新架构应该显式化的行为。

### 4.5 CLI 与多级运行

`src/main.cpp` 负责：

```text
parse_options
→ read_bookshelf
→ center Gaussian 初始化或 load --initial-placement
→ 单级 global_place，或固定 64→128→256→512 的 multilevel 循环
→ write global.pl + summary.txt
```

多级模式只传递坐标；每级重建 fixed occupancy 与 exact overlap。regional price 可以 prolongate，但 lambda control 每级重新平衡。当前没有通用 pipeline 描述，H375 主链由 `scripts/run_h375_transfer.ps1` 手工调度。

## 5. 当前代表性成果调用链（以 adaptec1 / H375 为例）

下列链条是当前最适合用来理解“阶段模块化”需求的真实案例：

```text
raw Bookshelf
→ H219：HPWL-only Adam seed
→ H219：128×128 DCT/Poisson 同伦
→ H221：512×512 DCT/Poisson 细化并选择低-overflow mu-zero checkpoint
→ H252：15% cap 的 2 轮 exact recovery（含 net block）
→ H253：中等 density scale 的 81 步 Adam bridge
→ H254：强 density scale 的 120 步 retighten，重新进入 7% 可行区
→ H255：4 轮 node-only exact recovery
→ H257：1 轮预算化 node-only exact recovery（line search 4）
→ H372：surplus-only 容量递归二分，不做 axis pre-map
→ H372：10 轮 exact breakpoint/compact/net-block recovery
→ H375：5 轮 net-aware equal-shape swap polish
```

这里的 `raw Bookshelf` 指从原始 benchmark 建立问题；H219 随后会用 seed 219、`sigma_ratio=0.001` 的中心高斯坐标覆盖 movable 初始位置，并不是直接从原始 `.pl` 中的 movable 坐标做 HPWL descent。固定节点位置仍来自 benchmark。

当前复现脚本：`scripts/run_h375_transfer.ps1`。

这条链用于说明现有算法和阶段交接，不表示新项目必须把它原样作为第五题合规方案。尤其 H219/H221 的正 `mu` DCT/Poisson 阶段应被标为历史路线；新方案需要给出默认启用的纯非光滑 pipeline，或对辅助项的合规性作出明确论证。

### 5.1 分阶段代码映射

| 阶段 | 输入/输出 | 实际程序与核心调用 | 关键参数/行为 |
|---|---|---|---|
| H219 HPWL seed | raw/center Gaussian → `global.pl` | `epsilon_active.exe` → `read_bookshelf()` → `initialize_center_gaussian()` → `global_place()` | 512 bins，200 Adam iterations，`density_weight_scale=0`，seed 219 |
| H219 coarse | H219 seed → `snapshots/iter_1100.pl` | 归档 `homotopy.exe` / `homotopy_main.cpp` | 128 grid，24×50 steps，DCT/Poisson electrostatic + exact overlap + exact HPWL，adaptive/gradient-balanced weights |
| H221 fine | H219 iter 1100 → `best.pl` | 同一归档 homotopy | 512 grid，8×50 steps；最终阶段强制 `mu=0`；选择低-overflow/feasible HPWL checkpoint |
| H252 | H221 `best.pl` → `global.pl` | `global_place()` → `recover_hpwl_under_overflow()` | cap 15%，2 sweeps，breakpoint/compact node + density-aware net blocks |
| H253 | H252 → `global.pl` | `global_place()` 主 GP loop | 81 Adam steps，density scale 0.5，net-batch weight 1 |
| H254 | H253 → feasible `global.pl` | `global_place()` 主 GP loop | 120 Adam steps，density scale 4.75，较小 step，7% selector |
| H255/H257 | 前序布局 → `global.pl` | `recover_hpwl_under_overflow()` | 7% cap 下 node-only breakpoint + compact directions；H255 四轮，H257 一轮/line-search 4 |
| H372 capacity cut | H257 → `global.pl` | `recursive_hypergraph_bisection()` | position-seeded，surplus-only，leaf 8 bins，nearest-capacity + HPWL-guided；禁用 axis pre-map |
| H372 recovery | cut 输出 → `global.pl` | `recover_hpwl_under_overflow()` | 10 sweeps，line search 8，node + net block density direction + contraction，7% cap |
| H375 polish | H372 recovery → `global.pl` | `recover_hpwl_with_equal_shape_swaps()` | 5 sweeps，radius 64，32 candidates，exact shortlist 16，net-aware，只允许 equal shape |

### 5.2 这条链的代表性指标

以下只用于让规划者理解阶段作用，不要求新项目把数值硬编码：

| checkpoint | HPWL (M) | exact overflow |
|---|---:|---:|
| raw fresh | 53.2148 | 99.9662% |
| H219 HPWL seed | 43.0348 | 96.0563% |
| H219 coarse iter 1100 | 60.6175 | 71.7078% |
| H221 selected | 109.4255 | 7.7515% |
| H252 after 2 relaxed sweeps | 94.5485 | 15.0000% |
| H253 bridge | 86.2420 | 13.8912% |
| H254 retighten | 89.5493 | 6.5567% |
| H257 strict polish | 87.1195 | 7.0000% |
| H372 surplus-only cut | 93.7055 | 2.3064% |
| H372 recovery ×10 | 86.1896 | 7.0000% |
| H375 equal-shape ×5 | 85.9993 | 7.0000% |

数值来源为 `report/h375/data/h375_stage_summary.csv`。历史记录中的 runtime 口径存在“阶段增量”“从某个 checkpoint 开始”“完整 raw-to-final”和不同线程/机器估计混用的情况；规划时应建立统一的运行 provenance，而不要直接把某个汇总值当作跨机器可比的绝对结论。

## 6. 当前架构中最值得解决的耦合

### 6.1 静态问题、可变布局、oracle cache 混在一起

`Database` 把 netlist/geometry 与 node coordinates 放在同一对象中，各阶段原地修改；`ExactOverlapDensity` 又缓存 occupancy。阶段接口不能从类型上说明自己读取了什么、输出了什么，也不容易安全地并行比较两个候选分支。

### 6.2 模块边界存在，但调度边界不存在

算法函数已有 `Config + Stats` 雏形，但被 `global_place()` 固定串接。实验只能通过“大配置关掉其他模块 + iterations=1”或外部脚本多次启动程序来组合。

### 6.3 指标、接受准则和记录逻辑重复/分散

多个 recovery/transport/swap 实现自行调用 exact oracle，`placer.cpp` 还负责写不同 CSV。相同概念在 config、CLI、summary 和研究脚本中重复。新的顶层核心应提供统一 evaluator/auditor 与统一 run record，模块只返回轨迹事件或候选统计。

### 6.4 参数空间被历史实验开关淹没

当前 `PlaceConfig` 和 CLI 同时暴露主 GP、所有早期算子、所有消融开关以及多种 post 阶段。它适合积累实验，却不适合回答“某个模块的最小参数是什么”。

### 6.5 代码、参数与结果物理混放

实验目录同时容纳 protocol、analysis、巨大 `.pl`、CSV、PNG/PDF，归档可执行文件也在其中。Git 工作树很容易被研究产物淹没，且很难只上传算法和参数。

### 6.6 归档 homotopy 的依赖不可自描述

H219/H221 依赖相邻 `dreamplace-cpp` 的 DCT/electric 实现，但仓库内只有 `homotopy_main.cpp` 和未跟踪二进制。新项目必须选择：提取最小必要依赖、显式链接一个可定位依赖，或把该模块作为可选外部 adapter；不能让核心复现依赖某台机器上的相对目录。

## 7. 新项目的明确要求

这一节是对网页版规划者的硬要求，而不是可选建议。

### 7.1 统一的优化模块协议

每个阶段都必须能被理解为：

```text
静态 Problem（nodes/nets/pins/fixed geometry/region）
+ 输入 Layout（当前或前序模块的节点坐标）
+ 该模块自己的参数
+ RunContext（seed、threads、输出位置、日志/快照策略等）
→ 执行数值优化或离散搜索
→ 输出新的 Layout
+ 模块统计/轨迹
+ 由顶层权威 evaluator 复算的 before/after metrics
```

概念接口可类似：

```cpp
ModuleResult run(
    const Problem& problem,
    const Layout& input,
    const ModuleConfig& config,
    RunContext& context);
```

不要求照抄这个 C++ 签名，但所有模块必须遵守相同的输入/输出语义。默认不应让模块依赖神秘的全局状态或前一个模块留下的隐式缓存。若某个连续模块确实需要传递 Adam moments、lambda state 或 price field，应把它声明为显式可选 state artifact，并允许“只传坐标、重新初始化状态”。

### 7.2 统一且权威的顶层核心

顶层 core 至少负责：

- benchmark/problem 读取；
- layout 读取、写出、复制与坐标 clamp；
- exact weighted pin-offset HPWL；
- exact rectangle/bin overlap、density energy、overflow、max density；
- fixed macro/`terminal_NI` 的统一语义；
- before/after/final 的独立复算审计；
- checkpoint 与模块 provenance；
- 必要的共享增量审计原语，例如 node/group move 的 exact density delta。

优化模块可以使用非光滑方向启发式、价格场和 exact-audited 候选生成器。DCT/Poisson 或其他平滑 surrogate 只能存在于默认关闭的历史复现/对照模块中；如果第五题合规解释不允许任何平滑辅助，它们不得进入正式主 pipeline。无论模块类型如何，都必须清楚区分：

1. 用于生成方向/候选的内部量；
2. 用于接受、选择和报告的权威 exact metrics。

### 7.3 极简 pipeline 调度

应能用短小的 YAML/TOML/JSON 或少量代码排列模块，例如：

```yaml
pipeline:
  - module: hpwl_adam
    config: modules/hpwl_adam/params/seed.yaml
  - module: exact_joint_gp
    config: modules/exact_joint_gp/params/coarse_capacity.yaml
  - module: surplus_bisection
    config: modules/surplus_bisection/params/capacity_cut.yaml
  - module: exact_recovery
    config: modules/exact_recovery/params/cap07.yaml
  - module: equal_shape_swap
    config: modules/equal_shape_swap/params/net_aware_polish.yaml
```

这是默认的纯非光滑示意；H219/H221 的历史 DCT/Poisson 重放可以另设一份明确标记为 `historical` 的 pipeline。调度器只做顺序执行、layout/state 交接、指标审计、失败停止条件和记录；不要把算法逻辑重新集中到调度器里。模块应可重复、删除、交换顺序，并允许从任意已审计 checkpoint 启动。

### 7.4 代码、参数、结果分离

物理目录采用一个顶层架构文件夹和一个模块文件夹，不建立更多复杂层级：

```text
new-project/
├─ framework/
│  ├─ code/                        # problem/layout、Bookshelf I/O、exact evaluator、薄 runner
│  ├─ params/                      # dataset root、case list、pipeline 配置、全局默认值
│  └─ results/                     # 八 case 汇总和 pipeline 总览
└─ modules/
   └─ <module>/
      ├─ code/                     # 模块实现与最小说明
      ├─ params/                   # 模块参数 preset
      └─ results/                  # 模块中间布局、轨迹、日志和 checkpoint
```

具体目录名只能做小幅调整，但每个模块的 `code/params/results` 三分必须一眼可见。只上传 `framework/code`、`framework/params`、各模块 `code/params` 和必要说明即可复现实验设计；`results/` 目录可以保留空结构，但大体积内容默认被 `.gitignore` 排除。结果中仍保留完整 config snapshot、版本/commit、输入 checkpoint hash、seed、threads、机器信息、开始/结束指标和 runtime。

### 7.5 首批应模块化的主链单元

新项目第一阶段不必迁移全部历史启发式。优先覆盖能复现代表性调用链的模块：

1. layout 初始化/加载；
2. HPWL-only nonsmooth Adam/subgradient；
3. exact-overlap nonsmooth joint GP / retighten；
4. exact recovery（node、breakpoint/compact、可选 net block）；
5. surplus-only recursive bisection；
6. net-aware equal-shape swap；
7. evaluator/auditor 与 pipeline runner；
8. 一条不依赖光滑 wirelength 或 smooth density warm start 的默认纯非光滑 pipeline。

DCT/Poisson electrostatic homotopy 只在需要复现 H219/H221 历史结果时作为可选模块迁移，不应阻塞上述轻量第一版，也不应成为默认合规路线。

`coarse_flow`、旧 transport、density-coordinate cluster、assignment、各种失败消融可以后续按价值迁移，不应为了“接口完整”一次性搬入新仓库。

### 7.6 必须保留的实验可复现信息

- 输入 benchmark 与输入 layout 的可验证标识；
- 完整模块 config，而不是只保存命令行；
- seed、有效线程数、数值精度/编译信息；
- 每个阶段独立 before/after exact metrics；
- module runtime 与 pipeline cumulative runtime 分开；
- selected checkpoint 与 trajectory-last checkpoint 分开；
- 若使用 surrogate/homotopy，最终必须 exact audit，并记录 surrogate state 是否被清零；
- 支持按模块保存少量轨迹 CSV 和可配置 snapshot，不默认产生大量文件。

### 7.7 八个 case 的统一验收入口

顶层 `framework` 应有一个很薄的 batch 入口，直接枚举指定数据根目录下的：

```text
adaptec1 adaptec2 adaptec3 adaptec4
bigblue1 bigblue2 bigblue3 bigblue4
```

它对每个 case 使用同一套指标口径，分别记录成功/失败、HPWL、overflow、CPU runtime、有效并行数和迭代次数，再生成八 case 汇总。汇总必须能直接判断第五题的六项技术诉求，其中第（6）项在题意澄清前以显式配置项表示。不要引入任务队列或通用实验管理平台，一个普通循环和 CSV/Markdown 汇总就足够。

## 8. 这不是一份传统软件工程项目

规划和实现时应显式采用以下优先级：

```text
数值正确性与指标口径一致
> 快速组合/替换优化阶段
> 实验可复现与 provenance
> 性能和内存成本可控
> 代码短小、容易修改
> 通用软件工程完备性
```

暂时不需要：

- 微服务、RPC、数据库、Web 后端；
- 大型插件框架、依赖注入框架或复杂 runtime registry；
- 为未来未知算法设计过度通用的 DSL；
- 企业级权限、服务治理、复杂兼容层；
- 把所有历史实验和消融都迁移；
- 为了“纯架构”拆出大量只有一两个函数的层；
- 过早追求稳定公共 API、ABI 或跨语言 SDK。

需要的是研究友好的最小结构：少量清楚的数据类型、一个权威 evaluator、一个很薄的 runner、若干可直接阅读和替换的模块、文本参数文件，以及各模块中结构明确但内容默认被 Git 忽略的 `results/` 目录。

测试也应以数值契约为中心：HPWL pin offset、固定宏容量、`terminal_NI`、exact move/group delta、一致的 before/after 审计、同形交换的 occupancy 不变性、checkpoint round-trip，以及少量小规模 pipeline smoke test。无需先建立庞大的软件工程测试矩阵。

## 9. 对重构方案输出的期待

请基于本文规划一个**新仓库**，不要默认在当前仓库上做渐进式大手术。方案应回答：

1. 如何严格落实 `framework/` 与 `modules/<module>/{code,params,results}/` 的最小目录结构和依赖关系；
2. `Problem`、`Layout`、`Metrics`、`ModuleResult`、显式可选 solver state 应如何划分；
3. evaluator 与模块如何共享高性能 exact oracle/增量审计，而不产生隐藏状态错误；
4. 模块的注册/选择机制如何保持极简；
5. pipeline 配置和模块 config 如何组织、校验与覆盖；
6. H219→H375 主链如何一对一映射到模块；
7. 面向第五题的纯非光滑默认 pipeline 是什么，为什么没有用光滑近似替代原始目标；
8. 历史 DCT/Poisson 依赖是否迁移、如何隔离，以及它能否进入合规路线；
9. 如何默认从 `D:\codex_project\HUAWEI_EDA\alg-electronic\ispd2005` 读取全部八个 case；
10. run directory、checkpoint manifest、metrics CSV 和 config snapshot 的最小格式是什么；
11. 从当前代码迁移的最小顺序是什么，每一步如何用数值回归确认没有改变评价口径；
12. 哪些历史模块暂不迁移，以及为什么。

方案应给出足够明确、可以马上照着创建代码的文件级实施步骤、接口草案、配置示例和最小运行命令，而不是只输出抽象架构图。避免把项目设计成通用工业框架。最终目标是：研究者能在很短时间内新增一个“读入布局 → 非光滑优化 → 输出布局”的模块，用一行 pipeline 配置插入、删除或重排它，并能可信地比较每个阶段以及八个 case 的 exact HPWL、overflow、runtime 与迭代预算。

## 10. 建议随本文一起提供给规划模型的核心文件

若上传大小允许，优先顺序如下：

1. `include/epsilon_active/types.hpp`
2. `include/epsilon_active/placer.hpp`
3. `src/placer.cpp`
4. `src/main.cpp`
5. `include/epsilon_active/hpwl.hpp` 与 `src/hpwl.cpp`
6. `include/epsilon_active/density.hpp` 与 `src/density.cpp`
7. `include/epsilon_active/recovery.hpp` 与 `src/recovery.cpp`
8. `include/epsilon_active/bisection.hpp` 与 `src/bisection.cpp`
9. `include/epsilon_active/swap_recovery.hpp` 与 `src/swap_recovery.cpp`
10. `include/epsilon_active/optimizer.hpp`、`src/optimizer.cpp`
11. `include/epsilon_active/lambda_controller.hpp`、`src/lambda_controller.cpp`
12. `experiments/h72_homotopy_electrostatic/homotopy_main.cpp`
13. `scripts/run_h375_transfer.ps1`
14. `report/h375/data/h375_stage_summary.csv`
15. `tests/core_tests.cpp`

若上传大小受限，至少提供 1–9、12、13；`bisection.cpp`、`recovery.cpp` 较长时，可先提供对应 header 与本文，待规划模型需要具体实现再补充源码。

## 11. 当前事实来源

- 非光滑题目与八 case 技术诉求：`D:\codex_project\HUAWEI_EDA\第148期 - EDA专题第四期.pdf` 第五题
- 指定 ISPD2005 数据：`D:\codex_project\HUAWEI_EDA\alg-electronic\ispd2005`
- 数学与模块边界：`README.md`、`docs/architecture.md`、`项目介绍.md`
- 当前研究状态：`research-state.yaml`、`findings.md`、`research-log.md`
- 当前 public API：`include/epsilon_active/*.hpp`
- 当前实现与调用顺序：`src/*.cpp`
- H219–H375 复现链：`scripts/run_h375_transfer.ps1`
- H252–H257、H372、H375 阶段意图：相应 `experiments/<experiment>/protocol.md` 与 `analysis.md`
- 阶段指标：`report/h375/data/h375_stage_summary.csv`
- DCT/Poisson homotopy：`experiments/h72_homotopy_electrostatic/homotopy_main.cpp`
