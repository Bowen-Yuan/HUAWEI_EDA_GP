# epsilon-active 新项目代码实施方案 V2（供本地 AI 直接执行）

> 目标：基于现有 `epsilon-active` 的数学契约与已验证阶段边界，新建一个**轻量、模块化、面向非光滑 Global Placement 数值实验**的 C++17 项目。
>
> 本文件不是抽象架构说明，而是可直接用于创建仓库、迁移代码、编写接口、建立测试和运行数值实验的实施蓝图。日常开发默认只运行 `adaptec1`；只有在用户明确要求完整实验时才运行全部八个 ISPD2005 case。
>
> 设计原则：**数值契约优先于架构优雅；exact evaluator 是唯一权威；默认主链必须是纯非光滑路线；模块可顺序组合；结果可独立复算；不要建设通用工业框架。**

> **V2 实验工作流约定**：日常开发只跑 `adaptec1`；不做 baseline 比较；pipeline 仅是模块组合示例而非最终方案；stage 间必须可靠交接；默认只保留结果、轻量 trajectory 和可继续使用的 audited `selected.pl`，其余中间内存/磁盘产物及时释放或删除。

---

## 1. 本地 AI 执行时必须遵守的最高优先级约束

### 1.1 原始问题不能被改写成平滑 GP

新项目正式路线直接处理：

```text
min W(x, y) + lambda * D(x, y)

W = exact weighted pin-offset HPWL (max-min)
D = exact rectangle/bin overlap 对应的分段非光滑 density penalty
```

必须满足：

1. 禁止用 weighted-average wirelength、log-sum-exp、Moreau envelope 等替代 exact HPWL。
2. 禁止用 smooth density surrogate 替代最终 exact rectangle/bin overlap 指标。
3. 允许 epsilon-active subgradient、普通 subgradient、bundle/proximal、exact breakpoint/coordinate search、exact-audited discrete move、capacity assignment 等非光滑方法。
4. 模块内部用于“生成方向/候选”的量，与“接受/选择/报告”的 exact metrics 必须分开。
5. 每个模块结束后，由 `framework` 从输出 `Layout` 重新构建 evaluator 并独立复算 HPWL / density / overflow；不得直接相信模块内部累计值。
6. H219/H221 的 DCT/Poisson homotopy 只能属于 `historical_dct_poisson` 可选模块；默认 challenge pipeline 禁止引用它。

### 1.2 项目不是企业软件

禁止首版引入：

- 微服务、RPC、数据库、Web 后端；
- 动态插件框架、依赖注入；
- 复杂 registry / workflow engine；
- 为未来未知算法设计 DSL；
- 大量 manager/service/adapter 空层；
- 为 ABI、SDK、跨语言接口做提前设计。

允许且推荐：

- C++17；
- OpenMP；
- CMake；
- 一个很小的 JSON 配置依赖；
- 静态函数表完成模块选择；
- 普通 `struct + vector + function`。

### 1.3 数据集位置

默认数据根目录必须是：

```text
D:\codex_project\HUAWEI_EDA\alg-electronic\ispd2005
```

数据目录包含八个 case：

```text
adaptec1 adaptec2 adaptec3 adaptec4
bigblue1 bigblue2 bigblue3 bigblue4
```

**日常运行默认 `adaptec1`。** 完整八 case 只在用户明确要求完整实验时触发。

每个 case 的 benchmark base：

```text
<dataset_root>/<case>/<case>
```

需要读取：

```text
.aux .nodes .nets .pl .scl
.wts   # 若存在
```

数据禁止复制进新仓库。

### 1.4 评价语义必须逐字保真

HPWL：

```text
pin_x = node_center_x + pin_offset_x
pin_y = node_center_y + pin_offset_y

HPWL_e = max(pin_x) - min(pin_x)
       + max(pin_y) - min(pin_y)

HPWL = sum_e net_weight_e * HPWL_e
```

Density：

```text
occupancy_b = movable rectangle 与 bin 的精确相交面积
            + physical fixed macro 与 bin 的精确相交面积

rho_b = occupancy_b / bin_area

density_energy = sum_b 0.5 * bin_area * max(rho_b - target_density, 0)^2

overflow = sum_b bin_area * max(rho_b - target_density, 0)
         / total_movable_area
```

必须保留：

- fixed macro 消耗容量且不能移动；
- `terminal_NI` 不消耗物理容量；
- 内存坐标为 node center；Bookshelf `.pl` 为左下角；
- epsilon-active 只改变方向选择，不改变 exact HPWL 值；
- overflow 达标不代表 legalized；
- 本项目不做 row/site snapping、overlap removal、legalization、detailed placement。

---

## 2. V1 项目范围

V1 的目标不是“迁移旧仓库所有功能”，而是建立一条可信、可修改、可复现实验的主干。

### 2.1 V1 必须有

1. Bookshelf I/O。
2. `Problem` / `Layout` 分离。
3. exact HPWL evaluator。
4. exact density/overflow evaluator。
5. node/group move 的 exact density incremental audit。
6. 顶层 fresh audit。
7. 极简 module protocol。
8. JSON pipeline runner。
9. `layout_init`。
10. `hpwl_adam`：epsilon-active HPWL subgradient + Adam/SGD 等。
11. `exact_joint_gp`：exact HPWL direction + exact overlap direction + lambda control + optimizer。
12. `exact_recovery`。
13. `surplus_bisection`。
14. `equal_shape_swap`。
15. 单 case（默认 `adaptec1`）实验入口；保留可选八 case batch 入口，但不作为日常开发默认路径。
16. 统一 metrics / manifest / config snapshot / runtime。
17. 数值契约测试。
18. 一个“参考 pipeline”配置，用来证明模块可顺序交接、重复、删除和重排；它不是最终算法方案。
19. 实验产物生命周期管理：只长期保留结果、关键轨迹数据和可继续衔接的 audited `.pl`，其余临时文件/缓存默认清理。

### 2.2 V1 暂不迁移

除非主链性能验证表明确显示必要，否则暂缓：

- coarse flow 全部变体；
- 旧 transport 全部变体；
- density-coordinate cluster 大量启发式；
- assignment 的历史分支；
- 大量失败消融开关；
- H219/H221 DCT/Poisson（仅在历史复现确有需求时单独迁移）。

原因：这些功能会显著扩大参数面和维护面积，但不是建立新的非光滑实验主干所必需。

---

## 3. 推荐仓库名称与物理目录

示例仓库名：

```text
nonsmooth-gp/
```

必须维持“一个 framework + modules”的两层逻辑：

```text
nonsmooth-gp/
├─ CMakeLists.txt
├─ README.md
├─ .gitignore
│
├─ framework/
│  ├─ code/
│  │  ├─ common.hpp
│  │  ├─ problem.hpp
│  │  ├─ problem.cpp
│  │  ├─ layout.hpp
│  │  ├─ layout.cpp
│  │  ├─ bookshelf.hpp
│  │  ├─ bookshelf.cpp
│  │  ├─ metrics.hpp
│  │  ├─ exact_hpwl.hpp
│  │  ├─ exact_hpwl.cpp
│  │  ├─ exact_density.hpp
│  │  ├─ exact_density.cpp
│  │  ├─ incremental_density.hpp
│  │  ├─ incremental_density.cpp
│  │  ├─ auditor.hpp
│  │  ├─ auditor.cpp
│  │  ├─ module_api.hpp
│  │  ├─ module_registry.hpp
│  │  ├─ module_registry.cpp
│  │  ├─ config.hpp
│  │  ├─ config.cpp
│  │  ├─ run_record.hpp
│  │  ├─ run_record.cpp
│  │  ├─ pipeline.hpp
│  │  ├─ pipeline.cpp
│  │  ├─ batch.hpp
│  │  ├─ batch.cpp
│  │  └─ main.cpp
│  │
│  ├─ params/
│  │  ├─ defaults.json
│  │  ├─ cases.json
│  │  └─ pipelines/
│  │     ├─ reference_nonsmooth_chain.json
│  │     ├─ smoke.json
│  │     └─ historical_h219_h375.json   # 仅历史参考
│  │
│  └─ results/
│     └─ .gitkeep
│
├─ modules/
│  ├─ layout_init/
│  │  ├─ code/
│  │  │  ├─ layout_init.hpp
│  │  │  └─ layout_init.cpp
│  │  ├─ params/
│  │  │  ├─ raw.json
│  │  │  └─ center_gaussian_219.json
│  │  └─ results/.gitkeep
│  │
│  ├─ hpwl_adam/
│  │  ├─ code/
│  │  │  ├─ optimizer.hpp
│  │  │  ├─ optimizer.cpp
│  │  │  ├─ hpwl_adam.hpp
│  │  │  └─ hpwl_adam.cpp
│  │  ├─ params/
│  │  │  ├─ seed_h219_like.json
│  │  │  └─ smoke.json
│  │  └─ results/.gitkeep
│  │
│  ├─ exact_joint_gp/
│  │  ├─ code/
│  │  │  ├─ lambda_controller.hpp
│  │  │  ├─ lambda_controller.cpp
│  │  │  ├─ batch_acceptance.hpp
│  │  │  ├─ batch_acceptance.cpp
│  │  │  ├─ exact_joint_gp.hpp
│  │  │  └─ exact_joint_gp.cpp
│  │  ├─ params/
│  │  │  ├─ bridge_h253_like.json
│  │  │  ├─ retighten_h254_like.json
│  │  │  ├─ pure_coarse.json
│  │  │  └─ smoke.json
│  │  └─ results/.gitkeep
│  │
│  ├─ exact_recovery/
│  │  ├─ code/
│  │  │  ├─ exact_recovery.hpp
│  │  │  └─ exact_recovery.cpp
│  │  ├─ params/
│  │  │  ├─ cap15_netblock.json
│  │  │  ├─ cap07_node.json
│  │  │  └─ cap07_full.json
│  │  └─ results/.gitkeep
│  │
│  ├─ surplus_bisection/
│  │  ├─ code/
│  │  │  ├─ surplus_bisection.hpp
│  │  │  └─ surplus_bisection.cpp
│  │  ├─ params/
│  │  │  └─ h372_like.json
│  │  └─ results/.gitkeep
│  │
│  ├─ equal_shape_swap/
│  │  ├─ code/
│  │  │  ├─ equal_shape_swap.hpp
│  │  │  └─ equal_shape_swap.cpp
│  │  ├─ params/
│  │  │  └─ h375_like.json
│  │  └─ results/.gitkeep
│  │
│  └─ historical_dct_poisson/           # 可选，V1 可先只有 README
│     ├─ code/
│     │  └─ README.md
│     ├─ params/
│     │  └─ disabled.json
│     └─ results/.gitkeep
│
└─ tests/
   ├─ test_bookshelf.cpp
   ├─ test_hpwl.cpp
   ├─ test_density.cpp
   ├─ test_incremental_density.cpp
   ├─ test_checkpoint_roundtrip.cpp
   ├─ test_equal_shape_swap.cpp
   └─ test_pipeline_smoke.cpp
```

### 3.1 依赖方向必须单向

```text
framework/code/common, problem, layout
        ↓
framework/code/bookshelf, exact_hpwl, exact_density, incremental_density
        ↓
framework/code/auditor, module_api, run_record
        ↓
modules/*
        ↓
framework/code/module_registry, pipeline, batch, main
```

关键规则：

- `framework` 的 evaluator 不依赖任何 module。
- module 可以依赖 `framework`。
- module 之间不直接 include 对方实现。
- pipeline 只通过统一 `ModuleRunner` 调用 module。
- 历史 DCT/Poisson 不得成为 core link dependency。

---

## 4. 编译与第三方依赖

### 4.1 编译标准

```text
C++17
CMake >= 3.20
OpenMP
```

建议只引入一个配置依赖：`nlohmann/json`。

理由：

- header-only；
- 配置结构明确；
- 避免自己实现 parser；
- 比首版引入完整 YAML 依赖更简单。

如果目标机器已经通过包管理器安装：

```cmake
find_package(nlohmann_json CONFIG REQUIRED)
```

若机器离线，可在仓库 `third_party/` 放单个 `json.hpp`；不要通过运行时联网下载依赖。

### 4.2 根 `CMakeLists.txt` 目标

推荐最终只生成一个主程序：

```text
nsgp.exe
```

子命令：

```text
nsgp run
nsgp batch
nsgp audit
nsgp list-modules
```

另生成一个测试目标：

```text
nsgp_tests.exe
```

不要为每个 module 生成单独 executable；它们作为静态库或直接 sources 链入主程序即可。

---

## 5. 核心数据结构设计

核心目标是彻底拆开“静态问题”和“可变布局”。

### 5.1 `common.hpp`

```cpp
#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace nsgp {

using Real = double;
using NodeId = std::int32_t;
using NetId  = std::int32_t;
using PinId  = std::int32_t;

struct Rect {
    Real xl = 0;
    Real yl = 0;
    Real xh = 0;
    Real yh = 0;
};

} // namespace nsgp
```

首版保持 `double`，不要同时支持 float/double 模板化，避免扩大验证面。

### 5.2 `Problem`

`Problem` 只保存静态信息，一旦读入 benchmark 后不再改变。

```cpp
struct NodeInfo {
    std::string name;
    Real width = 0;
    Real height = 0;
    bool movable = false;
    bool fixed = false;
    bool terminal_ni = false;
};

struct PinInfo {
    NodeId node = -1;
    NetId net = -1;
    Real offset_x = 0;
    Real offset_y = 0;
};

struct NetInfo {
    std::string name;
    Real weight = 1.0;
    std::vector<PinId> pins;
};

struct RowInfo {
    // 只保留当前 exact global placement 所需 Bookshelf row/site 信息。
};

struct Problem {
    std::string case_name;
    std::filesystem::path benchmark_base;

    std::vector<NodeInfo> nodes;
    std::vector<PinInfo> pins;
    std::vector<NetInfo> nets;
    std::vector<RowInfo> rows;

    std::vector<NodeId> movable_ids;
    std::vector<NodeId> fixed_ids;

    Rect region;
    Real movable_area = 0;
};
```

要求：

- `Problem` 不存 node 当前 `x/y`。
- fixed node 的位置也放进 `Layout`，因为 Bookshelf placement 是一个完整布局；但算法禁止修改 fixed 坐标。
- 可附加只读索引（如 node-to-pin），但它们必须由 benchmark 静态构建。

### 5.3 `Layout`

```cpp
struct Layout {
    std::vector<Real> x;   // center x
    std::vector<Real> y;   // center y
    std::vector<std::uint8_t> orientation;

    std::uint64_t revision = 0;

    void touch() { ++revision; }
};
```

约束：

1. `x.size() == y.size() == problem.nodes.size()`。
2. 坐标均为 center。
3. fixed 节点不得被 module 改写。
4. 每个批量坐标更新完成后必须 `layout.touch()`。
5. incremental cache 记录建立时的 `revision`，发现 revision 不匹配立即报错。

不要把 metrics、optimizer moments、lambda 值塞进 `Layout`。

### 5.4 `Metrics`

```cpp
struct Metrics {
    Real hpwl = 0;
    Real density_energy = 0;
    Real overflow = 0;
    Real max_density = 0;
};
```

必要时可增加：

```cpp
std::uint64_t evaluated_layout_revision;
```

但不要把非权威内部 loss 混进此结构。

### 5.5 Direction / internal objective 必须使用不同类型

例如：

```cpp
struct SearchDirection {
    std::vector<Real> dx;
    std::vector<Real> dy;
};

struct InternalObjectiveTrace {
    Real hpwl_direction_norm = 0;
    Real density_direction_norm = 0;
    Real lambda = 0;
};
```

任何 `smooth_loss`、`price_loss`、`surrogate_energy` 均禁止写入 `Metrics`。

---

## 6. Bookshelf I/O 设计

### 6.1 接口

```cpp
struct LoadedBenchmark {
    Problem problem;
    Layout raw_layout;
};

LoadedBenchmark read_bookshelf(
    const std::filesystem::path& benchmark_base);

Layout read_placement(
    const Problem& problem,
    const std::filesystem::path& pl_path);

void write_placement(
    const Problem& problem,
    const Layout& layout,
    const std::filesystem::path& pl_path);

Layout initialize_center_gaussian(
    const Problem& problem,
    const Layout& base,
    std::uint64_t seed,
    Real sigma_ratio);
```

### 6.2 坐标转换

读取 `.pl`：

```text
center_x = bookshelf_xl + width / 2
center_y = bookshelf_yl + height / 2
```

写 `.pl`：

```text
bookshelf_xl = center_x - width / 2
bookshelf_yl = center_y - height / 2
```

### 6.3 I/O 测试

必须有：

1. raw `.pl` 读取 → 写出 → 再读取；
2. movable/fixed 坐标在 tolerance 内一致；
3. orientation 保留；
4. pin offset 不受 placement I/O 影响。

---

## 7. 权威 exact evaluator

## 7.1 Exact HPWL

建议接口：

```cpp
struct HpwlDirectionConfig {
    Real epsilon = 0;
    Real active_power = 1;
    int degree_limit = -1;
};

class ExactHpwl {
public:
    explicit ExactHpwl(const Problem& problem);

    Real evaluate(const Layout& layout) const;

    Real evaluate_with_direction(
        const Layout& layout,
        const HpwlDirectionConfig& cfg,
        SearchDirection& direction) const;
};
```

必须保证：

- 两个接口返回的 HPWL 完全同一数学定义；
- `epsilon/active_power/degree_limit` 只能影响 direction；
- `degree_limit` 不得改变 exact HPWL 的计算 net 集合；
- pin offset 和 net weight 必须包含。

## 7.2 Exact density

```cpp
struct DensityGridConfig {
    int bins_x = 0;
    int bins_y = 0;
    Real target_density = 0;
};

struct DensityDirectionConfig {
    // epsilon-active overlap piece selection 等，只影响方向。
};

class ExactDensity {
public:
    ExactDensity(const Problem& problem, DensityGridConfig cfg);

    Metrics evaluate_density_only(const Layout& layout) const;

    void evaluate_with_direction(
        const Layout& layout,
        const DensityDirectionConfig& cfg,
        Metrics& metrics,
        SearchDirection& direction) const;
};
```

实现约束：

- fixed occupancy 在构造时可预计算；
- 每次 fresh audit 都从 fixed occupancy 出发重新累计 movable occupancy；
- `terminal_NI` 不进入 physical occupancy；
- `max_density` 与 overflow 来自同一 occupancy。

## 7.3 顶层 `ExactAuditor`

```cpp
class ExactAuditor {
public:
    Metrics evaluate(
        const Problem& problem,
        const Layout& layout,
        const DensityGridConfig& density_cfg) const;
};
```

`ExactAuditor` 是 pipeline 的唯一正式 metric 入口。

模块内部可以调用同一底层 exact 实现，但 pipeline 在模块前后必须重新调用一次 fresh audit。

### 7.4 Fresh audit 的强制时机

必须：

```text
benchmark load 后
pipeline 输入 checkpoint 后
每个 module 执行前
每个 module selected output 后
每个 module trajectory-last output 后（若不同）
pipeline final 后
batch summary 写出前
```

---

## 8. Incremental exact density 审计

recovery、swap、transport 类型 module 需要高性能增量审计，但不能让 cache 穿越 module 边界。

### 8.1 作用域规则

`IncrementalDensity` 必须是 module-local object：

```cpp
class IncrementalDensity {
public:
    IncrementalDensity(
        const Problem& problem,
        const Layout& layout,
        DensityGridConfig cfg);

    DensityMove evaluate_move(
        const Layout& layout,
        NodeId node,
        Real new_x,
        Real new_y) const;

    DensityGroupMove evaluate_group_move(
        const Layout& layout,
        const std::vector<NodeMove>& moves) const;

    void commit_move(Layout& layout, const DensityMove& move);
    void commit_group_move(Layout& layout, const DensityGroupMove& move);

    void assert_revision(const Layout& layout) const;
};
```

规则：

- constructor 记录 `layout.revision`；
- commit 同时更新 occupancy cache 和 layout，并更新 expected revision；
- 外部直接改 layout 后，旧 incremental object 必须失效；
- module 结束后销毁；
- pipeline 不传递它。

### 8.2 必须测试

随机小问题上：

```text
incremental move delta
vs
fresh full recomputation after applying move
```

对以下量逐一比：

- density energy；
- overflow；
- max density（允许因实现方式做 fresh compare）；
- group move。

---

## 9. Module 统一协议

### 9.1 通用上下文

```cpp
struct RunContext {
    std::string run_id;
    std::string case_name;
    int stage_index = 0;
    int threads = 1;
    std::uint64_t seed = 0;

    std::filesystem::path stage_result_dir;
    bool save_trajectory = false;
    int snapshot_every = 0;
};
```

### 9.2 显式可选 solver state artifact

不要让 `framework` 理解 Adam moments 的字段。

首版最简单的做法：state 以 module 私有序列化文件存在，framework 只传路径和元数据。

```cpp
struct StateArtifact {
    std::string kind;                 // e.g. "hpwl_adam_state_v1"
    std::filesystem::path path;
    std::string file_sha256;
};
```

如果 module 不需要 continuation，返回 `nullopt`。

pipeline stage 可指定：

```json
"state_input": null,
"state_policy": "reset"
```

允许值：

```text
reset       只传 Layout，重新初始化 optimizer/lambda/price
resume      要求读入显式 state artifact
```

V1 默认所有正式 pipeline stage 使用 `reset`，除非某个实验明确研究 continuation。

### 9.3 `ModuleResult`

```cpp
struct ModuleStats {
    int iterations = 0;
    int accepted_moves = 0;
    int rejected_moves = 0;
    double wall_seconds = 0;
};

struct ModuleResult {
    Layout selected_layout;
    Layout last_layout;
    bool selected_is_last = true;

    ModuleStats stats;
    std::optional<StateArtifact> state_out;

    std::string status;   // "ok" / "failed"
    std::string message;
};
```

注意：

- `ModuleResult` 不负责填写权威 before/after Metrics。
- pipeline 接到 result 后再 fresh audit。
- `selected_layout` 与 `last_layout` 必须显式分开，彻底消除旧 `global.pl`/`last.pl` 名称歧义。

### 9.4 Module 函数签名

使用 JSON 作为 module config 容器：

```cpp
using Json = nlohmann::json;

using ModuleRunner = ModuleResult (*)(
    const Problem& problem,
    const Layout& input,
    const Json& config,
    RunContext& context);
```

为什么不建立复杂 `ModuleConfig` 基类：

- 每个 module 参数完全不同；
- module 自己 `parse_and_validate_config()` 即可；
- framework 不需要理解算法参数。

---

## 10. 极简 module registry

不要动态加载 DLL。

```cpp
struct ModuleEntry {
    const char* name;
    ModuleRunner runner;
    const char* compliance_tag;
};

const std::unordered_map<std::string, ModuleEntry>& module_registry();
```

`module_registry.cpp` 直接静态列出：

```cpp
{
  {"layout_init",        {"layout_init",        &run_layout_init,        "exact_nonsmooth"}},
  {"hpwl_adam",         {"hpwl_adam",         &run_hpwl_adam,         "exact_nonsmooth"}},
  {"exact_joint_gp",    {"exact_joint_gp",    &run_exact_joint_gp,    "exact_nonsmooth"}},
  {"exact_recovery",    {"exact_recovery",    &run_exact_recovery,    "exact_nonsmooth"}},
  {"surplus_bisection", {"surplus_bisection", &run_surplus_bisection, "exact_nonsmooth"}},
  {"equal_shape_swap",  {"equal_shape_swap",  &run_equal_shape_swap,  "exact_nonsmooth"}},
  {"historical_dct_poisson", {"historical_dct_poisson", &run_historical_dct_poisson, "historical_surrogate"}}
};
```

增加新 module 的步骤必须只有：

1. 建 `modules/<name>/{code,params,results}`；
2. 写 `run_<name>()`；
3. 在 registry 加一行；
4. pipeline 加一行。

这就是“插件机制”的全部，不要继续抽象。

---

## 11. Challenge compliance gate

为防止历史 surrogate 混入正式路线，pipeline 顶层增加：

```json
"compliance_mode": "challenge_nonsmooth"
```

runner 规则：

```text
challenge_nonsmooth:
    只允许 compliance_tag == exact_nonsmooth

research:
    所有已注册 module 可运行
```

若 `historical_dct_poisson` 出现在 challenge pipeline：

```text
立即拒绝启动
打印明确错误
不允许仅靠改 module 名规避
```

历史 module 的 manifest 必须另外写：

```json
"uses_smooth_surrogate": true,
"surrogate_kind": "DCT/Poisson electrostatic",
"surrogate_zeroed_at_final_stage": true
```

即使最后 `mu=0`，也不得自动标记为 challenge-compliant。

---

## 12. Pipeline 配置格式

V1 用 JSON，不同时支持 YAML/TOML，避免多套解析逻辑。

### 12.1 顶层默认配置 `framework/params/defaults.json`

```json
{
  "dataset_root": "D:\\codex_project\\HUAWEI_EDA\\alg-electronic\\ispd2005",
  "threads": 1,
  "seed": 219,
  "density": {
    "target_density": 0.9,
    "bins_x": 512,
    "bins_y": 512
  },
  "output_root": "framework/results",
  "save_trajectory": true,
  "trajectory_every": 1,
  "snapshot_every": 0,
  "artifacts": {
    "policy": "minimal",
    "keep_selected_pl": true,
    "keep_last_pl": false,
    "keep_solver_state": false,
    "cleanup_on_stage_success": true
  }
}
```

注意：`target_density=0.9` 这里仅作为**配置示例**，实现时应从旧项目正式配置/题目要求确认后再决定默认值；不要因为示例存在就把它视为已从 brief 证明的挑战参数。

更稳妥的首版实现方式是：若未显式提供 `target_density`，直接报错，不自作默认。

### 12.2 默认 case 与完整 case 清单

日常数值开发默认只跑 `adaptec1`。建议：

`framework/params/cases.json`：

```json
{
  "default_case": "adaptec1",
  "all_cases": [
    "adaptec1", "adaptec2", "adaptec3", "adaptec4",
    "bigblue1", "bigblue2", "bigblue3", "bigblue4"
  ]
}
```

CLI 的 `run` 若没有给 `--case`，直接使用 `adaptec1`。不要因为完整列表存在就自动跑八个 case。只有显式执行 `batch --all-cases` 或等价命令时才跑完整实验。

### 12.3 pipeline schema

下面只是**能力参考配置**，用于验证模块顺序交接，不代表最终算法链。runner 绝不能对其中任何顺序做特殊处理。

```json
{
  "name": "reference_nonsmooth_chain",
  "compliance_mode": "challenge_nonsmooth",
  "stages": [
    {
      "module": "layout_init",
      "config": "modules/layout_init/params/center_gaussian_219.json",
      "state_policy": "reset"
    },
    {
      "module": "hpwl_adam",
      "config": "modules/hpwl_adam/params/seed_h219_like.json",
      "state_policy": "reset"
    },
    {
      "module": "exact_joint_gp",
      "config": "modules/exact_joint_gp/params/pure_coarse.json",
      "state_policy": "reset"
    },
    {
      "module": "exact_joint_gp",
      "config": "modules/exact_joint_gp/params/retighten_h254_like.json",
      "state_policy": "reset"
    },
    {
      "module": "surplus_bisection",
      "config": "modules/surplus_bisection/params/h372_like.json",
      "state_policy": "reset"
    },
    {
      "module": "exact_recovery",
      "config": "modules/exact_recovery/params/cap07_full.json",
      "state_policy": "reset"
    },
    {
      "module": "equal_shape_swap",
      "config": "modules/equal_shape_swap/params/h375_like.json",
      "state_policy": "reset"
    }
  ]
}
```

### 12.4 参数覆盖顺序

固定为：

```text
framework defaults
< pipeline-level overrides
< module preset file
< per-stage overrides
< CLI overrides
```

例如：

```json
{
  "module": "exact_recovery",
  "config": "modules/exact_recovery/params/cap07_full.json",
  "overrides": {
    "sweeps": 12
  }
}
```

### 12.5 配置校验规则

- 未知 module：报错；
- config 文件不存在：报错；
- module config 未知 key：默认报错；
- `threads < 1 || threads > 40`：报错；
- bins 非正：报错；
- target density 非法：报错；
- overflow cap 非法：报错；
- challenge pipeline 含 surrogate tag：报错；
- resume 但没有 state artifact：报错；
- pipeline stage 数量可以是 1 到 N；同一 module 可出现多次；
- runner 不校验“算法顺序是否符合某个参考链”，只校验输入/配置/state 契约。

研究代码中“静默忽略未知参数”非常危险，禁止。

---

## 13. 各 module 详细职责

## 13.1 `layout_init`

用途：把 raw benchmark placement 转为实验起点。

模式：

```text
raw
center_gaussian
load_external
```

配置示例：

```json
{
  "mode": "center_gaussian",
  "seed": 219,
  "sigma_ratio": 0.001,
  "preserve_fixed": true
}
```

必须：

- fixed 坐标保持原 benchmark；
- movable 使用 deterministic RNG；
- 输出后 fresh audit。

## 13.2 `hpwl_adam`

对应历史 H219 HPWL seed。

核心迭代：

```text
ExactHpwl.evaluate_with_direction()
→ epsilon-active wire direction
→ optimizer.compute_delta()
→ movable coordinate update
→ region clamp
→ optional checkpoint
```

不得引入 density surrogate。

配置建议字段：

```json
{
  "iterations": 200,
  "optimizer": "adam",
  "learning_rate": 0.0,
  "maximum_delta": 0.0,
  "epsilon": 0.0,
  "active_power": 1.0,
  "degree_limit": -1,
  "checkpoint_policy": "best_hpwl"
}
```

其中 `learning_rate/maximum_delta/epsilon` 的实际数值必须从旧代码/旧参数迁移或重新实验，不应由本方案臆造。

`checkpoint_policy`：

```text
last
best_hpwl
```

因为 HPWL-only 阶段不以 overflow 可行为 selector。

## 13.3 `exact_joint_gp`

这是默认纯非光滑 pipeline 中替代“必须依赖 DCT/Poisson warm start”的关键模块。

每轮严格按：

```text
1. exact overlap metrics + overlap active-piece direction
2. exact HPWL + epsilon-active wire direction
3. lambda initialize/update
4. combine directions
5. optional precondition
6. optimizer.compute_delta
7. optional exact batch backtracking
8. commit + clamp
9. trajectory record
```

重要：

- HPWL 数值来自 exact max-min；
- density 数值来自 exact overlap；
- “price field”只允许影响 direction；
- batch acceptance 必须重新审计 exact HPWL / overflow；
- grid resolution 改变时，lambda state 默认 reset/rebalance；
- 若要继承 price/lambda，必须显式 state artifact。

配置字段建议：

```json
{
  "iterations": 81,
  "density_grid": {
    "bins_x": 512,
    "bins_y": 512,
    "target_density": 0.0
  },
  "hpwl_direction": {
    "epsilon": 0.0,
    "active_power": 1.0,
    "degree_limit": -1
  },
  "density_direction": {
    "mode": "exact_overlap_active_piece"
  },
  "lambda": {
    "strategy": "ratio",
    "density_weight_scale": 0.5,
    "reset_on_start": true
  },
  "optimizer": {
    "name": "adam",
    "learning_rate": 0.0,
    "maximum_delta": 0.0
  },
  "exact_batch_acceptance": {
    "enabled": true,
    "max_backtracks": 0
  },
  "selector": {
    "mode": "lowest_hpwl_under_overflow_cap",
    "overflow_cap": 0.07
  }
}
```

`0.5`、`81`、`0.07` 等仅用于映射已有 H253/H254/H257 语义；新路线先在 `adaptec1` 上完成数值验证和调参，效果值得扩大验证后再运行完整八 case。

## 13.4 `exact_recovery`

模块内部使用 `IncrementalDensity`，候选必须 exact-audited。

首版需要覆盖：

- node move；
- axis-separated breakpoint；
- compact direction；
- optional net-rigid block；
- density direction；
- contraction；
- overflow cap；
- line search。

选择原则：

```text
candidate accepted only by exact delta/audit
```

最终 selected layout 再由 framework fresh audit。

配置 preset：

```text
cap15_netblock.json   # H252-like
cap07_node.json       # H255/H257-like
cap07_full.json       # H372 recovery-like
```

## 13.5 `surplus_bisection`

对应 H372 capacity cut。

V1 只迁移代表性路径：

```text
position seeded
surplus only
leaf capacity assignment
nearest capacity
HPWL guided leaf
axis pre-map disabled
```

不要同时迁移所有 bisection 历史开关。

模块应输出：

- moved node count；
- cut level stats；
- capacity before/after；
- iterations/sweeps；
- selected/last layout。

## 13.6 `equal_shape_swap`

对应 H375。

V1 必须支持：

```text
equal-shape only
net-aware neighborhood
candidate radius
candidate count
exact shortlist
sweeps
```

关键不变量：

对于完全同形节点交换：

```text
occupancy field 应保持不变
```

因此测试必须验证交换前后 fresh exact occupancy / overflow 在数值容差内一致。

HPWL 改善必须 exact evaluate。

## 13.7 `historical_dct_poisson`

V1 建议先不编译实现，只放：

```text
code/README.md
params/disabled.json
```

README 明确：

- 来源：历史 H219/H221；
- 依赖：DCT/electric / dreamplace-cpp 相邻实现；
- 属于 smooth surrogate search aid；
- challenge pipeline 禁止；
- 若未来迁移，必须显式定位依赖，不得依赖某台机器的相对目录或未跟踪 exe；
- final `mu=0` 只是一项记录，不能自动证明整个路线符合“不使用光滑化近似”。

---

## 14. H219 → H375 历史链到新 module 的一对一映射

| 历史阶段 | 新项目 module | preset / 说明 | challenge 默认允许 |
|---|---|---|---|
| raw Bookshelf | framework + `layout_init` | raw 或 center gaussian | 是 |
| H219 HPWL seed | `hpwl_adam` | `seed_h219_like.json` | 是 |
| H219 128 DCT/Poisson | `historical_dct_poisson` | coarse historical | 否 |
| H221 512 DCT/Poisson | `historical_dct_poisson` | fine historical, final mu=0 | 否 |
| H252 | `exact_recovery` | cap15 + net block | 是 |
| H253 | `exact_joint_gp` | bridge H253-like | 是 |
| H254 | `exact_joint_gp` | retighten H254-like | 是 |
| H255 | `exact_recovery` | cap07 node, 4 sweeps | 是 |
| H257 | `exact_recovery` | cap07 node, 1 sweep / LS4 | 是 |
| H372 capacity cut | `surplus_bisection` | surplus-only | 是 |
| H372 recovery | `exact_recovery` | cap07 full, 10 sweeps / LS8 | 是 |
| H375 polish | `equal_shape_swap` | net-aware, equal shape | 是 |

新的架构应做到：历史链只需改 pipeline 文件，不再通过 PowerShell 多次启动“全功能程序 + 关闭绝大多数开关”。

---

## 15. Pipeline 只是可组合模块序列，不是固定算法方案

本项目必须支持类似下面的顺序执行：

```text
layout_init
→ hpwl_adam
→ exact_joint_gp
→ exact_recovery
→ surplus_bisection
→ equal_shape_swap
```

但这只是**参考能力链**，不是最终方案，也不能写死在 C++ 中。实际研究时应允许：

```text
A → B → C
A → C
A → B → B → C
C → B
从某个 audited .pl → B → D
只运行单个 module
```

只要模块的输入语义满足，runner 就必须能够顺序交接。每个模块必须满足统一契约：

```text
Problem + input Layout + module config + RunContext
→ selected Layout + stats
```

模块不得通过“前一个模块一定是谁”来决定行为，也不得读取前一模块遗留的隐式全局变量、occupancy cache 或临时目录。若确实需要 Adam moments、lambda state、price field 等 continuation state，只能通过显式 `state_in/state_out` 使用，并且默认 `reset`。

### 15.1 模块交接的双保险

同一 `nsgp run` 进程内，stage 之间优先直接移动/传递内存中的 `Layout`：

```cpp
current = std::move(result.selected_layout);
```

同时，每个成功 stage 必须写一个 fresh-audited `selected.pl`。这个文件有三个用途：

1. 程序异常或人工停止后，可以从最近 stage 边界恢复；
2. 可以把某一阶段输出直接作为另一个实验/模块的输入；
3. 可以独立重放并复核该 stage 的 exact metrics。

因此，**内存中的 Layout 是正常顺序执行的快速交接方式；`selected.pl` 是持久化的可靠交接点**。

### 15.2 正式研究路线仍保持非光滑约束

任何被用于正式非光滑研究的 pipeline 都必须继续满足：

- wire objective/reporting 为 exact weighted pin-offset HPWL；
- density reporting 为 exact rectangle/bin overlap；
- 方向来自原函数 active pieces/subgradient 或 exact-audited discrete candidates；
- 不使用 WA/LSE wirelength 替代原始 HPWL；
- 不使用 smooth density surrogate 替代原始 density；
- 每个 stage 边界 fresh audit。

历史 DCT/Poisson 可以保留为单独 reference/historical module，但不改变上述默认研究边界。

## 16. Pipeline runner 详细流程

伪代码：

```cpp
RunSummary run_pipeline(case_name, pipeline_cfg, cli_overrides) {
    auto loaded = read_bookshelf(benchmark_base(case_name));
    Problem problem = std::move(loaded.problem);
    Layout current = std::move(loaded.raw_layout);

    RunContext base_ctx = build_context(...);

    auto input_metrics = fresh_audit(problem, current, default_density_cfg);
    write_run_header(...);

    for (stage_index, stage_cfg : pipeline.stages) {
        check_module_compliance(stage_cfg.module, pipeline.compliance_mode);

        Json module_cfg = load_and_merge(stage_cfg);
        auto entry = module_registry().at(stage_cfg.module);

        auto before = fresh_audit(problem, current, density_cfg_for_stage_or_global);
        auto input_ref = record_input_reference(...);   // 不复制大文件

        RunContext ctx = make_stage_context(...);
        auto start = steady_clock::now();
        ModuleResult r = entry.runner(problem, current, module_cfg, ctx);
        auto end = steady_clock::now();

        if (r.status != "ok") {
            write_failed_stage(...);
            stop_case();
        }

        auto selected_metrics = fresh_audit(problem, r.selected_layout, ...);
        auto last_metrics = r.selected_is_last
            ? selected_metrics
            : fresh_audit(problem, r.last_layout, ...);

        save audited selected.pl;   // 作为下一模块/后续实验的持久交接点
        if (artifact_policy.keep_last_pl && !selected_is_last) save last.pl;

        write stage_manifest.json;
        append framework stage_metrics.csv;
        append slim trajectory if requested;

        cleanup_stage_temporaries(ctx.stage_dir, artifact_policy);
        release module-local workspace/state not explicitly retained;

        current = std::move(r.selected_layout);
    }

    auto final_metrics = fresh_audit(...);
    write case_summary.json/csv;
    return summary;
}
```

### 16.1 module 异常处理

module 运行中若：

- NaN/Inf；
- vector size 错；
- fixed node 被移动；
- layout 超出 region 且无法 clamp；
- config 非法；
- incremental cache revision 不一致；

应抛出明确异常或返回 failed，case 记录失败，不要继续产出“看似成功”的 final checkpoint。

batch 需要继续执行其余 case，并在总表记录该 case failed，从而满足“无 case 崩溃”的审计需要。

---

## 17. Result directory、关键产物与自动清理

目标不是保存所有运行痕迹，而是只保留**可复现、可画图、可继续衔接**所需的最小结果集。

默认 `artifact_policy = minimal`。推荐目录：

```text
framework/results/<run_id>/
├─ run_manifest.json
├─ stage_metrics.csv
└─ case_summary.csv

modules/<module>/results/<run_id>/<case>/stage_<NN>/
├─ config_snapshot.json
├─ selected.pl              # 必留：fresh-audited，可供后续 stage/实验继续使用
├─ stage_manifest.json      # 必留：模块调用、输入来源、结果、runtime、hash
└─ trajectory.csv           # 按需保留：只保存画收敛曲线需要的标量
```

默认不长期保留：

```text
完整 input.pl 副本
last.pl（除非 selected != last 且研究需要）
每迭代 .pl snapshots
候选 move 列表
临时 occupancy/bin cache
临时 price/gradient/force 数组
调试 dump
solver_state.bin（除非显式要求 continuation）
模块内部 scratch 文件
```

### 17.1 产物生命周期策略

建议配置：

```json
{
  "artifacts": {
    "policy": "minimal",
    "keep_selected_pl": true,
    "keep_last_pl": false,
    "keep_solver_state": false,
    "keep_snapshots": false,
    "trajectory_every": 1,
    "trajectory_fields": [
      "iteration",
      "hpwl",
      "overflow",
      "density_energy",
      "max_density",
      "lambda",
      "accepted",
      "wall_seconds"
    ],
    "cleanup_on_stage_success": true,
    "cleanup_on_run_success": true
  }
}
```

`trajectory_every` 可以按算法成本改成 5/10/N。只有真正用于分析/画曲线的标量才写 CSV，不把大向量序列化。

### 17.2 内存释放规则

“清理中间产物”不仅指磁盘。每个 module 返回后：

- module-local `workspace`、gradient、candidate buffers、incremental occupancy cache 立即析构；
- `ModuleResult` 不得携带大 scratch buffer；
- pipeline 只保留 `current Layout`、必要 auditor 对象和小型 stats；
- `state_out` 只有 `state_policy=resume/emit` 时才保留，否则立刻释放；
- stage 使用不同 bin grid 时，旧 grid workspace 在 stage 结束后释放。

使用 RAII / 局部作用域实现，不依赖手工“最终清理”。

### 17.3 `run_id`

```text
YYYYMMDD-HHMMSS_<pipeline-name>_<short-git>
```

### 17.4 `run_manifest.json`

最少记录：

```json
{
  "run_id": "...",
  "pipeline": "reference_nonsmooth_chain",
  "git_commit": "...",
  "git_dirty": true,
  "compiler": "...",
  "build_type": "Release",
  "real_type": "double",
  "openmp": true,
  "requested_threads": 8,
  "effective_threads": 8,
  "machine": {"hostname": "...", "cpu": "...", "os": "..."},
  "dataset_root": "D:\\codex_project\\HUAWEI_EDA\\alg-electronic\\ispd2005",
  "cases": ["adaptec1"],
  "pipeline_snapshot": "..."
}
```

### 17.5 `stage_manifest.json`

```json
{
  "run_id": "...",
  "case": "adaptec1",
  "stage_index": 3,
  "module": "exact_recovery",
  "status": "ok",
  "seed": 219,
  "threads": 8,
  "input": {
    "source": "previous_stage:selected.pl",
    "sha256": "..."
  },
  "output": {
    "placement": "selected.pl",
    "sha256": "...",
    "metrics": {
      "hpwl": 0.0,
      "overflow": 0.0,
      "density_energy": 0.0,
      "max_density": 0.0
    }
  },
  "before_metrics": {},
  "iterations": 10,
  "accepted_moves": 0,
  "rejected_moves": 0,
  "wall_seconds": 0.0,
  "trajectory": "trajectory.csv",
  "state_out": null,
  "config_snapshot": "config_snapshot.json"
}
```

它必须能回答：“这次结果经过了哪些模块、按什么顺序、各用了什么参数、每个阶段 HPWL/overflow/runtime 如何变化”。

### 17.6 `stage_metrics.csv`

固定列建议：

```text
run_id
case
stage_index
module
status
threads
seed
hpwl_before
hpwl_after
overflow_before
overflow_after
density_energy_before
density_energy_after
max_density_before
max_density_after
iterations
accepted_moves
rejected_moves
module_wall_seconds
audit_wall_seconds
pipeline_cumulative_wall_seconds
input_sha256
output_sha256
config_sha256
```

这张表本身就是最主要的实验结果表。

### 17.7 `case_summary.csv`

单 case 日常运行也写 summary：

```text
case
status
module_sequence
final_hpwl
final_overflow
final_density_energy
final_max_density
total_iterations
module_wall_seconds_sum
audit_wall_seconds_sum
io_wall_seconds
total_wall_seconds
threads
final_pl
```

当前阶段不做任何对照判定；只忠实报告测得的结果。

### 17.8 Hash

对输入 placement（或 benchmark 原始 placement）、config snapshot、`selected.pl` 使用 SHA-256。输入如果就是上一 stage 的 `selected.pl`，记录引用和 hash 即可，不再复制一份文件。

---

## 18. Trajectory：只保存画曲线需要的数据

默认建议开启**轻量 trajectory**，但只写标量 CSV，不保存每次迭代布局。典型字段：

```text
iteration
hpwl
overflow
density_energy
max_density
lambda
step_size
accepted
wall_seconds
```

关键要求：

- trajectory row 和对应状态必须定义在 update 的同一侧；
- 可以每 N 轮记录一次，避免 exact audit 太贵；
- 若 row 中的 HPWL/overflow 标记为 `exact`，必须确实由 exact evaluator 得到；
- 内部 surrogate/direction norm 可以记录，但字段名必须明确，不能冒充 exact metrics；
- `.pl` 迭代 snapshot 默认关闭。

若为了调试临时打开 snapshots，stage 成功后默认只保留最终 `selected.pl`，除非配置显式要求保留调试快照。

---

## 19. 实验运行入口：默认 adaptec1，完整八 case 仅显式触发

日常开发：

```powershell
nsgp.exe run `
  --case adaptec1 `
  --pipeline framework/params/pipelines/reference_nonsmooth_chain.json `
  --threads 8
```

甚至可以让 `--case` 缺省为 `adaptec1`：

```powershell
nsgp.exe run --pipeline <pipeline.json> --threads 8
```

完整实验只在用户明确要求时执行：

```powershell
nsgp.exe batch `
  --all-cases `
  --pipeline <pipeline.json> `
  --threads 8
```

`batch` 只是普通循环，并最终拼接所有 case 的 `case_summary.csv`。当前不做对照比较或 PASS/FAIL gate。

---

## 20. 当前阶段只报告实验结果，不做对照比较

V1/V1.5 的结果系统只负责报告：

```text
执行了哪些 modules / 顺序
每个 stage 的配置与 seed/threads
每个 stage before/after HPWL
每个 stage before/after overflow
density_energy / max_density
iterations / accepted / rejected（适用时）
module runtime
audit runtime
I/O runtime
total runtime
最终 selected.pl
trajectory.csv（需要画收敛曲线时）
status / error（若失败）
```

**不实现对照结果 loader，不输出 degradation、runtime ratio、PASS/FAIL；优化 runner 只负责运行和记录。**

未来如果需要做最终对照比较，再在现有 `case_summary.csv` 之上增加独立 analysis 脚本即可，不要把比较逻辑耦合进优化 runner。

---

## 21. Runtime 计时规范

### 21.1 module runtime

从进入 `runner()` 前到返回后：

```cpp
std::chrono::steady_clock
```

### 21.2 pipeline runtime

记录：

```text
module wall time
fresh audit time
I/O/checkpoint time
pipeline total wall time
```

建议 `stage_metrics.csv` 的 module `wall_seconds` 只记录 module 计算主体；另外在 run manifest 记录 end-to-end total。

这样能避免旧历史中“阶段增量”“checkpoint 起点”“raw-to-final”混为同一个 runtime。

### 21.3 OpenMP

启动时：

```cpp
omp_set_num_threads(requested_threads);
```

并通过 parallel region 实测 effective thread count，写进 manifest。

如果实际线程数与请求不一致，验收 runtime 比较标 invalid。

---

## 22. 数值测试矩阵

测试优先级只围绕数学契约。

### 22.1 HPWL pin offset

构造 2–3 节点 net，手工计算：

```text
node center + pin offset
max-min
net weight
```

断言 exact。

### 22.2 fixed macro occupancy

固定宏覆盖若干 bin：

- movable 不动；
- fresh density 必须包含 fixed occupancy。

### 22.3 `terminal_NI`

将 `terminal_NI` 设置很大 width/height，验证它不增加 occupancy。

### 22.4 density exact rectangle intersection

构造 cell 横跨 bin boundary，手算四个 overlap 面积。

### 22.5 move delta

对 node 随机 move：

```text
incremental predicted after
== fresh full evaluator after commit
```

### 22.6 group delta

同上，至少 2–5 nodes。

### 22.7 checkpoint round trip

```text
Layout → .pl → Layout
```

fresh HPWL/overflow 一致。

### 22.8 selected vs last

构造 dummy module：

```text
selected != last
```

验证 manifest/metrics 不混淆。

### 22.9 equal-shape occupancy invariance

同形节点交换：

```text
fresh occupancy before == after
fresh overflow before == after
```

### 22.10 pipeline smoke

用一个很小 benchmark 或真实 case 的极少 iteration：

```text
layout_init → hpwl_adam → exact_joint_gp
```

确保 config、result dir、fresh audit、manifest 全链路可工作。

### 22.11 case 解析 smoke

不优化，只做：

```text
read → fresh audit → write summary
```

八个 case 必须全通过，作为正式优化前的第一道 gate。

---

## 23. 从旧代码迁移的最小顺序

这里是本地 AI 最应该按顺序执行的部分。

## Phase 0：建立空仓库

创建目录、CMake、JSON 依赖、`nsgp list-modules`。

Definition of Done：

```text
Release build 成功
nsgp list-modules 能输出模块名
ctest 能运行空测试
```

## Phase 1：迁移 `Problem/Layout/Bookshelf`

优先参考旧：

```text
include/epsilon_active/types.hpp
include/epsilon_active/bookshelf.hpp
对应 src
```

只迁移 benchmark 数据与 placement I/O，不迁移 placer config。

DoD：

- 8 cases 全部 parse；
- raw `.pl` roundtrip；
- node 数、net 数、movable/fixed 分类与旧程序一致。

## Phase 2：迁移 exact HPWL

优先保持旧 `ExactHpwl` 计算实现，不做“顺便重写优化”。

DoD：

- 同一 `.pl` 新旧 HPWL 在严格 tolerance 内一致；
- pin offset test 通过；
- degree_limit 不影响 exact value test 通过。

## Phase 3：迁移 exact density

迁移 fixed occupancy、rectangle/bin exact overlap、energy、overflow、max density。

DoD：

- 8 case raw placement 新旧 metrics 对齐；
- fixed macro / terminal_NI test；
- 512 grid 等主要 historical grid 可复算。

## Phase 4：迁移 incremental density

从旧 `evaluate_move/evaluate_group_move/commit_move` 迁移。

DoD：

- random node move delta vs full recompute；
- group move；
- cache revision guard。

## Phase 5：实现 auditor + result protocol

先把模块前后 fresh audit、selected/last、manifest 做出来。

此时先写一个 `identity` 测试 module（可只存在 tests 内）验证 pipeline。

DoD：

- 任意 checkpoint `nsgp audit` 输出统一 metrics；
- stage manifest 可复查；
- selected 与 last 分离。

## Phase 6：迁移 optimizer + `hpwl_adam`

迁移 Adam/AMSGrad/AdaGrad/HeavyBall/SGD；首版可实际只启用 Adam、SGD，其他保持原实现后再开放配置。

DoD：

- H219-like seed 可运行；
- 相同起点/参数/threads 与旧程序 trajectory 前若干步对齐或解释并记录差异；
- exact HPWL 始终由 auditor 复算。

## Phase 7：迁移 `exact_joint_gp`

从旧主 GP loop 提取，而不是迁移整个 `global_place()`。

依次迁移：

```text
exact density direction
exact HPWL direction
lambda controller
optimizer
exact batch acceptance
selector
```

DoD：

- 单模块可以从任意 `.pl` 启动；
- H253/H254-like preset 可以独立运行；
- module 不知道前后其他阶段。

## Phase 8：迁移 `exact_recovery`

先 node-only，再 net-block/compact/contraction。

DoD：

- H255/H257-like node-only 指标变化与旧程序对齐；
- H252/H372-like 复杂配置逐步回归；
- 每次 accepted move 的 overflow guard 使用 exact incremental audit。

## Phase 9：迁移 `surplus_bisection`

只迁移 H372 需要的配置路径。

DoD：

- H372-like 从相同输入 checkpoint 得到可比较输出；
- 不为了兼容历史全部实验迁移几十个开关。

## Phase 10：迁移 `equal_shape_swap`

DoD：

- equal-shape occupancy invariance；
- H375-like HPWL polish 可跑；
- exact shortlist 与最终 exact audit。

## Phase 11：结果汇总 + 可选 full batch

先完成单 case 结果汇总；保留显式 full batch 能力。

DoD：

- `adaptec1` run 能生成完整 `stage_metrics.csv` / `case_summary.csv` / final `selected.pl`；
- full batch 被显式调用时，单个 case 失败不影响其余 case 汇总；
- 不依赖任何外部对照结果；
- stage 成功后临时产物自动清理。

## Phase 12：仅在需要时考虑历史 DCT/Poisson

不要阻塞前 11 个 phase。

---

## 24. 每个迁移阶段的数值回归方法

任何“迁移完成”必须至少有以下一种 reference：

1. 旧 `epsilon_active.exe` 对同一 checkpoint 的 exact 输出；
2. 旧单元测试；
3. 历史 checkpoint + 已知 stage summary；
4. 手工构造的可精确计算 toy case。

优先比较 evaluator，再比较 optimizer trajectory。

### 24.1 推荐 tolerance

不要全局写一个随意 tolerance。

建议：

- HPWL：absolute + relative tolerance；
- overlap area：按 bin area scale；
- overflow：relative/absolute 双阈值；
- checkpoint roundtrip：坐标 tolerance；
- 多线程 reduction：若旧实现存在 reduction 顺序差异，应记录合理浮点 tolerance。

具体数值需在迁移旧实现时依据量级确定，而不是由本方案凭空设定。

---

## 25. `nsgp` CLI 设计

保持小而明确。

### 25.1 单 case run（默认开发方式）

```powershell
nsgp.exe run `
  --case adaptec1 `
  --pipeline framework/params/pipelines/reference_nonsmooth_chain.json `
  --threads 8
```

### 25.2 指定起点 checkpoint

```powershell
nsgp.exe run `
  --case adaptec1 `
  --pipeline framework/params/pipelines/reference_nonsmooth_chain.json `
  --initial-placement D:\path\checkpoint.pl `
  --start-stage 4 `
  --threads 8
```

语义：

- benchmark topology 仍从指定 dataset root 读取；
- initial placement 只覆盖坐标/方向；
- 起点立即 fresh audit；
- solver state 默认 reset。

### 25.3 audit

```powershell
nsgp.exe audit `
  --case adaptec1 `
  --placement D:\path\selected.pl `
  --bins 512 `
  --target-density <value>
```

### 25.4 batch

```powershell
nsgp.exe batch `
  --all-cases `
  --pipeline framework/params/pipelines/reference_nonsmooth_chain.json `
  --threads 8
```

### 25.5 覆盖 dataset root

```powershell
--dataset-root D:\other\ispd2005
```

但默认配置必须仍是项目指定绝对路径。

---

## 26. 代码实现细节：如何避免隐藏状态错误

### 26.1 evaluator 构造原则

- `ExactHpwl` 可长期绑定 `const Problem&`，不缓存 layout-dependent value。
- `ExactDensity` 可缓存 fixed occupancy，但 fresh evaluate 必须重新构造 movable occupancy。
- `IncrementalDensity` 才允许缓存 movable occupancy，而且仅 module-local。

### 26.2 layout revision

每次 module 大批量 update：

```cpp
for (...) {
    layout.x[id] += dx;
    layout.y[id] += dy;
}
layout.touch();
```

incremental object 若不是通过自己的 commit 更新，则应失效并重建。

### 26.3 fixed node 防护

debug build：

```text
module 前保存 fixed 坐标 hash
module 后 assert fixed hash 未变
```

release build 可保留轻量检查或在 fresh auditor 中验证。

### 26.4 NaN gate

每个 module 输出前：

```text
all x/y finite
metrics finite
```

否则 module failed。

---

## 27. 性能原则

### 27.1 不因“解耦”重复昂贵工作

分层不等于每次函数都重新解析 benchmark。

一个 case：

```text
Problem 只读取一次
ExactHpwl 拓扑索引只构造一次
每个 grid 的 fixed occupancy 可构造一次
module 内迭代重复使用 workspace
module 边界才做 fresh audit
```

### 27.2 不跨 module 携带危险 occupancy cache

性能损失只发生在阶段边界的一次 fresh rebuild；相比迭代总成本通常可控，并换来数值可信度。

### 27.3 OpenMP

优先迁移旧 exact core 已经验证的并行策略，不在第一轮迁移时顺便重写并行 reduction。

等 evaluator 数值完全对齐后，再做：

- cache locality；
- thread-local buffer；
- allocation reuse；
- schedule 调优。

---

## 28. Git 与 results 管理

`.gitignore`：

```gitignore
/build/
*.exe
*.pdb

framework/results/**
!framework/results/.gitkeep

modules/*/results/**
!modules/*/results/.gitkeep

# 可选：大型 checkpoint
*.pl.tmp
```

注意：如果需要提交极小 golden test `.pl`，使用 `tests/data/` 并通过 `.gitignore` 例外保留。

应该提交：

```text
framework/code
framework/params
modules/*/code
modules/*/params
tests
README
CMakeLists
```

不提交八个 benchmark 数据和大规模实验输出。

---

## 29. README 最少内容

新仓库 README 不写长篇历史，只写：

1. 研究目标：exact nonsmooth global placement。
2. 非光滑合规边界。
3. 默认 dataset root。
4. build 命令。
5. audit 单 case。
6. run 单 case。
7. 默认 `adaptec1` 的 run；完整八 case batch 作为显式可选命令。
8. 如何新增 module。
9. results 在哪里。
10. historical DCT/Poisson 为什么默认禁用。

---

## 30. 新增一个算法 module 的标准模板

本地 AI 后续创建新算法时，复制：

```text
modules/my_method/
├─ code/
│  ├─ my_method.hpp
│  └─ my_method.cpp
├─ params/
│  └─ default.json
└─ results/.gitkeep
```

`my_method.hpp`：

```cpp
#pragma once
#include "framework/code/module_api.hpp"

namespace nsgp {
ModuleResult run_my_method(
    const Problem& problem,
    const Layout& input,
    const Json& config,
    RunContext& context);
}
```

`my_method.cpp` 骨架：

```cpp
ModuleResult run_my_method(...) {
    auto cfg = parse_config(config);
    Layout work = input;

    // 1. 构造该模块自己的 workspace
    // 2. 生成非光滑方向或 exact-audited candidates
    // 3. 更新 work
    // 4. 记录 module-local stats
    // 5. 不伪造权威 Metrics

    ModuleResult out;
    out.selected_layout = work;
    out.last_layout = work;
    out.selected_is_last = true;
    out.status = "ok";
    return out;
}
```

然后只在 `module_registry.cpp` 加一行。

---

## 31. 建议的首批参数文件职责

参数文件名称要描述“算法意图”，不要继续只用 H 编号。

允许在注释/metadata 里记录历史来源：

```json
{
  "name": "strict_cap07_node_recovery",
  "origin": "H255/H257-like",
  "...": "..."
}
```

推荐：

```text
hpwl_adam/params/
  seed_center_gaussian.json
  smoke.json

exact_joint_gp/params/
  coarse_capacity_build.json
  bridge_medium_density.json
  retighten_strong_density.json
  smoke.json

exact_recovery/params/
  relaxed_cap15_netblock.json
  strict_cap07_node.json
  strict_cap07_full.json

surplus_bisection/params/
  surplus_capacity_cut.json

equal_shape_swap/params/
  net_aware_polish.json
```

H219/H253 等放 `origin`，不要把历史实验号变成新的 API。

---

## 32. 研究过程中建议新增的纯非光滑模块（V2 候选）

只有 V1 稳定后再加。

优先级建议：

1. `bundle_gp`：HPWL + overlap subgradient 的 small-memory bundle/proximal step。
2. `coordinate_breakpoint_density`：exact overlap breakpoint/coordinate search。
3. `capacity_assignment`：exact-audited bin capacity assignment。
4. `coarse_capacity_flow`：若后续多 case overflow 表明 bisection 不够。
5. `transport`：仅迁移能明显改善 overflow/runtime 的一个代表模式。

每个都必须仍遵守统一 `Layout -> Layout` module protocol。

---

## 33. 不要做的错误重构

本地 AI 若出现以下设计，应立即撤回：

### 错误 1：把旧 `PlaceConfig` 原样搬过来

正确：每个 module 只解析自己的参数。

### 错误 2：把 `global_place()` 改名为 PipelineRunner，但内部顺序仍写死

正确：pipeline 配置决定顺序，而且任何示例 pipeline 都只是配置数据；runner 不认识“标准顺序”。

### 错误 3：每个 module 自己重新实现 HPWL/overflow

正确：只调用 framework exact core。

### 错误 4：把 density cache 放进 `Problem`

正确：静态 fixed geometry 可缓存；layout-dependent occupancy 不属于 Problem。

### 错误 5：为了速度省掉 module boundary fresh audit

正确：module 内高性能，module 边界权威复算。

### 错误 6：历史 homotopy 改名后放入 challenge pipeline

正确：compliance tag 硬 gate。

### 错误 7：把结果重新堆进 `experiments/Hxxx` 风格目录

正确：代码/参数固定位置，结果按 run_id 隔离。

### 错误 8：先迁移 477 个实验目录

正确：只迁移能构成 V1 主链的核心实现。

---

## 34. 第一个可运行里程碑

在迁移 recovery/bisection/swap 之前，先做到：

```text
read Bookshelf
→ raw fresh audit
→ center Gaussian
→ hpwl_adam 10 iterations
→ exact_joint_gp 10 iterations
→ selected.pl
→ fresh audit
→ manifest
```

先只在 `adaptec1` 用极少 iterations 跑通；完整八 case smoke 留到用户明确要求或准备完整实验时。

这一步的意义是尽早验证：

- data root；
- I/O；
- exact metrics；
- module protocol；
- config；
- result paths；
- stage failure handling；
- OpenMP threads；
- provenance。

而不是追求 placement 质量。

---

## 35. 第二个里程碑：历史 exact 阶段可重放

加入：

```text
exact_recovery
surplus_bisection
equal_shape_swap
```

从旧 H221/H257/H372 checkpoint 作为 external initial placement 启动（如果这些 checkpoint 可获得），验证：

- H252-like recovery；
- H253/H254-like joint GP；
- H372-like capacity cut + recovery；
- H375-like polish。

目的是证明新 module 边界没有破坏已有 exact 算法。

注意：这是“从历史 checkpoint 重放 exact 阶段”，不是把历史 DCT/Poisson 纳入正式路线。

---

## 36. 第三个里程碑：纯非光滑 raw-to-final

移除任何 historical checkpoint 依赖：

```text
raw Bookshelf
→ reference_nonsmooth_chain
```

先在 `adaptec1` 完整运行并调试。只有效果值得扩大验证时，再显式运行八 case。

每次只输出实验结果表：

```text
case
status
final_hpwl
final_overflow
wall_seconds
threads
iterations
```


---

## 37. 调参策略也要模块化，而不是继续产生海量代码分支

调参只改 JSON：

```text
learning rate
epsilon active parameters
lambda rule
density weight scale
grid
recovery cap
sweeps
line search
swap radius/candidates
```

算法逻辑变化才改 C++。

若需要 sweep，写一个简单 PowerShell/Python 外部脚本生成多个 JSON 并调用 `nsgp.exe`，不要把 sweep engine 写进 core。

---

## 38. 关于历史代表性指标的使用方式

历史 adaptec1 阶段指标可作为 sanity reference，例如：

```text
H219 HPWL seed
H221 selected
H252
H253
H254
H257
H372
H375
```

但它们只能用于：

- 理解阶段作用；
- 检查迁移是否出现数量级错误；
- 对同一历史 checkpoint 做回归。

不能：

- 硬编码进 selector；
- 当成固定目标；
- 当成跨机器可直接比较的 runtime 结论；
- 据此声称 pure pipeline 已达题目指标。

---

## 39. 最终实验报告格式（当前阶段）

`framework/results/<run_id>/case_summary.csv` 是主结果；同时可以自动生成一个轻量 `run_summary.md`：

```markdown
# Run Summary

Case: adaptec1
Pipeline: reference_nonsmooth_chain
Threads: 8
Status: ok
Total runtime: ... s
Final PL: .../selected.pl

| stage | module | HPWL before | HPWL after | overflow before | overflow after | iterations | runtime(s) |
|---:|---|---:|---:|---:|---:|---:|---:|
| 0 | layout_init | ... | ... | ... | ... | ... | ... |
| 1 | hpwl_adam | ... | ... | ... | ... | ... | ... |
...

Final HPWL: ...
Final overflow: ...
Final density energy: ...
Final max density: ...
Total iterations: ...
Total wall time: ...
Trajectory files: ...
```

这里只陈述事实，不做 PASS/FAIL 对照判定。

---

## 40. 本地 AI 的推荐执行清单

按以下顺序逐项完成，不要跳着大规模迁移：

```text
[ ] 1. 新建仓库和目录
[ ] 2. 建 CMake / OpenMP / JSON
[ ] 3. common.hpp / Problem / Layout
[ ] 4. Bookshelf read/write
[ ] 5. adaptec1 parse smoke（完整 case 列表只验证路径存在/可按需读取）
[ ] 6. ExactHpwl + tests
[ ] 7. ExactDensity + tests
[ ] 8. IncrementalDensity + delta tests
[ ] 9. ExactAuditor
[ ] 10. Module API / registry
[ ] 11. Result/provenance writer
[ ] 12. Pipeline runner
[ ] 13. audit CLI
[ ] 14. layout_init
[ ] 15. hpwl_adam
[ ] 16. exact_joint_gp
[ ] 17. 单 case smoke pipeline
[ ] 18. adaptec1 smoke pipeline + stage 间 selected.pl 交接测试
[ ] 19. exact_recovery node-only
[ ] 20. exact_recovery net-block/compact
[ ] 21. surplus_bisection H372 path
[ ] 22. equal_shape_swap H375 path
[ ] 23. 用参考 pipeline 在 adaptec1 raw-to-final 运行
[ ] 24. minimal artifact cleanup + trajectory/result summary
[ ] 25. 从任意 stage selected.pl 重新启动后续 module
[ ] 26. 可选 full batch 入口（不默认执行）
[ ] 27. performance / memory profiling
[ ] 28. 仅在明确需要时研究 historical DCT/Poisson adapter
```

每完成一项，先跑对应数值回归；不要积累十几个未验证改动再一起调试。

---

## 41. 本地 AI 每次修改代码时的自检问题

提交/阶段完成前检查：

1. 我是否改变了 exact HPWL 数学定义？
2. 我是否改变了 fixed macro / `terminal_NI` 的 occupancy 语义？
3. 我是否让 direction 参数意外影响了 reported metric？
4. 我是否让 layout-dependent cache 跨 module 存活？
5. 我是否混淆 selected checkpoint 和 last checkpoint？
6. 我是否新增了一个本可放进 module JSON 的全局 CLI 开关？
7. 我是否把算法逻辑写进了 pipeline runner？
8. 我是否为一个实验新增了不必要的框架层？
9. 我是否用 smooth surrogate 替换了原目标？
10. 我是否在没有 fresh audit 的情况下写了最终 metrics？
11. 当前任务是否真的要求八 case？如果没有，我是否误跑了完整实验浪费时间？
12. stage 完成后是否释放了不需要的内存/scratch，并删除了不需要的中间磁盘产物？
13. 我是否保留了可供后续模块继续使用的 fresh-audited `selected.pl`？
14. 我是否把某个参考 pipeline 的顺序误写成了 runner 的固定逻辑？

任何一项答案不理想，先修再继续。

---

## 42. 需要从旧仓库继续提供给本地 AI 的源码优先级

为了真正执行迁移，而不是重新发明算法，建议按顺序给本地 AI：

```text
1. include/epsilon_active/types.hpp
2. include/epsilon_active/placer.hpp
3. src/placer.cpp
4. src/main.cpp
5. include/epsilon_active/hpwl.hpp + src/hpwl.cpp
6. include/epsilon_active/density.hpp + src/density.cpp
7. include/epsilon_active/recovery.hpp + src/recovery.cpp
8. include/epsilon_active/bisection.hpp + src/bisection.cpp
9. include/epsilon_active/swap_recovery.hpp + src/swap_recovery.cpp
10. optimizer.hpp/cpp
11. lambda_controller.hpp/cpp
12. homotopy_main.cpp（仅历史模块需要）
13. run_h375_transfer.ps1
14. h375_stage_summary.csv
15. tests/core_tests.cpp
```

迁移原则：**优先复制已经验证的数学实现，再改接口；不要同时“重构接口 + 重写算法公式 + 优化性能”。**

---

## 43. 最终项目应达到的开发体验

新研究者需要能做到：

### 新增算法

```text
创建 modules/my_algo/code/my_algo.cpp
创建 modules/my_algo/params/default.json
registry 加一行
pipeline 加一行
```

### 比较一个阶段

```text
输入同一个 audited checkpoint
→ pipeline A: exact_recovery
→ pipeline B: my_algo
→ framework 自动复算 before/after exact metrics
```

### 默认跑 adaptec1

```powershell
nsgp.exe run --case adaptec1 --pipeline <pipeline.json> --threads <N>
```

### 需要时才跑八 case

```powershell
nsgp.exe batch --all-cases --pipeline <pipeline.json> --threads <N>
```

### 一眼看到这次实验发生了什么

```text
framework/results/<run_id>/case_summary.csv
framework/results/<run_id>/run_summary.md
modules/<module>/results/<run_id>/adaptec1/stage_<NN>/selected.pl
```

这才是新项目重构成功的标准：模块可任意组合并可靠交接，实验结果清楚，中间垃圾不持续堆积，而不是“类更多、目录更复杂、接口更通用”。

---

# 附录 A：建议的最小 `module_api.hpp`

```cpp
#pragma once

#include "problem.hpp"
#include "layout.hpp"
#include <nlohmann/json.hpp>
#include <filesystem>
#include <optional>
#include <string>

namespace nsgp {

using Json = nlohmann::json;

struct StateArtifact {
    std::string kind;
    std::filesystem::path path;
    std::string file_sha256;
};

struct RunContext {
    std::string run_id;
    std::string case_name;
    int stage_index = 0;
    int threads = 1;
    std::uint64_t seed = 0;
    std::filesystem::path stage_result_dir;
    bool save_trajectory = false;
    int snapshot_every = 0;
};

struct ModuleStats {
    int iterations = 0;
    int accepted_moves = 0;
    int rejected_moves = 0;
    double wall_seconds = 0.0;
};

struct ModuleResult {
    Layout selected_layout;
    Layout last_layout;
    bool selected_is_last = true;
    ModuleStats stats;
    std::optional<StateArtifact> state_out;
    std::string status = "ok";
    std::string message;
};

using ModuleRunner = ModuleResult (*)(
    const Problem&,
    const Layout&,
    const Json&,
    RunContext&);

} // namespace nsgp
```

---

# 附录 B：建议的 `metrics.hpp`

```cpp
#pragma once
#include "common.hpp"

namespace nsgp {

struct Metrics {
    Real hpwl = 0;
    Real density_energy = 0;
    Real overflow = 0;
    Real max_density = 0;
};

} // namespace nsgp
```

---

# 附录 C：建议的 `reference_nonsmooth_chain.json` 初始模板

注意：此配置用于把代码链跑起来，并不声明参数已满足比赛门槛。涉及旧 brief 未明确给出的数值，应在迁移旧源码/参数后补齐。

```json
{
  "name": "reference_nonsmooth_chain",
  "compliance_mode": "challenge_nonsmooth",
  "dataset_root": "D:\\codex_project\\HUAWEI_EDA\\alg-electronic\\ispd2005",
  "stages": [
    {
      "module": "layout_init",
      "config": "modules/layout_init/params/center_gaussian_219.json"
    },
    {
      "module": "hpwl_adam",
      "config": "modules/hpwl_adam/params/seed_h219_like.json"
    },
    {
      "module": "exact_joint_gp",
      "config": "modules/exact_joint_gp/params/pure_coarse.json"
    },
    {
      "module": "exact_joint_gp",
      "config": "modules/exact_joint_gp/params/bridge_h253_like.json"
    },
    {
      "module": "exact_joint_gp",
      "config": "modules/exact_joint_gp/params/retighten_h254_like.json"
    },
    {
      "module": "surplus_bisection",
      "config": "modules/surplus_bisection/params/h372_like.json"
    },
    {
      "module": "exact_recovery",
      "config": "modules/exact_recovery/params/cap07_full.json"
    },
    {
      "module": "equal_shape_swap",
      "config": "modules/equal_shape_swap/params/h375_like.json"
    }
  ]
}
```

---

# 附录 D：实施决策总结

本方案做出的关键决策如下：

```text
1. 新仓库，不对旧仓库做渐进式大手术。
2. C++17 + OpenMP。
3. JSON 配置，只支持一种格式。
4. Problem / Layout 分离。
5. ExactAuditor 是唯一权威指标出口。
6. IncrementalDensity 只在 module 内存活。
7. module state 默认 reset；需要 continuation 时用显式 artifact。
8. 静态函数表 registry，不做动态插件。
9. selected 与 last 从类型层面分开。
10. challenge compliance 用硬 gate 排除 historical smooth surrogate。
11. 架构保持支持全部八个 case，但日常默认只运行 `adaptec1`；绝不自动触发完整实验。
12. 当前阶段只报告 HPWL、overflow、density、iterations、runtime、module sequence 等实测结果，不实现 baseline/acceptance 逻辑。
13. stage 间以内存 `Layout` 顺序交接，同时持久化 fresh-audited `selected.pl` 作为恢复/复用边界。
14. 只保留关键 trajectory 和交接 checkpoint；其他中间缓存、scratch、debug dumps 默认清理。
15. V1 只迁移主链 exact 模块。
16. 先保证 evaluator 数值一致，再追求 placement 质量和 runtime。
```

如果本地 AI 按此文件执行，第一阶段应优先产出一个**默认能读 `adaptec1`、能 fresh audit、能顺序运行两个非光滑模块、能保存可信 manifest/trajectory/`selected.pl`，并能在 stage 完成后释放无用中间资源的最小项目**；随后再逐个迁移 recovery / bisection / swap。完整八 case 只保留显式入口，不作为日常开发默认动作。
