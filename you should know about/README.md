# non-smooth 项目接手入口

这个目录是所有后续 AI、开发者和实验执行者的强制交接入口。开始工作前，按下面顺序完整阅读：

1. [`00_START_HERE.md`](00_START_HERE.md)：仓库位置、当前分支、最快上手步骤和不可破坏的约束。
2. [`01_PROJECT_REQUIREMENTS.md`](01_PROJECT_REQUIREMENTS.md)：项目目标、数学契约、开发原则和 Git 规则。
3. [`02_ARCHITECTURE_AND_CODE_MAP.md`](02_ARCHITECTURE_AND_CODE_MAP.md)：真实构建关系、顶层目录和核心代码地图。
4. [`03_MODULE_CATALOG.md`](03_MODULE_CATALOG.md)：每个算法模块的输入、输出、代码、参数和状态。
5. [`04_EXPERIMENT_RULES.md`](04_EXPERIMENT_RULES.md)：实验命名、调用链、指标、产物保留和复现规则。
6. [`05_CURRENT_STATE.md`](05_CURRENT_STATE.md)：当前提交、已有结果、已实现能力和已知缺口。
7. [`06_MAINTENANCE_CHECKLIST.md`](06_MAINTENANCE_CHECKLIST.md)：修改前后检查表，以及本目录的同步更新规则。
8. [`07_FILE_INDEX.md`](07_FILE_INDEX.md)：按目录列出的源码、配置、测试和方案索引。
9. [`CHANGELOG.md`](CHANGELOG.md)：面向接手者的架构与实验变更记录。

## 本目录的权威性

- 代码行为以当前工作树为最终事实；这里负责给出可验证的导航，不替代源码。
- `plan/` 描述目标和设计意图，不能据此假定功能已经实现。
- `README.md` 是简短用户说明；本目录是更完整的模型交接说明。
- 一旦源码、CLI、参数结构、模块边界、实验记录格式、权威 checkpoint 或结果发生变化，必须在同一个提交中更新本目录相关文件和 `CHANGELOG.md`。
- 如果发现本文档与源码冲突，先以源码和可复现实验为准，再立即修正文档；禁止让冲突长期存在。

## 一句话项目定义

这是一个面向 ISPD 2005 Bookshelf 数据的 C++17 非光滑全局布局研究项目：使用 exact HPWL 和 exact rectangle/bin overlap 作为权威评价，围绕可组合模块做轻量、可追溯的数值实验，而不是企业级布局软件或完整 legalization/detailed-placement 工具。
