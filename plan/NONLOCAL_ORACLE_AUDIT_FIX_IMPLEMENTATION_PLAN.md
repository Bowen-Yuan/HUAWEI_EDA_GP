# HUAWEI EDA GP Non-local Non-smooth Oracle 修复与重验方案

> 仓库：`Bowen-Yuan/HUAWEI_EDA_GP`  
> 分支：`nonsmooth-gp-v1`  
> 本方案审计基线 HEAD：`e2304255c9f9e9a8c51f695301bb840e7d60cccf`  
> 日期：2026-09-08  
> 用途：交给本地 AI 执行代码修复、契约测试、数值复验和收敛分析。  
>
> **重要：本文件是修改任务书，不是当前源码事实。执行前必须重新读取 `AGENTS.md` 与 `you should know about/` 全部 Markdown，检查当前 branch / HEAD / `git status --short`，避免覆盖用户已有修改。若 HEAD 已变化，以当前 CMake、源码、测试和 fresh exact audit 为最终权威。**

---

# 0. 本轮任务定位

本轮**不是继续增加新算法**，也不是继续大规模调 lambda。

本轮的唯一目标是：

> 修复当前五个 non-local pure-descent module 中已经审计出的共同工程/数学问题，使后续数值实验真正能够回答“这些 oracle 思想是否有效”。

当前五个模块：

```text
density_multiscale_active_gp
density_cut_pressure_gp
density_transport_gp
density_charge_gp
finite_radius_oracle_gp
```

已经满足以下工程目标：

- 独立 stage；
- 无逐步 acceptance；
- 无 backtracking；
- 无 rollback；
- 无 best-feasible restore；
- exact audit 只做 measurement / lambda feedback；
- stage 输出 last applied state。

这些原则继续保持。

但是当前 100 轮结果全部为负：

```text
H375 canonical 起点:
HPWL     85,999,318.311947
Overflow 6.99999030053%

100-round:
multiscale    139,548,754.76 / 14.8103%
cut           156,286,419.76 / 69.0618%
transport     224,641,758.21 / 24.2072%
charge         87,801,002.56 / 87.1644%
finite-radius 195,518,322.69 / 13.8906%
```

**不能用这组结果直接否定五种思想。**

原因是当前实现存在：

1. Cut / Transport / Charge 的 gradient sign 错误；
2. nonlocal common engine 步长尺度严重偏大；
3. `maximum_delta` 没有从 JSON 读入，实际一直使用默认 4.0；
4. auxiliary direction 无论质量和稀疏度都会被强制 RMS normalize；
5. lambda controller 与当前 `exact_joint_gp` 已验证实现发生重复和语义漂移；
6. lambda update interval=5 与最近 V6 证据冲突；
7. Multi-scale 实现没有按照“per-scale normalize 后 ensemble”；
8. Multi-scale auxiliary 重复包含 epsilon=0 exact density；
9. Finite-radius 只 probe 前 4 个 movable、只测 +x/+y；
10. Finite-radius auxiliary 内部已经混合 HPWL+density，外层又二次混合；
11. Transport 当前只是 nearest-deficit 指向，不是 mass-conserving transport；
12. 缺少最基本的 synthetic direction-contract tests。

因此本轮必须先修复“测量仪器”，再重新评价算法。

---

# 1. 不允许改变的项目契约

## 1.1 Exact HPWL

继续使用：

\[
W(x)=
\sum_{e\in E}w_e
\left[
\max X_p-\min X_p+
\max Y_p-\min Y_p
\right].
\]

pin offset 必须保留。

禁止 WA/LSE 替代 final HPWL。

## 1.2 Exact density / overflow

继续使用 exact rectangle/bin intersection：

\[
\rho_b=
\frac{\operatorname{occupancy}_b}{A_b},
\]

\[
D(x)=
\sum_b
\frac12 A_b[\rho_b-\rho_t]_+^2,
\]

\[
O(x)=
\frac{
\sum_b A_b[\rho_b-\rho_t]_+
}{
A_{\rm movable}
}.
\]

保持：

- fixed macro 消耗容量且不移动；
- `terminal_NI` 不消耗物理容量；
- canonical 正式 evaluator = 512×512 / target_density 1.0；
- exact audit 公式不得改变。

---

# 2. 继续坚持 Pure Descent：禁止重新引入 hard acceptance

本轮修复以后仍然必须：

\[
x_{k+1}
=
\Pi_{\rm die}
\left(
x_k-\Delta x_k
\right).
\]

禁止重新出现：

```text
overflow 上升 -> reject
HPWL 上升 -> reject
overflow > 15% -> reject
candidate merit 不降 -> reject
exact backtracking
rollback
restore best feasible
```

15% 只能作为 lambda controller 的参考信息，不是 feasibility gate。

允许：

\[
7\%\rightarrow10\%\rightarrow15\%\rightarrow12\%\rightarrow8\%\rightarrow7\%.
\]

---

# 3. 第一优先级：统一 gradient 与 desired displacement 的符号契约

当前最严重的问题是：

> 各 oracle 返回的是“希望 cell 移动的方向”，但 `run_nonlocal_descent()` 把它当成 gradient，而 optimizer 最终执行 `x -= delta`。

因此必须统一。

## 3.1 新的唯一契约

建议将：

```cpp
using AuxiliaryDirection = ...
```

改名为：

```cpp
using AuxiliaryGradient = ...
```

或者至少在 header 中写清楚：

```cpp
// Contract:
// callback returns a gradient-like vector g_aux.
// The live layout is updated by x <- x - optimizer(g_total).
// Therefore a desired displacement d must be encoded as g_aux = -d.
```

以后所有 oracle 都必须输出：

\[
g_{\rm aux}
\]

而不是：

\[
d_{\rm desired}.
\]

---

# 4. 修复 Cut Pressure 的符号

当前：

\[
S_x(c)=M_L(c)-C_L(c).
\]

若：

\[
S_x(c)>0,
\]

说明 left side overfull，cell 应向 \(+x\) 移动。

因为 engine 执行：

\[
x\leftarrow x-\eta g_x,
\]

因此正确 gradient 必须：

\[
g_x<0.
\]

## 4.1 当前错误

当前近似：

```cpp
ax[id] = prefix_surplus;
```

当 left surplus > 0：

\[
g_x>0
\]

会让 cell 向左，更拥塞。

## 4.2 最小修复

第一轮只修符号，不同时重新设计 cut model：

\[
g_{i,x}^{cut}
=
-\operatorname{prefixSurplus}_x(i),
\]

\[
g_{i,y}^{cut}
=
-\operatorname{prefixSurplus}_y(i).
\]

代码意义：

```cpp
ax[id] = -(bx == 0 ? 0.0 : px[bx - 1]);
ay[id] = -(by == 0 ? 0.0 : py[by - 1]);
```

本轮不要顺手加入 net criticality、macro routing、adaptive cut weight 或 multi-level cut。

---

# 5. 修复 Transport 的符号

当前：

\[
d_T=
p_{\rm sink}-p_{\rm source}
\]

是 desired displacement，但 engine 把它当 gradient。

因此应返回：

\[
g_T=
p_{\rm source}-p_{\rm sink}.
\]

即：

```cpp
ax[id] = source_x - sink_x;
ay[id] = source_y - sink_y;
```

从而：

\[
x\leftarrow x-\eta g_T
=
x+\eta(p_{\rm sink}-p_{\rm source}).
\]

---

# 6. 修复 Charge Field 的符号

定义：

\[
q_b=\operatorname{occ}_b-C_b.
\]

使用：

\[
K_R(r)=[R-r]_+,
\]

\[
H_C(x\mid q)
=
\sum_iA_i\sum_bq_bK_R(r_{ib}).
\]

在 active support 内：

\[
\frac{\partial H_C}{\partial x_i}
=
-q_b\operatorname{sign}(x_i-x_b).
\]

因此：

\[
g_x=
-q_b\operatorname{sign}(dx).
\]

当前代码使用了正号，必须改成：

```cpp
ax[id] -= q * sign(dx);
ay[id] -= q * sign(dy);
```

保证：

```text
positive overflow charge -> repel
negative free-capacity charge -> attract
```

---

# 7. 第二优先级：修复 common engine 的步长尺度

当前：

```cpp
optimizer->compute_delta(
    gradient,
    learning_rate * die_scale,
    maximum_delta * learning_rate * die_scale,
    delta
);
```

默认：

```text
learning_rate = 0.002
maximum_delta = 4.0
```

而 `nonlocal_common.hpp` 没有从 JSON 解析 `maximum_delta`。

因此五模块一直使用 4.0。

---

# 8. 必须让 maximum_delta 进入 JSON

修改：

```text
modules/nonlocal_common.hpp
```

解析：

```json
{
  "optimizer": {
    "name": "amsgrad",
    "learning_rate": 0.002,
    "maximum_delta": 0.25
  }
}
```

映射：

```cpp
c.maximum_delta =
    o.value("maximum_delta", c.maximum_delta);
```

增加参数校验。

---

# 9. 新默认采用已有 V6 证据

第一轮正式复验默认：

```text
optimizer = AMSGrad
learning_rate = 0.002
maximum_delta = 0.25
beta1 = 0.90
beta2 = 0.99
```

不要继续默认 beta2=0.999。

需要强调：

> 同样写 `lr=0.002`，在不同 engine 中并不代表相同物理步长。必须同时对齐 maximum_delta 与方向处理。

---

# 10. 第三优先级：禁止 auxiliary 无条件全局 RMS→1

当前 wire、exact density、auxiliary 都被分别全局 RMS normalize。

这会让：

- 很弱的 auxiliary；
- 很噪的 auxiliary；
- 极稀疏 auxiliary；

都获得和 exact density 接近的强度。

---

# 11. 改成 active-RMS auxiliary scaling

wire 与 exact density 可以继续：

\[
\hat g_W=
\operatorname{RMSNorm}(g_W),
\]

\[
\hat g_D=
\operatorname{RMSNorm}(g_D).
\]

对 auxiliary 计算：

```text
aux_global_rms
aux_active_rms
aux_nonzero_fraction
```

其中：

\[
\operatorname{activeRMS}(g_A)
=
\sqrt{
\frac{
\sum_{i:\|g_{A,i}\|>\epsilon}\|g_{A,i}\|^2
}{
\#\{i:\|g_{A,i}\|>\epsilon\}
}
}.
\]

然后：

\[
\tilde g_A
=
\frac{
g_A
}{
\operatorname{activeRMS}(g_A)+\epsilon
}.
\]

组合：

\[
g_D^{combined}
=
\hat g_D
+
\alpha_A\tilde g_A.
\]

默认：

\[
\alpha_A=0.10.
\]

后续只筛：

\[
\alpha_A
\in
\{0.05,0.10,0.25,0.50\}.
\]

---

# 12. AMSGrad 特殊提醒

AMSGrad/Adam 第一轮近似：

\[
\Delta_i
\approx
\eta\operatorname{sign}(g_i).
\]

因此简单缩放：

\[
g_A\leftarrow0.1g_A
\]

在“该坐标只有 auxiliary 非零”时并不会显著缩小第一步。

所以真正控制幅度的第一手段必须是：

```text
maximum_delta
```

其次才是 aux weight。

不能把“aux_weight 调小”当成步长校准替代品。

---

# 13. 增加 clipping telemetry

每轮记录：

```text
step_rms
step_max
fraction_coordinates_clipped
fraction_cells_moved
```

如果 clipping fraction 长期很高，说明当前算法依赖 clip 而不是自然尺度。

---

# 14. 第四优先级：统一 LambdaController，不继续维护三套

当前实际存在：

```text
exact_joint_gp -> LambdaController
adaptive_lambda_gp -> FunnelLambdaController
nonlocal_descent -> custom trajectory PI(D)
```

这已经产生语义漂移。

---

# 15. Lambda 架构决策

采用：

\[
\boxed{
\text{共享 controller 实现}
+
\text{stage 内持有 controller state}
}
\]

禁止把 per-iteration lambda 更新放进：

```text
main.cpp
microkernel.cpp
```

顶层只负责：

```text
stage dispatch
fresh stage-boundary audit
logging
provenance
```

---

# 16. Shared LambdaController 位置

建议把当前：

```text
modules/exact_joint_gp/code/lambda_controller.cpp
```

中通用实现迁移到：

```text
framework/kernel/src/lambda_controller.cpp
```

header 保持/确认位于：

```text
framework/kernel/include/epsilon_active/lambda_controller.hpp
```

CMake 中：

```text
NUMERIC_KERNEL_SOURCES
```

加入 shared source，并从 `MODULE_ALGORITHM_SOURCES` 移除旧 source path。

迁移时必须先保证 `exact_joint_gp` 数值行为不变。

---

# 17. Generalize trajectory policy，而不是复制 nonlocal controller

给 shared `LambdaConfig` 增加：

```cpp
Real trajectory_start_overflow = -1.0;
```

语义：

```text
<0:
    使用 initial overflow，保持 exact_joint_gp 旧行为

>=0:
    使用显式 soft reference start
```

nonlocal 设置：

```text
trajectory_start_overflow = 0.15
stop_overflow = 0.07
```

15%→7% 只是 soft reference，不是 hard bound。

---

# 18. Nonlocal 统一调用 shared LambdaController

删除 `nonlocal_descent.cpp` 中自己的：

```text
lambda
integral
previous_error
custom PI(D)
```

改为：

```cpp
LambdaController controller(config.lambda);
```

controller 只输出 lambda，不得控制 acceptance。

---

# 19. Lambda update interval

V6 已发现 interval=5 反应偏慢。

nonlocal 第一默认：

```text
update_interval = 2
```

并做小消融：

```text
1 / 2 / 5
```

---

# 20. 第五优先级：修复 Multi-scale 定义

当前 auxiliary 是 raw gradient 直接平均，而且包含 epsilon=0。

修复后：

exact density 已单独存在：

\[
g_D^{exact}=g_D^{(0)}.
\]

auxiliary 默认尺度：

```text
0.25 / 0.5 / 1 / 2 bin
```

不要包含 epsilon=0。

每个 scale：

\[
\hat g_s=
\operatorname{RMSNorm}
(g_D^{(\varepsilon_s)}).
\]

然后：

\[
g_A^{MS}
=
\frac1{N_{nonzero}}
\sum_{s:\|g_s\|>0}
\hat g_s.
\]

最后作为 auxiliary 输入 common engine。

---

# 21. Multi-scale telemetry

至少记录：

```text
nonzero_scale_count
scale_0.25_rms
scale_0.5_rms
scale_1_rms
scale_2_rms
```

可选记录 pairwise cosine，用于判断多尺度是否真的提供新方向。

---

# 22. 第六优先级：Transport 增加最小质量守恒

先修 sign。

如果 20-round 仍出现明显 sink oversubscription，再实现最小 mass-conserving greedy：

维护：

\[
r_s=remaining\ surplus,
\]

\[
r_d=remaining\ deficit.
\]

每次：

\[
f=\min(r_s,r_d),
\]

然后：

\[
r_s\leftarrow r_s-f,
\]

\[
r_d\leftarrow r_d-f.
\]

若 source 分到多个 sink：

\[
v_s=
\frac{
\sum_t f_{st}(p_t-p_s)
}{
\sum_t f_{st}+\epsilon
},
\]

gradient：

\[
g_s=-v_s.
\]

本轮不要直接写重型 min-cost-flow solver。

---

# 23. 第七优先级：Finite-radius 重做为真正 plateau oracle

当前问题：

1. 只 probe 前 4 个 movable；
2. 只 probe +x/+y；
3. 不选择 zero-gradient/stalled cells；
4. auxiliary 内部已经混合 HPWL+density。

---

# 24. 扩展 AuxiliaryContext

建议 common engine 向 oracle 提供：

```cpp
struct AuxiliaryContext {
    int iteration;
    Real lambda;

    const std::vector<Real>& wire_x;
    const std::vector<Real>& wire_y;

    const std::vector<Real>& exact_density_x;
    const std::vector<Real>& exact_density_y;

    Real hpwl;
    DensityMetrics density;
};
```

callback：

```cpp
using AuxiliaryGradient = std::function<void(
    const Database&,
    ExactOverlapDensity&,
    const AuxiliaryContext&,
    std::vector<Real>&,
    std::vector<Real>&)>;
```

---

# 25. Finite probe cell selection

对 cell \(i\)：

\[
r_i=
\sqrt{
g_{W,i,x}^2+g_{W,i,y}^2+
g_{D,i,x}^2+g_{D,i,y}^2
}.
\]

优先 probe \(r_i\) 最小的 movable cells。

默认：

```text
probe_count = 32
```

不再固定 movable ID。

---

# 26. 必须 probe 四方向

```text
+x
-x
+y
-y
```

---

# 27. Finite auxiliary 第一版只做 density secant

不要在 auxiliary 内部再次构造：

\[
W+\lambda D.
\]

第一版只做：

\[
s_{i,+x}^{D}
=
\frac{
D(x+\delta e_{i,x})-D(x)
}{
\delta
}.
\]

其他方向类似。

外层统一：

\[
g_W+\lambda g_D.
\]

如果未来要 finite HPWL oracle，单独设计，不混入当前 density auxiliary。

---

# 28. 第八优先级：新增 synthetic strategy contracts

新增：

```text
tests/test_nonlocal_oracle_contracts.cpp
```

CTest：

```text
nsgp_nonlocal_oracle_contracts
```

这是本轮最重要的测试。

---

# 29. Common engine contracts

必须验证：

### C1 no rollback
一步导致 overflow/HPWL 上升时，live layout 仍保留新位置。

### C2 max delta
每坐标满足：

\[
|\Delta_i|
\le
\eta L_{die}\cdot m.
\]

### C3 JSON maximum_delta
配置 0.25 必须真正进入 engine。

### C4 shared lambda
nonlocal 使用 shared LambdaController。

### C5 last-state handoff
stage 输出等于最后一次 applied state。

---

# 30. Cut contract

左侧 over-capacity、右侧 free：

\[
mean(g_x^{left})<0.
\]

镜像输入方向应镜像。

---

# 31. Transport contract

source 在左，sink 在右：

\[
g_x<0.
\]

一步后：

\[
x_{new}>x_{old}.
\]

若实现质量守恒：

\[
\sum_t f_{st}\le s_s,
\]

\[
\sum_s f_{st}\le d_t.
\]

---

# 32. Charge contract

positive charge 在 cell 左侧：

\[
q>0,\ x_i>x_b
\]

要求：

\[
g_x<0
\]

使 cell 远离。

negative free-capacity charge 在 cell 右侧时，cell 应被吸向右。

若：

\[
r\ge R,
\]

必须：

\[
g=0.
\]

---

# 33. Multi-scale contract

构造：

```text
epsilon=0 local density gradient 很弱/为0
epsilon>0 能感知附近边界
```

要求：

```text
aux multiscale != 0
```

且 epsilon=0 不进入 auxiliary scale list。

---

# 34. Finite-radius contract

构造：

\[
g_D^{exact}=0
\]

但某个 finite displacement：

\[
D(x+\delta e_i)<D(x)
\]

或：

\[
D(x-\delta e_i)<D(x).
\]

要求：

```text
cell 被选中
四方向可 probe
输出非零 finite gradient
```

---

# 35. 第九优先级：完善 telemetry

`trajectory.csv` 每轮建议记录：

```text
hpwl
overflow_percent
density_energy
max_density

lambda_effective
lambda_reference_percent
lambda_control

wire_rms
exact_density_rms

aux_global_rms
aux_active_rms
aux_nonzero_fraction

combined_gradient_rms
wire_density_cosine

learning_rate
maximum_delta
step_rms
step_max
clipped_coordinate_fraction
```

模块特有：

### Cut
```text
max_abs_cut_surplus
mean_abs_cut_surplus
```

### Transport
```text
total_surplus
total_deficit
transported_mass
mean_transport_distance
```

### Charge
```text
positive_charge_mass
negative_charge_mass
active_charge_pairs
```

### Multi-scale
```text
nonzero_scale_count
```

### Finite
```text
probe_count
selected_zero_gradient_count
nonzero_secant_count
```

---

# 36. 第十优先级：建立 exact-only common-engine control

五个 oracle 必须和**同一个 nonlocal common engine**下的 exact-only baseline 比较。

不能只拿 `exact_joint_gp` 当唯一 baseline，因为其 preconditioning、lambda initialization、selector 和 update semantics 不完全相同。

---

# 37. Exact-only control

推荐新增：

```text
framework/params/pipelines/nonlocal_exact_control.json
```

运行 common engine，但：

```json
"direction": {
  "aux_density_weight": 0.0
}
```

如果本地 AI 认为用某个 oracle module 但 aux weight=0 容易误导，可以新增极薄：

```text
nonlocal_exact_control_gp
```

它的 auxiliary 恒为 0，但必须复用同一 common engine。

不要创建第二套 descent implementation。

---

# 38. 修复后的数值实验顺序

## Phase 0 — 纯 contract

先运行：

```text
nsgp_numeric_contracts
nsgp_adaptive_lambda_contracts
nsgp_nonlocal_oracle_contracts
retention contract
```

任何 sign test 失败，禁止跑 H375 正式实验。

---

## Phase 1 — 一步尺度校准

H375 canonical：

```text
SHA:
ba0fd17e9111b5bcbd13ed746ab9857106d962652552d26076c6e82007551e6b

512 x 512
target=1.0
threads=1
seed=219
```

对：

```text
exact-only
5 oracle
```

各跑 1 iteration。

默认：

```text
AMSGrad
lr=0.002
maximum_delta=0.25
beta1=0.90
beta2=0.99
lambda interval=2
```

这一步只做诊断，不做 hard gate。

如果仍出现 7%→20%+，优先继续缩 maximum_delta，不允许加 acceptance。

---

## Phase 2 — step calibration

只用 exact-only common engine：

```text
maximum_delta:
0.125 / 0.25 / 0.5

lr=0.002
20 iterations
```

选择不过度 clip、又有可见移动的量级。

---

## Phase 3 — lambda interval

exact-only + 最佳 step：

```text
interval:
1 / 2 / 5

50 iterations
```

观察 lambda response、maximum overflow、recovery timing、HPWL。

---

## Phase 4 — 五 oracle 20-round corrected screen

固定：

```text
same AMSGrad
same lr
same maximum_delta
same shared lambda
same interval
same evaluator
same seed
```

只改变 oracle。

这是第一次真正有意义的横向比较。

---

## Phase 5 — auxiliary weight

只对健康 oracle：

```text
0.05 / 0.10 / 0.25 / 0.50
```

50 iterations。

同时看：

```text
sign flip fraction
clipped fraction
```

---

## Phase 6 — 100-round primary

晋级配置：

```text
100 iterations
threads=1
```

报告：

```text
last HPWL
last overflow
best-feasible-observed
best-any
first return <=7%
maximum overflow
lambda peak
runtime
```

---

## Phase 7 — optimizer 对比

只有 oracle 本身证明有效后才测试：

```text
AMSGrad
Adam
DualAveraging
```

不要立即做 5×7 大扫。

---

## Phase 8 — Fusion 暂缓

五个 standalone 未修好前禁止 fusion ranking。

至少两个 standalone 出现正信号后，再另行设计：

```text
nonlocal_joint_gp
```

整个 run 共享：

```text
one optimizer state
one LambdaController state
```

五个独立模块保持不动。

---

# 39. Lambda 放哪里：最终架构要求

不把 lambda 放到 microkernel 顶层。

推荐：

```text
framework/kernel
    shared LambdaController implementation
             ↑
             │
continuous optimization stage
    owns LambdaController state
```

也就是：

- controller 实现共享；
- controller state 属于当前连续 stage；
- microkernel 不参与 per-iteration lambda 更新。

---

# 40. 当前项目 lambda 状态的文档修正

修复前：

```text
exact_joint_gp:
    LambdaController

adaptive_lambda_gp:
    FunnelLambdaController

nonlocal_descent:
    custom trajectory PI(D)
```

修复后：

```text
framework/kernel:
    shared LambdaController primitive

exact_joint_gp:
    owns controller state

nonlocal descent:
    owns controller state

adaptive_lambda_gp:
    暂保留旧私有 FunnelLambdaController
    仅维持 V4 historical behavior
```

本轮不要顺手重写 `adaptive_lambda_gp`。

---

# 41. 推荐修改文件

## Shared kernel

```text
framework/kernel/include/epsilon_active/nonlocal_descent.hpp
framework/kernel/src/nonlocal_descent.cpp

framework/kernel/include/epsilon_active/lambda_controller.hpp
framework/kernel/src/lambda_controller.cpp

CMakeLists.txt
```

## Common adapter

```text
modules/nonlocal_common.hpp
```

## Five modules

```text
modules/density_multiscale_active_gp/code/module.cpp
modules/density_cut_pressure_gp/code/module.cpp
modules/density_transport_gp/code/module.cpp
modules/density_charge_gp/code/module.cpp
modules/finite_radius_oracle_gp/code/module.cpp
```

## Params

更新：

```text
modules/<five-module>/params/*.json
```

至少加入：

```text
maximum_delta
beta2=0.99
lambda update_interval
aux_density_weight
```

## Tests

```text
tests/test_nonlocal_oracle_contracts.cpp
```

## Experiment configs

```text
framework/params/pipelines/nonlocal_exact_control.json
framework/params/pipelines/nonlocal_primary_*.json
```

---

# 42. 文档同步

最终至少更新：

```text
you should know about/02_ARCHITECTURE_AND_CODE_MAP.md
you should know about/03_MODULE_CATALOG.md
you should know about/04_EXPERIMENT_RULES.md
you should know about/05_CURRENT_STATE.md
you should know about/07_FILE_INDEX.md
you should know about/CHANGELOG.md
```

如果 trajectory schema 改变，`04_EXPERIMENT_RULES.md` 必须更新。

---

# 43. 本轮禁止顺手做的事情

```text
新增第六/第七种 oracle
smooth HPWL
smooth density
Poisson
DCT
FFT electrostatic
bundle
hard overflow cap
acceptance
backtracking
best-feasible restore
大规模 7 optimizer sweep
修改 exact evaluator
重写 microkernel
第二套 optimizer hierarchy
第二套 Database/Layout
```

---

# 44. Git / 验证要求

执行完成后：

```powershell
git diff --check
git status --short
ctest --test-dir build -C Release --output-on-failure
```

确认：

```text
external checkpoint SHA 前后不变
无 dataset 进入 Git
无 build/
无 exe
无常规 experiment directory
无 .pl
无 snapshot
```

---

# 45. 建议 commit 拆分

## Commit 1

```text
fix: correct nonlocal oracle signs and descent scaling
```

包含：

```text
gradient contract
Cut sign
Transport sign
Charge sign
maximum_delta parsing
telemetry
synthetic contracts
```

## Commit 2

```text
refactor: share lambda controller across continuous GP stages
```

包含：

```text
lambda source relocation
trajectory_start_overflow generalization
nonlocal_descent integration
docs
```

如果无法保证第二步 exact_joint_gp 数值完全不变，必须拆 commit。

---

# 46. 最低成功标准

本轮不要求五个 oracle 立刻赢过 H375。

最低成功标准：

1. Cut / Transport / Charge sign contracts 全通过；
2. `maximum_delta` JSON 真正生效；
3. H375 一步不再默认 7%→20%+ 大爆炸；
4. nonlocal 使用 shared LambdaController；
5. update interval 可配置；
6. Multi-scale 不再重复 epsilon=0；
7. Multi-scale per-scale normalize；
8. Finite-radius 不再只看前 4 个 cell；
9. Finite-radius 支持 ±x/±y；
10. Finite auxiliary 不再内部二次混合 HPWL+density；
11. exact-only common-engine baseline 建立；
12. 20-round corrected oracle comparison 稳定完成；
13. trajectory.csv 足够解释方向尺度与 lambda；
14. 正负结果都写入 handoff。

---

# 47. 可直接交给本地 AI 的提示词

## 执行任务

你正在修改：

```text
Bowen-Yuan/HUAWEI_EDA_GP
branch: nonsmooth-gp-v1
```

首先完整阅读：

```text
AGENTS.md
you should know about/*.md
```

然后执行：

```powershell
git branch --show-current
git rev-parse HEAD
git status --short
git log -5 --oneline
```

不要覆盖用户未提交修改。

---

## 核心目标

修复五个 non-local pure-descent oracle 的实现和共享 engine，使当前负结果能够被重新公平验证。

五个模块：

```text
density_multiscale_active_gp
density_cut_pressure_gp
density_transport_gp
density_charge_gp
finite_radius_oracle_gp
```

必须继续保持独立。

---

## 绝对禁止

不得加入：

```text
overflow hard constraint
overflow decrease acceptance
HPWL decrease acceptance
candidate accept/reject
backtracking
rollback
best-feasible restore
bundle
smooth HPWL
smooth density
DCT/Poisson
```

每一步仍然：

\[
x_{k+1}=\Pi_{die}(x_k-\Delta_k).
\]

exact audit 只用于：

```text
measurement
lambda feedback
logging
```

---

## 必须修复 1：方向符号

统一：

```text
oracle callback 返回 gradient-like vector
layout update 使用 x <- x - optimizer(g)
```

因此：

### Cut
left overfull → desired move right → gradient x 必须为负。

### Transport
source→sink 是 desired move，因此 gradient = source - sink。

### Charge
必须：

```text
positive overflow charge -> repel
negative free-capacity charge -> attract
```

对应：

\[
g=-q\operatorname{sign}(x-x_b).
\]

必须有 synthetic sign tests。

---

## 必须修复 2：步长

`nonlocal_common.hpp` 必须读取：

```json
optimizer.maximum_delta
```

正式默认先采用：

```text
AMSGrad
lr=0.002
maximum_delta=0.25
beta1=0.90
beta2=0.99
```

不要继续隐式使用 maximum_delta=4。

---

## 必须修复 3：auxiliary scaling

禁止把极稀疏 auxiliary 用全局 RMS 强制放大到 1。

至少计算：

```text
aux_global_rms
aux_active_rms
aux_nonzero_fraction
```

auxiliary 使用 active-RMS scaling。

默认：

```text
aux_density_weight=0.10
```

后续小范围筛：

```text
0.05 / 0.10 / 0.25 / 0.50
```

注意 AMSGrad 第一轮对幅度缩放不敏感，真正步长控制必须靠 maximum_delta。

---

## 必须修复 4：lambda 架构

不要把 lambda 放进 main/microkernel 顶层。

把通用 `LambdaController` 作为 shared kernel primitive。

建议迁移：

```text
modules/exact_joint_gp/code/lambda_controller.cpp
```

到：

```text
framework/kernel/src/lambda_controller.cpp
```

保持 exact_joint_gp 旧行为不变。

给 trajectory policy 增加可选：

```text
trajectory_start_overflow
```

默认 `<0` 时保持 initial overflow→stop。

nonlocal 设置：

```text
15% -> 7%
```

这只是 soft reference，不是 hard bound。

nonlocal 默认：

```text
lambda update interval=2
```

并比较：

```text
1 / 2 / 5
```

---

## 必须修复 5：Multi-scale

exact density 已单独存在。

auxiliary 尺度：

```text
0.25 / 0.5 / 1 / 2 bin
```

不要包含 epsilon=0。

每个 scale 先 RMS normalize，再 ensemble。

---

## 必须修复 6：Transport

先修 sign。

若 20-round 仍出现 sink oversubscription，再实现最小 mass-conserving greedy：

\[
f=\min(remaining\ surplus,\ remaining\ deficit).
\]

不要直接引入重型 solver。

---

## 必须修复 7：Finite-radius

不能再固定：

```text
movable_ids[0..3]
```

根据当前 wire/exact-density gradient 选择 low-information / near-zero-gradient cells。

默认至少：

```text
probe_count=32
```

probe：

```text
+x
-x
+y
-y
```

finite auxiliary 第一版只构造 density secant。

不要在 auxiliary 内部再次混合 HPWL+density。

---

## 必须新增测试

新增：

```text
tests/test_nonlocal_oracle_contracts.cpp
```

覆盖：

```text
common no-rollback
maximum_delta JSON
stage output = last state
Cut sign
Transport sign
Charge repel/attract/support
Multi-scale plateau
Finite-radius zero-gradient escape
```

任何 synthetic sign test 不通过，不允许跑 H375。

---

## 正式实验顺序

### Phase 1
exact-only + 五 oracle：

```text
1 iteration
```

### Phase 2
exact-only：

```text
maximum_delta = 0.125 / 0.25 / 0.5
20 iterations
```

### Phase 3
lambda interval：

```text
1 / 2 / 5
50 iterations
```

### Phase 4
五 oracle corrected screen：

```text
20 iterations
```

### Phase 5
健康 oracle：

```text
aux weight = 0.05 / 0.10 / 0.25 / 0.50
```

### Phase 6
晋级配置：

```text
100 iterations
```

只有完成以上，才允许 optimizer/fusion 研究。

---

## 必须画图

至少：

```text
HPWL vs iteration
Overflow vs iteration
Lambda vs iteration
HPWL-Overflow phase trajectory
step_max / clipped_fraction
wire-density cosine
```

---

## 结果记录

无论正负，都更新：

```text
05_CURRENT_STATE.md
CHANGELOG.md
```

pure descent 必须同时报告：

```text
last
best-feasible-observed
best-any
```

不能把 best-any 冒充 final。

---

## 收尾

执行：

```powershell
git diff --check
git status --short
ctest --test-dir build -C Release --output-on-failure
```

确认 external checkpoint SHA 未改变，retention contract 通过。

最终回复必须包含：

```text
commit SHA
修改文件
测试结果
H375 数值结果
仍然存在的问题
```

---

# 48. 最终研究判断原则

只有同时满足：

1. synthetic direction contract 正确；
2. 步长尺度已校准；
3. lambda feedback 使用共享 controller；
4. exact-only common-engine baseline 已建立；
5. oracle 是唯一变量；
6. 无 hard acceptance / rollback；
7. 100-round exact trajectory 完整；

才可以真正评价某个 oracle。

在此之前：

\[
\boxed{
\text{“当前结果为负”}
\neq
\text{“该 oracle 思想无效”}
}
\]

本轮最重要的目标是：

\[
\boxed{
\text{先保证代码真正实现了想测试的数学方向。}
}
\]
