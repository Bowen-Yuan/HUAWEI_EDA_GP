# HUAWEI EDA GP：五类 Non-local Non-smooth Oracle 独立模块与数值实验方案

> 仓库：`Bowen-Yuan/HUAWEI_EDA_GP`  
> 分支：`nonsmooth-gp-v1`  
> 本方案依据的远端 HEAD：`0fbf1b9863420601e1442814e80eac242fc27b41`  
> 日期：2026-09-08  
> 用途：交给本地 AI 执行代码实现、单模块实验、收敛曲线分析与组合实验。  
>
> **重要：本文件是实施计划，不是当前源码事实。真正执行前，本地 AI 必须重新读取 `AGENTS.md` 和 `you should know about/` 全部 Markdown，并重新确认 branch / HEAD / dirty working tree。若仓库已经前进，以当前源码与 fresh exact audit 为最终权威。**

---

# 0. 研究目标

本轮要验证的核心假设不是“换一个 optimizer 就能解决最后优化”，而是：

> 当前 exact non-smooth HPWL / density 的局部次梯度信息过于贫乏，尤其在大量 plateau / zero-subgradient 区域中，optimizer 看不到有价值的非局部移动方向。  
> 如果保持 exact objective 完全不变，但引入 **non-local、non-smooth、可解释的辅助 direction oracle**，是否能够在无硬约束、无逐步接受/拒绝机制下，通过连续梯度式下降 + dynamic lambda，自发完成：
>
> `HPWL exploration → overflow 上升 → lambda 增大 → density recovery → 最终回到约 7%`

本轮实现五个互相独立、可以单独运行的 stage：

1. `density_multiscale_active_gp`
2. `density_cut_pressure_gp`
3. `density_transport_gp`
4. `density_charge_gp`
5. `finite_radius_oracle_gp`

五个模块必须：

- 独立注册；
- 独立 JSON 参数；
- 独立运行；
- 不修改彼此源码；
- 不依赖彼此私有 state；
- 不在模块之间共享 occupancy cache、optimizer moment 或 hidden global state；
- 不把任何一种方案偷偷混进另一种方案的 baseline。

组合实验优先通过 **pipeline 顺序组合** 完成，因此不会污染五个独立模块。

---

# 1. 当前项目约束：本轮必须遵守

## 1.1 exact evaluator 不改变

HPWL 仍为：

\[
W(x)=
\sum_{e\in E} w_e
\left[
\max_{p\in e} X_p-\min_{p\in e} X_p
+
\max_{p\in e} Y_p-\min_{p\in e} Y_p
\right].
\]

其中 pin 坐标必须包含 pin offset：

\[
X_p=x_{\mathrm{node}}+\Delta x_p,\qquad
Y_p=y_{\mathrm{node}}+\Delta y_p.
\]

Density / overflow 仍由 rectangle/bin exact intersection 得到：

\[
\rho_b(x)=
\frac{\operatorname{occ}_b(x)}{A_b},
\]

\[
D(x)=
\sum_b \frac{1}{2}A_b
[\rho_b(x)-\rho_t]_+^2,
\]

\[
O(x)=
\frac{
\sum_b A_b[\rho_b(x)-\rho_t]_+
}{
A_{\mathrm{movable}}
}.
\]

其中：

\[
[z]_+=\max(z,0).
\]

必须保持：

- fixed macro 位置不可改变且消耗容量；
- `terminal_NI` 不消耗物理容量；
- 内存 `Node.x/y` 是中心坐标；
- `.pl` 是左下角坐标；
- exact metric 必须 fresh audit；
- smooth surrogate 不能冒充 final metric。

---

# 2. 本轮最重要的新原则：完全取消逐步 acceptance / hard overflow constraint

这是本方案与当前 `global_view_gp`、当前 `adaptive_lambda_gp` 的最大区别。

## 2.1 禁止以下行为

五个新模块中禁止：

```text
candidate overflow > threshold => reject
candidate overflow must decrease => accept
candidate HPWL must decrease => accept
candidate merit must decrease => accept
exact backtracking until accepted
restore previous layout after a bad step
```

也就是说，不能存在：

\[
x_{k+1}=
\begin{cases}
\tilde x_{k+1}, & \text{if accepted}\\
x_k, & \text{otherwise}
\end{cases}
\]

这种离散 acceptance gate。

---

## 2.2 采用无条件连续更新

每轮都直接：

\[
g_k =
\hat g_{W,k}
+
\lambda_k \hat g_{D,k},
\]

optimizer 计算：

\[
\Delta x_k
=
\operatorname{Optimizer}(g_k),
\]

然后：

\[
x_{k+1}
=
\Pi_{\mathrm{die}}
\left(
x_k-\Delta x_k
\right).
\]

其中：

\[
\Pi_{\mathrm{die}}
\]

只负责 physical placement boundary clamp；它不是 overflow 可行性判断。

每一步移动之后才 fresh audit：

\[
W_{k+1},D_{k+1},O_{k+1}.
\]

这些 exact 指标：

- 用于日志；
- 用于下一轮 dynamic lambda；
- 用于最终实验评价；

**不用于撤销当前 step。**

---

## 2.3 允许 overflow 自由探索

允许轨迹出现：

\[
7\%
\rightarrow
9\%
\rightarrow
12\%
\rightarrow
15\%
\rightarrow
13\%
\rightarrow
10\%
\rightarrow
8\%
\rightarrow
7\%.
\]

甚至如果短时间超过 15%，也不因为一个 hard ceiling 直接 reject。

15% 在本轮的含义是：

> lambda controller 的“高压参考区”，而不是可行域边界。

如果 overflow 过高，算法应该通过：

\[
\lambda\uparrow
\]

自然增强 density direction，而不是通过拒绝布局解决。

---

# 3. 五个模块共同使用的 Pure-Descent 骨架

为了公平比较五种 oracle，五个 stage 必须使用同一套 optimizer / lambda / normalization 语义，只让 `aux_density_direction` 不同。

建议抽取极小的通用 helper，而不是创建新 runner：

```text
framework/kernel/include/epsilon_active/
    trajectory_lambda.hpp          # 可选，共享纯控制逻辑
    direction_utils.hpp            # RMS normalize / cosine / finite check

framework/kernel/src/
    trajectory_lambda.cpp          # 如非 header-only
```

这些 helper 只能提供：

- direction normalization；
- dynamic lambda；
- 数值 finite guard；
- 通用统计；

不能包含任何五种具体策略。

五个 module 各自保留自己的 oracle 实现。

---

## 3.1 共同方向形式

建议每个模块统一使用：

\[
g_D^{(m)}
=
\alpha_{\mathrm{exact}}
\hat g_D^{\mathrm{exact}}
+
\alpha_{\mathrm{aux}}
\hat g_D^{(m,\mathrm{aux})},
\]

再 normalize：

\[
\hat g_D^{(m)}
=
\operatorname{RMSNorm}
\left(
g_D^{(m)}
\right).
\]

wire direction：

\[
\hat g_W=
\operatorname{RMSNorm}(g_W).
\]

最终：

\[
g=
\operatorname{RMSNorm}
\left(
\hat g_W+\lambda\hat g_D^{(m)}
\right).
\]

原因：

1. exact density direction 保留局部“真实性”；
2. auxiliary oracle 只负责补充 non-local information；
3. 各 oracle 的原始量纲差异很大，先 normalize 才能公平讨论 \(\lambda\)；
4. \(\lambda\) 被解释为 wire-vs-density relative pressure，而不是依赖原始梯度量纲的神秘常数。

第一轮默认：

\[
\alpha_{\mathrm{exact}}=1,\qquad
\alpha_{\mathrm{aux}}=1.
\]

后续只对表现明显有价值的模块再做 \(\alpha\) 消融。

---

# 4. Dynamic Lambda：15%→7% 只是参考轨迹，不是约束

当前仓库 V5 已经验证了 `trajectory` lambda + 无 acceptance 机制能够工作，因此五个模块优先复用同一思想。

定义总迭代数 \(T\)，进度：

\[
p_k=\frac{k}{T-1}.
\]

定义 smoothstep：

\[
s(p)=p^2(3-2p).
\]

参考 overflow 轨迹：

\[
\bar O_k
=
O_{\mathrm{end}}
+
(O_{\mathrm{start}}-O_{\mathrm{end}})
[1-s(p_k)].
\]

建议第一轮：

\[
O_{\mathrm{start}}=15\%,\qquad
O_{\mathrm{end}}=7\%.
\]

注意：

\[
O_k>\bar O_k
\]

不会 reject。

只构造误差：

\[
e_k=
\frac{O_k-\bar O_k}{O_{\mathrm{scale}}},
\]

并在 log-domain 调 lambda：

\[
\log\lambda_{k+1}
=
\operatorname{clip}
\left[
\log\lambda_k
+
\eta_p e_k
+
\eta_d(e_k-e_{k-1})
+
\eta_i I_k
\right].
\]

其中：

\[
I_k=\beta_I I_{k-1}+e_k.
\]

直觉：

```text
O << reference
    -> lambda 缓慢下降
    -> HPWL 更自由地收缩

O > reference
    -> lambda 增大
    -> density/global oracle 权重增强

O 快速继续恶化
    -> derivative term 进一步增大 lambda

O 开始恢复
    -> lambda 不会瞬间消失
    -> 保持恢复惯性
```

本轮禁止让 controller 返回“是否接受”。

它只能返回：

```text
lambda_next
```

---

# 5. Optimizer 与步长：必须依据已有数值结果，而不是重新拍脑袋

当前已有两类关键证据。

## 5.1 H375 / V4 acceptance 型实验

在旧的 funnel-acceptance 框架下：

- `dual-averaging + constant` 100 轮曾得到约 `-0.644%` feasible HPWL；
- Adam 能探索但 recovery 较差；
- cosine step 太早缩步，recovery 失败；
- constant 比 cosine 更可靠。

这证明：

> 对最后 recovery 问题，不能默认“越平滑衰减越好”。

---

## 5.2 H221 / V5 无逐步 acceptance 纯下降实验

最新无 acceptance / 无 backtracking 的 trajectory-lambda 实验中：

| optimizer | best-feasible HPWL | overflow |
|---|---:|---:|
| AMSGrad | **99,136,915.82** | 5.4217% |
| Adam | 99,960,810.50 | 5.5886% |
| AdaGrad | 101,850,109.32 | 6.9879% |
| DualAveraging | 105,949,370.15 | 6.9907% |

因此本轮五模块的主 optimizer 顺序建议为：

```text
Primary:
1. AMSGrad
2. Adam

Secondary comparator:
3. DualAveraging

Later / only if needed:
4. AdaGrad
5. SGD
6. HeavyBall
7. NormalizedSGD
```

原因：

> 本轮机制是“无 acceptance 的纯下降”，与 V5 比 V4 更相似，所以 optimizer 的优先级应更重视 V5，而不是机械沿用 V4 的 dual-averaging 排名。

---

## 5.3 步长默认

第一轮不要使用 trust policy，因为 trust 的语义依赖 accept/reject 反馈，而本轮明确取消 acceptance。

默认：

```text
step_policy = constant
```

cosine 只保留为后续 ablation，不作为主实验。

V5 有效的学习率是：

```text
lr = 0.002
```

因此第一轮 pure-descent screen 以它为中心：

\[
\mathrm{lr}
\in
\{0.0005,\;0.001,\;0.002,\;0.004\}.
\]

如果新模块统一做 RMS direction normalization，则这一组学习率更有可比性。

不要在五个模块上各自随意设完全不同的 lr。

---

# 6. 模块一：`density_multiscale_active_gp`

## 6.1 研究问题

Exact rectangle/bin density 的 subgradient 只反映当前 active geometry。

某 cell 在 bin overlap plateau 内时可能：

\[
g_D^{\mathrm{exact}}=0,
\]

但移动有限距离后马上会跨越新的 bin edge。

研究假设：

> 如果 density oracle 同时观察当前位置附近多个 active-face 尺度，就能提前看到下一批潜在 density kink，从而缓解 zero-subgradient 与过度局部的问题。

---

## 6.2 数学思想

定义 enlarged subdifferential：

\[
\partial_\varepsilon D(x)
=
\operatorname{conv}
\left[
\bigcup_{\|z-x\|_\infty\le\varepsilon}
\partial D(z)
\right].
\]

实际不需要连续求整个集合，只使用有限尺度：

\[
\varepsilon_s
\in
\{0,\;0.25h,\;0.5h,\;1h,\;2h\},
\]

其中 \(h\) 是 bin characteristic size。

得到：

\[
g_D^{(\varepsilon_0)},
\ldots,
g_D^{(\varepsilon_S)}.
\]

每个方向先 RMS normalize：

\[
\hat g_D^{(s)}
=
\frac{
g_D^{(\varepsilon_s)}
}{
\operatorname{RMS}(g_D^{(\varepsilon_s)})+\epsilon
}.
\]

然后做简单 convex weighted ensemble：

\[
g_D^{\mathrm{MS}}
=
\sum_s w_s\hat g_D^{(s)},
\qquad
w_s\ge0,\quad
\sum_s w_s=1.
\]

第一版不使用 bundle，不保留 temporal direction history。

---

## 6.3 为什么不属于 smoothing

模块不能定义：

\[
\tilde D(x)
\]

并把它当目标函数。

它只是用多个 exact/non-smooth active-face 方向构造当前 search direction。

final / iteration metric 仍是 canonical exact \(D,O\)。

---

## 6.4 模块目录建议

```text
modules/density_multiscale_active_gp/
├─ code/
│  ├─ multiscale_active.hpp
│  ├─ multiscale_active.cpp
│  └─ module.cpp
├─ params/
│  ├─ smoke.json
│  ├─ amsgrad_lr002.json
│  └─ README.md
└─ results/.gitkeep
```

---

## 6.5 独立性要求

该模块不能调用：

```text
cut pressure
transport
charge field
finite-radius oracle
bundle
```

它的唯一 auxiliary information 是：

```text
multi-scale epsilon-active density
```

这样实验才能回答：

> 单纯扩大 non-smooth active set，究竟有没有价值？

---

# 7. 模块二：`density_cut_pressure_gp`

## 7.1 研究问题

局部 density gradient 不知道：

> 某一大片区域从全局面积守恒角度看是不是“无论如何都塞不下”。

Cut pressure 用全局 prefix capacity imbalance 提供方向。

---

## 7.2 竖直 cut surplus

对第 \(c\) 条竖直 cut：

\[
S_x(c)
=
M_L(c)-C_L(c),
\]

其中：

\[
M_L(c)=
\sum_{b:\,x_b<c}\operatorname{occ}_b,
\]

\[
C_L(c)=
\sum_{b:\,x_b<c}\rho_t A_b.
\]

如果：

\[
S_x(c)>0,
\]

说明 cut 左侧总面积超过左侧总容量。

无论局部 cells 怎么重排，最终都需要有面积跨 cut 向右。

---

## 7.3 cell-level direction

对 cell \(i\)，可定义概念方向：

\[
d^{\mathrm{cut}}_{i,x}
=
A_i
\left[
\sum_{c:x_i<c}[S_x(c)]_+
-
\sum_{c:x_i>c}[-S_x(c)]_+
\right].
\]

类似：

\[
d^{\mathrm{cut}}_{i,y}.
\]

含义：

- cell 位于 overfull cut 的 overfull 一侧 → 推向另一侧；
- cell 位于 free-capacity 一侧 → 不制造无意义排斥；
- 多条 cut 的压力自然形成多尺度 global signal。

最后：

\[
g_D^{\mathrm{cut}}
=
-\;d^{\mathrm{cut}}.
\]

---

## 7.4 特点

该方法不需要：

```text
FFT
Poisson
DCT
Gaussian
smooth density map
```

复杂度可以接近：

\[
O(B+N)
\]

级别的 prefix scan + cell lookup。

---

## 7.5 模块目录

```text
modules/density_cut_pressure_gp/
├─ code/
│  ├─ cut_pressure.hpp
│  ├─ cut_pressure.cpp
│  └─ module.cpp
├─ params/
│  ├─ smoke.json
│  ├─ amsgrad_lr002.json
│  └─ README.md
└─ results/.gitkeep
```

---

## 7.6 第一版不要做的事情

不要加入：

```text
transport sink matching
charge radius
net-aware cut weighting
macro-specific special heuristic
```

先纯粹验证：

\[
\boxed{
\text{global cumulative capacity imbalance 是否能改善 pure descent}
}
\]

---

# 8. 模块三：`density_transport_gp`

## 8.1 研究问题

Local density gradient 通常只告诉 cell：

> 离开这里。

但缺少：

> 应该去哪里。

Transport oracle 显式建立 surplus region 到 deficit region 的 mass flow。

---

## 8.2 surplus / deficit

对 bin \(b\)：

\[
s_b=
[\operatorname{occ}_b-C_b]_+,
\]

\[
d_b=
[C_b-\operatorname{occ}_b]_+,
\]

其中：

\[
C_b=\rho_tA_b.
\]

---

## 8.3 离散 transport

定义：

\[
\pi_{bc}\ge0
\]

表示从 surplus bin \(b\) 到 deficit bin \(c\) 的面积运输量。

第一版使用 Manhattan cost：

\[
c_{bc}
=
|x_b-x_c|+|y_b-y_c|.
\]

概念优化：

\[
\min_\pi
\sum_{b,c}\pi_{bc}c_{bc},
\]

满足：

\[
\sum_c\pi_{bc}\le s_b,
\]

\[
\sum_b\pi_{bc}\le d_c.
\]

不要求第一次就实现最重的 exact min-cost-flow solver。

可以先使用：

- coarse bin；
- greedy nearest deficit；
- auction；
- sparse neighbor graph；

但算法必须明确记录是哪一种。

---

## 8.4 flow → direction

对 source bin \(b\)：

\[
v_b
=
\frac{
\sum_c\pi_{bc}(p_c-p_b)
}{
\sum_c\pi_{bc}+\epsilon
}.
\]

cell \(i\) 位于 source bin \(b\) 时：

\[
d_i^{T}
=
\omega_i v_b.
\]

第一版 \(\omega_i\) 尽量简单，例如：

\[
\omega_i=1
\]

或与 cell area 做规范化。

不要第一版就加入 net criticality，否则无法知道效果来自 transport 还是 HPWL-aware cell selection。

---

## 8.5 与当前 `global_capacity_transport` 的区别

仓库已经有 relocation/transport 类 module，但它主要是：

```text
生成 relocation candidate
→ exact audit
→ 根据接受条件执行
```

新的 `density_transport_gp` 必须是：

```text
occupancy
→ transport vector field
→ gradient-like direction
→ optimizer unconditional step
```

因此不要直接复用旧模块的 candidate acceptance behavior。

可以复用纯数据结构 / cost 计算原语，但不能复用“move only if accepted”的逻辑。

---

## 8.6 模块目录

```text
modules/density_transport_gp/
├─ code/
│  ├─ transport_field.hpp
│  ├─ transport_field.cpp
│  └─ module.cpp
├─ params/
│  ├─ smoke.json
│  ├─ amsgrad_lr002.json
│  └─ README.md
└─ results/.gitkeep
```

---

# 9. 模块四：`density_charge_gp`

## 9.1 研究问题

希望得到类似 electrostatic placement 的远程 capacity awareness，但不引入 smooth electrostatic objective。

核心思想：

\[
\boxed{
\text{capacity imbalance = charge}
}
\]

---

## 9.2 定义 charge

每轮 fresh occupancy 后冻结：

\[
q_b
=
\operatorname{occ}_b-C_b.
\]

其中：

\[
q_b>0
\]

表示 overflow source；

\[
q_b<0
\]

表示 free-capacity sink。

注意：

> 一轮 direction 计算期间把 \(q_b\) 当作 frozen field，不对 \(q_b(x)\) 本身继续求导。

这样 oracle 的物理意义清晰，也避免把 exact overlap derivative 与 charge derivative 重复混合。

---

## 9.3 Non-smooth compact-support kernel

使用 Manhattan 距离：

\[
r_{ib}
=
|x_i-c_b^x|
+
|y_i-c_b^y|.
\]

定义 hinge kernel：

\[
K_R(r)
=
[R-r]_+.
\]

辅助 potential：

\[
H_C(x\mid q)
=
\sum_i A_i
\sum_b
q_b
K_R(r_{ib}).
\]

对 frozen \(q\) 求 subgradient：

\[
g_C
\in
\partial_x H_C(x\mid q).
\]

因此：

- \(q_b>0\)：产生排斥；
- \(q_b<0\)：产生吸引。

---

## 9.4 为什么它不是 smooth electrostatic

没有：

\[
1/r,\qquad
1/r^2,\qquad
\nabla^2\phi=\rho,
\]

没有 Gaussian，也不做 spectral convolution。

只有：

\[
|\cdot|,
\quad
[\cdot]_+,
\]

即 piecewise-linear non-smooth field。

---

## 9.5 半径

第一轮不要大规模多半径搜索。

只测试：

\[
R\in
\{2h,\;4h,\;8h\}.
\]

以 \(4h\) 作为默认。

原因：

- \(R\) 太小退化成局部 density；
- \(R\) 太大会让大量 cell 同时受到近似同质的远程力。

---

## 9.6 模块目录

```text
modules/density_charge_gp/
├─ code/
│  ├─ capacity_charge.hpp
│  ├─ capacity_charge.cpp
│  └─ module.cpp
├─ params/
│  ├─ smoke.json
│  ├─ radius4_amsgrad_lr002.json
│  └─ README.md
└─ results/.gitkeep
```

---

# 10. 模块五：`finite_radius_oracle_gp`

## 10.1 研究问题

如果：

\[
g_W=0,
\qquad
g_D=0,
\]

不代表 finite displacement 没有价值。

它只说明：

> infinitesimal neighborhood 暂时位于同一个 piecewise region。

这个模块要直接利用 **finite-radius exact one-sided secant information** 跨过 plateau。

---

## 10.2 不使用平滑 central finite difference

为了避免把方法解释成“finite-difference smoothing”，第一版不采用：

\[
\frac{F(x+\delta)-F(x-\delta)}{2\delta}
\]

作为主要定义。

而使用 one-sided exact secant：

\[
s_i^+(\delta)
=
\frac{
J_k(x+\delta e_i)-J_k(x)
}{
\delta
},
\]

\[
s_i^-(\delta)
=
\frac{
J_k(x)-J_k(x-\delta e_i)
}{
\delta
}.
\]

其中：

\[
J_k(x)
=
\hat W(x)
+
\lambda_k\hat D(x).
\]

这里 \(\hat W,\hat D\) 只是为了量纲归一化；真实 audit 仍报告原始 exact metrics。

---

## 10.3 breakpoint 优先

优先让 \(\delta\) 接近“下一次几何事件”：

HPWL：

```text
pin 成为 / 离开 net max/min
```

Density：

```text
cell edge 跨过 bin boundary
```

定义有限候选集合：

\[
\mathcal B_i
=
\{
\delta_{i,1},
\delta_{i,2},
\dots
\}.
\]

方向只选择最有利 one-sided secant 的符号与强度：

\[
d_i
=
-\operatorname{signed\_secant}(\mathcal B_i).
\]

**注意：这只是 direction construction。**

布局更新仍然：

```text
optimizer direction
→ unconditional move
```

不能：

```text
probe candidate better => commit
probe candidate worse => rollback
```

---

## 10.4 计算预算

该 oracle 很可能是五个里面最贵的。

第一版应只对：

- combined gradient RMS 很小的 cells；
- zero-gradient cells；
- top-density-risk region；
- 固定数量 sampled movable cells；

做 finite-radius probing。

这是计算预算筛选，不是 acceptance。

---

## 10.5 模块目录

```text
modules/finite_radius_oracle_gp/
├─ code/
│  ├─ finite_radius_oracle.hpp
│  ├─ finite_radius_oracle.cpp
│  └─ module.cpp
├─ params/
│  ├─ smoke.json
│  ├─ amsgrad_lr002.json
│  └─ README.md
└─ results/.gitkeep
```

---

# 11. 五模块统一 JSON 参数面

为了公平，建议共同字段尽量一致：

```json
{
  "iterations": 100,

  "optimizer": {
    "name": "amsgrad",
    "learning_rate": 0.002,
    "beta1": 0.9,
    "beta2": 0.999,
    "numerical_epsilon": 1e-8
  },

  "step": {
    "policy": "constant"
  },

  "lambda": {
    "policy": "trajectory",
    "initial": 0.25,
    "minimum": 0.01,
    "maximum": 256.0,
    "start_overflow_percent": 15.0,
    "stop_overflow_percent": 7.0,
    "update_interval": 5,
    "kp": 0.20,
    "ki": 0.02,
    "kd": 0.10
  },

  "direction": {
    "wire_epsilon_bin_scale": 0.5,
    "exact_density_weight": 1.0,
    "aux_density_weight": 1.0,
    "rms_normalize_components": true
  },

  "logging": {
    "iteration_metrics": true,
    "log_every": 1
  }
}
```

具体默认值应先对照当前 `exact_joint_gp` V5 参数，尽量继承已经验证的 trajectory-lambda 数值，而不是重新创造一套完全不同参数。

五个模块额外只增加自己的：

```text
strategy-specific parameters
```

---

# 12. trajectory.csv：为了画收敛曲线，需要补充通用逐轮 telemetry

当前实验仍保持 metrics-only 三文件：

```text
params.json
experiment.md
trajectory.csv
```

禁止为每个 module 单独产生：

```text
gradient_dump.csv
density_map.csv
debug.csv
```

如果当前 `ExperimentLog` 只能记录 stage boundary，则增加一个非常小的通用：

```cpp
log_iteration(...)
```

接口。

不要新增第四种长期日志文件。

---

## 12.1 建议每轮记录字段

```text
global_iteration
stage
stage_iteration

hpwl
overflow_percent
density_energy
max_density

lambda
overflow_reference_percent

wire_grad_rms
exact_density_grad_rms
aux_density_grad_rms
combined_grad_rms

wire_density_cosine

step_rms
step_max

learning_rate
optimizer
```

模块特有字段尽量少。

例如：

### multiscale active

```text
active_scale_count
```

### cut pressure

```text
cut_pressure_rms
max_abs_cut_surplus
```

### transport

```text
transported_area
mean_transport_distance
```

### charge

```text
charge_rms
active_charge_pairs
```

### finite radius

```text
probed_cells
nonzero_secant_cells
```

---

# 13. 模块输出语义：不 restore、不卡 7%、不回滚

这点必须明确。

每个新 stage 的输出：

\[
x_{\mathrm{out}}=x_T.
\]

即最后一个 unconditional descent state。

不能像现有 `adaptive_lambda_gp` 那样：

```text
search last state
→ restore best feasible
→ stage output
```

这属于隐含 selection，会破坏纯下降轨迹。

---

## 13.1 仍然可以记录 best-feasible

实验分析可以观察：

\[
\min_{k:\,O_k\le7\%} W_k.
\]

但这只是一个 **diagnostic metric**：

```text
best_feasible_observed
```

不能改变 stage 的实际输出布局。

因此最终必须同时报告：

```text
last HPWL / overflow
best-feasible-observed HPWL / overflow
best-any HPWL / overflow
```

不能把 best-feasible-observed 冒充 last。

---

# 14. 代码注册与 CMake

五个 stage 都通过当前 `ModuleRegistry` 注册。

建议新增：

```text
register_density_multiscale_active_gp(...)
register_density_cut_pressure_gp(...)
register_density_transport_gp(...)
register_density_charge_gp(...)
register_finite_radius_oracle_gp(...)
```

只在：

```text
microkernel.hpp
microkernel.cpp
CMakeLists.txt
```

做最小注册/构建改动。

不要在 `main.cpp` 增加针对五种算法的 `if/else`。

---

# 15. 单模块测试

## 15.1 所有模块共通 contract

synthetic 小例验证：

1. fixed node 不移动；
2. `terminal_NI` capacity 语义不变；
3. 每轮 optimizer step 都实际应用；
4. overflow 变高不会触发 rollback；
5. HPWL 变高不会触发 rollback；
6. 没有 backtracking loop；
7. lambda 只根据 trajectory feedback 更新；
8. stage 输出等于最后一次 applied state；
9. fresh exact audit 不改变布局；
10. NaN/Inf 显式报错，而不是静默 rollback。

---

## 15.2 每模块特有测试

### Multi-scale

验证：

\[
\varepsilon=0
\]

退化到局部 active oracle；

更大 \(\varepsilon\) 能让特定 plateau 小例出现非零辅助方向。

### Cut pressure

构造：

```text
左侧 over-capacity
右侧 under-capacity
```

要求大多数左侧 movable 得到向右 pressure。

镜像输入应镜像方向。

### Transport

一个 source + 一个 sink 小例应产生：

\[
source\rightarrow sink
\]

方向。

总 transported mass 不超过：

\[
\min(\sum s_b,\sum d_b).
\]

### Charge

正 charge 应排斥，负 charge 应吸引。

超出 \(R\)：

\[
K_R=0.
\]

### Finite radius

构造 infinitesimal subgradient = 0，但下一个 breakpoint 后 exact merit 改善的小例，要求 secant oracle 输出非零方向。

---

# 16. 数值实验统一起点

主实验从当前 final-optimization 问题的 canonical H375 checkpoint 开始：

```text
D:\codex_project\HUAWEI_EDA\epsilon-active\experiments\
h375_a1_surplus_recovery10_exchange\run\global.pl
```

要求 SHA-256：

```text
ba0fd17e9111b5bcbd13ed746ab9857106d962652552d26076c6e82007551e6b
```

canonical：

```text
case           = adaptec1
grid           = 512 x 512
target density = 1.0
threads        = 1
seed           = 219
```

起始 exact：

```text
HPWL     = 85,999,318.311947
Overflow = 6.99999030053%
```

实验前后必须检查 input SHA 不变。

---

# 17. Experiment 0：纯下降 control

在测试五个新 oracle 前，先得到同环境 control。

使用当前 `exact_joint_gp`：

```text
lambda.policy = trajectory
batch acceptance = disabled
optimizer = AMSGrad
lr = 0.002
iterations = 100
```

只使用：

```text
exact HPWL direction
+
exact density direction
```

不要 auxiliary oracle。

这个 control 是五种方案真正应该比较的 baseline。

原因：

> 当前问题已经变成“无 acceptance 的 pure descent 是否因增加 non-local oracle 而改善”，不能再用 strict-cap `global_view_gp` 当唯一 baseline。

---

# 18. Experiment 1：五模块 20-round health screen

五个模块全部：

```text
optimizer = AMSGrad
lr = 0.002
step = constant
iterations = 20
lambda = same trajectory policy
```

目的：

- finite；
- direction nonzero；
- 没有明显爆炸；
- 每步都应用；
- trajectory logging 正常；
- runtime 粗评估。

不在 20 轮阶段做算法排名。

---

# 19. Experiment 2：五模块 100-round primary comparison

全部固定：

```text
optimizer = AMSGrad
lr = 0.002
iterations = 100
threads = 1
seed = 219
same lambda schedule
same exact evaluator
```

比较：

```text
control exact-only
M1 multiscale-active
M2 cut-pressure
M3 transport
M4 charge
M5 finite-radius
```

Primary metrics：

```text
last exact HPWL
last exact overflow
best-feasible-observed HPWL
first return to <=7% (if trajectory went above)
maximum overflow excursion
runtime
```

辅助：

```text
HPWL reduction at 50 / 100
overflow at 25 / 50 / 75 / 100
lambda peak
gradient cosine
step RMS
```

---

# 20. Experiment 3：Optimizer screen 只对 top modules 做

不要对五模块一开始做 5×7 全排列。

从 Experiment 2 中选最多三个：

```text
top 3 promising oracle
```

测试：

```text
AMSGrad
Adam
DualAveraging
```

理由：

- AMSGrad：V5 pure descent 当前最好；
- Adam：V5 第二；
- DualAveraging：V4 H375 recovery 历史上有明显优势，值得验证这种优势在完全无 acceptance 时是否仍存在。

预算：

```text
100 iterations
same lr first
```

如果 optimizer 的 lr 语义差异明显，再进入 LR screen。

---

# 21. Experiment 4：Learning-rate screen

只对 top 2 oracle × top 2 optimizer。

统一：

\[
lr\in
\{0.0005,\;0.001,\;0.002,\;0.004\}.
\]

每组：

```text
100 iterations
```

不要同时改变：

```text
lambda gains
oracle weight
oracle radius
epsilon scales
```

避免混杂变量。

---

# 22. Experiment 5：Oracle-specific parameter screen

只有前面确认 oracle 本身有效后才做。

## Multi-scale

比较少量：

```text
{0,.25,.5,1}
{0,.5,1,2}
```

## Cut pressure

主要比较：

```text
full-resolution cuts
coarse cuts
```

不要一开始引入很多人工权重。

## Transport

比较：

```text
coarse resolution
transport solver
```

## Charge

比较：

\[
R=\{2h,4h,8h\}.
\]

## Finite-radius

比较：

```text
probe fraction
breakpoint radius
```

---

# 23. 收敛曲线必须绘制

绘图使用独立 analysis script：

```text
scripts/analyze_nonlocal_oracle_study.py
```

只读取 experiment 的：

```text
params.json
trajectory.csv
```

不读取 placement。

图保存到本地：

```text
framework/results/analysis/<study_id>/
```

该目录默认不进入 Git。

---

## 23.1 图 1：HPWL vs iteration

纵轴建议同时提供：

\[
W_k
\]

和 relative：

\[
100\cdot\frac{W_k-W_0}{W_0}.
\]

这样不同 run 易比较。

---

## 23.2 图 2：Overflow vs iteration

画：

\[
100O_k.
\]

并只作为视觉参考画两条虚线：

```text
7%
15%
```

注意：

> 这两条线只是 annotation，不是 algorithm hard bound。

---

## 23.3 图 3：Lambda vs iteration

使用 log y-axis。

观察：

```text
什么时候开始提高 density pressure
是否存在过冲
恢复后 lambda 是否能重新下降
```

---

## 23.4 图 4：HPWL–Overflow phase trajectory

横轴：

\[
O_k
\]

纵轴：

\[
W_k.
\]

这个图比单纯 iteration curve 更重要。

它可以直接显示：

```text
HPWL 探索是否通过更高 overflow 换来
后续 recovery 是否沿更优路径返回 7%
```

---

## 23.5 图 5：Direction conflict

画：

\[
\cos\theta_k
=
\frac{
g_W^Tg_D
}{
\|g_W\|\|g_D\|
}.
\]

如果长期：

\[
\cos\theta\approx-1,
\]

说明 HPWL shrink 与 density spread 几乎完全对抗。

比较五种 oracle 是否能改变这种冲突结构。

---

# 24. 组合实验：不污染五个独立模块

第一阶段组合只允许通过 pipeline sequence。

即：

```text
M_A
→
M_B
```

而不是进入 `M_A` 源码增加：

```text
if use_B ...
```

---

# 25. Fusion F1：Multi-scale → Cut Pressure

假设：

```text
Multi-scale:
先改善 local / near-active geometry

Cut Pressure:
再提供 global capacity redistribution
```

pipeline：

```text
density_multiscale_active_gp(50)
→
density_cut_pressure_gp(50)
```

反向也要测：

```text
cut_pressure(50)
→
multiscale(50)
```

比较 order effect。

---

# 26. Fusion F2：Cut Pressure → Charge

这是我最看好的组合之一。

```text
Cut Pressure:
提供非常低频 / global x-y mass balance

Charge:
提供二维中尺度 sink/source attraction
```

测试：

```text
cut(50) → charge(50)
charge(50) → cut(50)
```

---

# 27. Fusion F3：Cut Pressure → Transport

目的：

```text
cut:
先纠正大区域 mass imbalance

transport:
再做 source-to-destination refinement
```

这样 transport 不需要一开始承担全部 global redistribution。

---

# 28. Fusion F4：Global oracle → Finite-radius refinement

对表现最好的 global oracle：

```text
cut / charge / transport
```

之后：

```text
finite_radius_oracle_gp(25~50)
```

测试 plateau escape 是否适合作为 late refinement。

---

# 29. Fusion F5：Alternating pipeline

如果某两个模块互补，可测试：

```text
A(25)
→ B(25)
→ A(25)
→ B(25)
```

总预算仍为 100。

模块边界：

- fresh exact audit；
- optimizer state reset；
- lambda controller 默认 reset；
- layout state 直接传递。

不共享 hidden solver state。

---

# 30. 第一轮不要做 gradient-level fusion module

不要马上写：

```text
g = g_cut + g_charge + g_transport + ...
```

否则：

- 难以归因；
- 参数迅速爆炸；
- 五个 module 独立性失去意义。

只有当：

1. standalone 已证明至少两个 oracle 有正结果；
2. sequential fusion 证明两者有互补性；

才可以另开第六个独立 module：

```text
nonlocal_oracle_fusion_gp
```

该第六模块必须是新 stage，不能修改五个已有模块。

---

# 31. 公平比较的 budget

Standalone：

```text
100 optimizer updates
```

Fusion：

```text
总计同样 100 updates
```

例如：

```text
A50 + B50
```

而不是：

```text
A100 + B100
```

再与 standalone 100 比。

如果要跑 200，则所有 control/standalone/fusion 都应有同 budget 对照。

---

# 32. 关于 “best feasible” 的评价

由于本轮无 acceptance / 无 restore：

最重要的是：

```text
last state
```

正式比较顺序建议：

### 第一层：last 是否回到目标附近

\[
O_T\le7\%
\]

是理想结果，但不是算法中硬性 enforce 的条件。

### 第二层：last HPWL

如果：

\[
O_T\le7\%
\]

则直接比较 \(W_T\)。

### 第三层：best feasible observed

如果 last 没回 7%，仍报告：

\[
\min_{k:O_k\le7\%} W_k
\]

用于理解轨迹，但不能把它包装成 stage output。

### 第四层：best-any

用于判断 exploration capability：

\[
\min_k W_k
\]

以及对应 overflow。

---

# 33. 主要成功判据

一个 oracle 被认为值得继续，至少应满足以下之一：

## 类型 A：最终可行改进

\[
O_T\le7\%
\]

且：

\[
W_T<W_0
\]

并明显优于 exact-only pure descent control。

这是最强证据。

---

## 类型 B：更优 return trajectory

虽然中期：

\[
O_k>7\%
\]

但能够回到 7%，且返回点 HPWL 明显更好。

这直接验证 non-local exploration 的价值。

---

## 类型 C：显著改善全局方向质量

即使第一轮最终没有更优可行点，但表现为：

- zero-gradient cell 比例下降；
- auxiliary direction 非零覆盖率明显提高；
- wire-density cosine 不再长期接近 -1；
- recovery 速度明显快于 exact-only。

这种结果值得进入第二阶段。

---

# 34. 负结果同样必须记录

例如：

```text
charge field 把 HPWL 迅速推坏
transport 方向波动过大
multiscale active 与 exact direction 高度共线，几乎无新信息
finite-radius runtime 无法接受
cut pressure 只降低低频 imbalance，但 exact overflow 不改善
```

都应写进：

```text
05_CURRENT_STATE.md
CHANGELOG.md
```

不能只保存正结果。

---

# 35. 建议实施顺序

## Phase A — framework 最小准备

1. 重新读取 handoff；
2. 确认 HEAD / dirty；
3. 确认现有 V5 trajectory-lambda 的实际代码；
4. 提取/复用通用 lambda 与 direction-normalization helper；
5. 增加 iteration telemetry API；
6. 不改变现有 module 数值行为。

---

## Phase B — 五个独立 module

建议按风险从低到高：

```text
1. density_multiscale_active_gp
2. density_cut_pressure_gp
3. density_charge_gp
4. density_transport_gp
5. finite_radius_oracle_gp
```

每完成一个：

```text
compile
→ unit contract
→ 1 iter
→ 5 iter
→ 20 iter smoke
```

不要等五个全写完才测试。

---

## Phase C — standalone study

```text
pure exact control
+
5 standalone
```

统一 100 round AMSGrad/lr0.002。

---

## Phase D — optimizer / lr refinement

只筛 top modules。

---

## Phase E — pipeline fusion

只组合 standalone 已证明有价值的模块。

---

# 36. 需要新增/更新的文档

最终代码提交至少同步：

```text
you should know about/02_ARCHITECTURE_AND_CODE_MAP.md
you should know about/03_MODULE_CATALOG.md
you should know about/04_EXPERIMENT_RULES.md
you should know about/05_CURRENT_STATE.md
you should know about/07_FILE_INDEX.md
you should know about/CHANGELOG.md
```

如果 `iteration telemetry` 改变 trajectory schema，则 `04_EXPERIMENT_RULES.md` 必须明确更新。

如果 pure-descent 成为新的正式推荐协议，则 `01_PROJECT_REQUIREMENTS.md` 也应补充：

```text
exact audit 可以仅作为 measurement / feedback；
不要求所有算法都存在 candidate acceptance gate。
```

---

# 37. 不允许出现的实现

本轮禁止：

```text
bundle / temporal bundle
overflow hard cap
overflow decrease acceptance
HPWL decrease acceptance
merit decrease acceptance
exact backtracking acceptance
best-feasible restore
candidate rollback
smooth HPWL objective
smooth density objective
DCT/Poisson/electrostatic final metric
新的 optimizer hierarchy
第二套 Database/Layout
五模块互相 include 私有实现
在某个独立模块中偷偷加入另一个 oracle
```

---

# 38. 允许的数值安全措施

以下不属于“硬 overflow acceptance”：

```text
die boundary clamp
NaN/Inf abort
optimizer epsilon
lambda numeric min/max
finite per-step displacement safety
OpenMP deterministic settings
```

但这些只用于防数值崩溃，不能基于：

```text
overflow better/worse
HPWL better/worse
```

决定 rollback。

---

# 39. 最终希望回答的研究问题

完成本轮后，必须能独立回答：

## Q1

\[
\text{Multi-scale active information 是否能缓解 exact density plateau？}
\]

## Q2

\[
\text{Cumulative cut imbalance 是否提供真正有价值的 global capacity direction？}
\]

## Q3

\[
\text{Explicit surplus→deficit transport 是否优于单纯 local spreading？}
\]

## Q4

\[
\text{Non-smooth capacity charge 是否能复现 electrostatic 的远程 spreading 优势而不做 smoothing？}
\]

## Q5

\[
\text{Finite-radius exact oracle 是否能穿越 zero-subgradient plateau？}
\]

## Q6

\[
\text{哪些 oracle 是互补的，pipeline sequence 能否优于单模块？}
\]

## Q7

在完全无 acceptance / hard overflow constraint 的条件下：

\[
\boxed{
\text{能否仅靠 dynamic }\lambda
+
\text{non-local direction}
+
\text{optimizer}
\text{ 自发回到 }7\%
}
\]

这应当是本轮最重要的最终研究问题。

---

# 40. 最核心的整体结构

最终五个模块虽然 oracle 不同，但共同遵循：

```text
                 exact HPWL
                     │
                     ▼
               wire subgradient
                     │
                     │
exact occupancy ─────┼─────────────┐
                     │             │
                     ▼             ▼
             exact density    module-specific
                direction     non-local oracle
                     │             │
                     └──────┬──────┘
                            ▼
                    density direction
                            │
                            ▼
              trajectory dynamic lambda
                 (15% → 7% reference)
                            │
                            ▼
                RMS-normalized combined
                        direction
                            │
                            ▼
                 AMSGrad / Adam / ...
                            │
                            ▼
                  unconditional update
                            │
                            ▼
                     die clamp only
                            │
                            ▼
                   fresh exact audit
                            │
                 ┌──────────┴──────────┐
                 ▼                     ▼
            metrics/log             lambda feedback
                 │                     │
                 └──────── next iteration
```

这里不存在：

```text
accept?
reject?
rollback?
overflow cap?
```

只有：

\[
\boxed{
\text{measure}
\rightarrow
\text{feedback}
\rightarrow
\text{continuous update}
}
\]

这正是本轮应当保持的算法哲学。
