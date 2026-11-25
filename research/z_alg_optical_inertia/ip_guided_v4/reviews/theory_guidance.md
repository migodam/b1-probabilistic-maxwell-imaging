# 背景场、局部扇区与扰动界

原意见由独立 CLI 研究会话形成。下文保留其数学提案、当时可见数据与未复核判断。作者接受、补算和纠正的结论见[作者回应](RESPONSE_ZH.md)。

以下判断严格基于对工作区原稿/代码/数据的直接读取，以及已完成并落盘于 scratch 目录的复算结果（la_acc.py、exp1_audit_and_array.py、exp2_fullwave_and_multiplicity.py、exp3_deep_probe.py 及日志 task-70.log、task-76.log、task-82.log）。原始研究资料未作任何改动。

---

## 一、已核准的数学命题与符号/同余审计（全部解析证明 + 数值核验至 $10^{-15}$）

记材料分支盒 $X = M + \Delta$，$M = \operatorname{diag}(m_i)$，$D_h = \operatorname{diag}(h_i)$，$|\delta_i| \le h_i$，$q = J - ME = \Delta E$，$\Delta Y = Y - F(M)$，约定 $\operatorname{Herm}(A)=(A+A^*)/2$、$\operatorname{Im}_H(A)=(A-A^*)/(2i)$（与 model.tex 一致）。

### 1. 符号与物理结构核查（已核准）
当 $I - GM$（等价于 $I - MG$）可逆时：
- **精确算子关系**：$E_c = (I-GM)^{-1}B$，$G_c = (I-GM)^{-1}G = G(I-MG)^{-1}$，$S_c = S(I-MG)^{-1}$，$E = E_c + G_c q$，$F(X) - F(M) = S_c q$。
- **互易性与被动辐射保持**：若 $G^T = G$ 且 $M$ 为实对角，则 $G_c^T = G_c$；且
  $$R_c := \operatorname{Im}_H(G_c) = (I - GM)^{-1} R (I - GM)^{-*} \succeq 0 \quad (R = \operatorname{Im}_H G \succeq 0),$$
  即失真背景格林算子 $G_c$ 的辐射部分与自由空间 $R$ **Sylvester 同余**。
- **二次下界精确展开与易错符号**：将非正约束 $Q_{\mathrm{sec}} = q^* W q - E^* D_h^2 W E \preceq 0$（对角 $W \succeq 0$）与等式 $Q_{\mathrm{opt}} = \operatorname{Im}_H(q^* V E) = 0$（实对角 $V$，符号自由）加到 $\mathsf{D} = (S_c q - \Delta Y)^*(S_c q - \Delta Y)$ 上，得 $\mathsf{D} \succeq C_{0,c} + q^* K_c q - q^* L_c - L_c^* q$，其中：
  $$K_c = S_c^* S_c + W - G_c^* D_h^2 W G_c + \operatorname{Im}_H(V G_c), \qquad L_c = S_c^* \Delta Y + G_c^* D_h^2 W E_c + \frac{i}{2} V E_c,$$
  $$C_{0,c} = \Delta Y^* \Delta Y \;\mathbf{-\; E_c^* D_h^2 W E_c}, \qquad H_c = C_{0,c} - L_c^* K_c^{-1} L_c.$$
  **关键符号提醒**：① $C_{0,c}$ 中含负项 $-E_c^* D_h^2 W E_c$（与原 $[0,u]$ 模型 $C_0 = Y^*Y \succeq 0$ 不同）；② $\operatorname{Im}_H(q^* V E_c) = -\operatorname{Herm}(q^*(iVE_c))$，故 $L_c$ 中为 $+\frac{i}{2}VE_c$；③ 单像素处 $(q_i, e_i)$ 满足 $q_i^* q_i \preceq h_i^2 e_i^* e_i$ 且 $\operatorname{Im}_H(q_i^* e_i) = 0$ 充要等价于存在共同实数 $\delta_i \in [-h_i, h_i]$ 使 $q_i = \delta_i e_i$（因 $\operatorname{Im}_H(q_i^* e_i) = -(\operatorname{Im}\delta_i)e_i^* e_i$）。

### 2. 精确坐标同余定理（已核准）
令 $A_M = I - MG$，则 $q = A_M(J - J(M))$ 为可逆仿射变换，且对**任意** $J$（不必物理实现）恒有：
- $Q_{\mathrm{opt}}(q, E) = \operatorname{Im}_H((J - ME)^* V E) \equiv \operatorname{Im}_H(J^* V E)$（因 $E^* M V E$ 为 Hermitian）；
- $Q_{\mathrm{sec}}(q, E) \equiv J^* W J - 2\operatorname{Herm}(J^* M W E) + E^*(M^2 - D_h^2)W E$。

由此得**精确同余**：记未中心化一般区间盒 $[l_i, u_i] = [m_i - h_i, m_i + h_i]$ 的状态阵与源阵为 $K_J, H_J$，则
$$K_J = (I - MG)^* K_c (I - MG), \qquad n_-(K_c) \equiv n_-(K_J), \qquad H_c \equiv H_J \quad (\text{作为 } s\times s \text{ 矩阵严格相等}).$$
**实算核验**（study/array/results/physical_arrays.npz `disks08_low` 全 38 盘、非均匀 $h, W, V$，见 task-70.log）：
- $\|E-(E_c+G_c q)\|_F = 5.264\times 10^{-16}$，$\|(F-F(M))-S_c q\|_F = 6.120\times 10^{-17}$，$\|G_c - G_c^T\|_F = 5.291\times 10^{-16}$，$\|R_c - (I-GM)^{-1}R(I-GM)^{-*}\|_F = 4.850\times 10^{-16}$；
- 一般盒 $[m_i\pm h_i]$ 同余相对误差 $4.280\times 10^{-16}$，$H_c=H_J$ 相对误差 $4.933\times 10^{-15}$，$\lambda_{\min}(\mathsf{D}-\Phi(q)) = 6.655\times 10^{-3} > 0$；
- 特例 $[0, u]$（$m_i = h_i = u/2 \Rightarrow M^2 - D_h^2 = 0$）与原稿 model.tex 的同余误差 $3.746\times 10^{-16}$、$H_c=H_{\mathrm{old}}$ 误差 $3.441\times 10^{-16}$。

**战略推论（极重要）**：在原 $[0, u]$ 盒上改写成中心化失真 Born 坐标，对同一对 $(W, V)$ **数学上零增益**（$\nu$ 和 $H$ 严格不变）；中心化分支盒的真正数学红利全部来自 $l_i = m_i - h_i > 0$（即 $M^2 - D_h^2 = L_b U_b \succ 0$）或局部小半宽 $h_i \ll m_i$。

### 3. 为何单靠 $h\|G_c\|_2 \ge 1$ 不能推出负惯性，以及更锋利的充要刻画（已核准）
$h\|G_c\|_2 < 1$ 只是保证 $W=I, V=0$ 时 $\nu(K_c)=0$ 的**粗糙复非结构充分条件**；反向 $h\|G_c\|_2 \ge 1 \not\Rightarrow n_-(K_c) > 0$ 有四个独立机制，对应四条更锋利的定理（均已数值验证）：

1. **边界等号**：即使 $S_c=0, V=0, W=I$，$n_-(I - h^2 G_c^* G_c) = \#\{j : h\sigma_j(G_c) > 1\}$（严格大于），等号处 $\nu=0$。
2. **光学等式引起的精确移位圆定理（最关键物理机制）**：对均匀 $W=I$、$V = 2\alpha h^2 I$（$\alpha\in\mathbb{R}$），配完全平方恒得
   $$I - h^2 G_c^* G_c + 2\alpha h^2 R_c = (1+\alpha^2 h^2)I - h^2 (G_c - i\alpha I)^*(G_c - i\alpha I),$$
   故**精确恒等式**：
   $$n_-\big(I - h^2 G_c^* G_c + 2\alpha h^2 R_c\big) = \#\Big\{j : \sigma_j(G_c - i\alpha I) > \sqrt{h^{-2}+\alpha^2}\Big\}.$$
   由 Finsler 引理，存在标量光学权 $v \ge 0$ 使 $I - h^2 G_c^* G_c + v R_c \succ 0$ **当且仅当** $h\|G_c|_{\ker R_c}\|_2 < 1$。特别地，在离散开域 TM 网格上只要 $R_c \succ 0$，**无论 $h\|G_c\|_2$ 多大，都存在 $v$ 使 $\nu=0$**（代价是 $v$ 过大时 $L_c$ 会压垮 $H_c$）。
   - **实算核验**（真实 Helmholtz `disks08_high`，$n_{\mathrm{act}}=8, p=6, \|G_c\|_2=0.4335$，36 组测试 100% 吻合上式，见 task-70.log）：当 $h\|G_c\|_2 = 1.5, 2.5, 5.0$ 时，$\nu_{\mathrm{bare}} = 1, 2, 5$，而仅加标量光学移位 $\alpha = 0.2, 1.0, 50.0$ 即全部降为 **$\nu = 0$**。
3. **接收项稳定化定理**：记 $P_W = S_c^* S_c + W \succ 0$，则 $V=0$ 时精确恒有
   $$n_-(S_c^* S_c + W - G_c^* D_h^2 W G_c) = \#\Big\{j : \sigma_j\big(W^{1/2} D_h G_c P_W^{-1/2}\big) > 1\Big\}.$$
   取 $W = \mu I, \mu \downarrow 0$ 极限（叠加 $V = 2\alpha\mu h^2 I$）得精确暗空间判据：
   $$\lim_{\mu\downarrow 0} n_-(K_c) = \#\Big\{j : \sigma_j\big((G_c - i\alpha I)|_{\ker S_c}\big) > \sqrt{h^{-2}+\alpha^2}\Big\}.$$
   - **实算核验**：`disks08_high` 在 $h\|G_c\|_2=1.5$ 时 $\nu_{\mathrm{bare}}=1$，但 $h\|G_c|_{\ker S_c}\|_2 = 0.671 < 1$，不加光学项仅凭 $S_c^* S_c$ 即得 $\nu=0$；在 $24^2=576$ 未知支撑 `r10_twod_0001` 上，即便 $h\|G_c\|_2 = 1.84\text{–}4.41$（$\nu_{\mathrm{bare}}=2$），48 个接收器仍使 $K_c$ 保持 **$\nu=0$**（见 task-82.log）。
4. **对角 $D$-$G$ 缩放与实混合 $\mu$（IQC）刻画 + 精确逆定理**：
   - 无接收无光学时，$n_-(W - G_c^* D_h^2 W G_c) = \#\{j : \sigma_j(W^{1/2} D_h G_c W^{-1/2}) > 1\}$。对非均匀 $h_i$（分支盒常态），实测 $2\times 2$ 互易被动反例满足 $\|D_h G_c\|_2 = 1.5297 > 1$，取 $W = D_h^{-1}$ 即得 $\lambda = [0.2484, 4.2529] > 0 \Rightarrow \nu = 0$。
   - 存在对角 $(W\succ 0, V\text{实})$ 使 $K_c(W,V) \succ 0$，**充要等价于**经典 Fan–Tits–Doyle 实混合 $\mu$ 的 $D$-$G$ 上界（叠加接收稳定化）严格小于 $1$；其精确对偶障碍（仿 pro.txt 定理 4.1）是存在非零 $C \succeq 0$ 满足 $S_c C = 0$、$\operatorname{diag}(C - D_h^2 G_c C G_c^*) \le 0$、$\operatorname{Im}\operatorname{diag}(G_c C) = 0$。物理上，只要盒内存在**接收不可见的实内共振** $(I-\Delta G_c)v=0$（$v\in\ker S_c, |\delta_i|<h_i$），则对**一切** $(W\succ 0, V)$ 必有 $v^* K_c v < 0 \Rightarrow n_-(K_c) \ge 1$。

---

## 二、过度外推的危险（直接来自本轮复算与原稿红线）

1. **危险 1：把中心化 $(q, E)$ 误当 $[0, u]$ 瓶颈的数学突破。** 由同余定理，在 $[0, u]$ 上它与原 contrast.tex 完全等价，不会突破 $96^2$ 原生网格 $u=1$ 的 $d/r \approx 0.301$ 墙。
2. **危险 2：在未裁剪的 $m_i \pm h$ 盒上测试分支排除，导致 $-E_c^* D_h^2 W E_c$ 被虚假负背景吃穿。** 在 `r10_twod_0001`（$n=576$，未知支撑）的错误 CIE 盆地 $m_{\mathrm{CIE}}$ 处，**469/576 格（81.4%）紧贴下界 $0.02$**。未裁剪的 $m_{\mathrm{CIE}} \pm 0.50$ 允许背景取负到 $-0.48$，导致纯中心化尾谱仅 $0.551 r$（失败）；一旦做精确物理裁剪 $[0.02, 20]\cap[m_{\mathrm{CIE}}-h, m_{\mathrm{CIE}}+h]$（背景半宽减半为 $h/2$），同公式立刻从 **$0.551 r$ 跃升至 $1.900 r$**，直接整盒严格排除（task-82.log）。
3. **危险 3：以为 $K_c$ 正定（$\nu=0$）就意味着下界强，或以为只靠单向改进就能普适支配。**
   - 在未裁剪 `r10_twod_0001` $h=0.50$ 处，所有 $w\in[10^{-5}, 0.05]$ 均保持 $\nu=0$，但 $H_c$ 正特征值数随 $w$ 增大从 $7 \to 3 \to 1 \to 0$ 迅速坍缩（因为统一标量 $w$ 的 $-w h^2 E_c^* E_c$ 惩罚淹没了弱照明方向）。
   - 单纯强制把 $q$ 拆成 $q_1 + \eta$（即下节 $\theta=1$）在中/小 $h\|G_c\|$ 下大胜 `CentTail`（最多 24 倍），但在强多次散射（如 `disks16_high` $w=0.80, h\|G_c\|=0.883$ 或 $h\|G_c\|=1.5\text{–}2.0$）处反而会跌回 $0$，而 $\theta=0$ 仍保留正下界（`0.0090` 与 `0.0507`）。不能拿单一端点宣称处处支配。
4. **危险 4：在 $\nu > 0$ 时误把固定源子空间迹 $\operatorname{tr}(Q_s^* H_{\mathrm{ms}}(\delta) Q_s)$ 当作 $\operatorname{Tail}_\nu(H_{\mathrm{ms}}(\delta))$ 的下界。** $\operatorname{Tail}_\nu$ 对应的是**最小** $s-\nu$ 个正特征值之和，Ky Fan 极小原理方向相反，固定 $Q_s$ 给出的是上界而非下界；$\nu > 0$ 时必须对矩阵 $H_{\mathrm{ms}}(\delta)$ 本身作 Loewner 下包后再取 $\operatorname{Tail}_\nu$（或只在 $\nu=0$ 处使用标量迹/分方向定价）。
5. **危险 5：把“局部弱方向/假盆地排除”混同为“消除真实多解”。** 如 study/multiplicity/README.md 所示，同一局部复读数下可存在 $2^L$ 个正则真根；只要某分支盒内包含另一合法真根（或无支撑先验下的等数据材料），任何合法下界都必然 $\le r$，这是信息论底线而非松弛缺陷。

---

## 三、一个有效的物理创意（新下界 + 可证反例 + 判别实验）

### 1. 核心新下界：统一 $\theta$-同伦失真 Born–IQC 下界（Unified $\theta$-Homotopy DB-IQC）
**物理洞察**：原拟用约束 $q^* W q - E^* D_h^2 W E \preceq 0$ 把一阶失真 Born 源 $q_1 = \Delta E_c$ 与高阶多次散射源 $\eta = \Delta G_c q$ 混在一起，从而**丢弃了 $\delta_i \in [-h_i, h_i]$ 在跨 $s$ 束照明间共享同一实标量及其线性箱体几何**；而传统 Neumann 界又在 $h\|G_c\|_2 \ge 1$ 处彻底失效。

**精确恒等式（对任意 $\theta \in [0, 1]$，不需任何小对比度近似）**：由于 $q_i = \delta_i (E_{c,i} + (G_c q)_i)$，对每格恒有 $q_i - \theta\delta_i E_{c,i} = \delta_i((1-\theta)E_{c,i} + (G_c q)_i)$，因此对一切对角 $W \succeq 0$ 与实对角 $V$ 严格成立：
$$(q - \theta \Delta E_c)^* W (q - \theta \Delta E_c) - \big((1-\theta)E_c + G_c q\big)^* D_h^2 W \big((1-\theta)E_c + G_c q\big) \preceq 0,$$
$$\operatorname{Im}_H\!\Big((q - \theta\Delta E_c)^* V \big((1-\theta)E_c + G_c q\big)\Big) = 0.$$
**两个关键数学奇迹**：
1. **$\theta$ 无关的状态矩阵**：展开后关于 $q$ 的二次型矩阵**恒为同一个 $K_c(W, V) = S_c^* S_c + W - G_c^* D_h^2 W G_c + \operatorname{Im}_H(V G_c)$**，完全独立于 $\theta$！一次分解 $K_c$ 即可扫整条 $\theta \in [0, 1]$ 同伦族：
   - $\theta = 0$ 精确退化为原中心化尾谱 `CentTail`（保留在 $h\|G_c\| \gg 1$ 强多次散射区的生存能力）；
   - $\theta = 1$ 把线性项 $S_c \Delta E_c$ 精确锁在真实箱体 $\delta \in [-h, h]^n$ 内，仅对高阶多次散射源 $\eta$ 施以全阶 IQC（不含 Neumann 级数截断，当 $G_c = i\beta$ 时对任意 $h\beta \in (0, \infty)$ **零间隙精确**！）。
2. **两种可计算认证形式（已在 $\nu=0$ 实算验证）**：
   - **中小维已知支撑（DB-IQC-QP）**：对 $q$ 完成平方后得 $\delta$ 的实二次型 $c_0 + 2b^T\delta + \delta^T Q_{\mathrm{eff}}\delta$，用 Shor 对偶 $\max_{u\ge 0:\, Q_{\mathrm{eff}}+\operatorname{diag}(u)\succ 0} \big[c_0 - b^T(Q_{\mathrm{eff}}+\operatorname{diag}(u))^{-1}b - \sum_i u_i h_i^2\big]$ 给出严格下界；
   - **大维未知支撑（DB-IQC-FullK，与 $n$ 无关因子损失）**：把 $n\times n$ 矩阵 $M_{\mathrm{quad}} = \theta^2(W - A_{\mathrm{lin}}^* K_c^{-1} A_{\mathrm{lin}})$ 的负部用谱下界包裹，得一步可算闭式下界 $c_0 - 2\sum_i h_i |b_i| - \lambda_{\max}((-M_{\mathrm{quad}})_+)\|D_h E_c\|_F^2$（仅需与 `CentTail` 同阶的线性求解，无 $n\times n$ SDP）。

### 2. 已算成的核心证据与可证反例（全部通过真值上界一致性断言）

| 算例与设定 | $h\|G_c\|_2$ | 支撑/数据子空间 `Subsp` | `Neum` / `DB+Neum` | 原中心化 `CentTail` ($\theta=0$) | **新 `DB-IQC` ($\theta=1$)** | 真值上界 |
|---|---:|---:|---:|---:|---:|---:|
| **精确 $2\times 2$ 反例** ($G_c=0$, 见 pro.txt 命题 4.2 中心化) | $0.00$ | $0$ | $1.2553$ / $1.4142$ | $\le 1.3656$ ($d^2\le \frac{373}{200}$) | **$1.4142$ ($d^2=2.0000$, 100%精确)** | $1.4142$ |
| **精确 $2\times 2$ 强散射反例** ($G_c=2.5i I_2$, 超 Neumann 半径) | **$1.25 > 1$** | **$0$** | **$0$ / $0$** | $0.1819$ ($d^2=0.0331$, 全对角最优) | **$0.2055$ ($d^2=0.0422$, +27.7%)** | $0.2276$ |
| 已知支撑 `disks12_low` $w=0.80$（同支撑目标） | $0.156$ | $0$ | $0$ / $0.0028$ | $0.0003$ | **$0.0072$（24× `CentTail`）** | $0.0130$ |
| 已知支撑 `disks12_high` $w=0.80$（同支撑目标） | $0.431$ | $0$ | $0$ / $0$ | $0.0004$ | **$0.0023$（5.7× `CentTail`）** | $0.0260$ |
| 已知支撑 `disks16_high` $w=0.40 \to 0.80$（$\theta$-同伦互补） | $0.441 \to 0.883$ | $0$ | $0$ / $0$ | $0.0222 \to \mathbf{0.0090}$ | $\mathbf{0.0253} \to 0$（取 $\max_\theta$ 严格支配） | $0.0326 / 0.0659$ |
| 已知支撑 `disks16_high` 超 Neumann 半径 | **$1.5 \to 2.0$** | **$0$** | **$0$ / $0$** | **$0.0986(\nu=1) \to 0.0507(\nu=0)$** | $0.0147 \to 0$（$\theta=0$ 接管） | $\|\Delta Y\|=0.2139$ |
| **未知支撑全图 $24^2=576$** `r10_twod_0001` 绕错盆地 $m_{\mathrm{CIE}}$ ($h=0.15$, 未裁剪) | $0.552$ | $0$ | $0$ / — | $5.969\,r$ | **$10.432\,r$（+74.8%）** | $\le 13.04\,r$ |
| **未知支撑全图 $24^2=576$** `r10_twod_0001` 绕 $m_{\mathrm{CIE}}$ ($h=0.30$, **超 Neumann!**) | **$1.104 > 1$** | **$0$** | **$0$ / —** | $3.139\,r$ ($\nu=0$) | **$3.768\,r$（+20.0%，整盒排除!）** | $\le 13.04\,r$ |
| **未知支撑全图 $24^2=576$** `r10_twod_0001` **物理裁剪盒** ($h=0.50$) | $0.839$ | $0$ | $0$ / — | **$1.900\,r$（整盒排除!）** | 未定（见下节） | $\le 13.04\,r$ |

**为何能同时打赢审稿人关心的三大基线**：
- **打赢“精确支撑/数据子空间 Gram”**：同支撑材料变更或未知支撑（$n > p$）下 $\Delta Y \in \operatorname{ran}(S_c)$，`Subsp` 恒为 $0$，而 $\theta$-DB-IQC 给出严格正下界并完成整盒排除。
- **打赢“Neumann 源/收约束”**：当 $h\|G_c\|_2 \ge 1$（如 $2\times 2$ $\beta=2.5$、`disks16_high` $1.1\text{–}2.0$、`r10_twod_0001` $h=0.30$）所有 Neumann 界归零，而光学移位圆 + 接收稳定化 + $\theta$-DB-IQC 仍能整盒排除；即使 $h\|G_c\|_2 < 1$（如 `disks12_low/high`），`DB-IQC-QP` 也比一阶 `DB+Neum` 高 2.6 倍至无穷倍。
- **打赢“直接全波 Top-K”的正确定位**：不宣称单次求解比 Top-K 更快找解，而是**零漏杀（0% false pruning）的整簇剪枝认证**——在多根/折叠盒中，Top-K 未找到可行点无法证明盒内无根，而单次 $K_c$ 分解即可排除围绕错误局部极小 $m_{\mathrm{CIE}}$ 的整个 576 维分支盒（如 $h=0.30, 0.50$），并在分支定界中直接砍掉整棵子树。

### 3. 推荐的一条决定性判别实验（含已定与未定边界）
**实验设计（孪生多解 + 未知支撑假盆地双臂判别）**：
- **臂 1（已知支撑 + 真实双解/近双解防误杀与跨口整盒剪枝）**：取 study/array/results/physical_arrays.npz 冻结的 `null_self` 孪生对 $(m_-, m_+)$（以及 study/multiplicity/results/materials_and_ports.npz 的 $L=2$ 四根 `AA/AB/BA/BB`）。在仅用本地/3 通道 `base` 臂时验证证书**绝不误杀**含另一真根的分支盒 $B(m_+, h)$；切换至含新增空间角/跨组端口的 `rich` 臂、以 $Y = F(m_-)$ 为读数时，统一对比 `Subsp`、`DB+Neum`、`CentTail` 与 `θ-DB-IQC` 在同一 $h$ 阶梯上能够整盒淘汰的最大半径 $h_{\max}$，并与等算力预算的直接全波 Top-K 对比假淘汰率与剪枝成本。
- **臂 2（未知支撑大域假盆地整盒淘汰）**：在 `r10_twod_0001`（$24^2$）与 `austria_eps4_400MHz`（$32^2$）上，围绕历史错误收敛点（如 $m_{\mathrm{CIE}}$）在**物理裁剪盒** $[l_i, u_i]$ 上比较 `CentTail`、`θ-DB-IQC-FullK` 与**相干方向独立定价** $\sum_{j=1}^s \max_{w_j, v_j: \nu_j=0} (a_j^* H_c(w_j, v_j) a_j)_+$（延续 contrast.tex 命题 2 的思想，专治统一标量 $w$ 下弱照明拖垮强照明的问题）。

---

## 四、明确区分：未完成/未运行项（逐项标作“未定”）

以下项目已提出数学构造或脚本设计，但**尚未完成原生数值运行，现阶段一律标作未定（猜想/待验证）**：
1. **相干方向独立定价 + 物理裁剪 `DB-IQC-FullK` 在 `r10_twod_0001` 裁剪盒 $h=0.50, 0.80, 1.00$ 上的实际增益**（`exp4_decisive_tests.py` 的 EXP 4A，未运行）：能否把裁剪 $h=0.80$（$\|D_h G_c\|=1.159 > 1$）从当前 `CentTail` 的 $0.735 r$ 推过 $1.0 r$，**未定**。
2. **`study/array` 冻结 `null_self` 孪生对在 `base` vs `rich` 臂上的全盒剪枝最大半径**（`exp4_decisive_tests.py` 的 EXP 4C，未运行）及 **`study/multiplicity` $L=2$ 四根的跨口整盒剪枝**：**未定**。
3. **非正定情形 $\nu(K_c) > 0$ 下 $\theta > 0$ 的 Loewner 矩阵下包版本（$\operatorname{Tail}_\nu(H_\theta(\delta))$ 的可算紧下界）**：数学形式已明确，但数值紧度**未定**。
4. **原生 $96^2=9216$ 大域无损算例在局部候选中心 $M \neq 0$ 附近的分支盒排除收益与单侧 Arb 区间认证**：**未定**。
5. **与直接全波 Top-K 在完整分支定界树上的端到端墙钟时间与剪枝比实测**：**未定**。
