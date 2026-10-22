> 技术附录采用独立章节、公式与文献编号；主报告为最终研究判断与实验配额的权威版本。本附录提供扩展推导，不报告新的 Maxwell 实验。

# Research P3-L：几何、测度与材料风险独立审查稿

本稿为主报告的数学补充材料。日期：2026-10-19。没有运行新的 Maxwell 求解、训练或付费实验。

## 0. 阅读范围与新旧贡献边界

已完整阅读 P1 `MODEL_AND_PRIOR.md`、P1 Zero-Physics Audit `RESEARCH_REPORT_ZH.md`、P2 `THEORY_REFINEMENTS.md`，以及历史 Research L 的 `RESEARCH_L_FULL_REPORT_ZH.md` 全部数学/概率部分和完整 `THEORY_APPENDIX_ZH.md`。主报告负责四份项目完整证据及其他历史包。本稿不把这些已读文档的结论伪称为本轮新发现。

Research L 已经完成或提出：五类局部集合、先验内常秩 fiber、原先验横截消除 ambient nullspace 的反例、finite Maxwell resolvent 恒等式及二阶余项、decoder 与 Maxwell 两部分曲率、A17 inverse-precision 缺陷闭合、Jacobian 误差及 Davis–Kahan 扰动、coarea、normal Laplace 厚度、局部 posterior 不能给全局区域质量、pCN/LIS 和 manifold lifting。**这些都属于历史成果或成熟理论，P3-L 不能再包装成新定理。**

本轮相对历史 L 的实质推进是：

1. 针对真实 P1 quintic saturation 的 C^{2,1}、uniform property 边界与混合维度材料推前，建立适用的条件 ridge 正则性。
2. 把“弱方向探索”改为“强条件分布积分 + 弱边缘积分 + 材料条件均值”，明确哪些降维能保留 Bayes estimator、哪些会造成 singular posterior approximation。
3. 给出支持截断、物理近似、normal volume、条件材料均值及 Monte Carlo 的统一材料风险误差分解。
4. 证明固定少量零宽轨迹不能承担全维 posterior density，并给出可计算的 tube 质量与覆盖尺度。
5. 用 P1 有限网格、分段多项式 decoder 和有理 Maxwell forward 推出半代数分层/有限拓扑类型的特化结论。
6. 将曲率与谱旋转界落实成可验证的 geometry-refresh 条件，同时保留其可能过保守、没有实际加速的失败结论。

下面每项均标明“标准结论”“本轮条件推导”或“待实验假设”。标准定理在 P1 中的正确组合可以成为方法基础，但独立 TAP 贡献仍需要可计算界、正确概率及净成本收益。

## 1. 固定统计合同、坐标与材料度量

P1 原始坐标为 q=(z,u)∈R^{64}×[0,1]^2，z∼N(0,I)，两个 u 独立均匀。相对于 dq=dz du：

\[
p_0(q)=(2\pi)^{-32}\exp(-\|z\|^2/2)\,1_{[0,1]^2}(u).
\]

用已固定逐源实虚排序及白化，定义 f(q)=\mathcal W\mathcal R F(g(q))∈R^{512}，y_w=\mathcal W\mathcal R y，r=f-y_w。P1 的复噪声 CN(0,σ²) 对应实虚方差 σ²/2，因此

\[
\Phi(q)=\tfrac12\|r(q)\|^2=\|F(g(q))-y\|^2/\sigma^2,
\qquad U(q)=\Phi(q)+\tfrac12\|z\|^2.
\]

原后验是 Π(dq|y)=Z^{-1}e^{-U(q)}1_{[0,1]^2}(u)dq。Noise-ball 以及 ΔΦ gate 均不是这个 posterior 的硬支持，除非明确声明另一次条件化。

材料统计量可实化为 X(q)=W_χ^{1/2}\mathcal R_χg(q)。取 W_χ=(1/576)I 时，欧氏范数即全网格复材料 RMS。材料局部度量

\[
M(q)=Dg(q)^T W_\chi Dg(q)
\]

只是任务几何；它不会改变 latent prior，也不允许在原 posterior 中凭空乘或除 Gram determinant。只有把同一个推前测度改写为材料 Hausdorff/体积密度时才出现相应换测度项。

P1 的全局材料推前不是一个固定维度的平滑密度：所有 ROI cell 均落在 pure-host plateau 的 latent preimage 可含开集，因而纯宿主材料可以拥有正 prior mass；条件似然有限时也保留非零 posterior atom。审计 anchor 的 rank 66 只说明该局部规则片，没有消除全局平台、材料原子或低秩片。所有“点/曲线概率为零”的结论必须限定在相应绝对连续规则分量。

均匀 u 的内部负 log-prior Hessian 为零，不能把 12I 自动称为其 prior precision。12I 可以作为 Var(u)^{-1} 的尺度度量，但真实 log-prior 仍是 bounded uniform。将 u=Φ_N(s) 高斯化是合法先验 transport，但边界 MAP u=0 或 1 对应无穷 s，不能直接将原边界锚点当成有限高斯内点。

## 2. 定理 G1：满列秩意味着什么、不意味着什么

**性质：标准常秩/逆函数理论的直接应用；历史 L 已有 exact-fiber 版本。本轮补充有限噪声拓扑与反例。**

设 q_0 是 prior 内点，f:Q⊂R^n→R^m 为 C^1，n≤m，Df(q_0) 满列秩 n。

1. 存在邻域 N，使 f|_N 为嵌入；因此 f^{-1}(f(q_0))∩N={q_0}。
2. 若 Φ(q_0)<c，则 {q∈N:Φ(q)<c} 包含 q_0 的 n 维开球。它是全维统计区域，不是 n−rank(Df) 维的 fiber。
3. 若 ∇Φ(q)≠0，则 {Φ=c} 在该点附近是 n−1 维规则边界；先验边界会产生角、面及交叉 strata。
4. 满秩本身不保证 Φ 或 U 的 Hessian 正定，不保证 posterior unimodality、凸性、单连通或“单一解族”。

证明 1：从 Df 的 n 个独立行选出线性观测投影 P，使 D(Pf)(q_0) 可逆。逆函数定理给 Pf 的局部逆，因此 f 局部单射。证明 2 来自 Φ 连续性。证明 3 是 scalar level-set theorem。

若进一步 f∈C²、f(q_*)=y_w，则 ∇²Φ(q_*)=J_*^TJ_*≻0。由 Hessian 连续性，可取足够小的凸邻域使 Φ 严格凸，足够低的局部 sublevel 是凸的、可缩的，并与 n 维闭球同胚。这是附加的局部强条件，不是一般带噪拟合的结论。

一般情况

\[
\nabla^2U=J^TJ+P+\sum_i r_i\nabla^2f_i,
\quad P=\operatorname{diag}(I_{64},0_2)
\]

在原 uniform 内部成立。可用充分条件

\[
\lambda_{\min}(J^TJ+P)>\left\|\sum_i r_i\nabla^2f_i\right\|
\]

证明局部正曲率。高维白化 residual 的范数通常不是“小于 1”，因而不能不检验最后一项。

### 满秩却出现环形 likelihood 的反例

取 f(x_1,x_2)=(x_1,x_2,x_1²+x_2²)，y=(0,0,1)。其 Jacobian 对所有 x 都满列秩 2，且 f 全局单射，但令 s=x_1²+x_2²，

\[
\Phi=\tfrac12[s+(s-1)^2]=\tfrac12(s²-s+1).
\]

最小值出现在 s=1/2，较低 sublevel 是环域。数据位于 forward image 之外造成了这一几何；并不存在任何正维 exact fiber。加入足够弱的 Gaussian prior 或把第三个观测值增大，仍可保留环形 posterior。

### P1 实际含义

已审计 anchor 的 Dg 和 Df 数值秩均为 66，意味着不存在已证实的正维 local exact fiber。数值满秩仍需误差与阈值条件，不能升级为所有 q 的数学证明。P2 的有限噪声材料变化与上述事实完全相容。

## 3. 三种谱必须分开

**性质：成熟 LIS、任务导向 Bayes 降维的区分；本轮建议将其作为 P1 实验的首要消融。**

### 3.1 相对先验的信息谱

高斯参考内，H_L=J^TJ，prior precision P_0≻0，LIS 解 H_Lv=λP_0v。在线性 Gaussian 情形中，λ衡量该方向从 prior 到 posterior 的方差收缩。P1 的 bounded uniform 不能不加说明地纳入同一个 Gaussian precision 模型。[R1,R2]

### 3.2 单位材料变化的观测谱

在 M≻0 的规则点，解 H_Lv=λMv。这衡量同样材料 RMS 变化的可见性。P1 audit 中 material-normalized condition≈1.91×10^4，与 raw latent-J condition≈1.86×10^8 相差很大，恰说明 decoder 缩放不能当作物理弱可见性。这个谱仍未编码某方向达到任务位移需要支付多少原 Gaussian prior 代价。

### 3.3 后验材料不确定性谱

在有根据的局部 Gaussian Π_G=N(m,Σ) 下，令 K=Σ^{1/2}MΣ^{1/2}。若选 d 个线性坐标 a=V^TΣ^{-1/2}(q-m)，V^TV=I，则

\[
E[q\mid a]=m+\Sigma^{1/2}Va,
\]

未保留的材料条件方差为

\[
\mathcal E_d(V)=\operatorname{tr}[(I-VV^T)K].
\]

由标准 Ky Fan/PCA 原理，取 K 最大 d 个特征方向使这个量最小，最小值为其尾特征值之和。等价广义谱为 Mv=τΣ^{-1}v。它关注真正残留的材料不确定性；并不等于选 raw J 最小奇异向量，也不等于选最大 prior-to-posterior 收缩方向。任务导向 Bayesian 近似已有成熟先例，不能单凭此谱宣称原创。[R3]

**核心反例：**线性模型 y=εq+η、q∼N(0,I_n)。所有 likelihood eigenvalue ε² 很小，但 posterior 仍近似整个 n 维 prior。把 n−d 个方向冻结成一个点，得到的 measure 与原 posterior 相互奇异，TV=1。弱信息代表数据不能压缩这些方向，并不代表可以删除这些方向。

如果采用 LIS 的正确补 prior 近似，而非冻结补空间，令保留强方向、忽略余下 λ_i 对 likelihood 的影响。在线性 Gaussian prior-predictive 平均下，其 KL 误差为

\[
E_Y\operatorname{KL}(\Pi(\cdot\mid Y)\|\widetilde\Pi(\cdot\mid Y))
=\tfrac12\sum_{i\in\mathrm{omitted}}\log(1+\lambda_i).
\]

推导：各坐标 posterior variance s_i=(1+λ_i)^{-1}，posterior mean 的 prior-predictive 二阶矩为 1−s_i；将二者代入 KL(N(m_i,s_i)||N(0,1)) 即得。可见大量“小但非零”方向的累计作用仍可能大。非线性问题仅在一个 anchor 的谱小还不能给此结论，需要 posterior 范围内的梯度/曲率控制。Zahm 等的 certified reduction 基于积分的 log-likelihood 梯度二阶矩及 log-Sobolev 条件，而非一个点的最小奇异值。[R2]

## 4. 定理 G2：强条件 ridge 的正则性、尺度与材料曲率

**性质：IFT 与 Lipschitz 导数的条件推导；P1 正则性及边界是本轮新增适用性审查。**

先用固定可逆 V=[V_w,V_s] 建立 affine coordinates q=q_c+V_wa+V_sb。把常数 |detV| 吸收进归一化，U(a,b) 仍相对于 da db 定义。假定强条件 optimum 在 interior，

\[
U_b(a_0,b_0)=0,\qquad U_{bb}(a_0,b_0)\succ0.
\]

如果 U∈C²，则唯一局部 stationary graph b_*(a) 存在且为 C¹，

\[
Db_*=-U_{bb}^{-1}U_{ba},
\quad
\nabla^2U_{\rm prof}=U_{aa}-U_{ab}U_{bb}^{-1}U_{ba}.
\]

U_prof(a)=U(a,b_*(a)) 的 Hessian 已在 C² 假设下有意义；但 D²b_* 和连续曲率一般没有。

若在一个 graph 邻域内有 U_bb≽μI、||U_ba||≤γ、Hessian ∇²U 的 Lipschitz 常数≤L_H，则 b_*∈C^{1,1}，||Db_*||≤k=γ/μ，且几乎处处

\[
\|D^2b_*\|\le\frac{L_H}{\mu}(1+k^2).
\]

证明：对 U_b(a,b_*(a))=0 一次求导得到上式。Lipschitz Hessian 与有界 inverse 使 Db_* 局部 Lipschitz。令 Q(a)=(a,b_*(a))，在几乎处处可微的 a 上，对 U_bq(Q(a))Q'(a)=0 求导：

\[
U_{bb}D^2b_*[u,v]
=-D_a[U_{bq}\circ Q][v]Q'u,\quad Q'(a)u=(u,Db_*u).
\]

这里取的是 Hessian 沿 graph restriction 的导数，不能假定 ambient U''' 恰好在 graph 上存在。restriction 局部 Lipschitz，且 ||D_a[U_bq∘Q]||≤L_H||Q'||，由 ||Q'||≤sqrt(1+k²) 得界。若 U∈C³，右侧成为熟悉的 −D²(U_b)[Q'v,Q'u]，逐点成立并给连续二阶导；若仅 C^{2,1}，必须使用几乎处处或分片版本。

### 可验证的存在半径

设 U_bb≽μI 在 A_R×B_r(b_0) 上成立，U_b(a_0,b_0)=0，且 ||U_b(a,b_0)||≤γ||a−a_0||。若 γR<μr，则每个 a∈A_R 的条件最小点都在 B_r 内部，唯一且

\[
\|b_*(a)-b_0\|\le\frac\gamma\mu\|a-a_0\|.
\]

证明：在球边界 b=b_0+rv，v^TU_b(a,b)≥μr−||U_b(a,b_0)||>0，因此最小点不在边界；强凸性保证唯一。再对两个梯度作强单调性估计得到距离界。这个结果需要对整个产品邻域的 μ 和 γ 有效界，不能用单点 Hessian 代替。

### P1 不是一般 C²，而是局部 C^{2,1}

P1 饱和函数：s≤−1 时 h=0，s≥1 时 h=1，−1<s<1 时 t=(s+1)/2，h=6t^5−15t^4+10t³。

\[
h'(s)=15t²(1-t)²,
\quad h''(s)=15t(1-t)(1-2t).
\]

h'' 在接缝连续、全局 Lipschitz；h''' 的内部单侧极限在两端均为 7.5，外部为 0。准确的全局标量界为 ||h'||∞=15/16、||h''||∞=5sqrt(3)/6、Lip(h'')=7.5。因此 h∈C^{2,1}，不是 C³。固定非奇异 Maxwell 紧邻域内，g、f、U 也局部 C^{2,1}。条件 MAP ridge 可以拥有有界几乎处处曲率；固定每个 cell 的 saturation 状态的开 strata 内更光滑。跨接缝不能假定连续 principal curvature 或任意三阶 Taylor。已知 decoder 界可在不调用 Maxwell 的情况下由固定 field matrix 推出。

### C² forward 不保证 C² MAP ridge 的反例

取 f(a,b)=(εa,b,c+|a+b|³)，ε,c>0，y=0。f∈C²，Df 处处满列秩，U=||f||²/2。条件方程为

\[
b+3(c+|a+b|³)(a+b)|a+b|=0.
\]

在原点附近 b_*(a)=−3c a|a|+o(a²)，一阶导连续，而二阶导两侧极限不同。这是 finite residual 乘第三阶 forward 正则性缺口的实际影响。

### 材料曲率不是 latent 曲率

设 X(a)=W_χ^{1/2}g(q(a))，q(a)=q_c+V_wa+V_sb_*(a)。在欧氏正交 V、||Dg||_W≤G_1、||D²g||_W≤G_2 下，几乎处处

\[
\|D²X\|\le(1+k²)\left(G_2+G_1 L_H/\mu\right).
\]

若 σ_min(DX)≥s_*>0，则材料第二基本形式的保守界为

\[
\|\mathrm{II}\|\le
\frac{(1+k²)(G_2+G_1L_H/\mu)}{s_*²}.
\]

当 decoder 接近平台使 s_* 很小，材料 metric 会退化，即使 latent graph 很平整，材料曲率与 chart 条件数也可能爆炸。曲率界只控制局部弯曲；reach 还受到远处自接近/自交、边界的限制，不能只写 reach≥1/κ 而不处理全局分离。

## 5. 更平滑的 strong-response spine 与三角 tube chart

**性质：标准 IFT/三角换元；本轮建议利用其数值中心自由度，避免把 MAP ridge 当唯一“自然流形”。**

固定强观测 frame C_s，定义

\[
h(a,b)=C_s^T[f(q_c+V_wa+V_sb)-f(q_c)].
\]

若 h_b=C_s^TJV_s 可逆，则 h=0 给 C²（P1 中 C^{2,1}）graph b_h(a)。若 C_s 取强方向 J V_s 的正交输出 frame、JV_s 与 JV_w 正交，则 b_h'(0)=0 且

\[
b_h''(0)[u,v]
=-(C_s^TJV_s)^{-1}C_s^TD²f[V_wu,V_wv].
\]

它保持选定 strong response，不保持全部数据，也不是 posterior density ridge。

对任何足够规则的 spine c(a)，用

\[
q=T_c(a,z)=q_c+V_wa+V_s[c(a)+z]
\]

建立全维 tube coordinates。其导数为

\[
DT_c=[V_w,V_s]
\begin{bmatrix}I&0\\Dc&I\end{bmatrix},
\quad |\det DT_c|=|\det[V_w,V_s]|.
\]

**因此，对 z 的全部合法范围正确积分后，π_w(a) 与 spine c(a) 的选择无关。** 选择平滑 projected-response spine、MAP ridge 或条件均值中心，仅影响 numerical conditioning 和 strong proposal 的效率。只保留 z=0 时才改变支持与统计目标。

如果强范围来自原 uniform box，积分域 B(a) 必须随 recentering 正确平移。强方向不能因为“用来 correction”就被从 prior measure 中删除。

在 material surface 上人为定义 density w.r.t. dvol_G 与上述 weak marginal 不同。如果 weak coordinate density 为 ρ(a)，其推前在 surface 上的密度是 ρ(a)/sqrt(detG(a))；不能把 sqrt(detG) 再乘入原 marginal。

## 6. 固定强方向、profile 与 marginal：精确比较

写原坐标 joint posterior 为 π(a,b)∝e^{-U(a,b)}1_{B(a)}(b)。

| 操作 | 实际得到的对象 | 丢失的内容 |
|---|---|---|
| b=b_0 固定 | 选定版本的条件 slice π(a|b_0)，或人工 slice density | 强条件均值随 a 的变化、强方向体积、强方差 |
| b=b_*(a) 条件最小化后用 e^{-U_*} | profile-potential law | normal volume、边界体积、非Gaussian强条件形状 |
| 对 b 正确积分 | π_w(a)∝Z_s(a)=∫e^{-U(a,b)}db | 不丢失 weak marginal；材料估计仍需 E[g|a] |

条件在一个连续点 b_0 本身也是零概率事件；需要指定 regular conditional density 或 slice 极限，不能说该 slice 拥有原 posterior 的非零质量。

### 强积分的精确 Laplace 校正分解

设 b_* 为 interior mode，H_s=U_bb(a,b_*)≻0，R_a(z)=U(a,b_*+z)−U_*−z^TH_sz/2。则完全精确地

\[
Z_s(a)=e^{-U_*}(2\pi)^{s/2}|H_s|^{-1/2}\,C(a),
\]

\[
C(a)=E_{\xi\sim N(0,I_s)}\!left[
e^{-R_a(H_s^{-1/2}\xi)}
1_{\{b_*+H_s^{-1/2}\xi\in B(a)\}}
\right].
\]

普通 Laplace 相当于 C(a)≈1。真正导致 weak marginal bias 的是 log C(a) 随 a 的变化；若它是同一个常数，归一化后会抵消。

如果无边界且所有 strong b 上

\[
(1-δ)H_s\preceq U_{bb}(a,b)\preceq(1+δ)H_s,
\quad0<δ<1,
\]

则由二次上下界积分得

\[
(1+δ)^{-s/2}\le C(a)\le(1-δ)^{-s/2}.
\]

这显示误差随 strong 维数 s 累积；“每个方向都强”不等于高维 Laplace normalization 准。边界、多峰 strong 条件、nonquadratic tails 都可使此假设失效。

### P1 property active boundary 的解析反例

令 inward strong coordinate b≥0，U(a,b)=u(a)+c(a)b+h(a)b²/2，c≥0。条件 optimum 在 b_*=0，但

\[
Z_s(a)=e^{-u(a)}\sqrt{2\pi/h(a)}
\exp\!\left(\frac{c(a)^2}{2h(a)}\right)
\Phi_N\!\left(-\frac{c(a)}{\sqrt{h(a)}}\right).
\]

只有 c=0 时才是半个 Gaussian；c≫sqrt h 时宽度约为 1/c。直接加 −(1/2)logdetH 仍会算错 volume。P1 审计 anchor A 处在损耗上界，故该问题不是理论上的罕见情况。应显式做 truncated conditional integration，或选择有代表性的 interior posterior states；不能在 active MAP face 上赋予虚构的 surface mass。

## 7. 定理 G3：条件材料均值、降维与平方风险

**性质：标准条件期望/方差分解；在 P3-L 中是研究目标转向的核心，而非新概率定理。**

在声明局部区域 C 内，令 X=g(q)，A=a(q)，

\[
m(a)=E[X\mid A=a,y,C],
\quad C_s(a)=\operatorname{Cov}(X\mid A=a,y,C).
\]

则

\[
m_C=E[X\mid y,C]=E_{\pi_w}[m(A)],
\]

\[
\operatorname{Cov}(X\mid y,C)
=E[C_s(A)]+\operatorname{Cov}(m(A)).
\]

由此有三个直接结论：

1. 精确 strong marginalization 与 conditional material mean 可以完全保留 full local MMSE，即使 strong 方差非零。
2. 只平均 g(a,b_*(a)) 一般不等于 MMSE；decoder 曲率及 strong mean-mode 差会产生系统偏差。
3. 若输出只有 m(a) 的样本，即使最终均值正确，仍会丢失 E[C_s(A)]，其 credible intervals 会过窄。

若条件 strong distribution 对 b 强 log-concave、U_bb≽μ(a)I、无边界且具有足够尾衰减，则

\[
E[\|b-b_*\|²\mid a]\le s/\mu(a).
\]

证明：积分分部给 E[(b−b_*)^TU_b]=s，强单调性给 integrand≥μ||b−b_*||²。若 g 对 b 的材料 Lipschitz 常数为 L_s(a)，ridge collapse 的均值误差满足

\[
\|m(a)-g(a,b_*)\|_W\le L_s(a)\sqrt{s/\mu(a)}.
\]

这只是保守上界。若强条件关于 b_* 对称，decoder 有界二阶导且条件 mean=b_*，则可改为二阶量：

\[
\|m(a)-g(a,b_*)\|_W
\le\tfrac12\sup_b\|D^2_{bb}g\|_W\,E\|b-b_*\|².
\]

若不对称，则还要加 D_bg·(E[b]−b_*)。不能仅以“强特征值大”忽略这一项。

更精细的 conditional material variance 可由条件 Poincaré/Brascamp–Lieb 型界估计：在相应强 log-concavity 假设下，

\[
\operatorname{tr}W C_s(a)
\le E\left[\|W^{1/2}D_bg\,K(a)^{-1/2}\|_F²\mid a\right]
\]

当 U_bb≽K(a)≻0。此处的正则性、域边界及 tail 假设必须单列；这属于成熟概率不等式的应用，不是 Maxwell 新定理。[R13]

### 材料均值可以不在原 decoder image

P1 每个样本共享一个全局 a+ib，但 E[φ_i(a+ib)] 的实虚比可随像素 i 改变。故 Eχ 可能不再等于某个 g(q)，即使仍满足点wise passivity。无约束平方损失下它仍为 Bayes 最优。

若必须输出原 material law 中的单个合法材料，最优决策为

\[
\hat\chi_{\mathrm{legal}}
\in\arg\min_{\chi\in g(Q)}\|\chi-m_C\|_W²,
\]

而非 MAP。其相对于无约束 MMSE 的额外条件 Bayes risk 恰为 dist_W(m_C,g(Q))²。

若 complex MAE=Σ_i|\hatχ_i−χ_i|，其 Bayes decision 是每像素复平面的 geometric median。仅当损失定义为 Σ(|ΔRe|+|ΔIm|) 时才是分别的 real/imag median。

## 8. 定理 G4：支持、物理、体积、条件均值与 Monte Carlo 的统一风险界

**性质：本轮明确推导的组合界；各工具标准。它把 P3-L 的可实现误差预算直接连接到材料风险，具有优先验证价值。**

设 g(Q) 在材料 W 范数下有有限直径 D。P1 bounded passive material law 满足这一条件。记 full posterior Π，Π_C=Π(·|C)，ε_C=Π(C^c)。以近似 potential \widetildeΦ 得到同 prior、同 C 的 Π_C^F；其 weak marginal 为 ρ_F、条件材料均值为 m_F(a)。算法构造 weak density \widehatρ 和 conditional mean approximation \widehat m(a)，最终总体均值

\[
\widetilde m=\int\widehat m(a)\widehat\rho(a)da.
\]

假设这些均值输出属于 conv(g(Q))。定义

\[
e_{\rm cond}=\int\|\widehat m(a)-m_F(a)\|_W\widehat\rho(a)da,
\]

\[
e_{\rm marg}=\operatorname{TV}(\widehat\rho,\rho_F),
\quad\omega_F=\operatorname{osc}_{C}(\widetilde\Phi-\Phi).
\]

则

\[
\boxed{
\|\widetilde m-m_{\rm full}\|_W
\le D\epsilon_C+D\tanh(\omega_F/4)+D e_{\rm marg}+e_{\rm cond}.
}
\]

若 weak unnormalized density 的 log 近似误差振幅为 ω_V，则 e_marg≤tanh(ω_V/4)。该项可以来自 Laplace volume、learned density 或 trajectory/tube quadrature；应按实际误差来源分别报告，不能仅给总近似分数。

### 证明

第一步，m_full=(1−ε_C)m_C+ε_Cm_out，故 ||m_full−m_C||≤Dε_C。

第二步，对任意材料取值直径 D 的两个概率分布，最大耦合给 ||E_μg−E_νg||≤D TV(μ,ν)。

第三步，若 ν∝h μ，h∈[m,M]，令 c=E_μh，则

\[
\operatorname{TV}(\mu,\nu)=\frac{E|h-c|}{2c}
\le\frac{(M-c)(c-m)}{c(M-m)}.
\]

右侧在 c=sqrt(mM) 最大，等于 (sqrt M−sqrt m)/(sqrt M+sqrt m)=tanh(log(M/m)/4)。取 h=e^{-(\widetildeΦ−Φ)} 即得 physics TV bound。特别地，若 sup_C|\widetildeΦ−Φ|≤δ，则 ω_F≤2δ，故 TV≤tanh(δ/2)。

第四步，插入 ∫m_F\widehatρ，再用三角不等式与上面的有界材料耦合，即得 e_cond+D e_marg。

### 连接平方 Bayes risk

对正确 full posterior，任意材料估计 t 的条件平方风险满足精确恒等式

\[
R_y(t)-R_y(m_{\rm full})=\|t-m_{\rm full}\|_W².
\]

因此上述右侧平方直接界定 deterministic approximation 的 Bayes excess risk。若 Monte Carlo 估计 \widehat m_N 对 \widetilde m 无偏，则

\[
E_{\rm MC}[R_y(\widehat m_N)-R_y(m_{\rm full})]
=\|\widetilde m-m_{\rm full}\|_W²+
\operatorname{tr}(W\operatorname{Var}_{\rm MC}\widehat m_N).
\]

有限链 initialization 或 approximate stationary bias 另加到 deterministic error，不能靠 ESS 隐去。对相关链，其实际 covariance 是 lag-covariance sum；“一个 latent ESS”不能统一替代所有材料统计量的 ESS。

### Maxwell output error 如何进入这一界

若 whitened forward 的误差 ||f−\widetilde f||≤ε_f，则

\[
|\Phi-\widetilde\Phi|
\le\|\widetilde f-y_w\|\epsilon_f+\epsilon_f²/2.
\]

这是已有 L/P2 finite-response bound 与 posterior risk 的直接接口。任何只在几个端点拟合的 ε_f 都不是 C 上的 uniform bound；不能据此宣称 certificate。

### 解释及可检验预测

- 当 escape mass 主导时，延长局部链或扩大 Flow 无法修复支持偏差。
- 当 e_cond 主导时，增加 weak sample 数无效；应改进 strong conditional integration 或 decoder conditional mean。
- 当 e_marg 主导时，candidate novelty 与 Best-of-K 均不能说明问题；应做 normal volume 与 density 校正。
- 当 physics 误差已经远低于 MC 方差时，继续提高 ROM 精度未必改善材料风险；此时提升 ESS 更有价值。
- 同等预算下，以材料 risk/error 分量选 refresh 次数，比统一固定误差阈值更可能有增量。这是待实验假设，不是已证明性能优势。

## 9. GT coverage 与局部模型失配

GT 差量大量落在 anchor 的弱材料切空间，只是离线局部表达能力诊断，既不证明有限可达，也不证明这些方向拥有足够 prior mass，更不证明局部区域 C 拥有足够 posterior mass。

对正确模型的 prior-predictive 重复实验，若 A_y⊂C_y 且 Π(A_y|y,C_y)=1−α，则

\[
P(q^\dagger\in A_Y)
=E_Y[(1-\epsilon_C(Y))(1-\alpha)].
\]

若局部 escape mass 非零，局部“95% credible”没有全局 95% coverage。这个等式属于 Bayes prior-predictive coverage；并不保证每个固定 GT 在重复噪声下的 frequentist coverage。

仅观察被限制在 C 内的链，无法估计 ε_C。Posterior predictive checks 也可能无法发现材料弱方向的漏质量，因为错误材料可以产生相似数据。需要独立 full-target reference、SMC/tempering、允许逃逸的 control kernel，或明确标 UNKNOWN。

GT 到 surface S 的最短距离只下界“必须输出 S 中一个点”的误差；它不能下界 ensemble mean 的误差。均值属于 conv(g(S))，故正确的几何障碍是到这个凸包的距离。例：单位圆上所有点距离原点 1，均匀平均却正好是原点。报告必须同时保留 surface distance、conditional mean bias、strong variance 与 escape mass。

## 10. 定理 G5：固定轨迹、tube 与 surface 概率的覆盖障碍

**性质：面积公式与球覆盖的直接推导；本轮用于排除错误的 trajectory probability 解释。**

设 weak marginal 在 A⊂R^d（d>1）具有连续密度 ρ，固定有限 K 条 C¹ 或有限长度曲线 γ_k 的并 Γ 的 d 维 Lebesgue measure 为零。因此以这些固定曲线为全部生成支持的 law 与 ρ(a)da 相互奇异，TV=1。不能把所有采样点的 softmax 称为 surface posterior。

限定非常重要：随机初始点、随机轨迹族或全维噪声可能使单个样本的边缘 law 为全维。普通 MCMC 的有限经验 measure 也离散、与连续 target 的 TV 为 1，但它可以在统计量/弱收敛意义上一致。这里否定的是**固定低维生成支持**，不是轨迹形式的 MCMC。

若 trajectory generator 为 a=T(k,t,z)，z∈R^r，T 可微，则其单样本 law 需要有效 rank≥d 才可能有 d 维密度；仅一个 t 加 discrete k 不够。随机曲线的额外随机性本身并不会免除有效维数要求。

### 有限宽 tube 的质量上界

假设 ρ≤B，所有曲线总长度为 L，定义 Γ_ε={a:dist(a,Γ)≤ε}。沿每条曲线以间距≤ε取点，球覆盖给

\[
\operatorname{Vol}_d(\Gamma_\epsilon)
\le 2^d\omega_d(L\epsilon^{d-1}+K\epsilon^d),
\]

\[
\Pi_w(\Gamma_\epsilon)
\le B2^d\omega_d(L\epsilon^{d-1}+K\epsilon^d).
\]

证明：每条长度 L_k 曲线可被至多 L_k/ε+1 个半径 2ε 的球覆盖其 ε-neighborhood，求和即可。常数保守，但无需预设曲线不自交。要覆盖质量 1−δ，宽度很小而 d 大时，总长度必须相应增长。少量 trajectory 无法在一个真正宽的 16/24 维 posterior 内免费绕过维数问题。

### 如何给单条轨迹合法权重

若人为定义 surface arclength 条件 density ρ_S，轨迹 q(t) 的坐标密度为 ρ_S(q(t))||q'(t)||_G。改变 t 参数化必须同步变换密度。它仍是指定的曲线条件 law，不是 full posterior mass。

若要代表原后验一条轨迹邻域，则必须给出 normal tube chart T(t,z)、横截面范围、Jacobian 及积分

\[
w(t)\propto\int e^{-U(T(t,z))}|\det DT(t,z)|dz.
\]

多条 tube 重叠时，可用 partition of unity ψ_k≥0、Σ_kψ_k=1 在覆盖域上拆分目标质量。每个 chart 积分包含 ψ_k，不能按各轨迹访问次数分配权重。若用 mixture importance proposal r=Σ_kα_kr_k，应在 denominator 中使用完整 mixture density，不能仅除以生成该点的那条 r_k。

### 少量轨迹何时仍可能有用

1. 主要材料统计量的 conditional variation 很低，目标是积分少数均值/风险函数，而不是恢复全密度。
2. posterior 已经在测度意义上集中于窄 tube，并有独立 ε_C 证明或数值估计。
3. 轨迹只是建立低维 atlas 的 skeleton，实际 sampler 仍保留足够 normal randomness。
4. 用轨迹学习 preconditioner、方向场、refresh 或计算预算，不把 visitation 当作 posterior。

有限正权 quadrature 可以匹配有限个 moment，并不等价于估计 posterior density。其权重必须由合法积分/importance 推断，不能通过 held-out GT 倒推。

## 11. 定理 G6：P1 的半代数分层与有限拓扑类型

**性质：本轮针对真实 P1 的具体推导；结论使用成熟的半代数、Whitney 与 Hardt 定理。它不是 Maxwell 通用拓扑定理。**

假设：固定 finite grid、ROI、DCT field basis、Green matrix G、receivers S、incident fields、y 和噪声；保留原 q=(z,u)；在局部半代数区域 Q_R 内 L(q) 非奇异；local constraints 使用 polynomial/rational inequalities，例如 ||q−q_0||≤R、材料 RMS 半径、uniform bounds 或 det 下界。

**结论：**以下集合均为半代数：

- finite-noise likelihood sublevels {Φ≤c}∩Q_R；
- 原 dzdu reference 下的 latent posterior HPD region {U≤c}∩Q_R；
- Dg、Df 的各固定 rank 集合；
- 固定 frame 的 strong-response constraints h=0；
- interior strong stationarity U_b=0，以及原 uniform box 上的 KKT 集合；
- 通过在半代数 strong 域内最小化 U 所定义的 argmin graph（允许多值）。

它们具有有限个连通分量、有限 Betti numbers，并可取有限兼容 Whitney stratification。在紧致的 Q_R 上，sublevel family 随阈值 c 只有有限个半代数拓扑类型区间；跨临界阈值/边界可能发生变化。普通 Morse 描述需要额外非退化性；泛型 stratified Morse 不能自动视作 P1 已满足。

### 证明

1. 每个 cell 的 Gaussian field f_i(z) 对 z 线性；saturation 界 f_i=±0.5 是有限 hyperplane arrangement，给出有限 sign cells。
2. 每个 sign cell 内，φ_i 是零、壹或五次多项式；复材料 χ_i=χ_b+φ_i(a(u)+ib(u)) 是实变量 q 的多项式（总次数最多六）。
3. 离散 current Maxwell system L=I−diagχ G 的实化矩阵条目为实多项式。在 detL≠0 区域，L^{-1}=adj(L)/detL，forward 的实虚部为有理函数。也可使用 2N×2N real block matrix，其 determinant 为 |det_complexL|²。
4. Φ 是有理函数，U=Φ+||z||²/2 仍有理；固定边界跨有限 sign cells 取并得到半代数函数/集合。rank 条件由 minors 等式与不等式表达。导数在每个 cell 内有理，并在 C² 接缝相容。
5. argmin graph 可写为量词公式 U(a,b)≤U(a,b') 对所有合法 b' 成立；Tarski–Seidenberg 消元保证其半代数性。局部 argmin 也可用存在半径的量词表达。
6. 半代数集的分层/三角剖分给有限拓扑性质。对 E={(q,c):q∈Q_R,U(q)≤c} 到 c 的投影应用 Hardt local triviality，得到有限参数分割与区间内的拓扑平凡性。[R6–R9]

### 不能进一步声称什么

posterior density e^{-U}/Z 本身通常不是半代数；这里只用 HPD 与 U sublevel 等价。**积分 strong coordinates 后的 marginal density、material-volume HPD（含 log Gram）、以及 nonlinear prior transport 后的 HPD，不自动保留上述半代数结论。**

有限 strata 数可能巨大。256 个 ROI cell 产生多达 512 个 saturation hyperplanes，66 维有理约束的代数消元复杂度可能完全不可用。有限 Betti numbers 不等于低 Betti numbers，更不等于单连通。

本结论依赖 fixed discretization；未知几何导致 Green kernel 随连续位置变化、continuous Maxwell 的无限维极限、非多项式 learned decoder、未知 nuisance law，均需另证。fixed G 中的 Bessel/Hankel 数值只作为实常数，不破坏半代数性。

### 对算法的有限但明确价值

- saturation regime/event 可作为 atlas metadata，避免把接缝误判成 Maxwell 奇异性。
- 在一个已验证 regular stratum 中用光滑 predictor；到接缝切换 C^{1,1}/secant 法或刷新。
- 可以在 2D/4D slices 对阈值进行 continuation 与拓扑持久性检查；这些 slice 不能证明 66D 全局拓扑。
- 不值得在 66D 直接做完整符号分层或 Betti 枚举。该 lemma 提供理论边界与实验解释，不能独立保证算法加速。

## 12. 曲率、Maxwell 反馈与 geometry-refresh 半径

**性质：finite resolvent 与谱扰动在历史 L 已有。本节只给 P1 可实施的组合形式；其计算收益待验证。**

current 方程 j=X(e+Gj)，L=I−XG，R=L^{-1}。有限材料增量 ΔX 对应

\[
\Delta j=(I-R\Delta XG)^{-1}R\Delta X E.
\]

设 K=RΔXG、b=RΔX E，κ=||K||<1，则

\[
\|\Delta f-Sb\|\le\|S\|\frac\kappa{1-\kappa}\|b\|.
\]

对 latent step δq，先精确计算便宜的 Δχ=g(q+δq)−g(q)，用完整材料导数 J_χΔχ 作 predictor；其 Maxwell 余项与 decoder 余项必须分离：

\[
f(q+δq)-f(q)-J_qδq
=J_χ[\Deltaχ-Dg\,δq]+e_{\rm Maxwell}.
\]

P1 已读 audit 的一个 pair 中，decoder secant 是主要失真来源，故不能把全部有限步失真归因于高阶 Maxwell feedback。

二阶 current derivative

\[
j_{uv}=R\{D²X[u,v]E+DX[u]Gj_v+DX[v]Gj_u\}
\]

给出 ||D²f|| 的上界。若 ||R||≤B_R、||E||≤B_E、||DX||≤G_1、||D²X||≤G_2，则

\[
\|D²f\|
\le\|S\|B_RB_E\,[G_2+2B_R\|G\|G_1²].
\]

所有 whiten、实化、多源和 material norms 必须一致。这个上界可能很松；passivity 并不自动给出 raw Euclidean ||R|| 的实用小界。

### 12.1 预测方向何时不再能保证弱

若 ||D²f||≤L_f，则 ||J(q)−J(q_0)||≤L_f||q−q_0||。对固定单位 v、起始响应 α_0=||J(q_0)v||、弱阈值 α_*，距离

\[
s\le(\alpha_*-\alpha_0)/L_f
\]

保证该固定方向仍弱（右端须为正）。这只是安全半径下界，不能精确预言方向会在哪一点变强。若方向 v(t) 也旋转，需加入 ||J||·||v'||，其来源可能是 spine 曲率或谱平行输运。

### 12.2 moving material metric 的谱旋转

令 H=J^TJ、M=Dg^TW Dg，A=M^{-1/2}HM^{-1/2}，假设邻域 M≽mI。

\[
\|A(q)-A(q_0)\|
\le m^{-1}\|H(q)-H_0\|
+m^{-2}\|H_0\|\|M(q)-M_0\|.
\]

推导使用 ||M^{-1/2}−M_0^{-1/2}||≤||M−M_0||/(2m^{3/2})，可由 inverse-square-root 的 resolvent 积分证明。

若 J 的 Lipschitz 常数为 L_J、Dg 在材料范数下的 Lipschitz 常数为 L_D，s=||q−q_0||，则

\[
E_H(s)=(2\|J_0\|+L_Js)L_Js,
\quad E_M(s)=(2\|Dg_0\|_W+L_Ds)L_Ds.
\]

设所跟踪弱 cluster 与其补空间在 anchor 的 gap 为 γ，E_A(s)=m^{-1}E_H+m^{-2}||H_0||E_M<γ/2，则 Davis–Kahan 型界给

\[
\|P_w(q)-P_w(q_0)\|\le2E_A(s)/\gamma.
\]

可据此取 projector tolerance 下的最大复用半径。内部 eigenvalue crossing 不要求刷新单个 eigenvector 身份；cluster projector 才是稳定对象。gap→0 时该证书失效，应合并 cluster 或用 smooth spectral filter；这不是证明物理 geometry 发生奇点。[R5]

### 12.3 P1 的零 Maxwell saturation-event horizon

沿 affine latent ray z(t)=z_0+tv，s_i(t)=s_i(0)+t\dot s_i。下一次 saturation 接缝时间可由各正的 (±1−s_i(0))/\dot s_i 直接取最小值。这个计算不需要 Maxwell。

接缝时 g 的一、二阶导连续，所以不必每次强制 full refresh；但它终止了某个固定 polynomial stratum 的高阶证书。保留 C^{2,1} 全局界可以跨越接缝；若只用分片高阶 Taylor，则必须切换 regime。

### 成本与正确性

单点数值曲率估计不等于 uniform radius。若上界太松导致频繁刷新，该理论不能带来净加速。应比较：fixed interval、empirical secant trigger、可证曲率 trigger、概率决策误差 trigger，并统一计入 JVP、LU、adjoint、eigenproblem、fallback 成本。A17 只是构造或加速这些量的一个候选，普通 optimized Full、Krylov 和 goal-ROM 同样拥有复用权利。

### 12.4 动态弱谱分布未必可积：不能默认存在处处相切的 weak manifold

**性质：Frobenius 条件的本轮显式反例；不是新微分拓扑定理，却修正了原始研究设想中的一个重要缺口。**

在有 spectral gap 的区域，P_w(q) 的 range 给出一个秩 d 的 distribution E_q。要存在 d 维 surface，其每点 tangent 都等于 E_q，必须满足 involutivity：对所有取值于 E 的光滑局部向量场 v,w，Lie bracket [v,w] 也取值于 E。这是 Frobenius 可积性条件。[R12]

Jacobian 的 weak eigenvectors 一般不满足这个条件。取一个完全光滑、满列秩、没有 decoder gauge 的 forward

\[
f(x,y,z)=(\sqrt\varepsilon x,\sqrt\varepsilon y,\sqrt\varepsilon z,
z-xy,\sqrt\alpha y),\quad \varepsilon>0,\ 0<\alpha<1.
\]

其 GN metric 为

\[
H=\varepsilon I+\nabla(z-xy)\nabla(z-xy)^T+\alpha e_ye_y^T.
\]

原点 eigenvalues 为 ε、ε+α、ε+1。取最强单位 eigenvector n（其 z 分量为正），一阶扰动给

\[
n(x,y,z)=(-y,-x/(1-\alpha),1)+O(x²+y²).
\]

故

\[
n\cdot(\nabla\times n)(0)=-\alpha/(1-\alpha)\ne0.
\]

符号可直接核验：在原点，∂_xH=−e_ye_z^T−e_ze_y^T，∂_yH=−e_xe_z^T−e_ze_x^T。按简单最大 eigenvalue 的导数公式，∂_xn=−e_y/(1−α)、∂_yn=−e_x、∂_zn=0，故 curl 的 z 分量为 −1/(1−α)+1。最大 eigenvalue 与两个弱方向的 gap 分别为 1、1−α，εI 保证 Df 满列秩。

已做独立非 Maxwell 数值核验：由显式 5×3 J 构造 J^TJ、eigh 求最大 eigenvector、中心差分求 curl。α=0.2/0.6、ε=0.01、h=10^{-5} 时，分别得到 −0.249999999913281 和 −1.4999999983312495，与理论 −0.25/−1.5 一致；h 从 10^{-3} 缩至 10^{-5} 的误差按二阶下降。结果保存在 `p3l_frobenius_check.json`；新 Maxwell/LU/RHS 均为 0。用 [v_1,v_2]=Dv_2v_1−Dv_1v_2 的约定，normal bracket 为 +α/(1−α)e_z，与 one-form 外微分的负号一致。

三维中 plane distribution n^⊥ 的 Frobenius 条件正是 n·curl n=0，因此两个弱 eigen-directions 根本不能被拼成一张处处相切的二维 surface。数据/Jacobian 满秩、metric 光滑、谱隙非零都没有解决这一问题。

这不是 P1 Maxwell 的已测反例，而是对“满秩 + 有 gap 的弱谱 ⇒ 存在 weak manifold”这一一般推理的反例。P1 的具体弱分布是否近似可积，仍要用实际物理导数独立检验；不能把本代数例的 nonintegrability 数值移植成 P1 结论。

对一般 P_w，可测量 normal bracket leakage

\[
\mathcal B(v_i,v_j)=(I-P_w)
\{DP_w[v_i]v_j-DP_w[v_j]v_i\},\qquad v_i,v_j\in E_q.
\]

其为零是局部可积性的必要条件，在足够规则条件下对所有方向为零即 involutivity。不同 frame 的导数含 gauge 项，但 normal bracket 本身可不依赖 frame 标签地解释。

**算法含义：**

- Dynamic weak-frame trajectories 可以存在，路径次序甚至会产生二阶 normal 位移；这不构成某张 weak surface 的存在证明。
- 一个由 IFT 构造的 projected-response surface 的 tangent 是可积的，但不会一般等于所有位置重新计算的最弱 eigenspace。
- full-dimensional tube sampler 或正定 anisotropic diffusion 不需要 weak distribution 可积；这比把 P_w 直接解释为流形 Laplace–Beltrami operator 更稳健。
- 若要声称 genuine weak manifold，Stage I 应同时检查 curvature 和 normal bracket leakage；仅检查 spectrum rotation 不够。

### 12.5 可用的几何对象及正则性范围

对于真实 immersion X(a)，诱导 metric、tangent projector、normal projector 和第二基本形式分别为

\[
G=X_a^TX_a,\quad P_T=X_aG^{-1}X_a^T,
\quad P_N=I-P_T,\quad \mathrm{II}_{ij}=P_NX_{ij}.
\]

高 codimension 下 principal curvature 需要指定 normal direction n，再研究 G^{-1}(n^TX_{ij})；不能给出一个未声明 normal 的“主曲率列表”。

若 X∈C^{2,1} 且 G 一致正定，则 G∈C^{1,1}，Christoffel symbols

\[
\Gamma^k_{ij}=G^{k\ell}\langle X_\ell,X_{ij}\rangle
\]

局部 Lipschitz，geodesic ODE 可按通常存在唯一理论使用。P1 的 projected-response spine 可在规则区满足这个层次。仅 C^{1,1} 的 MAP spine 则只给 Lipschitz metric；其 connection 为几乎处处对象，不能未经额外论证套用同一 classical geodesic 唯一性结论。

对有 gap 的弱子空间 projector P(t)，可沿曲线通过

\[
\dot Q=[\dot P,P]Q
\]

运输 orthonormal frame；生成元 skew-symmetric 保持正交性。它避免逐 eigenvector 的 sign/permutation 跳变。闭合回路可能存在 frame holonomy，但这不是 posterior mode、材料形状拓扑或解集不连通的证据。需要将 frame gauge、distribution 可积性、intrinsic surface curvature 三者分开。

## 13. 论文可检验的核心预测与停止条件

1. **概率目标预测：**同一 weak frame 下，marginalized strong response 比 profile-only 更准确的程度，将与 log normal volume 的 variation、active boundary 和条件 skewness 相联系；与候选去重数未必相关。
2. **材料风险预测：**MMSE 的提升由 conditional mean 与正确 weak mass 决定，不由 Best-of-K 决定。线性近Gaussian对象上简单 blocked sampler/Rao–Blackwell 可能已足够。
3. **decoder 预测：**exact-decoder/material-linear surrogate 应在 P1 某些长 latent steps 明显优于 raw latent Taylor；是否转化为 posterior ESS 净收益需实测。
4. **geometry 预测：**refresh 可被拖延的对象应有较大 spectral gap、较小 metric variation、较低响应曲率及较远 saturation event；单看 weak eigenvalue 大小不能预测复用效率。
5. **轨迹预测：**固定 K 的窄 trajectory family 随 weak dimension 从 2 增至 8/16 时会迅速漏质量，除非 strong tube thickness 或随机 chart 维数同步增长。
6. **局部限制预测：**若 local escape mass 或 conditional normal material risk 大，增加局部 NN 大小不会解决校准与风险；应扩大 region 或保留 strong随机性。

最有论文价值的组合是：**原 prior 不变的 full tube target + 可计算的 conditional-mean/volume/physics risk budget + geometry-aware reuse**。其理论基础中的 IFT、coarea、LIS、PCA、Laplace、DK、TV 重加权界均不能单独宣称原创。真正需要新实验回答的是：P1 的材料/PDE 结构是否让这些误差分量可以低成本估计，且由此获得可重复的材料风险与 ESS/成本优势。

## 14. 已浏览的原始研究与出处

以下来源只用于支持对应成熟理论，不用于冒称本轮经验结果。主报告如使用网页引用，主报告需要按其要求再 open 对应来源。

- **[R1] Cui, Martin, Marzouk, Solonen & Spantini (2014).** Likelihood-informed dimension reduction for nonlinear inverse problems. Inverse Problems 30, 114015. https://arxiv.org/abs/1403.4680 ，DOI https://doi.org/10.1088/0266-5611/30/11/114015 。已浏览 arXiv abstract，工具 ref ``。支持 LIS 相对 prior 的低维 likelihood 变化，以及补空间 prior 与 Rao–Blackwell 思路。
- **[R2] Zahm, Cui, Law, Spantini & Marzouk (2022).** Certified dimension reduction in nonlinear Bayesian inverse problems. Mathematics of Computation 91, 1789–1835. https://arxiv.org/abs/1807.03712 ，DOI https://doi.org/10.1090/mcom/3737 。已浏览 arXiv abstract ``。支持 posterior KL、log-Sobolev 与积分 score 二阶矩的 certified reduction，不能简化成 single-anchor SVD。
- **[R3] Spantini, Cui, Willcox, Tenorio & Marzouk (2017).** Goal-Oriented Optimal Approximations of Bayesian Linear Inverse Problems. SIAM J. Sci. Comput. 39(5), S167–S196. https://arxiv.org/abs/1607.01881 ，DOI https://doi.org/10.1137/16M1082123 。已浏览期刊原论文 abstract ``，作者机构版本 https://dspace.mit.edu/entities/publication/6701f08c-9911-4e86-97c0-d4d35a8ce68f 。支持任务/prior/data 联合决定的 Bayes 近似，不支持将本稿简单 PCA 结果说成新发现。
- **[R4] Cotter, Dashti & Stuart (2010).** Approximation of Bayesian Inverse Problems for PDEs. SIAM J. Numer. Anal. 48, 322–345. https://arxiv.org/abs/0909.2126 ，DOI https://doi.org/10.1137/090770734 。已浏览 ``。成熟工作已把 forward approximation error 传到 posterior Hellinger 和 moments。本稿有界材料 TV/risk 分解是针对本任务的直接组合推导。
- **[R5] Davis & Kahan (1970).** The Rotation of Eigenvectors by a Perturbation. III. SIAM J. Numer. Anal. 7(1), 1–46. https://doi.org/10.1137/0707001 。已浏览原期刊 ``。支持 gap 控制的 subspace perturbation，不能把 spectral refresh 本身当全新矩阵理论。
- **[R6] Hardt (1980).** Semi-Algebraic Local-Triviality in Semi-Algebraic Mappings. American Journal of Mathematics 102(2), 291–302. https://www.jstor.org/stable/2374240 ，DOI https://doi.org/10.2307/2374240 。原论文检索 ``；JSTOR 全文页面未完整解析。Hardt 本人课程页 https://math.rice.edu/~hardt/501F07/ (``) 直接解释 finite strata、triangulation 及该文出处。
- **[R7] Whitney (1965).** Tangents to an Analytic Variety. Annals of Mathematics 81(3), 496–549. https://www.jstor.org/stable/1970400 ，DOI https://doi.org/10.2307/1970400 。已检索原论文 ``，全文页面未完整解析。此文原始论题主要为 analytic varieties，不能只靠标题把所有 semialgebraic 扩展归给这篇原文。
- **[R8] Kaloshin (2000 preprint).** A geometric proof of the existence of Whitney stratifications. https://arxiv.org/abs/math/0010144 ，PDF https://arxiv.org/pdf/math/0010144 。已读 PDF ``，Theorem 1 在印刷第 2 页，明确涵盖 semialgebraic / semianalytic sets，并说明 Whitney、Thom、Łojasiewicz 的历史归属。
- **[R9] Diaconis, Holmes & Shahshahani (2013).** Sampling From A Manifold. IMS Collections 10, 102–125. https://arxiv.org/abs/1206.6913 ，DOI https://doi.org/10.1214/12-IMSCOLL1006 。已浏览原研究入口 ``, ``。面积/coarea 和流形上合法测度是成熟方法。
- **[R10] Graham & Storkey (2017).** Asymptotically exact inference in differentiable generative models. AISTATS, PMLR 54. https://proceedings.mlr.press/v54/graham17a.html ，扩展版 DOI https://doi.org/10.1214/17-EJS1340SI 。已浏览原入口 ``。生成器约束下的 coarea 与 constrained inference 不是 P3-L 首创。
- **[R11] Weidling & Hohage (2017).** Variational source conditions and stability estimates for inverse electromagnetic medium scattering problems. Inverse Problems and Imaging 11(1), 203–220. https://www.aimsciences.org/article/doi/10.3934/ipi.2017010 。已读原论文摘要/定理入口 ``。支持连续 Maxwell inverse stability 的条件性与困难，但不能把其 full near/far field 条件直接套到 P1 有限 8×32 acquisition 上。
- **[R12] Frobenius, Ueber das Pfaffsche Problem.** Journal für die reine und angewandte Mathematik 82, 230–315，publisher DOI https://doi.org/10.1515/crll.1877.82.230 ，原期刊卷目录 https://www.degruyterbrill.com/journal/key/crll/1877/82/html (``)。本文的历史年份在其他目录有 1875/1877 的差异，按 publisher DOI/volume 记录引用；此处使用的是成熟 Frobenius involutivity theorem，显式 forward 反例为本轮推导。
- **[R13] Brascamp & Lieb (1976).** On extensions of the Brunn–Minkowski and Prékopa–Leindler theorems, including inequalities for log concave functions, and with an application to the diffusion equation. Journal of Functional Analysis 22(4), 366–389. https://www.sciencedirect.com/science/article/pii/0022123676900045 ，作者机构入口 https://collaborate.princeton.edu/en/publications/on-extensions-of-the-brunn-minkowski-and-pr%C3%A9kopa-leindler-theorem/ 。原论文摘要检索 ``，机构作者记录 ``；正文获取 403，未宣称阅读完整证明。本稿 conditional variance bound 属于该成熟不等式在材料函数上的应用。

## 15. 最终审查结论

P3-L 最稳健的数学对象是**原 prior 与 finite-noise likelihood 下的局部全维统计区域，配备低维 weak coordinates 与可积分的 strong conditional distributions**。局部 surface 可以是数值中心、可视化结构或受控的近似支持，但必须声明是哪一种。

对 P1 而言，认真研究 geometry 有理由：decoder 的 C^{2,1} 分片结构、材料/观测/先验三套尺度、active property boundary、及固定离散 Maxwell 的有理响应都可被利用。它们并不保证低维 posterior、GT 可达、简单拓扑或几何采样优于普通 blocked MCMC。

本轮应优先验证“概率正确的 strong integration + 材料 conditional mean + 计算误差预算”，然后检验 geometry 是否降低其总成本。若这些简单而正确的 Bayes 基线已足够，就不需要把 RL、Diffusion 或复杂流形 PDE 当成主线。
