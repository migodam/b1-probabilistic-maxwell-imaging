# Agent B/D 独立报告：A17、反馈导数与辅助电流的材料曲率

日期：2026-10-31。本文是独立数学推导、原文审计与小型复矩阵核验；没有执行新的高对比度成像实验。

## 0. 核心判断

最值得继续发展的结构不是任意选择一个指数或 tanh 权重，而是**区分物理电流切向与偏离共享材料约束的辅助电流法向**，再研究后者怎样与 CIE 的有限材料更新、当前原始数据残差及电流空间释放共同作用。

两种“增加电流方向”不等价：

- A17 的 Galerkin tangent 修改，改变对同一个物理材料导数的近似；材料法方程曲率差通常不定，且必须保留两个 pullback channel。
- CSI/CIE 中释放自由辅助 current，改变一个松弛问题；在冻结线性化、精确剖面化条件下，新增自由度从 material Hessian 中减去半正定项。这有精确的 step 稳定界，但减掉的可能是原 reduced model 的虚假刚度，不能解释成必然丢失真实测量信息。

与 Agent A 的 CIE 有限材料读出恒等式、主研究者推导的二阶/四阶激活规律相结合，可以形成一条比“有界反馈更稳定”具体得多的研究主线。以下证明并不建立历史首次性；Schur 消元、可逆变量变换、最小二乘剖面化本身都是已有工具。

## 1. 实际原文及版本

A17 仓库锁定 SHA：2ac81b989a2f13d62df9d5a7abd690a33c2189c7。

已取得并读取：

1. [CORE_THEORY.tex](https://github.com/migodam/a17-maxwell-certified-exchange/blob/2ac81b989a2f13d62df9d5a7abd690a33c2189c7/closure_final/manuscript/CORE_THEORY.tex)：物理 tangent、Galerkin Schur、双 material pullback。
2. [CONDITIONAL_EXCHANGE_THEORY.tex](https://github.com/migodam/a17-maxwell-certified-exchange/blob/2ac81b989a2f13d62df9d5a7abd690a33c2189c7/closure_final/manuscript/CONDITIONAL_EXCHANGE_THEORY.tex)：同维 exchange、有限二次 gain 及误差界。
3. PAPER_TAP_FINAL.tex、PAPER_SUPPLEMENT_FINAL.tex，以及其 FINAL_EXCHANGE_METHOD、FINAL_PHYSICAL_WITNESSES、CLOSURE_RESULTS、FINAL_DESIGN_DISCUSSION、FINAL_DDA_MODEL、CLOSURE_SUPPLEMENT_COST 输入。
4. A17-X 正式原报告 [RESEARCH_A17X_FULL_REPORT_ZH.md](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/6659aa49dff2535996d40962791e9a5a6e974adc/research/a17x_f3_competition_audit/pro_source/report/RESEARCH_A17X_FULL_REPORT_ZH.md)。
5. A17-X 后续审计 THEORY_AND_SOURCE_COMPATIBILITY、MECHANISM_REVISION、RESEARCH_REPORT_ZH。

主仓库 SHA：6659aa49dff2535996d40962791e9a5a6e974adc。A17X 原文已包含一般多源 current 可实现性、minors、row material readout、精确 defect-output 公式、有限双侧 exchange、有限响应尾界及探索任务驱动 current selector 的建议。不能将这些列作本轮首次推导。

A17 的物理实验为三维 vector DDA，真实本构是 polarizability a(χ)，注入必须保留 a'(χ)。本报告以下简化到任务书二维线性 X=diag(m) 时明确如此；不能把它冒称原 DDA 的逐式同一个实现。

## 2. 固定记号与一般反馈精确等价

取 n 个标量材料单元、p 个 source、r 个 receiver：

\[
j=XQ,\quad Q=E+Gj,\quad L=I-XG,\quad F=Sj,\quad R=F-Y.
\]

材料 m 为实变量，复数矩阵的实微分内积为 Re tr(A*B)。L 可逆。若 X=diag(m)，沿实方向 h，

\[
T_h=Dj[h]=L^{-1}X_hQ,\quad DQ[h]=GT_h.
\]

### 2.1 非交换反馈的正确顺序

令 B 可依赖 m，但在每个状态满足 I+BX 可逆，定义

\[
V=Q+Bj=E+(G+B)j,\qquad
A_B=X(I+BX)^{-1}=(I+XB)^{-1}X.
\]

则

\[
j=A_B\{E+(G+B)j\}.
\]

并有精确身份

\[
I-A_B(G+B)=(I+XB)^{-1}(I-XG).
\]

证明只需左乘 I+XB。这个身份表明反馈表示的状态方程与原 Maxwell 方程等价，而且原物理极点没有因换表示消失。变换自身产生的奇点必须另行排除。

B=βI 时 A_B=m/(1+βm)，乘 β 得 CIE modified contrast。**这说明 CIE 和 A17 均可放在 resolvent/Schur 框架内，不等于已经证明 CIE 是 A17 特定 current-exchange 定理的特殊情形。** A17 的 Schur 矩阵由选定子空间决定；CIE 的 β 是另一个自由选择。

非局部 B 一般使 A_B 稠密。若在 joint-CIE 中把这个稠密 A_B 当任意矩阵优化，就把局部材料的 n 个自由度扩大为非局部响应的 n² 个自由度，必须重新施加可实现性。不能保留 scalar modified-contrast 的便宜材料更新，同时偷偷放松其非局部结构。

### 2.2 “等价”对新权重几乎没有约束

若 X 可逆，任意可逆 M 都可写成 M=I+BX，只需 B=(M-I)X^{-1}。因此精确反馈等价性自身不能证明 M 是有物理意义的设计，更不能证明其搜索性能。

若 B 随 m 改变，

\[
DV[h]=(G+B)T_h+DB[h]j.
\]

遗漏第二项会改变优化目标。若 B 由 residual、current 或子空间间接计算，DB 必须包括这些链式项；也可声明外循环冻结 B、内循环优化一个固定 surrogate，二者不是同一算法。

### 2.3 Reciprocity、passivity 与 Cayley 的适用边界

若 X 和 B 复对称，且逆存在，则 A_B^{-1}=X^{-1}+B 复对称，保持 reciprocal 代数结构。Reciprocity 使用转置，不能与 Hermitian 条件混淆。

若采用 e^{-iωt}，Im X=(X-X*)/(2i)⪰0，并且 B=B*，则

\[
\operatorname{Im}A_B
=A_B^*X^{-*}(\operatorname{Im}X)X^{-1}A_B\succeq0.
\]

这是纯矩阵计算；仍要求 X 与变换可逆，并且实际离散量 X 真正满足相应 passivity convention。对实正 X 和 B⪰0，I+BX 在 X 诱导度量下有良性谱，但欧氏奇异值可很小：取
X=diag(10^6,10^{-6})，
B=[[10^{-6},1],[1,10^6]]，
则 I+BX=[[2,10^{-6}],[10^6,2]]，
特征值为1、3，而最小奇异值约3×10^{-6}。仅看 eigenvalues 会漏掉非正规放大。

若一个明确的算子 Z 在指定内积中 accretive，C=(Z-I)(Z+I)^{-1} 才满足

\[
I-C^*C=4(Z+I)^{-*}\operatorname{Re}Z\,(Z+I)^{-1}\succeq0.
\]

Maxwell 的 X^{-1}-G 在高对比度和近共振下并不自动 accretive。即使这个 Cayley 变换 contractive，也不等于材料 inverse iteration contractive。

## 3. Gram 反馈梯度：完整重推导

令 V=M(m)Q，A=V*V≻0，c 为固定正数，

\[
f_M=\frac c2\operatorname{tr}(RA^{-1}R^*),\quad
W=RA^{-1},\quad D=W^*W.
\]

由 d(A^{-1})=-A^{-1}(dA)A^{-1}，

\[
Df_M[h]
=c\operatorname{Re}\langle W,DR[h]\rangle_F
-c\operatorname{Re}\langle VD,DV[h]\rangle_F.
\]

所有材料导数均为实方向；这已处理复数共轭，不能再额外乘2。

### 3.1 题中的 -cKᵀd 成立条件

若 M=diag(exp φ_i)、φ_i 实值、M(m0)=I，且 K_ij=∂φ_i/∂m_j，

\[
DV[h]=DQ[h]+\operatorname{diag}(Kh)Q
\]

在 m0 成立。两个目标的 DR、DQ 相同，故

\[
D(f_M-f_0)[h]
=-c\sum_i(Kh)_i(QDQ^*)_{ii}.
\]

由于 Q†=(Q*Q)^{-1}Q*，

\[
(QDQ^*)_{ii}=\|(RQ^\dagger)_{:,i}\|_2^2=d_i,
\quad
\boxed{\nabla f_M-\nabla f_0=-cK^\top d.}
\]

若 c=c(m)，需再加 (∇c)(f_M/c-f_0/c)；在两目标值参考点相同，该额外差为零，但各自梯度仍含 ∇c。相同 material regularizer 相减后消失，不同正则则不能忽略。

使用 A=Q*Q+λH_s 时，结论仍有同形表达，但 d_i 应改为 \|(RA^{-1}Q*)_{:,i}\|²，不能称其 Q†。

### 3.2 标量 CIE 归一化问题

若参考材料为常数 b，归一化 M_i=(1+βm_i)/(1+βb)，则 K=β/(1+βb)I。

原始 M_i=1+βm_i 未归一化时，f_raw(b)=(1+βb)^{-2} f0(b)，且

\[
\nabla f_{\rm raw}(b)
=(1+\beta b)^{-2}\{\nabla f_0(b)-c\kappa d\},
\quad \kappa=\beta/(1+\beta b).
\]

不能把未归一化 raw 目标直接代入参考点 M=I 的恒等式。若 b 空间变化，K为相应对角而非一个标量乘 I。

### 3.3 固定 B 与动态 B 的导数

固定 B 时 V=E+(G+B)j，故

\[
\boxed{Df_B[h]
=c\operatorname{Re}\langle S^*W-(G+B)^*VD,T_h\rangle_F.}
\]

动态 B 还需

\[
-c\operatorname{Re}\langle VD,DB[h]j\rangle_F.
\]

这个固定目标的导数式支持题中的 A17 梯度预测差，但仅在相同 full anchor、相同 V/W/D、指定 reduced tangent 比较时成立。

### 3.4 可逆 source 混合

对任意可逆 C，若 E→EC、Y→YC，则 Q→QC、R→RC，

\[
(RC)(QC)^\dagger=RQ^\dagger,\quad
f_M(RC,QC)=f_M(R,Q).
\]

满列秩是关键。固定 λI ridge 破坏一般 GL(p) 不变性；只有 unitary source mixing 仍无问题。若需要一般可逆源坐标不变，ridge metric 必须按 H_s→C*H_sC 协变。

这个不变性只说明改换已有源的基。它既不增加观测信息，也不能代替新增物理 illumination。

## 4. K_field 的准确含义

设 H0=Q0*Q0，Z=Q0H0^{-1/2}，题中 K_field 的 numerator 为

\[
\operatorname{Re}\langle Z_{i,:},(GT_{e_j}H_0^{-1/2})_{i,:}\rangle.
\]

### 4.1 源基不变性成立

在 Q0→Q0C 后，新的白化 Z'=ZU，其中
U=H0^{1/2}C(C*H0C)^{-1/2} 为 unitary。
同一 GT_ej 白化后也右乘 U，故逐行内积和行范数不变。分母使用固定相同 ε 时，K_field 也不变。

### 4.2 它不是 projector/leverage 导数

上式冻结了 H0^{-1/2}。真正 P_Q=Q(Q*Q)^{-1}Q* 的导数为

\[
DP_Q[h]=(I-P_Q)DQ[h]Q^\dagger+\text{adjoint}.
\]

因此

\[
\frac12D(P_Q)_{ii}[h]
=\operatorname{Re}\langle Z_{i,:},
[(I-P_Q)GT_hH^{-1/2}]_{i,:}\rangle.
\]

若 DQ=QA，span 完全不变，DP=0，但 K_field 通常非零。它测的是在当前 source 白化规范下的局部场强变化，包含 source-span 内增益；不能用它代替 illumination-subspace rotation。

ε 只防止数值除零，不是严格处理总场节点。节点处应使用带阈值的未归一化 feature 或保留该材料变量，明确读出不稳定性。3D DDA 还必须按每单元3×p block 和物理体积度量处理。

### 4.3 物理生成不等于优化有益

这个 K 来自 Maxwell response，因此不是随便造的一组数字。但它不包含 residual 的改善符号，不能由此推出 -Kᵀd 是良好方向。任意 K 在 d≠0 时能产生任意梯度修正 g：取 K=-d gᵀ/(c||d||²)。这是 gradient-shaping 的明确 no-go。

## 5. 有界图册与 jet 对照的真实结论

若 |φ_i|≤α，则

\[
e^{-2\alpha}Q^*Q\preceq Q^*M^*MQ\preceq e^{2\alpha}Q^*Q,
\quad
e^{-2\alpha}f_0\le f_M\le e^{2\alpha}f_0.
\]

它保证目标值尺度等价，**不保证梯度、驻点、basin 或迭代收缩**。有界光滑函数可振荡得很快。

在参考点 b，w_CIE=(1+βm)/(1+βb)、w_exp=e^{κ(m-b)}、w_bound=exp[αtanh(κ(m-b)/α)]，κ=β/(1+βb)，三者一阶相同。w_CIE''(b)=0，而 w_exp''(b)=w_bound''(b)=κ²。所以 exp 与 bounded 已经匹配二阶 jet；它们的差从三阶开始。用“exp/bounded二阶不同”解释实验会错误。

一单元 Maxwell-like 模型提供反例。取 E=S=1，G=g，物理 F=m/(1-gm)，观测对应 m*=y/(1+gy)。Gram 反馈目标化为

\[
f_M(m)=\tfrac12(1+gy)^2(m-m_*)^2/M(m)^2.
\]

导数符号额外含 1-(m-m*)φ'(m)。正 β 的 CIE 在正材料域中该因子=(1+βm*)/(1+βm)>0；指数权重则在 m=m*+1/κ 产生非真实驻点，并在其右边将下降推向更大材料。有限 box 上可形成错误边界极小值。这个例子不声称同一锚点的连续梯度流一定走到不同解；它证明相同初始梯度不够推断远处景观更好。

可逆重参数化 m=ψ(z) 保留 interior critical points、极小值类型和 exact fibers；它可以改变离散迭代动力学，但没有消灭原目标坏极小值。

## 6. A17 的两个材料通道与辅助目标增益

原文给

\[
\Delta T_h=\bar C\Gamma_C^{-1}N_C(h),
\quad \Delta J=Y_C\Gamma_C^{-1}N_C.
\]

在同一个 frozen real material quadratic 中，

\[
s_C-s_0
=-M_0J_0^*Y_C(\Gamma_C^{-1}N_Cs_C)
-M_0N_C^*\Gamma_C^{-*}Y_C^*(r+J_Cs_C).
\]

第一个通道反映旧敏感度与新增输出的 Hessian 交叉，第二个通道反映当前 residual 对新材料注入的激活。只用 receiver sensitivity、current energy 或 ΔJ* r 都不够。

固定 feedback 目标时，令 Z_B=S*W-(G+B)*VD，则

\[
\Delta Df_B[h]
=c\operatorname{Re}\langle Z_B,\bar C\Gamma_C^{-1}N_C(h)\rangle_F.
\]

若 reduced state 本身改变，V/W/D 也变，上式不再是总梯度差。应重新计算或加相应误差项。

原 full−old tangent error分成已解析 exchange 与 full−new tail；Schur 更新不消除 tail。给定 tail action E_h 的范数上界 ε_h，有

\[
|\Delta Df_{\rm full-reduced}[h]|
\le c\|Z_B\|_F\epsilon_h.
\]

这只是针对该方向该目标的界；不覆盖遗漏 DB、状态误差或任意未测试方向。

原文有限 quadratic gain

\[
Q=-g_F(x)^Td-\tfrac12d^TH_Fd
\]

仍是比“最大场变化减小”更直接的任务判据。一次 full Jd 支付 p 个 source RHS；“rank很低”不等于这些验证免费。A17 原报告还明确：同空间 ordinary Galerkin 与 Schur 更新精度相同，且原 Krylov 对照在部分历史 solver 任务上更强。

## 7. 本轮候选命题：辅助电流释放的精确 profile 更新

### 7.1 假设和模型

在一个物理一致 anchor (m0,j0)，实化所有复输出与 current 系数，固定正定 material damping Λ，固定 box 活动面。若活动面变动，以下公式只在新的自由坐标重新建立后适用。

对 current basis U 和 material step h，使用局部 joint objective

\[
\frac12\|r+SUc\|^2
+\frac\rho2\|LUc-B_mh\|^2
+\frac12h^T\Lambda h+\ell^Th,
\quad B_mh=X_hQ_0.
\]

写成
A_U=[SU;\sqrt\rho LU]，
b(h)=b_0+\mathcal B h=[r;-\sqrt\rho B_mh]。
令 P_U 为 A_U 的欧氏正交 projector，profile current 后

\[
\psi_U(h)=\tfrac12\|(I-P_U)(b_0+\mathcal Bh)\|^2
+\tfrac12h^T\Lambda h+\ell^Th.
\]

CIE 中可用冻结状态权重 Wβ 替代 √ρI，即 A_U=[SU;WβLU]，mathcal B=[0;−Wβ B_m]。在物理 anchor，D Wβ 的项乘 state residual=0，可以在一阶线性化中消失；离开物理 manifold 后必须保留它。

### 7.2 增广公式

新增 current basis C，令

\[
Z=(I-P_U)A_C,\quad Z^TZ\succ0,
\]
\[
U_m=(Z^TZ)^{-1/2}Z^T\mathcal B,\quad
v=(Z^TZ)^{-1/2}Z^Tb_0.
\]

则

\[
\boxed{\psi_{U+C}(h)=\psi_U(h)-\tfrac12\|v+U_mh\|^2.}
\]

证明：A_[U,C]的 column space 是 Range(A_U) 与 Range(Z) 的正交直和，故 P_new=P_U+Z(ZᵀZ)^{-1}Zᵀ，代入即可。

于是

\[
\boxed{H_n=H_o-U_m^TU_m,\qquad g_n=g_o-U_m^Tv.}
\]

这是自由 nuisance current 的 elimination，和物理 tangent approximation 的 A17 不定 ΔH 完全不同。

### 7.3 精确 step 差及稳定界

令 p=-H_o^{-1}g_o、p_n=-H_n^{-1}g_n，若 H_n≻0，

\[
\boxed{p_n-p=H_n^{-1}U_m^T(v+U_mp).}
\]

由两 normal equations 相减即可证明。令

\[
\eta=\|U_mH_o^{-1/2}\|_2^2<1.
\]

则 H_n⪰(1−η)H_o，并有

\[
\boxed{\|p_n-p\|_{H_o}
\le\frac{\sqrt\eta}{1-\eta}\|v+U_mp\|_2.}
\]

证明写 A=U_mH_o^{-1/2}，
H_n=H_o^{1/2}(I−AᵀA)H_o^{1/2}，
左端为 ||(I−AᵀA)^{-1}Aᵀ(v+U_mp)||，用奇异值上界即可。

它给出 β/ρ/current-mode-count 联合决策的一个实际尺度：候选不是按 current energy排序，而是按其造成的 step变动、residual relief及原数据finite gain排序。η阈值不能单独充当最优selector。

### 7.4 不能回避的反例：虚假刚度

取一单元 G=0、E=1、S=ε、L=B_m=1。无自由current时 H_o=ρ+λ；释放current并消去后

\[
H_n=\frac{\rho\epsilon^2}{\rho+\epsilon^2}+\lambda.
\]

当 ε→0，物理数据本来就不辨识材料。释放current只是去掉了“把current冻结”的虚假state-penalty刚度。死守 η≪1 会阻止正确模型，巨大的 λ 也可以掩盖问题。

所以本命题只称 **profile step 稳定控制**，不称 Fisher 信息保留定理。需要允许逐步释放、调节 ρ/β、重锚和独立物理corrector。原数据 J 的信息与最终真值误差另行评价。

### 7.5 精确物理 tangent/normal 分解

定义实线性 B_m:h↦diag(h)Q。每行 Q_i≠0 时

\[
(B_m^\dagger Z)_i=
\frac{\operatorname{Re}\langle Q_i,Z_i\rangle}{\|Q_i\|^2}.
\]

令 P_B=B_mB_m†，则 current-space projector

\[
P_{\rm phys}=L^{-1}P_B L
\]

满足 P_phys²=P_phys，Range(P_phys)=Range(L^{-1}B_m)，而它在内积
\(\langle V,W\rangle_L=\operatorname{Re}\langle LV,LW\rangle_F\)
中正交。相应 normal defect为

\[
N(V)=(I-P_B)LV.
\]

计算 candidate 的 normal signature 只需 L action 和逐行投影，无需构造完整 material Jacobian；要真正构造 P_phys V 则还需 L solve。多源必须共同投影到同一个**实**标量行，而不是每source独立恢复材料。

## 8. 与 CIE 有限材料偏置的二阶/四阶激活

本节有限读出恒等式由文献 Agent 独立提出，二阶/四阶规律由 root 提出，我重新证明并独立数值核验。

**实现边界：以下主公式针对 real-material / real-R 受限最小二乘 block。2016 原文 Eq.14 是 complex-R 更新，不能把此主公式不加限定地称为原文完整算法的实际轨迹。** 若具体代码先做complex-R LS再投影实材料，必须使用下述复版本，不得混用实/复行投影。

对一行 q,j，a=||j||²、b=Re⟨q,j⟩、c=||q||²>0。CSI real-material readout m_CSI=b/c。CIE fixed-current modified-contrast LS 经逆变换得到

\[
m_\beta=\frac{b+\beta a}{c+\beta b},
\quad
\boxed{m_\beta-m_{\rm CSI}
=\frac{\beta\|j-m_{\rm CSI}q\|^2}{c+\beta b}.}
\]

β>0且分母正时，偏置向高材料方向；只有 off-manifold shared-material defect 非零才激活。这不是对完整CIE迭代每一步都单调增大的断言，且受material box/regularization修改。

在物理anchor j0=m0q0，沿linear current路径 j(t)=j0+tV、q(t)=q0+tGV，

\[
e_{\rm CSI}(t)=t\,\Pi_{q_0}^{\perp,\mathbb R}(LV)+O(t^2).
\]

故
\[
m_\beta(t)-m_{\rm CSI}(t)
=t^2\frac{\beta\|\Pi_{q_0}^{\perp,\mathbb R}LV\|^2}
{c_0(1+\beta m_0)}+O(t^3).
\]

若 V=T_h 为真实 physical tangent，LV=hq0。记 g=GV，则有**精确**

\[
j(t)=(m_0+th)q(t)-t^2h g,
\]
\[
e_{\rm CSI}(t)=-t^2h\,\Pi_{q(t)}^{\perp,\mathbb R}g.
\]

因此

\[
\boxed{m_\beta(t)-m_{\rm CSI}(t)
=t^4\frac{\beta h^2\|\Pi_{q_0}^{\perp,\mathbb R}GT_h\|^2}
{c_0(1+\beta m_0)}+O(t^5).}
\]

各行分别成立；需 q0非零、实material interior、β有界且分母远离零。沿一条精确物理材料路径，e_CSI恒零，所有β恢复同一材料，区别完全消失。

意义：CIE的有限读出机制在自由normal current上是二阶，在linear physical tangent上是四阶。这为“为何A17 tangent attenuation不能直接解释CIE basin收益”提供清晰的局部数学解释；它不推导全局CIE优越性。

若允许复R，改记复数b=⟨q,j⟩，则

\[
m_\beta^{\mathbb C}=\frac{b+\beta a}{c+\beta\bar b},\quad
m_\beta^{\mathbb C}-m_{\rm CSI}^{\mathbb C}
=\frac{\beta\|j-m_{\rm CSI}^{\mathbb C}q\|^2}{c+\beta\bar b}.
\]

若随后只取实部，其差为

\[
\operatorname{Re}(m_\beta^{\mathbb C}-m_{\rm CSI}^{\mathbb C})
=\frac{\beta(c+\beta\operatorname{Re}b)}{|c+\beta\bar b|^2}
\|j-m_{\rm CSI}^{\mathbb C}q\|^2.
\]

二阶/四阶结论相应使用complex span的正交投影。一个纯虚材料方向 V=i h q 在G=0时违反真实材料约束，却属于complex-material tangent，因此不激活这个complex版本的defect bias。不能把“real正常方向必有bias”推广到所有读出实现。

## 9. 有限响应候选：可精确枚举射线，但不再包装成新 A17-X

令合法 material ray m(t)=m0+t h，h 支持于 s 个单元，ΔX=t UDU*。定义

\[
P=L_0^{-1}U,\ H=U^*GP,\ O=SP,\ B=U^*Q_0.
\]

则

\[
R(t)=R_0+tO(I-tDH)^{-1}DB.
\]

这是已存在A17-X/T-matrix公式。新报告可补充“线内全局搜索”的明确代数形式：d(t)=det(I−tDH)，
q(t)=R0 d(t)+t O adj(I−tDH)DB，
N(t)=||q(t)||²，
D(t)=|d(t)|²。
在无极点实紧区间，

\[
\Phi(t)=N(t)/(2D(t)).
\]

N,D 的次数≤2s，导数分子 N'D−ND' 次数≤4s−2（最高同次项相消）。全部实驻点与可行端点给出该射线的全局最小值；不等于整个材料问题全局解。多项式数值求根未经过interval validation时，也不能称严格认证。

成本为 s 个传播RHS建立P，后续每t求s×s反馈系统。一个低维latent仍可能使ΔX满秩，旧A17X审计已有7维latent但575个非零材料算子项的反例。

approximate current ktilde(t) 的状态残差

\[
\rho_t=L_t\tilde k-tX_hQ_0
\]

给
\[
F_{\rm true}(t)-F_{\rm approx}(t)=-SL_t^{-1}\rho_t.
\]

若有真output error bound ε(t)，residual norm在
[max(0,||Rapprox||−ε),||Rapprox||+ε] 内。可据候选原始损失区间是否重叠决定是否付费增广current。这个任务型selector的总体想法旧A17-X已经提出，本轮应检验它能否和normal-current分支机制形成独立增量，而非重新命名。

## 10. 两个完整候选算法六元组

### 候选 I：按物理/法向分解释放辅助电流，并控制 material profile step

- 优化变量：完整实material m、共享current basis中的复系数 a；材料不受任意深latent限制。
- 物理约束：box/passivity；state constraint最终由完整Maxwell校验；中间允许带预算的state defect。
- 目标：原始data residual + 固定CIE/CSI state residual + material prior；可用AL multiplier推进可行性。
- 更新：profile或近似profilecurrent，material GN/LS；候选current增广用第7节精确profile差与step界预测，top candidates用原data finite gain比较。
- 子空间：FFT-TSOM backbone；A17生成的dressed atoms及普通Krylov/dual atoms共同竞争；normal signature N(C)区分物理与松弛角色。
- Globalization：β/ρ和mode-count联合continuation；有限normal分支保留多个候选，局部物理corrector最后重投影；不要求每个探索态原data单调。

伪代码：

1. 从相同数据和成熟CIE/FFT初始化，保留原full material维度。
2. 在物理anchor上冻结当前chart，获得state/data residual和candidate current atoms。
3. 对每个atom或小组计算N(C)与其material读出signature；排除重复及不稳定Schur核心。
4. 对当前joint线性子问题计算profile step、candidate U_m,v,η；记录实际step-change bound，不以η单一排名。
5. 同时比较少量 (β/ρ,current-set) actions；若增广只是揭去假刚度，允许必要释放并重锚。
6. 用独立full data评价有限material动作或物理corrector后的候选；保留有数据质量及结构差异的少量分支。
7. 若normal defect已小到容许值，切回物理局部求解；继续原完整data验收。

成本：current residual和G actions；每新增模式形成S C、L C及小Gram；material solve维度仍n；full Jd validation每方向支付p RHS。只有模式数远小于n且多个候选摊销同一anchor时才可能有优势。原ordinary goal-ROM同样允许factorization与Schur复用。

失败：normal freedom过多只拟合data却不物理；虚假刚度门阻止必要自由度；FFT外atoms提高成本但不改善；β固定图册导数遗漏；branch count扩大而无同预算增量。

最小决定实验：固定小型complex VIE、同FFT基础和同支付actions，比较原progressive FFT、无门的normal enrichment、profile-step门、A17 tangent enrichment；同一初始gradient和same-pool条件下，记录normal defect、material误差、原data residual、current释放造成的阶次效应。尚不能直接跳到576参数全面campaign。

### 候选 II：小支持有限物理响应的多候选搜索与任务型current计算

- 优化变量：完整m中的小patch update D，或标量ray t；多个合法候选分支。
- 物理约束：真实m box、真实本构；反馈resolvent invertible；source共享自动由完整material模型保持。
- 目标：原始Maxwell data residual及声明过的regularizer。
- 更新：exact有限有理响应枚举多个驻点；或small constrained D优化；候选原物理corrector。
- 子空间：material小patch不是全局latent；内部current由exact-port/full solver/ordinary ROM/A17 selectable ROM分别计算。
- Globalization：保留多个分离有限低谷、不同patch/方向及branch重访；同预算普通多起点作为强对照。

这个算法与候选I真正不同：它始终在物理material上搜索，不依靠off-manifold current freedom。但其核心已有A17-X/R4先例，R3独立贡献不能仅来自把射线算得更准。适合作为机制对照和次路线。

最小实验：s=2或3的精确射线全部驻点，验证true objective与surrogate tail区间；同候选同current空间比较A17与ordinary ROM；若精度相同则只能比较setup/selection成本。主任务是区分finite-physics收益与更多搜索次数。

## 11. 已完成代数检查

脚本：check_a17_math.py。固定seed=20261031，n=6、p=3、r=4，随机复对称G和实正m；是 Maxwell-like VIE 合成矩阵，不是真实天线/Green网格、不是成像benchmark。

| 检查 | 结果 |
|---|---:|
| -cKᵀd 与中心差分梯度 | 1.33×10^-10 相对/尺度误差 |
| K_field 可逆source混合 | 3.87×10^-16 |
| Gram objective source混合 | 2.22×10^-16 |
| 固定B导数与中心差分 | 2.00×10^-10 |
| profile Hessian低秩差 | 1.22×10^-16 |
| profile gradient差 | 4.37×10^-16 |
| 精确profile step差 | 1.65×10^-16 |
| η | 0.1790654 |
| 实际step H-norm / 上界 | 0.2628524 / 0.3539512 |
| two-support finite response与直接求解 | 1.44×10^-15 |
| s=2射线导数多项式次数 | 6，即4s−2 |
| 任意current的CIEbias最后二分阶次 | 2.00041 |
| physical tangent的CIEbias最后二分阶次 | 4.01631 |
| finite CIEbias identity误差 | 约10^-16 |

四阶项较小时直接相减容易浮点消失，所以阶次用代数等价的非负defect表达计算，同时独立核对两readout之差，不能把floating underflow当零机制。

## 12. 建议给主报告的主张边界与优先级

优先核心：**共享材料约束定义的current normal geometry + CIE finite defect bias + 二阶/四阶激活分离**。它把原A17物理tangent角色与成熟CIE的辅助current作用真正连接起来，同时避免强行证明前者“包含”后者。

第二核心：**joint-current profile step更新与可控自由度释放**。它提供可执行的β/current-subspace联合控制，但只是局部条件定理，不能称全局收敛或真实信息保持。

有限响应、多源minors、Schur双侧修饰、任意current读出都应归还A17-X与既有T-matrix/CSI背景。任意exp/tanh与K_field建议降为E0机制对照；在没有normal-current机制增量前，不作为主论文算法。

FBG-R2的8/12应承担“表明辅助几何值得研究的正面动机”，不能承担“已优于成熟CIE”或“已证实A17独立收益”的结论。

最终建议继续**辅助电流与材料实现性之间的非线性几何**，保留反馈作为控制离物理流形探索的图册，保留A17作为物理切向/有限响应计算的任务接口。不要让论文主线回到寻找一个更好看的M(m)函数。

