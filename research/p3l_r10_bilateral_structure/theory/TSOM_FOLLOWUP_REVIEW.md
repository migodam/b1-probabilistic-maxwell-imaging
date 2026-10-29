# TSOM 跟进：能否绕开中间完整逆算子

独立零物理数学审查，2026-10-30。只核对代数、原材料约束与已有方法边界；无新增 forward/J、无新文献检索、无核心代码修改。主线程保留推导与最终判断。沿用 q=256、8 个源，所有源必须共享同一个原 Gaussianized66 decoder 产生的 D(x)。

**短结论：可以用 SOM/TSOM/CSI 风格的电流变量分裂，在反演迭代中避开每步 q 阶 LU；但不能把一个截断电流空间的最小二乘消元直接称为原正演的精确求逆。** 下面的消元是标准可分离最小二乘/variable projection，材料约化目标仍然非线性。实际节约、覆盖或概率资格均为 `NOT_RUN`。

## 1. 系数消元成立，但要说明 rank 与所消去的目标

固定正交列 Q∈C^(q×r)，令所有源电流组成矩阵

\[
h=h_d+Q\alpha,\quad K_D=I-DC,\quad T_D=K_DQ,\quad B_D=DE_R-K_Dh_d.
\]

若只最小化状态残差 \(\frac12\|T_D\alpha-B_D\|_F^2\)，则

\[
\alpha_*(D)=T_D^\dagger B_D,
\qquad
\alpha_*=(T_D^*T_D)^{-1}T_D^*B_D
\quad\text{仅当 }\operatorname{rank}T_D=r.
\]

Q 满列秩且 K_D 可逆是后一条件的充分条件；一般仍须检查 rank/条件数。实现应使用 QR/SVD 或正则化小问题，不显式求逆 normal matrix，因为 normal equations 会平方条件数。秩变化处伪逆和约化目标可能不光滑。

\[
\phi(D)=\frac12\|(I-T_DT_D^\dagger)B_D\|_F^2
\]

仍通过 D 依赖的投影而非线性。它并不是 \(\frac12\|F(D)-y\|^2\)。只有 \(B_D\in\operatorname{range}T_D\)，等价于精确状态 \(h(D)-h_d\in\operatorname{range}Q\)，才得到零状态残差并恢复该 D 的精确物理电流。正交 Q 不意味着 T_D 正交，也不能把 \(T_D^*T_D\) 换成 I。

在唯一内层最优、Q 与 h_d 固定、状态权重不随 D 变化时，envelope 给出一阶导数

\[
d\phi=-\operatorname{Re}\operatorname{tr}\big[r_s^*\,dD(E_R+Ch)\big],
\qquad r_s=K_Dh-DE_R.
\]

一阶导数不必显式求 dα，但精确 Hessian 或不同的近似 GN 仍涉及系数响应；不能把这称为把所有材料导数都消除了。

## 2. 若 Q 只是弱辐射空间，必须把数据项一起消元

记 Y=y−F_h。只有 **OQ=0 的严格接收零空间** 才使数据项对 α 不变。截断接收 SVD 留下的“弱方向”一般满足 OQ≠0，不能丢掉其数据影响。对联合目标

\[
\tfrac12\|W_d(Oh_d+OQ\alpha-Y)\|_F^2
+\tfrac\lambda2\|W_s(T_D\alpha-B_D)\|_F^2
+\tfrac\gamma2\|\alpha\|_F^2+\mathcal R(x),
\]

正确系数来自 stacked least squares：

\[
\mathcal A_D=\begin{bmatrix}W_dOQ\\\sqrt\lambda W_sT_D\\\sqrt\gamma I\end{bmatrix},
\quad
\mathcal b_D=\begin{bmatrix}W_d(Y-Oh_d)\\\sqrt\lambda W_sB_D\\0\end{bmatrix},
\quad\alpha_*=\mathcal A_D^\dagger\mathcal b_D.
\]

W_d 必须对应原复噪声的实化/白化约定；状态权重 λ 是算法惩罚，不能当新观测噪声。若权重或正则项随 D 变化，导数需包含这些依赖。

即使 OQ=0，把含噪数据的强奇异分量通过 h_d 精确钉死，也不等价于保留原有限噪声 likelihood：它把这些数据分量改成硬约束。可以优化可见电流系数，或保留完整数据误差并明确称为 current-assisted proposal。以数据构造的 h_d 导致这里的近似响应依赖 y，不能未经说明当作独立于数据的原 forward map F(D)。

## 3. Q 是投影后的像，且可缓存的是 CQ，不是变化的 B_D

\[
Z=(I-P_{visible})V_{internal},\qquad Q=\operatorname{orth}_{rank}(Z).
\]

此空间一般不是两个原空间的交集。应作一次 rank-revealing QR/SVD，去掉消失和相关列；通常 ZZ* 不是正交投影，不能拿来代替 QQ*。内部作用更直接的对象是 C(I−P_visible)，其谱基与“先对 C 取谱再投影”也不必相同。Fourier 基是一种可声明的选择，仍需投影、正交化和完整残差验证；不能仅靠名字宣布与原 TSOM 完全一致。

固定 C,Q 和本次观测的 h_d 时，可缓存

\[
G_Q=CQ,\quad E_d=E_R+Ch_d,\quad OQ,\quad Oh_d,
\]

并以

\[
T_D=Q-DG_Q,\qquad B_D=DE_d-h_d
\]

逐材料构造，无须重新做 C 的稠密乘法。**符号纠错：若 B 指本节右端 B_D，CB_D 并不固定；只有 B 被另定义为固定基时才可说“CB 固定缓存”。** 建议统一使用 CQ，避免与 R10 的 B=PCP 及本节 B_D 冲突。h_d 或 Q 若动态更新，相应缓存必须更新并计费。

每次主要成本包括 O(qr²) 的小列数 QR、O(qr n_s) 的多源操作、原 decoder/梯度与外层迭代；不是仅 O(r³)。这能避免 q 阶 LU，但在 q=256 的小规模上不自动更快，尤其 r 较大或迭代增多时。

## 4. 精确块 Schur 显示被删电流如何反馈

取完整酉基 U=[Q,Z_⊥]，将校正 h−h_d 写成 Q a+Z_⊥ b。定义 K_ab 为 U* K_D U 的对应块，b_Q=Q*B_D、b_⊥=Z_⊥*B_D。若 K_⊥⊥ 可逆，精确消去被删坐标得到

\[
(K_{QQ}-K_{Q\perp}K_{\perp\perp}^{-1}K_{\perp Q})a
=b_Q-K_{Q\perp}K_{\perp\perp}^{-1}b_\perp,
\]
\[
b=K_{\perp\perp}^{-1}(b_\perp-K_{\perp Q}a).
\]

输出同时变成

\[
F=F_h+Oh_d+OZ_\perp K_{\perp\perp}^{-1}b_\perp
+(OQ-OZ_\perp K_{\perp\perp}^{-1}K_{\perp Q})a.
\]

删掉 b 不只是少了若干直接辐射电流，还删掉了状态反馈、右端修正和可能的端口修正。即使 OZ_⊥=0，状态 Schur 反馈项仍通常非零；接收不可见不等于构成方程封闭。完整 K_D 可逆也不保证选定的主块 K_⊥⊥ 可逆，Schur 推导需要这个额外假设。

上面精确消元仍含通常随 D 变化的 (q−r) 阶逆，不自动便宜。反之，第1点的全状态 least squares 令 b=0 后解的是 Q*K_D*K_DQ 的 normal system，不是这个精确 Schur 补，也不等于一般的 Q*K_DQ Galerkin 系统。

## 5. 保留原66维材料与所有源共享，是电流辅助法的合法边界

应始终令 D=D(x)，x∈R^66 经过原 decoder；α 可逐源不同，但同一个 D 必须同时满足八列构成关系

\[
h_{i,:}=d_i(E_R+Ch)_{i,:}.
\]

不能逐源自由反推 d_is，也不能将逐格自由复数材料拟合后当作原 Gaussian64×Uniform2 prior。状态 residual 小只是当前离散方程的一种近似；原材料与最终数据资格需以原物理误差门或精确回退确认。

对原精确解，有

\[
F(D)-(F_h+Oh)=OK_D^{-1}(DE_R-K_Dh).
\]

因此，仅优化状态残差而没有 resolvent/端口稳定性控制，不保证原噪声尺度的数据误差小。若用于最终材料可能性集合，保留原噪声门、decoder 合法性及所有失败/回退成本。若用于准确 posterior 或原 MAP，还需证明消元/惩罚没有更换原 likelihood；有限 λ 的增广目标一般不具此保证。

## 6. 传统方法归属与可能的独立增量

电流作为辅助变量、测量方程与构成方程分开、接收 SVD 的确定/模糊分量、内部传播谱或 Fourier 的进一步筛选，属于传统 **SOM/TSOM/CSI 家族的变量分裂思想**。固定材料时线性消去电流系数是标准 separable least squares / variable projection。这里对 α 是线性的，和 R10 尚未实现的“完整非线性复幅度 profile”不是同一个问题。没有逐式核对 TSOM 原篇目标函数，因此应称 **TSOM-style 具体变体**，不能称上述每个 λ、regularizer、投影规则就是原算法原样。

可能检验的独立算法增量，是把该成熟结构嵌入原 66 维合法材料生成：共享精确 ROI 后端，按原白化 likelihood 的误差控制扩充 Q 或重建，保留被删电流反馈的校正/预条件，并比较原完整 LU/J 基线的整个 anchor＋候选＋最终资格成本和材料覆盖。必须同时面对成熟 TSOM/CSI 或相同变量分裂强基线。单纯“对电流做 TSOM”“消去 α”“固定 CQ 缓存”都不足以形成新颖性。全部这些跟进性能为 `NOT_RUN`。

已有来源与阅读资格：[公式集02 F-OP-004—007](../../maxwell_formula_atlas/02_FIXED_OPERATORS_AND_PROJECTIONS.md)、[公式集05 F-ERR-003—008](../../maxwell_formula_atlas/05_STATIC_OPERATOR_ERROR_CONTROL.md)、[R10 先例账本](../PRIOR_ART_AND_NOVELTY.md)。原 TSOM DOI [10.1088/0266-5611/25/8/085003](https://doi.org/10.1088/0266-5611/25/8/085003)，FFT-TSOM DOI [10.1109/TAP.2010.2103027](https://doi.org/10.1109/TAP.2010.2103027)；沿用此前身份核验与作者综述/摘要证据，不新增全文核验声明。
