# 05 · 从固定算子直接推到可校正的压缩算法

本章继续推进“固定载体”的数学后果。目标不是把材料固定在某个弱方向，而是让固定算子承担离线计算，让变化材料只进入必要的小系统，并让误差控制面向实际观测/任务。以下推导不需要再跑一个恒等式数值测试。尚未实施的算法与已运行的精确ROI消元分开。

## F-ERR-001 · 一次对角化为何不能消除所有非线性求解

精确固定宿主消元后写为

$$L_Dh=b_D,\quad L_D=I-DC,\quad b_D=DE_R,\quad F=F_h+Oh,$$

其中 $C,E_R,O,F_h$ 固定，$D=\operatorname{diag}(m_{ROI}-m_h)$ 变化。对于任意两格，

$$[D,C]_{ij}=(d_i-d_j)C_{ij}.$$

因此，在物理 Green 有跨格耦合时，一般非均匀材料不与 $C$ 对易。若一个固定子空间要对**所有对角材料矩阵**都不变，它必须由一组坐标轴张成：证明是让 $D$ 依次取每个坐标投影，子空间内向量的每个非零坐标分量均须留在子空间。若还要求对连接这些格点的稠密 $C$ 不变，除空集/全空间外通常没有这样的公共子空间。

**限制范围：** 上述是对“全部对角材料族”的结构结果，不是证明 P1受限decoder像中不存在特殊小不变子空间。原P1未给出这样一个共同不变性，因此不能无条件沿用均匀幅度的对角化。

均匀材料 $D=vI$ 则不同。对 $C=QTQ^*$ 做一次复Schur后，

$$F(v)=F_h+(OQ)v(I-vT)^{-1}(Q^*E_R).$$

每次仅解三角系统，保留全部256维，不是rank-one更新。此公式已用于本轮R6条件模型；一般Full66材料只能复用固定 $C$，不能复用同一个三角resolvent。

## F-ERR-002 · 固定低秩 Green 近似与小系统

取固定离线近似 $\widetilde C=UV^*$，$U,V\in\mathbb C^{q\times k}$。可以来自截断SVD、物理端口Krylov或其它明确的拟合；此处不预设哪种最好。

由 Woodbury，若 $I-V^*DU$ 可逆，

$$\widetilde h=(I-DUV^*)^{-1}DE_R
=DE_R+DU(I-V^*DU)^{-1}V^*DE_R.$$

这是真实近似算子的精确解。固定 $U,V$；只有 $k\times k$ 系统依赖材料。其复杂度包含 $O(qk^2+k^3+qkn_s)$，不能只报 $k^3$。低秩近似支持仍是全256格对角材料；它没有改变 decoder 或 prior，但改变了forward，故未经校正不能作为原likelihood。

## F-ERR-003 · 原 Maxwell 方程的精确缺陷

定义实际原算子 residual：

$$r=b_D-L_D\widetilde h=D(C-\widetilde C)\widetilde h.$$

直接代回得

$$h-\widetilde h=L_D^{-1}r,\qquad F-\widetilde F=OL_D^{-1}r.$$

因此，低秩近似的价值不由 $\|C-\widetilde C\|$ 单独决定，而由 **材料加权、当前场激发、最终接收端口**共同决定。左端 $D$ 与右端 $\widetilde h$ 均不可丢掉。计算真实 $r$ 要付出一次 $C\widetilde h$ 的费用；这仍可能远低于重新LU，但不能记为免费。

## F-ERR-004 · 可收敛的缺陷校正与精确 fallback

把 $\widetilde L=I-D\widetilde C$ 当预条件器，

$$h_{j+1}=h_j+\widetilde L^{-1}(b_D-L_Dh_j),\qquad
 e_{j+1}=T_De_j,\quad T_D=\widetilde L^{-1}D(C-\widetilde C).$$

若有可计算上界 $\|T_D\|\le\gamma<1$，则

$$\|e_j\|\le\frac{\|h_{j+1}-h_j\|}{1-\gamma}.$$

这是由 $e_j=(I-T_D)^{-1}(h_{j+1}-h_j)$ 的Neumann级数得到。非正规系统中只看 $\rho(T_D)<1$ 不能把该范数界替换成同样常数；无法得到有用 $\gamma$ 时，可以使用预条件Krylov及真实residual，或回退完整256LU。**本章没有假设这种迭代一定更快。** 原P1规模小，直接256LU很强，只有整个求解链成本比较才能决定是否值得实施压缩。

## F-ERR-005 · 接收任务的双边误差表示

对线性观测任务 $\ell(F)=\operatorname{Re}\operatorname{tr}(W^*F)$，令

$$L_D^*Z=O^*W.$$

则严格恒等式为

$$\ell(F)-\ell(\widetilde F)=\operatorname{Re}\operatorname{tr}(Z^*r).$$

若使用近似dual $\widetilde Z$，并有 $\|L_D^{-1}\|\le\beta$，令 $r_d=O^*W-L_D^*\widetilde Z$，则

$$\left|\ell(F)-\ell(\widetilde F)-\operatorname{Re}\operatorname{tr}(\widetilde Z^*r)\right|
\le\beta\|r_d\|_F\|r\|_F.$$

证明来自 $Z-\widetilde Z=L_D^{-*}r_d$ 和Cauchy–Schwarz。dual也随材料变化；固定离线dual只可作近似或预条件。要节省成本，应只为少量实际任务求dual，不为每个材料格点建立巨大免费证书。没有稳定性上界时，该式仍是误差表示，但不是可认证误差区间。

## F-ERR-006 · 从物理误差传到 likelihood，而非追求机器精度

设白化预测误差 $\|f_w-\widetilde f_w\|\le\varepsilon$，近似 residual $\widetilde r=\widetilde f_w-d_w$。由三角不等式，

$$\tfrac12\max(0,\|\widetilde r\|-\varepsilon)^2\le\Phi\le
\tfrac12(\|\widetilde r\|+\varepsilon)^2.$$

也有

$$|\Phi-\widetilde\Phi|\le\|\widetilde r\|\varepsilon+\tfrac12\varepsilon^2.$$

似然为 $L=e^{-\Phi}$，直接反转两端指数即可得 $L$ 包络。与接受阈值或任务积分区间比较后，才能决定是否继续解物理方程。注意 $\varepsilon$ 已包含 $\sqrt2/\sigma$；原始复预测误差小不等于高SNR下likelihood误差小。

## F-ERR-007 · 固定 receiver 投影遗漏的信息

设 $S=U\Sigma V^*$，保留前 $p$ 个左奇异方向，$P_U=U_pU_p^*$。

$$\|(I-P_U)F\|_F=\|\Sigma_{tail}V_{tail}^*J_c\|_F
\le\sigma_{p+1}(S)\|J_c\|_F.$$

这是**输出截断**界，不是材料误差界。令复 residual $R=(F-Y)/\sigma\in\mathbb C^{32\times8}$，原512维实白化 residual 为 $r=\sqrt2\operatorname{pack}R$，则

$$\|r\|^2=2\|P_UR\|_F^2+2\|(I-P_U)R\|_F^2.$$

丢掉第二项就改变了likelihood，除非它与材料无关；即使尾部预测很小，尾部数据噪声与变化预测的交叉项也须处理。对于严格receiver-null current，数据的确不变，但 $J_c=D_m(E^{inc}+KJ_c)$ 仍是所有源共享的硬构成关系，不能省略。

## F-ERR-008 · 投影式 TSOM 的正确压缩对象

固定接收弱投影 $P_w$ 后，真正影响内部场的是 $KP_w$。对其作rank-$k$近似或选固定Fourier作用基，控制的是

$$\|K P_w-KP_w Q_kQ_k^*\|,$$

而不是计算两个子空间的几何交集。若先选 $V_k$ 再投影，$P_wV_k$ 通常不正交；应一次QR/SVD得到 $Q$，保留其实际rank。串联投影 $P_wP_b$ 一般不是正交投影，因为 $(P_wP_b)^*=P_bP_w$；只有两者对易时才是交集投影。使用投影复合是合法算法构造，但须明确改变了哪些方向及其metric。

## F-ERR-009 · NN 的数学位置：预测状态，不替换物理目标

如果未来用网络预测 $\widehat h_\theta(D,E_R)$ 或低维系数，其资格与上述低秩近似相同：

$$r_\theta=b_D-L_D\widehat h_\theta,\qquad
h=\widehat h_\theta+L_D^{-1}r_\theta.$$

可以先用便宜预测，再用真实residual决定继续Krylov或完整回退；也可以由 F-ERR-005/006 控制任务。网络训练误差不是 $\beta$ 或全域置信界，不能把经验误差直接代入严格证书。这样网络承担难以解析表示的状态映射，原likelihood、prior和材料支持保持明确。当前未训练此网络。

## F-ERR-010 · 计算是否值得做的代数门

设一次固定setup成本 $C_0$，原每状态成本 $c_f$，替代路径含fallback的平均成本 $c_r$。对同一固定采集的 $n$ 个实际状态：

$$C_{new}=C_0+nc_r,\qquad C_{old}=nc_f.$$

盈亏平衡 $n>C_0/(c_f-c_r)$，只在 $c_f>c_r$ 时成立。达到20%费用门要求

$$n\ge\frac{C_0}{0.8c_f-c_r},\qquad c_r<0.8c_f.$$

材料优化产生的状态数也可能变化，所以完整系统需要实际 $n_{new},n_{old}$；不能用同state核函数加速比冒称整个重建加速。R9当前精确ROI核已有对照，但完整重建和新对象泛化未运行。

## 由推导直接得到的决策

- 优先保留精确固定宿主消元；它不改变材料/概率目标，适用于完整原ROI支持。
- 均匀幅度可继续固定Schur；一般非均匀材料不可强行沿用标量可交换性。
- 固定SOM/TSOM投影用于组织current或预条件，不把非法电流当材料样本。
- 暂不新增压缩rank sweep或NN训练：只有在256LU确为剩余瓶颈，且上述缺陷校正有实际余量时，才做一项完整端到端裁决实验。
- 若做该实验，强基线须为精确256ROI solver，而不是旧576系统；固定与动态探索均共享此优化。比较完整anchor+探索+任务成本，保留失败与fallback。

**来源与新颖性：** 本章为Woodbury、residual correction、dual-weighted residual与正交投影性质在本合同下的显式推导。它们是标准工具；尚无独立Maxwell新定理主张。源码对应精确分支 `src/stable_roi_operator.py`，压缩/NN分支仅推导未实验。
