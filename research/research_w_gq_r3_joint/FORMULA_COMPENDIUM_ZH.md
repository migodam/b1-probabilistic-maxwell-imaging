# Maxwell 联合成像公式集：从 Green 积分到固定电流坐标和材料更新

本文件由 Astra 维护。目标是记录能推出实现选择的公式、推导与适用条件，不用重新命名算子冒充推导。公式编号在本轮保持稳定。结论分为已有标准公式、在明确条件下推导的界、以及待实测算法假设。公式正确性主要通过推导审阅，不反复运行机器精度恒等式测试。

来源：Chen, *Computational Methods for Electromagnetic Inverse Scattering* (2018), DOI [10.1002/9781119311997](https://doi.org/10.1002/9781119311997)，正文页150–164、233–234；具体页码与原论文见 [来源表](research/delegated/sources/FORMULA_SOURCE_LEDGER.md)。本轮源码以实际离散弱式为准；连续公式不是连续 PDE 误差证书。

## F01 — 先固定方程、维度与 Green 符号

采用时间因子 $e^{-i\omega t}$、非磁二维 $E_z$，令背景标量算子 $L_b=-\Delta-k_0^2\epsilon_b$。Green 满足 $L_bG_b=I$。在均匀背景：

$$g_2(r,r')=\frac{i}{4}H_0^{(1)}(k_b|r-r'|),\qquad
g_3(r,r')=\frac{e^{ik_b|r-r'|}}{4\pi|r-r'|}.\tag{F01}$$

二维核不能用三维 $1/r$ 核代替。三维矢量 Maxwell 还需要 dyadic Green；本轮不据标量公式声称矢量结果。

令 $\chi=\epsilon-\epsilon_b$，$V=k_0^2\chi$ 为乘法算子。由 $(L_b-V)u=f$ 得

$$u=u_b+G_bVu,\qquad u_b=G_bf.\tag{F02}$$

只要逆存在，F02 保留全部多次散射；不要求 Neumann/Born 级数收敛。

## F03 — Born 积分为什么能变成固定矩阵乘法

弱散射只在这里近似 $u\simeq u_b$。单位线源/点源的固定归一化下：

$$y_{rs}^{\rm B}=k_0^2q_s\int_D
G_b(r_r,r)\,\chi(r)\,G_b(r,r_s)\,dr.\tag{F03}$$

取单元积分权 $m_j$，定义 $[G_S]_{rj}=k_0^2G_b(r_r,r_j)m_j$、$[E_b]_{js}=u_{b,s}(r_j)$：

$$Y^{\rm B}=G_S\operatorname{diag}(\chi)E_b,\qquad
\operatorname{vec}Y^{\rm B}=(E_b^T\odot G_S)\chi.\tag{F04}$$

$\odot$ 是逐列 Kronecker 积：第 $j$ 列是 $E_b(j,:)^T\otimes G_S(:,j)$。背景与采集固定时，F04 的矩阵可一次构造；$\chi$ 是唯一变化量。精确单元积分、FEM 质量矩阵不能无声明地换成中心采样。

**可用优化：**固定 Born 初始化矩阵、FFT/Fourier 实现、一次 SVD、批量源运算。**不能推出：**完整非线性材料 Jacobian 也固定。

## F05 — 远场近似是振幅与相位的两项近似

本段远场条件采用实数 $k_b>0$ 的均匀无损背景；有损复波数需另控制相应振幅误差。令接收点 $r_r=R\hat s$，目标坐标 $|r|\leq a$。展开

$$|R\hat s-r|=R-\hat s\cdot r+
\frac{|r|^2-(\hat s\cdot r)^2}{2R}+O(a^3/R^2).\tag{F05}$$

三维振幅误差受 $a/R$ 控制，省略二次相位须 $k_ba^2/R\ll1$：

$$g_3(R\hat s,r)\simeq \frac{e^{ik_bR}}{4\pi R}
e^{-ik_b\hat s\cdot r};\quad
g_2(R\hat s,r)\simeq \frac{e^{i(k_bR+\pi/4)}}{\sqrt{8\pi k_bR}}
e^{-ik_b\hat s\cdot r}.\tag{F06}$$

二维还使用 Hankel 大参数条件 $k_bR\gg1$。平面入射 $u_b=e^{ik_b\hat d\cdot r}$ 时，Born 数据采样 $\widehat\chi(k_b(\hat s-\hat d))$；双球面源若以源位置方向定义 $\hat s_s$，入射方向是 $-\hat s_s$，对应空间频率变为 $k_b(\hat s_r+\hat s_s)$。因此“把分母提出”只处理振幅，不能省略传播相位的条件。

**本轮实现边界：**F06 可指导初始化和固定 Fourier chart；主接受函数仍使用完整 FEM。

## F07 — 固定 Green 与完整反馈可以同时成立

定义对比源 $j_s=\chi e_s$，把 $k_0^2$ 与离散积分吸收到 $G_D,G_S$：

$$e_s=e_{b,s}+G_Dj_s,\quad
j_s=\operatorname{diag}(\chi)(e_{b,s}+G_Dj_s),\quad
y_s=G_Sj_s.\tag{F07}$$

若将**整个未知墙和目标**都保留在 $\chi$ 中，参考背景可以一直取真空，$G_D,G_S$ 即使墙未知也固定。把估计墙吸收入背景是一种可选的因子化；只有采用这种因子化时背景 Green 随墙改变。不能把这两种表示混为一谈。

对 F07 求微分，共同材料方向 $h$ 在所有源下相同：

$$ (I-\operatorname{diag}(\chi)G_D)\,Dj_s[h]
=\operatorname{diag}(e_s)h,\qquad
J_s=G_S(I-\operatorname{diag}(\chi)G_D)^{-1}\operatorname{diag}(e_s).
\tag{F08}$$

推导只需对 $j=\chi(e_b+G_Dj)$ 使用乘积法则。F08 说明：固定接收谱不等于固定材料切向；内部反馈和局部场同时改变材料更新的作用。

## F09 — SOM 可观测电流坐标与投影式 TSOM

在声明的数据白化和电流质量坐标中取 $G_S=U_S\Sigma_SV_S^*$。保留 $L$ 个稳定奇异方向：

$$j_s^+=V_{S,L}\Sigma_{S,L}^{-1}U_{S,L}^*y_s,
\quad P_-=I-V_{S,L}V_{S,L}^*.\tag{F09}$$

这将独立电流分为可由稳定数据估计的部分与剩余部分。它没有证明 $P_-v$ 是某个合法共同材料变化产生的电流。

令 $V_{D,k}$ 是内部传播 $G_D$ 的保留右奇异向量。教材 Eq.6.71 的实用投影式为

$$j_s=j_s^++B_{\rm TSOM}\beta_s,\qquad
B_{\rm TSOM}=P_-V_{D,k}.\tag{F10}$$

不需要计算空间交集。$B_{\rm TSOM}$ 一般不正交、甚至可能降秩；真正投影到其列空间的是

$$\Pi_B=B(B^*B)^\dagger B^*,\tag{F11}$$

或用秩揭示 QR 得到的正交基。$P_-P_D$ 一般不是正交投影，除非相应投影交换。FFT-TSOM 使用低空间频率 Fourier **子空间**近似内部保留子空间，并非声称每根奇异向量都是一个 Fourier mode。

原 SOM/TSOM 同时约束数据与状态残差。采用教材 Eq.6.64 的非零确定电流归一化，例如

$$\sum_s\frac{\|G_Sj_s-y_s\|^2}{\|y_s\|^2}
+\sum_s\frac{\|j_s-\chi(e_{b,s}+G_Dj_s)\|^2}{\|j_s^+\|^2}.\tag{F12}$$

原算法归一化、更新顺序和截断必须按对应实现声明。只把 F10 的 basis 放入 GN，并不等于复现完整 TSOM；内部反馈项不能作为“重复 seed”删除。

## F13 — 跨照明共同材料，而不是每源各一张材料图

固定当前 $j_s,e_s$，无约束复材料的逐点最小二乘为

$$\chi_{\rm LS}(r)=
\frac{\sum_s w_s\overline{e_s(r)}j_s(r)}{\sum_s w_s|e_s(r)|^2},\tag{F13}$$

仅当分母非零；弱场点不具有稳定恢复保证。若跨频约束为 $\chi(f)=a+i\ell/f$，以同一个 $a,\ell$ 最小化所有频率源的复残差，则两列实设计向量正交，得到

$$a=\frac{\sum_{f,s}w_{fs}\Re(\bar e_{fs}j_{fs})}{\sum_{f,s}w_{fs}|e_{fs}|^2},\quad
\ell=\frac{\sum_{f,s}w_{fs}f^{-1}\Im(\bar e_{fs}j_{fs})}{\sum_{f,s}w_{fs}f^{-2}|e_{fs}|^2}.\tag{F14}$$

假设权重为正实标量、参考材料已扣除。一般相关噪声/耦合正则需解相应实正规方程。F13/F14 是条件 LS 更新，不是材料后验或唯一性证书。

## F15 — FEM 的正确材料切向含右端导数

实际离散散射场方程为 $L(z)u_{\rm sc}=b(z)$，其中 $L=L_0-V(z)$，$b=b_W+b_X$。对实参数 $z_p$：

$$q_p=\partial_pb+(\partial_pV)u_{\rm sc},\qquad
\partial_pu_{\rm sc}=L^{-1}q_p,\qquad
J_p=RL^{-1}q_p.\tag{F15}$$

由 $L_p=-V_p$ 直接微分得到。漏掉 $b_p$ 会漏掉入射场材料激发。残差加权目标 $\phi=\tfrac12\|\Sigma^{-1/2}(F-y)\|^2$，复伴随满足

$$L^*\lambda=R^*\Sigma^{-1}(F-y),\qquad
\partial_p\phi=\Re\sum_s\lambda_s^*q_{p,s}.\tag{F16}$$

每源一次加权伴随可形成标量目标梯度，不必强迫基线计算全部接收 Jacobian。完整 GN 需要更多信息，本轮明确计量两者成本。

形状变化由边界速度 $v_p\cdot n$ 给出：

$$\partial_pV_{ij}=k_0^2\alpha\int_{\partial\Omega_X}
(v_p\cdot n)\varphi_i\varphi_j\,ds,
\quad
\partial_pb_{i,s}=k_0^2\alpha\int_{\partial\Omega_X}
(v_p\cdot n)\varphi_i u_{{\rm inc},s}\,ds.\tag{F17}$$

材料幅度另有体积分导数。椭圆边界 $r=p+R_\theta\operatorname{diag}(r_x,r_y)(\cos t,\sin t)^T$；对位置、两半轴、角度分别求导，再与几何法向点乘。所有源使用相同形状速度。本轮采用相同 cut-polygon 的体/边界积分；这是有限多边形几何导数。

## F18 — Green 迁移作用于场响应，不是任意材料或 SVD 坐标

设 $L_1=L_0-\Delta V$，两者可逆。相减得到

$$G_1-G_0=G_1\Delta VG_0=G_0\Delta VG_1,
\quad G_1g=(I-G_0\Delta V)^{-1}G_0g.\tag{F18}$$

不要求 $\|G_0\Delta V\|<1$；该条件只足以保证 Neumann 展开与简单范数界。若 $\eta=\|G_0\Delta V\|<1$，则

$$\|G_1g-G_0g\|\leq\frac{\eta}{1-\eta}\|G_0g\|.\tag{F19}$$

没有稳定裕度时，代数等式仍可成立但小墙变化可能被强反馈放大。把 F18 直接乘到材料系数、右电流 SVD 向量或 posterior fiber 上没有依据，必须先说明对象属于载荷还是场空间。

## F20 — 实际 R3 chart：有限元载荷质量坐标

固定公共网格，在公共 wall∪target 区域所激发的载荷自由度上选正对角 $D$（实现为材料质量矩阵的正对角，**不是**未经检查的高阶行和 lumping）。其余自由度排除；代码中对应逆根为零。弱载荷 $q$ 写成 $g=D^{-1/2}q$。本轮独立同方差噪声以 $\sigma=.01$ 白化：

$$S_z=\sigma^{-1}RL(z)^{-1}D^{1/2},\qquad
K_z h=D^{-1/2}q[h],\qquad J_z=S_zK_z.\tag{F20}$$

多源、多频按实际顺序堆叠，实材料参数共享。F20–F22 使用复白化数据范数（与堆叠实/虚范数相同）；限制到实材料参数后仍适用。若使用参数尺度矩阵 $T$，令 $K_z$ 作用于 $h=T\bar h$，后续所有范数使用同一尺度。相同标量白化不改变参考右奇异向量；若噪声相关则必须重新白化后构建谱。取公开参考状态 $z_0$ 的 $S_0$ 右奇异基 $U_r$，$P_r=U_rU_r^*$：

$$\widetilde J_z=S_zP_rK_z,
\quad E_z=J_z-\widetilde J_z=S_z(I-P_r)K_z.\tag{F21}$$

这是 **SOM-inspired material-forcing chart**；它投影微分载荷，并通过当前完整 $L(z)$ 传播保留模式，不是原 SOM 的独立对比源优化，更不是完整 TSOM。每频一组接收伴随（16 个 RHS）构建参考 SVD；当前传播需要 $r$ 个 RHS。比较完整接收伴随的 RHS 数及总 CPU，不能因为 $r$ 小就宣称加速。

参考谱尾与背景失配可分开：

$$\|E_z\|\leq
\big[\sigma_{r+1}(S_0)+\|S_z-S_0\|\big]\,
\|(I-P_r)K_z\|,\tag{F22}$$

其中 $\|S_z-S_0\|\leq\sigma^{-1}\|R\|\|G_z\Delta VG_0D^{1/2}\|$ 可再用乘积范数上界。推导用 $S_z=S_0+(S_z-S_0)$ 及截断 SVD 性质。

F22 有三个不同来源：参考谱尾、传播改变、当前材料激发落在遗漏方向。单看 $\sigma_{r+1}$ 不足。有限方向 $v$ 的实际 $\|(J-\widetilde J)v\|$ 只检测该方向；不是全空间误差证书。

## F23 — 真正的联合 GN 与墙 nuisance 消元

下面全部为实白化、按公开参数尺度归一化的坐标。令 $J=[A,B]$，其中 $A$ 为墙、$B$ 为目标。对固定、墙与目标分块独立的二次 prior $\Lambda=\operatorname{diag}(\Lambda_W,\Lambda_X)$ 的局部模型：

$$H=J^TJ+\Lambda,\qquad g=J^Tr+\Lambda(z-\mu),
\qquad H\delta=-g.\tag{F23}$$

若 prior 非二次，使用其实际局部曲率并保留目标函数。数值 LM damping 可加入求步矩阵，但不属于 posterior precision。相关 prior 的交叉精度必须加进 $C$，以下 F24–F27 的简式不直接适用于相关 prior。分块 $H_w=A^TA+\Lambda_W$、$C=A^TB$、$H_x=B^TB+\Lambda_X$：
$$S_X=H_x-C^TH_w^{-1}C,\quad
\delta x=-S_X^{-1}(g_x-C^TH_w^{-1}g_w),\quad
\delta w=-H_w^{-1}(g_w+C\delta x).\tag{F24}$$

这是标准 Schur 消元。带 box 的问题须在共同 active face 上约化，或者直接解联合有界最小二乘；不能无视活跃约束套用 F24。

若 $\Lambda_W\succ0$，定义 $\bar A=A\Lambda_W^{-1/2}$，则由 Woodbury：

$$Q_A=I-A(A^TA+\Lambda_W)^{-1}A^T
=(I+\bar A\bar A^T)^{-1},\quad S_X=B^TQ_AB+\Lambda_X.\tag{F25}$$

因此 $0\prec Q_A\preceq I$。该对象是 proper prior 下的正则消元；不再等于裸正交投影，也不能把所有广义特征值都叫 $\sin^2\theta$。

## F26 — A17-X 基础上的目标 Schur 专门化：chart 误差的更新后果

固定相同 prior，令 $\widehat A=A+E_A$、$\widehat B=B+E_B$，$\bar E_A=E_A\Lambda_W^{-1/2}$。逆矩阵差恒等式给出

$$\|\widehat Q_A-Q_A\|
\leq (2\|\bar A\|+\|\bar E_A\|)\|\bar E_A\|.\tag{F26}$$

证明：$\widehat Q-Q=\widehat Q(\bar A\bar A^T-\widehat{\bar A}\widehat{\bar A}^T)Q$，两侧逆的范数均不超过1，再展开中项。进而

$$\|\widehat S_X-S_X\|\leq
2\|B\|\|E_B\|+\|E_B\|^2+
\|B\|^2(2\|\bar A\|+\|\bar E_A\|)\|\bar E_A\|.
\tag{F27}$$

证明：先在同一 $\widehat Q$ 下展开 $\widehat B^T\widehat Q\widehat B-B^T\widehat QB$，再加 $B^T(\widehat Q-Q)B$。此界不需要变化中的奇异子空间有固定谱隙，但依赖墙 prior 的正定尺度；若 prior 退化，不能直接套用。

**状态：PROVED_UNDER_ASSUMPTIONS。** 推导由标准逆差公式得到，单独不是已证明的新论文贡献。其研究用途是把 F21 的材料切向遗漏与目标 Schur 系统联系起来，检验哪些投影错误会改变目标更新。

## F28 — 不必先形成完整 Jacobian 的更新后果检查

chart 求得联合步 $\widehat\delta$。计算它在完整局部二次模型中的残差

$$e=H\widehat\delta+g
=J^T(J\widehat\delta+r)+\Lambda(z+\widehat\delta-\mu).
\tag{F28}$$

一次 native JVP 与一次加权 VJP 可求 F28，不必形成全部 $J$；实际多源 RHS 仍全部计费。若完整解位于相同无约束空间/共同活跃面，以下块全部在共同自由变量上约化，且 $H_w\succ0,S_X\succ0$：

$$\widehat\delta x-\delta x
=S_X^{-1}(e_x-C^TH_w^{-1}e_w),\tag{F29}$$

$$\|\widehat\delta x-\delta x\|
\leq\frac{\|e_x\|+\|C^TH_w^{-1}\|\|e_w\|}
{\lambda_{\min}(S_X)}.\tag{F30}$$

符号由 $H(\widehat\delta-\delta)=e$ 得到。目标弱信息使很小的切向误差仍可能造成很大的目标步变化。实际采样 JVP 误差不能提供未计算的 $S_X$ 下界；没有廉价下界时应只称经验触发或离线审计。

以下比较共享相同的当前完整残差、prior、damping 与可行步集合；如果前向残差也近似，须另加残差误差项。另一完全显式界：若 $H\succeq\mu_H I$、相同凸可行步集合，$\delta,\widehat\delta$ 分别为完整/近似 QP 的精确最小点，则变分不等式相加得

$$\|\widehat\delta-\delta\|
\leq\frac{\|(H-\widehat H)\widehat\delta+g-\widehat g\|}{\mu_H},
\quad
\|H-\widehat H\|\leq2\|J\|\|E\|+\|E\|^2,
\quad\|g-\widehat g\|\leq\|E\|\|r\|.\tag{F31}$$

F31 允许活跃面变化，但只给联合步界；非精确 QP 还要加优化残差。若目标无 proper prior，$\mu_H$ 可能只来自算法 damping，结论只适用于该阻尼更新而非后验不确定性。

## F32 — 为什么全波接受不能替代材料正确性

候选步通过真实完整目标函数的下降与 trust-region 比率接受：

$$\rho=\frac{\mathcal L(z)-\mathcal L(z+\delta)}
{-g^T\delta-\tfrac12\delta^TH\delta}.\tag{F32}$$

它能防止近似 chart 导致接受一个上升步，却不能证明材料/目标更接近真值；联合补偿和错误局部极小仍存在。真值只能在输出封存后评分。多个分支按相同完整目标函数选择，不能按目标误差选分支。

## F33 — 有限材料补偿代价与真实任务分辨率

对指定目标方向 $v$，合法墙步集合 $\mathcal C_W$ 下的线性残余：

$$c(v)=\min_{h\in\mathcal C_W}\|Ah+Bv\|^2.\tag{F33}$$

无约束、加二次墙惩罚时最优 $h=-(A^TA+\Lambda_W)^{-1}A^TBv$，最小总惩罚为 $v^TB^TQ_ABv$。有限 box 则需有界最小二乘，不能用伪逆的无限幅度补偿代替物理合法性。

若 $h$ 实际达到 F33 的最小值（否则用 $\|Ah+Bv\|$ 替代 $\sqrt{c(v)}$），$t\geq0$ 且整条线段合法，沿 $(h,v)$ 存在统一二阶余项界 $\|F(z+t(h,v))-F(z)-tJ(h,v)\|\leq Mt^2/2$，则

$$\|F(z+t(h,v))-F(z)\|\leq t\sqrt{c(v)}+Mt^2/2.\tag{F34}$$

只有实际合法端点给出真实有限歧义见证；局部小奇异值与随机抽样曲率最大值均不自动提供全局结论。本轮优先成像和关键路径对照，不重启大规模补偿 campaign。

## F35 — 哪些 chart 错误真正改变目标信息？

F27 的全算子范数界可能很保守。还可以直接对“墙体已作最佳局部补偿”的目标信息求导。固定正定 $\Lambda_W$、固定 $\Lambda_X$，设

$$N=(A^TA+\Lambda_W)^{-1}A^TB,\qquad
R=B-AN=Q_AB. \tag{F35}$$

这里 $-Nv$ 不是人为锁定的弱坐标，而是目标方向 $v$ 对应的最优正则墙体补偿。由二次最小化，

$$v^TS_Xv=\min_h\{\|Ah+Bv\|^2+h^T\Lambda_Wh+v^T\Lambda_Xv\}.
$$

最优点满足 $h_*=-Nv$。对 $A,B$ 的光滑变化使用包络定理，最优 $h_*$ 的一阶变化项因驻点条件抵消，因此

$$dS_X=(dB-dA\,N)^TR+R^T(dB-dA\,N),
\qquad
v^T(dS_X)v=2\langle Rv,(dB-dA\,N)v\rangle. \tag{F36}$$

**推导后果：**目标信息的一阶误差由“目标切向误差减去配套墙补偿的切向误差”决定，而非单独由 current 奇异值、$\|dA\|$ 或 $\|dB\|$ 决定。两种很大的误差可以在这个特定目标方向上抵消，也可以共同放大。对 chart 的有限误差，可用 $dA=E_A,dB=E_B$ 得到一阶预测；实际有限 Schur 差还包含高阶项，不能把此预测当成严格误差上界。

F36 是标准二次消元与包络求导的明确推论，不单独声称原创。本轮使用已付费 Jacobian 进行离线归因；在线若为获得 $N,R$ 又计算完整 Jacobian，费用不能隐去。它提示一种更有依据的后续模式优先级：检查合法共同材料方向 $(h_*,v)$ 的实际载荷误差，再按它与剩余目标数据方向 $Rv$ 的作用排序，而非把所有接收弱方向锁死。当前固定主方法不会据此事后改动盲测。

## F37 — 墙变化与目标自身变化必须分开归因

令 $G_{ij}=G(W_i,X_j)$，$\Delta V_W=V(W_1)-V(W_0)$，$\Delta V_X=V(X_1)-V(X_0)$，所有逆均存在。沿一个明确路径使用 F18：

$$G_{11}-G_{00}
=G_{11}\Delta V_XG_{10}
+G_{10}\Delta V_WG_{00}. \tag{F37}$$

第一项改变目标、固定 $W_1$，第二项改变墙、固定 $X_0$。换路径会改变两项归因，完整差不变；墙–目标反馈使它们不是两个独立可相加的材料主效应。故本轮安排公共非真值的墙/目标 $2\times2$ 状态对照，保留 $G_{11}-G_{10}-G_{01}+G_{00}$ 的交互含义。仅在“参考无目标→当前强目标”处发现固定 chart 失配，不能把全部失配归因于未知墙。

## F38 — 终止前真正需要检查的是加权伴随梯度

当前完整前向残差相同、prior相同，令 $E=J-\widetilde J$。那么

$$g-\widetilde g=E^Tr,\qquad
\mathcal G_\eta(z;g)=\frac{z-\Pi_{\mathcal C}(z-\eta g)}{\eta},\qquad
\|\mathcal G_\eta(z;g)-\mathcal G_\eta(z;\widetilde g)\|
\leq\|E^Tr\|.\tag{F38}$$

证明只用同一闭凸box的欧氏投影非扩张性；$\eta>0$且坐标尺度固定。故小 $\|E\delta\|$ 不控制完整驻点，甚至 $E\delta=0$ 时也可有 $E^Tr\ne0$。本轮adaptive正是只按候选方向审核，末端完整梯度的失败不能被前向下降或近似gradient stop掩盖。

对于FEM的复响应导数 $J_p^c=R L^{-1}q_p$，定义原始复残差 $e_s=F_s-y_s$；本轮独立等方差实/虚噪声下，以每源加权接收残差求 $L^*\lambda_s=R^*(e_s/\sigma^2)$，再取 $\operatorname{Re}\sum_s\lambda_s^*q_{p,s}$，可得到全部实材料坐标的完整数据梯度。这里不把实白化向量$r$直接当作复接收负载；一般相关噪声应先按实化协方差做权重，再映回共轭伴随负载。$q_p$仍含材料依赖右端。当前4源×3频需要12个加权adjoint RHS；完整接收Jacobian用16接收×3频需要48个adjoint RHS。**这是同一已缓存状态下的RHS比较，实际重分解、组装、shape积分、I/O仍需计费，并非实测端到端加速。**

**实际终止量边界：**当前冻结程序把投影梯度逐列除以 $\max(\sqrt{\operatorname{diag}(J_{\rm scaled}^TJ_{\rm scaled})_j},1)$ 后才使用0.005门。12个加权伴随RHS只给完整梯度，不免费给出全部Jacobian列范数。因而不能直接宣称12RHS无改动复现本轮精确驻点门；未来若用固定公开尺度的gradient mapping，需要独立冻结其门，或仍计费取得原曲率对角。F38的非扩张不等式针对明示的欧氏投影映射，不把当前数据依赖对角量隐去。

这条标准投影/伴随结论给出具体的下一步：用可支付的完整加权梯度决定是否真的终止，再讨论如何补足造成目标后果的chart方向。它不是新的A17X之外的原创定理；本轮科学方法已经冻结，没有依据盲测事后替换终止策略，也没有声称该改动已经改善成像。

## F39 — k-space 中面积、材料与形状为什么可能相互补偿

对连续椭圆 $\Omega=p+R_\theta\operatorname{diag}(a,b)\mathbb D$，把积分变量换到单位圆并作角向积分，得到精确的指示函数 Fourier 变换：

$$\widehat{1_\Omega}(q)=\pi ab\,e^{-iq\cdot p}\frac{2J_1(Q)}{Q},\quad
Q^2=q^TR_\theta\operatorname{diag}(a^2,b^2)R_\theta^Tq.\tag{F39}$$

$Q=0$取连续极限。用 $2J_1(Q)/Q=1-Q^2/8+O(Q^4)$ 可知，在适用的 Born/低空间频率条件下，首先被约束的是面积乘复对比度 $\pi ab(c+i\ell/f)$，位置进入相位，形状二阶矩进入下一阶。路径 $a,b\mapsto e^t(a,b)$、$c,\ell\mapsto e^{-2t}(c,\ell)$ 精确保留积分幅度，但不保留二阶形状信息，也不保留一般完整 Maxwell 响应。这是具体的可检验耦合，不是全波精确 gauge。当前近场和高对比度未必满足低 $Q$；不能用此展开替代已有完整 Jacobian 或宣称形状失败原因已经证明。

可用 $u=\log(A/A_*)$、$v=\tfrac12\log(a/b)$、$m_R=Ac$、$m_I=A\ell$ 重新表达同一合法模型。改变信赖域坐标可能改善有限预算路径，但原参数box在新坐标中是耦合约束，不能偷换先验与可行集。只用Schur精确消去同一个QP的幅度变量，不会改变该QP的步。

## F40 — FFT加速固定Green作用，不等于一次FFT求完非均匀全波

同一规则格点上的平移不变参考Green离散可写为卷积。适当零填充并正确处理自项后：

$$G_0j=\operatorname{crop}\,\mathcal F^{-1}\!\left[\widehat g_0\odot\mathcal F(\operatorname{pad}j)\right],\qquad
j=D(e^{inc}+G_0j).\tag{F40}$$

第一式是线性卷积的实现恒等式；第二式仍有完整材料反馈。非均匀材料乘法 $D$ 通常不与 $G_0$ 交换，因此 $(I-DG_0)^{-1}$ 不能一般化为一次频域逐点除法。FFT可服务于矩阵作用、Krylov或预条件，也可以加速近场Born初值。当前4个TX与有限16RX不支持把补零FFT点当作新增测量；P3非结构FEM、材料弱式及PML也不能无声明地换成周期卷积。

## F41 — 固定形状的幅度有理族：完整反馈与材料依赖右端同时保留

固定墙和目标形状，令 $L_\delta=L_W-\delta M_S$、$b_\delta=b_W+\delta b_S$，$G_W=L_W^{-1}$、$u_W=G_Wb_W$。从两套散射场方程相减：

$$u_\delta-u_W
=\delta(I-\delta G_WM_S)^{-1}G_W(b_S+M_Su_W).\tag{F41}$$

这里 $M_S$ 是当前FEM形状的非对角质量块，不是随意的rank1材料参数。设 $B=G_WM_S$、$a=G_W(b_S+M_Su_W)$，则输出增量为 $\delta R(I-\delta B)^{-1}a$，幅度导数为 $R(I-\delta B)^{-2}a$，因为 $B$ 与自身resolvent交换。跨频 $\delta_f=c+i\ell/f$ 仍只有两个共同实材料参数。该条件有理族可预付setup后多次profile幅度；改变墙或形状通常使setup失效。它继承标准resolvent/T-matrix及相邻R10结构，不是本轮新理论或已实测加速。

三条路线的具体接入条件、近场角谱、有限孔径与噪声边界见 [Astra FFT/k-space审阅](research/delegated/theory_breakthrough/FFT_KSPACE_REVIEW.md)。本轮已冻结的盲测没有因此事后换方法；这些推导直接指导下一步计算，而不另跑机器精度公式测试。

## 与已有 A17-X 的继承关系

用户提醒后的归属核对确认：A17-X §1.1 已明确共同材料激发—内部反馈—接收映射，§1.3 已给残差相关的 paired GN 更新后果及活跃约束边界，§3.1 已有有限 resolvent 公式，§9.1–9.3 已有 Fisher 缺陷、材料移动附加误差与提案加权 action 选择。因此本文件相关式子属于**已有结果复用、标准数学或向当前 FEM/墙–目标变量的专门化**，不作为本轮新理论。F35/F36 的墙 nuisance 消元表达是标准包络与 Schur 求导；即使在所查 A17-X 文本中没有同样写法，也不能据此声称论文原创性。逐项对应见 [A17-X 归属表](research/delegated/sources/A17X_REUSE_CROSSWALK.md)。本轮可能新增的只有独立成像及方法对照证据，是否有增量由真实结果决定。

## 公式到实现的取舍

| 关系 | 可直接采用的优化 | 本轮限制 |
|---|---|---|
| F03–F06 | 固定 Born/Fourier 初值与批量矩阵操作 | 不替换 full-wave 目标函数 |
| F07–F10 | 固定公共参考 current chart | 必须更新局部场/材料激发，TSOM 第二层不能删除 |
| F15–F17 | 复用 LU，接收伴随或加权目标伴随，连续几何导数 | 保留右端导数与公共材料参数 |
| F18–F22 | 固定载荷 chart、当前场传播、失配触发修复 | 回退完整 Jacobian，所有 action 计费 |
| F24–F31 | nuisance-profiled 步及目标后果审计 | 标准 Schur；无谱下界不得假称证书 |
| F32–F34 | 全波接受、多分支保留、有限补偿解释 | 残差下降不替代独立目标成像质量 |
| F35–F36 | 依据补偿后的目标切向误差安排模式审计 | 一阶、局部；本轮离线分析，不宣称廉价在线证书 |
| F37 | 区分墙变化、目标变化及交互反馈 | 用固定公共状态作归因，不根据 blind 真值选例 |

后续版本将把真实代码模块、实测成本、有效/失败场景链接到公式。若新数学没有改变实现选择，就不为它额外增加实验。
