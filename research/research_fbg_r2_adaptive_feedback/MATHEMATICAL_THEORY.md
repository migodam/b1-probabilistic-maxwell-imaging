# FBG-R2：固定反馈如何改变辅助几何

本文区分三件事：真实 Maxwell 方程、同一方程的等价表示，以及用于搜索的另一个残差度量。反馈不增加测量，也不改变合法物理解；它可能改变优化轨迹，是否改善材料恢复由实验决定。

## F01　精确反馈表示（标准代数恒等式）

假设固定离散 Green 矩阵 G、接收 S、入射 E；实对角材料 X=diag(m)，L=I−XG 可逆。源为列，j=L⁻¹XE，Q=E+Gj，F=Sj。固定 B 且 I+XB 可逆时，

$$
\widetilde X=(I+XB)^{-1}X=X(I+BX)^{-1},\quad
\widetilde G=G+B,\quad \widetilde Q=E+(G+B)j=(I+BX)Q.
$$

因为 I−X̃G̃=(I+XB)⁻¹L，所以变换后在**原入射 E**下的电流仍为 j。Q̃ 是总场，不能把它再当作新入射场代入。非对角 B 通常使 X̃ 非局部；本轮始终优化原来的 576 个独立实材料，绝不把 X̃ 自由当作材料。

只改变表示并同步运输残差度量时，是 exact-equivalent no-op。若使用 Q̃ 的 Euclidean Gram，则是主动改变辅助目标，不能再说原目标被等价凸化。

更明确地，Mf=I+BX、V=MfQ，运输后的场度量 C=Mf^{−*}Mf^{−1} 满足 V*CV=Q*Q，故相应Gram目标逐点严格相同。C随材料的导数也要一起运输，不能将它冻结后仍称no-op。本轮证明这个控制并在Stage0验证其状态层前提，没有为同一目标再跑一份重复优化。

## F02　固定反馈切向的消去（本轮关键算法连接）

材料方向 h 的完整电流切向 T 满足

$$
LT=X_hQ,\qquad DQ[h]=GT,\qquad DF[h]=ST.
$$

固定 B 时，

$$
D\widetilde Q[h]=BX_hQ+(I+BX)GT
=B(I-XG)T+(I+BX)GT
=\boxed{(G+B)T}.
$$

这是标准切向方程推导出的消去，不是新的 Maxwell 定理。算法价值在于：拟合 (G+B)T 恰好拟合了变换后总场的方向变化，同时将完整反馈梯度压缩为一次原物理伴随。若外层更新 B，内层必须冻结；若每步令 B=B(m)，还要加入 (DB[h])j，本轮不使用后一种方式。

## F03　反馈目标与单伴随梯度

设 V=Q̃ 满列秩，H=V*V，R=F−Y，固定正数 c。定义

$$
f_B=\tfrac12c\operatorname{tr}(RH^{-1}R^*),\quad
W=RH^{-1},\quad D=W^*W.
$$

对实材料求导：

$$
Df_B[h]=c\Re\{\langle W,ST\rangle_F-
\langle VD,(G+B)T\rangle_F\}.
$$

令

$$
a=L^{-*}[S^*W-(G+B)^*VD],\qquad
(\nabla f_B)_i=c\Re\sum_s\overline{a_{is}}Q_{is}.
$$

仅需共享 LU、多源正演和一次多源伴随，再加 p×p Gram；没有逐材料有限差分。代码 Cholesky 求解 H，不显式求逆。固定参考归一化和普通材料正则另行加入；不能让 c 随当前 m 改变后仍用此梯度。

## F04　源归一化的电流拟合及其边界

令 H₀=Q*Q 正定，Zₗ=TₗH₀⁻¹ᐟ²。源基变换 E→EC、Y→YC，C 可逆时，所有 Zₗ 同时右乘同一酉矩阵。因此用同一组方向、相同 ridge 定义的 Frobenius 回归不依赖任意源坐标。

对角反馈最小化 Σₗwₗ‖(G+B)Zₗ‖²+λ‖b−bref‖²，有

$$
b_i={\lambda b_{ref,i}-\sum_{\ell,s}w_\ell\bar Z_{\ell,is}(GZ_\ell)_{is}
\over\lambda+\sum_{\ell,s}w_\ell|Z_{\ell,is}|^2}.
$$

实约束取相应实部的正规方程；b≥0 则投影到非负半轴。不是把任意复数解取实部后保留原复数最优性。

局部块将方向/源拼接为 Zb、Vb=(GZ)b，有

$$
B_b=(\lambda B_{ref,b}-V_bZ_b^*)
(Z_bZ_b^*+\lambda I)^{-1}.
$$

这是已知 ridge 最小二乘闭式解。新增的研究问题是这种**精确变换后场切向衰减**能否预测有用的材料搜索，而不是回归本身是否新颖。生产实现使用 2×2 空间块，即 4×4 矩阵，批量求解。

它不最小化真实数据曲率，也不最小化投影曲率。即使 T 精确，压低变换后场的变化也可能压低有用的信息。尤其 B=−G 完全消去这一变化，却没有把材料反演变成线性问题。

## F05　投影导数与 A17 尾项

对满列秩 V，P=VV†，固定 B，

$$
DP[h]=(I-P)(G+B)TV^\dagger+\mathrm{h.c.}
$$

A17 对指定 retained/candidate Galerkin 电流空间给出条件 Schur 增量 ΔT。使用**共同精确 anchor V**时，

$$
\Delta DP[h]=(I-P)(G+B)\Delta T V^\dagger+\mathrm{h.c.}
$$

两个 Galerkin 近似之间此式精确；它不是 Full-J 误差证书。若 T̂ 与完整切向的残差 r=XhQ−LT̂，则

$$
T-\widehat T=L^{-1}r,\quad
\|DP-D\widehat P\|_F\le
2\|(I-P)(G+B)L^{-1}r\|_F\|V^\dagger\|_2.
$$

这给出 Schur 近似怎样进入 illumination-projector 误差的明确通路。尚不能据此推出跨盆地成功；尾项须支付完整作用或其他可靠界。本轮主成像如使用完整切向，其收益不能归给 A17 Schur。src/schur_current.py 仅移植 scalar-TM 电流空间代数，不移植 A17 的向量 DDA 材料约定。

## F06　两个反例与尺度防护

**反例一：** B=−G，V=E、DV=0，而 F=S(I−XG)⁻¹XE 仍为非线性。照明投影不转动不是成像成功证书。

**反例二：** 正材料下 B=tI、t→∞，V=Q+t j。若 j 满列秩，固定非零 R 的 raw ψB=O(t⁻²)，原始物理残差完全不变。仅按 ψB(m₀) 选 B 是退化设计。

生产使用公共 mref=.1 的 sB=‖(I+BXref)Qref‖²/‖Qref‖²，再乘 ψB。它消除整体常数增益，不消除所有材料相关条件数投机。最终在线选择始终是原始数据损失加同一正则；候选 pilot 也直接比较该物理分数。

## F07　可逆性与限定 box

材料 box [.02,20]。对角 B 的 |1+bi mi| 在实区间上的极小值可直接求出，最大值在端点。块反馈一次缩放，使 maxb‖Bb‖₂·20≤.8，则对 box 内所有材料，

$$
\sigma_{min}(I+BX)\ge .2,\quad
\sigma_{max}(I+BX)\le1.8,\quad\kappa_2\le9.
$$

同一范数条件也适用于 I+XB。对角情况两者相同。此界很保守，可能限制反馈效用，应在结果中报告实际缩放。谱半径小不能替代奇异值条件，非正规矩阵的反例保留在 Stage 0。

## F08　信息、能量与不可声称的结论

H 正定时 fB=0 ⇔ R=0。精确根附近、没有正则时，反馈加权 Hessian 与原始 Jacobian 的零空间一致；反馈没有增加可辨识信息。远离根时 Gram 的材料导数可改变梯度与吸引域，但没有全局可达保证。

可逆非酉源混合保持 raw ψB 不变，却通常改变 ‖RC‖²、重新计算的 R4 归一化与正则相对权重。源混合审计固定原参考常数并运输噪声度量，不能用新归一化的不同目标冒充不变性。

清洁数据下，所有 E/Y 同乘 α，归一化数据目标、归一化反馈目标和材料正则相同，因此两种整体照明能量协议不必重复跑确定性反演。噪声绝对方差固定时此等价一般失效；该噪声实验须另行测量，不能从 clean 推出。

Stage 0 是小矩阵代数审计，数值结果在 audit/STAGE0_MATH.json；它不构成 576 维成像证据。当前没有全域收敛、全根遍历、完整概率解集或自动 GT 选择定理。

## F09　均匀参考下的残差驱动对比度项（本轮条件命题）

**结论。** 在公共均匀参考 m=b·1、b>0，取实标量反馈 B=βI、a=1+βb>0，并按同一参考使用 sβ=a²。设 f₀=(c/2)tr[R(Q*Q)⁻¹R*]，c>0 在该材料微分中固定。W=R(Q*Q)⁻¹，

$$
d_i=(QW^*WQ^*)_{ii}=\|WQ_{i,:}^*\|_2^2\ge0.
$$

则归一化反馈目标在该参考点满足

$$
\boxed{\nabla\widehat f_\beta(b\mathbf1)
=\nabla f_0(b\mathbf1)-c\frac{\beta}{1+\beta b}d.}
$$

**假设。** 原 Maxwell resolvent 可逆、Q 满列秩、源/接收及数据固定、材料方向为独立实像素。反馈与归一化常数在材料微分时冻结。公共正则相同，故梯度差中抵消。若参考缩放来自另一个均匀材料，式右侧另有一个正的全局因子，但分量符号阈值不变。

**证明。** 均匀材料使 L*=I−bG*，反馈场为 aQ。缩放后的伴随 RHS 可写成 S*W−(G*+βI)QD/a。与原 RHS 比较，

$$
\frac{G^*+\beta I}{a}-G^*
=\frac{\beta}{a}(I-bG^*)=\frac{\beta}{a}L^*.
$$

原物理伴随解消去 L*，恰得上式。特别地，两个目标在公共参考**值相同**，但梯度相差非负空间场 d 的负倍数。因此该效应不是整体缩放目标值。正 β 会在下降方向中加入非负的、由当前观测残差和总场共同决定的对比度提升项。

**算法后果。** 可以在公共 b=.1 的一个已缓存完整场上，计算 g₀/c 与 d，无 GT、无额外入射数据。对 g₀,i>0、dᵢ>0，令 tᵢ=g₀,i/(c dᵢ)。若 0<tᵢ<1/b，

$$
\beta_i=\frac{t_i}{1-bt_i}
$$

是该像素梯度变号的阈值；β=βi 时为零，严格变负需 β>βi。Stage 2 开发规则取 d 高于其全图中位数的合法正阈值之中位数，截到 [.05,.4]；空集合取 0。每个目标仅从公共参考校准一次，内层固定 B。此规则的有效性必须另做盲测，不能由梯度符号定理推出。

**边界。** d 不是 GT 支撑概率。对比度提升可能加重过估计，也可能毫无帮助。局部梯度变号不保证非局部盆地可达、真实材料恢复或更低物理费用。该条件命题及其校准器不使用 A17 Schur；若最终有效，应归因于反馈残差几何而非 Galerkin 加速。

**源基不变性推论。** 对同一入射集合的可逆混合 C，有 Q′=QC、W′=WC^{−*}，因此 Q′W′*W′Q′*=QW*WQ*。raw Gram 梯度本身也不变，故 g₀/c、d 及其变号阈值、eligible 集合和中位数校准均不依赖此源坐标。这不消除归一化物理排名和噪声协方差的运输要求，也不把 24 个实际新源等同于 12 源的可逆混合。


## F10　粗尺度、latent 与 characteristic subspace

对任意可微材料表示 m=decoder(z)，条件F09在一个映射到均匀参考的 z 上经链式法则变成

$$
\nabla_z\widehat f_\beta-\nabla_zf_0
=-c\frac{\beta}{1+\beta b}(D_zm)^T d.
$$

本轮双线性粗到细插值的权重非负，所以粗系数也得到非负的下降偏置；这不是把576阶物理矩阵降到64阶。一般NN decoder或SVD characteristic basis列可正可负，投影后的latent分量不再全部同号，不能把像素空间的“对比度提升”直接搬到latent上。还要检查GT表示能力、材料box和decoderJacobian。

## F11　真正带符号的辅助方向曲率

固定B和实材料方向h，T₁=L⁻¹XhQ，T₂=2L⁻¹XhGT₁。令V=Q̃，V₁=(G+B)T₁，V₂=(G+B)T₂，R₁=ST₁，R₂=ST₂；

$$
H_1=V_1^*V+V^*V_1,\quad
H_2=V_2^*V+V^*V_2+2V_1^*V_1.
$$

K=H⁻¹满足K₁=−KH₁K、K₂=2KH₁KH₁K−KH₂K。由乘积法则，raw ψB 的二阶导为

$$
\psi_B''=2\Re\operatorname{tr}(R_2KR^*)
+2\operatorname{tr}(R_1KR_1^*)
+4\Re\operatorname{tr}(R_1K_1R^*)
+\operatorname{tr}(RK_2R^*).
$$

除第二项外都可能带符号；最后两项是动态Gram的贡献，不能被非负的projector旋转范数替代。src/objective_geometry.py通过原LU切向和小型Cholesky作用直接计算这些量，无有限差分Hessian。本轮使用一条公共数据梯度方向作开发机制诊断，未把它称为全材料Hessian谱或盆地边界。


## F12　与原 GlobalDual 可观测残差的连接

F09 的非负场无需重新定义一个没有计算意义的算子。令已有最小范数残差 lift 为 C_R=RQ†，则

$$
d_i=\|(C_R)_{:,i}\|_2^2,\qquad
\sum_i d_i=\|RQ^\dagger\|_F^2=\psi_0.
$$

如果沿用 R1 数据消元的 W₀E=Y、C_Y=(S+W₀G)X−W₀，则 C_YQ=R，因此 C_R=C_YP_Q。反馈的对比度驱动恰好来自原有 illumination-observable compatibility residual 的逐单元列能量，而不是与观测无关的任意材料抬升。

沿全1材料方向，两个梯度之差满足

$$
\mathbf1^T(\nabla\widehat f_\beta-\nabla f_0)
=-c\frac{\beta}{1+\beta b}\psi_0.
$$

在同一均匀参考固定物理场，若 Y=Fref+εZ，则 d=O(ε²)，该梯度修正**精确二阶**激活。它不会在任意弱数据扰动下自动成为主要收益；这与低/中对比度普通Full-J已可靠时反馈不保证更优相容，但不是把任意材料对比度A直接等同数据ε的Born近似。

在均匀真根且参考一致、无额外正则时，R=0，归一化反馈Gram与原Gram对一阶数据扰动的加权Hessian相同。因而该条件下的差异主要来自远离拟合根的高阶残差几何，不能仅描述为一个更好的根附近Jacobian条件数。非均匀真材料不享有这个标量抵消，需重新分析。

C_R 的列能量不是概率密度，也不代表GT支撑。它可因错误残差方向或模型失配而激发错误位置，所有成像结论仍需真正的新鲜数据与原物理验收。
