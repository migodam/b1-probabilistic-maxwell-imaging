# 一轮实际修订：从全端口表示转向提案条件的有限RHS

## 原假设与失败证据

原端口路线希望复用低维P=R0 U、H/O进行大量有限提案。17份F3缓存显示DeltaX rank575，七端口截断误差大。因此“低latent即低端口”不能作为默认计算假设，全部575列setup也未必比Full-J便宜。

## 本轮单一机制改变

**不改材料proposal、不改prior/noise，只改有限响应的内部计算表示。** 对一个完整合法材料增量，直接解

\[
L_D k=b_D,\qquad L_D=L_0-\Delta XG,\quad b_D=\Delta X E_0.
\]

同一个current空间V中求k_V；增广C时对L_D使用原Galerkin Schur更新。RHS仅有实际source列数，不先构造s列R0 U或s×s H。剩余误差由rho=b_D−L_D k_V经过full/dual控制，普通goal-ROM使用同样的current pool和复用权利。

这沿用已有Galerkin恒等式，不是新散射定理，也不是证明source数小就一定便宜：basis、L_D作用、预条件、刷新与尾界仍会付费。它也不提供新的材料可实现性；材料先由原prior产生，评价器始终固定材料，避免current比值离开支持或破坏proposal measure。

## 冻结合成对照与结果

seed20261019，12个dense full-support合成增量，矩阵N18、3sources、4outputs；保持相同V/C，与full、exact T-matrix及ordinary同空间Galerkin比较。原结果未覆盖，保存为checks/SOURCE_RHS_REVISION.json及带SYNTHETIC标签的数组。

- exact-port与full最大相对差6.33e−16。
- 直接finite-state exchange与ordinary同空间Galerkin最大差4.61e−16。
- 浮点残差输出界12/12覆盖，但不标certified，未验证真实后端速度。
- 新增physical states/RHS/LU均0；本次合成约.0340s。

**决定：保留此实现作为未来固定提案求值竞争默认臂；不主张A17增量。** 若ordinary goal-ROM相同或更低成本，就保留普通方法。无P1/新预算不进入真实测试，亦不根据合成结果训练NN。
