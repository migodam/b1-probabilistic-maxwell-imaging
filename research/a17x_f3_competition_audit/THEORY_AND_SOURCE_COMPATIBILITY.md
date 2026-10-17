# 理论独立复核与新增可证伪诊断

## 保留的条件恒等式

对精确支持DeltaX=U D U*、固定full anchor，T=D(B+HT)、delta j=P T、delta F=O T及Q=(I−DH)^{-1}D成立。本构增量必须是真正DeltaX；不能把3D DDA的极化率增量换成任意chi差。

原A17对端口U作Galerkin增广，得到P_new−P_V=barC Gamma^{-1}N_U；再用resolvent两次修饰，有限差为

\[
\Delta F_{new}-\Delta F_V=(S+O_VQ_VU^*G)\bar C\Gamma^{-1}N_UQ_{new}B.
\]

推导正确，但左边仅是指定两个current模型之间的差。full−new尾项没有消失。R_U=U−L0 P_V时，有限状态残差rho_D=R_U Q_V B成立。若有效eps_H/eps_O和theta=||Q_V||eps_H<1，包内尾界成立；浮点SVD所得常数仍是NUMERICAL_BOUND_ESTIMATE，不是区间认证。

独立合成反例theta=2，而完整feedback矩阵为−1、最小奇异值1；此时界不可用，却并不奇异。不能将bound-unavailable解释为材料不合法。

## 相同current空间下不能凭恒等式获得精度增益

普通Galerkin ROM若采用同一个[V,C]，其P_direct与A17 Schur更新P_V+barC Gamma^{-1}N_U完全相同。因此对任意共同D，两者有限输出也相同。独立合成检查差为2.75e−16。

这不是否定A17的计算潜力：它可通过更便宜的增量更新、不同的付费action选择或更紧的输出尾界产生价值。但不能让两臂使用同一空间后仍宣称A17拥有不同物理精度，也不能故意禁止普通ROM复用Schur更新。此处是标准投影代数的排除性结论，不宣称新原创定理。

## 多源和vector块是实质约束

包内人为共线的两源例有非零有理joint chart，重放通过；独立一般两源控制则有3个复receiver-dark方向，而加入feedback行相容后joint kernel在1e−7/1e−10/1e−12三个阈值均为0。两者都是合成矩阵，不是真实天线场景，更不证明真实P1是否有解。

独立3×2块反例允许每个极化分量各自精确拟合，但共同各向同性标量拟合仍有相对残差.38125。它验证了三分量及所有source必须共享同一个材料，不能逐分量读出不同介电常数。零场chart明确拒绝，没有分母epsilon或隐藏clip。

## 直接用F3已付费状态检查直线current路径

设k=j1−j0、T=L0 k、B=E0、Z=Gk。这些量都可从缓存相乘取得，无需R0。若选择k(t)=t k，某单元存在共同标量材料增量的充要行条件为

\[
tT_{i,:}=d_i(t)(B_{i,:}+tZ_{i,:}).
\]

对全部17个缓存端点，原完整endpoint关系最大相对误差4.82e−13，接收差关系最大7.83e−13；说明本次读到的数组遵循同一物理合同。

对三个合格weak端点，t=.25/.5/.75的最优共同标量拟合仍有相对行残差约5.77e−4–2.06e−3，显著高于原端点数值误差。其读出的材料增量还产生约2.91e−4–1.34e−3的固定相位相对缺陷。**这些具体的直线current插值不对应精确共享标量、固定相位P0材料。**

这里没有调用Maxwell验证读出材料，没有产生新的物理候选；只是一项明确的可实现性反证。它不排除弯曲current路径、其他有理方向、finite-noise近似或新的P1，也不能和旧F3的q直线中点full-forward结果混为同一实验。

## 本轮不扩大的主张

原A17冻结GN定理未修改；F3 D不是A17；正确source体积公式不等于已运行正确posterior。DCTMC/CSI/ROM等相关工作的判断沿用包内来源与访问限制，本轮未进行新原创性排查，不宣称理论首创。无需以论文叙事取代P1机制实验。
