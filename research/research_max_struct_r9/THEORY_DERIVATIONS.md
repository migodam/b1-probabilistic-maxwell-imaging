# R9：非凸结构、A/A+与可执行修正

## 1. A+解决了哪一层问题

本轮沿用二维scalar-TM有限VIE，时间因子exp(−iωt)，实材料contrast m=εr−1，固定μ；三维与连续PDE不由本轮自动覆盖。令

\[
L(m)=I-\operatorname{diag}(m)G,\quad
J_m=L(m)^{-1}\operatorname{diag}(m)E,\quad F(m)=SJ_m.
\]

所有最终候选以原始完整forward及同一先验评价
\(\Phi(m)=\|F(m)-Y\|_F^2/(2\|Y\|_F^2)+R(m)\)。这里沿用历史归一化目标，不将它重新称作未经定义的白化似然。材料576维，box与先验不变。

A在固定current时改变本构残差的空间metric。用既有记号，d_i=β−Gii，z_i=m_i/(1+d_i m_i)，H_i=Q_i+d_iJ_i，
\[
R_A=\beta(I+D_zG_{\rm off})(J-D_zH).
\]
固定其余坐标后，对单个z_i的条件结构可被解析利用；但多坐标同时变化仍有z_i z_j耦合。A+的Möbius候选与MM接受减少这个材料读出的内优化成本。它保留既有current更新，既不增加测量信息，也不消除完整forward的resolvent。

**DERIVED：**A与A+应先在相同global current更新下比较，材料读出的增量应对照global-dual；sourcewise A+应对照sourcewise Dual-CIE。仅比较global A+与sourcewise Dual-CIE不能将差别归因于A材料metric。Möbius内更新做得精确，不意味着整个联合问题或最终Φ已变凸。

## 2. 单调变换与可逆坐标为什么不能消灭局部极小

**PROVED_UNDER_ASSUMPTIONS：**若h严格递增，则m*是Φ的局部极小当且仅当它是h∘Φ的局部极小。证明直接由局部函数值序关系得到。若h可微且h′>0，在驻点处
\(\nabla^2(h\circ\Phi)=h′\nabla^2\Phi\)，因此惯性不变。

若m=ψ(z)为局部微分同胚，在驻点处
\(\nabla_z^2(\Phi\circ\psi)=D\psi^T\nabla_m^2\Phi D\psi\)。Sylvester惯性定律保持正、负、零方向计数。边界处使用相对可行集的局部极小定义，可逆连续映射仍保留局部极小。坐标可以改善条件数和有限步轨迹，却不能据此宣称去掉原目标的错误极小。

这不是说A的辅助metric不能改变轨迹：它确实改变了辅助目标。这里排除的是“单靠可逆Möbius改写已经凸化原Φ”的说法。 noisy/prior问题中GT也未必是驻点；历史若干成功例不足以证明“GT吸引域半径扩大”。

## 3. 噪声进入DPCS的具体通道

S=UΣV*，冻结强奇异子空间的读出为J_D=V_rΣ_r⁻¹U_r*Y。若观测噪声独立圆复高斯，E|N_rs|²=σ²，则
\[
\mathbb E\|\Delta J_D\|_F^2=s\sigma^2\sum_{i\le r}\sigma_i(S)^{-2}.
\]
这是预期噪声放大，不是对实际GT电流误差的观测。本轮只复用此结构诊断，不重复其机器精度验证。固定相对奇异阈值不随测量噪声变化，可能硬锁弱方向中的噪声；R3已有clean/release失败表明这不是所有错误材料终点的统一解释。

测试的简单修正：以公开σ定义一个完整数据99%噪声球
\(\epsilon^2=\sigma^2\chi^2_{2ns,0.99}/2\)，n为接收数、s为源数。取最小连续奇异前缀r，使\(\|(I-U_rU_r^*)Y\|_F\le\epsilon\)。因所有正交投影收缩，在同一事件\(\|N\|_F\le\epsilon\)上均有相同上界，选择r不需要为每个r另用一个点态概率界。

**边界：**这只给出数据近似的公开discrepancy规则，不证明r是最佳可恢复秩，也不证明被保留模式没有噪声，更不证明CIE的Fourier可行集合包含真解。与原rTol=.001算法分开列账。额外receiver SVD收费。非IID或缺失σ的任务不套用该规则。

## 4. 真正改变中间目标：完整current松弛

这是已有WRI/CSI思想的控制，不是新发明。固定λ>0，令
\[
\Psi_\lambda(m,J)=\frac{\|SJ-Y\|_F^2}{2y_d}
+\frac{\lambda\|L(m)J-D_mE\|_F^2}{2e_d}+R(m),
\quad y_d=\|Y\|_F^2,\ e_d=\|E\|_F^2.
\]
完整current内解满足
\[
KJ=S^*Y/y_d+\lambda L^*D_mE/e_d,\quad
K=S^*S/y_d+\lambda L^*L/e_d.
\]
在L可逆、λ>0时K正定。设e=J−D_m(E+GJ)，精确内解的包络梯度为
\[
\partial_{m_i}\min_J\Psi_\lambda
=-\lambda\,\operatorname{Re}\sum_s\overline{e_{is}}(E+GJ)_{is}/e_d+\partial_iR.
\]
这里没有R2固定CIE profile中的材料相关p(m)权重；不能漏掉p导数后冒称同一算法。实现复用既有G/E/S与物理corrector。

定义A(m)=S L(m)⁻¹、r=F(m)−Y，消去current的等价表达为
\[
\min_J\Psi_\lambda
=\tfrac12 r^*\big(y_d I+(e_d/\lambda)AA^*\big)^{-1}r+R(m),
\]
对多源逐列求和。它改变了材料相关的残差metric，弱罚权可压低某些端口方向；因此低松弛目标不等于低Φ，且没有全局凸性保证。λ→∞在可逆有界区域恢复原Φ，但λ→0一般会丧失材料信息，而非自动得到有用凸逆问题。

固定初值的比例尺度为\(\lambda_0=\|A(m_0)\|_2^2e_d/y_d\)。本轮开发协议仅测试λ/λ0=.01,1,100各100次，与固定λ0共300次对照。最后统一原物理corrector3000次。中间延拓传递最后iterate，不以原Φ逐阶段否决；同时保留所有已付费原Φ候选，最终同规则选最小Φ。这一差别被明确记账，不声称与历史逐阶段选择完全同轨迹。

## 5. 同一数据的目标延拓

定义\(Y_t=(1-t)F(m_0)+tY\)，固定原y_d和R，仅改变辅助target；t=.1,.3,.6各250次，最后t=1原corrector3000次。t=0时m0的数据残差为零，但先验可使其非驻点；连续路径可分叉、撞边界或丢失稳定性，不能预言一定连到正确解。中间Yt是算法构造，不是新观测，更不是合成低频或加损耗的真实观测。

这个控制能区分“原目标上的局部步太短”和“另一目标路径能否到达不同合法终点”。等迭代不能替代等CPU或等物理LU/RHS，须读完整成本曲线。

## 6. 交叉熵与有条件的凸化

交叉熵随机搜索改变候选分布；KL/CE数据loss改变误差几何；熵正则运输凸的是运输计划；相对熵prior引入先验。这四者都不自动凸化m↦F(m)。归一化复场幅度来套概率loss可能丢掉总幅度和相位，必须另论物理合同。此轮不无依据训练网络或大规模CE采样。

近先例van Leeuwen–Yang的特殊二次材料目标要求整个状态空间就是配对source/receiver Riesz基的有限span、固定内积与仿射PDE系数；方形可逆数据还可给直接线性解。当前48×6/12外部矩形端口与576个内部材料自由度不满足此资格。把内部状态强行截到端口span可以做出凸小模型，但不能称其为保真完整Maxwell材料逆问题。

近先例与原文限制详见delegated/nonconvex_prior_art/PRIOR_ART_ZH.md。最终成像判断等待实际封存结果；上述推导不预先担保收益。

## 7. 错误端点的实际曲率，而非按停止标签猜测

为了区分尚有负曲率的终点和受box约束的稳定错误状态，本轮对两个已冻结公开端点计算原Φ的有限模型Hessian。令r=SJ−Y、a=L^{-*}S*r，则对实方向v：
\[
\dot J=L^{-1}D_vQ,\quad \dot Q=G\dot J,\qquad
\dot a=L^{-*}(S^*S\dot J+G^*D_va),
\]
\[
H_\Phi v=\operatorname{Re}\sum_s(\overline{\dot a_s}\odot Q_s+
\overline{a_s}\odot\dot Q_s)/y_d+H_Rv.
\]
smooth-TV的**精确**二阶差分权重是
\(w=(2\,smooth+tv\,\epsilon^2/(d^2+\epsilon^2)^{3/2})/N_{edges}\)。原Regularizer返回的\((2smooth+tv/\sqrt{d^2+\epsilon^2})/N_{edges}\)是IRLS/MM权重，不能冒称此处的精确Hessian。没有因此指控原求解器错误；用途不同。

只对距box两端超过1e−5的自由坐标构造矩阵；保留全部特征值、对称缺陷、特征对残差和一次必要的新HVP接线检查。它不是活跃约束打开后的全临界锥证明，也不是连续Maxwell或全局最优证书。

**EMPIRICAL，公开开发：**40ac的A+错误终点有55自由/521活动坐标，原投影梯度约9.43e−9，自由Hessian最小特征值+4.27e−9，没有发现该子空间的负曲率。b20d的普通Global终点有49自由/527活动坐标，投影梯度约9.67e−7，自由Hessian有一个−1.62e−8特征值；不能把它写成已经认证的局部极小。固定8个合法有限步中，最佳原Φ下降约3.17e−8，随后另作相同3000步物理续跑与不扰动续跑的配对。存在一个微弱负方向不等于它能恢复材料。

## 8. 新确认与二分支控制

在新六材料生成/评分前冻结噪声规则、A+版本、松弛延拓与预算。两项附加组合仅为已知方法的付费对照：noise-aware A+＋profile，及noise-aware Global-Dual＋profile；均按最终原Φ选一项，始终支付两次完整运行。确认封存时先保存选择，再读取GT。这样可检验互补是否在新材料上成立、A材料metric是否仍提供普通Global不能解释的增量。

这是改变搜索路径/保留候选的机制研究。两分支本身不是新数学，也不保证材料误差单调变好；小Φ可能选错材料。不能用事后GT oracle替代真实选择，也不能把两次求解的费用漏掉。
