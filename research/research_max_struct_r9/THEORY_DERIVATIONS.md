# R9：非凸结构、A/A+与可执行修正

## 1. A+解决了哪一层问题

本轮沿用二维scalar-TM有限VIE，时间因子exp(−iωt)，实材料contrast m=εr−1，固定μ；三维与连续PDE不由本轮自动覆盖。令

\[
L(m)=I-\operatorname{diag}(m)G,\quad
J_m=L(m)^{-1}\operatorname{diag}(m)E,\quad F(m)=SJ_m.
\]

所有最终候选以原始完整forward及同一先验评价
\(\Phi(m)=\|F(m)-Y\|_F^2/(2\|Y\|_F^2)+R(m)\)。这里沿用历史归一化目标，不将它重新称作未经定义的白化似然。材料576维，box与先验不变。

A在固定current时改变本构残差的空间metric。完整总电场为Q=E+GJ。用既有记号，d_i=β−Gii，z_i=m_i/(1+d_i m_i)，H_i=Q_i+d_iJ_i。这里的G_off仅保留不相交2×2空间块内的非对角耦合，即G_off=⊕_B(G_BB−diag G_BB)，不是完整G−diag G；完整物理forward中的G仍保留所有单元间的作用。
\[
R_A=\beta(I+D_zG_{\rm off})(J-D_zH).
\]
固定其余坐标后，对单个z_i的条件结构可被解析利用；但多坐标同时变化仍有z_i z_j耦合。A+的Möbius候选与MM接受减少这个材料读出的内优化成本。它保留既有current更新，既不增加测量信息，也不消除完整forward的resolvent。

**DERIVED：**A与A+应先在相同global current更新下比较，材料读出的增量应对照global-dual；sourcewise A+应对照sourcewise Dual-CIE。仅比较global A+与sourcewise Dual-CIE不能将差别归因于A材料metric。Möbius内更新做得精确，不意味着整个联合问题或最终Φ已变凸。

**固定current时的梯度解释（DERIVED，非新定理）：**令e=J−D_mQ、a_i=1+d_i m_i、T(m)=β(I+D_zG_off)D_(1/a)，因此R_A=T(m)e。在固定J、Q、β及Hermitian源metric Ω时，材料辅助目标为f_A=ρ‖T(m)e‖²_Ω/2+R_prior。令C_i为只有第i行等于Q_i的矩阵，⟨X,Z⟩_F=tr(X*Z)，则
\[
\partial_i f_A=-\rho\operatorname{Re}\langle TC_i,R_A\Omega\rangle_F
+\rho\operatorname{Re}\langle(\partial_iT)e,R_A\Omega\rangle_F
+\partial_iR_{\rm prior}.
\]
第一项是传播后本构切向的加权读出，第二项来自metric本身随材料变化。故A不是一个在整个材料内循环保持固定的预条件矩阵；离开e=0时不能漏掉第二项而仍声称优化同一目标。R4的onehop_material_value_gradient与R7/R8的A+候选接受都保留该目标。在e=0处，固定current辅助目标的二阶变分是ρ‖T D_(δm)Q‖²_Ω加先验二阶变分；这不删除物理forward的resolvent，也不赋予全局材料凸性。A+缩短材料内优化与A改变辅助路径是不同层次的作用。

## 2. 单调变换与可逆坐标为什么不能消灭局部极小

**PROVED_UNDER_ASSUMPTIONS：**若h严格递增，则m*是Φ的局部极小当且仅当它是h∘Φ的局部极小。证明直接由局部函数值序关系得到。若h与Φ均为C²且h′>0，在驻点处
\(\nabla^2(h\circ\Phi)=h′\nabla^2\Phi\)，因此惯性不变。

若Φ为C²且m=ψ(z)为C²局部微分同胚，在驻点处
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

近先例与原文限制详见delegated/nonconvex_prior_art/PRIOR_ART_ZH.md。第一阶段的实际封存结果已经完成，见RESEARCH_REPORT_ZH.md；上述推导不预先担保成像收益。

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

时间边界：上句的生成前冻结仅指主方法；二分支组合、single-noise补充与weak-handoff补充是在数据生成后、首次确认GT评分前加入，不能全部称作预生成注册。实际时间与对应配置保留在audit/FINAL_SCIENTIFIC_REVIEW.md。

这是改变搜索路径/保留候选的机制研究。两分支本身不是新数学，也不保证材料误差单调变好；小Φ可能选错材料。不能用事后GT oracle替代真实选择，也不能把两次求解的费用漏掉。

## 9. 从已付范数识别未解决的端口误差传播

原尺寸三维的matched-grid控制移除了粗细网格差异，CIE系仍出现小辅助损失、大完整物理残差。这里复用已有Maxwell defect-transport恒等式，不将它记作新理论。固定m，设辅助current为J，e=L(m)J−D_mE，d=SJ−Y，物理current为J_m，则
\[
J_m=J-L(m)^{-1}e,\qquad F(m)-Y=d-SL(m)^{-1}e.
\]
这说明辅助current对数据拟合得很好，仍可能由未消除的本构缺陷经完整resolvent传播而形成较大的真实端口误差。局部材料metric和系数变化停止量都不自动控制这个传播项。

用实际实现的源归一化范数
\(\|Z\|_d^2=s^{-1}\sum_j\|Z_j\|^2/\|Y_j\|^2\)、
\(\|Z\|_e^2=s^{-1}\sum_j\|Z_j\|^2/\|E_j\|^2\)，
以及当前CIE残差缩放p(m)，令a=‖d‖_d、b=‖F(m)−Y‖_d、c=‖p(m)e‖_e。由正反三角不等式，c>0时，这个**实际残差方向**的增益满足
\[
\frac{\max(0,b-a)}{c}\le
\frac{\|SL(m)^{-1}e\|_d}{\|p(m)e\|_e}
\le\frac{b+a}{c}.
\]
右端不是该算子的operator-norm上界，不能用来给所有未知材料的停止证书。以下只使用已保存的有限精度范数，没有严格线性解/离散误差包络，因此数值区间标记为EMPIRICAL_DIAGNOSTIC，而非认证界。

|matched-grid方法|辅助data RMS a|完整物理RMS b|普通current状态RMS c|实际方向增益区间估计|
|---|---:|---:|---:|---:|
|source CIE|.01676|.36712|.02892|[12.11,13.27]|
|Global-Dual|.02743|.37576|.05000|[6.97,8.06]|
|A+ global|.02742|.37734|.05000|[7.00,8.09]|
|LSIE|.01807|.70295|.06129|[11.17,11.76]|

不同β的c定义不同，不能单看增益大小排序算法；Global-Dual/A+这一配对使用相同缩放，可直接观察到A+并未明显降低该端口缺陷。A+行的c是实际日志中普通current残差，不能误称材料one-hop损失。无需再运行正演或读取GT就能得到这些数，代码与原始输入在`experiments/analyze_original_residual_transport.py`和`RESULTS/original_residual_transport/`。

上述范数不能回答梯度/材料恢复问题，因此另冻结matched-grid空气、Global-Dual和A+三个初值的共同八步原物理校正，保留原九臂不变。三臂已完成：Global/A+的χ相对L2分别由.42275/.42347降至.24899/.24950，空气仍.97501。这是已公开开发对象上的机制控制，不是盲测，也不是全局收敛证明；原辅助端点仍有有效物理方向，而A+没有独有收益。全部新增费用与同数据材料评分见短校正报告。

这次使用的精确离散梯度为：Q=E+GJ，w_s=1/||Y_s||²，P=L^{-*}S*((SJ−Y)w)，则对原尺寸目标Φ=mean_s||SJ_s−Y_s||²/||Y_s||²，g_i=(2/s)Re Σ_{分量,源}conj(P_i)Q_i。实各向同性χ在三个分量共享，不把三分量独立优化成新材料类别。公开空气点仅做一次必要方向接线，解析方向导数−.61365626与h=1e−5差分−.61365319相符；数值接线不作为新恒等式贡献。
