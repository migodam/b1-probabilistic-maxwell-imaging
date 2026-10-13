# 材料事件概率为何必须带条件分布一起交付

这里给出精确Bayesian恒等式与有限近似的误差界，说明实验中哪些表现来自信息保留。它们是条件概率与KL分解的直接结果，不作为新原创定理申报；实验贡献必须另外由Maxwell机制、概率参考和计算成本支持。

设同一声明后验为 \(\pi(d\chi\mid y)\)，互不重叠材料事件 \(A_j\) 覆盖全部已建模支持。令 \(p_j=\pi(A_j\mid y)\)，\(\pi_j\) 为事件内条件分布。一个给定材料后与散射观测独立的传感器给出似然 \(l_c(\chi)=p(c\mid\chi)\)。则

\[
m_j(c)=E_{\pi_j}[l_c(\chi)],\qquad
p_j(c)=\frac{p_jm_j(c)}{\sum_kp_km_k(c)}.
\]

证明是把 \(\int_{A_j}l_c\,d\pi\) 写为 \(p_j\int l_c\,d\pi_j\)，再对整个支持归一化。分母包括未展示的已建模事件。模型外没有声明的prior时，不能为它凭空赋予数字质量。

只有事件质量和代表图一般无法求出 \(m_j\)。用代表点 \(\widehat\chi_j\) 的 \(l_c(\widehat\chi_j)\) 代替该期望，隐含了事件内部全部材料具有同一context响应的假设。一个事件可以同时包含两个不同空腔体积；窄的独立几何传感器在它们处的似然差异很大，即使代表图和原始事件质量完全准确，该替代也会失败。对于任意事件内条件分布，这种替代恒成立的条件是context似然在事件内几乎处处为常量；否则改变内部质量位置就会改变更新。

因此图中实心与大中心空腔的事件概率是有用的第一层输出，但完整后续推断还需要事件内加权材料样本、或可积分的条件密度。`TopK_conditional_samples_plus_remainder`与whole-particle更新相同，是保留相同支持/权重后的表示恒等关系，不是另一个算法战胜了采样器。

若各事件的条件似然期望近似满足 \(|\widehat m_j-m_j|\le\delta_j\)，记 \(D=\sum_jp_j\delta_j\)、\(Z_c=\sum_jp_jm_j>0\)，且近似总证据也为正，则

\[
\|p(c)-\widehat p(c)\|_1\le2D/Z_c,
\qquad \operatorname{TV}\le D/Z_c.
\]

证明：令 \(a_j=p_jm_j\)、\(\widehat a_j=p_j\widehat m_j\)。对 \(a/Z_c-\widehat a/\widehat Z_c\) 加减 \(\widehat a/Z_c\)，三角不等式给 \(\ell_1\le(\|a-\widehat a\|_1+|Z_c-\widehat Z_c|)/Z_c\le2D/Z_c\)。对于取值在 \([0,L_{max}]\) 的事件决策损失，由近似后验选出的Bayes动作在真实后验下的regret不超过 \(2L_{max}\operatorname{TV}\)。这需要实际的误差上界；本实验的有限积分细化/样本bootstrap不提供对所有context的uniform界，不能将该公式当成本阶段的已认证保证。

还有一个重要的坐标问题。decoder在inactive coefficients上不是一一映射。对可测decoder \(D\) 及适当的条件分解，KL链式法则给

\[
\mathrm{KL}(P_x\Vert Q_x)
=\mathrm{KL}(D_\#P_x\Vert D_\#Q_x)
+E_{D_\#P_x}\mathrm{KL}(P_x(\cdot\mid\chi)\Vert Q_x(\cdot\mid\chi)).
\]

第二项衡量同一物理材料的latent冗余坐标。原prior中inactive coefficient条件分布为uniform，Maxwell likelihood对它不变。因此C1改善latent KL可能部分来自这些不可观测的冗余坐标；不能直接宣称material-field KL改善同样数值。本阶段只测量latent KL与明确的材料事件/CDF，没有估计完整pushforward material KL。

补偿projection也有类似边界：无约束白化Jacobian投影是局部least-squares诊断。proper nuisance prior有有限尺度且gain有边界，Bayesian marginal还包含logdet、prior体积和截断常数。非线性有限secant给出一个可行补偿见证，最小残差仍不能替代上面的 \(p_j\) 或evidence。实际两频小孔记录的投影后斜率约2.975，未显示普遍的leading cubic cancellation；不能从局部投影形式推断所有结构都变成五次量级。
