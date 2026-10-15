# 理论修订：材料注入、局部评分与有限成核

原离散Maxwell方程为 Lχ c = diag(χ)e_inc，Lχ=I−diag(χ)G；e=e_inc+Gc，y=Sc。c/e均沿用原质量归一化，χ是raw susceptibility。所有照明同一χ。

## 精确有限注入不是任意dark current

对于两个已求得的合法材料χ0、χ1，令Δχ=χ1−χ0，则

`L1(c1−c0)=diag(Δχ)e0`，故 `F1−F0=S L1^{-1} diag(Δχ)e0`。

这来自原方程相减，不是原创恒等式。相较Jacobian `J0d=S L0^{-1}diag(d)e0`，有限编辑使用新材料resolvent L1。40个实际端点的恒等式缺陷最大4.10e-15；这验证本轮injection与场/材料合同，不是continuum认证。没有将任意receiver-dark current当成材料候选。

## 有界线性补偿的正确范围

对当前状态及单个trust box，`min_{l<=d<=h} ||r+Ad||`。任意||v||≤1给合法下界

`L(v)=max(0, v^T r + Σ_j min(l_j (A^T v)_j, h_j (A^T v)_j))`。

任一合法d给上界U。上下界对这个凸线性问题成立。只有额外取得对同一整个box的uniform nonlinear remainder ε，才可在该box用L−ε>tau排除或用某个实际试点加误差界支持接受。多步更新改变box/resolvent，当前L不能排除后续路径；U>tau只是该attained点失败，甚至不是线性最优下界。

实际反例：v1有13个U>tau但三步非线性通过。不得把单步upper、未收敛primal残差、局部dual或profile minimum升级为全材料不可能性或posterior odds。

## 成核结构必须是约束，而非展示掩膜

core固定χ=0但不保护周围材料，不能保证封闭孔；guard保持Reχ≥.30且闭合，使所有core在8连通空气定义下被包围。最低6cell及2cell半径排除单像素阈值收益。仍有额外非core相变化，v2保护使其减少但损失成功率，说明补偿自由度有实质作用。

## 数据相容和成像质量不同

减小||Wr||未必减小材料误差，尤其1152实材料自由度对应512观测并有弱方向。观测数本身不是有效秩定理，但本轮full vs Haar的实证明确暴露过拟合。minimum-change规则只在已付费合格端点集合中选最靠近N3者；没有声称求得全局最小改动解。

## A17的作用边界

A17 paired law是在固定材料注入、残差、先验和指定Galerkin差异下刻画局部二次步差。全材料exact adjoint的S Lχ^{-1}本身包含完整feedback。增设空气core/guard会改变活动约束；只有固定free face后重新限制H/J才有相应局部线性结论，不能沿用无约束公式认证有限拓扑变化。当前没有相同成本paired/Krylov付费比较，独立增量NOT_RUN。未改写A17冻结稿。

## 能声明与不能声明

可以声明：特定exposed离散TM目标下，实际构造了多个材料上显著不同、满足同一数据球和被动盒的端点。不能声明：统计不可区分、后验显著、完整覆盖、连续可行路径、全局不唯一性解决、网格收敛、一般三维效果或算法原创定理。标准GN、有限resolvent恒等式与约束LS本身不构成新颖性。
