# 后验目标与概率边界

固定共同target：原A30 conditional smooth P0、原复Gaussian观测、8sources/32receivers、固定gain1、原背景和mesh。q=(六个单位uniform shape,t)，rho=(2t−1)L，L=log(9/4)。

\[
\pi(q\mid y)\propto \exp[-\Phi(g(q))]\,2L\,p_\rho((2t-1)L),\qquad \Phi=\tfrac12\|W(F(g(q))-y)\|^2.
\]

前六坐标密度1；2L只乘一次。共同weight尺度已正确边缘化，不把7维设均匀，也不额外乘材料Gram determinant。若改报告材料Hausdorff体积密度，应另处理Jacobian及所有前像；本轮没有这样改目标。

rho的端点属于合法材料支持闭包，但密度为0，无边界概率原子。报告“positive prior density”是点密度的有限性标记，不是把零密度边界删出支持。任意边界邻域仍有非零概率，尾质量按.72*delta²衰减。不能只因点密度0宣布物理不合法；反过来，点密度正也不意味着该材料区域有显著mass。

Gaussian标签交换仍是离散同材料重复；碰合Gaussian还会额外退化。共同尺度quotient不等于全局唯一atlas。保持原有带标签推前测度，没有擅自排序或乘除2。

固定对角反射MH满足对称q提案，所以接受比为full posterior密度比。初始边界测度零点的扩展规则明确；所有正常连续内部点使用标准规则。拒绝重复没有丢弃，也没有把proposal行伪装成posterior样本。独立重放了16次proposal、接受随机数和保留链，全部一致。

Data-feasible和ΔΦ≤2只用于候选评价，不参与MH目标截断；尤其不能用noise-ball通过率给候选“posterior probability”。当前没有独立7Dposterior reference、evidence、event partition mass、region probability、credible interval或模式发现认证，全部NOT_ESTIMATED/null。

每个连续图像单点质量仍为0；图库仅展示candidate-only端点。未展示质量、未探索质量及模型外质量均无法量化。原P0条件支持不是一般二维未知结构或开放三维支持。
