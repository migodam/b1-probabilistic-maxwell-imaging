# 理论与数值核验：哪些成立，哪些尚未成立

## 响应坐标

固定正交坐标x=x_ref+Va+Wb，选择不依赖待求s的正交U(a)，s=U(a)^T(f(x)-d)/sigma。在声明域内存在一一可微逆且B=U^T Df W非奇异时，换元体积为sigma^r/|det B|。完整目标还包含原prior、未投影残差、分支和重叠权重。该公式是成熟换元工具，不作为新定理。

数学单元在条件体积、x²双分支与线性高条件数模型上通过。[原始检查](stage_a_response_coordinates/math/A1_MATH_RESULTS.json)对比了删除体积的错误控制。真实Maxwell图的有限求根成功不能代替全域单射/覆盖证明。离散求积点集对连续posterior的TV为1；本轮没有将两者直接比较成TV证书。

## 实际decoder的混合秩

原配置给出field=.65*x0、tau=.5。第二个二维子模型因此存在rank0/2/1三层；上下平台各有.22086048先验质量。原latent后验仍相对于二维Lebesgue测度定义，退化发生在材料/响应推送。分层积分须保留真实区域体积。[推导](audit/MIXED_RANK_DECODER_DERIVATION.md)

H(-1+e)=5e³/4-15e⁴/16+3e⁵/16，支持在正则物理材料作用条件下的sigma^(1/3)接缝尺度。真实Full-Maxwell有限路径斜率2.996758，latent seam J为0而材料作用非零。它只验证一个预定子模型方向，不是Full-66局部fiber存在性证明。[实际数组与指标](audit/seam_maxwell/RESULTS.json)

## 局部稳定与权重界

固定强方向上，利用有界Dg、D²g及Neumann逆扰动，评估了两个中心的充分半径。六个R=3测试均不满足噪声半径条件，界很保守。这不是不存在逆或后验不确定性的证据。`L_on_rho_ball`遗留机器字段在本实验只表示**固定W、固定弱坐标的一维线段**，不是全二维球。

范数计算、残差与曲率公式均为float64，未包含定向舍入证书；中心不一定满足选中响应残差为零。有限节点的权重振幅不称全域Delta，不据此宣称TV/KL上界实际成立。[独立审查](audit/A5_INDEPENDENT_REVIEW.json)

## Gaussian decoder矩

对联合Gaussian (S,U)，先条件化U|S=s，再用E Phi_N(U)|S=Phi_N(ell/sqrt(1+v_c))及相应二元CDF，积分delta(U)H(S)与平方。实际恢复的Gaussian checkpoint场—属性交叉块为零，才采用独立特例。任意相关/退化/尾部例另有独立求积分验证。

四通道是实/虚均值及各自方差，不把Hermitian covariance、pseudo-covariance混为一谈。原K32的ddof=0与population variance有限K差异保留。积分正确只描述近似q的矩，不描述未经校准的真实posterior。

## 任务低秩修正与功能蒸馏

对冻结基准m0和仅依赖y的特征psi，损失加权交叉矩C与输入二阶矩B给出已知RRR修正A_k=[CB^(-1/2)]_k B^(-1/2)。它筛选可预测残差，不能用GT残差PCA替代。实际嵌套OOF修正均无收益，有限类的负结果保留。

新的小型provider压缩的是原Gaussian的四通道功能，不是完整概率密度。其收益不依赖该Gaussian已校准，也不说明输出统计对应同一个Gaussian。若teacher head对统计误差为Lipschitz，可由输出差d和teacher风险R给出风险变化的保守上界2sqrt(R)d+d²；本轮未为深网给出紧全局常数，实际风险由验证样本测量。

## 原创性边界

IFT、Gaussian矩、RRR、RTO、Laplace、LIS和条件期望均有成熟先例。[已有工作与本轮主张](PRIOR_ART_AND_CLAIMS.md)。可能值得进一步研究的是实际decoder分层响应概率、有限预算下可验证的概率覆盖，以及独立确认后的材料统计压缩；目前尚无新Maxwell定理或强TAP性能主张。
