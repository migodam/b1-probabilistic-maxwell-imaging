# 受控墙–目标家族的真实全波概率参考

本阶段取得了可靠到预注册检查精度的低维数值posterior，而非原高维任意未知墙体的完整posterior。两个公开墙场在旧物理存在性搜索完成后，先冻结成一个新的可归一化材料家族，再产生16个独立参数parent、各2份独立噪声观测。墙体分支与目标位置区间在这个新prior中相关；不能将它解释成独立墙prior乘任意目标位置prior。

新观测生成前的prior与协议见 `configs/W_CONTROLLED_PRIOR.json`、`P1_PROBABILITY_PROTOCOL.json`。原宽材料类未定义的最小特征合同不被追认PASS。这里的Research W P1，也不是此前64维Gaussian材料先验的P1。

## 测度与实际积分

分支b具有质量0.5/0.5；分支内u、目标介电对比c、共享gain幅度a和相位φ均具有明示的有界uniform prior。对所有五频、四源、十二收的240个复数据，使用480维实Gaussian likelihood，Phi=一半完整白化残差平方，保留其归一化常数。

五频/四源共享同一gamma=a exp(iφ)。先形成完整likelihood，再积gamma一次，不能分别给各source或frequency配置自由gain，也不能用profile最小值冒充证据。

材料积分采用Full3、Full5、Full9规则。CDF在每个声明截点重新积分原密度的截断区间，分母为同阶完整证据；不把节点铺成bin，不将截断prior重新归一化。所有规则由实际高阶FEM响应组成，未用response surrogate构造参考。去重后共678个全波材料向量，13560 complex RHS、3390 LU；16个新parent均值生成另付320 RHS、80 LU。总物理成本13880 RHS、3470 LU。

## 已执行的数值资格

32/32观测通过Full5/9及直接截断CDF的冻结门槛：

| 核验量 | 32观测最大差 |
|---|---:|
| 分支事件质量：Full5 vs Full9 | 0.000321111 |
| 声明CDF检查点：截断5 vs 截断9 | 0.000775065 |
| 预注册容差 | 0.02 |

共享gain的11/21阶检查已经满足2×10⁻⁴的log-evidence差门槛；各CDF均经过范围与单调性检查。每个观测的完整阶次记录、权重、条件证据、均值和SD均保留在 `results/P1A/reference/`，所有Full响应规则在 `results/P1A/rules/`。

这是一组独立Full节点集上的经验求积收敛证据，不是连续PDE严格认证、全prior误差上界或模型外质量上界。原未知墙体类仍可能有未建模分支；其质量为null/unquantified。没有用一条短MCMC或ESS充当概率真值。

新P1-B观测使用另一个种子、另16个parent。其连续Full5/9参考与实际650材料节点有限规则分别核对，不能将P1-A的资格自动赋给新观测。
