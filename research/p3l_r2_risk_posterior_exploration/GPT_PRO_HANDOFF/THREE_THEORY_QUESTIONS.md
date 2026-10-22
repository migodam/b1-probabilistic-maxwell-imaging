# 三个具体的新理论问题

## 1. 材料风险方向与Maxwell信息方向如何统一，而不是事后冠名

已知真实risk改善、等容量/空间置乱控制，以及第二批的Maxwell训练谱与随机decoder-envelope控制。请推导固定观测函数下，何种可计算的局部风险算子能够优于仅按J最小奇异值分割；明确哪些信息只是decoder导数/位置编码，哪些需要真实Maxwell谱。应给出可证伪的必要对照，而不是声称LIS首次用于逆散射。还需解释为何局部材料输出优于受限latent MLP、而大部分真实风险收益又位于prior线性切空间。

## 2. 饱和decoder与强弱耦合能否形成可归一化的非Gaussian conditional transport

10维proper NLL改善不能补齐余下56维；多中心/full-support/高weightESS都不足以保证探索。请利用实际被拒提案的decoder secant/Maxwell余项、sourceSMC祖先丢失、field2连续CDF与材料矩，构造有明确密度/Jacobian及合法Full校正的提案。需要说明它解决什么条件耦合，如何避免强方向Gaussian体积失配，并预言一个比已有控制更高的material ESS/RHS或更低的MCSE，不能只加flow层数。

## 3. 非对合弱分布下，什么探索或概率积分值得支付物理费用

真实P1的rank4弱场bracket有显著normal leakage；有限噪声路径仍能走，但弧长平均无风险收益。请研究以风险为目标的横向补偿、tube/多chart测度或随机核：什么时候必须加入弱分布以外的方向，什么时候局部采样的统计收益小于一次监督风险动作？提出含decoder余项、Maxwell余项、谱gap/rotation及统计误差的可检验成本判据。任何不变分布或probability claim均须给出参考测度、反向概率/重叠权重；零宽trajectory occupation不是后验。

以上都需要新的实质数学或明确反例，不要求维护A17、Gaussian、单流形或某个NN架构。未来实验与训练必须有新阶段授权，不能自动复用已封账余额。
