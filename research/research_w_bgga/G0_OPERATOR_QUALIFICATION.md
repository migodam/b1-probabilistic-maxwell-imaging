# G0 — 公共空间与精确吸收资格

状态：**QUALIFIED_DISCRETE_FULLWAVE_FOR_TESTED_CONTRACT**；机制stress的部分局部rank仍UNRESOLVED；连续与统一prior界未建立。

以固定P3空间、公共完整target ROI，取代依目标位置/半径重网格的旧接口。目标圆盘采用384边形cut-cell积分，解析incident进入材料右端。散射场的增量右端为b_X+V_Xu_W，不能机械套用固定f的total-field公式。完整target feedback由(I−C V)求解保留，full/woodbury/Schur的场、响应与方程残差同时审核。

12配置fixed-v-fitted白化差.006925–.020328；28精确对照场相对误差最大2.04e−13、响应5.17e−14、backward5.51e−17。实际JVP/实伴随约4e−15，shape有限差分单列。无jitter、伪逆、按结果增秩或full fallback。低损耗高对比inverse-norm估计不是严格稳定下界或共振证明。

公共2,257端口需要7.64 CPU秒/2257RHS setup；完整Green images1.026GB。7案例平均query CPU：full .33988、woodbury .57170、Schur .87109、背景预条件Krylov .47185。ordinary goal-ROM .01114，但held-out误差.237–6.517不合格。没有独立完整冷部署速度PASS。

原始：[资格结果](results/g0_qualification/RESULTS.json)、[算子结果](results/g0_operators/RESULTS.json)、[setup](results/g0_operators/SETUP.json)、[导数](results/g0_derivatives/RESULTS.json)、[stress表](results/G0_STRESS_TABLE.csv)。NPZ在本轮Release，按阶段目录恢复。

四类各3fixture，共9独立物理parent，known-support条件研究。C2/C3是匹配观测合同，不能当独立重复。finite15先验与tangent Gaussian信息不混用；C4_1只建立单端口当前观测与future25不同的经验数值见证，主观测能够区分，目标位置没改变。完整解释见中文报告。
