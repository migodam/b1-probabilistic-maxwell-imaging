# P0-R 可表达性：回放、投影与非线性拟合分开读

状态：**WORK_IN_PROGRESS / 保存结果说明**。本文件只读现有结果、审计与源码文字；新增 Maxwell、LU、RHS 均为0。科学判决与后续实验由主线程维护，本页不修改门或账本。

首先问“同一个连续材料能否解释原观测”，再问“某个逆表示能否表达相容材料”。直接把真值栅格化后的失配，不是该逆表示全部可行参数的最小失配。

## 1. 原连续真值的同合同回放

[results/R0/RESULTS.json](results/R0/RESULTS.json)保存44行：36个材料–频率状态和8个完整数据向量；两seed、两个band、cal与固定probe8的8个`d_mean`均为0。其记录的独立噪声长度20.68480395–22.77434031，属于`d_observed_truth`，不能与均值距离混为一谈。原有跨网格资格被复用而未整轮重跑，文件对此明示。

这个事实支持同一已保存前向合同可以复现原均值，不等于该材料已经满足全部几何先验，也不是所有新候选的数值资格。

## 2. 两个wide对象的Oracle-near非线性拟合

|seed|band|d_mean|projected_gradient|iterations|success|stop|
|---|---|---|---|---|---|---|
|2026102201|wide|0.967572217|0.467013702|4|False|objective_floor_reached (not gradient convergence)|
|2026102202|wide|0.987976116|0.305986096|7|False|objective_floor_reached (not gradient convergence)|

原始路径：[seed2201](results/R3/oracle_2026102201_wide_TR/RESULTS.json)、[seed2202](results/R3/oracle_2026102202_wide_TR/RESULTS.json)。两者使用真值子像素投影派生支持及95%投影初值，身份都是`ORACLE_DIAGNOSTIC`。这与仅凭测量找到材料是不同任务。

两个保存端点的均值距离小于1，说明该**Oracle派生受限像素表示中**找到了数据拟合端点。它们都以objective floor结束，`success=false`，投影梯度约0.4670/0.3060；并未宣称KKT驻点、全局最优或表示误差下界。

## 3. 为什么完整几何合法性仍未证明

[几何审计](AUDIT/geometry/REPORT.md)与[COMPUTED.json](AUDIT/geometry/COMPUTED.json)指出：历史`minimum_feature=0.1`只有标量，没有宽度、滚动球、拓扑持续性、两相、尖角或误差规则。它的缺口同时影响连续真值和候选，不能只检查像素孤岛而默认真值已获认证。

审计解析包络支持真值的材料盒范围、被动性、跨频材料律、面积上界1.7628<3，以及距query的竖直间隔至少0.1996。seed2201的孔腔与上边界相交，形成开口凹槽和约35°尖楔；较强双侧滚动圆盘定义与容许角点的侵蚀拓扑持续性是不同的新合同，不能事后任选一个解释原先验。栅格诊断也不是连续几何证书。

因此，`d_mean<1`与“完整原几何类的合格材料”必须分开。当前证据不支持把Oracle端点升级为完整合法物理见证，也不支持把直接投影误差大解释为全像素类不可表达。

## 4. 资源边界

[R_SEALED.json](R_SEALED.json)记录R阶段已封账：12000 RHS、1195次因子分解、2251.515528 process CPU秒，额外授权0；disk与smooth_block controls未完成。本文件没有借用其他阶段预算补R实验。campaign仍在继续展示，未在此作最终判决。

直接投影明细见[P0_R_PROJECTION_ANALYSIS.md](P0_R_PROJECTION_ANALYSIS.md)，在线差异见[P0_R_ONLINE_VS_ORACLE.md](P0_R_ONLINE_VS_ORACLE.md)。


## R阶段封账时点说明

上述保存说明中的WORK_IN_PROGRESS与“后续未运行”描述R封账时点，不是整轮当前状态。R原始数字与阶段Gate保持原样；最终已完成的B2/C/受控P1-A/P1-B及其限制见[整轮裁决](SCIENTIFIC_DECISION.md)和[阶段汇总](GATE_SUMMARY.json)。本补记没有重开R预算或将Online失败改成PASS。
