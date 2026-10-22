# 两批独立新对象与冻结合同

## 第一批：主要成像比较

原旧P3-L128个locked parent已暴露，仅作历史证据。本轮第一批R2_locked_00000–00255在第一份新观测产生前冻结，source revision `ccf3987fd3910c656fa515fc076aa636fb89c47e`，77检查点/主融合/轨迹/失败fallback/额外16-state普通优化/统计规则均写入METHODS_FROZEN_FOR_LOCKED.json。256个在线点、全部神经输出和轨迹闭合后，独立evaluator才打开GT。

两项主比较以parent配对bootstrap2000次、seed2026102292，给出97.5%双侧区间。主风险降低25.27%，区间[22.34%, 27.99%]；弧长动作没有改善。无parent删除，无point/context/model fallback。三训练seed不能被当三倍独立parent。

## 第二批：机制的独立确认

根据第一批无法拆分通道的明确缺口，以及P1灵敏度恒等式，另建两份protocol设计记录和MECHANISM_CONFIRMATION_FROZEN.json。21个checkpoint及source revision `5d5332a144b6665a9206d073f51ad4cdc572d852` 在第二批第一份观测前登记。新IDs00256–00511使用原master seed规则下新的独立parent/noise，不是旧数据的新命名。

只使用旧512train/128validation训练9个新控制，未使用第一批GT作新标签，也未按第二批选winner。三项主比较用固定seed2026102393、2000次parent bootstrap、98.3333%区间。每family是三个固定训练模型的平均**损失**，不是ensemble图。

| 第二批事前冻结的主比较 | 相对风险降低 | 98.333% 配对区间 | 实际裁决 |
|---|---|---|---|
| 完整灵敏度＋回投影 / 仅残差回投影图 | 3.24% | [2.13%, 4.34%] | CANDIDATE_LOWER_RISK |
| 完整灵敏度＋回投影 / 仅训练谱灵敏度图 | 0.26% | [-0.23%, 0.75%] | UNRESOLVED |
| 仅训练谱灵敏度图 / 仅随机谱 decoder 灵敏度图 | 0.70% | [0.17%, 1.24%] | CANDIDATE_LOWER_RISK |

第二批估计均为raw/no-blend，不改变第一批投影融合主结果；两批不混作一个事前512人主检验。确认对象完成全部在线估计与物理子作业后才读取GT，全部256个对象及fallback进入分母。成本仍在原L3m/总20m预算中，不追加授权。

## 外推范围

仅验证同一合成P1 prior、原噪声和acquisition下的新parent。没有新噪声/新频率/新prior泛化，没有真实测量或三维vector-Maxwell外推。独立测试成功不消除模型误差或模型外材料质量。
