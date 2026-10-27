# R6：给下一轮GPT Pro的实证交接

先读[完整报告](../RESEARCH_REPORT_ZH.md)、[逐路线裁决](../GO_NO_GO.md)、[数学与条件](../THEORY_AND_COUNTEREXAMPLES.md)、[机器结果](../RESULTS_SUMMARY.json)。本轮真实消耗6,860,984 RHS，全部物理/训练关闭；不是未来计划。

## 最值得推理的新现象

1. **同材料gauge可变成伪材料不确定性。** S3固定材料V=QD的正确纤维上，材料完全不变；即使精确知道Q、D边缘，独立组合仍产生均值偏差V(EQ·E1/Q−1)>0。单次32对置乱的两例Full Maxwell归一化复观测差均值15.53/117.58，实际448 RHS。边缘乘积积分与单次置乱分开，非Full66 posterior样本。
2. **零切向Jacobian没有涵盖所有Bayes任务响应。** host平台anchor的latent J为零，但有限分层后验的材料均值响应非零，between-source项约占导数范数28%。三档固定创新支持条件尺度；没有全域tail证书，不能升级为Full66普适定律。
3. **优化阶段的目标不相同。** Native的Uniform属性点目标常推向box边界；Probit另有Gaussianized坐标代价。花更多Native费用没有减少后续Probit工作。跳过Native的普通Probit很有竞争力，复杂stage allocation并无已证明的独立价值。
4. **KKT闭合在低噪下有独立风险收益，但不是普遍胜利。** 残差反投影用正号q代替的规则无需新权重。512新对象的原global8均支持成本/风险门；KKT机制仅sigma/2成立。强化Full-q参考后，低噪1% NI仍HOLD。
5. **计时漂移能夸大“简单基线更快”。** 相同端点重放零差，但不同窗口绝对wall明显不同。同窗口Direct对Low没有确定优势。原窗口通过不等于普遍跨窗口20%保证。

## 已经充分尝试而应停止扩大

N1024四目标×3seed、12provider+12head，加C的9个matched head，共33个真实训练阶段。新的action/GT目标未击败历史N32768部署基准。任务交叉矩12个ridge试验暴露高SNR回归一致性问题。不要把它改写成NN无法学习，也不要重跑旧R5-B512来证明同一结果。

## 资格

A只有2D/3D截断条件模型、局部Q1/Q2；任务成本分离和4D没有建立。Full66材料posterior仍未认证。新D只验证两档noise，不验证新layout。旧≥80%retention仍HOLD。当前TAP判定REDIRECT，不能仅凭省时或已知恒等式宣布新论文成立。

详见[A证据](A_BRANCH_EVIDENCE_NOTE.md)、[下一轮三个问题](NEXT_THEORY_QUESTIONS.md)。
