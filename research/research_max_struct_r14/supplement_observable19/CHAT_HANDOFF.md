# G19 Observable19 补充交接

**本次已实际完成144次，不是计划。** 24个已曝光R14对象×G/S/O0三种19维右基×PGP/GP；120整进程CPU/72proposal、原G19/current/prior/box/corrector/原Phi选根，576材料不减，真实forward一直保留完整G。无新NN、材料、照明、三维、Fresnel或rank扫描。所有预测先封存，后GT评分，selector720/720一致、issues0。

|19维选基|PGP primary恢复 / mean NRMSE|GP primary恢复 / mean NRMSE|
|---|---:|---:|
|内部G|10/24 / .198433|10/24 / .197667|
|接收S|12/24 / .158919|12/24 / .161074|
|公开m0=.1的O0=S(I−Dm0G)^−1|12/24 / .159954|12/24 / .162209|

共同已付费CPU与LU/RHS候选前缀均为G10对S/O0的11，不能把所有净增2例都称排除预算时序后的增量。S-PGP新增0003/0012/0021/0022、丢失0009/0020；S-GP保留0020但不能恢复0021；O0没有再增加集合。全ID均r14_new_*。旧r13_fresh_0021另一个材料，旧R14新0021原赢家是11，不是G19。

新0021：G-PGP .499289/.201646，S-PGP .051057/.932540，O0-PGP .050784/.932540（NRMSE/IoU）。S-PGP raw已.102341/.851449，corrector进一步改善；GP仍失败。不能把全部效果归于corrector。原G的Phi下降却材料变差。

21共同G-PGP状态×六条件读出，126条no-GT OFFLINE_DIAGNOSTIC全完成；原读出复现m_after最大差0。早期观测端口遗漏更低，但后期G坏轨迹中不再如此，局部coverage/Phi下降不能预测完整恢复。21共同状态来自一条基准轨迹，不是完整PR checkpoint/所有方法因果动力学。

理论澄清：传播或一阶残差谱界≠双块材料闭环谱界；不同metric可产生扩张但不改变噪声可行集合。算法不需要全局加法性，不可加仅限制统一线性Jacobian解释；不要据此否定Möbius/PR。真正材料Jacobian含S、内部resolvent和D_Q源场，不能直接把current右基当材料可辨识基。

主CPU4.778小时、并行26.90分钟，98387物理LU/1913532复RHS；setup、小19系统、FFT均收费。诊断v1新增输出字段曾误引用state，21条G-PGP后续验收不完整，已修复并以v2重放；两版原数组/费用都保留，主144预测完全不变。见FAILURES与COST_LEDGER。

阅读顺序：RESEARCH_REPORT_ZH → THEORY_CLARIFICATION_ZH → EXPERIMENT_CONTRACT → analysis/observable19_v1逐例表图 → analysis/fixed_state_summary → audit/fixed_state_v2。复现和完整full-field分卷在REPRODUCE/evidence.html；历史R14保持冻结。

判断：观测选基是值得保留的简单正结果，尚不是新盲确认或普适最优算法；O0未证明额外复杂度值得。下一轮最值得解释PGP/GP为何交换0020与0021的可达性，以及早期材料切向/端口结构如何累计改变盆地。本次到此收尾，不自动开展新算法循环。
