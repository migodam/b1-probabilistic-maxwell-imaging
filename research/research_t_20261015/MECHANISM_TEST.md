# 有限结构编辑与补偿

|轮次|基础规则／主要变化|raw通过|补偿后通过|独立parent有成功|
|---|---|---:|---:|---:|
|1|最低残差base＋普通全材料补偿|0/16|15/16|4/4|
|2|同提案＋非core正材料相保护|0/16|10/16|4/4|
|3|最小N3改动的合格base；固定两个编辑|0/8|7/8|4/4|


三轮raw均0通过，所有成功均满足原观测tau、全球被动材料盒、共享gain1、本构/场检查及core/guard。core13cell，物理半径2/24；守护层1cell。材料距离为体积RMS/固定1.2，phase差≥6cell、primary absχ≤.1，另报.075/.125敏感性。这些是离散几何定义，不是连续persistent homology。

v1 raw残差265.99–394.47，经补偿后通过端点为20.24–23.79内；单个失败27.691。v2有6失败；v3有1失败。所有失败保留完整数组和line-search/KKT记录。

[逐编辑成本与物理距离](PER_EDIT_METRICS.csv) · [same-observation配对材料/观测距离](PAIRWISE_VALID_MATERIALS.csv) · [图与离线误差](offline/index.html)。每条成功core的13个格点在原base都不是primary近空气；纠正器不能用填孔或开通气道冒充成功。

成功伴随全域补偿与额外相变化；材料质量没有提高。噪声球相容与低Bayesian evidence可能并存，不报告任何posterior probability。没有验证全部中间路径相容或证明disconnected branch。
