# Research A — 最终科研裁决

**保留径向多候选成像benchmark与测量交接；停止当前feedback/自适应复杂化与三维扩展。**

|问题|裁决|证据|
|---|---|---|
|包内代数与零物理接口|PASS|E0_PACKAGE_ALGEBRA、E0_INTERFACE、16项必要contract检查；只证明相应有限数值/接口|
|实际full-Mie与A30物理接口|PASS_WITH_RETAINED_NUMERICAL_AUDIT_FAILURE|Mie独立递推和n14/18；A30实际finite/Schur/dual；原微小残差相对误差失败保留，改可分辨审计残差后通过|
|径向找到不同相容材料世界|PASS_IN_DECLARED_SUPPORT|32新对象、193个固定编辑full验证材料候选；每对象不同，不能说同一观测有193个|
|同预算Best-of-4对普通多起点|POSITIVE_LIMITED_RESULT|归一化MAE.112949 vs.136339，下降17.16%；parent paired95%差[-.03578,-.01245]|
|冻结强算法验收的联合条件|FAIL|Best4条件通过，但Coverage@4仅+.0625<+.10；实际搜索成本下降约10.9%<20%；二者均未满足第二条件|
|在线单图选择|NO_IMPROVEMENT|Selected .152732 vs共同base .151770；低残差不能稳定选到更接近truth的图|
|历史自适应next-solution增量|NOT_SUPPORTED|Best4 .113789，未优于固定.112949，区间跨零且耗时更多|
|Maxwell-specific feedback机制增量|NOT_ESTABLISHED|相同保护子域的独立起点重现端点，最大材料差2.84e−5；未有TSOM/ROM公平成本胜出|
|A30二维困难条件|NO_GO_CURRENT_MECHANISM|新增≥6cell结构族0；5/8无合法解；其余3条件archive Best4全部方法相同；唯一额外材料已在无补偿臂出现|
|一般三维E3/GPU/大型训练|NOT_RUN_NOT_AUTHORIZED|本次没有相应授权；E2机制失败亦不支持提交启动|
|候选概率/后验校准|NOT_ASSIGNED|candidate_only，所有weights/posterior_weight为null；不以残差伪造质量|
|A→B交接|PASS_VERSIONED_COMMON_EXCHANGE|32真实原生execution-profile与commonexchange；不冒称已适配未读取的B原生schema|
|强IEEE TAP算法论文|INSUFFICIENT_EVIDENCE|有限径向改进主要来自结构分层，A30未越过反证，缺一般3D与独有机制|

E1进入E2的事前门槛是物理/接口一致、发现真实结构族及材料误差改善信号。它与最终强算法验收不同；本轮确实按该较弱探索入口进入，并保留了E2失败。没有事后放宽噪声、材料支持、最终强条件或历史Gate。
