# R6 证据、合同与理论覆盖独立审计

2026-10-27。范围仅 `audit/r6_evidence/`；新 Maxwell states/RHS=0，训练=0。不作最终理论有效性或原创性裁决。

## 主要发现

1. **R5T 的 B0 状态已过时，但其历史原文必须保留。** R5T Source §§1.5、12.5、末段与 NEXT 把新256+256记为 NOT_RUN；实际 R5B 科学提交 `d6fb550001169569b829d3f4b96cf10f625068e5` 已完成512 independent parents。当前 sibling HEAD `3e61f963` 含后续机械QA。R6只能将其列为 `reused_evidence`，不能再算新确认、重新生成同512或把其真值当未暴露。
2. **独立风险增益不解除80% HOLD。** R5B 原冻结四比较 family、20,000 parent bootstrap、Bonferroni4、98.75% CI，seed先平均loss。retention=80.69595%，二侧区间[60.28254%,113.12494%]，80%线性contrast的单侧下界−5.430356e−5，故 HOLD。它不同于R5开发合同三比较/4000抽样，也不同于R5T建议三比较/20,000抽样；不能互换置信水平重新宣布成功。
3. **所有原 SELECTION_FROZEN 的绝对checkpoint/输入引用现均存在。** 共清点118,207个R5B `.npz/.pt/.jsonl`文件名与大小（118,183在raw），包括原始状态、在线contexts、封存预测和PARENT_LOSSES；未打开GT/重评分、未重算大文件身份。原checkpoint有自带BLAKE2b清单；存在不等于内容本次重新认证。
4. **R5T decoder三文件与ad71固定对象并非逐字节相同：各多一个末尾换行。** Python AST和两份JSON语义完全相等；附BLAKE2b与明确差异。不得把“exact committed”描述为原字节完全一致；这不是统计模型变化。
5. **解析与decoder回放成功，但仅为重复非物理证据。** 解析单元数值断言通过；decoder最高档host/pure对称估计0.0006880817773100669，1,146,880 Rao–Blackwell点，与输入结果一致。host derivative=0、pure rank=2。无posterior、Maxwell likelihood或normalizer计算。脚本首跑缺少threadpoolctl，错误保留；仅audit副本改为线程环境请求，未安装依赖/改历史。
6. **先验原子计算已限内存。** 原脚本单个65536×256 b数组约128MiB，超过本任务100MB批限制。audit副本按4096行计算b（8MiB）；Sobol u与z_rest单个最大各31.5MiB，其他小矩阵另计。未声称进程总RSS≤100MB或实测硬件峰值。
7. **文献查询受429阻断，全球新颖性仍未审。** 使用 scholarqa-research 的collect，三查询snippet+paper均429；间隔后仅重试一个失败paper请求（no snippets）仍429。无candidate ID、无可verify论文，无泛用web替代。不能将未检出当成无先例。
8. **旧负结果与费用边界保持。** R5 A5六门全false、覆盖NOT_CERTIFIED；宿主两档均值未过1%SD；RRR物理专属增量未建立；原CPU/MPS 3e−6门FAIL未抹去。R5B 1,999,504历史RHS单列，R6新增0；provider压缩不代表全系统明显加速（R5B约0.0117% anchor-inclusive节省）。

## 原统计与求解合同

- 66维Gaussianized坐标对应物理 prior Gaussian64×Uniform2，原quintic精确平台，24×24 TM/VIE、256ROI、8×32复观测。sigma=0.00119211798514519；每源32real后32imag；f=√2 pack(F)，Phi=||F−y||²/sigma²。不可重新缩放sigma或加材料Jacobian prior。
- R5 A三个二维模型(0,1)/(0,64)/(12,27)，其余64固定0，活动标准正态截断[-4,4]²；6 parents各三档共享noise，不是18独立parent，也不是full66边缘。
- R5 A门：TV≤.05、可算双向KL≤.1、branch≤.02、遗漏≤.01、均值≤.01 posterior SD；退化绝对floor1e−8；reference误差<比较效应1/5；扩大需要同资格成本至少降20%。当前失败未更改这些门。
- R5B Native `reconstruct(full_TRF,max_evaluations=96,polish=False)`，随后 `probit_optimize(cap=64)`；原ftol=xtol=1e−9、gtol=1e−7、exact TRF、x_scale=jac。512 Native finite endpoints；Probit success338，174保留原最佳有限点。success不等于驻点、全局逆或posterior资格。
- R5B先封两批全部512预测再offline GT评分。TEST已暴露；未来开发须新数据/版本。raw complex MSE分母576；clipped、单seed、ensemble为secondary。并未在本审计重新打开truth数组。

## P0的本次覆盖与缺口

|单元|本次证据|资格与剩余工作|
|---|---|---|
|标量cubic平台/seam|独立路径回放原analytical_units.py，断言通过|Q1数值重复；不是新定理或Maxwell|
|Full66先验原子|静态源码审阅+有界分块回放|decoder-only Q2诊断；不是概率包络或posterior质量|
|h^6、z+h^3、(hu)^3、固定负residual|Source §4.7全读，未在本任务另编程实现|这四个R6 P0反例仍需负责P0的worker完成，不能用本表记PASS|
|3D QD gauge、source积分与合并density|§3.3推导与范围限制已读|未运行新积分对比，NOT_RUN|
|branch normalizer/有符号分子/task遗漏界|§§5.6–5.7已读|未运行区间传播新单元，NOT_RUN|
|连续latent与材料原子一致性|§§2.2–2.4、8.2已读；host decoder witness复核|未运行完整posterior sampler；纯prior证据不可升级|

## 供主Codex复核的可疑或条件性主张

- §4.2幂律依赖统一非退化、可积支配、远分支一二阶矩控制与正平台质量；源报告已写这些条件。Full66小p_H不证明当前噪声已进入平台主导区，不能把标量指数当真实Maxwell定律。
- §4.4必须验证P_A⊥c；二维固定虚属性制造的新增秩不能外推到3D或Full66。§3.3的正面积overlap仅属于新3D模型。
- §4.3 Hausdorff密度还需导数余项；§4.6 log阈值是scalar regular-variation条件结果。报告已限范围；下一轮需避免用命名算法重述成普遍结论。
- §5.2 receiver-task候选界和§5.7任务遗漏界的公式意义不等于常数足够紧或证书廉价；点检查、区间包络和全局覆盖必须分列。§8.5经典成本加权Neyman分配仅在真实均值中心化的一阶ratio近似下最优；pilot plug-in并非精确有限样本最优。
- §6.5 action扰动风险恒等式包含未知teacher Bayes excess，不能用D²单独保证风险提升。§6.6任务加权谱压缩是线性head模型结果；非线性head须实测独立风险。
- §7.1/§8.4 Gaussian-channel covariance identity不自动适用于错误q。三矩/中心化标签须保留同fold baseline、无噪f或原innovation，误差÷sigma²另立门；小identity residual不能认证正确。
- “task-weighted Bayesian quadrature”应区分GP概率数值积分与当前分层adaptive/IS任务误差设计；不能仅凭都估计QoI而声称同算法或新Bayesian quadrature贡献。

## 交付物

`SOURCE_AND_ASSET_AUDIT.json`保存全部Source章节标题/行号覆盖、三源码差异、原绝对输入引用存在性及R5B raw目录统计；`THEORY_SECTION_COVERAGE.md`补充逐节用途；`CITATION_LEDGER.json`保存检索缺口与输入文献未验证状态；两份ScholarQA原始bundle保留429；`analytical_replay/`与`decoder_replay/`保存隔离副本及结果；`EXECUTION_FAILURES.json`保留首次运行依赖失败。
