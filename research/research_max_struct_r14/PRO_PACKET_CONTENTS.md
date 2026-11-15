# 给Pro的最小curated证据包建议

这是curated审阅包的内容规范，最终实际成员、bytes与生成状态以公开PRO_PACKET_MANIFEST为准。它不替代完整raw release，也不把打包/上传当科学验收。目录名、任务数和二进制文件大小来自实际文件；正文/论文以主线程定稿为准。后续版本另纳入新0012/0021的01/11各128/512/1500完整state共12个原档，供研究早期metric为何选择相反材料区域。

## 1. 首先放入的阅读顺序

1. RESEARCH_REPORT_ZH.md、THEORY_DERIVATIONS.md、THEORY_STATUS.md：给出主问题、实际数学与已证明/未证明层级。主报告若尚未完成，不制作假的完成版。
2. METHODS_AND_ABLATIONS.md、BASELINE_AUDIT.md、PRIOR_ART_AND_NOVELTY.md、CONTROL_DERIVATION.md：区分四度量、readout、PR memory、stage与统一预算；先例以已核原文/证据账为准。
3. TRAJECTORY_ANALYSIS.md、BENCHMARKS.md、COST_LEDGER.md、FAILURES.md、REPRODUCE.md：存在时包含；不以缺失报告作隐含结果。PRO_HANDOFF/CHAT_HANDOFF只收作者整理稿，不放原聊天缓存。
4. 本文件及一页主线程自己写的提问：建议要求Pro逐条指出哪个命题由哪个证据支持、什么替代解释尚未排除，不要求先肯定论文价值。

## 2. 恢复指标只放紧凑且足够核对的表

|证据组|最小文件集合|
|---|---|
|主暴露48|analysis/exposed48_v2/{RECEIPT.json,SUMMARY.csv,ALL_VIEWS.csv,PAIRED_METRIC_DIFFERENCES.csv,RECOVERY_GAIN_LOSS.csv,COMMON_CPU_PREFIXES.csv,CORRECTOR_EFFECTS.csv,RAW_INITIAL_AUDIT.csv,STRICT48_PAID_POOL_BRACKETS.csv,COST_LEDGER.csv}|
|独立24/264任务|analysis/independent24_v1/{RECEIPT.json,SUMMARY.csv,ALL_VIEWS.csv,POSTHOC_GAIN_LOSS.csv,COMMON_CPU_PREFIXES.csv,DIRECT_PHYSICAL_SEARCH_VIEWS.csv,PUBLIC_SELECTOR_REPLAY.csv,PUBLIC_SELECTOR_COUNTS.json,RAW_INITIAL_AUDIT.csv,COST_LEDGER.csv,ISSUES.json}|
|机制控制|analysis/self_only16_v1与stage_control16_v1的RECEIPT.json、SUMMARY.csv、ALL_VIEWS.csv、COST_LEDGER.csv；shared_readout4_v1的SCORING_RECEIPT.json及实际summary/paired表|
|3D类比|analysis/vector4_v1/SCORING_RECEIPT.json、METHOD_VIEW_SUMMARY.csv、SCORES.csv、CORRECTOR_DECOMPOSITION.csv、COST_LEDGER.csv、audit/VECTOR_CONTRACT.md|

对应冻结configs、公共选择规则、授权seal/authorization的**公开整理版本**一起提供，保留其预测先于评分的出处。原Windows seal不要改成Mac路径伪装原件；需要本地映射时附独立迁移说明。264是最终11方法×24，不使用早期216草案。暴露main192、self48、stage32、shared16、vector24分别列，不混成独立样本总数。

不放NATIVE_STAGE_HISTORY约1.8GB JSON/1.1GB CSV、全量NATIVE_STEPS、PAID_CANDIDATE_HISTORY等大表进curated包；这些原件留完整release，并在内容索引给出真实asset/member路径。不要将丢弃原大表描述成“只需要小结果”的科学简化。

## 3. 七个完整checkpoint：不能换成m-only

固定原始路径位于RESULTS/dualmetric_pilot12_v1，均为pair01历史完整状态：

|case|checkpoint路径尾部|实际bytes|
|---|---|---:|
|r13_fresh_0002|r13_fresh_0002__pair01/CHECKPOINT_00128.npz|870274|
|r13_fresh_0002|r13_fresh_0002__pair01/CHECKPOINT_00512.npz|1284208|
|r13_fresh_0017|r13_fresh_0017__pair01/CHECKPOINT_00128.npz|851391|
|r13_fresh_0017|r13_fresh_0017__pair01/CHECKPOINT_00512.npz|1121377|
|r13_fresh_0017|r13_fresh_0017__pair01/CHECKPOINT_02500.npz|2213283|
|r13_fresh_0037|r13_fresh_0037__pair01/CHECKPOINT_00128.npz|879061|
|r13_fresh_0037|r13_fresh_0037__pair01/CHECKPOINT_00512.npz|1293101|

七个合计8,512,695bytes。保存原完整NPZ，不压精度、不删J/Q/A/old_grad/old_direction/表示或stage cache。再加入3个对应data_public/exposed48/r13_fresh_XXXX.npz（只有公开G/E/S/Y等），以及src/dualmetric_core.py、intervene_checkpoints.py和它们实际依赖的历史源版本说明。完整checkpoint是科学输入，不能拿export_scoring产生的selected-m NPZ替代。

七个干预各自RECEIPT.json、BRANCHES.json、CONFIGURATION.json、COST_LEDGER.json、FAILURES.json，以及audit/dynamics_all192_v1/INTERVENTION_BRANCHES.json和INTERVENTION_RECEIPTS.csv一并收录，便于核56分支与1/8/32步行为。想最小化时不重复收整套7×8×全部step字段；完整分支状态仍在fullstate release。若Pro要独立重算force，必须再取得下面对应56对step0/step1完整状态；原7checkpoint只允许重新执行干预，不等于已有全部force输入。

## 4. force数组与来源

必须收：audit/material_force_decomposition_v1/FIRST_STEP_FORCE_ARRAYS.npz，实际20,501,012bytes。源码表明它对每个7×8分支保存m0、m1、J1、Q1、e、Pe、Os、delta_m、g_const、g_metric、g_prior、g_full，足以直接查看已有分解和实际步；不读取任何GT。

同时放入该目录的：FIRST_STEP_MATERIAL_FORCES.csv、MECHANICAL_GROUP_SUMMARY.csv、MECHANICAL_BRIEF.json、PROTOCOL_BEFORE_ACTIONS.json、COUNTER_AND_SOURCE_SCOPE.json、VALIDATION_COST_EVENTS.csv、RECEIPT.json、FAILURES.json、SOURCE_AT_ACTIONS.py、REPORT_ZH.md。保留原g=gradient约定，下降force是−g；点积符号不改写成恢复收益。

七原checkpoint与这一force NPZ合计29,013,707bytes，另加公开输入、紧凑表和源代码仍是有限小包；这只是基于已有字节数的规模估计，不冒充已经生成的ZIP大小。PNG优先MATERIAL_FORCE_COMPONENTS和主线程最后选定的实际图；不要重复放FIRST_LAYOUT旧布局当最终展示。

## 5. 根与有限程序的最小数学反例材料

若此次给Pro的问题包含“理想T能否解释finitePR”，追加audit/local_closed_loop_v1/{SUMMARY.json,TABLES.csv,ANALYSIS_NOTES.md,ROOT_FULL_STATE.npz,ROOT_PUBLIC_AND_PHYSICAL.npz,DIRECTIONS.npz,IDEAL_00.npz,IDEAL_01.npz,IDEAL_10.npz,IDEAL_11.npz}，以及实际RESPONSE/ADDITIVITY数组和src/local_closed_loop.py。

加finite_mode_trajectory_v1与coupled_mode_trajectory_v1各自的FROZEN_CONFIG.json、MODE_SELECTION.json、TRAJECTORIES.csv、FIRST_STEP_POLICY_EQUALITY.json、RECEIPT.json、COST_LEDGER.json、FAILURES.json、RESULTS_BRIEF.md、INTERPRETATION_BOUNDARIES.md；完整逐步轨迹fullstate仍放大release。用源程序声明的同一mode/epsilon/PR策略，不让Pro把方向差分自动当统一Fréchet Jacobian，也不把辅助current memory当新增物理可观测信息。

## 6. 图、引用与不应混入的材料

主材料图采用FIGURES/exposed48/materials_01–08.png及all48_readable.pdf，保留完整48对象。独立24的实际FIGURES需要完整阅读时由主线程决定是否全部纳入，不只挑成功case。3D说明只配真实central-slice图，承认不是每voxel视觉审查。旧analysis/exposed48_v1 bugged导出不进主包，仅保留FAILURES中的出处。

文献只放audit/prior_art/REVIEW.md、CLAIM_EVIDENCE.csv、SEARCH_LOG.md、已核书目及作者/DOI链接，不复制第三方全文PDF/HTML。不要包含private_operational、SSH路径、token、sources/CHAT_READBACK、原聊天缓存、密封真值或生成器私有材料标签。这里的scientific query不需要暴露权限/机器操作日志。

最终包由主线程逐项核现存文件、定稿与来源隔离后制作；另附成员路径/bytes清单和真实full-release恢复入口。本文是内容规范，实际发布以DOWNLOADS及匿名readback为准；不因制作一个包新增科学结论。
