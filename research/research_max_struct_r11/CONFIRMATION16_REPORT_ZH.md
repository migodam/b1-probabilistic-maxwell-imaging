# 16个独立高对比/混合对象：修正的确认结果

**判断：两块Green修正不应替代原Green材料分母。** 它在旧桥上的正例没有在本组形成增量，甚至破坏了一个原Green能恢复的新环。原材料-only Green在本组相对CIE保留了恢复增量，但有实际成本代价；这也不是要求它在所有对象上都胜出。

本组16材料×5臂80次全部完成；预测、组合选根、共同LU预算选择在GT读取前封存。0执行失败、0重试。成功与失败均展示；不能把“执行成功”当“材料恢复成功”。48对象主比较与真三维研究仍在运行，本页只代表这一组。

|方法|材料成功|数据CN99合格|平均幅度NRMSE|中位NRMSE|平均实际CPU秒/对象|
|---|---:|---:|---:|---:|---:|
|CIE|4/16|12/16|.23134|.11219|47.96|
|Dual-CIE|5/16|13/16|.21068|.10889|48.57|
|Global-Dual控制|6/16|13/16|.21344|.10855|47.39|
|原Green材料分母|7/16|14/16|.14485|.10358|60.85|
|两块Green修正|5/16|13/16|.17357|.10891|62.84|

这批所有方法实际自然结束都在240 CPU秒以内，480秒只是统一允许上限。原Green实际CPU总和为CIE的1.269倍；两块修正为1.310倍。共同已达到的全波LU预算下，成功数仍为4/5/6/7/5，但物理LU相同不代表辅助运算免费。

## 哪些变化有实质意义

高对比环0000是明确反例：CIE NRMSE .9364，原Green .1003且IoU .7237，两块修正却变为.5495和.2968。既有开发环、开发桥的正例不能用来遮掉这个新失败。

高对比开缝壳0002中，原Green/两块的NRMSE .1068/.1089很接近，IoU .7022/.6649分居冻结门两侧；不能把这次成功标签变化夸大成和环同等严重的崩溃。混合开缝壳0003则CIE失败，Dual、Global、原Green及两块均恢复，说明这一例收益不是Green独有。

两块对原Green的paired平均NRMSE差为+.02873，10000次high/mixed分层对象bootstrap 95%区间[+.000176,+.084961]，5胜11负。主比较支持拒绝“普遍匹配两块就能修复”的说法。多个方法/子组的探索比较不作未经校正的普遍显著性宣传；仅16对象的区间也不代表所有几何总体。

原Green相对CIE增加3个成功，相对Dual增加2个；增量中包含大幅环恢复，也包含边界IoU变化。应同时看逐对象连续误差和完整图，不能只把7对4当作全面优势。原Green仍有高对比包覆体、桥及细结构共同失败。

同global current的普通材料控制已6/16；G19对它唯一新增的成功是高对比环0000，无反向成功回退。CIE→Global解释了另两个成功增量，不能归给分母。该环的实际CPU为CIE73.70、Dual64.79、Global47.40、G1954.48秒：相对CIE/Dual的环恢复并非额外CPU换来，相对Global仍有约15%的成本增量。

## 两种付费组合也没有解决问题

预声明CIE+Dual与CIE+BothBlocks都实际支付两个候选，再按同一原Phi选根，结果均5/16。前者部署CPU总和1544.51秒，后者1772.86秒；后者约多14.8%，未带来成功数增量。两组合共用了实验中的CIE运行，研究账本不会将其重复相加。

由于全臂自然结束早于240秒，240+240的同允许预算视图等于完整终点；它不是同实际花费证明。组合也没有额外信息或GT选根。

## 补齐四格消融后的解释

三个已曝光对象的Green-current/ordinary-material控制已另行完成。环/壳/桥NRMSE为.4220/.5151/.5546，全部失败；其中桥的两块Green则.1009。因此旧桥修复不是“current单独更好”，而是两个非线性更新块的路径交互。完整表见[四格开发报告](delegated/diagnostics/current_only_dev/REVIEW_ZH.md)。

这是一条有内容的负结论：即使两块使用同一Green残差族，也不产生全局恢复保证；当前双β仍无共同目标。单块凸化、较好current条件数及更多Green能量都不足以预测材料恢复。它没有否定原Green在特定环上的实际价值。

另外，本组12–14个数据CN99合格而只有4–7个材料成功，且多种方法留下很相似的平滑/缺失细结构。全16×5的零新FW审计已经完成：65个候选和GT都在实际y的CN99球内，其中38个材料FAIL、涉及8个独立材料。这65个候选同时比GT有更低数据项和更低prior，不能只责怪prior。另15个候选在球外而GT在球内，仍需解释算法路径/预算问题。38对预测的sigma距离为4.30–6.09，并非接近零；噪声球共同容纳不能称Gaussian统计不可区分或精确多解。详见[完整审计](delegated/diagnostics/bothblocks_confirmation/finite_noise_prior_audit/FACTUAL_SUMMARY_ZH.md)。

## 已完成的预算与停止补充诊断

公开成本/停止状态规则在32个CIE或Dual对象对中触发7次续算，用同对象Green实际CPU减去原基线CPU作为新增上限。新增62.857739 CPU秒、7370个576阶LU、176880 complex RHS；预测先封存后评分。CIE/Dual仍4/16和5/16，mean NRMSE仅从.231343/.210683变为.231240/.210673。4次触cap、3次提前PGTOL结束；不空转填满预算，也不把晚于cap的候选倒填。原子及序列化尾部.789021秒照常计费。它支持该范围内的增量并非仅由这些基线少用的corrector时间解释，但不证明任意预算下的优势。

这不是无缝延长：L-BFGS历史重置，原已PGTOL或花费已更高的对保留；它是POSTHOC_BUDGET_CONTROL，主80结果不改。完整选择及成本见[预算审计](delegated/diagnostics/bothblocks_confirmation/posthoc_budget_control/FACTUAL_SUMMARY_ZH.md)。

首次公开CN99命中的65个保留状态，材料NRMSE全部较最终差，材料成功0改善、19回退。读取已评分轨迹的GT oracle边界检查也未发现任何先满足材料成功后在最终丢失的运行。两者只针对已保存Phi-incumbents，不能排除未保留状态或其他停止规则，但没有支持用当前噪声门早停来修复。见[停止诊断](delegated/diagnostics/bothblocks_confirmation/first_cn99_hit/FACTUAL_SUMMARY_ZH.md)。

## 全部图像与原始事实

![全部8个high材料](FIGURES/bothblocks_confirmation/MATERIALS_HIGH.png)

![全部8个mixed材料](FIGURES/bothblocks_confirmation/MATERIALS_MIXED.png)

- [全部80条评分](delegated/diagnostics/bothblocks_confirmation/SINGLE_ARM_SCORES.csv)
- [事实与成本汇总](delegated/diagnostics/bothblocks_confirmation/FACTUAL_SUMMARY_ZH.md)
- [配对不确定性](delegated/diagnostics/bothblocks_confirmation/PAIRED_UNCERTAINTY.csv)
- [实际CPU轨迹](FIGURES/bothblocks_confirmation/PAID_CPU_FRONTIERS.png)
- [全波LU轨迹](FIGURES/bothblocks_confirmation/PAID_LU_FRONTIERS.png)

80次反演精确child CPU4281.798299秒，队列wall2182.972411秒；生成1.372750 CPU秒、controller .448456秒，评分/绘图77.417秒另计。current-only三例DEV128.297977秒另计，不混入确认成功数。
