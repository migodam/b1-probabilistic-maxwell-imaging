# Stack tables independent arithmetic and provenance audit

检查既有表格、公开候选索引和运行回执，未调用评分器/物理函数，未读NPZ或GT，未改算法或既有评分代码。检查项148，未通过0。

## 分母、索引与版本

- Uniform48为48对象×6方法=288运行；固定every-third16为16×6×uniform/NN，其中uniform96已属于288。新增NN96加optional NNdirect16，总唯一400运行；不是400个材料。
- 五视图2000行，五prefix2000行；576共同成本行=uniform_six288+NN_six96+factorial_twelve192。所有行available，0科学失败。4576候选ID与先落盘公开索引一致；2000视图ID另与源receipt/initial/corrected交叉核对。
- 实际回执128份A+核对：uniform48+NN16的64Arc均arc roots和gtol checks为正、无quintic roots；64MM均quintic roots为正、无arc/guard计数。精确逐回执来源见VERSION_COUNTER_AUDIT。Arc是R7最终readout选项放进R13协议，不能称完整R7实验原样重跑。

## Uniform48实际端点与新六方法共同成本

|solver|primary NRMSE均值|primary IoU均值|球内/48|平均whole-child CPU秒|common-six NRMSE均值|
|---|---:|---:|---:|---:|---:|
|cie|0.19453592|0.69455730|35|106.147135|0.19468446|
|dual_cie|0.17984533|0.71591743|37|106.263997|0.17978262|
|g19|0.16758204|0.68576885|39|119.083333|0.18221556|
|dr|0.17671934|0.68726032|39|84.496094|0.17671944|
|a_plus_mm_r13|0.15210017|0.75062659|40|119.367188|0.15745965|
|a_plus_arc_r7|0.15214123|0.75152546|40|113.222656|0.15281554|

原输出Arc−MM均值差4.106108725e-05，描述bootstrap95区间[-0.0003739628386,0.0004885523854]。这不证明未来等价或一般优越；当前均值近似而费用不同。球内数不是材料成功数。

## NN同16及corrector（端点描述）

|solver|uniform16 primary|NN16 primary|NN−uniform|NN raw→corrected|NN实际CPU均值|
|---|---:|---:|---:|---:|---:|
|cie|0.14919834|0.14980489|+0.00060655|0.15955411→0.14980489|101.773438|
|dual_cie|0.12786475|0.12818136|+0.00031660|0.14100461→0.12818136|102.695312|
|g19|0.15007514|0.12655605|-0.02351909|0.18831921→0.12655605|120.596680|
|dr|0.12312718|0.11935740|-0.00376978|0.14472286→0.11935740|83.851562|
|a_plus_mm_r13|0.10966243|0.10943134|-0.00023109|0.12686486→0.10943134|120.579102|
|a_plus_arc_r7|0.10970313|0.10996435|+0.00026121|0.12675307→0.10996435|108.638672|

## 必须保留的derived表边界

1. **NN_six与uniform_six预算不同，不能跨组相配。** 每个对象各自取其6个NN任务、6个uniform任务的实际费用最小值。当前16/16对象的这两个limit实际均不同（由PUBLIC_SIX_AND_TWELVE_CPU_LIMITS逐对象核对）。NN_MINUS_UNIFORM中的96个NN_six paired_available=false是结构性“不作比较”，不是缺失结果或科学失败。只有factorial_twelve使用该对象全部12运行的共同limit，提供同成本NN−uniform配对（96个可用配对）。标准primary/prefix的同对象比较另外按其原预算语义读。
2. Optional NNdirect16没有uniform-direct纳入stack400，因此5视图80+5prefix80的160个配对均明确不可用。原数据并非缺失：本审计另从core48 direct限制同16，列于DIRECT16_REFERENCE，完全不改stack评分表。
3. 已记录的cohort label correction确实仅从旧36条共同成本summary中删除6条错误uniform48标签；当前30条均在备份中逐字段原样存在。被删的是factorial_twelve16的额外uniform48副本；预测、候选、误差和其余汇总数字未变。当前所有uniform48标签均为48。
4. Corrector差表的whole-job CPU/LU差为同一run视图相减，可为0；不是corrector阶段免费。阶段候选累计CPU/LU/RHS另列。raw→corrected指标有好有坏；不得将Phi下降等同材料恢复改善。

## 费用与复用

新增stack112实际whole-child CPU=12373.671875秒、物理LU96749、complexRHS1773918；旧288任务是已付运行引用，CPU30590.515625秒、LU253913、RHS4818594。五视图、prefix、六/十二方法表重复引用同一run，不能重复加总。

R8主physics NN训练加生成历史费用402.558889 Mac CPU秒已保留；本次未训练或推断，既有R13预处理另计，不能当作零成本部署。

## 历史global-PR Dual与direct控制

GLOBAL_DUAL16_REFERENCE.csv/.md包含原mechanism16 global_dual和A+的96条raw/corrected/primary行；DIRECT16_REFERENCE.csv包含同16 uniform/NNdirect32条原primary行。每行原索引和candidate均与原receipt核对。0045 Global校正材料恶化与A+不同raw起点的改善均保留。此控制不是新算法、重选候选或共同实际CPU比较。

## 审计结论范围

本审计确认表格分母、候选索引、版本分支、配对差算术和唯一费用复用的内部一致性；没有核验新成像主张、数学因果、新颖性或普遍性能。完整逐项结果见TABLE_AUDIT_CHECKS.csv；无算法修改。
