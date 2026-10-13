# B1 实验结果

2026-10-13。原始、派生与经验近似结果分列，不是严格概率认证。

|案例|实心%|单空腔%|双空腔%|原生5/7参考|经验细化|
|---|---:|---:|---:|---|---|
|distinguishable|100|7.98633e-10|4.83672e-07|REFERENCE_UNRESOLVED|EMPIRICALLY_REFINED_REFERENCE|
|competitive|8.93166|50.6418|40.4265|NUMERICAL_REFERENCE_PASS|EMPIRICALLY_REFINED_REFERENCE|
|information_limited|34.5821|32.3524|33.0655|NUMERICAL_REFERENCE_PASS|EMPIRICALLY_REFINED_REFERENCE|
|gain_mismatch|4.83998|2.66395|92.4961|REFERENCE_UNRESOLVED|EMPIRICALLY_REFINED_REFERENCE|
|gain_marginalized|3.3315e-07|49.7319|50.2681|REFERENCE_UNRESOLVED|EMPIRICALLY_REFINED_REFERENCE|
|representation_mismatch|1.21302e-26|83.5325|16.4675|REFERENCE_UNRESOLVED|EMPIRICALLY_REFINED_REFERENCE|

表中概率来自第二轮，百分比条件于三个模板世界；原生P数据逐项在CASE_PROBABILITIES.csv。Top2不重归一化，第三event剩余保留；模板之外的更大物理世界未量化。

## 同数据/名义冷成本对照

|方法|名义冷全波状态|Brier|log score|平均最大P差|完整P/logBF门槛比例|
|---|---:|---:|---:|---:|---:|
|gauss3|27|0.518486|0.852173|0.00314929|0.988889|
|gauss5|75|0.518211|0.851616|3.43196e-05|1|
|gauss7|147|0.518204|0.851607|0|1|
|frozen_importance48|75|0.530183|0.871099|0.0387663|0.172222|
|grid7_min_residual_softmax_HEURISTIC|147|0.530219|0.878695|0.0733014|∞ / 非概率点决策|
|grid7_hard_MAP_ONE_HOT_DECISION|147|0.883333|∞ / 非概率点决策|0.460127|∞ / 非概率点决策|

参考Gauss7，在24parent×30noise的已暴露开发panel上；固定importance的proposal来自一个y0，没有逐数据适应或重抽样。对照setup名义状态从已支付结果重算，实际本次新增全波0；Gauss5冷75，最初pipeline新72+复用3中心。不可把该表写成一般sampler优势或冷wall实测。

## 校准与材料不确定性

Brier=0.518204，log=0.851607，Top2真H覆盖=0.865278，平均报告mass=0.866522。

独立review的parent分层bootstrap、边际/真实H条件coverage与分箱偏差见CALIBRATION_INTERVAL_REANALYSIS_NATIVE.json。部分分箱区间不含0，不作校准PASS；rankSBC未运行。quantile仍是prior-mass-bin近似。

## 物理与成本

|任务|新材料–频率状态|
|---|---:|
|AUDIT_A13_HOMOGENEOUS|1|
|AUDIT_LOSSLESS|1|
|AUDIT_MAIN|3|
|AUDIT_RAYLEIGH|2|
|AUDIT_RECURSIVE|3|
|AUDIT_REDUNDANT_INTERFACE|1|
|AUDIT_TRUNCATION|2|
|AUDIT_VACUUM|3|
|DEFENSIVE_IMPORTANCE|48|
|EXISTING_A15_TREAMS_INDEPENDENT|4|
|FRESH_CONTINUOUS_CALIBRATION_PARENT|24|
|FRESH_OBSERVATION_TRUTH|3|
|MULTISTART_SELECTED_FULL_VERIFICATION|9|
|OUTSIDE_TEMPLATE_NEGATIVE_CONTROL|1|
|PAID_LINEAR_DERIVATIVE|12|
|QUADRATURE_ORDER_3|27|
|QUADRATURE_ORDER_5|72|
|QUADRATURE_ORDER_7|144|
|ROUND2_FORWARD_INTERPOLATION_AUDIT|24|
|合计|384/384|

真实forward包围计时之和2.266938s，包含treams子进程；它嵌套于job wall，不可重复相加。GPU0，NN0，旧locked/truth0。历史原包2D74+本阶段3D384=458；A31既有微型物理另列。

最终11个wrapper jobs计时sum wall=14.467812s，sumCPU=10.262424s。最终文件检查另记stdlib代码计时（不含解释器启动），未运行新物理或数值实验。native视觉阶段wall19.845079s/CPU18.896490s；早期renderer内部wall5.3861s/CPU5.2974s另保留；native数值review3.152s/CPU2.802s按自身receipt。网络、模型推理、交互读写与工具等待未完整测量，不称免费或全部研究wall。

## 失败/未运行

原生四案例reference不足、单例importance误差.02565、冻结panelimportance两容差仅124/720、gain/模板错配PPC异常都保留。两个review问题（coverage标签、gain NPZ元数据）已修复；第一次修正误用系统解释器的失败与错误CLI重放亦保留。不存在新物理失败，但实现失败不是零。

原生A17反馈、未知位置/边界、复杂一般形状、uniform概率认证、rankSBC、锁定测试、NN、公开发布与HTML现场浏览器验收均未运行。
