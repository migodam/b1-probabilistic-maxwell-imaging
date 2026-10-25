# Evidence ledger

|证据|路径|资格|
|---|---|---|
|数学三单元/19扫描/7检查|math/A1_MATH_RESULTS.json|PASS，解析模型数值验证|
|积分/尾部/box三测试|math/UNIT_TEST_STDERR.log|PASS，不是物理验收|
|六parent/paired innovations|datasets/PARENTS_AND_NOISE.npz|StageA自生成；不涉及B truth|
|18实际观测|v2/datasets/|sigma每份明确，与parent聚类|
|18推断尝试汇总|OBSERVATION_SUMMARY.csv|逐case参考资格；全局HOLD|
|coupled mixed strata|v3_strata/*RESULT.json|正常区数值强、host仍有限制|
|host直接已付费节点readback|v3_strata/HOST_DIRECT_READBACK.json|zero new physics；不声称strict omitted-mass bound|
|regular GN/Laplace密度|v3_strata/*GAUSSIAN_BASELINES.json|仅收敛参考上的数值TV/KL|
|higher-field normal rank2|rank2_higher_field/*RESULT.json|P0三噪声；P1partial保留|
|实际QR-RTO|v2/,v3_strata/,remaining_spaces/*QR_RTO.json|local ideal map density；global normalization NOT_QUALIFIED|
|全成本|COST_LEDGER_FINAL.json,COST_BY_JOB.json|32705 states/291784RHS，含失败与父线程probe|
|独立审核|../audit/A_INDEPENDENT_REVIEW.json,../audit/A_V3_MIXED_RANK_REVIEW.json,../audit/A_V3_CLOSEOUT_QA.json|不替代数值包络|
|父线程曲率充分条件|../audit/A5_INDEPENDENT_REVIEW.json,../audit/curvature_bounds/RESULTS.json|六条件FAIL；scope fixed weak center|
|父线程接缝物理probe|../audit/seam_maxwell/RESULTS.json|实际Maxwell，slope2.996758；非连续PDE证书|

没有用离散quadrature measure的TV冒充连续posterior TV；没有用模型响应通过代替branch/full66校准。
