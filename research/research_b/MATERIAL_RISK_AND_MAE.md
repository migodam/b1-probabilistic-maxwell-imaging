# 材料风险、实际MAE与费用

| 策略 | E1 固定625池选择 MAE | E2 固定625池选择 MAE | E2 有限参考辅助选择 MAE† | E2 新增blocks均值 |
| --- | --- | --- | --- | --- |
| none | 0.162279 | 0.162279 | 0.104046 | 0.000000 |
| fixed_1.2 | 0.149763 | 0.140054 | 0.085696 | 3.000000 |
| fixed_1.6 | 0.143048 | 0.133089 | 0.066354 | 3.000000 |
| fixed_2 | 0.150269 | 0.114810 | 0.062181 | 3.000000 |
| fixed_2.4 | 0.122534 | 0.106367 | 0.059581 | 3.000000 |
| fixed_2.8 | 0.123681 | 0.107033 | 0.052725 | 3.000000 |
| fixed_3.2 | 0.101982 | 0.090885 | 0.039812 | 3.000000 |
| random | 0.133483 | 0.105370 | 0.055309 | 3.000000 |
| uniform | 0.149763 | 0.097503 | 0.053684 | 3.000000 |
| raw | 0.102027 | 0.093068 | 0.043414 | 3.000000 |
| classical_Fisher | 0.101982 | 0.090885 | 0.039812 | 3.000000 |
| joint_projected_Fisher | 0.101982 | 0.089865 | 0.040140 | 3.000000 |
| event_EIG | 0.104847 | 0.090107 | 0.042696 | 3.000000 |
| expected_material_risk | 0.101982 | 0.093775 | 0.040103 | 3.000000 |
| robust | 0.103714 | 0.093842 | 0.044747 | 3.000000 |
| conflict_diagnostic | 0.141574 | 0.113446 | 0.049028 | 3.000000 |
| calibration_only | 0.141574 | 0.134254 | 0.084007 | 3.000000 |
| calibrate_then_risk | 0.141574 | 0.094999 | 0.038684 | 3.000000 |
| geometry_only | 0.130434 | 0.125114 | 0.070882 | 3.000000 |
| risk_cost_stop | 未运行 | 0.093687 | 0.040103 | 2.984375 |
| robust_abstain | 未运行 | 0.094046 | 0.045372 | 2.984375 |
| two_step_bundle | 未运行 | 0.095791 | 0.047853 | 2.000000 |

| 选择池/信息来源 | 策略 | MAE差：策略−固定3.2 | parent配对95%区间 |
| --- | --- | --- | --- |
| mode1_candidate_only | event_EIG | -0.000778 | [-0.005217, +0.003528] |
| mode1_candidate_only | expected_material_risk | +0.002890 | [-0.000581, +0.007754] |
| mode1_candidate_only | robust | +0.002957 | [-0.002879, +0.010341] |
| mode2_reference_selection | event_EIG | +0.002884 | [+0.000838, +0.005081] |
| mode2_reference_selection | expected_material_risk | +0.000291 | [-0.001564, +0.002216] |
| mode2_reference_selection | robust | +0.004935 | [+0.000635, +0.010134] |

Mode1字段指固定625池的profile选择；event EIG/material-risk行的选测仍依赖外部参考。Mode2有限连续proxy参考未通过新history的完整积分收敛，不能称upper bound。raw full-Mie truth评分与参考概率可靠性分别判断。

Material MAE是精确体积加权Re/Im ε分量L1和、s_ref1。Selected-MAE不是Best-of-K；固定池后者均值.0622242099、Coverage@K≤.05为.40625，选测不能改变它。结构/事件正确率另存CSV。材料loss是决策任务，本轮没有引入其他工业下游损失。

32parent是统计单位；先平均两noise，再固定2000次bootstrap，不把64episodes当64对象。E2同E1人群且已曝光；四例A不与之合并。block成本不是秒/能量，当前没有硬件成本最优性主张。

## 真实A候选的独立四例

| 真实A交接方法 | Selected-MAE | Best-of-K MAE | Coverage@K ≤.05 | 平均K |
| --- | --- | --- | --- | --- |
| A_base_only_no_measurement | 0.125480 | 0.125480 | 0% | 1.00 |
| A_fixed_candidate_profile_no_measurement | 0.125480 | 0.110881 | 0% | 6.75 |
| A_robust_action_then_profile | 0.119000 | 0.110881 | 0% | 6.75 |
| A_fixed3.2_then_profile | 0.146238 | 0.110881 | 0% | 6.75 |
| A_reference_action_then_profile | 0.146238 | 0.110881 | 0% | 6.75 |
| A_robust_then_protected_eps_corrector | 0.064836 | 0.057371 | 50% | 11.75 |
| A_robust_then_independent_prior_multistart_corrector | 0.036170 | 0.031589 | 100% | 9.75 |
| declared_continuous_P2_reference_comparator | 0.027106 | 0.019557 | 100% | 96.00 |
| A_base_only_same_action_eps_corrector | 0.089090 | 0.078574 | 0% | 2.00 |

protected/prior-corrector只安装双频fullvalidated且数据兼容的新端点，才允许改变Best-of-K。原池始终保留。四例不足以证明policy优势，也缺同容量无新测量corrector对照。参考comparator未验证完整history，不安装材料。

逐状态原表为results/E1_POLICY_RESULTS.json、E2_POLICY_RESULTS.json、E3/RESULTS.json；逐parent与配对表见results/PARENT_METRICS.csv和PAIRED_COMPARISONS.csv。全部随Release提供。
