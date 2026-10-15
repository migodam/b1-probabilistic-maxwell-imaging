# 实际行动排名

## E0：缓存full端点与连续proxy参考

| 旧暴露案例 | repeat1.6 预测风险下降 | 2.4 预测风险下降 | 校准 预测风险下降 |
| --- | --- | --- | --- |
| solid | 0.016637 | 0.032561 | 0.019930 |
| small_cavity | 0.019053 | 0.029510 | 0.022591 |
| large_cavity | 0.014245 | 0.030506 | 0.013998 |
| annular_cavity | 0.045243 | 0.065187 | 0.019764 |
| combined | 0.040926 | 0.063276 | 0.015613 |
| prior_boundary | 0.017355 | 0.027318 | 0.013688 |
| gain_negative_control | 0.033116 | 0.054580 | 0.013884 |

七例的最高材料风险收益都为2.4；raw/Fisher/joint Fisher/EIG/robust相同，conflict都选calibration。预测风险来自有限proxy posterior，不是实测MAE。全部动作含分通道控制的数值见results/E0_ACTION_RANKINGS.json与E0_CRITERION_COMPARISON.json。所有0新增Maxwell。

## E1固定第一个parent、第一noise、第一block的完整数值

这是预先固定首例的展示，未按表现挑选；所有64episode×全部动作均在raw rankings。

| 行动 | raw | robust | conflict | joint Fisher | event EIG† | 预测材料风险下降† |
| --- | --- | --- | --- | --- | --- | --- |
| frequency_1.2 | 23.296486 | 0.039978 | 0.119405 | 0.880918 | 0.194476 | 0.044513 |
| repeat_1.6 | 31.832327 | 0.133678 | 0.000000 | 1.870911 | 0.318041 | 0.063697 |
| frequency_2 | 35.752711 | 0.073071 | 0.098621 | 2.894720 | 0.404171 | 0.062894 |
| frequency_2.4 | 39.124116 | 0.182168 | 0.187036 | 4.777085 | 0.738609 | 0.110420 |
| frequency_2.8 | 56.540814 | 0.455941 | 0.557554 | 6.588345 | 1.053202 | 0.132411 |
| frequency_3.2 | 66.568655 | 0.532758 | 0.895695 | 7.478498 | 1.171405 | 0.143469 |
| reference_calibration | 0.000000 | 0.000000 | 4.300551 | 1.733021 | 0.191696 | 0.036408 |

分数单位/尺度不同，不能跨列比大小。每列内部越大越偏好；risk reduction/EIG属于有限参考†。本例raw、robust、Fisher、EIG、risk均选3.2，conflict选calibration。校准direct material sensitivity为0，但joint Fisher非0，因为约束gain。

E1中fixed3.2与Fisher/risk首步完全同选；E2后续会分开，但没有稳定材料MAE优势。零额外measurement是baseline，停止准则与成本假设见configs/E2_FROZEN.json。

## 真实A交接

confirm_00..03 robust动作依次3.2、2.8、3.2、2.8；reference-assisted全部3.2。每个导入包有完整rankings、新观测与updated_data_history。没有A posterior weights；未把任意等权视为真实posterior。

[主材料结果](MATERIAL_RISK_AND_MAE.md) · [完整库存](E0_ACTION_INVENTORY.json)
