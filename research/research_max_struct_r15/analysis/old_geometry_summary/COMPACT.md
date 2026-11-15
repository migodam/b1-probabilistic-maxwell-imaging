# 旧端点几何：机械离线汇总

范围固定为24个旧案例 × G19/DR两种旧primary端点，共48/48 COMPLETE。只读取原授权CSV、rank_all24实材料JSON及原费用回执；未打开任何预测NPZ、GT数组或public物理数组，未增加physics/优化或更改算法。所有结果为已曝光端点的事后描述，科学解释由父线程完成。

## 定义、分母与适用范围

- h = m_hat − m_coarse_reference；材料误差总能量为 ||h||²。表中份额是实材料正交投影 ||P_b h||²/||h||²，四区可加。这里不是NRMSE，也不是逐坐标支持率。
- real算子以 [Re A; Im A] 实化；先按S的完整右基以tau=.1分强/弱，再在各S区内对G*G压缩矩阵分条件G强/弱。U11/U10为S强内的G强/弱，U01/U00为S弱内的G强/弱；不是全局S强/G强空间的四交集、共同SVD或四个可交换投影。
- 背景S/G strong: sigma/sigma_max >= .1；weak_observable: 1e-12 < sigma/sigma_max < .1；numerical_null: <=1e-12。weak_including_null是后两者之和；背景G weak不同于条件U00。null名称始终按声明数值容差。
- L_hat=I−diag(m_hat)G。有限chord分量为 deltaJ_b=L_hat^-1 diag(P_b h)Q_reference；local分量为 dJ_b=L_hat^-1 diag(P_b h)Q_hat。端口能量是 ||S deltaJ_b||_F² 或 ||S dJ_b||_F²；whitened列再除以原公开标量sigma²。source维按Frobenius求和，未额外除以source数。
- 完整chord分母为 ||S(J_hat−J_reference)||_F²；完整local分母为 ||sum_b S dJ_b||_F²。分别列出原始能量及分母。孤立分量传播后不再正交：signed cross = 完整能量 − 各孤立能量之和 = 2 sum_{b<c} Re<v_b,v_c>。因此传播比率不强迫sum=100%，单分量比率可>1，signed cross可正可负。
- exact_real_Snull是背景实S数值null内的h投影，与上述四区分割独立。S h_null约零不等于材料Jacobian作用 S L_hat^-1 diag(h_null)Q_hat为零，也不等于finite-chord端口响应为零；此处不作材料不可辨识性或恢复保证。
- coarse_reference仅为原反演网格上的Maxwell参考。跨网格观测保留原acquisition mismatch；physical Phi/残差也不构成noise或GT证书。

## real rank元数据

| 组 | 案例数 | U11 | U10 | U01 | U00 | S rank | S null | G strong | G null | global commutator F |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| operator_group_00 | 16 | 11 | 0 | 11 | 554 | 48 | 528 | 19 | 0 | 1.01953 |
| operator_group_01 | 8 | 9 | 0 | 13 | 554 | 40 | 536 | 19 | 0 | 0.692337 |

## 全24案例的材料能量汇总

机械描述：G19/DR的U00材料误差份额逐例中位数为 98.93370451% / 99.08862568%；pooled份额为 44.72028547% / 44.02549751%。两种算法的24个端点均有负的finite-chord S-port signed cross。指定0021、0023的总材料误差能量远高于0009、0020，故本表不以单一pooled数替代逐例数；不从这些描述推出恢复机制或预测器。

pooled_ratio明确是跨24案例的分子能量总和/分母能量总和；case_median是24个案例比率的中位数。大误差案例会影响pooled，因此并列保留，未把48个端点或每个网格点当独立统计检验样本。

| 方法 | 指标 | pooled ratio | case median | min | max |
| --- | --- | --- | --- | --- | --- |
| g19 | U11_material_error_fraction | 0.43394 | 0.00209076 | 0.000478009 | 0.615181 |
| g19 | U10_material_error_fraction | 0 | 0 | 0 | 0 |
| g19 | U01_material_error_fraction | 0.118857 | 0.00719562 | 0.00187044 | 0.208284 |
| g19 | U00_material_error_fraction | 0.447203 | 0.989337 | 0.212967 | 0.997555 |
| g19 | background_S_weak_including_null_error_fraction | 0.56606 | 0.997909 | 0.384819 | 0.999522 |
| g19 | background_S_numerical_null_error_fraction | 0.554485 | 0.982176 | 0.377375 | 0.997484 |
| g19 | background_G_weak_including_null_error_fraction | 0.473383 | 0.991597 | 0.238033 | 0.997853 |
| g19 | Snull_isolated_chord_over_total_chord_port_energy | 70.598 | 42.8597 | 0.631719 | 6379.29 |
| g19 | Snull_isolated_local_over_total_local_port_energy | 0.111766 | 1.29965 | 0.0716976 | 8.79171 |
| g19 | S_port_cross_over_total | -151.111 | -219.448 | -14589.3 | -3.43665 |
| g19 | local_S_port_cross_over_total | -0.0505671 | -2.77738 | -16.1276 | 0.21586 |
| dr | U11_material_error_fraction | 0.42674 | 0.00195281 | 0.000336714 | 0.57683 |
| dr | U10_material_error_fraction | 0 | 0 | 0 | 0 |
| dr | U01_material_error_fraction | 0.133006 | 0.00596627 | 0.00147335 | 0.207524 |
| dr | U00_material_error_fraction | 0.440255 | 0.990886 | 0.215646 | 0.99819 |
| dr | background_S_weak_including_null_error_fraction | 0.57326 | 0.998047 | 0.42317 | 0.999663 |
| dr | background_S_numerical_null_error_fraction | 0.562631 | 0.990714 | 0.422269 | 0.998052 |
| dr | background_G_weak_including_null_error_fraction | 0.468885 | 0.993465 | 0.240633 | 0.998309 |
| dr | Snull_isolated_chord_over_total_chord_port_energy | 255.295 | 72.8172 | 0.612013 | 1013.63 |
| dr | Snull_isolated_local_over_total_local_port_energy | 0.0502279 | 1.42832 | 0.0220532 | 5.57144 |
| dr | S_port_cross_over_total | -477.175 | -230.953 | -3616.07 | -3.84715 |
| dr | local_S_port_cross_over_total | -0.220062 | -1.67985 | -10.0769 | -0.0920417 |

## 指定逐例：0009 / 0020 / 0021 / 0023

下列份额以0–1表示；完整48行、原始分子与分母在 ENDPOINT_COMPACT.csv。

| 案例 | 方法 | ||h||² | U11 | U10 | U01 | U00 | S weak+null | S null | G weak+null |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| r14_new_0009 | g19 | 641.284 | 0.00187912 | 0 | 0.00450543 | 0.993615 | 0.998121 | 0.996289 | 0.992263 |
| r14_new_0009 | dr | 510.919 | 0.00144162 | 0 | 0.00457977 | 0.993979 | 0.998558 | 0.996451 | 0.992741 |
| r14_new_0020 | g19 | 712.786 | 0.00163162 | 0 | 0.00377666 | 0.994592 | 0.998368 | 0.99532 | 0.994732 |
| r14_new_0020 | dr | 703.644 | 0.00111333 | 0 | 0.00418929 | 0.994697 | 0.998887 | 0.99575 | 0.994743 |
| r14_new_0021 | g19 | 38678.3 | 0.615181 | 0 | 0.071548 | 0.313271 | 0.384819 | 0.377375 | 0.331781 |
| r14_new_0021 | dr | 39304.6 | 0.489881 | 0 | 0.0583562 | 0.451763 | 0.510119 | 0.491885 | 0.455559 |
| r14_new_0023 | g19 | 51525.3 | 0.37028 | 0 | 0.163081 | 0.466639 | 0.62972 | 0.616514 | 0.524105 |
| r14_new_0023 | dr | 51638.6 | 0.369621 | 0 | 0.162951 | 0.467428 | 0.630379 | 0.617192 | 0.524787 |

同一指定端点的S-null传播与signed cross（原始端口能量，不white）：

| 案例 | 方法 | ||S h_null||² | null chord energy | total chord denominator | null local energy | total local denominator | chord cross | local cross |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| r14_new_0009 | g19 | 5.90661e-29 | 0.00125365 | 2.35893e-06 | 0.00036642 | 0.00111773 | -0.00335511 | -0.000966734 |
| r14_new_0009 | dr | 7.4801e-29 | 0.000719972 | 2.38022e-06 | 0.000325833 | 0.000998483 | -0.00194954 | -0.000305143 |
| r14_new_0020 | g19 | 7.64424e-29 | 0.000318034 | 2.71558e-06 | 0.00117305 | 0.000431258 | -0.000978893 | -0.00151041 |
| r14_new_0020 | dr | 1.08277e-28 | 0.000250637 | 2.67001e-06 | 0.000562684 | 0.000440033 | -0.000739318 | -0.000732721 |
| r14_new_0021 | g19 | 1.21395e-25 | 2.44938 | 0.0486616 | 3.42706 | 42.8663 | -5.20215 | 9.25312 |
| r14_new_0021 | dr | 9.83517e-26 | 4.34018 | 0.0285641 | 1.39629 | 63.3147 | -8.01863 | -15.8617 |
| r14_new_0023 | g19 | 2.11988e-27 | 2.24571 | 0.00909647 | 1.81011 | 15.0218 | -5.41486 | -1.34357 |
| r14_new_0023 | dr | 2.13924e-27 | 2.24259 | 0.00909589 | 1.80258 | 15.0042 | -5.41319 | -1.38102 |

各band的材料误差、finite-chord/local端口能量（全部192行在 CONDITIONAL_BANDS_WITH_DENOMINATORS.csv；此处固定32行）：

| 案例 | 方法 | band | 材料误差份额 | chord port energy | local port energy |
| --- | --- | --- | --- | --- | --- |
| 0009 | g19 | U11 | 0.00187912 | 0.00134078 | 0.00133237 |
| 0009 | g19 | U10 | 0 | 0 | 0 |
| 0009 | g19 | U01 | 0.00450543 | 0.000535928 | 0.000526387 |
| 0009 | g19 | U00 | 0.993615 | 0.00148076 | 0.000225701 |
| 0009 | dr | U11 | 0.00144162 | 0.000800151 | 0.000790822 |
| 0009 | dr | U10 | 0 | 0 | 0 |
| 0009 | dr | U01 | 0.00457977 | 0.000289833 | 0.000286492 |
| 0009 | dr | U00 | 0.993979 | 0.000861932 | 0.000226312 |
| 0020 | g19 | U11 | 0.00163162 | 0.000261109 | 0.000261687 |
| 0020 | g19 | U10 | 0 | 0 | 0 |
| 0020 | g19 | U01 | 0.00377666 | 0.000297558 | 0.000301005 |
| 0020 | g19 | U00 | 0.994592 | 0.000422941 | 0.00137897 |
| 0020 | dr | U11 | 0.00111333 | 0.000238939 | 0.000239049 |
| 0020 | dr | U10 | 0 | 0 | 0 |
| 0020 | dr | U01 | 0.00418929 | 0.000284463 | 0.000283834 |
| 0020 | dr | U00 | 0.994697 | 0.000218586 | 0.000649871 |
| 0021 | g19 | U11 | 0.615181 | 2.66406 | 26.5339 |
| 0021 | g19 | U10 | 0 | 0 | 0 |
| 0021 | g19 | U01 | 0.071548 | 1.89747 | 5.15794 |
| 0021 | g19 | U00 | 0.313271 | 0.689281 | 1.92135 |
| 0021 | dr | U11 | 0.489881 | 4.31094 | 70.42 |
| 0021 | dr | U10 | 0 | 0 | 0 |
| 0021 | dr | U01 | 0.0583562 | 0.699465 | 2.46858 |
| 0021 | dr | U00 | 0.451763 | 3.03679 | 6.28779 |
| 0023 | g19 | U11 | 0.37028 | 2.42918 | 8.93924 |
| 0023 | g19 | U10 | 0 | 0 | 0 |
| 0023 | g19 | U01 | 0.163081 | 2.39555 | 5.60547 |
| 0023 | g19 | U00 | 0.466639 | 0.599224 | 1.82062 |
| 0023 | dr | U11 | 0.369621 | 2.42577 | 8.93971 |
| 0023 | dr | U10 | 0 | 0 | 0 |
| 0023 | dr | U01 | 0.162951 | 2.39887 | 5.61401 |
| 0023 | dr | U00 | 0.467428 | 0.597643 | 1.83155 |

## 文件与费用

- REAL_FOUR_SECTOR_MATERIAL_ERROR.png/pdf：全24×2四区材料误差堆叠，合法可加的材料投影能量。
- SNULL_TRANSPORTED_PORT_VISIBILITY.png/pdf：全24×2背景S-null份额与finite-chord/local端口能量比；对数纵轴，非信息百分比。
- SIGNED_CROSS_ENERGIES.csv保留384条 current/G/port/whitened、chord/local原始能量和signed cross；BACKGROUND_MATERIAL_BANDS.csv保留288条背景谱分区。
- METHOD_SUMMARY.csv保留每项pooled分子/分母与case median/min/max；REAL_RANK_CASES.csv保留24案例分组和rank定义。
- 原离线端点几何计算费用（不重复计费）：CPU 4.177991000s，原whole-process receipt 4.208414000s，72 LU / 6000 complex RHS。原物理counter子项与总RHS重叠，不能相加。
- 本次新增仅CSV/JSON汇总和静态绘图，0 LU / 0 RHS / 0 SVD / 0优化；实际新增whole-process CPU与wall见 ANALYSIS_RECEIPT.json，外层进程结束后补记OS child CPU。无新科学通过阈值或恢复判断。
