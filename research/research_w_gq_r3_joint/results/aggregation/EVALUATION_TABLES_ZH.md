# 独立评估数据表（自动生成，科学解释由主线程负责）

更新：2026-10-29T18:12:56.828501+00:00。当前有 3 个独立盲测 parent；各 parent 噪声编号为 {'blind_cavity': [0], 'blind_normal': [0, 1], 'blind_strong': [0]}。不能把同一 parent 的两次噪声当成两个独立世界。

主工程门严格为：位置≤0.1、IoU≥0.7、目标实部介电常数相对误差≤0.2，三者同时满足。损耗、f=1复介电常数相对误差、排序后两个连续半轴的绝对/相对误差与图像L2单独报告；排序只移除椭圆参数的轴交换对称，不使用像素分类标签。所有端点先按完整测量目标选择、封存，再打开私有真值评分。

历史ID `blind_cavity` 是连边界空气缺口，不是封闭内腔。非线性前向 F 在各臂均为完整 Maxwell；forcing-current tangent chart 只近似材料切向 J。

| Parent / 历史ID | 噪声 | 方法 | 运行状态 | 位置误差 | IoU | ε实部相对误差 | 损耗绝对误差 | 复ε相对误差f1 | 短半轴相对误差 | 长半轴相对误差 | 材料图相对L2 | 三门合取 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| blind_cavity | 0 | ADAPTIVE_CURRENT_CHART | POST_SEAL_SCORED | 0.01083 | 0.8416 | 0.01279 | 0.1252 | 0.02896 | 0.09845 | 0.01444 | 0.4058 | PASS |
| blind_cavity | 0 | FIXED_CURRENT_CHART | POST_SEAL_SCORED | 0.01753 | 0.7472 | 0.03424 | 0.1248 | 0.04293 | 0.1184 | 0.135 | 0.5079 | PASS |
| blind_cavity | 0 | FULL_GN | POST_SEAL_SCORED | 0.0116 | 0.8407 | 0.006562 | 0.1252 | 0.0268 | 0.09288 | 0.004582 | 0.408 | PASS |
| blind_cavity | 0 | KNOWN_WALL_ORACLE | POST_SEAL_SCORED | 0.01388 | 0.8304 | 0.09481 | 0.1252 | 0.09827 | 0.114 | 0.06196 | 0.4294 | PASS |
| blind_cavity | 0 | PRIOR_ONLY | POST_SEAL_SCORED | 0.2172 | 0.02783 | 0.09017 | 0.0001766 | 0.09014 | 0.01006 | 0.1571 | 1.394 | FAIL |
| blind_cavity | 0 | WALL_FIRST_FIXED | POST_SEAL_SCORED | 0.003988 | 0.7795 | 0.02899 | 0.1252 | 0.03893 | 0.2002 | 0.02874 | 0.4716 | PASS |
| blind_cavity | 1 | ADAPTIVE_CURRENT_CHART | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME |
| blind_cavity | 1 | FIXED_CURRENT_CHART | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME |
| blind_cavity | 1 | FULL_GN | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME |
| blind_normal | 0 | ADAPTIVE_CURRENT_CHART | POST_SEAL_SCORED | 0.03921 | 0.5041 | 0.001679 | 0.03542 | 0.008748 | 0.2827 | 0.3206 | 0.8028 | FAIL |
| blind_normal | 0 | FIXED_CURRENT_CHART | POST_SEAL_SCORED | 0.04368 | 0.4552 | 0.1129 | 0.03542 | 0.1132 | 0.2109 | 0.5104 | 0.8337 | FAIL |
| blind_normal | 0 | FULL_GN | POST_SEAL_SCORED | 0.03917 | 0.5025 | 0.003568 | 0.03542 | 0.009297 | 0.2827 | 0.3155 | 0.8051 | FAIL |
| blind_normal | 0 | KNOWN_WALL_ORACLE | POST_SEAL_SCORED | 0.03763 | 0.5199 | 0.003634 | 0.03542 | 0.009323 | 0.2827 | 0.3341 | 0.7876 | FAIL |
| blind_normal | 0 | PRIOR_ONLY | POST_SEAL_SCORED | 0.5539 | 0 | 0.2728 | 0.08958 | 0.2736 | 0.2808 | 0.2808 | 2.005 | FAIL |
| blind_normal | 0 | WALL_FIRST_FIXED | POST_SEAL_SCORED | 0.03093 | 0.5731 | 0.1202 | 0.03542 | 0.1205 | 0.2827 | 0.4238 | 0.6958 | FAIL |
| blind_normal | 1 | ADAPTIVE_CURRENT_CHART | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME |
| blind_normal | 1 | FIXED_CURRENT_CHART | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME |
| blind_normal | 1 | FULL_GN | POST_SEAL_SCORED | 0.02599 | 0.6183 | 0.01678 | 0.03542 | 0.01884 | 0.2827 | 0.2972 | 0.6687 | FAIL |
| blind_strong | 0 | ADAPTIVE_CURRENT_CHART | POST_SEAL_SCORED | 0.03861 | 0.6492 | 0.2996 | 0.04256 | 0.2996 | 0.1625 | 0.2311 | 0.6002 | FAIL |
| blind_strong | 0 | FIXED_CURRENT_CHART | POST_SEAL_SCORED | 0.1086 | 0.3122 | 0.2252 | 0.2074 | 0.2268 | 0.06883 | 0.1425 | 0.9266 | FAIL |
| blind_strong | 0 | FULL_GN | POST_SEAL_SCORED | 0.03834 | 0.6511 | 0.2982 | 0.04256 | 0.2982 | 0.162 | 0.2302 | 0.5985 | FAIL |
| blind_strong | 0 | KNOWN_WALL_ORACLE | POST_SEAL_SCORED | 0.04874 | 0.5968 | 0.3298 | 0.04256 | 0.3298 | 0.1773 | 0.2634 | 0.6471 | FAIL |
| blind_strong | 0 | PRIOR_ONLY | POST_SEAL_SCORED | 0.2666 | 0 | 0.3311 | 0.08244 | 0.3312 | 0.05656 | 0.0871 | 1.153 | FAIL |
| blind_strong | 0 | WALL_FIRST_FIXED | POST_SEAL_SCORED | 0.07722 | 0.4768 | 0.4603 | 0.1003 | 0.4605 | 0.3042 | 0.3146 | 0.7452 | FAIL |
| blind_strong | 1 | ADAPTIVE_CURRENT_CHART | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME |
| blind_strong | 1 | FIXED_CURRENT_CHART | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME |
| blind_strong | 1 | FULL_GN | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME |

完整目标、收敛与费用：

| Parent | 噪声 | 方法 | 完整J投影梯度 | 活动坐标范围 | 停止原因 | 收敛 | 完整目标 | 数据χ² | FE | RHS | CPU秒 | CPU缺测 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| blind_cavity | 0 | ADAPTIVE_CURRENT_CHART | 0.02926 | all 19 joint coordinates | PROJECTED_GRADIENT | 否 | 210.5 | 417.7 | 765 | 12384 | 366.9 | 否 |
| blind_cavity | 0 | FIXED_CURRENT_CHART | 0.8597 | all 19 joint coordinates | TRUST_REGION_STALLED | 否 | 214.2 | 424.6 | 729 | 10548 | 326.6 | 否 |
| blind_cavity | 0 | FULL_GN | 0.0002628 | all 19 joint coordinates | RELATIVE_OBJECTIVE | 是 | 210.5 | 417.7 | 492 | 8400 | 289.8 | 否 |
| blind_cavity | 0 | KNOWN_WALL_ORACLE | 0.0002486 | seven target coordinates conditional on fixed wall | RELATIVE_OBJECTIVE | 是 | 212.3 | 421.6 | 498 | 8376 | 519.5 | 否 |
| blind_cavity | 0 | PRIOR_ONLY | NOT_RECORDED | not optimized | NOT_RECORDED | NOT_RECORDED | 667 | 1333 | 3 | 12 | 3.901 | 否 |
| blind_cavity | 0 | WALL_FIRST_FIXED | 0.0002341 | seven target coordinates conditional on fixed wall | RELATIVE_OBJECTIVE | 是 | 216.6 | 429.9 | 420 | 6960 | 259.1 | 否 |
| blind_cavity | 1 | ADAPTIVE_CURRENT_CHART | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | 0 | 0 | NOT_RUN_TIME | 否 |
| blind_cavity | 1 | FIXED_CURRENT_CHART | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | 0 | 0 | NOT_RUN_TIME | 否 |
| blind_cavity | 1 | FULL_GN | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | 0 | 0 | NOT_RUN_TIME | 否 |
| blind_normal | 0 | ADAPTIVE_CURRENT_CHART | 0.005553 | all 19 joint coordinates | TRUST_REGION_STALLED | 否 | 165 | 328 | 993 | 17412 | 577 | 否 |
| blind_normal | 0 | FIXED_CURRENT_CHART | 0.5206 | all 19 joint coordinates | RELATIVE_OBJECTIVE | 否 | 166.3 | 330.4 | 393 | 4596 | NOT_RECORDED | 是 |
| blind_normal | 0 | FULL_GN | 0.0001734 | all 19 joint coordinates | PROJECTED_GRADIENT | 是 | 165 | 328 | 726 | 11736 | NOT_RECORDED | 是 |
| blind_normal | 0 | KNOWN_WALL_ORACLE | 6.279e-05 | seven target coordinates conditional on fixed wall | PROJECTED_GRADIENT | 是 | 169.2 | 336.4 | 393 | 6276 | 295.2 | 否 |
| blind_normal | 0 | PRIOR_ONLY | NOT_RECORDED | not optimized | NOT_RECORDED | NOT_RECORDED | 413.6 | 826.4 | 3 | 12 | 3.79 | 否 |
| blind_normal | 0 | WALL_FIRST_FIXED | 6.331e-05 | seven target coordinates conditional on fixed wall | PROJECTED_GRADIENT | 是 | 170.4 | 338.8 | 819 | 13356 | 509.7 | 否 |
| blind_normal | 1 | ADAPTIVE_CURRENT_CHART | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | 0 | 0 | NOT_RUN_TIME | 否 |
| blind_normal | 1 | FIXED_CURRENT_CHART | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | 0 | 0 | NOT_RUN_TIME | 否 |
| blind_normal | 1 | FULL_GN | 0.000511 | all 19 joint coordinates | RELATIVE_OBJECTIVE | 是 | 209.9 | 417.7 | 483 | 7644 | 262.6 | 否 |
| blind_strong | 0 | ADAPTIVE_CURRENT_CHART | 0.01565 | all 19 joint coordinates | RELATIVE_OBJECTIVE | 否 | 184 | 364.8 | 759 | 11868 | 400.5 | 否 |
| blind_strong | 0 | FIXED_CURRENT_CHART | 3.906 | all 19 joint coordinates | TRUST_REGION_STALLED | 否 | 212.7 | 421.8 | 729 | 10740 | 328.5 | 否 |
| blind_strong | 0 | FULL_GN | 0.0002718 | all 19 joint coordinates | RELATIVE_OBJECTIVE | 是 | 184 | 364.8 | 483 | 8220 | 271.7 | 否 |
| blind_strong | 0 | KNOWN_WALL_ORACLE | 0.0006093 | seven target coordinates conditional on fixed wall | RELATIVE_OBJECTIVE | 是 | 191.5 | 381.2 | 414 | 6840 | 438.9 | 否 |
| blind_strong | 0 | PRIOR_ONLY | NOT_RECORDED | not optimized | NOT_RECORDED | NOT_RECORDED | 532.9 | 1065 | 3 | 12 | 3.573 | 否 |
| blind_strong | 0 | WALL_FIRST_FIXED | 0.000346 | seven target coordinates conditional on fixed wall | RELATIVE_OBJECTIVE | 是 | 201.3 | 399 | 528 | 8496 | 308 | 否 |
| blind_strong | 1 | ADAPTIVE_CURRENT_CHART | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | 0 | 0 | NOT_RUN_TIME | 否 |
| blind_strong | 1 | FIXED_CURRENT_CHART | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | 0 | 0 | NOT_RUN_TIME | 否 |
| blind_strong | 1 | FULL_GN | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | NOT_RUN_TIME | 0 | 0 | NOT_RUN_TIME | 否 |

只有三主方法都完成的同一 r0 parent 才进入以下配对描述汇总；当前 parent 集为 ['blind_cavity', 'blind_normal', 'blind_strong']。不作总体置信推断。

| 方法 | 独立parent数 | 位置RMSE | 平均IoU | 平均ε误差 | 合取通过数 | FE | RHS | 完整CPU秒 | CPU缺测行 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| FULL_GN | 3 | 0.03234 | 0.6648 | 0.1028 | 1 | 1701 | 28356 | NOT_RECORDED | 1 |
| FIXED_CURRENT_CHART | 3 | 0.06834 | 0.5049 | 0.1241 | 1 | 1851 | 25884 | NOT_RECORDED | 1 |
| ADAPTIVE_CURRENT_CHART | 3 | 0.03238 | 0.665 | 0.1047 | 1 | 2517 | 41664 | 1344 | 0 |

连续轮廓/分辨率复核（0 PDE，不改冻结门）：

| Parent | 原评分IoU | 细栅格IoU | 中心距/真值等效半径 | 平均边界距λ₀ | 最大边界距λ₀ |
| --- | --- | --- | --- | --- | --- |
| blind_cavity | 0.8407 | 0.8396 | 0.08561 | 0.01148 | 0.02361 |
| blind_normal | 0.5025 | 0.5056 | 0.4013 | 0.03046 | 0.06674 |
| blind_strong | 0.6511 | 0.6498 | 0.2846 | 0.03146 | 0.06353 |

轮廓放大图使用同一0.5λ₀物理视窗和等比例。细化栅格只复核显示/采样，不提高测量信息；FEM独立前向核验另列，不能从高分辨率图推出统计分辨率或连续PDE证书。

开发案例（不计为盲测）：

| Parent / 历史ID | 噪声 | 方法 | 运行状态 | 位置误差 | IoU | ε实部相对误差 | 损耗绝对误差 | 复ε相对误差f1 | 短半轴相对误差 | 长半轴相对误差 | 材料图相对L2 | 三门合取 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| dev | 0 | ADAPTIVE_CURRENT_CHART | POST_SEAL_SCORED | 0.04171 | 0.4463 | 0.05372 | 0.1845 | 0.06846 | 0.2101 | 0.6532 | 0.9862 | FAIL |
| dev | 0 | FIXED_CURRENT_CHART | POST_SEAL_SCORED | 0.05947 | 0.4391 | 0.1998 | 0.1845 | 0.2041 | 0.2101 | 0.6042 | 1.128 | FAIL |
| dev | 0 | FULL_GN | POST_SEAL_SCORED | 0.04162 | 0.4461 | 0.05976 | 0.1845 | 0.07328 | 0.2101 | 0.6404 | 0.9899 | FAIL |

保留的非主开发尝试：

| Parent / 历史ID | 噪声 | 方法 | 运行状态 | 位置误差 | IoU | ε实部相对误差 | 损耗绝对误差 | 复ε相对误差f1 | 短半轴相对误差 | 长半轴相对误差 | 材料图相对L2 | 三门合取 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| dev | 0 | FULL_GN | POST_SEAL_SCORED | 0.04163 | 0.4468 | 0.05981 | 0.1845 | 0.07332 | 0.2101 | 0.6403 | 0.9885 | FAIL |
| dev | 0 | FIXED_CURRENT_CHART | POST_SEAL_SCORED | 0.04753 | 0.4208 | 0.2211 | 0.1845 | 0.225 | 0.21 | 0.4793 | 1.131 | FAIL |

全部 R3 已预约费用（含失败和拒绝，不只是表中端点）：11684 FE / 177928 RHS。主表及控制已记录 CPU 小计 5165.202568 秒；2 行 CPU 缺测，缺失值保持 NOT_RECORDED，不能将小计称为完整 CPU。

四起点的逐分支指标见 BRANCH_RESULTS.csv；所有已封存但未评分项目见 SEALED_NOT_SCORED.csv。未运行或未完成部分不能以现有成功端点替代。没有全局最优、任意未知墙、连续 PDE 或全局 posterior 主张。
