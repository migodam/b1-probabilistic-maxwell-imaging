# R7 独立评分交付：原30对象与已冻结附加诊断

状态：COMPLETE。明确收到父线程 ALL_PRIMARY_SEALED 后才打开新30对象 sealed_truth。父线程公开封存为 audit/PRIMARY_PREDICTION_SEAL.json。主180端点、rational12端点、extra47端点全部实际评分；没有新物理调用，没有改算法、预测、预算、阈值或按GT选替代结果。后续家族后选的新20对象不在本交付范围，不混原30分母。

## 原30对象终点表

|臂|材料成功/30|NRMSE中位数|总CPU秒|总LU|总复RHS|
|---|---:|---:|---:|---:|---:|
|cie_single|21|0.087223|1543.596|54873|1154196|
|cie_dual|20|0.085844|1539.148|53408|1137060|
|global_dual|20|0.086526|1512.969|53923|1137684|
|a_material|21|0.081773|2495.143|54384|1148820|
|a_plus|21|0.082334|1711.336|52156|1103040|
|physical_lbfgs|10|0.314888|838.182|47754|1009764|

single CIE同样21/30，不能宣称A/A+赢过所有成熟基线。A/globalDual两大阳性为obj12、obj20；一个阈值损失为obj09（IoU .696629<.7，而globalDual .712644）。obj17只是在与single比较时的另一阈值边界，不能移植到global比较。

## 原有配对对象bootstrap95%CI

固定10000次对象级成对bootstrap，n=30；六臂不是180独立样本，条件不是额外样本。差值均左减右。

|比较|指标|均值差|95%CI|
|---|---|---:|---|
|a_material vs global_dual|material_nrmse|-0.0438157|[-0.100877, -0.00567822]|
|a_material vs global_dual|material_image_success|0.0333333|[-0.0666667, 0.133333]|
|a_material vs global_dual|cpu_seconds|32.7391|[24.9837, 40.7549]|
|a_plus vs a_material|material_nrmse|0.00060392|[-0.000923286, 0.00287233]|
|a_plus vs a_material|material_image_success|0|[0, 0]|
|a_plus vs a_material|cpu_seconds|-26.1269|[-33.8199, -18.6234]|

A-global的均值误差CI低于0，但配对中位差CI、成功率差CI包含0；改善集中在少数大阳性，不是普遍逐对象优势。A较global更耗CPU。A+与A的21个成功对象完全相同，[0,0]成功差CI只是经验bootstrap无法创造未观察分歧的退化，不能证明未来等价。没有事后非劣性margin。

A+总CPU比A低 31.41%；这支持本批次工程成本降低及观察到的恢复保持，不能将未做的非劣性检验补成通过。

## A+对成熟sourcewiseDual：追加SECONDARY_USER_REQUESTED

|指标|A+减Dual均值|95%CI|
|---|---:|---|
|material_nrmse|-0.0377787|[-0.0929187, -0.00308697]|
|material_image_success|0.0333333|[-0.0666667, 0.133333]|
|cpu_seconds|5.73959|[1.91586, 9.54081]|
|LU|-41.7333|[-207.806, 101.335]|
|actual_complex_RHS|-1134|[-4440.85, 1665.22]|

A+终点总CPU比sourcewiseDual高 11.19%；所以对A的加速不能偷换为对成熟Dual的加速。

成对共同最大CPU/LU/RHS窗口下A+22/30、Dual20/30；平均NRMSE差CI分别[-.09143,-.00269]、[-.09308,-.00376]、[-.09320,-.00381]。对应成功率差CI：CPU[-.0667,.2]，LU/RHS[0,.1667]。仍不得夸成稳健成功率优势。CPU预算早期obj17暂过IoU，LU/RHS预算早期obj09暂过IoU，最终更低Phi可丢失材料通过；预算规则而非GT决定检查点。

## 同网格与跨网格分层

|臂|同网格材料/24|同网格材料+CN99/24|跨网格材料/6|跨网格原始CN99/6|
|---|---:|---:|---:|---:|
|cie_single|18|15|3|3|
|cie_dual|16|14|4|4|
|global_dual|16|14|4|4|
|a_material|17|14|4|4|
|a_plus|17|14|4|4|
|physical_lbfgs|10|9|0|0|

跨36→24的sigma只表述噪声，不能充当含物理失配的完整误差预算。其材料恢复与原始数据/噪声判据分开报告。

24同网格GT可行Phi见证使用已存m/Y_clean及同一prior，无新forward：144同网格端点中27个有更低已知可行GT Phi且报告数值收敛，11个有更低可行Phi且预算停止。六跨网格例固定NOT_AVAILABLE_DIFFERENT_OPERATOR，未把fineY_clean冒充coarse m_ref可行预测。完整逐臂、stop、PG诊断在 GT_FEASIBLE_WITNESS.json。

## 额外诊断

- init18：同六个现有对象，host.1与init1成功数分别Dual4/6→4/6，A5/6→5/6，A+5/6→5/6。两端点均保留，未做GT最佳多初始化替代，不新增n。
- digit9：冻结A+8/9材料成功，中位NRMSE .08324。旧R6 dense结果仅附已公开引用，没有重跑；与新FFT A+存在backend/roundoff/轨迹及readout变化混杂，不能纯因果归因。
- 开发20：四披露对象×五readout。originalA、mobius、guarded、arc各3/4成功，projected2/4。arc总CPU197.78秒vs originalA284.21秒；恢复分数事后评分，不据此重选已冻结A+。
- rational：三旧失败锚点×四已选并校正端点，0/12材料成功。rank2/rank8与继续锚点均未获得合格材料恢复，所有负结果保留。

## 两个新大阳性的完整物理验证

详表 INDEPENDENT_AUDIT.md/json：每个终点576未知材料变量，完整576×576原始Maxwell验证，已计真实LU/所有source和adjoint复RHS。obj12是6source，obj20是12source；不是低秩代价图上的示意结果。

|case|A+NRMSE|A+IoU|A+数据残差|A+线性残差|A+LU|A+复RHS|A+CPU秒|Dual NRMSE/IoU|
|---|---:|---:|---:|---:|---:|---:|---:|---|
|obj12|.0391449|.953125|9.4195e-4|1.51e-15|1745|20916|51.529|.735001/.159763|
|obj20|.0687210|.964286|.00105732|2.52e-15|2005|48060|63.258|.296361/.622093|

single在obj20也成功。direct全部30报告数值收敛（28PG、2相对F），但仅10/30材料成功，低CPU可含很早的坏驻点；不能外推已停止方法的未来表现。

## 核查、成本与图

180主端点与47额外端点均为SEALED且合法，失败分母没有删除。保存Y的残差复算与回执差0；Phi最大差1.11e-16。180终点均有已付费梯度。预算float32材料与最终float64最大差9.52e-7，对NRMSE最大影响1.99e-9；最终快照之后序列化CPU尾部最大.04241秒。

主评分5.419秒wall/4.993秒CPU，主分析1.892/1.833，rational评分.410/.372，extra评分1.069/1.012；全部0新增LU/RHS。后续图排版/secondary统计另有回执，不与物理实验成本混加。

六组主图覆盖全部30×6。逐例视觉检查了全部主图、extra47与rational12。旧digit长标题发生裁切，保留原图并在 FIGURES/primary/readable_extra/ 增补短标题版本；三对成本frontier增补 readable_analysis/；secondary推荐 secondary_aplus_vs_dual_cost_frontiers_readable.png。没有改变任何分数或预测。

## 入口文件

- RESULTS/scoring/primary/SCORES.json, SCORE_RECEIPT.json, GT_FEASIBLE_WITNESS.json, BUDGET_AUDIT.json
- RESULTS/scoring/analysis/ANALYSIS.md, PAIRED_BOOTSTRAP.json, PAIRWISE_COST_FRONTIERS.json, INDEPENDENT_AUDIT.md/json, SECONDARY_APLUS_VS_DUAL.md/json, INFERENCE_LIMITS.md
- RESULTS/scoring/extra/EXTRA_ANALYSIS.md, EXTRA_SCORES.json, SCORE_RECEIPT.json, ANALYSIS_RECEIPT.json
- RESULTS/scoring/rational/RATIONAL_SCORES.json, RATIONAL_COSTS.json, SCORE_RECEIPT.json
- experiments/score_primary.py, analyse_primary.py, score_extra.py, score_rational.py

没有读取或评分后续targeted新20对象，也未对它们提出算法调整。父线程继续负责最终科学解读与发布。
