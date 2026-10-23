# P1B 最终冻结统计结果独立审计

状态：**STATISTICAL_RECOMPUTATION_MATCH / SCIENTIFIC_DECISION_NOT_ASSIGNED**。

完成标记首次检查时已经存在，因此等待0秒。初始队列为RUNNING；标记状态为OFFLINE_SCORES_READY_FOR_PARENT_REVIEW，32 observations、16 parents、automatic_GO=false。本审计只复算冻结结果，不运行原main、不推断科学GO。

## 允许读取的证据与复核边界

阶段根目录为 `/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-w-p0r-to-p1/research/research_w_p0r_to_p1`。以下路径均相对该目录。统计复算只读取以下8个文件；初始等待/完成检查另查看了P1B_QUEUE_CHECKPOINT.json状态。文件size/mtime在复算进程内精确核对前后相同；报告写出时再次核对字节数与传输可表示的mtime一致，纳秒元数据最终以字符串保存避免JSON/JavaScript整数舍入，没有计算hash：

- `results/P1B/offline_scores/COMPLETE.json`（171 bytes）。
- `results/P1B/offline_scores/SCORES.json`（337823 bytes）。
- `results/P1B/offline_scores/PARENT_SCORES.json`（109123 bytes）。
- `results/P1B/offline_scores/SUMMARY.json`（26293 bytes）。
- `results/P1B/offline_scores/PAIRED_CI.json`（16097 bytes）。
- `results/P1B/inference/RESULTS.json`（1202803 bytes）。
- `results/P1B/point_validation/RESULTS.json`（399828 bytes）。
- `results/P1B/continuous_point_control/RESULTS.json`（350155 bytes）。

位置、分支真值标签只使用允许读取的冻结SCORES.json内已有评分字段，用于核对损失；未读取OFFLINE_TRUTH、GENERATION_PLAN、观测、规则、原场或任何GT数组。报告不复制逐对象真值标签。本轮不读取生产源码，不导入项目模块；所有复算是独立写出的纯JSON/NumPy统计操作。

## 人口、失败与未合格对象

32个观测在原inference、point_validation、continuous_point_control中全部齐全，parent ID覆盖0–15；每parent对应2个noise，每观测有7臂，SCORES共224条、PARENT_SCORES共112条，均无重复/缺失。全部224条scoring_status=SCORED，64个原point validation记录均FULL_CHECK_PASSED。

当前32/32 continuous reference及32/32 finite reference均qualified，其资格标记在各评分行一致保留；没有通过丢弃unqualified子集获得16-parent均值。不过本批次没有实际unqualified-reference对象，不能从“未发生”外推异常路径已被本批数据测试。此前中止/重试日志不在本轮允许阅读清单，本审计不宣称重新检查了全部历史尝试。

256个SMC seed runs（32观测×2 schedule×4 seed）完整保留于原inference输出。审计发现86个run至少有一项branch/CDF误差超过.02。全体64个ensemble评分对象均保留，这86个超差run对应的评分对象也全部存在；并未因不够精确而删掉该观测。

## 独立复算方法与数值一致性

先分别复算224条位置平方/绝对损失及已有的branch Brier/log、CDF和interval评分，再在每parent内平均两个noise，最后对16个parent平均。对非概率臂，branch/CDF/interval指标原样为null，不补造概率分数。对每一项已保存对照，重建16对parent损失差，独立初始化default_rng(2026102317)，生成(2000,16)配对parent重采样索引，取均值差的线性percentile2.5/97.5分位数。

因此统计单位是**16个parent**，既不是32个noise、256个seed runs，也不是224条arm记录。SMC四seed先在每观测内构成混合后验/均值行动，不能将其另算作独立parent。父代生成的概率独立性来自原实验设计；本轮只检查冻结文件的分组实现，不读取生成计划来重新验证。

共完成5840项有限数值比较，另核对预期null及完整性字段。224条位置损失、112条parent记录的全部指标、7臂全部均值、16个paired区间及所有相关gain/coverage计数均与保存值完全相同；这些主要复算的最大绝对差为**0**。唯一非零读回差为SMC branch质量的roundoff归一化差1.11022302463e-16；没有超过容差的异常。

### 七臂均值（16 parent，noise先平均）

|arm|平方位置损失|绝对位置损失|branch Brier|branch log loss|
|---|---:|---:|---:|---:|
|direct_point|0.00442938219378|0.0487376022699|—|—|
|progressive_point|0.00442916199491|0.0487404648551|—|—|
|direct_SMC_mean|0.00254510752399|0.0473951635064|0.265515983818|0.730480081418|
|progressive_SMC_mean|0.00254660913941|0.0473683723213|0.265619179805|0.730143100661|
|continuous_quadrature_Bayes_mean|0.00255661112534|0.0474942268464|—|—|
|zero_observation_prior_position_mean|0.00236615893363|0.0483023344109|0.25|0.69314718056|
|full_validated_continuous_point_control|0.00404181235048|0.045678575103|—|—|

全部9个summary指标、完整父代覆盖和各项误差上界见REVIEW.json的seven_arm_means、numeric_checks与maximum_recompute_errors。

### 所有已保存paired区间

以下gain为loss(a)-loss(b)，正值有利于b。与保存值逐项复算误差均为0；SUMMARY.contrasts与PAIRED_CI.json完全一致。区间为描述性paired parent percentile区间，未作多重比较校正，本审计不将它自动转为假设或科学门结论。

|family|a → b|loss|gain|95% paired CI|
|---|---|---|---:|---|
|multi_hypothesis|direct_point → direct_SMC_mean|position_squared_error|0.00188427466978|[0.000434178551745, 0.00349843951454]|
|multi_hypothesis|direct_point → direct_SMC_mean|position_absolute_error|0.00134243876345|[-0.0115635558981, 0.0154360184839]|
|multi_hypothesis|progressive_point → progressive_SMC_mean|position_squared_error|0.00188255285549|[0.000433505629858, 0.00350238910309]|
|multi_hypothesis|progressive_point → progressive_SMC_mean|position_absolute_error|0.0013720925338|[-0.0115365893151, 0.0154614125803]|
|progressive_schedule|direct_SMC_mean → progressive_SMC_mean|position_squared_error|-1.50161541889e-06|[-1.78415316377e-05, 1.59221376366e-05]|
|progressive_schedule|direct_SMC_mean → progressive_SMC_mean|position_absolute_error|2.67911850798e-05|[-0.000152001580452, 0.000219545608468]|
|progressive_schedule|direct_point → progressive_point|position_squared_error|2.20198871603e-07|[-2.25978960778e-06, 2.92038622259e-06]|
|progressive_schedule|direct_point → progressive_point|position_absolute_error|-2.86258526288e-06|[-0.00013473760834, 0.000126149852551]|
|observation_information|zero_observation_prior_position_mean → direct_SMC_mean|position_squared_error|-0.000178948590365|[-0.000739362424697, 0.000385324341894]|
|observation_information|zero_observation_prior_position_mean → direct_SMC_mean|position_absolute_error|0.000907170904521|[-0.00501630087884, 0.00736563267421]|
|observation_information|zero_observation_prior_position_mean → continuous_quadrature_Bayes_mean|position_squared_error|-0.000190452191708|[-0.000761380969291, 0.000374075284658]|
|observation_information|zero_observation_prior_position_mean → continuous_quadrature_Bayes_mean|position_absolute_error|0.000808107564496|[-0.00517473251593, 0.00727584996068]|
|strong_point_control|direct_point → full_validated_continuous_point_control|position_squared_error|0.000387569843297|[-8.87960281164e-06, 0.00117438692676]|
|strong_point_control|direct_point → full_validated_continuous_point_control|position_absolute_error|0.00305902716681|[-0.000145478331453, 0.00934943093188]|
|strong_point_control|full_validated_continuous_point_control → direct_SMC_mean|position_squared_error|0.00149670482649|[0.000298337109955, 0.00279413756095]|
|strong_point_control|full_validated_continuous_point_control → direct_SMC_mean|position_absolute_error|-0.00171658840336|[-0.0129602172466, 0.0101420211979]|

就已保存数值而言，零观测prior平方风险为0.00236615893363，direct SMC为0.00254510752399，continuous Bayes mean为0.00255661112534；两项prior-vs-data平方风险区间均跨0。direct/progressive schedule平方风险对照也跨0。如何形成科学表述由父线程决定，不能只列point-vs-SMC对照而隐藏这些对照。

## 强连续点控制：actual Phi选择与KKT范围

32个best均与保存的两个连续branch候选加原已验证grid点的**actual Full Phi最小值**一致，argmin及actual_Phi_improvement复算最大差均为0。30个best来自连续branch候选，2个（observation25、27）保留原validated grid点；改进范围[0, 0.0371566875286]。

该核验使用Phi字段及候选身份，没有使用材料误差或material GT作选择。需要准确限定证据强度：保存结果证明其与actual-Phi选择规则一致；本轮没有读控制实现，并且64个连续候选在32例中的proxy/actual分支排序恰好都一致，不能据此单独验证“若两者逆序时实现仍选择actual”的反事实路径，也不能认证没有记录的隐藏路径。

64个proxy均标记proxy_KKT_only=true，最大proxy_projected_gradient_scaled=7.74935821624e-08。所有best、candidate及顶层actual_Full_KKT_certified均false；32个global_continuous_MAP_certified均false。**proxy KKT不能写成Full KKT或global continuous MAP**。本轮验证的是已有候选的Full-response目标值，不是Full目标导数或全局最优证明。

64个候选actual_full_response标记均true；保存的response_error_white最大0.00613665113942，potential_error最大0.0105219144365。32行control成本字段之和new_rhs=760，这是已有控制实验的来源记录，并非本审计新增费用；本轮不读取账本或宣称等成本。所有control truth_access=false，独立的原始响应未再次打开或求解。

## SMC branch/CDF误差与保留情况

branch_MISS是显著分支完全缺失诊断，与branch质量误差>.02不同。两种schedule的branch_MISS run数均为0，但不能由此写成所有概率数值都合格。

|schedule|runs|branch-vs-finite>.02|CDF-vs-finite>.02|CDF-vs-continuous>.02|最大CDF-vs-continuous|ensemble CDF-vs-continuous>.02|
|---|---:|---:|---:|---:|---:|---:|
|direct|128|17|24|36|0.0355472968215|0|
|progressive|128|18|32|39|0.0385178793397|1|

显著的集成例外为 **observation21 / progressive**：CDF_max_error_vs_continuous=0.0203742640445 >.02；同例CDF_max_error_vs_finite=0.0147619795941，branch误差=0.00814597931266。该ensemble完整留在SCORES和16-parent风险中；reference qualified不等同于该SMC ensemble跨对象精度全合格。

所有256个run的branch/CDF误差均从其保存的质量/CDF与对应finite/continuous reference重新计算，差为0。每观测ensemble的branch、CDF、x_mean从四seed重新平均，差为0；全部64个ensemble的分支质量、CDF误差及MISS率都追踪回评分。86个至少一项超差的run逐项列于REVIEW.json addendum.SMC_error_exceptions.run_records，未静默裁剪或仅挑最好seed。

## Finite、continuous与原先验范围

实际finite grid统一为[2,65,5,21,21]，即650个材料节点×441个共享gain节点；interpolation_used=false且physical_prior_changed=false。有限原子概率与连续受控prior参考分开保存，并未将有限枚举叫作连续精确posterior。

32例finite-vs-continuous分支最大差=0.000207355820543，五个CDF checkpoint最大差=0.00992655677325；continuous5/9分支最大差=0.000333555815837，截断CDF最大差=0.00123821115455。这些只是保存参考及独立复算的数值关系，不是连续PDE或全部CDF位置的严格认证。

SMC interval在源文件中明确属于FINITE_ATOMIC_APPROXIMATION，continuous coverage未被五个CDF checkpoint资格认证。冻结评分的经验覆盖为1不等于完整连续分布的已证明95%校准；本轮只验证已保存区间与评分算术，没有读取原粒子数组来重算分位点。

family_id为W_FIXED_SUPPORT_JOINT_V1，scope是new restricted public fixed-support family, not general unknown-wall online。原minimum_feature=UNRESOLVED；outside_family_mass=null。受控新prior的数值资格不能转写成原任意墙体prior或原minimum_feature条件已qual。

## 需父线程裁决与执行成本

没有发现需修复的算术、人口删除、配对单位或bootstrap实现结果错误。需保留并由父线程解释的具体事项为：一个progressive ensemble的CDF连续参考超差；prior-vs-data与schedule的跨0区间；强point的proxy/Full KKT边界与选择证据范围；新受控prior不能替代原先验资格。无最终科学GO。

本审计两段独立纯分析窗口CPU合计0.149055秒，wall合计1.6716785829秒（NumPy 2.3.5，包含JSON加载/数组统计，排除交互阅读与写报告时间）。PDE/RHS/LU/Book/生产main调用均0；未修改结果、源码、协议、输入或ledger，只写本目录REVIEW.md、REVIEW.json。

生成UTC：2026-10-23T17:19:12.264505+00:00。
