# 原生状态证据摘录（仅既有表格汇总）

对象为事后指定的8个既有对象、40个case/method任务；原生分析保留561个snapshot、475个CIE系readout transition、475个current-force行、1208个field-block行。汇总未读取GT/NPZ、未新增forward/adjoint/optimizer/评分。原native analyzer此前已付844次LU、18924 complex RHS；不将这些当成本次新开销。

## 一、回放对齐不是全状态等价

39个已请求raw target中，36个是材料数组逐单元完全相等，3个UNREACHED。另有1任务未请求raw target：r13_fresh_0005__dr，不能把它计为通过或失败。40任务AVAILABLE仅代表轨迹文件可读。

未到达：r13_fresh_0032__g19（1883）、r13_fresh_0032__a_plus（2060）、r13_fresh_0038__g19（1863）。完整原候选ID、iteration、snapshot和源行见RAW_FIDELITY。即使材料逐位相等，也没有核验辅助J/Q、PR历史或corrector端点相等。

## 二、实际readout的物理变化与反事实梯度

475个保留readout上，实际原物理DeltaPhi严格负/零/正为336/0/139。其中反事实onehop loss下降但物理Phi上升有111行；ordinary反事实loss下降但物理Phi上升有114行。这些是保留点逐行描述，不是全部native迭代率，也不是独立样本统计。

`actual_delta_Phi=Phi(m_after)-Phi(m_before)`来自既有同原物理评估；`physical_first_order_step=gPhi_before·(m_after−m_before)`只是局部一阶量。`counterfactual_ordinary/onehop_loss_change`与其梯度在同一个真实postcurrent J/Q处另算，G19轨迹上的ordinary/onehop梯度也不是G19原生材料导数。不能把其cosine当作“实际使用了这一梯度”。

ACTUAL_READOUT_DETAILS逐行保留方向量、物理评估ID与counterfactual值；READOUT_EXTREME_EXAMPLES给每case/method最大上升及最大下降（仅按公开DeltaPhi，非GT）。

## 三、原生PR方向与readout后的current force

`CURRENT_FORCE_DETAILS`的方向是已保存native current方向。before/after force都在同一postcurrent J/Q/A上用beta_current计算，tau=gA_after−gA_before。saved native gradient来自pre-current位置，三者不能混称同一时刻。负gA·direction表示该已存方向对该重算force的一阶下降；并不表示已经执行下一步PR，也不证明下次line search结果。

CASE_METHOD_SUMMARY给严格符号矩阵及normalized dot、gradient rotation、tau/||gA_before||的范围。严格符号计数不设容差；before dot常极接近零，其浮点正负不能直接解释为稳定的下降/上升切换，必须连同具体量级读取。DR没有CIE readout/PR transport行，保持结构性缺失。

## 四、u/v及端口差

`u=J−mQ`，`v=Q−Einc−GJ`，forcing为`u+m v`。current discrepancy为`||J−Jphys(m)||`；port discrepancy为`||S J−Yphys(m)||/delta`；auxiliary residual为`||S J−Y||/delta`，physical residual为`||Yphys(m)−Y||/delta`。`directional_amplification=||J−Jphys||/||u+m v||`是此状态方向上的观察比值，不是resolvent算子范数或普遍上界。

FIELD_BLOCK_SUMMARY按case/method/field block分列min/median/max与最大值源行。CIE fixed_current_before_readout和fixed_current_after_readout使用同J/Q，不同m；DR分别保留shadow/reflected-A/pre-iterator语义，不把不同field block池化。

## 五、native auxiliary history、原选择与corrector

DENSE_HISTORY_SUMMARY总结104754条原native scalar/audit记录：current_step_objective_change是原辅助current步目标变化，并不是DeltaPhi；state_term_current与state_term_material使用的beta可不同。CIE逐迭代只有wall，缺失CPU保持null；不能用wall画成实际CPU。DR sparse shadow gap只对已付审计记录成立。

SELECTION_CORRECTOR_PUBLIC取自此前comparison48/RUN_INDEX.json中原冻结initial/raw/corrected/primary/secondary的候选ID和公开D/R/Phi/ball/CPU，不附新GT分数。raw→corrected、raw→primary等差值只在同原运行内部比较。Replay的paid audit相位数另列，不把计时回放的corrector端点冒充原corrector。所有40任务保留，不按Phi改善筛任务。

## 来源、定位与范围

JSON源行ID为TABLE[零起索引]，对应CSV物理行号index+2（首行为header）。表内source_row_ids或压缩连续range可回到精确输入行；SNAPSHOT字段与physical_evaluation_id进一步定位。`SOURCE_SCHEMA`列完整实际schema。

定义源码：supplement_aplus48/src/analyze_native_states.py:171–193（field defects）、219–264（current force）、546–580（readout/counterfactual）。详细合同见src/README_NATIVE_STATES.md。以上没有显著性检验、成像成功判定、新算法优越性或因果解释；科学判断留给父线程。

|case|method|readout下降/不变/上升|raw target状态|source coverage|
|---|---|---|---|---|
|r13_fresh_0005|cie|13/0/2|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[0]|
|r13_fresh_0005|dual_cie|13/0/2|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[1]|
|r13_fresh_0005|g19|8/0/6|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[2]|
|r13_fresh_0005|a_plus|13/0/2|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[3]|
|r13_fresh_0005|dr|不适用DR|{}|TASK_COVERAGE[4]|
|r13_fresh_0032|cie|8/0/7|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[5]|
|r13_fresh_0032|dual_cie|10/0/5|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[6]|
|r13_fresh_0032|g19|6/0/7|{"UNREACHED": 1}|TASK_COVERAGE[7]|
|r13_fresh_0032|a_plus|11/0/3|{"UNREACHED": 1}|TASK_COVERAGE[8]|
|r13_fresh_0032|dr|不适用DR|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[9]|
|r13_fresh_0038|cie|12/0/3|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[10]|
|r13_fresh_0038|dual_cie|11/0/5|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[11]|
|r13_fresh_0038|g19|10/0/3|{"UNREACHED": 1}|TASK_COVERAGE[12]|
|r13_fresh_0038|a_plus|10/0/5|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[13]|
|r13_fresh_0038|dr|不适用DR|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[14]|
|r13_fresh_0002|cie|12/0/3|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[15]|
|r13_fresh_0002|dual_cie|14/0/1|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[16]|
|r13_fresh_0002|g19|4/0/10|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[17]|
|r13_fresh_0002|a_plus|11/0/4|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[18]|
|r13_fresh_0002|dr|不适用DR|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[19]|
|r13_fresh_0037|cie|11/0/4|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[20]|
|r13_fresh_0037|dual_cie|11/0/4|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[21]|
|r13_fresh_0037|g19|9/0/5|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[22]|
|r13_fresh_0037|a_plus|10/0/5|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[23]|
|r13_fresh_0037|dr|不适用DR|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[24]|
|r13_fresh_0016|cie|7/0/9|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[25]|
|r13_fresh_0016|dual_cie|9/0/7|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[26]|
|r13_fresh_0016|g19|7/0/7|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[27]|
|r13_fresh_0016|a_plus|9/0/6|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[28]|
|r13_fresh_0016|dr|不适用DR|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[29]|
|r13_fresh_0023|cie|11/0/4|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[30]|
|r13_fresh_0023|dual_cie|12/0/3|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[31]|
|r13_fresh_0023|g19|13/0/3|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[32]|
|r13_fresh_0023|a_plus|13/0/3|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[33]|
|r13_fresh_0023|dr|不适用DR|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[34]|
|r13_fresh_0012|cie|10/0/5|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[35]|
|r13_fresh_0012|dual_cie|13/0/2|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[36]|
|r13_fresh_0012|g19|12/0/2|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[37]|
|r13_fresh_0012|a_plus|13/0/2|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[38]|
|r13_fresh_0012|dr|不适用DR|{"EXACT_RECORDED_MATERIAL_MATCH": 1}|TASK_COVERAGE[39]|

## 数值定位示例（按实际DeltaPhi极值，非GT挑选）

|method|规则|case / snapshot|DeltaPhi|物理一阶项|ordinary反事实Delta loss|onehop反事实Delta loss|源行|
|---|---|---|---:|---:|---:|---:|---|
|cie|最小|r13_fresh_0032 / 500|-0.25464445|-0.30099449|-0.00074543067|-0.0010096096|READOUT_TRANSITIONS[64]|
|cie|最大|r13_fresh_0023 / 0|0.76641481|0.33878321|-0.30038957|-0.28478193|READOUT_TRANSITIONS[354]|
|dual_cie|最小|r13_fresh_0032 / 500|-0.23988763|-0.28360276|-0.00047692342|-0.00067995815|READOUT_TRANSITIONS[79]|
|dual_cie|最大|r13_fresh_0023 / 0|0.76641481|0.33878321|-0.30038957|-0.28478193|READOUT_TRANSITIONS[369]|
|g19|最小|r13_fresh_0023 / 7|-0.36894425|-0.44241752|-0.028856525|-0.02887366|READOUT_TRANSITIONS[385]|
|g19|最大|r13_fresh_0023 / 0|0.9344127|0.55888898|-0.26903581|-0.24710787|READOUT_TRANSITIONS[384]|
|a_plus|最小|r13_fresh_0032 / 500|-0.24534514|-0.29929408|-0.00041267252|-0.00063635129|READOUT_TRANSITIONS[107]|
|a_plus|最大|r13_fresh_0023 / 0|0.68118468|0.29158755|-0.29965512|-0.28566729|READOUT_TRANSITIONS[400]|

这些极值仅定位已保存轨迹片段。原生辅助目标、反事实目标和原Phi各列独立；不能把最大/最小行当作代表性分布。
