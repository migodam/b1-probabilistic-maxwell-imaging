# R10 已记录费用账本

当前已知CPU下界：**42169.780829 秒**。这是成本记录汇总，不是科学结论，也不是包含全部未知开销的精确总额。

|分支|平台|类别|已知CPU秒|未知行|
|---|---|---|---:|---:|
|campaign_bookkeeping|Mac|metadata|0.419486 (1行)|0|
|conditional_confirmation|Mac|control_or_metadata|0.005676 (2行)|0|
|conditional_confirmation|Mac|generation|1.117061 (1行)|0|
|conditional_confirmation|Mac|packaging_transport|1.769202 (7行)|0|
|conditional_confirmation|XINAN|control_or_metadata|1.031250 (5行)|0|
|conditional_confirmation|XINAN|inversion|1336.140625 (9行)|0|
|conditional_confirmation|XINAN|packaging_transport|5.125000 (1行)|0|
|delivery_metadata|UNKNOWN|metadata|0.140309 (1行)|0|
|external_retrieval|UNKNOWN|literature|0.000000 (0行)|1|
|historical_DEV|Mac|DEV|376.015011 (10行)|0|
|historical_DEV|Mac|metadata|0.022816 (1行)|0|
|historical_DEV|UNKNOWN|DEV|120.474626 (8行)|11|
|historical_DEV|UNKNOWN|energy_or_scoring|174.988561 (10行)|0|
|postseal_analysis|UNKNOWN|offline_postprocessing|75.932700 (8行)|0|
|publication_delivery|Mac|packaging_transport|17.690543 (2行)|0|
|publication_packaging|UNKNOWN|packaging_transport|136.994389 (1行)|0|
|right_exposed_diagnostic|Mac|packaging_transport|1.335368 (9行)|0|
|right_exposed_diagnostic|XINAN|control_or_metadata|0.546875 (3行)|0|
|right_exposed_diagnostic|XINAN|inversion|701.640625 (3行)|0|
|right_exposed_diagnostic|XINAN|packaging_transport|0.515625 (1行)|0|
|scoring_energy|UNKNOWN|offline_postprocessing|21.974009 (5行)|1|
|twod_baseline|Mac|generation|0.857409 (1行)|0|
|twod_baseline|Mac|packaging_transport|13.704980 (9行)|0|
|twod_baseline|XINAN|control_or_metadata|7.937500 (2行)|0|
|twod_baseline|XINAN|inversion|19351.046875 (120行)|0|
|twod_baseline|XINAN|packaging_transport|89.062500 (1行)|0|
|twod_matrix_candidates|Mac|packaging_transport|5.921477 (9行)|0|
|twod_matrix_candidates|XINAN|control_or_metadata|4.453125 (2行)|0|
|twod_matrix_candidates|XINAN|inversion|9124.375000 (48行)|0|
|twod_matrix_candidates|XINAN|packaging_transport|36.328125 (1行)|0|
|vector|Mac|packaging_transport|7.139364 (13行)|0|
|vector|Mac|vector_process|4088.168261 (16行)|1|
|vector|XINAN|vector_process|5962.796875 (36行)|0|
|vector_exposed_continuation|Mac|vector_process|1.359581 (5行)|0|
|vector_exposed_continuation|XINAN|vector_process|502.750000 (5行)|0|

whole-child/wait4/GetProcessTimes与controller各计一次；R8/R10内部、auxiliary、corrector、部署视图均不再次收费。原24例数据生成仅在baseline计一次；conditional4为独立新数据，另计一次。两种条件部署共享CIE前缀的费用视图不能相加为生产总账。

打包/传输/归档在表中单列类别。DEV成本审计已包含其历史energy/评分及失败尝试，不能再加其raw receipt或energy总和。vector只读取ACCOUNTING_FINAL（若有）或CURRENT的additive行，compiler/transport已由该producer收费，不再加每臂receipt。

原Maxwell物理LU/RHS、移位辅助矩阵、低秩分母小LU/RHS保留在JSON独立计数视图，维度不合并；fresh Twofold实际rank30，不能沿用历史rank28标签。GPU事件时间不加CPU，且不是active-kernel时间。并行wall不相加为日历工期。

显式unknown记录：52条。缺失/进行中producer见CAMPAIGN_COSTS.json；无记录不填0。

最终刷新：先更新vector与各分支cost ledger，再运行 `python research/research_max_struct_r10/experiments/compile_campaign_costs.py`。本脚本只读取已知cost/setup JSON，不加载GT/NPZ、不运行物理或评分。
