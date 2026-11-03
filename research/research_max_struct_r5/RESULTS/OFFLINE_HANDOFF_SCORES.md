# R5 natural-stage handoff：离线评分

两次handoff均来自同一个已公开旧弯条开发对象；独立新盲测对象0。自然分界为原mf3/stage0之后，续跑原mf6/mf10阶段；阶段边界PR重置和Fourier零填充保持原规则。

独立封存：2026-11-03T12:13:54.670469+00:00；本次首次GT读取：2026-11-03T12:13:54.670868+00:00。本次评分新增物理LU/RHS为0，原封存未修改。

两个corrector均达到1000迭代预算上限，未认证收敛。候选和选中seed仅作轨迹参考，不计为额外独立运行；不按GT改选。

| handoff | 状态 | RMSE/12 | IoU | 相对残差 | Phi | 历史联合门 |
|---|---|---:|---:|---:|---:|---|
|global_stage0_to_material_only|auxiliary_candidate: reused_host|0.355104290|0.000000000|1.023252923|0.523523272|FAIL|
|global_stage0_to_material_only|auxiliary_candidate: handoff_initial_original_stage0|0.442158707|0.318725100|0.744002638|0.276774597|FAIL|
|global_stage0_to_material_only|auxiliary_candidate: remaining_original_stage1|0.304449847|0.559440559|0.281629073|0.0396642675|FAIL|
|global_stage0_to_material_only|auxiliary_candidate: remaining_original_stage2|0.217386202|0.750000000|0.167645897|0.014060555|FAIL|
|global_stage0_to_material_only|selected_before_corrector: Frozen original-Phi selection; not GT-selected|0.217386202|0.750000000|0.167645897|0.014060555|FAIL|
|global_stage0_to_material_only|final_prediction: Sealed final original-Phi incumbent|0.388613370|0.414141414|0.107658905|0.00581528123|FAIL|
|material_only_stage0_to_global|auxiliary_candidate: reused_host|0.355104290|0.000000000|1.023252923|0.523523272|FAIL|
|material_only_stage0_to_global|auxiliary_candidate: handoff_initial_original_stage0|0.420694824|0.329218107|0.726518803|0.263919054|FAIL|
|material_only_stage0_to_global|auxiliary_candidate: remaining_original_stage1|0.378421280|0.536912752|0.630330087|0.198665313|FAIL|
|material_only_stage0_to_global|auxiliary_candidate: remaining_original_stage2|0.346448468|0.628099174|0.700267817|0.245196713|FAIL|
|material_only_stage0_to_global|selected_before_corrector: Frozen original-Phi selection; not GT-selected|0.378421280|0.536912752|0.630330087|0.198665313|FAIL|
|material_only_stage0_to_global|final_prediction: Sealed final original-Phi incumbent|0.558737600|0.257812500|0.097692719|0.00479704463|FAIL|

沿用历史RMSE/12≤.15、预测自身Otsu256 IoU≥.7、公共99%噪声球与材料盒，未增加门。GT仅用于封存后描述，禁止在线择根或参数调整。
