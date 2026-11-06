# R8成本与发布机械检查（最终已记录成本，保留缺口）

本地最终主预测回执 750/750。仅扫描回执与文件元数据；未读GT，未评分，未运行反演。

任务wall求和是工作量，不能称项目历时。不同硬件、物理模型、矩阵维度分组，不给统一LU等价费用。已知CPU和包含单独标记的失败进程下界；缺口不当零。

|类别|分支|硬件|行数|已知CPU秒|未知CPU行|任务wall和秒|
|---|---|---|---:|---:|---:|---:|
|data_generation|paper_benchmarks|SOURCE_RECEIPT_SCOPE|8|12.607625|0|12.931157|
|data_generation|warm_data|Mac|5|57.956547|0|76.223817|
|data_generation|warm_validation_observations|Mac|1|0.174864|0|0.189323|
|failed|fresnel_engineering|SOURCE_RECEIPT_SCOPE|4|405.520000|1|804.000000|
|failed|paper_benchmarks|SOURCE_RECEIPT_SCOPE|3|320.950000|1|505.957248|
|failed|vector3d_engineering|SOURCE_RECEIPT_SCOPE|1|344.460000|0|0.000000|
|failed|vector3d_serialization|SOURCE_RECEIPT_SCOPE|1|0.000000|1|0.000000|
|failed|warm_remote_dependency_failure|XINAN|8|0.593750|0|1.149570|
|inference|frozen_initializer|SOURCE_RECEIPT_SCOPE|3|0.337369|0|0.453613|
|inference|history_warm_preprocess|SOURCE_RECEIPT_SCOPE|1|0.077402|0|0.142219|
|inference|paper_benchmarks|SOURCE_RECEIPT_SCOPE|2|0.038816|0|0.080571|
|online_inversion|fresnel|SOURCE_RECEIPT_SCOPE|20|7893.046791|0|9667.139365|
|online_inversion|fresnel_air0_development|SOURCE_RECEIPT_SCOPE|3|1886.573247|0|2156.520402|
|online_inversion|fresnel_multifrequency|SOURCE_RECEIPT_SCOPE|2|824.120663|0|1014.098677|
|online_inversion|history_warm|Mac|21|1274.261502|0|1583.139183|
|online_inversion|paper_benchmarks|SOURCE_RECEIPT_SCOPE|59|3214.530186|0|3908.639960|
|online_inversion|vector3d|SOURCE_RECEIPT_SCOPE|52|1199.898052|0|1448.923002|
|online_inversion|vector3d_original_size|SOURCE_RECEIPT_SCOPE|1|379.592774|0|395.867609|
|online_inversion|warm_anchor_confirmation|Mac|66|3927.724018|0|4507.093370|
|online_inversion|warm_anchor_confirmation|XINAN|60|11137.515625|0|11405.216463|
|online_inversion|warm_anchor_development|Mac|8|529.067704|0|721.548892|
|online_inversion|warm_anchor_development_extended|Mac|24|1499.015794|0|1883.124693|
|online_inversion|warm_development|Mac|48|2984.810872|0|3692.829137|
|online_inversion|warm_primary|Mac|405|24744.172129|0|31126.976429|
|online_inversion|warm_primary|XINAN|345|53921.562500|0|55395.747132|
|publication_preparation|anonymous_publication_previous_check|SOURCE_RECEIPT_SCOPE|1|0.144744|0|8.057443|
|publication_preparation|anonymous_publication_verification|SOURCE_RECEIPT_SCOPE|1|0.016693|0|4.013028|
|publication_preparation|fetch_control_transport|Mac|5|0.000000|5|155.183650|
|publication_preparation|final_release_packaging|SOURCE_RECEIPT_SCOPE|1|106.938462|0|131.197553|
|publication_preparation|github_asset_upload|SOURCE_RECEIPT_SCOPE|1|15.974058|0|140.168166|
|publication_preparation|publication_ignore_compaction|SOURCE_RECEIPT_SCOPE|1|0.952206|0|0.997929|
|publication_preparation|publication_ignore_review|SOURCE_RECEIPT_SCOPE|1|0.249919|0|2.607718|
|publication_preparation|remote_incremental_pack|XINAN|5|263.515625|0|593.720487|
|scoring|conditional_deployment_amortization|SOURCE_RECEIPT_SCOPE|1|0.000588|0|0.000897|
|scoring|direct_physics_budget_extension|SOURCE_RECEIPT_SCOPE|1|1.555670|0|1.580075|
|scoring|history_warm_scoring|SOURCE_RECEIPT_SCOPE|1|1.697163|0|2.501586|
|scoring|neural_initializer_endpoint_visualization|SOURCE_RECEIPT_SCOPE|1|3.479344|0|3.638917|
|scoring|warm_anchor_confirmation_analysis|Mac|1|22.480945|0|23.142188|
|scoring|warm_confirmation_tables_figures|Mac|1|31.739264|0|32.456771|
|scoring|warm_initializer_survival_analysis|Mac|1|0.286685|0|0.653295|
|scoring|warm_offline_scoring|Mac|1|10.282360|0|12.176573|
|setup|anchor_remote_transfer|SOURCE_RECEIPT_SCOPE|1|0.000000|1|4.543605|
|setup|hardware_relocation_control|SOURCE_RECEIPT_SCOPE|1|0.000000|1|4.125208|
|setup|history_warm|Mac|21|0.199605|0|0.320874|
|setup|paper_shared_32grid|SOURCE_RECEIPT_SCOPE|4|2.052770|0|2.073108|
|setup|warm_anchor_confirmation|Mac|66|1.009627|0|1.288889|
|setup|warm_anchor_confirmation|XINAN|60|1.296875|0|2.107280|
|setup|warm_anchor_development|Mac|8|0.118714|0|0.241070|
|setup|warm_anchor_development_extended|Mac|24|0.335523|0|0.467767|
|setup|warm_development|Mac|48|0.733645|0|1.017786|
|setup|warm_primary|Mac|405|5.998646|0|8.573847|
|setup|warm_primary|XINAN|345|7.171875|0|9.546239|
|training|neural_training|SOURCE_RECEIPT_SCOPE|3|1406.920852|0|1431.360504|
|validation|deployment_measurement|Mac|2|7.545237|0|12.003913|
|validation|deployment_measurement|XINAN|2|24.093750|0|27.035218|
|validation|final_document_link_QA|SOURCE_RECEIPT_SCOPE|1|0.003273|0|0.007420|
|validation|fresnel_artifact_qa|SOURCE_RECEIPT_SCOPE|2|5.457461|0|0.000000|
|validation|paper_benchmarks|SOURCE_RECEIPT_SCOPE|6|18.262049|0|20.933327|
|validation|public_projection_transport|SOURCE_RECEIPT_SCOPE|2|0.443860|0|0.556263|
|validation|warm_anchor_development_analysis|Mac|1|0.009375|0|0.016943|
|validation|warm_anchor_development_analysis_extended|Mac|1|12.170881|0|14.524229|
|validation|warm_anchor_pilot_plotting|Mac|1|0.000000|1|0.000000|
|validation|warm_validation_analysis|Mac|1|0.329804|0|0.693594|

## 额外开发队列

- warm_anchor_development: 8/8外层回执；DEVELOPMENT_ONLY; not part of750primary; old baseline costs not re-added
- warm_anchor_development_extended: 24/24外层回执；DEVELOPMENT_ONLY; not part of750primary; old baseline costs not re-added
- warm_anchor_confirmation: 126/126外层回执；SEPARATE_FROZEN126_CONFIRMATION; not original750; all complete raw states retained
- history_warm: 21/21外层回执；DEVELOPMENT_ONLY; not part of750primary; old baseline costs not re-added

## 明确缺口

- audit/warm_anchor_xinan_FETCH_HISTORY.jsonl: CPU not metered/recorded; not zero
- audit/warm_anchor_xinan_FETCH_HISTORY.jsonl: CPU not metered/recorded; not zero
- audit/warm_xinan_FETCH_HISTORY.jsonl: CPU not metered/recorded; not zero
- audit/warm_xinan_FETCH_HISTORY.jsonl: CPU not metered/recorded; not zero
- audit/warm_xinan_FETCH_HISTORY.jsonl: CPU not metered/recorded; not zero
- audit/HARDWARE_REALLOCATION.json: CPU not metered/recorded; not zero
- audit/HARDWARE_REALLOCATION_POSTCHECK.json: Postcheck control CPU unmetered; expected final Mac405/66 and XINAN345/60 are allocation targets, not actual completion totals.
- audit/ANCHOR_XINAN_PREPARED.json: CPU not metered/recorded; not zero
- RESULTS/paper_benchmarks/COST_LEDGER.json: CPU not metered/recorded; not zero
- RESULTS/paper_benchmarks/COST_LEDGER.json: Do not sum comparative_algorithm_CPU when reused checkpoints already appear in interruptedattempt processCPU. Longcorrector fees incremental fromearly80. CPUlowerbounds/gaps and sharedoperatorassembly exclusions prevent claiming exactgrandtotal.
- delegated/warm_anchor_analysis/PILOT_RECEIPT.json: CPU not metered/recorded; not zero
- RESULTS/vector3d/ONLINE_COST_LEDGER.json: online successful arms only; sequential CPU can be summed; shared setup/simulation/qualification/original3lambda and failed attempts separately documented; no inferred overhead
- RESULTS/vector3d/RETAINED_FAILURE.json: CPU not metered/recorded; not zero
- RESULTS/fresnel/dielTM_dec8f_cie_dual_multipole_full/INTERRUPTED.json: CPU not metered/recorded; not zero
- delegated/cost_package/PACKAGING_EXECUTION_RECEIPT.json: Packaging timer excludes imports/public plan+gates/preflight/post-receipt index refresh; later independent52-archive verification CPU not separately metered. Recorded compact+ZIP CPU is a paid-phase measurement, not complete orchestration CPU.
- delegated/warm_direct_budget_analysis/EXECUTION_NOTES.md: First syntax/startup failure CPU NOT_METERED; no input data or physics accessed. Successful COST.json is separate; later read-only verification CPU not included/not metered.
- audit/ORCHESTRATION_FAILURES.json: Orchestration/setup/repair overhead is not independently metered; no invented zero.
- audit/XINAN_DETACHED_LAUNCH_FAILURE.json: Orchestration/setup/repair overhead is not independently metered; no invented zero.
- audit/failed_threadpoolctl_XINAN_CONNECTED_RUN.json: Orchestration/setup/repair overhead is not independently metered; no invented zero.
- audit/XINAN_DEPENDENCY_REPAIR.json: Orchestration/setup/repair overhead is not independently metered; no invented zero.
- global: No exact full experiment grand total: shared setup, report/plot/scoring/publication and failed first serialization are partly unmetered.
- global: Money/energy not recorded. No currency or zero-cost claim.
- setup: Standalone public/operator setup not universally metered; warm embedded setup partition is subtracted from online; paper_benchmark shared setup added once per case/timer; vector shared process scopes remain contextual.

## 失败清单

- audit/failed_missing_threadpoolctl/r8_test_0006__a_material__uniform/RECEIPT.json
- audit/failed_missing_threadpoolctl/r8_test_0006__a_plus__uniform/RECEIPT.json
- audit/failed_missing_threadpoolctl/r8_test_0006__a_plus_sourcewise__uniform/RECEIPT.json
- audit/failed_missing_threadpoolctl/r8_test_0006__cie_dual__uniform/RECEIPT.json
- audit/failed_missing_threadpoolctl/r8_test_0006__cie_single__nn/RECEIPT.json
- audit/failed_missing_threadpoolctl/r8_test_0006__cie_single__uniform/RECEIPT.json
- audit/failed_missing_threadpoolctl/r8_test_0006__global_dual__uniform/RECEIPT.json
- audit/failed_missing_threadpoolctl/r8_test_0006__physical_lbfgs__uniform/RECEIPT.json
- RESULTS/vector3d/refinement8_uncached_attempt/ENGINEERING_INTERRUPTION.json
- RESULTS/vector3d/RETAINED_FAILURE.json
- RESULTS/fresnel/FoamTwinDielTM_8GHz_n60_a_plus_v2/INTERRUPTED.json
- RESULTS/fresnel/FoamTwinDielTM_8GHz_n60_cie_dual_v1/INTERRUPTED.json
- RESULTS/fresnel/FoamTwinDielTM_8GHz_n60_cie_dual_v2/INTERRUPTED.json
- RESULTS/fresnel/dielTM_dec8f_cie_dual_multipole_full/INTERRUPTED.json
- audit/ORCHESTRATION_FAILURES.json
- audit/XINAN_DETACHED_LAUNCH_FAILURE.json
- audit/failed_threadpoolctl_XINAN_CONNECTED_RUN.json
- audit/XINAN_DEPENDENCY_REPAIR.json

## 实际归档快照

|去向|文件数|字节|
|---|---:|---:|
|exclude|59|340027068|
|git|18187|87628690|
|release|29823|16334522276|

最终PUBLICATION_INVENTORY来自实际ARCHIVE_INDEX快照，取代早期启发式27包提案。真实52包的成员、尺寸及核对见ARCHIVE_VERIFICATION与PACKAGING_EXECUTION_RECEIPT；上传和匿名访问见audit/PUBLICATION_RECEIPT。评分后合成truth已经进入独立release包。晚增Git文档/发布回执不改变原包。第三方原测量、原Chat、下载全文和本地依赖仍排除。

最终成本状态保留全部未知费用，发布完成不等于科学算法成功。
