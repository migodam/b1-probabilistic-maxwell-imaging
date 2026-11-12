# R13 已知费用清单（元数据快照）

这是按host/阶段分列的已知消耗下界，不是最终全费用或并行历时。当前无终态回执的队列保留为不完整；缺失/中断/评分计时为空，不补零。

|Host / 阶段|已知CPU秒|未知CPU记录|
|---|---:|---:|
|i9_remote / online|6449.437500|0|
|i9_remote / online_controller|2.171875|0|
|i9_remote / online_reference|0.000000|0|
|Apple_local / online|2842.706364|0|
|Apple_local / online_controller|0.725367|0|
|Apple_local / online_reference|0.000000|0|
|Apple_local / online_incomplete|0.000000|4|
|i9_remote / online_incomplete|0.000000|4|
|i9_remote / online_failed|0.000000|4|
|Apple_local / validation_failed|0.000000|2|
|Apple_local / setup_generation|5.125994|0|
|Apple_local / setup_NN|1.285002|0|
|UNMEASURED_LEGACY / inherited_prior|0.000000|0|
|Apple_local / setup_packing|3.511289|0|
|Apple_local / validation_exposed|9.146260|0|
|Apple_local / validation_wiring|3.436228|0|
|Apple_local / validation_physics|4.933219|0|
|Apple_local / validation_reference|0.000000|0|
|Apple_local / setup_validation|2.431455|3|
|i9_remote / release_compression|18.203125|0|
|Apple_local / release_download|0.000000|1|
|Apple_local / scoring_plotting|0.463953|5|
|Apple_local / release_site|0.000000|1|

外部child CPU只从权威队列计一次；内部计时、分析COST_TOTALS与恢复子队列为引用。队列wall不与child wall求和。不同矩阵维数的LU/RHS分列，projector/current/Gram动作不能当作物理LU RHS。

## 明确未知或重叠

- `queue:development_v1_recovery`：Merged parent queue already includes these jobs/controller window; reference only
- `controller:development_v1_recovery`：Included in merged parent controller total
- `job:development_v1_recovery:C__r8_test_0001__dr`：Already counted in merged parent
- `job:development_v1_recovery:C__r8_test_0001__ap`：Already counted in merged parent
- `job:development_v1_recovery:C__r8_test_0001__ap_dr`：Already counted in merged parent
- `job:development_v1_recovery:C__r8_test_0001__event_integer`：Already counted in merged parent
- `job:development_v1_recovery:C__r8_test_0001__event_fractional`：Already counted in merged parent
- `job:development_v1_recovery:C__r8_test_0001__random_event`：Already counted in merged parent
- `job:development_v1_recovery:C__r8_test_0001__normal_double`：Already counted in merged parent
- `job:development_v1_recovery:C__r8_test_0001__anderson`：Already counted in merged parent
- `job:development_v1_recovery:A__r10_twod_0001__ap`：Already counted in merged parent
- `job:development_v1_recovery:C__r8_test_0001__reset_fractional`：Already counted in merged parent
- `job:development_v1_recovery:A__r10_twod_0001__dr`：Already counted in merged parent
- `job:development_v1_recovery:A__r10_twod_0001__ap_dr`：Already counted in merged parent
- `job:development_v1_recovery:A__r10_twod_0001__event_integer`：Already counted in merged parent
- `job:development_v1_recovery:A__r10_twod_0001__event_fractional`：Already counted in merged parent
- `job:development_v1_recovery:A__r10_twod_0001__anderson`：Already counted in merged parent
- `job:development_v1_recovery:A__r10_twod_0001__random_event`：Already counted in merged parent
- `job:development_v1_recovery:A__r10_twod_0001__normal_double`：Already counted in merged parent
- `job:development_v1_recovery:A__r10_twod_0001__reset_fractional`：Already counted in merged parent
- `queue:representation_v1`：Queue elapsed is an independent interval; never sum child wall times.
- `internal:representation_v1:A__r10_twod_0004__global_dual`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:A__r10_twod_0004__dual_cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:A__r10_twod_0004__cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:A__r10_twod_0004__a_plus`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:A__r10_twod_0004__direct`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:A__r10_twod_0004__g19`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:A__r10_twod_0004__g19_full_apcs`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:A__r10_twod_0004__g19_full_current`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:B__r12_confirm24_0010__cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:B__r12_confirm24_0010__a_plus`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:B__r12_confirm24_0010__dual_cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:B__r12_confirm24_0010__global_dual`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:B__r12_confirm24_0010__direct`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:B__r12_confirm24_0010__g19`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:B__r12_confirm24_0010__g19_full_apcs`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:B__r12_confirm24_0010__g19_full_current`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:C__r8_test_0001__cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:C__r8_test_0001__dual_cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:C__r8_test_0001__global_dual`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:C__r8_test_0001__a_plus`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:C__r8_test_0001__g19`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:C__r8_test_0001__g19_full_apcs`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:C__r8_test_0001__g19_full_current`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:A__r10_twod_0001__cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:A__r10_twod_0001__dual_cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:C__r8_test_0001__direct`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:A__r10_twod_0001__global_dual`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:A__r10_twod_0001__a_plus`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:A__r10_twod_0001__direct`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:A__r10_twod_0001__g19_full_apcs`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:A__r10_twod_0001__g19_full_current`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:representation_v1:A__r10_twod_0001__g19`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `queue:development_v2`：Queue elapsed is an independent interval; never sum child wall times.
- `internal:development_v2:B__r12_confirm24_0010__random_integer`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v2:A__r10_twod_0004__online_random`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v2:A__r10_twod_0004__random_integer`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v2:A__r10_twod_0004__online_event`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v2:B__r12_confirm24_0010__online_random`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v2:B__r12_confirm24_0010__online_event`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v2:C__r8_test_0001__online_event`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v2:C__r8_test_0001__random_integer`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v2:C__r8_test_0001__online_random`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v2:A__r10_twod_0001__random_integer`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v2:A__r10_twod_0001__online_event`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v2:A__r10_twod_0001__online_random`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `queue:confirmation_local_reference16`：Queue elapsed is an independent interval; never sum child wall times.
- `internal:confirmation_local_reference16:r13_fresh_0000__dr`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0000__dual_cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0000__cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0000__direct`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0000__g19`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0000__physical_globaldual`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0003__cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0003__dual_cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0003__dr`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0003__g19`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0003__physical_globaldual`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0006__cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0006__dual_cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0003__direct`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0006__g19`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0006__dr`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0006__physical_globaldual`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0006__direct`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0009__cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0009__direct`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0009__physical_globaldual`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0009__dual_cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0009__dr`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0009__g19`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0012__cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0012__dual_cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0012__direct`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0012__physical_globaldual`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0012__dr`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0012__g19`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0015__cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0015__dual_cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0015__dr`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0015__physical_globaldual`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0015__direct`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0015__g19`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0018__dual_cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0018__cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0018__physical_globaldual`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0018__direct`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0018__dr`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0018__g19`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0021__cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0021__dual_cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0021__dr`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0021__g19`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0021__direct`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0021__physical_globaldual`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0024__g19`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0024__cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0024__direct`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0024__physical_globaldual`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0024__dual_cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0024__dr`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0027__cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0027__dr`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0027__dual_cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0027__direct`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0027__physical_globaldual`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0027__g19`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0030__dual_cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0030__direct`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0030__cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0030__physical_globaldual`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0030__dr`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0030__g19`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0033__cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0033__dual_cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0033__direct`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0033__physical_globaldual`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0033__dr`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0033__g19`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0036__cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0036__dual_cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0036__direct`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0036__dr`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0036__physical_globaldual`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0036__g19`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:confirmation_local_reference16:r13_fresh_0039__dual_cie`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `running:confirmation_local_reference16:r13_fresh_0039__cie`：Unfinished child exit CPU/physical totals not available; not zero
- `running:confirmation_local_reference16:r13_fresh_0039__g19`：Unfinished child exit CPU/physical totals not available; not zero
- `running:confirmation_local_reference16:r13_fresh_0039__dr`：Unfinished child exit CPU/physical totals not available; not zero
- `running:confirmation_local_reference16:r13_fresh_0039__direct`：Unfinished child exit CPU/physical totals not available; not zero
- `queue:development_v4_local`：Queue elapsed is an independent interval; never sum child wall times.
- `internal:development_v4_local:B__r12_confirm24_0010__online_midpoint`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v4_local:B__r12_confirm24_0010__online_midpoint_random`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `queue:development_v1`：Queue elapsed is an independent interval; never sum child wall times.
- `internal:development_v1:A__r10_twod_0004__ap`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:A__r10_twod_0004__event_integer`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:A__r10_twod_0004__dr`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:A__r10_twod_0004__ap_dr`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:A__r10_twod_0004__event_fractional`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:A__r10_twod_0004__anderson`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:A__r10_twod_0004__random_event`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:A__r10_twod_0004__normal_double`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:B__r12_confirm24_0010__dr`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:A__r10_twod_0004__reset_fractional`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:B__r12_confirm24_0010__ap`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:B__r12_confirm24_0010__ap_dr`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:B__r12_confirm24_0010__event_integer`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:B__r12_confirm24_0010__event_fractional`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:B__r12_confirm24_0010__random_event`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:B__r12_confirm24_0010__normal_double`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:B__r12_confirm24_0010__anderson`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:B__r12_confirm24_0010__reset_fractional`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:C__r8_test_0001__dr`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:C__r8_test_0001__ap`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:C__r8_test_0001__ap_dr`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:C__r8_test_0001__event_integer`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:C__r8_test_0001__event_fractional`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:C__r8_test_0001__random_event`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:C__r8_test_0001__normal_double`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:C__r8_test_0001__anderson`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:A__r10_twod_0001__ap`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:C__r8_test_0001__reset_fractional`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:A__r10_twod_0001__dr`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:A__r10_twod_0001__ap_dr`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:A__r10_twod_0001__event_integer`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:A__r10_twod_0001__event_fractional`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:A__r10_twod_0001__anderson`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:A__r10_twod_0001__random_event`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:A__r10_twod_0001__normal_double`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v1:A__r10_twod_0001__reset_fractional`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `queue:development_v3_local`：Queue elapsed is an independent interval; never sum child wall times.
- `internal:development_v3_local:B__r12_confirm24_0010__online_event`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v3_local:B__r12_confirm24_0010__online_random`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v3_local:B__r12_confirm24_0010__online_double`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v3_local:B__r12_confirm24_0010__online_noevent`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v3_local:B__r12_confirm24_0010__online_reset`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `internal:development_v3_local:B__r12_confirm24_0010__dr`：Internal pre-exit timer/counters overlap authoritative job row; do not add
- `unobserved_queue:confirmation_mechanism16`：Expected176tasks; local ledger cannot certify remote/current completion
- `unobserved_queue:physics_v1`：Expected60tasks; local ledger cannot certify remote/current completion
- `unobserved_queue:confirmation_core48`：Expected432tasks; local ledger cannot certify remote/current completion
- `unobserved_queue:confirmation_online_controls16`：Expected48tasks; local ledger cannot certify remote/current completion
- `analysis_cost_reference:representation_v1`：Same jobs as queue; analysis totals cross-check, not additional work
- `analysis_cost_reference:development_v2`：Same jobs as queue; analysis totals cross-check, not additional work
- `analysis_cost_reference:development_v4_local`：Same jobs as queue; analysis totals cross-check, not additional work
- `analysis_cost_reference:development_v1`：Same jobs as queue; analysis totals cross-check, not additional work
- `analysis_cost_reference:development_v3_local`：Same jobs as queue; analysis totals cross-check, not additional work
- `aborted_C1:C__r8_test_0001__dr`：child terminated on SSH/controller exit before its final receipt; no recoverable process handle
- `aborted_C1:C__r8_test_0001__ap`：child terminated on SSH/controller exit before its final receipt; no recoverable process handle
- `aborted_C1:C__r8_test_0001__ap_dr`：child terminated on SSH/controller exit before its final receipt; no recoverable process handle
- `aborted_C1:C__r8_test_0001__event_integer`：child terminated on SSH/controller exit before its final receipt; no recoverable process handle
- `local_psutil_import_failure`：Original CPU/wall was not captured; no recoverable timing record
- `dual_pressure_failed_serialization`：Not captured before failure; unknown, not imputed from wall time
- `inherited_NN_training`：Historical training CPU/GPU/forward cost unknown; not zero
- `Mie_first_summary_overlap`：Same volume-solve receipts retained unchanged in final combined refinement summary; overlap unresolved, not counted twice
- `new3d_generation_solve:0`：Per-solve CPU not separately timed; included in outer generation interval
- `new3d_generation_solve:1`：Per-solve CPU not separately timed; included in outer generation interval
- `new3d_generation_solve:2`：Per-solve CPU not separately timed; included in outer generation interval
- `new3d_generation_solve:3`：Per-solve CPU not separately timed; included in outer generation interval
- `physical_pack_check:INTEGRATION_CHECK.json`：Outer CPU or full solve counters absent where null; not reconstructed from successful comparison values
- `physical_pack_check:new3d/runner_public/PACKAGING_RECEIPT.json`：Outer CPU or full solve counters absent where null; not reconstructed from successful comparison values
- `physical_pack_check:fresnel_public/PACKAGING_RECEIPT.json`：Outer CPU or full solve counters absent where null; not reconstructed from successful comparison values
- `artifact_download_timing_unknown`：Artifact bytes alone are not download elapsed time; prior downloads not timed
- `scoring_plot:representation_v1`：Scorer/figure CPU and wall were not recorded where null; not zero
- `scoring_plot:development_v2`：Scorer/figure CPU and wall were not recorded where null; not zero
- `scoring_plot:development_v4_local`：Scorer/figure CPU and wall were not recorded where null; not zero
- `scoring_plot:development_v1`：Scorer/figure CPU and wall were not recorded where null; not zero
- `scoring_plot:development_v3_local`：Scorer/figure CPU and wall were not recorded where null; not zero
- `static_site_build`：Builder recorded wall only; historic builder reruns were not separately timed
