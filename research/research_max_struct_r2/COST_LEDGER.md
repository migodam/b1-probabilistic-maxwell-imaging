# R2 paid cost ledger

Current numerical receipts: **10,261 physical LU / 245,096 complex RHS**; **772 current Cholesky / 9,246 current RHS** are separate.

Numerical process sums: 102.699278s wall / 103.537034s CPU. Measured rendering/scoring-inclusive wall sum: 103.807448s. Known CPU component sum: 104.344143s; full CPU total is unknown because some components were not timed.

|Work|Wall/s|CPU/s|Physical LU|Source RHS|Adjoint RHS|All complex RHS|Current Cholesky|Current RHS|
|---|---:|---:|---:|---:|---:|---:|---:|---:|
|BLIND_GENERATION|0.150064|0.143569|3|36|0|36|0|0|
|PROFILE_GRADIENT_WIRING|0.004372|0.002186|1|6|0|6|3|18|
|obs_1c40f624b9334a6d_direct|0.239663|0.237332|41|492|480|972|0|0|
|obs_1c40f624b9334a6d_global|11.088903|10.958647|1023|12276|12192|24468|0|0|
|obs_1c40f624b9334a6d_native|10.578785|10.501823|1038|12456|12372|24828|0|0|
|obs_1c40f624b9334a6d_profile|8.055382|8.722192|1025|12300|12264|24564|260|3120|
|obs_b20d5a55677d4eca_direct|0.658564|0.647890|107|1284|1272|2556|0|0|
|obs_b20d5a55677d4eca_global|9.220951|9.109286|689|8268|8184|16452|0|0|
|obs_b20d5a55677d4eca_native|11.409796|11.113529|1020|12240|12156|24396|0|0|
|obs_b20d5a55677d4eca_profile|4.552350|5.223185|407|4884|4848|9732|255|3060|
|obs_ef4ed8a3806e4ad6_direct|4.891556|4.833050|856|10272|10260|20532|0|0|
|obs_ef4ed8a3806e4ad6_global|11.130057|10.959073|1031|12372|12288|24660|0|0|
|obs_ef4ed8a3806e4ad6_native|10.968202|10.708919|1025|12300|12216|24516|0|0|
|obs_ef4ed8a3806e4ad6_profile|7.985332|8.644281|1017|12204|12168|24372|254|3048|
|ACTIVE_plan|0.005769|0.005526|2|12|0|12|0|0|
|ACTIVE_acquire|0.003832|0.003813|1|6|0|6|0|0|
|ACTIVE_fit|0.020661|0.020566|10|60|0|60|0|0|
|ACTIVE_score|0.671965|0.786137|2|12|0|12|0|0|
|ACTIVE_coverage|0.003463|0.003434|1|2|0|2|0|0|
|ENDPOINTS|0.295090|0.280210|6|72|72|144|0|0|
|SOURCE_WIRING|0.034905|0.031273|5|30|0|30|0|0|
|DIAGNOSTIC_SUMMARY_ASSEMBLY|0.025044|0.014618|0|0|0|0|0|0|
|PREDICTION_SEAL|0.096409|0.071511|0|0|0|0|0|0|
|POST_CLOSURE_SCORING|0.074902|0.067045|0|0|0|0|0|0|
|BLIND_PLOT|0.707035|0.648707|0|0|0|0|0|0|
|ACTIVE_PLOT|0.199400|null|0|0|0|0|0|0|
|POST_HOC_CAPACITY|10.729617|10.601113|951|11412|11328|22740|0|0|
|POSTHOC_CAPACITY_SCORING|0.005380|0.005228|0|0|0|0|0|0|

Counts sum measured work across different matrix dimensions and decompositions. They are not a FLOP-equivalence certificate. Profile current Cholesky matrices are576×576; the gradient wiring check uses16×16. They are not free, nor are they additional physical Maxwell forward LU.

Outer inverse task, generation, gradient-check, active stage,diagnostic and scoring/rendering timers each added once; nested current/corrector/setup/raw state timers not added again. Process sums are not elapsed research time.

No online post-closure reselection changes the12 main predictions. Post-hoc controls, if included, are labelled development and excluded from blind success rates.

Unknown overhead, incomplete phase partitions and monetary fees remain null. Existing R1 costs are referenced historical work and never added as new R2 executions. See [structured totals](audit/COST_TOTALS.json) and [CSV ledger](audit/COST_LEDGER.csv).
