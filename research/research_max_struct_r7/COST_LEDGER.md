# R7 完整成本账本

Status: **ALL EXECUTED STUDIES SEALED AND SCORED; PACKAGED; PUBLIC READBACK FOLLOWS**. Primary180, targeted80, initialization18, discloseddigit9, Aplusdevelopment20, three rational/coherent anchors, vector and measured studies completed. No queued scientific work remains. The final machine ledger is [FINAL_COSTS.json](audit/FINAL_COSTS.json); known physical dimensions are [RECORDED_PHYSICS_BY_DIMENSION.json](audit/RECORDED_PHYSICS_BY_DIMENSION.json).

Accounting uses each process's final cumulative physical counters once. Primary, initialization and A-plus development receipts already include auxiliary validation and physical correction: nested counters and timers are not added. The three rational and three coherent top-level study receipts include serial comparator correctors; their nested endpoint receipts and cumulative budget curves are not added as independent expenditures. Those curves cannot establish isolated method costs. A-plus's 16 development and four arc-candidate runs remain paid development expenditure, separate from the primary180; targeted80 and historicaldigit9 are separately charged.

The numerical ledger preserves actual matrix dimensions. Scalar inversion is 576; generation has 24 solves at 576 and six at 1296. The accumulated `full_factorization_dimension` field is a sum over LU calls, not a matrix dimension. Vector-DDA grids have Maxwell dimensions 375 and 648, distinct from their 125 and 216 material coordinates. Scalar divisions, 4-by-4 spatial block work, source-Gram solves, and finite-action systems of dimensions 2, 8 and 32 are separate operation categories; their counts must not be presented as full Maxwell LU counts. The mixed rational query total is left unallocated if per-dimension counts are unavailable.

The initial vector benchmark includes its 12 reconstructions, three data-generation calls, initial derivative qualification and independent CM-coordinate check. Its whole-process timing is counted once; per-method/setup times are not summed again. The independent CM solve adds one 375-dimensional LU and six source RHS outside the adapter counters. The later off-root check is separate. Continuation includes four original 216-DOF reconstructions, one data generation, 16 new continuations and two actual-grid Mie forward checks; `previous_counts` are excluded from continuation deltas. Fractional-cell Mie validation is separate. The inherited R1 Mie receipt is reused evidence and contributes no new R7 work.

Six Fresnel inversion receipts contain disjoint training and held-out prediction counters; their nested auxiliary counters overlap training. Inversion CPU is about 229.436 s in the existing audit, excluding retrieval and later audit work. Two additional host-objective audit evaluations have no complete operation/time receipt; they remain explicitly unmetered. Literature/data retrieval, failed access attempts, figures, reporting, code/audit work and packaging are also not fully metered. Missing measurements are not zero-cost work.

Summed task wall time is paid task duration, not elapsed global wall time. Concurrent queues overlap; their individual intervals are retained without summing them into a global duration. CPU totals are sums of observed process CPU within declared scopes. RSS is a process peak, not additive; only reported whole-process peaks are retained. Missing RSS, CPU, wall, energy and money are null with reasons. No hardware-energy or monetary telemetry was collected.

Evidence packaging is performed after all primary and targeted predictions seal and GT scoring completes; its final inventory/timer/asset sizes are recorded in audit/PACKAGING_MANIFEST.json and audit/PACKAGING_RECEIPT.json. No extra physics is performed. `audit/PUBLICATION_SCOPE.json` defines the release boundary. The packager defaults to path-only inventory and requires an explicit parent seal with `status=ALL_PRIMARY_SEALED` and `release_authorized=true` before packaging. It includes ignored raw generated traces/states, excludes raw Fresnel measurements and measured-Y endpoints, excludes third-party source archives/full texts, and writes individually usable ZIP parts strictly below 2 GB with explicit inclusion/exclusion manifests. Packaging does not itself publish or establish scientific acceptance.

## 已记录的独立开销窗

|Study|Receipts|Observed CPU s|Sum task/process wall s|
|---|---:|---:|---:|
|primary|180|9640.373|10125.391|
|targeted_confirmation|80|4034.064|4148.705|
|initialization|18|1060.281|1097.412|
|historical_digit_regression|9|345.310|352.675|
|aplus_development|20|1420.475|1495.912|
|rational|3|446.094|491.449|
|rational_coherent|3|120.711|137.139|
|vector_initial|12|null:未测|null:未测|
|vector_216_initial|4|null:未测|null:未测|
|vector_continuation|16|null:未测|null:未测|
|vector_data_generation|4|null:未测|null:未测|
|vector_derivative_qualification|2|null:未测|0.018|
|vector_independent_CM_check|1|null:未测|null:未测|
|vector_actual_grid_Mie|2|null:未测|null:未测|
|vector_fractional_Mie|2|null:未测|null:未测|
|vector_process_timing|3|38.150|40.274|
|fresnel_inversion|6|229.436|236.110|
|generation30|30|null:未测|null:未测|
|generation_process_timing|1|3.786|3.881|
|targeted_generation20|20|null:未测|null:未测|
|targeted_generation_process_timing|1|0.159|0.170|
|aplus_wiring_audit|1|null:未测|null:未测|

这些行的operation/timing scope详见JSON。vector process timing已覆盖其计数行，不能把缺失per-run timing当0，也不能重复加入process内部时间。scoring/figure/report/retrieval未全部计时，已测回执另保留；金额、能耗、API账单未知，未虚构固定费用。

## 物理求解的已记录下界，按实际维度

|Complex Maxwell dimension|LU|Complex physical RHS|
|---|---:|---:|
|576|612737|12882540|
|375|4265|49962|
|648|884|10200|
|1296|6|72|

不同维度不可当相同算量；auxiliary Cholesky、FFT、coordinate roots、小resolvent候选与dense full Maxwell分别收费。Fresnel两次host审计没有完整操作回执，因此上表是已记录下界，不是伪造的全任务精确费用。targeted生成前语法失败0物理调用、CPU未测；修复不改科学协议。两次机制绘图布局的首次timer未单独保留，非零计时缺口明确保留。
