# R4 six-object confirmation scores

Post-seal independent scoring; six materials, three frozen arms per material. FAILED/UNKNOWN remain in each original six-object denominator. Descriptive counts only; no significance claim.

| Object | Family | Method | Disposition | RMSE/12 | IoU | Relative residual | Joint |
|---|---|---|---|---:|---:|---:|---|
|obs_100bd65498644322|ring|global|SEALED_PREDICTION|0.0733692801|0.924369748|0.000999205072|PASS|
|obs_100bd65498644322|ring|native_extended|SEALED_PREDICTION|0.0859883776|0.875|0.00103760748|PASS|
|obs_100bd65498644322|ring|onehop|SEALED_PREDICTION|0.0716205918|0.917355372|0.00101086421|PASS|
|obs_42cb233e4f224bb3|curved_bar|global|SEALED_PREDICTION|0.741539073|0.121827411|0.0903226157|FAIL/UNAVAILABLE|
|obs_42cb233e4f224bb3|curved_bar|native_extended|SEALED_PREDICTION|0.743118035|0.130653266|0.0908357667|FAIL/UNAVAILABLE|
|obs_42cb233e4f224bb3|curved_bar|onehop|SEALED_PREDICTION|0.761367699|0.110047847|0.0657570571|FAIL/UNAVAILABLE|
|obs_4751796c78134194|two_islands|global|SEALED_PREDICTION|0.0720526094|0.923076923|0.00102421428|PASS|
|obs_4751796c78134194|two_islands|native_extended|SEALED_PREDICTION|0.069851579|0.895522388|0.00102806193|PASS|
|obs_4751796c78134194|two_islands|onehop|SEALED_PREDICTION|0.0734396818|0.923076923|0.00102980336|PASS|
|obs_63abafa750fc483c|curved_bar|global|SEALED_PREDICTION|0.575754689|0.225352113|0.097670798|FAIL/UNAVAILABLE|
|obs_63abafa750fc483c|curved_bar|native_extended|SEALED_PREDICTION|0.577783509|0.251798561|0.0957809891|FAIL/UNAVAILABLE|
|obs_63abafa750fc483c|curved_bar|onehop|SEALED_PREDICTION|0.576876052|0.23943662|0.095880909|FAIL/UNAVAILABLE|
|obs_911efaf2fbac472d|two_islands|global|SEALED_PREDICTION|0.0728496202|0.855072464|0.000996460691|PASS|
|obs_911efaf2fbac472d|two_islands|native_extended|SEALED_PREDICTION|0.0712732109|0.855072464|0.000993287361|PASS|
|obs_911efaf2fbac472d|two_islands|onehop|SEALED_PREDICTION|0.0735597827|0.869565217|0.000996791342|PASS|
|obs_bff1b013e5aa42ae|ring|global|SEALED_PREDICTION|0.0634764306|0.911764706|0.000966834752|PASS|
|obs_bff1b013e5aa42ae|ring|native_extended|SEALED_PREDICTION|0.0673956814|0.889705882|0.000975034866|PASS|
|obs_bff1b013e5aa42ae|ring|onehop|SEALED_PREDICTION|0.066619268|0.926470588|0.000968545456|PASS|

Secondary selection was sealed before GT. Exact-Phi ties are all retained, without GT tiebreak. A partial failed/UNKNOWN pool has unknown full-pool minimum; available endpoints remain visible and costs remain charged.

| Object | Retained available methods | Complete3 | Joint all retained |
|---|---|---|---|
|obs_100bd65498644322|onehop|True|True|
|obs_42cb233e4f224bb3|onehop|True|False|
|obs_4751796c78134194|global|True|True|
|obs_63abafa750fc483c|native_extended|True|False|
|obs_911efaf2fbac472d|global|True|True|
|obs_bff1b013e5aa42ae|onehop|True|True|

| Method | Joint pass / original objects | Unavailable/incomplete |
|---|---:|---:|
|onehop|4 / 6|0|
|global|4 / 6|0|
|native_extended|4 / 6|0|
|secondary_min_original_Phi|4 / 6|0|

## Cost table

Top-level task physical totals only; nested setup/auxiliary/corrector are not added again. Each physical source solve has12 source RHS; small-block RHS are separate. Actual GPU tasks/workers=0 despite legacy runner metadata.

| Method task work | LU | Physical complex RHS | Worker wall sum s | Adapter CPU sum s |
|---|---:|---:|---:|---:|
|onehop|6202|148332|1147.093|1125.20312|
|global|6216|148680|577.952|564.421875|
|native_extended|6229|148848|2174.969|2112.42188|

Actual total is all18 task work plus shared generation/setup **once**. The costed secondary pool is exactly that same total, not an additional run. Unknown costs retain null totals and known subtotals in JSON. Worker wall sums are not concurrent elapsed. Solver import charges are counted once per worker. Peak RSS has worker-lifetime scope.

Measured remote elapsed used for time cost: 1179.03183 s, including foreground SSH/process startup and shutdown; bootstrap cannot be separated. Recorded remote adapter+solver-import CPU: 3804.34375 s; this excludes unmeasured bootstrap. Currency tariff/cost remains unknown.

Shared setup physical counts: `{"cached_operator_imports": 1, "FULL_FORWARD": 6, "actual_complex_RHS": 72, "LU": 6, "source_RHS": 72, "full_factorization_dimension": 3456}`. Campaign physical totals including shared once: `{"LU": {"known_subtotal": 18653.0, "complete": true, "missing_count": 0, "total": 18653.0}, "source_RHS": {"known_subtotal": 223836.0, "complete": true, "missing_count": 0, "total": 223836.0}, "adjoint_RHS": {"known_subtotal": 222096.0, "complete": true, "missing_count": 0, "total": 222096.0}, "actual_complex_RHS": {"known_subtotal": 445932.0, "complete": true, "missing_count": 0, "total": 445932.0}}`.

See JSON task rows for all failures, status pairs, per-task counters, local-operation counts, memory scope, source-RHS counts and incomplete costs. Scoring adds0 LU/0 RHS.
