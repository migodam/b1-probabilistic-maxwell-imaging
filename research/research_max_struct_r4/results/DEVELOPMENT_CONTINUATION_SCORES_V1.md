# Post-confirmation continuation scoring V1

**POST_CONFIRMATION_BUDGET_DIAGNOSTIC / exposed case42cb, not blind.** Both endpoints were sealed before this scoring pass reopened the already exposed GT. Original13/18 score tables and all original confirmation predictions/4-of-6 denominators remain unchanged.

Each assigned endpoint had an allowance of **at most4000** extra original corrector iterations; original FTOL/PGTOL rules could stop it earlier. No added parameter search or further allowance.

| Endpoint | RMSE/12 | IoU | Relative residual | White residual | Original Phi | Joint |
|---|---:|---:|---:|---:|---:|---|
|onehop|0.760862298|0.104761905|0.065583743213|1574.047260|0.00218159599768|FAIL|
|native_extended|0.750890266|0.119402985|0.088671833638|2128.174605|0.00395869725471|FAIL|

| Endpoint | m min / max | Box violation | Extra iterations / allowance | Stop | Projected gradient inf | Original target met |
|---|---:|---:|---:|---|---:|---|
|onehop|0.020000000 / 20.000000000|0|251 / 4000|FTOL relative reduction|3.86901418636e-06|False (1e-8)|
|native_extended|0.020000000 / 20.000000000|0|3361 / 4000|PGTOL projected gradient|7.45734496377e-09|True (1e-8)|

| Endpoint | Original Phi before → after | NRMSE before → after | Extra LU | Extra physical RHS | Wall / CPU s |
|---|---:|---:|---:|---:|---:|
|onehop|0.0021930144605 → 0.00218159599768|0.761367699 → 0.760862298|267|6432|3.504502 / 3.427509|
|native_extended|0.00415147686824 → 0.00395869725471|0.743118035 → 0.750890266|3385|81264|31.140772 / 30.314812|

onehop stops on relative function reduction before attaining the original1e-8 projected-gradient precision; that stop does not qualify KKT accuracy. native_extended attains that numerical bound-projected first-order tolerance. This is a precision statement at a single endpoint, not a local-minimum, uniqueness, or general algorithm qualification. All material/data/image gates remain independently reported above.

Projected gradients and stop reasons come from the existing numerical receipts; the scorer performs no new forward or derivative. Extra-run physical totals use each receipt top-level physical_cost once, including its endpoint diagnostics; nested corrector costs and previous6-LU stationarity diagnostics are not added again. The scoring itself adds0 LU/0 RHS.

No new result replaces any of the original18 confirmation endpoints or secondary choices. This exposed budget diagnostic cannot be counted as independent confirmation. The assigned scientific experiment sequence ends here; final scientific interpretation belongs to parent Codex.
