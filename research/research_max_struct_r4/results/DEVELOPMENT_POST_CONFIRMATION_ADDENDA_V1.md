# Separate post-confirmation development addenda V1

**Not blind; not part of the frozen confirmation denominator.** Both partition-average runs are POST_CONFIRMATION_DEVELOPMENT. Windows onehop is a disclosed original-R2 cross-platform replay. All3 predictions were sealed before this pass reopened their already exposed GT. Original13-row development and18-row confirmation tables remain unchanged; all original method/secondary counts remain4/6.

| Endpoint | Object | RMSE/12 | Otsu256 IoU | Relative residual | White residual | Original Phi | Joint |
|---|---|---:|---:|---:|---:|---:|---|
|partition_average_R2|obs_b20d5a55677d4eca|0.389718611|0.408163265|0.107660170671|2583.850640|0.00581543424693|FAIL|
|partition_average_42cb|obs_42cb233e4f224bb3|0.641395037|0.174418605|0.143490952543|3443.864739|0.010320544208|FAIL|
|onehop_Windows_R2|obs_b20d5a55677d4eca|0.041603262|0.951807229|0.000968350517|23.240471|6.4930492446e-06|PASS|

| Endpoint | Delta NRMSE from original onehop | Delta IoU | Delta relative residual |
|---|---:|---:|---:|
|partition_average_R2|0.348046130|-0.531595771|0.106692283416|
|partition_average_42cb|-0.119972663|0.064370758|0.077733895416|
|onehop_Windows_R2|-0.000069218|0.012048193|0.000000463262|

| Endpoint | Physical LU | Physical complex RHS | Receipt wall s | Receipt CPU s |
|---|---:|---:|---:|---:|
|partition_average_R2|1043|24948|42.490061|40.675342|
|partition_average_42cb|1019|24372|49.207173|47.853363|
|onehop_Windows_R2|1032|24684|376.196616|356.687500|

Each original run top-level physical total is charged once; no nested stage cost is added again. Local process timers and Windows adapter timers have different scopes, retained in JSON; these rows support no cross-platform speed claim. Windows import charge/caller wall are separately recorded. All correctors retain their1000-iteration cap stop. Scoring adds0 LU/0 RHS.

The new partition-average version does not replace frozen onehop, any of the original18 predictions, or the secondary choices. The42cb input was selected by public data failure and smallest opaque ID. Future claims about this revision need independent validation; this disclosed subset supplies no new blind rate. Further4000-step budget endpoints are not included here. Final scientific interpretation belongs to parent Codex.
