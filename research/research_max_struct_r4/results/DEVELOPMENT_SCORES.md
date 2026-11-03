# R4 independent development scores

**POST_SEAL_PUBLIC_DEVELOPMENT**: All 13 outputs were closed before this scoring pass. Previous nine-arm scores and scorer remain in `scoring_history/previous_nine/`; original four-arm scores remain in `scoring_history/initial_four/`. Parent authorized already-published R2 truth for scoring only. No new physics, GT-ranked selection, or algorithm change.

Unchanged R2 gates: RMSE/12≤.15; prediction-only256-bin Otsu IoU≥.7; original noisy whitened residual≤25.165422197; public material box.

| Endpoint | NRMSE/12 | Otsu IoU | Relative residual | White residual | Joint |
|---|---:|---:|---:|---:|---|
|R2_native (published reference)|0.510695259|0.245901639|0.100734862|in inherited report|FAIL|
|R3_scaled_original_restart (published reference)|0.510717418|0.254098361|0.100733687|in inherited report|FAIL|
|cooperative_clipped|0.510717418|0.254098361|0.100733687|2417.615|FAIL|
|weak_gn_clipped|0.510717418|0.254098361|0.100733687|2417.615|FAIL|
|box_cooperative|0.510717418|0.254098361|0.100733687|2417.615|FAIL|
|tmatrix_completion|0.524685470|0.250000000|0.095863919|2300.740|FAIL|
|A_block_absorption|0.042653872|0.928571429|0.000967014|23.208|PASS|
|global_common_corrector|0.389182533|0.422680412|0.107820546|2587.700|FAIL|
|C_illuminated_localization|0.383985005|0.082352941|0.636792344|15283.055|FAIL|
|B_lambda0_control|0.376259307|0.292134831|0.378129915|9075.141|FAIL|
|B_lambda1_port|0.533540174|0.244094488|0.093254845|2238.122|FAIL|
|self_absorption|0.389417644|0.422680412|0.107820652|2587.702|FAIL|
|block_absorption_cached|0.041691092|0.928571429|0.000969537|23.269|PASS|
|native_extended|0.557673743|0.255813953|0.097706559|2344.963|FAIL|
|onehop_absorption|0.041672481|0.939759036|0.000967887|23.229|PASS|

T-matrix endpoint versus R3 scaled-original: residual 0.100733687 → 0.095863919; NRMSE 0.510717418 → 0.524685470; IoU 0.254098361 → 0.250000000.
Residual decreased: True; material NRMSE decreased: False; IoU increased: False. Original joint gate: FAIL.

These metrics describe one disclosed development object. They do not establish a new-blind recovery rate, superiority over converged CIE, or a converged optimum. Original R2 blind results remain unchanged. Scoring adds0 LU and0 RHS.

## A versus the common-corrector global control

Both use the mature global CIE mf3/6/10 schedule with500updates per stage, unchanged host/box/prior/source metrics, then the same R4 scaled original corrector with1000iterations. The global control reuses its saved auxiliary endpoint and also retains its historical final incumbent; it does not rerun its auxiliary trajectory. This removes the historical unscaled-corrector mismatch but does not equalize actual work.

A measured new physical work: 1033LU/24708RHS, wall59.4237s. Common-global new work: 1013LU/24288RHS, wall10.6079s. Its reused auxiliary separately recorded 5.1811s and3LU/36RHS; these historical costs are identified, not charged again or called a new end-to-end timing.

A auxiliary material-inner evaluations=61440 versus the stored global control's7500. A also incurred18342720small4x4 factorizations and261872640small complexRHS. These are separate from full Maxwell LU/RHS and prevent an equal-compute or speed claim.

Full receipt-only prior/stopping/cost comparison: [DEVELOPMENT_AUDIT.md](../audit/DEVELOPMENT_AUDIT.md); machine-readable costs: [DEVELOPMENT_COST_AUDIT.json](../audit/DEVELOPMENT_COST_AUDIT.json).
