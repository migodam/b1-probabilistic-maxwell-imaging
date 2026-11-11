# Observed vector readout for parent scientific review

All18 original endpoints plus18 identical longer physical continuations were sealed before first offline GT read. There were no failed attempts or scientific reruns. Original common-corrector max120 stops were retained, followed by a parent-requested common max1000/60CPU continuation, registered before any score. All longer continuations stopped by relative-function reduction within205–471 accepted iterations; this does not certify stationarity or a local minimum.

| Observation | Original120 CIE / Dual / AP / DR NRMSE | Longer common continuation, all six methods |
|---|---|---|
| matched mixed coated cube | .10565 / .10410 / .09252 / .09456 | .09763–.09766 |
| matched connected inclusions | .13614 / .13234 / .08948 / .08146 | .08231–.08243 |
| coated cube generated8³/inverted6³ | .22367 / .22309 / .21922 / .22462 | .22554–.22558 |

The smaller-budget AP/DR differences do not persist after the longer common corrector. On matched cube AP, decreasing physical Phi from.0001282967 to.0001282785 worsened material NRMSE from.09252 to.09763. Cross-grid material bias remains about.2256 across all six methods despite residual about1.03%. Final novelty and claim interpretation belong to parent Codex; these are numerical observations only.

Actual scope:216 real chi unknowns/648 vector field components,60 dipole sources,180 receiver components; real isotropic chi bounds[0,4], epsilon[1,5], actual coated-cube chi≤3.2. Two nonspherical public DEV shapes, three observations,1% declared noise. No fresh independent confirmation, high-epsilon8 test, anisotropic recovery or literal2019 reproduction.

Per-arm initial-stage CPU lower bound totals402.267099s; longer common continuations add87.623642s; observation generation0.541862s, total490.432603s. Full audit-state archive serialization is included in arm timing and can dominate. Initial imports and final metadata writes are outside some near-exit timings; external child CPU is preserved for projection queue and all continuation children, and must not be added again to overlapping arm CPU. All fixed H Cholesky/setup RHS/port eigensystems and physical audit/corrector counts are separate in COST_LEDGER.json. This is fixed declared work, not equal CPU:120 inherited PR cycles versus3000 AP/DR projections. Prefix selection is public-Phi-only and was sealed before scoring; its inner callback clocks exclude later archive serialization, so it is explicitly an exploratory compute-prefix view.

Source/artifacts: run_vector.py, continue_physical.py, compile_report.py; CONTRACT.json, PROJECTION_CONTRACT.json, PREFIX_CONTRACT.json, CONTINUATION_CONTRACT.json; GENERATION_RECEIPT.json; PREDICTION_SEAL.json(36 paths), PREFIX_SELECTION_SEAL.json(108 selections); SCORES.json, LONG_SCORES.json, PREFIX_SCORES.json; REPORT.md; MATERIAL_SLICES.png and MATERIAL_SLICES_LONG.png. The latter was visually inspected and has readable labels, consistent0–4 color scale and all three observation rows. RESULTS/ contains2.2GB of full AP/DR audit states plus compact predictions; preserve states locally and decide public packaging deliberately.
